"""divide_conquer.py — Enumeração de tours via partição + matching de fronteira.

Estratégia:
  1. Corta o tabuleiro n×n horizontalmente em duas metades A e B.
  2. Enumera todos os caminhos Hamiltonianos abertos DENTRO de cada
     metade, com pontas (entry, exit) sobre a fronteira do corte.
  3. Combina fragmentos de A com fragmentos de B fechando-os por
     duas arestas de fronteira (exit_A→entry_B e exit_B→entry_A).

Ground truth de validação: n=6 deve produzir exatamente 9.862 tours.
"""

import sys
import time
import numpy as np

from .core import (
    State, OK, CONTRADICTION, SUBTOUR, COMPLETE, ACTIVE, FREE,
    KNIGHT_MOVES,
    fix_and_propagate, propagate_initial, uf_find,
)

sys.setrecursionlimit(20000)


# ---------------------------------------------------------------------------
# T1 — Geometria do corte
# ---------------------------------------------------------------------------

def make_cut(n, cut_row=None):
    """Particiona o tabuleiro n×n em A (linhas 0..cut_row-1) e B (linhas cut_row..n-1).

    Retorna dict com:
      'vertices_A', 'vertices_B'     : arrays int32 dos vértices de cada metade
      'boundary_edges'               : array (M, 2) int32, cada linha = (v_A, v_B)
      'boundary_A', 'boundary_B'     : vértices de fronteira (tocam arestas cruzando)
      'cut_row'                      : linha de corte usada
    """
    if cut_row is None:
        cut_row = n // 2
    if not (1 <= cut_row <= n - 1):
        raise ValueError(f"cut_row deve estar em [1, n-1], obteve {cut_row}")

    vertices_A, vertices_B = [], []
    for r in range(n):
        for c in range(n):
            v = r * n + c
            (vertices_A if r < cut_row else vertices_B).append(v)
    vertices_A = np.asarray(vertices_A, dtype=np.int32)
    vertices_B = np.asarray(vertices_B, dtype=np.int32)

    set_A = set(int(v) for v in vertices_A)
    set_B = set(int(v) for v in vertices_B)

    boundary_edges = []
    b_A_set, b_B_set = set(), set()
    for r in range(n):
        for c in range(n):
            u = r * n + c
            for dr, dc in KNIGHT_MOVES:
                rr, cc = r + dr, c + dc
                if not (0 <= rr < n and 0 <= cc < n):
                    continue
                w = rr * n + cc
                if u < w and ((u in set_A and w in set_B) or
                              (u in set_B and w in set_A)):
                    if u in set_A:
                        boundary_edges.append((u, w))
                        b_A_set.add(u)
                        b_B_set.add(w)
                    else:
                        boundary_edges.append((w, u))
                        b_A_set.add(w)
                        b_B_set.add(u)

    boundary_edges = np.asarray(boundary_edges, dtype=np.int32)
    boundary_A = np.asarray(sorted(b_A_set), dtype=np.int32)
    boundary_B = np.asarray(sorted(b_B_set), dtype=np.int32)

    # Verificações
    assert len(vertices_A) + len(vertices_B) == n * n
    for v_A, v_B in boundary_edges:
        assert int(v_A) in set_A and int(v_B) in set_B
    assert set(int(x) for x in boundary_A).issubset(set_A)
    assert set(int(x) for x in boundary_B).issubset(set_B)

    return {
        'n': n,
        'cut_row': cut_row,
        'vertices_A': vertices_A,
        'vertices_B': vertices_B,
        'boundary_edges': boundary_edges,
        'boundary_A': boundary_A,
        'boundary_B': boundary_B,
    }


# ---------------------------------------------------------------------------
# T2 — Grafo restrito a uma metade
# ---------------------------------------------------------------------------

def build_restricted_graph(n, vertices_half):
    """Constrói o grafo do cavalo restrito aos vértices de vertices_half.

    Arestas só são incluídas se AMBOS os endpoints estão em vertices_half
    (arestas cruzando a fronteira ficam de fora).

    Vértices são re-indexados 0..V_half-1.
    Endpoints em índices LOCAIS são armazenados em edge_endpoints;
    endpoints GLOBAIS em edge_endpoints_global (para inspeção/debug).
    """
    half_set = set(int(v) for v in vertices_half)
    vertex_map = {int(v): i for i, v in enumerate(vertices_half)}
    vertex_rmap = np.asarray(vertices_half, dtype=np.int32)
    V = len(vertices_half)

    edges_local, edges_global = [], []
    for v_g in vertices_half:
        v_g = int(v_g)
        r, c = divmod(v_g, n)
        for dr, dc in KNIGHT_MOVES:
            rr, cc = r + dr, c + dc
            if not (0 <= rr < n and 0 <= cc < n):
                continue
            w_g = rr * n + cc
            if w_g not in half_set:
                continue
            if v_g < w_g:
                edges_local.append((vertex_map[v_g], vertex_map[w_g]))
                edges_global.append((v_g, w_g))

    E = len(edges_local)
    edge_endpoints = np.asarray(edges_local, dtype=np.int32) if E > 0 \
                     else np.zeros((0, 2), dtype=np.int32)
    edge_endpoints_global = np.asarray(edges_global, dtype=np.int32) if E > 0 \
                            else np.zeros((0, 2), dtype=np.int32)

    adj_lists = [[] for _ in range(V)]
    for ei, (a, b) in enumerate(edges_local):
        adj_lists[a].append(ei)
        adj_lists[b].append(ei)
    adj_edges = [np.asarray(lst, dtype=np.int32) for lst in adj_lists]
    total_incident = np.asarray([len(lst) for lst in adj_lists], dtype=np.int32)

    return {
        'n': n,
        'V': V,
        'E': E,
        'edge_endpoints': edge_endpoints,
        'edge_endpoints_global': edge_endpoints_global,
        'adj_edges': adj_edges,
        'total_incident': total_incident,
        'vertex_map': vertex_map,
        'vertex_rmap': vertex_rmap,
    }


def is_connected(ctx):
    """BFS para verificar conectividade do grafo restrito."""
    V = ctx['V']
    if V == 0:
        return True
    adj_edges = ctx['adj_edges']
    edge_endpoints = ctx['edge_endpoints']
    seen = np.zeros(V, dtype=bool)
    stack = [0]
    seen[0] = True
    count = 1
    while stack:
        u = stack.pop()
        for ei in adj_edges[u]:
            a = int(edge_endpoints[int(ei), 0])
            b = int(edge_endpoints[int(ei), 1])
            w = b if a == u else a
            if not seen[w]:
                seen[w] = True
                count += 1
                stack.append(w)
    return count == V


# ---------------------------------------------------------------------------
# T3 — Enumeração de fragmentos
# ---------------------------------------------------------------------------

def _vertex_color_global(v_global, n):
    r, c = divmod(int(v_global), n)
    return (r + c) % 2


def _choose_next_edge_simple(state, ctx, rng):
    """Seleção sem heurística f∞ (a fase local do tabuleiro completo não
    aplica a sub-grafos). Usa apenas pressão de vértice."""
    V = ctx['V']
    total_inc = ctx['total_incident']
    adj_edges = ctx['adj_edges']
    target = state.degree_target

    best_score = -(10 ** 9)
    candidates = []
    for v in range(V):
        free_v = int(total_inc[v]) - int(state.degree[v]) - int(state.n_inactive[v])
        if free_v == 0:
            continue
        score = (int(state.degree[v]) * 100 +
                 int(total_inc[v] - state.n_inactive[v] - target[v]))
        if score > best_score:
            best_score = score
            candidates = [v]
        elif score == best_score:
            candidates.append(v)

    if not candidates:
        return -1
    if len(candidates) == 1:
        v_star = candidates[0]
    else:
        v_star = int(candidates[int(rng.integers(len(candidates)))])

    for ei in adj_edges[v_star]:
        if state.fixed[int(ei)] == FREE:
            return int(ei)
    return -1


def _is_complete_path_restricted(state, ctx, start_local):
    V = ctx['V']
    if not np.array_equal(state.degree, state.degree_target):
        return False
    root_s = uf_find(state, start_local)
    return int(state.uf_size[root_s]) == V


def _extract_path_restricted(state, ctx, start_local, end_local):
    V = ctx['V']
    E = ctx['E']
    edge_endpoints = ctx['edge_endpoints']
    nbrs = [[] for _ in range(V)]
    for ei in range(E):
        if state.fixed[ei] == ACTIVE:
            u = int(edge_endpoints[ei, 0])
            w = int(edge_endpoints[ei, 1])
            nbrs[u].append(w)
            nbrs[w].append(u)
    path = np.zeros(V, dtype=np.int32)
    path[0] = start_local
    prev = -1
    cur = start_local
    for i in range(1, V):
        cands = [w for w in nbrs[cur] if w != prev]
        if len(cands) != 1:
            raise RuntimeError(
                f"posição {i}: vértice {cur} tem {len(cands)} candidatos")
        nxt = cands[0]
        path[i] = nxt
        prev = cur
        cur = nxt
    if cur != end_local:
        raise RuntimeError(f"caminho terminou em {cur}, esperava {end_local}")
    return path


def _backtrack_restricted(state, ctx, rng, paths, start_local, end_local,
                          K_max):
    """Enumera caminhos Hamiltonianos start→end no grafo restrito.
    K_max é o teto absoluto (use um número grande para enumeração)."""
    if len(paths) >= K_max:
        return True
    e = _choose_next_edge_simple(state, ctx, rng)
    if e == -1:
        if _is_complete_path_restricted(state, ctx, start_local):
            paths.append(_extract_path_restricted(state, ctx, start_local, end_local))
        return len(paths) >= K_max
    for val in (1, 0):
        snap = state.snapshot()
        status = fix_and_propagate(state, ctx, e, val, mode='path')
        if status == OK:
            if _backtrack_restricted(state, ctx, rng, paths,
                                      start_local, end_local, K_max):
                state.restore(snap)
                return True
        state.restore(snap)
    return len(paths) >= K_max


def knight_path_restricted(ctx, start_local, end_local, K=None, seed=None):
    """Como knight_path mas opera em ctx restrito e suporta K=None (exaustivo)."""
    V = ctx['V']
    if start_local == end_local or not (0 <= start_local < V) or not (0 <= end_local < V):
        return []
    K_max = K if (K is not None and K > 0) else 10 ** 18

    vertex_rmap = ctx['vertex_rmap']
    n = ctx['n']
    s_g = int(vertex_rmap[start_local])
    e_g = int(vertex_rmap[end_local])
    cs = _vertex_color_global(s_g, n)
    ce = _vertex_color_global(e_g, n)
    if V % 2 == 0:
        if cs == ce:
            return []
    else:
        if cs != ce:
            return []

    degree_target = np.full(V, 2, dtype=np.int32)
    degree_target[start_local] = 1
    degree_target[end_local] = 1
    state = State(ctx['V'], ctx['E'], degree_target=degree_target)

    status = propagate_initial(state, ctx, mode='path')
    if status in (CONTRADICTION, SUBTOUR):
        return []
    paths = []
    if status == COMPLETE:
        if _is_complete_path_restricted(state, ctx, start_local):
            paths.append(_extract_path_restricted(state, ctx, start_local, end_local))
        return paths[:K_max]
    rng = np.random.default_rng(seed)
    _backtrack_restricted(state, ctx, rng, paths, start_local, end_local, K_max)
    return paths


def enumerate_fragments(n, cut_info, which='A', seed=None, verbose=False):
    """Enumera todos os fragmentos Hamiltonianos dentro de uma metade.

    Retorna dict {(entry_global, exit_global): [paths em índices globais]}.
    Cada path tem comprimento len(vertices_half).
    Apenas pares (entry, exit) com cores compatíveis e ambos em boundary
    são considerados.
    """
    vertices = cut_info['vertices_A'] if which == 'A' else cut_info['vertices_B']
    boundary = cut_info['boundary_A'] if which == 'A' else cut_info['boundary_B']

    ctx = build_restricted_graph(n, vertices)
    if not is_connected(ctx):
        if verbose:
            print(f"  AVISO: grafo restrito da metade {which} é DESCONEXO")

    vertex_map = ctx['vertex_map']
    vertex_rmap = ctx['vertex_rmap']
    V_half = ctx['V']

    fragments = {}
    profiles_tried = 0
    profiles_with = 0

    boundary_list = list(int(b) for b in boundary)
    for entry_g in boundary_list:
        for exit_g in boundary_list:
            if entry_g == exit_g:
                continue
            # Cores: V par → opostas; V ímpar → iguais
            cs = _vertex_color_global(entry_g, n)
            ce = _vertex_color_global(exit_g, n)
            if V_half % 2 == 0:
                if cs == ce:
                    continue
            else:
                if cs != ce:
                    continue
            profiles_tried += 1
            entry_l = vertex_map[entry_g]
            exit_l = vertex_map[exit_g]
            paths = knight_path_restricted(ctx, entry_l, exit_l, K=None, seed=seed)
            if paths:
                profiles_with += 1
                fragments[(entry_g, exit_g)] = [vertex_rmap[p].astype(np.int32)
                                                 for p in paths]

    if verbose:
        n_frags = sum(len(v) for v in fragments.values())
        if fragments:
            counts = [len(v) for v in fragments.values()]
            print(f"  Metade {which}: {profiles_tried} perfis tentados, "
                  f"{profiles_with} com fragmentos, total {n_frags} fragmentos "
                  f"(min={min(counts)} max={max(counts)} média={n_frags/len(counts):.1f})")
        else:
            print(f"  Metade {which}: {profiles_tried} perfis tentados, NENHUM com fragmentos")

    return fragments


# ---------------------------------------------------------------------------
# T4 — Matching de fronteira
# ---------------------------------------------------------------------------

def _assemble_tour(frag_A, frag_B):
    """Concatena frag_A (entry_A → ... → exit_A) + frag_B (entry_B → ... → exit_B).
    O tour é fechado: a aresta exit_A→entry_B fecha A↔B e exit_B→entry_A fecha
    o ciclo (verificadas externamente)."""
    return np.concatenate([frag_A, frag_B])


def match_fragments(fragments_A, fragments_B, cut_info, n,
                    count_only=False, verbose=False):
    """Combina fragmentos compatíveis para formar tours fechados.

    Compatibilidade do par (frag_A, frag_B):
      (exit_A, entry_B) ∈ boundary_edges  AND  (exit_B, entry_A) ∈ boundary_edges
    """
    boundary_set = set()
    for v_A, v_B in cut_info['boundary_edges']:
        boundary_set.add((int(v_A), int(v_B)))
        boundary_set.add((int(v_B), int(v_A)))

    count = 0
    tours = []

    # Indexar fragments_B por (entry_B, exit_B) para iteração direta
    for (entry_A, exit_A), frags_A_list in fragments_A.items():
        for (entry_B, exit_B), frags_B_list in fragments_B.items():
            if (exit_A, entry_B) not in boundary_set:
                continue
            if (exit_B, entry_A) not in boundary_set:
                continue
            n_pair = len(frags_A_list) * len(frags_B_list)
            count += n_pair
            if not count_only:
                for fA in frags_A_list:
                    for fB in frags_B_list:
                        tours.append(_assemble_tour(fA, fB))

    if verbose:
        print(f"  Matching: {count} tours encontrados "
              f"({'count_only' if count_only else f'{len(tours)} montados'})")

    return count, tours


# ---------------------------------------------------------------------------
# T5 — Interface principal
# ---------------------------------------------------------------------------

def knight_tours_dc(n, cut_row=None, count_only=False, seed=None, verbose=True):
    """Enumera tours fechados do cavalo n×n via divide-and-conquer.

    Args:
      n          : tamanho do tabuleiro (≥ 6)
      cut_row    : linha de corte (default: n//2)
      count_only : retorna apenas contagem
      seed       : semente para knight_path
      verbose    : imprimir progresso

    Returns:
      int (count_only=True) ou list[np.array] (count_only=False)
    """
    if n < 6:
        raise ValueError("n deve ser ≥ 6")

    if verbose:
        print(f"Divide-and-conquer n={n}, cut_row={cut_row or n//2}")

    t0 = time.perf_counter()
    cut = make_cut(n, cut_row)
    t_cut = time.perf_counter() - t0
    if verbose:
        print(f"  Corte: |A|={len(cut['vertices_A'])} |B|={len(cut['vertices_B'])} "
              f"boundary_edges={len(cut['boundary_edges'])} "
              f"|bA|={len(cut['boundary_A'])} |bB|={len(cut['boundary_B'])} "
              f"({t_cut:.4f}s)")

    t0 = time.perf_counter()
    frags_A = enumerate_fragments(n, cut, which='A', seed=seed, verbose=verbose)
    t_A = time.perf_counter() - t0
    if verbose:
        n_fa = sum(len(v) for v in frags_A.values())
        print(f"  Metade A: {n_fa} fragmentos em {t_A:.3f}s")

    t0 = time.perf_counter()
    frags_B = enumerate_fragments(n, cut, which='B', seed=seed, verbose=verbose)
    t_B = time.perf_counter() - t0
    if verbose:
        n_fb = sum(len(v) for v in frags_B.values())
        print(f"  Metade B: {n_fb} fragmentos em {t_B:.3f}s")

    t0 = time.perf_counter()
    count, tours = match_fragments(frags_A, frags_B, cut, n,
                                    count_only=count_only, verbose=verbose)
    t_match = time.perf_counter() - t0
    if verbose:
        print(f"  Matching: {count} tours em {t_match:.3f}s")
        print(f"  Total: {t_cut + t_A + t_B + t_match:.3f}s")

    return count if count_only else tours


# ---------------------------------------------------------------------------
# T6 — Validação e benchmark
# ---------------------------------------------------------------------------

def validate_dc():
    print("=== Validação divide-and-conquer ===\n")
    from .tours import verify_tour

    # n=6: ground truth
    count = knight_tours_dc(6, count_only=True)
    print(f"\n→ n=6: {count} tours (GT=9862)")
    if count != 9862:
        print(f"  ✗ ERRO: esperava 9862, obteve {count}")
        return False
    print("  ✓ Contagem correta")

    tours = knight_tours_dc(6, count_only=False, verbose=False)
    if len(tours) != 9862:
        print(f"  ✗ ERRO: lista tem {len(tours)} elementos, esperava 9862")
        return False
    n_check = min(200, len(tours))
    bad = [i for i in range(n_check) if not verify_tour(tours[i], 6)]
    if bad:
        print(f"  ✗ ERRO: {len(bad)}/{n_check} tours inválidos (primeiros: {bad[:5]})")
        return False
    print(f"  ✓ {n_check} primeiros tours verificados (knight moves + sem repetição)")

    return True


def benchmark_dc():
    print("\n=== Benchmark DC vs DFS direto ===\n")
    from .tours import knight_tours

    for n in (6,):
        print(f"--- n={n} ---")
        t0 = time.perf_counter()
        tours_dfs = knight_tours(n, 9999999, seed=0)
        t_dfs = time.perf_counter() - t0
        print(f"  DFS direto:           {len(tours_dfs)} tours em {t_dfs:.3f}s")

        t0 = time.perf_counter()
        count_dc = knight_tours_dc(n, count_only=True, verbose=False)
        t_dc = time.perf_counter() - t0
        print(f"  Divide-and-conquer:   {count_dc} tours em {t_dc:.3f}s")
        print(f"  Speedup DC/DFS:       {t_dfs/max(t_dc,1e-9):.2f}×")
        print()


if __name__ == '__main__':
    ok = validate_dc()
    if ok:
        benchmark_dc()
