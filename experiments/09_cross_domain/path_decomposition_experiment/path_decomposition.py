"""
Path-finding via DFS loop decomposition (XOR) vs. heurísticas clássicas.

Hipótese: o conjunto de caminhos hamiltonianos s->t pode ser enumerado
gerando-se um caminho-base e aplicando XOR (diferença simétrica de arestas)
com combinações de ciclos fundamentais do grafo.

Grafo de teste: passeio do cavalo em tabuleiros 6x6 e 8x8.
Resultados reportados honestamente, favoráveis ou não à hipótese.
"""

from __future__ import annotations

import itertools
import json
import random
import time
from pathlib import Path

import networkx as nx

RESULTS_DIR = Path(__file__).parent / "results"
RESULTS_DIR.mkdir(exist_ok=True)

# Representação de aresta: frozenset({u, v}) com u, v índices de vértice.
Edge = frozenset


# ----------------------------------------------------------------------------
# STEP 1 — Construção do grafo
# ----------------------------------------------------------------------------
KNIGHT_MOVES = [(1, 2), (2, 1), (-1, 2), (-2, 1),
                (1, -2), (2, -1), (-1, -2), (-2, -1)]


def build_knight_graph(n: int) -> nx.Graph:
    """Grafo do cavalo n x n. Vértice (r, c) -> índice r*n + c."""
    G = nx.Graph()
    G.add_nodes_from(range(n * n))
    for r in range(n):
        for c in range(n):
            u = r * n + c
            for dr, dc in KNIGHT_MOVES:
                rr, cc = r + dr, c + dc
                if 0 <= rr < n and 0 <= cc < n:
                    v = rr * n + cc
                    if u < v:
                        G.add_edge(u, v)
    return G


def color(idx: int, n: int) -> int:
    return ((idx // n) + (idx % n)) & 1


def feasible_parity(s: int, t: int, n: int) -> bool:
    """Caminho hamiltoniano tem n*n-1 arestas; paridade de cor exige
    extremos de cores opostas sse (n*n-1) é ímpar (sempre, p/ n>=1)."""
    n_edges = n * n - 1
    needs_opposite = (n_edges % 2 == 1)
    same = (color(s, n) == color(t, n))
    return (not needs_opposite) == same  # opostas exigidas <-> cores diferentes


# ----------------------------------------------------------------------------
# STEP 2 — Caminho-base via Warnsdorff
# ----------------------------------------------------------------------------
def warnsdorff_path(G: nx.Graph, s: int, t: int, n_vertices: int,
                    rng: random.Random | None = None,
                    tie_random: bool = False) -> list[int] | None:
    """Heurística de Warnsdorff a partir de s; aceita só se cobrir todos os
    vértices e terminar exatamente em t."""
    rng = rng or random
    visited = {s}
    path = [s]
    cur = s
    while len(path) < n_vertices:
        nbrs = [w for w in G.neighbors(cur) if w not in visited]
        if not nbrs:
            break
        # grau de Warnsdorff: menor nº de vizinhos não-visitados (excl. cur).
        def deg(w):
            return sum(1 for x in G.neighbors(w) if x not in visited and x != cur)
        best = min(deg(w) for w in nbrs)
        cands = [w for w in nbrs if deg(w) == best]
        cur = rng.choice(cands) if (tie_random and len(cands) > 1) else cands[0]
        visited.add(cur)
        path.append(cur)
    if len(path) == n_vertices and path[-1] == t:
        return path
    return None


def find_base_path(G: nx.Graph, s: int, t: int, n_vertices: int,
                   attempts: int = 5000) -> list[int] | None:
    """Tenta Warnsdorff (com desempate aleatório); se falhar, backtracking."""
    rng = random.Random(12345)
    p = warnsdorff_path(G, s, t, n_vertices, rng, tie_random=False)
    if p:
        return p
    for _ in range(attempts):
        p = warnsdorff_path(G, s, t, n_vertices, rng, tie_random=True)
        if p:
            return p
    # fallback determinístico
    return backtrack_one_path(G, s, t, n_vertices)


def backtrack_one_path(G: nx.Graph, s: int, t: int, n_vertices: int):
    """Acha UM caminho hamiltoniano s->t por backtracking (Warnsdorff-ordered)."""
    adj = {v: set(G.neighbors(v)) for v in G.nodes}
    visited = [False] * n_vertices
    path = [s]
    visited[s] = True

    def rec(cur):
        if len(path) == n_vertices:
            return cur == t
        nbrs = [w for w in adj[cur] if not visited[w]]
        nbrs.sort(key=lambda w: sum(1 for x in adj[w] if not visited[x]))
        for w in nbrs:
            visited[w] = True
            path.append(w)
            if rec(w):
                return True
            path.pop()
            visited[w] = False
        return False

    return list(path) if rec(s) else None


# ----------------------------------------------------------------------------
# Utilidades de arestas
# ----------------------------------------------------------------------------
def path_edges(path: list[int]) -> set[Edge]:
    return {frozenset((path[i], path[i + 1])) for i in range(len(path) - 1)}


def cycle_edges(cycle_nodes: list[int]) -> set[Edge]:
    es = set()
    k = len(cycle_nodes)
    for i in range(k):
        es.add(frozenset((cycle_nodes[i], cycle_nodes[(i + 1) % k])))
    return es


def is_valid_hamiltonian_path(edges: set[Edge], s: int, t: int,
                              n_vertices: int) -> bool:
    """edges é caminho hamiltoniano s->t? (graus, span, conexidade)."""
    if len(edges) != n_vertices - 1:
        return False
    deg: dict[int, int] = {}
    for e in edges:
        for v in e:
            deg[v] = deg.get(v, 0) + 1
    if len(deg) != n_vertices:
        return False  # algum vértice ficou isolado (grau 0)
    for v, d in deg.items():
        target = 1 if v in (s, t) else 2
        if d != target:
            return False
    # graus corretos + |E| = V-1 + s,t grau 1  =>  precisa ser conexo (1 componente)
    # verificação explícita de conexidade:
    adj: dict[int, list[int]] = {}
    for e in edges:
        u, w = tuple(e)
        adj.setdefault(u, []).append(w)
        adj.setdefault(w, []).append(u)
    seen = {s}
    stack = [s]
    while stack:
        x = stack.pop()
        for y in adj[x]:
            if y not in seen:
                seen.add(y)
                stack.append(y)
    return len(seen) == n_vertices


def classify_subgraph(edges: set[Edge], s: int, t: int,
                      n_vertices: int) -> str:
    """Classifica o resultado de um XOR: 'valid' | 'bad_count' |
    'bad_degree' | 'disconnected'. (Para diagnosticar o modo de falha.)"""
    if len(edges) != n_vertices - 1:
        return "bad_count"
    deg: dict[int, int] = {}
    for e in edges:
        for v in e:
            deg[v] = deg.get(v, 0) + 1
    if len(deg) != n_vertices:
        return "bad_degree"  # vértice isolado
    for v, d in deg.items():
        if d != (1 if v in (s, t) else 2):
            return "bad_degree"
    # graus + contagem ok => união de 1 caminho + ciclos; checar conexidade
    adj: dict[int, list[int]] = {}
    for e in edges:
        u, w = tuple(e)
        adj.setdefault(u, []).append(w)
        adj.setdefault(w, []).append(u)
    seen = {s}
    stack = [s]
    while stack:
        x = stack.pop()
        for y in adj[x]:
            if y not in seen:
                seen.add(y)
                stack.append(y)
    return "valid" if len(seen) == n_vertices else "disconnected"


# ----------------------------------------------------------------------------
# STEP 3 — Ciclos fundamentais via DFS
# ----------------------------------------------------------------------------
def fundamental_cycles(G: nx.Graph, root: int) -> list[set[Edge]]:
    """Base de ciclos fundamentais (uma árvore DFS-equivalente).
    networkx.cycle_basis devolve, para grafo conexo, |E|-|V|+1 ciclos."""
    cycles_nodes = nx.cycle_basis(G, root)
    return [cycle_edges(c) for c in cycles_nodes]


# ----------------------------------------------------------------------------
# STEP 4 — Geração de caminhos via Loop-XOR
# ----------------------------------------------------------------------------
def xor_enumerate(base_edges: set[Edge], cycles: list[set[Edge]],
                  s: int, t: int, n_vertices: int,
                  max_exhaustive_k: int = 20, n_samples: int = 10_000,
                  seed: int = 7):
    """Gera candidatos P XOR (uniao de subconjunto de ciclos), valida.

    Devolve (paths_validos, n_subsets_testados, compat_por_ciclo, modo)."""
    rng = random.Random(seed)
    k = len(cycles)
    valid: list[frozenset[Edge]] = []
    seen_keys: set[frozenset[Edge]] = set()
    compat_single = 0  # ciclos que sozinhos (XOR só com ele) geram caminho válido
    fail_modes = {"valid": 0, "bad_count": 0, "bad_degree": 0,
                  "disconnected": 0}

    # compatibilidade individual de cada ciclo
    for C in cycles:
        cand = base_edges ^ C
        if is_valid_hamiltonian_path(cand, s, t, n_vertices):
            compat_single += 1

    def consider(subset_indices):
        cand = set(base_edges)
        for i in subset_indices:
            cand ^= cycles[i]
        fail_modes[classify_subgraph(cand, s, t, n_vertices)] += 1
        if is_valid_hamiltonian_path(cand, s, t, n_vertices):
            key = frozenset(cand)
            if key not in seen_keys:
                seen_keys.add(key)
                valid.append(key)

    if k <= max_exhaustive_k:
        mode = f"exhaustive(2^{k})"
        n_tested = 0
        for r in range(0, k + 1):
            for combo in itertools.combinations(range(k), r):
                consider(combo)
                n_tested += 1
    else:
        mode = f"sampled({n_samples} of 2^{k})"
        n_tested = n_samples
        # inclui o caminho-base (subconjunto vazio) e cada ciclo isolado
        consider(())
        for i in range(k):
            consider((i,))
        for _ in range(n_samples):
            # subconjunto aleatório (tamanho variado, viés p/ pequenos)
            size = rng.randint(1, min(6, k))
            combo = rng.sample(range(k), size)
            consider(combo)

    return valid, n_tested, compat_single, mode, fail_modes


# ----------------------------------------------------------------------------
# STEP 5 — Baselines
# ----------------------------------------------------------------------------
def backtracking_all_paths(G: nx.Graph, s: int, t: int, n_vertices: int,
                           time_limit: float | None = None,
                           max_samples: int = 3000):
    """Conta/coleta caminhos hamiltonianos s->t por DFS com poda de
    conectividade (vértices não-visitados + atual devem permanecer conexos).
    Usa bitmasks (inteiros) p/ velocidade. Exaustivo se time_limit None.
    Devolve (count, sample_paths_edges, hit_limit, elapsed)."""
    import sys
    sys.setrecursionlimit(max(10000, n_vertices * 50))
    adjm = [0] * n_vertices
    for u, v in G.edges():
        adjm[u] |= 1 << v
        adjm[v] |= 1 << u
    FULL = (1 << n_vertices) - 1
    count = 0
    samples: list[frozenset[Edge]] = []
    path = [s]
    t0 = time.perf_counter()
    hit_limit = [False]

    def rec(cur, visited, depth):
        nonlocal count
        if depth == n_vertices:
            if cur == t:
                count += 1
                if len(samples) < max_samples:
                    samples.append(frozenset(path_edges(path)))
            return
        if time_limit is not None and (time.perf_counter() - t0) > time_limit:
            hit_limit[0] = True
            return
        free_un = FULL & ~visited
        # poda de grau: todo vértice não-visitado que será INTERIOR do trecho
        # restante (cur->...->t) precisa de >=2 vizinhos disponíveis (livres
        # ou o próprio cur); t precisa de >=1. Caso contrário, beco sem saída.
        avail_base = free_un | (1 << cur)
        f = free_un
        while f:
            b = f & -f
            v = b.bit_length() - 1
            f ^= b
            a = bin(adjm[v] & avail_base).count("1")
            if (a < 1) if v == t else (a < 2):
                return
        # poda de conectividade: flood-fill sobre não-visitados a partir de cur
        reach = adjm[cur] & free_un
        frontier = reach
        while frontier:
            nf = 0
            f = frontier
            while f:
                b = f & -f
                nf |= adjm[b.bit_length() - 1]
                f ^= b
            nf &= free_un & ~reach
            reach |= nf
            frontier = nf
        if free_un & ~reach:
            return  # algum não-visitado ficou inalcançável
        # ordem Warnsdorff (menos vizinhos livres primeiro)
        nb = []
        m = adjm[cur] & free_un
        while m:
            b = m & -m
            w = b.bit_length() - 1
            m ^= b
            nb.append((bin(adjm[w] & free_un).count("1"), w))
        nb.sort()
        for _, w in nb:
            path.append(w)
            rec(w, visited | (1 << w), depth + 1)
            path.pop()
            if hit_limit[0]:
                return

    rec(s, 1 << s, 1)
    return count, samples, hit_limit[0], time.perf_counter() - t0


def warnsdorff_repeated(G: nx.Graph, s: int, t: int, n_vertices: int,
                        n_runs: int = 100, seed: int = 99):
    """Warnsdorff com desempate aleatório, n_runs partidas a partir de s,
    exigindo terminar exatamente em t. Devolve (caminhos distintos, tempo,
    n_runs, n_sucessos)."""
    rng = random.Random(seed)
    found: dict[frozenset[Edge], None] = {}
    successes = 0
    t0 = time.perf_counter()
    for _ in range(n_runs):
        p = warnsdorff_path(G, s, t, n_vertices, rng, tie_random=True)
        if p:
            successes += 1
            found[frozenset(path_edges(p))] = None
    return list(found.keys()), time.perf_counter() - t0, n_runs, successes


# ----------------------------------------------------------------------------
# STEP 6 — Diversidade
# ----------------------------------------------------------------------------
def diversity_score(paths: list[frozenset[Edge]], max_pairs: int = 5000,
                    seed: int = 1) -> float:
    """Média da diferença de arestas (|simétrica diff|) entre pares de caminhos."""
    m = len(paths)
    if m < 2:
        return 0.0
    rng = random.Random(seed)
    total_pairs = m * (m - 1) // 2
    if total_pairs <= max_pairs:
        pairs = itertools.combinations(range(m), 2)
    else:
        pairs = ((rng.randrange(m), rng.randrange(m)) for _ in range(max_pairs))
    acc = 0
    cnt = 0
    for i, j in pairs:
        if i == j:
            continue
        acc += len(paths[i] ^ paths[j])
        cnt += 1
    return acc / cnt if cnt else 0.0


# ----------------------------------------------------------------------------
# Runner por tabuleiro
# ----------------------------------------------------------------------------
def run_board(n: int, exhaustive: bool, bt_time_limit: float | None):
    print(f"\n{'='*60}\nTABULEIRO {n}x{n}\n{'='*60}")
    G = build_knight_graph(n)
    V = n * n
    E = G.number_of_edges()
    beta1 = E - V + 1  # nº de ciclos fundamentais (conexo)
    print(f"V={V}  E={E}  ciclos fundamentais (beta1)={beta1}")

    # --- escolha de s, t com paridade viável ---
    s = 0
    t_pref = V - 1
    if feasible_parity(s, t_pref, n):
        t = t_pref
        st_note = f"(0,0)->({n-1},{n-1})"
    else:
        t = 1  # (0,1): cor oposta a (0,0)
        st_note = (f"(0,0)->({n-1},{n-1}) INVIÁVEL por paridade "
                   f"(mesma cor, {V-1} arestas ímpar) -> usando (0,0)->(0,1)")
    print(f"s={s} t={t}  {st_note}")

    # STEP 2
    base = find_base_path(G, s, t, V)
    base_found = base is not None
    print(f"Caminho-base: {'SIM' if base_found else 'NÃO'}")
    if not base_found:
        return {"n": n, "base_path_found": False, "st_note": st_note}
    base_edges = path_edges(base)
    assert is_valid_hamiltonian_path(base_edges, s, t, V), "base inválido?!"

    # STEP 3
    cycles = fundamental_cycles(G, s)
    lens = sorted(len(c) for c in cycles)
    from collections import Counter
    len_dist = dict(Counter(lens))
    print(f"Ciclos fundamentais coletados: {len(cycles)}  "
          f"(comprimentos min={lens[0]} max={lens[-1]})")
    print(f"Distribuição de comprimento: {len_dist}")

    # STEP 4 — XOR
    t0 = time.perf_counter()
    xor_paths, n_subsets, compat_single, xor_mode, fail_modes = xor_enumerate(
        base_edges, cycles, s, t, V)
    xor_time = time.perf_counter() - t0
    # o caminho-base sempre conta como "encontrado"
    base_key = frozenset(base_edges)
    xor_set = set(xor_paths)
    xor_set.add(base_key)
    xor_paths = list(xor_set)
    print(f"XOR modo={xor_mode}  subsets testados={n_subsets}")
    print(f"XOR caminhos válidos distintos: {len(xor_paths)}  "
          f"(tempo {xor_time*1000:.1f} ms)")
    print(f"Ciclos individualmente compatíveis: {compat_single}/{len(cycles)}"
          f" ({100*compat_single/len(cycles):.2f}%)")
    tot_fm = sum(fail_modes.values())
    print(f"Modos do XOR ({tot_fm} subsets): " + ", ".join(
        f"{k}={v} ({100*v/tot_fm:.1f}%)" for k, v in fail_modes.items()))

    # STEP 5A — backtracking
    bt_count, bt_samples, bt_hit, bt_time = backtracking_all_paths(
        G, s, t, V, time_limit=bt_time_limit)
    print(f"Backtracking: count={bt_count}"
          f"{' (LIMITE atingido, parcial)' if bt_hit else ' (exaustivo)'}"
          f"  tempo {bt_time*1000:.1f} ms")

    # STEP 5B — Warnsdorff repetido (exige terminar em t)
    ws_runs = 2000
    ws_paths, ws_time, ws_n, ws_succ = warnsdorff_repeated(
        G, s, t, V, n_runs=ws_runs)
    print(f"Warnsdorff {ws_n} partidas: {ws_succ} sucessos "
          f"({100*ws_succ/ws_n:.2f}%), {len(ws_paths)} caminhos distintos"
          f"  tempo {ws_time*1000:.1f} ms")

    # STEP 6 — métricas
    div_xor = diversity_score(xor_paths)
    div_ws = diversity_score(ws_paths)
    # diversidade de referência: amostra da população VERDADEIRA (backtracking)
    div_truth = diversity_score(bt_samples)
    print(f"Diversidade  XOR={div_xor:.2f}  Warnsdorff={div_ws:.2f}  "
          f"população-real(bt amostra n={len(bt_samples)})={div_truth:.2f}")

    result = {
        "n": n,
        "V": V, "E": E, "beta1": beta1,
        "s": s, "t": t, "st_note": st_note,
        "base_path_found": True,
        "fundamental_cycles": len(cycles),
        "cycle_len_min": lens[0], "cycle_len_max": lens[-1],
        "cycle_len_dist": len_dist,
        "xor_mode": xor_mode,
        "xor_subsets_tested": n_subsets,
        "xor_valid_paths": len(xor_paths),
        "xor_time_ms": round(xor_time * 1000, 3),
        "cycles_compatible_single": compat_single,
        "cycles_compatible_frac": round(compat_single / len(cycles), 6),
        "xor_failure_modes": fail_modes,
        "backtracking_count": bt_count,
        "backtracking_exhaustive": not bt_hit,
        "backtracking_time_ms": round(bt_time * 1000, 3),
        "warnsdorff_runs": ws_n,
        "warnsdorff_successes": ws_succ,
        "warnsdorff_success_rate": round(ws_succ / ws_n, 6),
        "warnsdorff_distinct_paths": len(ws_paths),
        "warnsdorff_time_ms": round(ws_time * 1000, 3),
        "diversity_xor": round(div_xor, 4),
        "diversity_warnsdorff": round(div_ws, 4),
        "diversity_truth_sample": round(div_truth, 4),
        "truth_sample_size": len(bt_samples),
    }

    if exhaustive and not bt_hit and bt_count > 0:
        cov = len(xor_paths) / bt_count
        result["coverage_xor"] = round(cov, 6)
        result["total_paths_st"] = bt_count
        print(f"COBERTURA XOR: {len(xor_paths)} / {bt_count} "
              f"= {100*cov:.4f}%")
    else:
        # 8x8: paths/sec
        result["xor_paths_per_sec"] = round(
            len(xor_paths) / xor_time if xor_time > 0 else 0, 2)
        result["warnsdorff_paths_per_sec"] = round(
            len(ws_paths) / ws_time if ws_time > 0 else 0, 2)
        print(f"Paths/s  XOR={result['xor_paths_per_sec']}  "
              f"Warnsdorff={result['warnsdorff_paths_per_sec']}")

    return result


def main():
    print("PATH FINDING s->t via DFS LOOP DECOMPOSITION (XOR)")
    results = {}
    results["6x6"] = run_board(6, exhaustive=True, bt_time_limit=None)
    results["8x8"] = run_board(8, exhaustive=False, bt_time_limit=60.0)

    out = RESULTS_DIR / "path_decomposition_experiment.json"
    out.write_text(json.dumps(results, indent=2, ensure_ascii=False))
    print(f"\nResultados salvos em {out}")

    # ----- relatório formatado -----
    r6 = results["6x6"]
    print("\n=== PATH FINDING s->t : 6x6 KNIGHT ===")
    print(f"Base path found: {'YES' if r6['base_path_found'] else 'NO'}")
    print(f"Fundamental cycles: {r6['fundamental_cycles']}")
    if "coverage_xor" in r6:
        print(f"XOR method coverage: {r6['xor_valid_paths']} / "
              f"{r6['total_paths_st']} ({100*r6['coverage_xor']:.4f}%)")
    print(f"XOR time: {r6['xor_time_ms']}ms | "
          f"Backtracking time: {r6['backtracking_time_ms']}ms")
    print(f"Diversity score XOR: {r6['diversity_xor']} | "
          f"Diversity score Warnsdorff: {r6['diversity_warnsdorff']} | "
          f"true-population: {r6['diversity_truth_sample']}")
    print(f"Warnsdorff success rate (s->t fixo): "
          f"{100*r6['warnsdorff_success_rate']:.2f}% "
          f"({r6['warnsdorff_successes']}/{r6['warnsdorff_runs']})")

    r8 = results["8x8"]
    print("\n=== PATH FINDING s->t : 8x8 KNIGHT ===")
    print(f"Fundamental cycles: {r8['fundamental_cycles']} "
          f"({r8['xor_mode']})")
    print(f"Valid paths found: {r8['xor_valid_paths']} in "
          f"{r8['xor_time_ms']/1000:.3f}s")
    print(f"Paths/sec XOR: {r8.get('xor_paths_per_sec','-')} | "
          f"Paths/sec Warnsdorff: {r8.get('warnsdorff_paths_per_sec','-')}")
    print(f"Diversity score: XOR={r8['diversity_xor']} "
          f"Warnsdorff={r8['diversity_warnsdorff']} "
          f"true-population={r8['diversity_truth_sample']}")
    print(f"Warnsdorff success rate (s->t fixo): "
          f"{100*r8['warnsdorff_success_rate']:.2f}% "
          f"({r8['warnsdorff_successes']}/{r8['warnsdorff_runs']})")


if __name__ == "__main__":
    main()
