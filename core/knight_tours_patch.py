"""knight_tours_patch.py — knight_tours.py + LUT de deformação GF(2).

Igual ao knight_tours.py, mas com:
  1. patch_LUT pré-computada em build_graph() — todos os 4-ciclos
     (diamantes) do grafo do cavalo.
  2. _apply_patch() helper — aplica Estado ⊕ C_patch revertendo e1
     e destruindo e2 (e3, e4 são deixadas para a propagação normal).
  3. _process_queue() consulta a LUT quando um sub-ciclo é detectado
     (ru == rv) e tenta deformar antes de retornar SUBTOUR.

Estatísticas de patch são contadas em ctx['patch_count'] / ctx['patch_miss']
e retornadas via knight_tours_patch(..., return_stats=True).
"""

# --- bootstrap: torna core/ importavel a partir de qualquer subpasta ---
import sys as _sys, pathlib as _pathlib
for _anc in _pathlib.Path(__file__).resolve().parents:
    if (_anc / "core").is_dir() and (_anc / "pyproject.toml").exists():
        _sys.path.insert(0, str(_anc / "core"))
        break
# --- fim do bootstrap ---

import sys
import numpy as np

sys.setrecursionlimit(10000)


# ---------------------------------------------------------------------------
# Constantes — heurística de fase local f∞(L)
# ---------------------------------------------------------------------------

F_INF = {0: 0.528, 1: 0.193, 2: 0.198, 3: 0.294, 4: 0.247, 5: 0.261}
F_INF_DEFAULT = 0.25

FREE = 0
ACTIVE = 1
INACTIVE = 2

OK = 0
CONTRADICTION = 1
SUBTOUR = 2
COMPLETE_TOUR = 3

_KNIGHT_MOVES = ((-2, -1), (-2, 1), (-1, -2), (-1, 2),
                 (1, -2), (1, 2), (2, -1), (2, 1))


# ---------------------------------------------------------------------------
# Componente 1 — Grafo do cavalo + LUT de diamantes (4-ciclos)
# ---------------------------------------------------------------------------

def build_patch_LUT(edge_endpoints, adj_edges, E):
    """Pré-computa, para cada aresta e1, todas as triplas (e2, e3, e4) tais
    que e1-e3-e2-e4 forma um 4-ciclo no grafo do cavalo.

    Geometria do diamante:
            u --e1-- v
            |        |
           e3       e4
            |        |
            x --e2-- y

    e1, e2 = paredes a DESTRUIR
    e3, e4 = pontes a CRIAR

    Os 4 vértices {u,v,x,y} têm grau LÍQUIDO invariante após o patch:
    cada um perde 1 (parede) e ganha 1 (ponte).
    """
    edge_to_idx = {}
    for ei in range(E):
        u = int(edge_endpoints[ei, 0])
        v = int(edge_endpoints[ei, 1])
        edge_to_idx[(u, v)] = ei
        edge_to_idx[(v, u)] = ei

    LUT = [[] for _ in range(E)]

    for e1 in range(E):
        u = int(edge_endpoints[e1, 0])
        v = int(edge_endpoints[e1, 1])

        for e3 in adj_edges[u]:
            e3 = int(e3)
            if e3 == e1:
                continue
            ep = edge_endpoints[e3]
            x = int(ep[1]) if int(ep[0]) == u else int(ep[0])
            if x == v:
                continue

            for e4 in adj_edges[v]:
                e4 = int(e4)
                if e4 == e1 or e4 == e3:
                    continue
                ep4 = edge_endpoints[e4]
                y = int(ep4[1]) if int(ep4[0]) == v else int(ep4[0])
                if y == u or y == x:
                    continue

                e2 = edge_to_idx.get((x, y), -1)
                if e2 == -1 or e2 == e1:
                    continue

                LUT[e1].append((e2, e3, e4))

    return LUT


def build_graph(n):
    """Igual ao knight_tours.py, mas inclui patch_LUT e n2 no ctx."""
    V = n * n
    edges = []
    for r in range(n):
        for c in range(n):
            u = r * n + c
            for dr, dc in _KNIGHT_MOVES:
                rr, cc = r + dr, c + dc
                if 0 <= rr < n and 0 <= cc < n:
                    w = rr * n + cc
                    if u < w:
                        edges.append((u, w))

    E = len(edges)
    edge_endpoints = np.asarray(edges, dtype=np.int32)

    adj_lists = [[] for _ in range(V)]
    for ei, (u, w) in enumerate(edges):
        adj_lists[u].append(ei)
        adj_lists[w].append(ei)
    adj_edges = [np.asarray(lst, dtype=np.int32) for lst in adj_lists]
    total_incident = np.asarray([len(lst) for lst in adj_lists], dtype=np.int32)

    patch_LUT = build_patch_LUT(edge_endpoints, adj_edges, E)

    return {
        'n': n, 'V': V, 'E': E, 'n2': V,
        'edge_endpoints': edge_endpoints,
        'adj_edges': adj_edges,
        'total_incident': total_incident,
        'patch_LUT': patch_LUT,
        'patch_count': 0,
        'patch_miss': 0,
    }


def vertex_level(v, n):
    r, c = divmod(int(v), n)
    return min(r, c, n - 1 - r, n - 1 - c)


def edge_level(u, v, n):
    return min(vertex_level(u, n), vertex_level(v, n))


# ---------------------------------------------------------------------------
# Heurística de fase local
# ---------------------------------------------------------------------------

def edge_priority(L):
    f = F_INF.get(L, F_INF_DEFAULT)
    return abs(f - 0.5), (1 if f > 0.5 else 0)


# ---------------------------------------------------------------------------
# Estado da busca
# ---------------------------------------------------------------------------

class State:
    __slots__ = ('fixed', 'degree', 'n_inactive',
                 'uf_parent', 'uf_size', 'uf_deg2', 'n_free')

    def __init__(self, V, E):
        self.fixed = np.zeros(E, dtype=np.uint8)
        self.degree = np.zeros(V, dtype=np.int32)
        self.n_inactive = np.zeros(V, dtype=np.int32)
        self.uf_parent = np.arange(V, dtype=np.int32)
        self.uf_size = np.ones(V, dtype=np.int32)
        self.uf_deg2 = np.zeros(V, dtype=np.int32)
        self.n_free = E

    def snapshot(self):
        return (self.fixed.copy(), self.degree.copy(), self.n_inactive.copy(),
                self.uf_parent.copy(), self.uf_size.copy(),
                self.uf_deg2.copy(), self.n_free)

    def restore(self, snap):
        np.copyto(self.fixed, snap[0])
        np.copyto(self.degree, snap[1])
        np.copyto(self.n_inactive, snap[2])
        np.copyto(self.uf_parent, snap[3])
        np.copyto(self.uf_size, snap[4])
        np.copyto(self.uf_deg2, snap[5])
        self.n_free = snap[6]


def uf_find(state, v):
    while state.uf_parent[v] != v:
        v = int(state.uf_parent[v])
    return v


# ---------------------------------------------------------------------------
# Aplicação do patch (XOR em GF(2))
# ---------------------------------------------------------------------------

def _apply_patch(state, e1, e2, edge_endpoints, new_deg2, ru):
    """Reverte e1 (acabou de virar ACTIVE em _process_queue) e destrói e2
    (era ACTIVE em outro componente). Ambas voltam a FREE.

    e3 e e4 NÃO são tocadas aqui — serão enfileiradas como (e, 1) e o
    fluxo normal de _process_queue aplica + faz o union do UF.

    O UF (uf_parent/uf_size) NÃO é split — o snapshot/restore do backtrack
    é quem cuida do rollback. Aqui só ajustamos uf_deg2.
    """
    # e1: ACTIVE → FREE (revert do incremento feito imediatamente antes)
    state.fixed[e1] = FREE
    u1 = int(edge_endpoints[e1, 0])
    v1 = int(edge_endpoints[e1, 1])
    state.degree[u1] -= 1
    state.degree[v1] -= 1
    state.n_free += 1
    state.uf_deg2[ru] -= new_deg2  # desfaz incremento recém-aplicado

    # e2: ACTIVE → FREE (destroy)
    state.fixed[e2] = FREE
    u2 = int(edge_endpoints[e2, 0])
    v2 = int(edge_endpoints[e2, 1])
    # se algum endpoint perde grau-2 (era 2, agora 1), uf_deg2 decrementa
    lost = (1 if int(state.degree[u2]) == 2 else 0) + \
           (1 if int(state.degree[v2]) == 2 else 0)
    state.degree[u2] -= 1
    state.degree[v2] -= 1
    state.n_free += 1
    if lost > 0:
        root2 = uf_find(state, u2)
        state.uf_deg2[root2] -= lost


# ---------------------------------------------------------------------------
# Componente 3a — Propagação R2 + detector de sub-ciclos + LUT de deformação
# ---------------------------------------------------------------------------

def _process_queue(state, ctx, queue):
    n2 = ctx['n2']
    total_inc = ctx['total_incident']
    adj_edges = ctx['adj_edges']
    edge_endpoints = ctx['edge_endpoints']
    patch_LUT = ctx['patch_LUT']

    while queue:
        e, val = queue.pop()
        cur = state.fixed[e]
        if cur != FREE:
            if (cur == ACTIVE and val == 1) or (cur == INACTIVE and val == 0):
                continue
            return CONTRADICTION

        u = int(edge_endpoints[e, 0])
        v = int(edge_endpoints[e, 1])

        if val == 1:
            state.fixed[e] = ACTIVE
            state.n_free -= 1
            state.degree[u] += 1
            state.degree[v] += 1
            if state.degree[u] > 2 or state.degree[v] > 2:
                return CONTRADICTION

            ru = uf_find(state, u)
            rv = uf_find(state, v)
            new_deg2 = (1 if state.degree[u] == 2 else 0) + \
                       (1 if state.degree[v] == 2 else 0)

            if ru == rv:
                state.uf_deg2[ru] += new_deg2
                if state.uf_size[ru] == n2:
                    return COMPLETE_TOUR

                # LEMA DE DEFORMAÇÃO GF(2) — consulta O(1) na LUT
                patch_applied = False
                for (e2, e3, e4) in patch_LUT[e]:
                    if state.fixed[e2] != ACTIVE:
                        continue
                    if state.fixed[e3] != FREE:
                        continue
                    if state.fixed[e4] != FREE:
                        continue

                    # e2 deve estar em componente diferente de e1
                    u2 = int(edge_endpoints[e2, 0])
                    if uf_find(state, u2) == ru:
                        continue

                    # PATCH VÁLIDO — aplicar XOR em GF(2)
                    _apply_patch(state, e, e2, edge_endpoints, new_deg2, ru)

                    # Enfileira as pontes para serem processadas normalmente
                    # (a propagação faz o union no UF via fluxo padrão)
                    queue.append((e3, 1))
                    queue.append((e4, 1))

                    ctx['patch_count'] += 1
                    patch_applied = True
                    break

                if not patch_applied:
                    ctx['patch_miss'] += 1
                    return SUBTOUR

                # Patch aplicado — pula o R2 deste e (que foi revertido)
                continue

            # Union por tamanho
            if state.uf_size[ru] < state.uf_size[rv]:
                ru, rv = rv, ru
            state.uf_parent[rv] = ru
            state.uf_size[ru] += state.uf_size[rv]
            state.uf_deg2[ru] += state.uf_deg2[rv] + new_deg2
        else:
            state.fixed[e] = INACTIVE
            state.n_free -= 1
            state.n_inactive[u] += 1
            state.n_inactive[v] += 1
            if total_inc[u] - state.n_inactive[u] < 2:
                return CONTRADICTION
            if total_inc[v] - state.n_inactive[v] < 2:
                return CONTRADICTION

        # R2 nos dois endpoints
        for w in (u, v):
            dw = int(state.degree[w])
            avail = int(total_inc[w] - state.n_inactive[w])
            if dw == 2:
                for ei in adj_edges[w]:
                    if state.fixed[ei] == FREE:
                        queue.append((int(ei), 0))
            elif avail == 2 and dw < 2:
                for ei in adj_edges[w]:
                    if state.fixed[ei] == FREE:
                        queue.append((int(ei), 1))

    return OK


def fix_and_propagate(state, ctx, e, val):
    return _process_queue(state, ctx, [(e, val)])


def _propagate_initial(state, ctx):
    V = ctx['V']
    total_inc = ctx['total_incident']
    adj_edges = ctx['adj_edges']
    queue = []
    for v in range(V):
        if total_inc[v] - state.n_inactive[v] == 2 and state.degree[v] < 2:
            for ei in adj_edges[v]:
                if state.fixed[ei] == FREE:
                    queue.append((int(ei), 1))
    return _process_queue(state, ctx, queue)


# ---------------------------------------------------------------------------
# Componente 3b — Seleção de variável
# ---------------------------------------------------------------------------

def _choose_next_edge(state, ctx, rng):
    V = ctx['V']
    n = ctx['n']
    total_inc = ctx['total_incident']
    adj_edges = ctx['adj_edges']
    edge_endpoints = ctx['edge_endpoints']

    best_score = -1
    candidates = []
    for v in range(V):
        free_v = int(total_inc[v]) - int(state.degree[v]) - int(state.n_inactive[v])
        if free_v == 0:
            continue
        score = int(state.degree[v]) * 100 + \
                int(total_inc[v] - state.n_inactive[v] - 2)
        if score > best_score:
            best_score = score
            candidates = [v]
        elif score == best_score:
            candidates.append(v)

    if not candidates:
        return -1, 0

    if len(candidates) == 1:
        v_star = candidates[0]
    else:
        v_star = int(candidates[int(rng.integers(len(candidates)))])

    best_pri = -1.0
    best_edge = -1
    best_first = 0
    for ei in adj_edges[v_star]:
        ei_i = int(ei)
        if state.fixed[ei_i] != FREE:
            continue
        a = int(edge_endpoints[ei_i, 0])
        b = int(edge_endpoints[ei_i, 1])
        L = edge_level(a, b, n)
        pri, fv = edge_priority(L)
        if pri > best_pri:
            best_pri = pri
            best_edge = ei_i
            best_first = fv

    return best_edge, best_first


# ---------------------------------------------------------------------------
# Componente 3c — Verificação / extração
# ---------------------------------------------------------------------------

def _is_complete_tour(state, ctx):
    V = ctx['V']
    if not np.all(state.degree == 2):
        return False
    root0 = uf_find(state, 0)
    return int(state.uf_size[root0]) == V


def _extract_tour(state, ctx):
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

    for v in range(V):
        if len(nbrs[v]) != 2:
            raise RuntimeError(
                f"vértice {v} tem {len(nbrs[v])} arestas ativas (esperava 2)")

    tour = np.zeros(V, dtype=np.int32)
    tour[0] = 0
    nxt = min(nbrs[0])
    tour[1] = nxt
    prev = 0
    cur = nxt
    for i in range(2, V):
        a, b = nbrs[cur]
        nxt = a if a != prev else b
        tour[i] = nxt
        prev = cur
        cur = nxt
    return tour


# ---------------------------------------------------------------------------
# Componente 3d — Backtracking
# ---------------------------------------------------------------------------

def _backtrack(state, ctx, rng, tours, K):
    if len(tours) >= K:
        return True

    e, first_val = _choose_next_edge(state, ctx, rng)

    if e == -1:
        if _is_complete_tour(state, ctx):
            tours.append(_extract_tour(state, ctx))
        return len(tours) >= K

    for val in (first_val, 1 - first_val):
        snap = state.snapshot()
        status = fix_and_propagate(state, ctx, e, val)
        if status == OK:
            if _backtrack(state, ctx, rng, tours, K):
                state.restore(snap)
                return True
        elif status == COMPLETE_TOUR:
            if _is_complete_tour(state, ctx):
                tours.append(_extract_tour(state, ctx))
            if len(tours) >= K:
                state.restore(snap)
                return True
        state.restore(snap)

    return len(tours) >= K


# ---------------------------------------------------------------------------
# Interface pública
# ---------------------------------------------------------------------------

def knight_tours_patch(n, K, seed=None, return_stats=False):
    """Igual a knight_tours, mas com tentativa de deformação GF(2) por LUT.

    Se return_stats=True, retorna (tours, stats) onde
    stats = {'patch_count', 'patch_miss', 'lut_size'}.
    """
    if n < 6:
        raise ValueError("n deve ser ≥ 6")
    if K < 0:
        raise ValueError("K deve ser não-negativo")
    if K == 0:
        return ([], {'patch_count': 0, 'patch_miss': 0, 'lut_size': 0}) if return_stats else []

    ctx = build_graph(n)
    rng = np.random.default_rng(seed)
    state = State(ctx['V'], ctx['E'])

    status = _propagate_initial(state, ctx)
    tours = []
    if status in (CONTRADICTION, SUBTOUR):
        pass
    elif status == COMPLETE_TOUR:
        if _is_complete_tour(state, ctx):
            tours.append(_extract_tour(state, ctx))
    else:
        _backtrack(state, ctx, rng, tours, K)

    if return_stats:
        stats = {
            'patch_count': ctx['patch_count'],
            'patch_miss': ctx['patch_miss'],
            'lut_size': sum(len(lst) for lst in ctx['patch_LUT']),
        }
        return tours, stats
    return tours


# alias para clareza no benchmark
knight_tours = knight_tours_patch


def verify_tour(tour, n):
    V = n * n
    if len(tour) != V:
        return False
    seen = set()
    for x in tour:
        xi = int(x)
        if xi < 0 or xi >= V or xi in seen:
            return False
        seen.add(xi)
    moves = set(_KNIGHT_MOVES)
    for i in range(V):
        u = int(tour[i])
        w = int(tour[(i + 1) % V])
        ru, cu = divmod(u, n)
        rw, cw = divmod(w, n)
        if (rw - ru, cw - cu) not in moves:
            return False
    return True


# ---------------------------------------------------------------------------
# Validação / benchmark
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    import time

    # ============================================================
    # TAREFA 1 — validação n=6 (GT = 9862)
    # ============================================================
    print("=" * 60)
    print("TAREFA 1 — Validação n=6 (GT = 9862)")
    print("=" * 60)

    import knight_tours as kt_orig

    # Original — enumeração completa
    t0 = time.time()
    tours_orig = kt_orig.knight_tours(6, 10**6, seed=0)
    t1 = time.time()
    print(f"\nOriginal  : {len(tours_orig):>6d} tours  {t1-t0:.3f}s")

    # Patch LUT — enumeração completa
    t0 = time.time()
    tours_patch, stats = knight_tours_patch(6, 10**6, seed=0, return_stats=True)
    t1 = time.time()
    print(f"Patch LUT : {len(tours_patch):>6d} tours  {t1-t0:.3f}s")
    print(f"  LUT size              : {stats['lut_size']}")
    print(f"  patches aplicados     : {stats['patch_count']}")
    print(f"  patches não encontrados: {stats['patch_miss']}")
    denom = stats['patch_count'] + stats['patch_miss']
    rescue = stats['patch_count'] / denom if denom else 0.0
    print(f"  taxa de resgate       : {rescue:.1%}")

    assert len(tours_patch) == 9862, \
        f"ERRO: {len(tours_patch)} != 9862 — abortando antes do benchmark"
    assert all(verify_tour(t, 6) for t in tours_patch), \
        "ERRO: tour inválido encontrado"
    print("\n✓ Contagem n=6 correta (9862) e todos tours válidos")

    # ============================================================
    # TAREFA 2 — Benchmark n=10, K=200
    # ============================================================
    print("\n" + "=" * 60)
    print("TAREFA 2 — Benchmark 10×10, K=200")
    print("=" * 60)

    K = 200

    t0 = time.time()
    tours_a = kt_orig.knight_tours(10, K, seed=0)
    ta = time.time() - t0
    print(f"\nA) Original  : {len(tours_a):>4d} tours  {ta:.3f}s "
          f"({ta/len(tours_a)*1000:.2f}ms/tour)")

    t0 = time.time()
    tours_b, stats_b = knight_tours_patch(10, K, seed=0, return_stats=True)
    tb = time.time() - t0
    print(f"B) Patch LUT : {len(tours_b):>4d} tours  {tb:.3f}s "
          f"({tb/len(tours_b)*1000:.2f}ms/tour)")
    print(f"   LUT size              : {stats_b['lut_size']}")
    print(f"   patches aplicados     : {stats_b['patch_count']}")
    print(f"   patches não encontrados: {stats_b['patch_miss']}")
    denom = stats_b['patch_count'] + stats_b['patch_miss']
    rescue = stats_b['patch_count'] / denom if denom else 0.0
    print(f"   taxa de resgate       : {rescue:.1%}")
    print(f"\n   Speedup patch vs original: {ta/tb:.2f}×")

    # ============================================================
    # TAREFA 3 — Análise da LUT por nível L
    # ============================================================
    print("\n" + "=" * 60)
    print("TAREFA 3 — Estrutura da LUT por nível L")
    print("=" * 60)

    for n in (6, 8, 10, 12):
        t0 = time.time()
        ctx = build_graph(n)
        t_build = time.time() - t0
        LUT = ctx['patch_LUT']
        E = ctx['E']
        edge_endpoints = ctx['edge_endpoints']
        sizes = [len(LUT[e]) for e in range(E)]

        print(f"\nn={n}: V={ctx['V']}, E={E}, t_build={t_build:.3f}s")
        print(f"  Arestas com ≥1 diamante  : {sum(s>0 for s in sizes):>4d}/{E}")
        print(f"  Arestas sem diamante     : {sum(s==0 for s in sizes):>4d}/{E}")
        print(f"  Média diamantes/aresta   : {np.mean(sizes):.2f}")
        print(f"  Máx diamantes/aresta     : {max(sizes)}")
        print(f"  Total diamantes na LUT   : {sum(sizes)}")

        # Por nível
        from collections import defaultdict
        by_L = defaultdict(list)
        for e in range(E):
            a = int(edge_endpoints[e, 0])
            b = int(edge_endpoints[e, 1])
            L = edge_level(a, b, n)
            by_L[L].append(sizes[e])
        print(f"  Por nível L:")
        for L in sorted(by_L.keys()):
            arr = by_L[L]
            print(f"    L={L}: {len(arr):>3d} arestas, "
                  f"média {np.mean(arr):>5.2f}, "
                  f"máx {max(arr):>3d}")
