"""knight_tours_klein.py — Gerador de tours fechados do cavalo na
"garrafa de Klein" (banda de Möbius com Y rígido).

Identificação: (x, y) ~ (x + m, n-1-y) — glide reflection no eixo x.
X é modular; Y é rígido com flip quando o movimento atravessa a fronteira X.

Algoritmo:
    Backtracking + propagação R2 + Union-Find incremental.
    Ancoragem em v=0; first_val sorteado aleatoriamente (corrige viés
    visto no toro).
    Heurística f∞(L) com L = min(y, n-1-y) (Y ainda é rígido).
"""

import sys
import numpy as np

sys.setrecursionlimit(10000)


F_INF = {0: 0.528, 1: 0.193, 2: 0.198, 3: 0.294, 4: 0.247, 5: 0.261}
F_INF_DEFAULT = 0.25

FREE, ACTIVE, INACTIVE = 0, 1, 2
OK, CONTRADICTION, SUBTOUR, COMPLETE_TOUR = 0, 1, 2, 3

_KNIGHT_MOVES = ((-2, -1), (-2, 1), (-1, -2), (-1, 2),
                 (1, -2), (1, 2), (2, -1), (2, 1))


# ---------------------------------------------------------------------------
# Aplicação do movimento de cavalo na Klein
# ---------------------------------------------------------------------------

def klein_step(y, x, dy, dx, n, m):
    """Aplica (dy, dx) a (y, x) na Klein. Retorna (new_y, new_x, wrap) ou None.
    `wrap` ∈ {0, 1}: 1 se atravessou a fronteira X (com flip Y)."""
    naive_x = x + dx
    naive_y = y + dy
    if naive_x < 0:
        wrap = 1
        new_x = naive_x + m
    elif naive_x >= m:
        wrap = 1
        new_x = naive_x - m
    else:
        wrap = 0
        new_x = naive_x
    new_y = (n - 1 - naive_y) if wrap else naive_y
    if new_y < 0 or new_y >= n:
        return None
    return new_y, new_x, wrap


# ---------------------------------------------------------------------------
# Componente 1 — Grafo
# ---------------------------------------------------------------------------

def build_graph(n, m):
    V = n * m
    edge_set = set()
    for y in range(n):
        for x in range(m):
            u = y * m + x
            for dy, dx in _KNIGHT_MOVES:
                step = klein_step(y, x, dy, dx, n, m)
                if step is None:
                    continue
                new_y, new_x, _wrap = step
                w = new_y * m + new_x
                if u == w:
                    continue
                a, b = (u, w) if u < w else (w, u)
                edge_set.add((a, b))
    edges = sorted(edge_set)
    E = len(edges)
    edge_endpoints = np.asarray(edges, dtype=np.int32) if E > 0 \
        else np.zeros((0, 2), dtype=np.int32)

    adj_lists = [[] for _ in range(V)]
    for ei, (u, w) in enumerate(edges):
        adj_lists[u].append(ei)
        adj_lists[w].append(ei)
    adj_edges = [np.asarray(lst, dtype=np.int32) for lst in adj_lists]
    total_incident = np.asarray([len(lst) for lst in adj_lists], dtype=np.int32)

    return {
        'n': n, 'm': m, 'V': V, 'E': E,
        'edge_endpoints': edge_endpoints,
        'adj_edges': adj_edges,
        'total_incident': total_incident,
    }


def vertex_level(v, n, m):
    """L(v) = min(y, n-1-y) — Y rígido."""
    y = int(v) // m
    return min(y, n - 1 - y)


def edge_level(u, v, n, m):
    return min(vertex_level(u, n, m), vertex_level(v, n, m))


def edge_priority(L):
    f = F_INF.get(L, F_INF_DEFAULT)
    return abs(f - 0.5), (1 if f > 0.5 else 0)


# ---------------------------------------------------------------------------
# Estado
# ---------------------------------------------------------------------------

class State:
    __slots__ = ('fixed', 'degree', 'n_inactive',
                 'uf_parent', 'uf_size', 'uf_deg2', 'n_free',
                 'nodes_explored')

    def __init__(self, V, E):
        self.fixed = np.zeros(E, dtype=np.uint8)
        self.degree = np.zeros(V, dtype=np.int32)
        self.n_inactive = np.zeros(V, dtype=np.int32)
        self.uf_parent = np.arange(V, dtype=np.int32)
        self.uf_size = np.ones(V, dtype=np.int32)
        self.uf_deg2 = np.zeros(V, dtype=np.int32)
        self.n_free = E
        self.nodes_explored = 0

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
# Propagação R2 + detector (idêntico ao cilindro/toro)
# ---------------------------------------------------------------------------

def _process_queue(state, ctx, queue):
    V = ctx['V']
    total_inc = ctx['total_incident']
    adj_edges = ctx['adj_edges']
    edge_endpoints = ctx['edge_endpoints']

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
                if state.uf_size[ru] == V:
                    return COMPLETE_TOUR
                return SUBTOUR

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
# Seleção (ancora em v=0, randomiza first_val)
# ---------------------------------------------------------------------------

def _choose_next_edge(state, ctx, rng, anchor_v0=False):
    V = ctx['V']
    n = ctx['n']
    m = ctx['m']
    total_inc = ctx['total_incident']
    adj_edges = ctx['adj_edges']
    edge_endpoints = ctx['edge_endpoints']

    if anchor_v0:
        free_0 = int(total_inc[0]) - int(state.degree[0]) - int(state.n_inactive[0])
        candidates = [0] if free_0 > 0 else []
    else:
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
        if anchor_v0:
            return _choose_next_edge(state, ctx, rng, anchor_v0=False)
        return -1, 0

    if len(candidates) == 1:
        v_star = candidates[0]
    else:
        v_star = int(candidates[int(rng.integers(len(candidates)))])

    best_pri = -1.0
    top_edges = []
    for ei in adj_edges[v_star]:
        ei_i = int(ei)
        if state.fixed[ei_i] != FREE:
            continue
        a = int(edge_endpoints[ei_i, 0])
        b = int(edge_endpoints[ei_i, 1])
        L = edge_level(a, b, n, m)
        pri, _fv = edge_priority(L)
        if pri > best_pri:
            best_pri = pri
            top_edges = [ei_i]
        elif pri == best_pri:
            top_edges.append(ei_i)

    if not top_edges:
        return -1, 0
    best_edge = top_edges[int(rng.integers(len(top_edges)))]
    # first_val SORTEADO (corrige viés de marginal visto no toro)
    first_val = int(rng.integers(2))
    return best_edge, first_val


# ---------------------------------------------------------------------------
# Verify / extract / backtrack
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


def _backtrack(state, ctx, rng, tours, K, depth=0):
    if len(tours) >= K:
        return True
    state.nodes_explored += 1
    e, first_val = _choose_next_edge(state, ctx, rng, anchor_v0=(depth == 0))
    if e == -1:
        if _is_complete_tour(state, ctx):
            tours.append(_extract_tour(state, ctx))
        return len(tours) >= K
    for val in (first_val, 1 - first_val):
        snap = state.snapshot()
        status = fix_and_propagate(state, ctx, e, val)
        if status == OK:
            if _backtrack(state, ctx, rng, tours, K, depth + 1):
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
# Interface
# ---------------------------------------------------------------------------

def knight_tours(n, m, K, seed=None, log_nodes=False):
    if n < 3 or m < 3:
        raise ValueError("n, m devem ser ≥ 3")
    if K < 0:
        raise ValueError("K não-negativo")
    if K == 0:
        return []

    ctx = build_graph(n, m)
    rng = np.random.default_rng(seed)
    state = State(ctx['V'], ctx['E'])

    status = _propagate_initial(state, ctx)
    if status in (CONTRADICTION, SUBTOUR):
        if log_nodes:
            print(f"[klein {n}x{m}] inviável (status={status})")
        return []

    tours = []
    if status == COMPLETE_TOUR:
        if _is_complete_tour(state, ctx):
            tours.append(_extract_tour(state, ctx))
    else:
        _backtrack(state, ctx, rng, tours, K)

    if log_nodes:
        per_tour = state.nodes_explored / max(len(tours), 1)
        print(f"[klein {n}x{m}] tours={len(tours)}/{K} "
              f"nós={state.nodes_explored} ({per_tour:.2f}/tour)")

    return tours[:K]


def verify_tour(tour, n, m):
    """Verifica tour fechado na Klein."""
    V = n * m
    if len(tour) != V:
        return False
    seen = set()
    for x in tour:
        xi = int(x)
        if xi < 0 or xi >= V or xi in seen:
            return False
        seen.add(xi)

    for i in range(V):
        u = int(tour[i])
        w = int(tour[(i + 1) % V])
        y_u, x_u = divmod(u, m)
        y_v, x_v = divmod(w, m)
        ok = False
        for dy, dx in _KNIGHT_MOVES:
            step = klein_step(y_u, x_u, dy, dx, n, m)
            if step is None:
                continue
            ny, nx, _wrap = step
            if (ny, nx) == (y_v, x_v):
                ok = True
                break
        if not ok:
            return False
    return True


# ---------------------------------------------------------------------------
# Testes
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    import time

    ctx = build_graph(6, 6)
    print(f"[grafo klein 6x6] V={ctx['V']} E={ctx['E']}")
    grau_min = int(ctx['total_incident'].min())
    grau_max = int(ctx['total_incident'].max())
    print(f"  grau ∈ [{grau_min}, {grau_max}]  (esperado: mín ≥ 4)")
    assert grau_min >= 4, f"grau mín {grau_min} < 4"

    print()
    print("=== n=m=6 klein, K=10 ===")
    t0 = time.time()
    tours = knight_tours(6, 6, 10, seed=0, log_nodes=True)
    t1 = time.time()
    print(f"  tempo: {t1-t0:.3f}s")
    assert len(tours) == 10
    for t in tours:
        assert verify_tour(t, 6, 6), "tour Klein inválido"
    print(f"  validados {len(tours)}")

    print()
    print("=== Reprodutibilidade ===")
    a = knight_tours(6, 6, 5, seed=42)
    b = knight_tours(6, 6, 5, seed=42)
    assert all(np.array_equal(x, y) for x, y in zip(a, b))
    print("  ok")

    print()
    print("=== Diversidade ===")
    tours = knight_tours(6, 6, 50, seed=0)
    unique = {tuple(int(x) for x in t) for t in tours}
    print(f"  únicos: {len(unique)}/{len(tours)}")

    print()
    print("=== Escalonamento ===")
    for (n, m) in [(6, 6), (6, 8), (8, 6), (8, 8), (6, 10), (10, 10)]:
        t0 = time.time()
        tours = knight_tours(n, m, 20, seed=0, log_nodes=True)
        t1 = time.time()
        ok = all(verify_tour(t, n, m) for t in tours)
        print(f"  {n}x{m}: {len(tours)} tours em {t1-t0:.2f}s  verify={ok}")

    print()
    print("OK")
