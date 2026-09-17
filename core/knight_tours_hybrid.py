"""knight_tours_hybrid.py — Gerador autocontido de tours fechados do cavalo.

Algoritmo:
    Backtracking + propagação R2 (vértice de grau 2) + detector de
    sub-ciclos por Union-Find incremental + heurística de fase local contínua.
"""

import sys
import math
import numpy as np

sys.setrecursionlimit(10000)

# Estados de aresta no array `fixed`
FREE = 0
ACTIVE = 1
INACTIVE = 2

# Códigos de retorno
OK = 0
CONTRADICTION = 1
SUBTOUR = 2
COMPLETE_TOUR = 3

_KNIGHT_MOVES = ((-2, -1), (-2, 1), (-1, -2), (-1, 2),
                 (1, -2), (1, 2), (2, -1), (2, 1))

# ---------------------------------------------------------------------------
# Componente 1 — Grafo do cavalo
# ---------------------------------------------------------------------------

def build_graph(n):
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

    return {
        'n': n, 'V': V, 'E': E,
        'edge_endpoints': edge_endpoints,
        'adj_edges': adj_edges,
        'total_incident': total_incident,
    }


def vertex_level(v, n):
    r, c = divmod(int(v), n)
    return min(r, c, n - 1 - r, n - 1 - c)


def edge_level(u, v, n):
    return min(vertex_level(u, n), vertex_level(v, n))


# ---------------------------------------------------------------------------
# Componente 2 — Heurística de fase local (Modelo Híbrido com LUT)
# ---------------------------------------------------------------------------

def _calc_f_inf_raw(L: int) -> float:
    if L == 0:
        return 0.528 
    A = -0.1985
    lambda_val = 1.1604
    return 0.25 + A * math.exp(-lambda_val * L)

# Pré-calcula a curva termodinâmica para todos os níveis possíveis até L=200 (n=400)
_F_INF_LUT = [_calc_f_inf_raw(L) for L in range(200)]

def edge_priority(L):
    f = _F_INF_LUT[L] if L < 200 else _calc_f_inf_raw(L)
    return abs(f - 0.5), (1 if f > 0.5 else 0)


# ---------------------------------------------------------------------------
# Estado da busca
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
# Componente 3a — Propagação R2 + detector de sub-ciclos
# ---------------------------------------------------------------------------

def _process_queue(state, ctx, queue):
    n2 = ctx['n'] * ctx['n']
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
                if state.uf_size[ru] == n2:
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
# Componente 3c — Verificação / extração de tour
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

    state.nodes_explored += 1

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

def knight_tours(n, K, seed=None):
    if n < 6:
        raise ValueError("n deve ser ≥ 6")
    if K < 0:
        raise ValueError("K deve ser não-negativo")
    if K == 0:
        return [], 0

    ctx = build_graph(n)
    rng = np.random.default_rng(seed)
    state = State(ctx['V'], ctx['E'])

    status = _propagate_initial(state, ctx)
    if status in (CONTRADICTION, SUBTOUR):
        return [], state.nodes_explored

    tours = []
    if status == COMPLETE_TOUR:
        if _is_complete_tour(state, ctx):
            tours.append(_extract_tour(state, ctx))
        return tours[:K], state.nodes_explored

    _backtrack(state, ctx, rng, tours, K)
    return tours, state.nodes_explored


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
