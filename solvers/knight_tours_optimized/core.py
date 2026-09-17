"""core.py — Estado, grafo do cavalo, e propagação R-target compartilhados.

Generaliza knight_tours.py:
  - degree_target heterogêneo (suporta start/end com target=1 para caminhos abertos)
  - modo 'cycle' para tours fechados; modo 'path' para caminhos abertos
    (em path, fechar qualquer ciclo é sempre contradição)
"""

import numpy as np


# Heurística de fase local f∞(L)
F_INF = {0: 0.528, 1: 0.193, 2: 0.198, 3: 0.294, 4: 0.247, 5: 0.261}
F_INF_DEFAULT = 0.25

# Estados de aresta
FREE = 0
ACTIVE = 1
INACTIVE = 2

# Códigos de retorno
OK = 0
CONTRADICTION = 1
SUBTOUR = 2
COMPLETE = 3

KNIGHT_MOVES = ((-2, -1), (-2, 1), (-1, -2), (-1, 2),
                (1, -2), (1, 2), (2, -1), (2, 1))


def build_graph(n):
    """Constrói o grafo do cavalo n×n."""
    V = n * n
    edges = []
    for r in range(n):
        for c in range(n):
            u = r * n + c
            for dr, dc in KNIGHT_MOVES:
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


def edge_priority(L):
    f = F_INF.get(L, F_INF_DEFAULT)
    return abs(f - 0.5), (1 if f > 0.5 else 0)


class State:
    """Estado mutável da busca. degree_target padrão = 2 (tour fechado)."""

    __slots__ = ('fixed', 'degree', 'n_inactive',
                 'uf_parent', 'uf_size', 'uf_deg2', 'n_free',
                 'degree_target')

    def __init__(self, V, E, degree_target=None):
        self.fixed = np.zeros(E, dtype=np.uint8)
        self.degree = np.zeros(V, dtype=np.int32)
        self.n_inactive = np.zeros(V, dtype=np.int32)
        self.uf_parent = np.arange(V, dtype=np.int32)
        self.uf_size = np.ones(V, dtype=np.int32)
        self.uf_deg2 = np.zeros(V, dtype=np.int32)
        self.n_free = E
        if degree_target is None:
            self.degree_target = np.full(V, 2, dtype=np.int32)
        else:
            self.degree_target = np.asarray(degree_target, dtype=np.int32)

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


def _process_queue(state, ctx, queue, mode='cycle'):
    """Processa fila de fixações com propagação R-target completa.

    mode='cycle': fechar ciclo = SUBTOUR (ou COMPLETE se cobrir todos)
    mode='path' : fechar ciclo = CONTRADICTION sempre (caminho aberto é acíclico)
    """
    n2 = ctx['n'] * ctx['n']
    total_inc = ctx['total_incident']
    adj_edges = ctx['adj_edges']
    edge_endpoints = ctx['edge_endpoints']
    target = state.degree_target

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
            if state.degree[u] > target[u] or state.degree[v] > target[v]:
                return CONTRADICTION

            ru = uf_find(state, u)
            rv = uf_find(state, v)
            new_at_target = (1 if state.degree[u] == target[u] else 0) + \
                            (1 if state.degree[v] == target[v] else 0)

            if ru == rv:
                if mode == 'path':
                    return CONTRADICTION
                state.uf_deg2[ru] += new_at_target
                if state.uf_size[ru] == n2:
                    return COMPLETE
                return SUBTOUR

            if state.uf_size[ru] < state.uf_size[rv]:
                ru, rv = rv, ru
            state.uf_parent[rv] = ru
            state.uf_size[ru] += state.uf_size[rv]
            state.uf_deg2[ru] += state.uf_deg2[rv] + new_at_target
        else:
            state.fixed[e] = INACTIVE
            state.n_free -= 1
            state.n_inactive[u] += 1
            state.n_inactive[v] += 1
            avail_u = total_inc[u] - state.n_inactive[u]
            avail_v = total_inc[v] - state.n_inactive[v]
            need_u = target[u] - state.degree[u]
            need_v = target[v] - state.degree[v]
            if avail_u < need_u or avail_v < need_v:
                return CONTRADICTION

        # R-target nos dois endpoints
        for w in (u, v):
            dw = int(state.degree[w])
            tw = int(target[w])
            avail = int(total_inc[w] - state.n_inactive[w])
            if dw == tw:
                for ei in adj_edges[w]:
                    if state.fixed[ei] == FREE:
                        queue.append((int(ei), 0))
            elif avail == tw and dw < tw:
                for ei in adj_edges[w]:
                    if state.fixed[ei] == FREE:
                        queue.append((int(ei), 1))

    return OK


def fix_and_propagate(state, ctx, e, val, mode='cycle'):
    return _process_queue(state, ctx, [(e, val)], mode=mode)


def propagate_initial(state, ctx, mode='cycle'):
    V = ctx['V']
    total_inc = ctx['total_incident']
    adj_edges = ctx['adj_edges']
    target = state.degree_target
    queue = []
    for v in range(V):
        avail = int(total_inc[v]) - int(state.n_inactive[v])
        if avail == int(target[v]) and int(state.degree[v]) < int(target[v]):
            for ei in adj_edges[v]:
                if state.fixed[ei] == FREE:
                    queue.append((int(ei), 1))
    return _process_queue(state, ctx, queue, mode=mode)


def choose_next_edge(state, ctx, rng):
    """Pressão de vértice: degree alto + maior 'slack' (avail - target).
    Tie-break aleatório entre candidatos com mesmo score."""
    V = ctx['V']
    n = ctx['n']
    total_inc = ctx['total_incident']
    adj_edges = ctx['adj_edges']
    edge_endpoints = ctx['edge_endpoints']
    target = state.degree_target

    best_score = -10**9
    candidates = []
    for v in range(V):
        free_v = int(total_inc[v]) - int(state.degree[v]) - int(state.n_inactive[v])
        if free_v == 0:
            continue
        score = int(state.degree[v]) * 100 + \
                int(total_inc[v] - state.n_inactive[v] - target[v])
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
