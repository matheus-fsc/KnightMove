#!/usr/bin/env python3
import sys
import time
import numpy as np

sys.setrecursionlimit(20000)

_KNIGHT_MOVES = ((-2, -1), (-2, 1), (-1, -2), (-1, 2),
                 (1, -2), (1, 2), (2, -1), (2, 1))

FREE = 0
ACTIVE = 1
INACTIVE = 2

OK = 0
CONTRADICTION = 1
SUBTOUR = 2
COMPLETE_TOUR = 3


def build_graph(n, m):
    V = n * m
    edge_set = set()
    for y in range(n):
        for x in range(m):
            u = y * m + x
            for dy, dx in _KNIGHT_MOVES:
                yy = (y + dy) % n
                xx = (x + dx) % m
                w = yy * m + xx
                if u == w:
                    continue
                a, b = (u, w) if u < w else (w, u)
                edge_set.add((a, b))
    edges = sorted(edge_set)
    E = len(edges)
    edge_endpoints = np.asarray(edges, dtype=np.int32)
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


def _choose_next_edge(state, ctx, anchor_v0=False):
    V = ctx['V']
    total_inc = ctx['total_incident']
    adj_edges = ctx['adj_edges']
    if anchor_v0:
        free_0 = int(total_inc[0]) - int(state.degree[0]) - int(state.n_inactive[0])
        if free_0 > 0:
            candidates = [0]
        else:
            candidates = []
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
            return _choose_next_edge(state, ctx, anchor_v0=False)
        return -1, 0
    # Deterministic choice: take the first candidate
    v_star = candidates[0]
    livres = [int(ei) for ei in adj_edges[v_star] if state.fixed[ei] == FREE]
    if not livres:
        return -1, 0
    # Deterministic choice of edge
    return livres[0], 0


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


def _backtrack_exhaustive(state, ctx, tours, depth=0):
    e, first_val = _choose_next_edge(state, ctx, anchor_v0=(depth == 0))
    if e == -1:
        if _is_complete_tour(state, ctx):
            tours.append(_extract_tour(state, ctx))
            if len(tours) % 5000 == 0:
                print(f"Encontrados {len(tours)} tours...")
        return

    # Try both values deterministically
    for val in (first_val, 1 - first_val):
        snap = state.snapshot()
        status = fix_and_propagate(state, ctx, e, val)
        if status == OK:
            _backtrack_exhaustive(state, ctx, tours, depth + 1)
        elif status == COMPLETE_TOUR:
            if _is_complete_tour(state, ctx):
                tours.append(_extract_tour(state, ctx))
                if len(tours) % 5000 == 0:
                    print(f"Encontrados {len(tours)} tours...")
        state.restore(snap)


def gf2_rank(M):
    if M.size == 0 or M.shape[0] == 0 or M.shape[1] == 0:
        return 0
    M = M.copy().astype(np.int64) % 2
    rows, cols = M.shape
    pivot_row = 0
    for col in range(cols):
        found = -1
        for row in range(pivot_row, rows):
            if M[row, col] & 1:
                found = row
                break
        if found == -1:
            continue
        M[[pivot_row, found]] = M[[found, pivot_row]]
        for row in range(rows):
            if row != pivot_row and M[row, col] & 1:
                M[row] ^= M[pivot_row]
        pivot_row += 1
    return pivot_row


def build_tour_matrix(tours, ctx):
    K = len(tours)
    E = ctx['E']
    V = ctx['V']
    edge_to_index = {}
    for ei in range(E):
        u, w = ctx['edge_endpoints'][ei]
        edge_to_index[(min(u,w), max(u,w))] = ei
    M = np.zeros((K, E), dtype=np.uint8)
    for k in range(K):
        tour = tours[k]
        for i in range(V):
            u = tour[i]
            w = tour[(i+1)%V]
            ei = edge_to_index[(min(u,w), max(u,w))]
            M[k, ei] = 1
    return M


if __name__ == "__main__":
    n = 6
    ctx = build_graph(n, n)
    state = State(ctx['V'], ctx['E'])
    _propagate_initial(state, ctx)
    
    t0 = time.time()
    tours = []
    print("Iniciando busca exaustiva no toro 6x6...")
    _backtrack_exhaustive(state, ctx, tours)
    elapsed = time.time() - t0
    
    num_tours = len(tours)
    print(f"Concluído! Total de tours distintos (com ancoragem): {num_tours} em {elapsed:.3f} segundos.")
    
    if num_tours > 0:
        M = build_tour_matrix(tours, ctx)
        t0_rank = time.time()
        rank = gf2_rank(M)
        print(f"Rank GF(2) exato de todos os tours: {rank} (tempo rank: {time.time()-t0_rank:.3f}s)")
        print(f"beta1(G_T(6)) = {ctx['E'] - ctx['V'] + 1} = 109")
        print(f"Deficit exato = 109 - {rank} = {109 - rank}")
