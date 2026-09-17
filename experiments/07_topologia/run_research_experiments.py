#!/usr/bin/env python3
"""
Knight's Tour Toroidal Topology and Wrap Descent Analysis
=========================================================
This script runs Experiments A, B, and C to verify the Toroidal Tightness Conjecture
and the Descent Theorem (Theorem U).
"""

import sys
import time
import numpy as np

# Set recursion limit for backtracking
sys.setrecursionlimit(10000)

# Knight's Moves
_KNIGHT_MOVES = ((-2, -1), (-2, 1), (-1, -2), (-1, 2),
                 (1, -2), (1, 2), (2, -1), (2, 1))

FREE = 0
ACTIVE = 1
INACTIVE = 2

OK = 0
CONTRADICTION = 1
SUBTOUR = 2
COMPLETE_TOUR = 3


# ---------------------------------------------------------------------------
# Graph Construction & Helpers
# ---------------------------------------------------------------------------

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


def classify_wrap(u, v, n, m):
    uy, ux = u // m, u % m
    wy, wx = v // m, v % m
    diff_x = abs(wx - ux)
    diff_y = abs(wy - uy)
    
    wrap_x = diff_x not in (1, 2)
    wrap_y = diff_y not in (1, 2)
    
    if wrap_x and wrap_y:
        return 'diag'
    elif wrap_x:
        return 'x-only'
    elif wrap_y:
        return 'y-only'
    else:
        return 'internal'


def get_corner_wraps(n, m):
    ctx = build_graph(n, m)
    edges = ctx['edge_endpoints']
    corner_wraps = []
    for ei in range(ctx['E']):
        u, v = edges[ei]
        if u == 0 or v == 0:
            other = v if u == 0 else u
            oy, ox = other // m, other % m
            is_internal = (abs(oy) == 1 and abs(ox) == 2) or (abs(oy) == 2 and abs(ox) == 1)
            if not is_internal:
                corner_wraps.append(ei)
    return corner_wraps


# ---------------------------------------------------------------------------
# Solver State & Propagation (Modified from knight_tours_torus.py)
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


def _choose_next_edge(state, ctx, rng, anchor_v0=False):
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
            return _choose_next_edge(state, ctx, rng, anchor_v0=False)
        return -1, 0

    if len(candidates) == 1:
        v_star = candidates[0]
    else:
        v_star = int(candidates[int(rng.integers(len(candidates)))])

    livres = [int(ei) for ei in adj_edges[v_star] if state.fixed[ei] == FREE]
    if not livres:
        return -1, 0
    ei = livres[int(rng.integers(len(livres)))]
    # Use INACTIVE first because torus has a high density of non-tour edges
    first_val = 0
    return ei, first_val


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
            raise RuntimeError(f"Vertex {v} has {len(nbrs[v])} active edges (expected 2)")

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


def set_edges_inactive(state, ctx, inactive_edges):
    for e in inactive_edges:
        u = int(ctx['edge_endpoints'][e, 0])
        v = int(ctx['edge_endpoints'][e, 1])
        if state.fixed[e] == FREE:
            state.fixed[e] = INACTIVE
            state.n_free -= 1
            state.n_inactive[u] += 1
            state.n_inactive[v] += 1


def knight_tours_with_exclusions(n, m, K, inactive_edges=None, seed=None, log_nodes=False):
    ctx = build_graph(n, m)
    rng = np.random.default_rng(seed)
    state = State(ctx['V'], ctx['E'])

    if inactive_edges is not None:
        set_edges_inactive(state, ctx, inactive_edges)

    status = _propagate_initial(state, ctx)
    if status in (CONTRADICTION, SUBTOUR):
        if log_nodes:
            print(f"[toro {n}x{m}] inviável na inicialização (status={status})")
        return []

    tours = []
    if status == COMPLETE_TOUR:
        if _is_complete_tour(state, ctx):
            tours.append(_extract_tour(state, ctx))
    else:
        _backtrack(state, ctx, rng, tours, K)

    return tours[:K]


# ---------------------------------------------------------------------------
# Algebraic GF(2) Rank and Q calculations
# ---------------------------------------------------------------------------

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


def compute_Q_for_experiment(n, m, inactive_edges):
    ctx = build_graph(n, m)
    E = ctx['E']
    V = ctx['V']
    
    active_mask = np.ones(E, dtype=bool)
    if inactive_edges is not None:
        for e in inactive_edges:
            active_mask[e] = False
            
    active_idx = np.where(active_mask)[0]
    num_active_edges = len(active_idx)
    if num_active_edges == 0:
        return 0
        
    amap = {old: new for new, old in enumerate(active_idx)}
    degree = np.zeros(V, dtype=int)
    adj = {v: [] for v in range(V)}
    for idx in active_idx:
        u, w = ctx['edge_endpoints'][idx]
        degree[u] += 1
        degree[w] += 1
        adj[u].append(idx)
        adj[w].append(idx)
        
    deg2 = [v for v in range(V) if degree[v] == 2]
    
    inc = np.zeros((V, num_active_edges), dtype=np.int64)
    for new, old in enumerate(active_idx):
        u, w = ctx['edge_endpoints'][old]
        inc[u, new] = 1
        inc[w, new] = 1
        
    rank_inc = gf2_rank(inc)
    
    if len(deg2) < 2:
        return 0
        
    forced = [adj[v][0] for v in deg2]
    xvecs = []
    for i in range(len(forced)):
        for j in range(i+1, len(forced)):
            vec = np.zeros(num_active_edges, dtype=np.int64)
            vec[amap[forced[i]]] = 1
            vec[amap[forced[j]]] = 1
            xvecs.append(vec)
            
    if not xvecs:
        return 0
        
    combined = np.vstack([inc, np.array(xvecs, dtype=np.int64)])
    Q = gf2_rank(combined) - rank_inc
    return Q


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


# ---------------------------------------------------------------------------
# EXPERIMENT A: Tightness no toro 6×6
# ---------------------------------------------------------------------------

def run_experiment_a():
    print("\n" + "="*80)
    print(" EXPERIMENTO A: Tightness no toro 6×6")
    print("="*80)
    
    n = 6
    K = 500
    seed = 42
    
    t0 = time.time()
    tours = knight_tours_with_exclusions(n, n, K, seed=seed)
    elapsed = time.time() - t0
    
    num_tours = len(tours)
    print(f"1. Coletados K = {num_tours} tours em {elapsed:.3f} segundos.")
    
    ctx = build_graph(n, n)
    M = build_tour_matrix(tours, ctx)
    
    t0_rank = time.time()
    rank = gf2_rank(M)
    elapsed_rank = time.time() - t0_rank
    
    beta1 = ctx['E'] - ctx['V'] + 1 # 144 - 36 + 1 = 109
    deficit = beta1 - rank
    
    print(f"2. Rank GF(2) calculado = {rank} (tempo rank: {elapsed_rank:.3f}s).")
    print(f"3. beta1(G_T(6)) = {beta1}")
    print(f"4. deficit_toro = beta1 - rank = {deficit}")
    print(f"5. deficit_plano = 3")
    
    is_tight = (deficit == 0)
    print(f"\nRESPOSTA DA PERGUNTA:")
    if is_tight:
        print("deficit_toro = 0. A conjectura de Tightness Toroidal é VERDADEIRA para n=6!")
    else:
        print(f"deficit_toro = {deficit} > 0. Existem obstruções adicionais no toro para n=6!")
        
    return deficit, elapsed


# ---------------------------------------------------------------------------
# EXPERIMENT B: Descida wrap-a-wrap
# ---------------------------------------------------------------------------

def run_experiment_b():
    print("\n" + "="*80)
    print(" EXPERIMENTO B: Descida wrap-a-wrap")
    print("="*80)
    
    n = 6
    K = 200
    seed = 100
    
    # 1. Identify and sort wrap edges of (0,0)
    corner_wraps = get_corner_wraps(n, n)
    ctx = build_graph(n, n)
    
    diag_wraps = []
    x_only_wraps = []
    y_only_wraps = []
    
    for ei in corner_wraps:
        u, v = ctx['edge_endpoints'][ei]
        cat = classify_wrap(u, v, n, n)
        if cat == 'diag':
            diag_wraps.append(ei)
        elif cat == 'x-only':
            x_only_wraps.append(ei)
        elif cat == 'y-only':
            y_only_wraps.append(ei)
            
    # Sort: 2 diag + 2 x-only + 2 y-only
    sorted_corner_wraps = diag_wraps + x_only_wraps + y_only_wraps
    
    print("Wraps do canto (0,0) identificadas e ordenadas:")
    for idx, ei in enumerate(sorted_corner_wraps):
        u, v = ctx['edge_endpoints'][ei]
        cat = classify_wrap(u, v, n, n)
        print(f"  e{idx+1}: {u // n, u % n} -> {v // n, v % n} ({cat})")
        
    # Table header
    print("\nTabela de Descida:")
    print(f"| k | |E_k| | β₁(k) | Q(k) | rank(k) | deficit(k) |")
    print(f"|---|---|---|---|---|---|")
    
    results = []
    t0_all = time.time()
    
    for k in range(7):
        removed_edges = sorted_corner_wraps[:k]
        num_edges = ctx['E'] - k
        beta1_k = num_edges - ctx['V'] + 1
        
        # Sample tours
        tours = knight_tours_with_exclusions(n, n, K, inactive_edges=removed_edges, seed=seed)
        
        if len(tours) == 0:
            rank_k = 0
        else:
            M = build_tour_matrix(tours, ctx)
            rank_k = gf2_rank(M)
            
        deficit_k = beta1_k - rank_k
        Q_k = compute_Q_for_experiment(n, n, removed_edges)
        
        print(f"| {k} | {num_edges} | {beta1_k} | {Q_k} | {rank_k} | {deficit_k} |")
        results.append((k, num_edges, beta1_k, Q_k, rank_k, deficit_k))
        
    elapsed = time.time() - t0_all
    print(f"\nTempo total do Experimento B: {elapsed:.3f} segundos.")
    
    # Check if deficit(k) = Q(k) for all steps
    match = all(r[3] == r[5] for r in results)
    
    print("\nRESPOSTA DA PERGUNTA:")
    if match:
        print("deficit(k) = Q(k) ao longo de TODA a descida! A tightness é perfeitamente preservada a cada passo de remoção de wraps.")
    else:
        print("deficit(k) != Q(k) em algum passo. Reportar a discrepância!")
        
    return results, elapsed


# ---------------------------------------------------------------------------
# EXPERIMENT C: Tightness toroidal n=8
# ---------------------------------------------------------------------------

def run_experiment_c():
    print("\n" + "="*80)
    print(" EXPERIMENTO C: Tightness toroidal n=8")
    print("="*80)
    
    n = 8
    K = 300
    seed = 88
    
    # Test generation speed and scale if needed
    t0 = time.time()
    tours = knight_tours_with_exclusions(n, n, K, seed=seed)
    elapsed = time.time() - t0
    
    num_tours = len(tours)
    print(f"1. Coletados K = {num_tours} tours em {elapsed:.3f} segundos.")
    
    ctx = build_graph(n, n)
    M = build_tour_matrix(tours, ctx)
    
    t0_rank = time.time()
    rank = gf2_rank(M)
    elapsed_rank = time.time() - t0_rank
    
    beta1 = ctx['E'] - ctx['V'] + 1 # 256 - 64 + 1 = 193
    deficit = beta1 - rank
    
    print(f"2. Rank GF(2) calculado = {rank} (tempo rank: {elapsed_rank:.3f}s).")
    print(f"3. beta1(G_T(8)) = {beta1}")
    print(f"4. deficit_toro(8) = {deficit}")
    
    return deficit, elapsed


# ---------------------------------------------------------------------------
# Main Execution
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    t_start = time.time()
    
    deficit6, tA = run_experiment_a()
    resB, tB = run_experiment_b()
    deficit8, tC = run_experiment_c()
    
    print("\n" + "="*80)
    print(" TABELA RESUMO DE TIGHTNESS")
    print("="*80)
    print(f"| n | β₁(plano) | deficit(plano) | β₁(toro) | deficit(toro) |")
    print(f"|---|---|---|---|---|")
    print(f"| 6 |     45    |       3        |    109    |      {deficit6}      |")
    print(f"| 8 |    105    |       3        |    193    |      {deficit8}      |")
    
    print(f"\nTempo total de execução dos 3 experimentos: {time.time() - t_start:.3f} segundos.")
