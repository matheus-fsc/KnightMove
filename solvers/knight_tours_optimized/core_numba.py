"""core_numba.py — Numba JIT do inner loop _process_queue.

Estratégia:
  - process_queue_nb compilado com @njit, recebe arrays planos (CSR)
  - StateNumba é wrapper Python que só guarda os buffers
  - Backtracking permanece em Python (snapshot/restore via numpy.copyto)
  - choose_next_edge reaproveita o do core.py (já passa StateNumba sem mudança)
"""

import sys
import numpy as np

try:
    from numba import njit
    NUMBA_AVAILABLE = True
except ImportError:
    NUMBA_AVAILABLE = False
    def njit(*args, **kwargs):
        if len(args) == 1 and callable(args[0]):
            return args[0]
        def deco(f):
            return f
        return deco

from .core import (
    build_graph, choose_next_edge,
    FREE, ACTIVE, INACTIVE,
    OK, CONTRADICTION, SUBTOUR, COMPLETE,
)

sys.setrecursionlimit(10000)

# Tamanho da pilha LIFO do _process_queue.
# Bound conservador: cada aresta pode ser empurrada várias vezes
# (até 2 × max_deg ~ 16) antes de ser processada. Para n≤14 (E~625),
# 8192 dá margem 8× sobre o pior caso realista.
MAX_QUEUE = 8192


# ---------------------------------------------------------------------------
# T1 — CSR: adj_flat + adj_ptr para Numba
# ---------------------------------------------------------------------------

def build_graph_numba(n):
    """build_graph + representação CSR das adjacências."""
    ctx = build_graph(n)
    V = ctx['V']
    adj_lists = ctx['adj_edges']

    adj_ptr = np.zeros(V + 1, dtype=np.int32)
    for v in range(V):
        adj_ptr[v + 1] = adj_ptr[v] + len(adj_lists[v])

    adj_flat = np.zeros(int(adj_ptr[V]), dtype=np.int32)
    for v in range(V):
        start = int(adj_ptr[v])
        for i, ei in enumerate(adj_lists[v]):
            adj_flat[start + i] = int(ei)

    ctx['adj_flat'] = adj_flat
    ctx['adj_ptr'] = adj_ptr
    return ctx


# ---------------------------------------------------------------------------
# T2 — process_queue_nb (kernel JIT)
# ---------------------------------------------------------------------------

@njit(cache=True)
def uf_find_nb(uf_parent, v):
    while uf_parent[v] != v:
        v = uf_parent[v]
    return v


@njit(cache=True)
def process_queue_nb(
    fixed, degree, n_inactive,
    uf_parent, uf_size, uf_deg2,
    n_free,
    degree_target,
    edge_endpoints, total_inc,
    adj_flat, adj_ptr,
    queue_e, queue_v, q_top,
    n2, mode,
):
    """Numba LIFO de _process_queue.

    mode: 0 = cycle, 1 = path
    Retorna 0=OK 1=CONTRADICTION 2=SUBTOUR 3=COMPLETE
    """
    FREE_V = 0
    ACTIVE_V = 1
    INACTIVE_V = 2

    while q_top[0] > 0:
        q_top[0] -= 1
        e = queue_e[q_top[0]]
        val = queue_v[q_top[0]]
        cur = fixed[e]
        if cur != FREE_V:
            if (cur == ACTIVE_V and val == 1) or (cur == INACTIVE_V and val == 0):
                continue
            return 1

        u = edge_endpoints[e, 0]
        v = edge_endpoints[e, 1]

        if val == 1:
            fixed[e] = ACTIVE_V
            n_free[0] -= 1
            degree[u] += 1
            degree[v] += 1
            if degree[u] > degree_target[u] or degree[v] > degree_target[v]:
                return 1

            ru = uf_find_nb(uf_parent, u)
            rv = uf_find_nb(uf_parent, v)
            new_at_target = 0
            if degree[u] == degree_target[u]:
                new_at_target += 1
            if degree[v] == degree_target[v]:
                new_at_target += 1

            if ru == rv:
                if mode == 1:
                    return 1
                uf_deg2[ru] += new_at_target
                if uf_size[ru] == n2:
                    return 3
                return 2

            if uf_size[ru] < uf_size[rv]:
                tmp = ru
                ru = rv
                rv = tmp
            uf_parent[rv] = ru
            uf_size[ru] += uf_size[rv]
            uf_deg2[ru] += uf_deg2[rv] + new_at_target
        else:
            fixed[e] = INACTIVE_V
            n_free[0] -= 1
            n_inactive[u] += 1
            n_inactive[v] += 1
            avail_u = total_inc[u] - n_inactive[u]
            avail_v = total_inc[v] - n_inactive[v]
            need_u = degree_target[u] - degree[u]
            need_v = degree_target[v] - degree[v]
            if avail_u < need_u or avail_v < need_v:
                return 1

        # Propagação R-target em u e v
        for w_iter in range(2):
            if w_iter == 0:
                w = u
            else:
                w = v
            dw = degree[w]
            tw = degree_target[w]
            avail = total_inc[w] - n_inactive[w]
            if dw == tw:
                for i in range(adj_ptr[w], adj_ptr[w + 1]):
                    ei = adj_flat[i]
                    if fixed[ei] == FREE_V:
                        queue_e[q_top[0]] = ei
                        queue_v[q_top[0]] = 0
                        q_top[0] += 1
            elif avail == tw and dw < tw:
                for i in range(adj_ptr[w], adj_ptr[w + 1]):
                    ei = adj_flat[i]
                    if fixed[ei] == FREE_V:
                        queue_e[q_top[0]] = ei
                        queue_v[q_top[0]] = 1
                        q_top[0] += 1

    return 0


# ---------------------------------------------------------------------------
# T3 — StateNumba + wrapper Python
# ---------------------------------------------------------------------------

class StateNumba:
    __slots__ = ('fixed', 'degree', 'n_inactive',
                 'uf_parent', 'uf_size', 'uf_deg2', 'n_free',
                 'degree_target',
                 'queue_e', 'queue_v', 'q_top')

    def __init__(self, V, E, degree_target=None):
        self.fixed      = np.zeros(E, dtype=np.uint8)
        self.degree     = np.zeros(V, dtype=np.int32)
        self.n_inactive = np.zeros(V, dtype=np.int32)
        self.uf_parent  = np.arange(V, dtype=np.int32)
        self.uf_size    = np.ones(V, dtype=np.int32)
        self.uf_deg2    = np.zeros(V, dtype=np.int32)
        self.n_free     = np.array([E], dtype=np.int32)
        if degree_target is None:
            self.degree_target = np.full(V, 2, dtype=np.int32)
        else:
            self.degree_target = np.asarray(degree_target, dtype=np.int32)
        self.queue_e = np.zeros(MAX_QUEUE, dtype=np.int32)
        self.queue_v = np.zeros(MAX_QUEUE, dtype=np.int32)
        self.q_top   = np.zeros(1, dtype=np.int32)

    def snapshot(self):
        return (self.fixed.copy(), self.degree.copy(),
                self.n_inactive.copy(), self.uf_parent.copy(),
                self.uf_size.copy(), self.uf_deg2.copy(),
                int(self.n_free[0]))

    def restore(self, snap):
        np.copyto(self.fixed, snap[0])
        np.copyto(self.degree, snap[1])
        np.copyto(self.n_inactive, snap[2])
        np.copyto(self.uf_parent, snap[3])
        np.copyto(self.uf_size, snap[4])
        np.copyto(self.uf_deg2, snap[5])
        self.n_free[0] = snap[6]
        self.q_top[0] = 0


def _run_pq(state, ctx, mode):
    return int(process_queue_nb(
        state.fixed, state.degree, state.n_inactive,
        state.uf_parent, state.uf_size, state.uf_deg2,
        state.n_free, state.degree_target,
        ctx['edge_endpoints'], ctx['total_incident'],
        ctx['adj_flat'], ctx['adj_ptr'],
        state.queue_e, state.queue_v, state.q_top,
        np.int32(ctx['V']), np.int32(mode),
    ))


def fix_and_propagate_nb(state, ctx, e, val, mode='cycle'):
    state.q_top[0] = 0
    state.queue_e[0] = int(e)
    state.queue_v[0] = int(val)
    state.q_top[0] = 1
    return _run_pq(state, ctx, 0 if mode == 'cycle' else 1)


def propagate_initial_nb(state, ctx, mode='cycle'):
    V = ctx['V']
    total_inc = ctx['total_incident']
    adj_edges = ctx['adj_edges']
    target = state.degree_target
    state.q_top[0] = 0
    for v in range(V):
        avail = int(total_inc[v]) - int(state.n_inactive[v])
        if avail == int(target[v]) and int(state.degree[v]) < int(target[v]):
            for ei in adj_edges[v]:
                ei_i = int(ei)
                if state.fixed[ei_i] == FREE:
                    state.queue_e[state.q_top[0]] = ei_i
                    state.queue_v[state.q_top[0]] = 1
                    state.q_top[0] += 1
    return _run_pq(state, ctx, 0 if mode == 'cycle' else 1)


def _uf_find_py(state, v):
    while state.uf_parent[v] != v:
        v = int(state.uf_parent[v])
    return v


def _is_complete_tour_nb(state, ctx):
    V = ctx['V']
    if not np.all(state.degree == 2):
        return False
    root0 = _uf_find_py(state, 0)
    return int(state.uf_size[root0]) == V


def _extract_tour_nb(state, ctx):
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


def _backtrack_nb(state, ctx, rng, tours, K):
    if len(tours) >= K:
        return True
    e, first_val = choose_next_edge(state, ctx, rng)
    if e == -1:
        if _is_complete_tour_nb(state, ctx):
            tours.append(_extract_tour_nb(state, ctx))
        return len(tours) >= K
    for val in (first_val, 1 - first_val):
        snap = state.snapshot()
        status = fix_and_propagate_nb(state, ctx, e, val, mode='cycle')
        if status == OK:
            if _backtrack_nb(state, ctx, rng, tours, K):
                state.restore(snap)
                return True
        elif status == COMPLETE:
            if _is_complete_tour_nb(state, ctx):
                tours.append(_extract_tour_nb(state, ctx))
            if len(tours) >= K:
                state.restore(snap)
                return True
        state.restore(snap)
    return len(tours) >= K


def knight_tours_numba(n, K, seed=None):
    """API compatível com knight_tours mas usando process_queue_nb."""
    if n < 6:
        raise ValueError("n deve ser ≥ 6")
    if K <= 0:
        return []
    ctx = build_graph_numba(n)
    rng = np.random.default_rng(seed)
    state = StateNumba(ctx['V'], ctx['E'])
    status = propagate_initial_nb(state, ctx, mode='cycle')
    if status in (CONTRADICTION, SUBTOUR):
        return []
    tours = []
    if status == COMPLETE:
        if _is_complete_tour_nb(state, ctx):
            tours.append(_extract_tour_nb(state, ctx))
        return tours[:K]
    _backtrack_nb(state, ctx, rng, tours, K)
    return tours


# ---------------------------------------------------------------------------
# T4 — microbenchmark
# ---------------------------------------------------------------------------

def benchmark_numba_vs_numpy(n=10, K=500):
    import time
    from .tours import knight_tours as kt_numpy

    print(f"Compilando Numba (warmup) n={n}...")
    knight_tours_numba(n, 5, seed=0)

    t0 = time.perf_counter()
    kt_numpy(n, K, seed=0)
    t_np = time.perf_counter() - t0

    t0 = time.perf_counter()
    knight_tours_numba(n, K, seed=0)
    t_nb = time.perf_counter() - t0

    print(f"  NumPy:   {t_np:.3f}s  ({K / max(t_np, 1e-9):.0f} tours/s)")
    print(f"  Numba:   {t_nb:.3f}s  ({K / max(t_nb, 1e-9):.0f} tours/s)")
    print(f"  Speedup: {t_np / max(t_nb, 1e-9):.2f}×")
    return t_np, t_nb
