"""tours.py — knight_tours: tours fechados (refator de knight_tours.py)."""

import sys
import numpy as np

from .core import (
    State, OK, CONTRADICTION, SUBTOUR, COMPLETE, ACTIVE,
    KNIGHT_MOVES,
    build_graph, fix_and_propagate, propagate_initial, choose_next_edge,
    uf_find,
)

sys.setrecursionlimit(10000)


def _is_complete_tour(state, ctx):
    V = ctx['V']
    if not np.all(state.degree == 2):
        return False
    root0 = uf_find(state, 0)
    return int(state.uf_size[root0]) == V


def _extract_tour(state, ctx):
    """Sequência canônica iniciando em 0, indo ao menor vizinho ativo."""
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


def _backtrack(state, ctx, rng, tours, K):
    if len(tours) >= K:
        return True
    e, first_val = choose_next_edge(state, ctx, rng)
    if e == -1:
        if _is_complete_tour(state, ctx):
            tours.append(_extract_tour(state, ctx))
        return len(tours) >= K
    for val in (first_val, 1 - first_val):
        snap = state.snapshot()
        status = fix_and_propagate(state, ctx, e, val, mode='cycle')
        if status == OK:
            if _backtrack(state, ctx, rng, tours, K):
                state.restore(snap)
                return True
        elif status == COMPLETE:
            if _is_complete_tour(state, ctx):
                tours.append(_extract_tour(state, ctx))
            if len(tours) >= K:
                state.restore(snap)
                return True
        state.restore(snap)
    return len(tours) >= K


def knight_tours(n, K, seed=None):
    """Gera até K tours fechados do cavalo no tabuleiro n×n."""
    if n < 6:
        raise ValueError("n deve ser ≥ 6")
    if K < 0:
        raise ValueError("K deve ser não-negativo")
    if K == 0:
        return []
    ctx = build_graph(n)
    rng = np.random.default_rng(seed)
    state = State(ctx['V'], ctx['E'])
    status = propagate_initial(state, ctx, mode='cycle')
    if status in (CONTRADICTION, SUBTOUR):
        return []
    tours = []
    if status == COMPLETE:
        if _is_complete_tour(state, ctx):
            tours.append(_extract_tour(state, ctx))
        return tours[:K]
    _backtrack(state, ctx, rng, tours, K)
    return tours


def verify_tour(tour, n):
    """Verifica se tour é um tour fechado válido do cavalo em n×n."""
    V = n * n
    if len(tour) != V:
        return False
    seen = set()
    for x in tour:
        xi = int(x)
        if xi < 0 or xi >= V or xi in seen:
            return False
        seen.add(xi)
    moves = set(KNIGHT_MOVES)
    for i in range(V):
        u = int(tour[i])
        w = int(tour[(i + 1) % V])
        ru, cu = divmod(u, n)
        rw, cw = divmod(w, n)
        if (rw - ru, cw - cu) not in moves:
            return False
    return True
