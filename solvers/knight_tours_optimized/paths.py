"""paths.py — knight_path: caminhos Hamiltonianos abertos start→end."""

import sys
import numpy as np

from .core import (
    State, OK, CONTRADICTION, SUBTOUR, COMPLETE, ACTIVE,
    KNIGHT_MOVES,
    build_graph, fix_and_propagate, propagate_initial, choose_next_edge,
    uf_find,
)

sys.setrecursionlimit(10000)


def _vertex_color(v, n):
    r, c = divmod(int(v), n)
    return (r + c) % 2


def _is_complete_path(state, ctx, start, end):
    V = ctx['V']
    if not np.array_equal(state.degree, state.degree_target):
        return False
    root_s = uf_find(state, start)
    return int(state.uf_size[root_s]) == V


def _extract_path(state, ctx, start, end):
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
    if len(nbrs[start]) != 1:
        raise RuntimeError(
            f"start={start} tem grau {len(nbrs[start])}, esperava 1")
    if len(nbrs[end]) != 1:
        raise RuntimeError(
            f"end={end} tem grau {len(nbrs[end])}, esperava 1")
    path = np.zeros(V, dtype=np.int32)
    path[0] = start
    prev = -1
    cur = start
    for i in range(1, V):
        cands = [w for w in nbrs[cur] if w != prev]
        if len(cands) != 1:
            raise RuntimeError(
                f"posição {i}: vértice {cur} tem {len(cands)} candidatos")
        nxt = cands[0]
        path[i] = nxt
        prev = cur
        cur = nxt
    if cur != end:
        raise RuntimeError(f"caminho terminou em {cur}, esperava {end}")
    return path


def _backtrack_path(state, ctx, rng, paths, K, start, end):
    if len(paths) >= K:
        return True
    e, first_val = choose_next_edge(state, ctx, rng)
    if e == -1:
        if _is_complete_path(state, ctx, start, end):
            paths.append(_extract_path(state, ctx, start, end))
        return len(paths) >= K
    for val in (first_val, 1 - first_val):
        snap = state.snapshot()
        status = fix_and_propagate(state, ctx, e, val, mode='path')
        if status == OK:
            if _backtrack_path(state, ctx, rng, paths, K, start, end):
                state.restore(snap)
                return True
        # path mode: COMPLETE não vem da propagação (não fecha ciclo)
        state.restore(snap)
    return len(paths) >= K


def knight_path(n, start, end, K=1, seed=None):
    """Gera até K caminhos Hamiltonianos abertos de start a end no tabuleiro n×n.

    start e end devem estar em [0, n²) e devem ter cores compatíveis
    com a bipartição do grafo do cavalo (V par → cores opostas;
    V ímpar → mesma cor). Caso contrário, retorna [].

    Raises:
      ValueError se start == end ou se start/end fora de [0, n²).
    """
    if n < 4:
        raise ValueError("n deve ser ≥ 4")
    V = n * n
    if start == end:
        raise ValueError("start deve ser diferente de end")
    if not (0 <= start < V):
        raise ValueError(f"start deve estar em [0, {V}), obteve {start}")
    if not (0 <= end < V):
        raise ValueError(f"end deve estar em [0, {V}), obteve {end}")
    if K < 0:
        raise ValueError("K deve ser não-negativo")
    if K == 0:
        return []

    cs = _vertex_color(start, n)
    ce = _vertex_color(end, n)
    if V % 2 == 0:
        if cs == ce:
            print(f"Aviso: start({start}) e end({end}) têm a mesma cor — "
                  f"nenhum caminho aberto Hamiltoniano existe.")
            return []
    else:
        if cs != ce:
            print(f"Aviso: start({start}) e end({end}) têm cores diferentes — "
                  f"nenhum caminho aberto Hamiltoniano existe.")
            return []

    ctx = build_graph(n)
    rng = np.random.default_rng(seed)
    degree_target = np.full(V, 2, dtype=np.int32)
    degree_target[start] = 1
    degree_target[end] = 1
    state = State(ctx['V'], ctx['E'], degree_target=degree_target)

    status = propagate_initial(state, ctx, mode='path')
    if status in (CONTRADICTION, SUBTOUR):
        return []
    paths = []
    if status == COMPLETE:
        if _is_complete_path(state, ctx, start, end):
            paths.append(_extract_path(state, ctx, start, end))
        return paths[:K]
    _backtrack_path(state, ctx, rng, paths, K, start, end)
    return paths


def verify_path(path, n, start, end):
    """Verifica se path é Hamiltoniano válido start→end no tabuleiro n×n."""
    V = n * n
    if len(path) != V:
        return False
    if int(path[0]) != int(start):
        return False
    if int(path[-1]) != int(end):
        return False
    seen = set()
    for x in path:
        xi = int(x)
        if xi < 0 or xi >= V or xi in seen:
            return False
        seen.add(xi)
    moves = set(KNIGHT_MOVES)
    for i in range(V - 1):
        u = int(path[i])
        w = int(path[i + 1])
        ru, cu = divmod(u, n)
        rw, cw = divmod(w, n)
        if (rw - ru, cw - cu) not in moves:
            return False
    return True
