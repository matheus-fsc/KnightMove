"""knight_tours_sheared.py — Gerador de tours fechados do cavalo em um
Toro com Cisalhamento (Sheared Torus).

Tabuleiro N (linhas, eixo Y modular) × M (colunas, eixo X modular).
Colagem: (x, y) ~ (x + m, y + s*wrap_x mod n), onde wrap_x ∈ {0, 1}
documenta quantas vezes cruzou a fronteira X.

Sem fronteira rígida: todo vértice tem grau 8 (para N,M ≥ 5).

Algoritmo:
    Backtracking + propagação R2 + detector Union-Find incremental.
    A heurística f∞(L) degenera (não há borda) — priority uniforme.
    Ancoragem em v=0 para quebrar simetria de translação C_n × C_m.

Dependência única: numpy.

Interface pública:
    knight_tours(n, m, s, K, seed=None, log_nodes=False) -> list[np.ndarray]
    verify_tour(tour, n, m, s) -> bool
"""

import sys
import numpy as np

sys.setrecursionlimit(10000)

F_INF_DEFAULT = 0.25  # Toro com cisalhamento: tudo é "interior"

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
# Componente 1 — Grafo do cavalo no toro com cisalhamento
# ---------------------------------------------------------------------------

def sheared_step(y, x, dy, dx, n, m, s):
    """Aplica movimento (dy, dx) a (y, x) em um toro com cisalhamento.
    
    Computação:
      1. wrap_x = (x + dx) // m   (quantas vezes cruzou fronteira X)
      2. new_x = (x + dx) % m
      3. new_y = (y + dy + wrap_x * s) % n
    
    Retorna (new_y, new_x, wrap_x).
    """
    wrap_x = (x + dx) // m
    new_x = (x + dx) % m
    new_y = (y + dy + wrap_x * s) % n
    return new_y, new_x, wrap_x


def build_graph(n, m, s):
    """Constrói grafo de cavalo no toro com cisalhamento s.
    
    Parâmetros:
        n: número de linhas (Y modular)
        m: número de colunas (X modular)
        s: parâmetro de cisalhamento (deslocamento em Y por cruzamento X)
    
    Retorna dicionário com metadados do grafo.
    """
    V = n * m
    edge_set = set()
    
    for y in range(n):
        for x in range(m):
            u = y * m + x
            for dy, dx in _KNIGHT_MOVES:
                new_y, new_x, _wrap_x = sheared_step(y, x, dy, dx, n, m, s)
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
        'n': n, 'm': m, 's': s, 'V': V, 'E': E,
        'edge_endpoints': edge_endpoints,
        'adj_edges': adj_edges,
        'total_incident': total_incident,
    }


# ---------------------------------------------------------------------------
# Componente 2 — Estado da busca
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
# Componente 3 — Propagação R2 + detector
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
    """Normalmente vazia em toro com cisalhamento (sem grau 2 forçado)."""
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
# Componente 4 — Seleção (ancoragem em v=0 no root)
# ---------------------------------------------------------------------------

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
    first_val = 0 if F_INF_DEFAULT < 0.5 else 1
    return ei, first_val


# ---------------------------------------------------------------------------
# Componente 5 — verify / extract / backtrack
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
# Interface pública
# ---------------------------------------------------------------------------

def knight_tours(n, m, s, K, seed=None, log_nodes=False):
    """Gera até K tours fechados no toro N×M com cisalhamento s.
    
    Parâmetros:
        n: número de linhas
        m: número de colunas
        s: parâmetro de cisalhamento (mod n)
        K: número máximo de tours para gerar
        seed: seed do RNG
        log_nodes: se True, imprime stats de exploração
    
    Retorna:
        lista de tours (cada um é um np.ndarray de tamanho n*m)
    """
    if n < 3 or m < 3:
        raise ValueError("n, m devem ser ≥ 3")
    if K < 0:
        raise ValueError("K não-negativo")
    if K == 0:
        return []

    ctx = build_graph(n, m, s)
    rng = np.random.default_rng(seed)
    state = State(ctx['V'], ctx['E'])

    status = _propagate_initial(state, ctx)
    if status in (CONTRADICTION, SUBTOUR):
        if log_nodes:
            print(f"[toro {n}x{m}, s={s}] inviável (status={status})")
        return []

    tours = []
    if status == COMPLETE_TOUR:
        if _is_complete_tour(state, ctx):
            tours.append(_extract_tour(state, ctx))
    else:
        _backtrack(state, ctx, rng, tours, K)

    if log_nodes:
        per_tour = state.nodes_explored / max(len(tours), 1)
        print(f"[toro {n}x{m}, s={s}] tours={len(tours)}/{K} "
              f"nós={state.nodes_explored} ({per_tour:.2f}/tour)")

    return tours[:K]


def verify_tour(tour, n, m, s):
    """Verifica se tour é um passeio do cavalo válido no toro com cisalhamento.
    
    Parâmetros:
        tour: lista/array de vértices
        n: número de linhas
        m: número de colunas
        s: parâmetro de cisalhamento
    
    Retorna:
        True se o tour é válido, False caso contrário.
    """
    V = n * m
    if len(tour) != V:
        return False
    seen = set()
    for x in tour:
        xi = int(x)
        if xi < 0 or xi >= V or xi in seen:
            return False
        seen.add(xi)

    moves = set(_KNIGHT_MOVES)
    sum_dy = 0
    sum_dx = 0
    for i in range(V):
        u = int(tour[i])
        w = int(tour[(i + 1) % V])
        uy, ux = divmod(u, m)
        wy, wx = divmod(w, m)
        
        found = False
        for dy, dx in moves:
            ny, nx, _wr = sheared_step(uy, ux, dy, dx, n, m, s)
            if (ny, nx) == (wy, wx):
                sum_dy += dy
                sum_dx += dx
                found = True
                break
        if not found:
            return False
    
    # Sanity checks: must close on sheared torus
    if sum_dx % m != 0:
        return False
    adjusted_y = sum_dy + s * (sum_dx // m)
    if adjusted_y % n != 0:
        return False


    return True
