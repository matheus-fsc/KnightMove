"""knight_tours_cylinder.py — Gerador de tours fechados do cavalo em
topologia cilíndrica (b1 = 1).

Tabuleiro N (linhas, eixo Y rígido) × M (colunas, eixo X modular).
A borda direita conecta-se à esquerda; teto e chão permanecem rígidos.

Algoritmo:
    Backtracking + propagação R2 + detector de sub-ciclos por Union-Find
    incremental + heurística de fase local f∞(L) (agora L depende só do Y)
    + ancoragem da primeira aresta na coluna x=0 para quebrar a simetria
    de translação C_m.

Dependência única: numpy.

Interface pública:
    knight_tours(n, m, K, seed=None, log_nodes=False) -> list[np.ndarray]
    verify_tour(tour, n, m) -> bool
"""

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
# Componente 1 — Grafo do cavalo no cilindro
# ---------------------------------------------------------------------------

def build_graph(n, m):
    """Grafo do cavalo no cilindro N×M.

    Y ∈ [0, N) rígido; X ∈ [0, M) modular.
    Indexação 2D→1D:  v = y*m + x.

    Para M pequeno (M ≤ 4) dois deslocamentos distintos podem resultar
    no mesmo destino — usamos `set` para garantir grafo simples.
    """
    V = n * m
    edge_set = set()
    for y in range(n):
        for x in range(m):
            u = y * m + x
            for dy, dx in _KNIGHT_MOVES:
                yy = y + dy
                if not (0 <= yy < n):
                    continue
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


def vertex_level(v, n, m):
    """L(v) = min(y, n-1-y) — interioridade só pelo eixo Y (X é cíclico)."""
    y = int(v) // m
    return min(y, n - 1 - y)


def edge_level(u, v, n, m):
    return min(vertex_level(u, n, m), vertex_level(v, n, m))


# ---------------------------------------------------------------------------
# Componente 2 — Heurística de fase local
# ---------------------------------------------------------------------------

def edge_priority(L):
    f = F_INF.get(L, F_INF_DEFAULT)
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
    """Propaga R2 inicial. No cilindro, normalmente vazio (sem cantos
    de grau 2), mas mantemos a chamada por uniformidade."""
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

def _choose_next_edge(state, ctx, rng, restrict_column0=False):
    """Escolhe próxima aresta. Se `restrict_column0`, só considera vértices
    com x=0 (ancoragem da simetria de translação C_m). Faz fallback para
    todos os vértices se a coluna 0 já estiver saturada."""
    V = ctx['V']
    n = ctx['n']
    m = ctx['m']
    total_inc = ctx['total_incident']
    adj_edges = ctx['adj_edges']
    edge_endpoints = ctx['edge_endpoints']

    best_score = -1
    candidates = []
    for v in range(V):
        if restrict_column0 and (v % m) != 0:
            continue
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
        if restrict_column0:
            return _choose_next_edge(state, ctx, rng, restrict_column0=False)
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
        L = edge_level(a, b, n, m)
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

def _backtrack(state, ctx, rng, tours, K, depth=0):
    if len(tours) >= K:
        return True

    state.nodes_explored += 1

    e, first_val = _choose_next_edge(
        state, ctx, rng, restrict_column0=(depth == 0))

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

def knight_tours(n, m, K, seed=None, log_nodes=False):
    """Gera até K tours fechados no cilindro N×M.

    Parâmetros:
      n:         linhas (eixo Y rígido)
      m:         colunas (eixo X modular)
      K:         número de tours desejado
      seed:      semente RNG (opcional)
      log_nodes: se True, imprime contagem de nós explorados ao final

    Retorna:
      lista de até K arrays numpy (n*m,) dtype int32 — cada um é a
      sequência fechada de vértices (tour[-1] → tour[0] também é movimento
      de cavalo cilíndrico).
    """
    if n < 3:
        raise ValueError("n deve ser ≥ 3")
    if m < 3:
        raise ValueError("m deve ser ≥ 3")
    if K < 0:
        raise ValueError("K deve ser não-negativo")
    if K == 0:
        return []

    ctx = build_graph(n, m)
    rng = np.random.default_rng(seed)
    state = State(ctx['V'], ctx['E'])

    status = _propagate_initial(state, ctx)
    if status in (CONTRADICTION, SUBTOUR):
        if log_nodes:
            print(f"[cilindro {n}x{m}] inviável na propagação inicial "
                  f"(status={status})")
        return []

    tours = []
    if status == COMPLETE_TOUR:
        if _is_complete_tour(state, ctx):
            tours.append(_extract_tour(state, ctx))
    else:
        _backtrack(state, ctx, rng, tours, K)

    if log_nodes:
        per_tour = state.nodes_explored / max(len(tours), 1)
        print(f"[cilindro {n}x{m}] tours={len(tours)}/{K} "
              f"nós={state.nodes_explored} ({per_tour:.2f}/tour)")

    return tours[:K]


def verify_tour(tour, n, m):
    """Verifica tour fechado do cavalo no cilindro N×M.

    Critérios:
      - cobre todos os n*m vértices exatamente uma vez
      - cada passo consecutivo é movimento de cavalo válido com Y rígido
        e X modular (m)
      - é fechado
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
    for i in range(V):
        u = int(tour[i])
        w = int(tour[(i + 1) % V])
        uy, ux = divmod(u, m)
        wy, wx = divmod(w, m)
        dy = wy - uy
        diff = (wx - ux) % m
        # No cilindro a mesma aresta pode corresponder a dx = diff ou diff - m
        ok = False
        for dx in (diff, diff - m):
            if (dy, dx) in moves:
                ok = True
                break
        if not ok:
            return False
    return True


# ---------------------------------------------------------------------------
# Validador estrutural (critério de aceite 2)
# ---------------------------------------------------------------------------

def validate_2_factor_connected(state, ctx):
    """Confirma grau exatamente 2 em todos os vértices e que o 2-fator
    forma um único componente (tour fechado, não união de ciclos)."""
    V = ctx['V']
    if not np.all(state.degree == 2):
        return False, "algum vértice não tem grau 2"
    root0 = uf_find(state, 0)
    if int(state.uf_size[root0]) != V:
        return False, "2-fator tem mais de um componente"
    return True, "ok"


# ---------------------------------------------------------------------------
# Testes
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    import time

    # Teste 0 — sanidade do grafo cilíndrico 6×6
    ctx = build_graph(6, 6)
    print(f"[grafo cilindro 6x6] V={ctx['V']} E={ctx['E']}")
    grau_min = int(ctx['total_incident'].min())
    grau_max = int(ctx['total_incident'].max())
    print(f"  grau ∈ [{grau_min}, {grau_max}]  "
          f"(esperado: sem grau 2; lateral interior = 8)")
    assert grau_min >= 4, f"grau mínimo {grau_min} < 4 — cilindro deveria eliminar cantos"

    # Teste 1 — critério de aceite: n=m=6 cilíndrico, log de nós
    print()
    print("=== Critério 1: n=m=6 cilíndrico, K=10, log nós ===")
    t0 = time.time()
    tours = knight_tours(6, 6, 10, seed=0, log_nodes=True)
    t1 = time.time()
    print(f"  tempo: {t1 - t0:.3f}s")
    assert len(tours) == 10, f"esperava 10, obteve {len(tours)}"
    for t in tours:
        assert verify_tour(t, 6, 6), "tour cilíndrico inválido"
    print("  todos os 10 tours validados (cobertura + adjacência + fechamento)")

    # Teste 2 — critério de aceite: validador estrutural (grau=2, 1 componente)
    # Reconstroi o 2-fator a partir do tour extraído e valida.
    print()
    print("=== Critério 2: validador estrutural (grau=2, 1 componente) ===")
    ctx6 = build_graph(6, 6)
    tour6 = knight_tours(6, 6, 1, seed=0)[0]
    V6 = ctx6['V']
    deg = np.zeros(V6, dtype=np.int32)
    parent = np.arange(V6, dtype=np.int32)
    size = np.ones(V6, dtype=np.int32)

    def _find(p, v):
        while p[v] != v:
            v = int(p[v])
        return v

    for i in range(V6):
        u, w = int(tour6[i]), int(tour6[(i + 1) % V6])
        deg[u] += 1
        deg[w] += 1
        ru, rw = _find(parent, u), _find(parent, w)
        if ru != rw:
            if size[ru] < size[rw]:
                ru, rw = rw, ru
            parent[rw] = ru
            size[ru] += size[rw]

    grau_ok = bool(np.all(deg == 2))
    conex_ok = int(size[_find(parent, 0)]) == V6
    print(f"  grau exatamente 2 em todos vértices: {grau_ok}")
    print(f"  componente único (tour fechado válido): {conex_ok}")
    assert grau_ok and conex_ok

    # Teste 3 — reprodutibilidade
    print()
    print("=== Reprodutibilidade ===")
    a = knight_tours(6, 6, 5, seed=42)
    b = knight_tours(6, 6, 5, seed=42)
    assert len(a) == len(b)
    assert all(np.array_equal(x, y) for x, y in zip(a, b))
    print("  ok")

    # Teste 4 — diversidade
    print()
    print("=== Diversidade (n=m=6, K=20) ===")
    tours = knight_tours(6, 6, 20, seed=0)
    unique = {tuple(int(x) for x in t) for t in tours}
    print(f"  {len(unique)} únicos de {len(tours)}")
    assert len(unique) == len(tours)

    # Teste 5 — escalonamento curto
    print()
    print("=== Escalonamento ===")
    for (n, m) in [(6, 6), (6, 8), (8, 6), (8, 8), (6, 10)]:
        t0 = time.time()
        tours = knight_tours(n, m, 20, seed=0, log_nodes=True)
        t1 = time.time()
        ok = all(verify_tour(t, n, m) for t in tours)
        print(f"  {n}x{m}: {len(tours)} tours em {t1 - t0:.2f}s  verify={ok}")

    print()
    print("OK")
