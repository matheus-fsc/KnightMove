"""Busca exaustiva de ciclos Hamiltonianos no toro cisalhado 6x6, s=3.

O estado e um problema de 2-fator sobre as arestas do grafo do cavalo. A
propagacao local fixa arestas quando um vertice fica com exatamente duas
opcoes restantes, e invalida arestas incidentes quando o grau chega a 2.

O Union-Find guarda componentes que sao caminhos: extremos e soma orientada
de dy real entre eles. Quando uma aresta fechar uma componente, o detector
calcula a paridade de W_y e conta uma poda topologica se ela for impar.
"""

from __future__ import annotations

import argparse
import time
from dataclasses import dataclass


FREE = 0
ACTIVE = 1
INACTIVE = 2

OK = 0
CONTRADICTION = 1
SUBTOUR = 2
COMPLETE_TOUR = 3
TOPO_PRUNE = 4

KNIGHT_MOVES = (
    (-2, -1), (-2, 1), (-1, -2), (-1, 2),
    (1, -2), (1, 2), (2, -1), (2, 1),
)


@dataclass
class Metrics:
    nodes: int = 0
    tours: int = 0
    canonical_tours: int = 0
    topo_prunes: int = 0
    subtour_prunes: int = 0
    contradictions: int = 0


def vertex_id(y: int, x: int, m: int) -> int:
    return y * m + x


def coords(v: int, m: int) -> tuple[int, int]:
    return divmod(v, m)


def sheared_step(y: int, x: int, dy: int, dx: int, n: int, m: int, s: int) -> tuple[int, int]:
    wrap_x = (x + dx) // m
    return (y + dy + wrap_x * s) % n, (x + dx) % m


def build_graph(n: int, m: int, s: int) -> dict[str, object]:
    edge_set = set()
    directed_dy = {}
    directed_dx = {}

    for y in range(n):
        for x in range(m):
            u = vertex_id(y, x, m)
            for dy, dx in KNIGHT_MOVES:
                ny, nx = sheared_step(y, x, dy, dx, n, m, s)
                v = vertex_id(ny, nx, m)
                if u == v:
                    continue
                edge_set.add((u, v) if u < v else (v, u))
                directed_dy[(u, v)] = dy
                directed_dx[(u, v)] = dx

    edges = sorted(edge_set)
    edge_index = {edge: i for i, edge in enumerate(edges)}
    adj_edges = [[] for _ in range(n * m)]
    for i, (u, v) in enumerate(edges):
        adj_edges[u].append(i)
        adj_edges[v].append(i)

    for lst in adj_edges:
        lst.sort()

    return {
        "n": n,
        "m": m,
        "s": s,
        "V": n * m,
        "E": len(edges),
        "edges": edges,
        "edge_index": edge_index,
        "adj_edges": adj_edges,
        "total_incident": [len(lst) for lst in adj_edges],
        "directed_dy": directed_dy,
        "directed_dx": directed_dx,
    }


class State:
    __slots__ = (
        "fixed", "degree", "inactive", "n_free",
        "uf_parent", "uf_size", "end_a", "end_b", "path_dy", "path_dx",
    )

    def __init__(self, V: int, E: int):
        self.fixed = bytearray(E)
        self.degree = [0] * V
        self.inactive = [0] * V
        self.n_free = E
        self.uf_parent = list(range(V))
        self.uf_size = [1] * V
        self.end_a = list(range(V))
        self.end_b = list(range(V))
        self.path_dy = [0] * V
        self.path_dx = [0] * V

    def snapshot(self):
        return (
            self.fixed[:],
            self.degree[:],
            self.inactive[:],
            self.n_free,
            self.uf_parent[:],
            self.uf_size[:],
            self.end_a[:],
            self.end_b[:],
            self.path_dy[:],
            self.path_dx[:],
        )

    def restore(self, snap) -> None:
        (
            self.fixed,
            self.degree,
            self.inactive,
            self.n_free,
            self.uf_parent,
            self.uf_size,
            self.end_a,
            self.end_b,
            self.path_dy,
            self.path_dx,
        ) = snap


def uf_find(state: State, v: int) -> int:
    while state.uf_parent[v] != v:
        v = state.uf_parent[v]
    return v


def endpoint_path_to(state: State, root: int, endpoint: int) -> tuple[int, int, int]:
    """Retorna (outro_extremo, soma_dy, soma_dx) outro_extremo -> endpoint."""
    a = state.end_a[root]
    b = state.end_b[root]
    total_dy = state.path_dy[root]
    total_dx = state.path_dx[root]
    if endpoint == a:
        return b, -total_dy, -total_dx
    if endpoint == b:
        return a, total_dy, total_dx
    raise RuntimeError(f"vertice {endpoint} nao e extremo da componente {root}")


def path_sum_between(state: State, root: int, start: int, end: int) -> tuple[int, int]:
    a = state.end_a[root]
    b = state.end_b[root]
    total_dy = state.path_dy[root]
    total_dx = state.path_dx[root]
    if start == a and end == b:
        return total_dy, total_dx
    if start == b and end == a:
        return -total_dy, -total_dx
    if start == end and state.uf_size[root] == 1:
        return 0, 0
    raise RuntimeError(f"extremos incompativeis: {start}->{end} na raiz {root}")


def directed_dy(ctx: dict[str, object], u: int, v: int) -> int:
    return ctx["directed_dy"][(u, v)]  # type: ignore[index]


def directed_dx(ctx: dict[str, object], u: int, v: int) -> int:
    return ctx["directed_dx"][(u, v)]  # type: ignore[index]


def cycle_has_odd_wy(ctx: dict[str, object], cycle_dy: int, cycle_dx: int) -> bool | None:
    n = int(ctx["n"])
    m = int(ctx["m"])
    s = int(ctx["s"])
    if cycle_dx % m != 0:
        return None
    adjusted_y = cycle_dy + s * (cycle_dx // m)
    if adjusted_y % n != 0:
        return None
    return (adjusted_y // n) % 2 != 0


def union_or_close_path(state: State, ctx: dict[str, object], u: int, v: int) -> int:
    V = int(ctx["V"])
    ru = uf_find(state, u)
    rv = uf_find(state, v)
    dy_uv = directed_dy(ctx, u, v)
    dx_uv = directed_dx(ctx, u, v)

    if ru == rv:
        if state.uf_size[ru] < V:
            path_dy, path_dx = path_sum_between(state, ru, u, v)
            cycle_dy = path_dy + directed_dy(ctx, v, u)
            cycle_dx = path_dx + directed_dx(ctx, v, u)
            if cycle_has_odd_wy(ctx, cycle_dy, cycle_dx):
                return TOPO_PRUNE
            return SUBTOUR

        path_dy, path_dx = path_sum_between(state, ru, u, v)
        cycle_dy = path_dy + directed_dy(ctx, v, u)
        cycle_dx = path_dx + directed_dx(ctx, v, u)
        if cycle_has_odd_wy(ctx, cycle_dy, cycle_dx):
            return TOPO_PRUNE
        return COMPLETE_TOUR

    if u not in (state.end_a[ru], state.end_b[ru]):
        return CONTRADICTION
    if v not in (state.end_a[rv], state.end_b[rv]):
        return CONTRADICTION

    other_u, sum_dy_other_u_to_u, sum_dx_other_u_to_u = endpoint_path_to(state, ru, u)
    other_v, sum_dy_v_to_other_v, sum_dx_v_to_other_v = endpoint_path_to(state, rv, v)
    new_dy = sum_dy_other_u_to_u + dy_uv + sum_dy_v_to_other_v
    new_dx = sum_dx_other_u_to_u + dx_uv + sum_dx_v_to_other_v

    if state.uf_size[ru] < state.uf_size[rv]:
        ru, rv = rv, ru
    state.uf_parent[rv] = ru
    state.uf_size[ru] += state.uf_size[rv]
    state.end_a[ru] = other_u
    state.end_b[ru] = other_v
    state.path_dy[ru] = new_dy
    state.path_dx[ru] = new_dx
    return OK


def process_queue(state: State, ctx: dict[str, object], metrics: Metrics, queue: list[tuple[int, int]]) -> int:
    edges = ctx["edges"]  # type: ignore[assignment]
    adj_edges = ctx["adj_edges"]  # type: ignore[assignment]
    total_incident = ctx["total_incident"]  # type: ignore[assignment]

    while queue:
        e, val = queue.pop()
        current = state.fixed[e]
        if current != FREE:
            if (current == ACTIVE and val == 1) or (current == INACTIVE and val == 0):
                continue
            metrics.contradictions += 1
            return CONTRADICTION

        u, v = edges[e]
        state.fixed[e] = ACTIVE if val == 1 else INACTIVE
        state.n_free -= 1

        if val == 1:
            state.degree[u] += 1
            state.degree[v] += 1
            if state.degree[u] > 2 or state.degree[v] > 2:
                metrics.contradictions += 1
                return CONTRADICTION

            status = union_or_close_path(state, ctx, u, v)
            if status == TOPO_PRUNE:
                metrics.topo_prunes += 1
                return TOPO_PRUNE
            if status == SUBTOUR:
                metrics.subtour_prunes += 1
                return SUBTOUR
            if status == CONTRADICTION:
                metrics.contradictions += 1
                return CONTRADICTION
            if status == COMPLETE_TOUR:
                return COMPLETE_TOUR
        else:
            state.inactive[u] += 1
            state.inactive[v] += 1
            if total_incident[u] - state.inactive[u] < 2:
                metrics.contradictions += 1
                return CONTRADICTION
            if total_incident[v] - state.inactive[v] < 2:
                metrics.contradictions += 1
                return CONTRADICTION

        for w in (u, v):
            if state.degree[w] == 2:
                for ei in adj_edges[w]:
                    if state.fixed[ei] == FREE:
                        queue.append((ei, 0))
            elif total_incident[w] - state.inactive[w] == 2 and state.degree[w] < 2:
                for ei in adj_edges[w]:
                    if state.fixed[ei] == FREE:
                        queue.append((ei, 1))

    return OK


def choose_next_edge(state: State, ctx: dict[str, object], depth: int) -> int:
    adj_edges = ctx["adj_edges"]  # type: ignore[assignment]
    total_incident = ctx["total_incident"]  # type: ignore[assignment]
    V = int(ctx["V"])

    if depth == 0:
        free_anchor = [ei for ei in adj_edges[0] if state.fixed[ei] == FREE]
        if free_anchor:
            return free_anchor[0]

    best_score = -1
    best_edge = -1
    for v in range(V):
        free = total_incident[v] - state.inactive[v] - state.degree[v]
        if free <= 0 or state.degree[v] >= 2:
            continue
        score = state.degree[v] * 100 + (total_incident[v] - state.inactive[v] - 2)
        if score > best_score:
            for ei in adj_edges[v]:
                if state.fixed[ei] == FREE:
                    best_score = score
                    best_edge = ei
                    break
    return best_edge


def active_edge_key(state: State, ctx: dict[str, object]) -> tuple[tuple[int, int], ...]:
    edges = ctx["edges"]  # type: ignore[assignment]
    return tuple(edges[i] for i, value in enumerate(state.fixed) if value == ACTIVE)


def translate_edge_key(edge_key: tuple[tuple[int, int], ...], n: int, m: int, dy: int, dx: int):
    translated = []
    for u, v in edge_key:
        uy, ux = coords(u, m)
        vy, vx = coords(v, m)
        tu = vertex_id((uy + dy) % n, (ux + dx) % m, m)
        tv = vertex_id((vy + dy) % n, (vx + dx) % m, m)
        translated.append((tu, tv) if tu < tv else (tv, tu))
    return tuple(sorted(translated))


def canonical_translation_key(state: State, ctx: dict[str, object]) -> tuple[tuple[int, int], ...]:
    n = int(ctx["n"])
    m = int(ctx["m"])
    key = active_edge_key(state, ctx)
    return min(translate_edge_key(key, n, m, dy, dx) for dy in range(n) for dx in range(m))


def is_complete_tour(state: State, ctx: dict[str, object]) -> bool:
    V = int(ctx["V"])
    return all(deg == 2 for deg in state.degree) and state.uf_size[uf_find(state, 0)] == V


def search(
    state: State,
    ctx: dict[str, object],
    metrics: Metrics,
    canonical_seen: set[tuple[tuple[int, int], ...]],
    depth: int,
    max_tours: int | None,
    node_limit: int | None,
) -> bool:
    if node_limit is not None and metrics.nodes >= node_limit:
        return True
    if max_tours is not None and metrics.canonical_tours >= max_tours:
        return True

    metrics.nodes += 1
    if metrics.nodes % 100000 == 0:
        print(f"[Progresso] Nos visitados: {metrics.nodes:,} | Tours canonicos: {metrics.canonical_tours:,}", flush=True)
    e = choose_next_edge(state, ctx, depth)
    if e < 0:
        if is_complete_tour(state, ctx):
            metrics.tours += 1
            key = canonical_translation_key(state, ctx)
            if key not in canonical_seen:
                canonical_seen.add(key)
                metrics.canonical_tours += 1
        return False

    for val in (1, 0):
        snap = state.snapshot()
        status = process_queue(state, ctx, metrics, [(e, val)])
        if status == OK:
            stop = search(state, ctx, metrics, canonical_seen, depth + 1, max_tours, node_limit)
        elif status == COMPLETE_TOUR and is_complete_tour(state, ctx):
            metrics.tours += 1
            key = canonical_translation_key(state, ctx)
            if key not in canonical_seen:
                canonical_seen.add(key)
                metrics.canonical_tours += 1
            stop = max_tours is not None and metrics.canonical_tours >= max_tours
        else:
            stop = False
        state.restore(snap)
        if stop:
            return True

    return False


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Busca completa no toro cisalhado com poda topologica.")
    parser.add_argument("--n", type=int, default=6)
    parser.add_argument("--m", type=int, default=6)
    parser.add_argument("--s", type=int, default=3)
    parser.add_argument("--max-tours", type=int, default=None, help="limite opcional para smoke tests")
    parser.add_argument("--node-limit", type=int, default=None, help="limite opcional de nos visitados")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    t0 = time.perf_counter()
    ctx = build_graph(args.n, args.m, args.s)
    state = State(int(ctx["V"]), int(ctx["E"]))
    metrics = Metrics()
    canonical_seen: set[tuple[tuple[int, int], ...]] = set()

    stopped = search(state, ctx, metrics, canonical_seen, 0, args.max_tours, args.node_limit)
    elapsed = time.perf_counter() - t0

    print("=" * 72)
    print(f"Busca completa - toro cisalhado {args.n}x{args.m}, s={args.s}")
    print("=" * 72)
    print(f"Vertices (V)                         : {ctx['V']}")
    print(f"Arestas (E)                          : {ctx['E']}")
    print(f"Tempo total de execucao              : {elapsed:.3f}s")
    print(f"Nos visitados no backtracking         : {metrics.nodes}")
    print(f"Tours rotulados validos encontrados   : {metrics.tours}")
    print(f"Tours canonicos validos encontrados   : {metrics.canonical_tours}")
    print(f"Galhos podados por topologia W_y impar: {metrics.topo_prunes}")
    print(f"Subciclos rejeitados pelo UF           : {metrics.subtour_prunes}")
    print(f"Contradicoes locais                    : {metrics.contradictions}")
    if stopped:
        print("Status                               : interrompido por limite configurado")
    else:
        print("Status                               : busca exaustiva concluida")


if __name__ == "__main__":
    main()
