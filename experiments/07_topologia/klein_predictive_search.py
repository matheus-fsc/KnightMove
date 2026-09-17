"""Busca exaustiva com Poda Preditiva e Heurística Topológica na Garrafa de Klein.

Implementação das otimizações:
1. Poda Preditiva (Look-ahead Bounding): aborta caminhos abertos que possuem 
   Winding Number ímpar mas não têm passos suficientes para alcançar a borda X.
2. Heurística Topológica: ordena/prioriza arestas candidatas no desempate da
   escolha de Warnsdorff para acelerar a correção ou preservação da topologia.
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
    predictive_prunes: int = 0
    subtour_prunes: int = 0
    contradictions: int = 0


def vertex_id(y: int, x: int, m: int) -> int:
    return y * m + x


def coords(v: int, m: int) -> tuple[int, int]:
    return divmod(v, m)


def klein_step(y: int, x: int, dy: int, dx: int, n: int, m: int) -> tuple[int, int, int] | None:
    wrap_x = (x + dx) // m
    nx = (x + dx) % m
    
    if wrap_x != 0:
        if abs(wrap_x) % 2 != 0:
            ny = (n - 1) - (y + dy)
        else:
            ny = y + dy
    else:
        ny = y + dy
        
    if ny < 0 or ny >= n:
        return None
        
    return ny, nx, wrap_x


def build_graph(n: int, m: int) -> dict[str, object]:
    valid_edges = set()
    dir_wx = {}
    
    for y in range(n):
        for x in range(m):
            u = vertex_id(y, x, m)
            for dy, dx in KNIGHT_MOVES:
                step_res = klein_step(y, x, dy, dx, n, m)
                if step_res is not None:
                    ny, nx, wrap_x = step_res
                    v = vertex_id(ny, nx, m)
                    if u != v:
                        valid_edges.add((u, v))
                        dir_wx[(u, v)] = wrap_x

    edges_sym = set()
    for (u, v) in valid_edges:
        if (v, u) in valid_edges:
            edges_sym.add(tuple(sorted([u, v])))
            
    edges = sorted(edges_sym)
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
        "V": n * m,
        "E": len(edges),
        "edges": edges,
        "edge_index": edge_index,
        "adj_edges": adj_edges,
        "total_incident": [len(lst) for lst in adj_edges],
        "directed_wx": dir_wx,
    }


class State:
    __slots__ = (
        "fixed", "degree", "inactive", "n_free",
        "uf_parent", "uf_size", "end_a", "end_b", "path_wx"
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
        self.path_wx = [0] * V

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
            self.path_wx[:],
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
            self.path_wx,
        ) = snap


def uf_find(state: State, v: int) -> int:
    while state.uf_parent[v] != v:
        v = state.uf_parent[v]
    return v


def endpoint_path_to(state: State, root: int, endpoint: int) -> tuple[int, int]:
    a = state.end_a[root]
    b = state.end_b[root]
    total_wx = state.path_wx[root]
    if endpoint == a:
        return b, -total_wx
    if endpoint == b:
        return a, total_wx
    raise RuntimeError(f"vertice {endpoint} nao e extremo da componente {root}")


def path_sum_between(state: State, root: int, start: int, end: int) -> int:
    a = state.end_a[root]
    b = state.end_b[root]
    total_wx = state.path_wx[root]
    if start == a and end == b:
        return total_wx
    if start == b and end == a:
        return -total_wx
    if start == end and state.uf_size[root] == 1:
        return 0
    raise RuntimeError(f"extremos incompativeis: {start}->{end} na raiz {root}")


def directed_wx(ctx: dict[str, object], u: int, v: int) -> int:
    return ctx["directed_wx"][(u, v)]  # type: ignore[index]


def union_or_close_path(state: State, ctx: dict[str, object], u: int, v: int) -> int:
    V = int(ctx["V"])
    ru = uf_find(state, u)
    rv = uf_find(state, v)
    wx_uv = directed_wx(ctx, u, v)

    if ru == rv:
        if state.uf_size[ru] < V:
            path_wx = path_sum_between(state, ru, u, v)
            cycle_wx = path_wx + directed_wx(ctx, v, u)
            if cycle_wx % 2 != 0 and ctx.get("enable_topo_prune", True):
                return TOPO_PRUNE
            return SUBTOUR

        path_wx = path_sum_between(state, ru, u, v)
        cycle_wx = path_wx + directed_wx(ctx, v, u)
        if cycle_wx % 2 != 0 and ctx.get("enable_topo_prune", True):
            return TOPO_PRUNE
        return COMPLETE_TOUR

    if u not in (state.end_a[ru], state.end_b[ru]):
        return CONTRADICTION
    if v not in (state.end_a[rv], state.end_b[rv]):
        return CONTRADICTION

    other_u, sum_wx_other_u_to_u = endpoint_path_to(state, ru, u)
    other_v, sum_wx_v_to_other_v = endpoint_path_to(state, rv, v)
    new_wx = sum_wx_other_u_to_u + wx_uv + sum_wx_v_to_other_v

    if state.uf_size[ru] < state.uf_size[rv]:
        ru, rv = rv, ru
    state.uf_parent[rv] = ru
    state.uf_size[ru] += state.uf_size[rv]
    state.end_a[ru] = other_u
    state.end_b[ru] = other_v
    state.path_wx[ru] = new_wx
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

    # Look-ahead Bounding Topológico (Poda Preditiva)
    if ctx.get("enable_lookahead", True):
        V = int(ctx["V"])
        m = int(ctx["m"])
        for r in range(V):
            if state.uf_parent[r] == r:
                if state.path_wx[r] % 2 != 0:
                    d = state.uf_size[r]
                    xa = state.end_a[r] % m
                    xb = state.end_b[r] % m
                    
                    # Distância mínima em passos de cavalo (eixo X) até a borda
                    jumps_a = min(xa // 2 + 1, (m - 1 - xa) // 2 + 1)
                    jumps_b = min(xb // 2 + 1, (m - 1 - xb) // 2 + 1)
                    min_jumps = min(jumps_a, jumps_b)
                    
                    # Condição de poda preditiva
                    limit = V - d
                    if ctx.get("use_strict_limit", False):
                        # Permite o wrap na última aresta de fechamento se limit + 1
                        limit = V - d + 1
                        
                    if min_jumps > limit:
                        metrics.predictive_prunes += 1
                        return TOPO_PRUNE

    return OK


def choose_next_edge(state: State, ctx: dict[str, object], depth: int) -> int:
    adj_edges = ctx["adj_edges"]  # type: ignore[assignment]
    total_incident = ctx["total_incident"]  # type: ignore[assignment]
    edges = ctx["edges"]  # type: ignore[assignment]
    V = int(ctx["V"])
    m = int(ctx["m"])

    if depth == 0:
        free_anchor = [ei for ei in adj_edges[0] if state.fixed[ei] == FREE]
        if free_anchor:
            return free_anchor[0]

    best_score = -1
    best_vertex = -1
    for v in range(V):
        free = total_incident[v] - state.inactive[v] - state.degree[v]
        if free <= 0 or state.degree[v] >= 2:
            continue
        score = state.degree[v] * 100 + (total_incident[v] - state.inactive[v] - 2)
        if score > best_score:
            best_score = score
            best_vertex = v

    if best_vertex == -1:
        return -1

    # Heurística Topológica para ordenar/desempatar arestas incidentes a best_vertex
    if ctx.get("enable_heuristic", True):
        best_edge = -1
        best_edge_score = -999999
        rv = uf_find(state, best_vertex)
        
        for ei in adj_edges[best_vertex]:
            if state.fixed[ei] == FREE:
                u, w = edges[ei]
                nbr = w if u == best_vertex else u
                
                is_wrap = (ctx["directed_wx"].get((best_vertex, nbr), 0) != 0)
                xv = best_vertex % m
                x_nbr = nbr % m
                dist_v = min(xv, m - 1 - xv)
                dist_nbr = min(x_nbr, m - 1 - x_nbr)
                
                topo_weight = 0
                if state.path_wx[rv] % 2 != 0:
                    # Paridade ímpar: prioriza wrap ou aproximar da borda
                    if is_wrap:
                        topo_weight += 20
                    elif dist_nbr < dist_v:
                        topo_weight += 10
                    elif dist_nbr > dist_v:
                        topo_weight -= 10
                else:
                    # Paridade par: penaliza wrap desnecessário
                    if is_wrap:
                        topo_weight -= 20
                
                if topo_weight > best_edge_score:
                    best_edge_score = topo_weight
                    best_edge = ei
        return best_edge
    else:
        # Padrão: escolhe a primeira aresta livre
        for ei in adj_edges[best_vertex]:
            if state.fixed[ei] == FREE:
                return ei
        return -1


def active_edge_key(state: State, ctx: dict[str, object]) -> tuple[tuple[int, int], ...]:
    edges = ctx["edges"]  # type: ignore[assignment]
    return tuple(edges[i] for i, value in enumerate(state.fixed) if value == ACTIVE)


def translate_edge_key(edge_key: tuple[tuple[int, int], ...], n: int, m: int):
    translated = []
    for dy in range(n):
        cur_t = []
        for u, v in edge_key:
            uy, ux = coords(u, m)
            vy, vx = coords(v, m)
            tu = vertex_id((uy + dy) % n, ux, m)
            tv = vertex_id((vy + dy) % n, vx, m)
            cur_t.append((tu, tv) if tu < tv else (tv, tu))
        translated.append(tuple(sorted(cur_t)))
    return min(translated)


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
    e = choose_next_edge(state, ctx, depth)
    if e < 0:
        if is_complete_tour(state, ctx):
            metrics.tours += 1
            key = translate_edge_key(active_edge_key(state, ctx), int(ctx["n"]), int(ctx["m"]))
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
            key = translate_edge_key(active_edge_key(state, ctx), int(ctx["n"]), int(ctx["m"]))
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
    parser = argparse.ArgumentParser(description="Busca na Garrafa de Klein com Poda Preditiva e Heurística Topológica.")
    parser.add_argument("--n", type=int, default=6)
    parser.add_argument("--m", type=int, default=6)
    parser.add_argument("--max-tours", type=int, default=None)
    parser.add_argument("--node-limit", type=int, default=None)
    parser.add_argument("--disable-lookahead", action="store_true", help="Desativa Poda Preditiva")
    parser.add_argument("--disable-heuristic", action="store_true", help="Desativa Heurística Topológica")
    parser.add_argument("--use-strict-limit", action="store_true", help="Usa V - d + 1 em vez de V - d para limite estrito")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    t0 = time.perf_counter()
    ctx = build_graph(args.n, args.m)
    
    ctx["enable_lookahead"] = not args.disable_lookahead
    ctx["enable_heuristic"] = not args.disable_heuristic
    ctx["use_strict_limit"] = args.use_strict_limit
    
    state = State(int(ctx["V"]), int(ctx["E"]))
    metrics = Metrics()
    canonical_seen: set[tuple[tuple[int, int], ...]] = set()

    stopped = search(state, ctx, metrics, canonical_seen, 0, args.max_tours, args.node_limit)
    elapsed = time.perf_counter() - t0

    print("=" * 72)
    print(f"=== Otimizações Topológicas na Garrafa de Klein ===")
    print(f"Tabuleiro: {args.n}x{args.m}")
    print(f"Poda Preditiva (Look-ahead)           : {'ATIVADA' if ctx['enable_lookahead'] else 'DESATIVADA'}")
    print(f"Limite estrito (V - d + 1)           : {'ATIVADA' if ctx['use_strict_limit'] else 'DESATIVADA'}")
    print(f"Heurística Topológica                : {'ATIVADA' if ctx['enable_heuristic'] else 'DESATIVADA'}")
    print("=" * 72)
    print(f"Tempo total de execucao              : {elapsed:.3f}s")
    print(f"Nos visitados no backtracking        : {metrics.nodes}")
    print(f"Tours rotulados validos encontrados  : {metrics.tours}")
    print(f"Tours canonicos validos encontrados  : {metrics.canonical_tours}")
    print(f"Poda Topológica Clássica (Ciclo)     : {metrics.topo_prunes}")
    print(f"Poda Preditiva Ativa (Look-ahead)    : {metrics.predictive_prunes}")
    print(f"Subciclos rejeitados pelo UF         : {metrics.subtour_prunes}")
    print(f"Contradicoes locais                  : {metrics.contradictions}")
    if stopped:
        print("Status                               : interrompido por limite configurado")
    else:
        print("Status                               : busca exaustiva concluida")


if __name__ == "__main__":
    main()
