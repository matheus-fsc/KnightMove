"""Script to analyze translational symmetries of tours in the 6x6, s=3 sheared torus.

It runs a backtracking search up to a configurable node limit, collects the
unique canonical tours, and computes their orbit size under the 36 translations.
Any tour with an orbit size < 36 has non-trivial translational symmetry.
"""

from __future__ import annotations
import sys
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
    return ctx["directed_dy"][(u, v)]

def directed_dx(ctx: dict[str, object], u: int, v: int) -> int:
    return ctx["directed_dx"][(u, v)]

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

def process_queue(state: State, ctx: dict[str, object], queue: list[tuple[int, int]]) -> int:
    edges = ctx["edges"]
    adj_edges = ctx["adj_edges"]
    total_incident = ctx["total_incident"]

    while queue:
        e, val = queue.pop()
        current = state.fixed[e]
        if current != FREE:
            if (current == ACTIVE and val == 1) or (current == INACTIVE and val == 0):
                continue
            return CONTRADICTION

        u, v = edges[e]
        state.fixed[e] = ACTIVE if val == 1 else INACTIVE
        state.n_free -= 1

        if val == 1:
            state.degree[u] += 1
            state.degree[v] += 1
            if state.degree[u] > 2 or state.degree[v] > 2:
                return CONTRADICTION

            status = union_or_close_path(state, ctx, u, v)
            if status in (TOPO_PRUNE, SUBTOUR, CONTRADICTION, COMPLETE_TOUR):
                return status
        else:
            state.inactive[u] += 1
            state.inactive[v] += 1
            if total_incident[u] - state.inactive[u] < 2:
                return CONTRADICTION
            if total_incident[v] - state.inactive[v] < 2:
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
    adj_edges = ctx["adj_edges"]
    total_incident = ctx["total_incident"]
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
    edges = ctx["edges"]
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

@dataclass
class SearchData:
    nodes: int = 0
    tours_collected: list[tuple[tuple[int, int], ...]] = None

def search(
    state: State,
    ctx: dict[str, object],
    data: SearchData,
    canonical_seen: set[tuple[tuple[int, int], ...]],
    depth: int,
    node_limit: int,
) -> bool:
    if data.nodes >= node_limit:
        return True

    data.nodes += 1
    e = choose_next_edge(state, ctx, depth)
    if e < 0:
        if is_complete_tour(state, ctx):
            key = canonical_translation_key(state, ctx)
            if key not in canonical_seen:
                canonical_seen.add(key)
                data.tours_collected.append(key)
        return False

    for val in (1, 0):
        snap = state.snapshot()
        status = process_queue(state, ctx, [(e, val)])
        if status == OK:
            stop = search(state, ctx, data, canonical_seen, depth + 1, node_limit)
        elif status == COMPLETE_TOUR and is_complete_tour(state, ctx):
            key = canonical_translation_key(state, ctx)
            if key not in canonical_seen:
                canonical_seen.add(key)
                data.tours_collected.append(key)
            stop = data.nodes >= node_limit
        else:
            stop = False
        state.restore(snap)
        if stop:
            return True

    return False

def main():
    n, m, s = 6, 6, 3
    node_limit = 100000 # we will search up to 100,000 nodes
    
    print("========================================================================")
    print("ANÁLISE DE SIMETRIA DE TRANSLACÃO NO TORO CISALHADO 6x6, s=3")
    print("========================================================================")
    print(f"Executando amostragem de busca até {node_limit:,} nós para coletar tours...")
    
    ctx = build_graph(n, m, s)
    state = State(int(ctx["V"]), int(ctx["E"]))
    data = SearchData(tours_collected=[])
    canonical_seen = set()
    
    t0 = time.perf_counter()
    search(state, ctx, data, canonical_seen, 0, node_limit)
    elapsed = time.perf_counter() - t0
    
    total_found = len(data.tours_collected)
    print(f"Busca concluída em {elapsed:.2f}s.")
    print(f"Total de tours canônicos coletados: {total_found:,}")
    
    if total_found == 0:
        print("Nenhum tour encontrado na amostragem.")
        return
        
    print("\nAnalisando tamanhos de órbitas de translação...")
    
    orbit_sizes = {}
    symmetric_tours = []
    
    for tour in data.tours_collected:
        # Compute size of the orbit under translations
        translations = set()
        for dy in range(n):
            for dx in range(m):
                trans_key = translate_edge_key(tour, n, m, dy, dx)
                translations.add(trans_key)
        
        sz = len(translations)
        orbit_sizes[sz] = orbit_sizes.get(sz, 0) + 1
        if sz < 36:
            symmetric_tours.append((sz, tour))
            
    print("\nResultados da Órbita de Translação:")
    print("------------------------------------------------------------------------")
    print(f"  Tamanho da Órbita |  Qtd de Tours  |  Proporção  | Significado")
    print("--------------------+----------------+-------------+--------------------")
    for sz in sorted(orbit_sizes.keys()):
        count = orbit_sizes[sz]
        prop = count / total_found * 100
        signif = "Assimétrico (Típico)" if sz == 36 else f"Simetria Translacional (Período {36 // sz})"
        print(f"  {sz:>16d} |  {count:>14d} |  {prop:>9.2f}% | {signif}")
    print("------------------------------------------------------------------------")
    
    if symmetric_tours:
        print(f"\n✔ SUCESSO: Foram encontrados {len(symmetric_tours)} tours com simetria translacional!")
        print("Exemplo de um tour simétrico (órbita menor):")
        sz, tour = min(symmetric_tours, key=lambda x: x[0])
        print(f"  Tamanho da órbita: {sz} (Invariante sob translação de ordem {36 // sz})")
        print("  Arestas do tour:")
        for u, v in tour[:10]:
            uy, ux = coords(u, m)
            vy, vx = coords(v, m)
            print(f"    ({uy},{ux}) - ({vy},{vx})")
        print("    ...")
    else:
        print("\nℹ Nenhum tour com simetria translacional estrita foi encontrado na amostragem de 100k nós.")
        print("Todos os tours encontrados pertencem a órbitas completas de tamanho 36 (assimétricos).")

if __name__ == "__main__":
    main()
