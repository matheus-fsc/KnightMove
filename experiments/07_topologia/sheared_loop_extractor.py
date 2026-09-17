"""Extrai ciclos fundamentais do grafo do cavalo no toro cisalhado.

O grafo usa a colagem (x, y) ~ (x + m, y + s mod n). Uma BFS a partir de
(0, 0) define a arvore geradora; cada aresta fora da arvore fecha um ciclo
fundamental, que e classificado pelo winding total (W_x, W_y).
"""

from __future__ import annotations

import argparse
from collections import Counter, deque
from dataclasses import dataclass


Vertex = tuple[int, int]  # (y, x)
Edge = tuple[Vertex, Vertex]

KNIGHT_MOVES: tuple[Vertex, ...] = (
    (-2, -1), (-2, 1), (-1, -2), (-1, 2),
    (1, -2), (1, 2), (2, -1), (2, 1),
)


@dataclass(frozen=True)
class FundamentalCycle:
    back_edge: Edge
    vertices: list[Vertex]
    winding: tuple[int, int]

    @property
    def size(self) -> int:
        return len(self.vertices) - 1


def canonical_edge(u: Vertex, v: Vertex) -> Edge:
    return (u, v) if u <= v else (v, u)


def sheared_step(y: int, x: int, dy: int, dx: int, n: int, m: int, s: int) -> Vertex:
    """Aplica um movimento do cavalo no toro cisalhado."""
    wrap_x = (x + dx) // m
    nx = (x + dx) % m
    ny = (y + dy + wrap_x * s) % n
    return ny, nx


def build_sheared_graph(n: int, m: int, s: int) -> dict[Vertex, list[Vertex]]:
    """Gera a adjacencia simples do grafo do cavalo no toro cisalhado."""
    adjacency: dict[Vertex, set[Vertex]] = {
        (y, x): set()
        for y in range(n)
        for x in range(m)
    }

    for y in range(n):
        for x in range(m):
            u = (y, x)
            for dy, dx in KNIGHT_MOVES:
                v = sheared_step(y, x, dy, dx, n, m, s)
                if u == v:
                    continue
                adjacency[u].add(v)
                adjacency[v].add(u)

    return {u: sorted(vs) for u, vs in sorted(adjacency.items())}


def graph_edges(adjacency: dict[Vertex, list[Vertex]]) -> set[Edge]:
    edges = set()
    for u, neighbors in adjacency.items():
        for v in neighbors:
            edges.add(canonical_edge(u, v))
    return edges


def bfs_spanning_tree(
    adjacency: dict[Vertex, list[Vertex]],
    root: Vertex = (0, 0),
) -> tuple[dict[Vertex, Vertex | None], dict[Vertex, int], set[Edge], list[Edge]]:
    """Constroi a BFS e devolve parent, depth, arestas da arvore e back-edges."""
    parent: dict[Vertex, Vertex | None] = {root: None}
    depth: dict[Vertex, int] = {root: 0}
    tree_edges: set[Edge] = set()
    back_edges: list[Edge] = []
    seen_back_edges: set[Edge] = set()

    queue: deque[Vertex] = deque([root])
    while queue:
        u = queue.popleft()
        for v in adjacency[u]:
            edge = canonical_edge(u, v)
            if v not in parent:
                parent[v] = u
                depth[v] = depth[u] + 1
                tree_edges.add(edge)
                queue.append(v)
            elif v != parent[u] and edge not in tree_edges and edge not in seen_back_edges:
                seen_back_edges.add(edge)
                back_edges.append(edge)

    if len(parent) != len(adjacency):
        missing = len(adjacency) - len(parent)
        raise RuntimeError(f"grafo desconexo: {missing} vertices nao alcancados pela BFS")

    return parent, depth, tree_edges, back_edges


def reconstruct_cycle(u: Vertex, v: Vertex, parent: dict[Vertex, Vertex | None]) -> list[Vertex]:
    """Reconstrui LCA -> ... -> u -> v -> ... -> LCA para a back-edge (u, v)."""
    path_u: list[Vertex] = []
    cur: Vertex | None = u
    while cur is not None:
        path_u.append(cur)
        cur = parent[cur]

    ancestors_u = set(path_u)
    path_v: list[Vertex] = []
    cur = v
    while cur not in ancestors_u:
        path_v.append(cur)
        cur = parent[cur]
        if cur is None:
            raise RuntimeError("falha ao encontrar LCA; parent tree inconsistente")

    lca = cur
    path_v.append(lca)
    path_u_to_lca = path_u[:path_u.index(lca) + 1]

    return list(reversed(path_u_to_lca)) + path_v


def decode_delta_x(dx_apparent: int, m: int) -> int:
    candidates = [dx for dx in (-2, -1, 1, 2) if dx % m == dx_apparent % m]
    if len(candidates) != 1:
        raise ValueError(f"dx aparente {dx_apparent} tem candidatos {candidates}")
    return candidates[0]


def decode_delta_y(dy_base: int, n: int) -> int:
    candidates = [dy for dy in (-2, -1, 1, 2) if dy % n == dy_base % n]
    if len(candidates) != 1:
        raise ValueError(f"dy base {dy_base} tem candidatos {candidates}")
    return candidates[0]


def decode_winding(cycle: list[Vertex], n: int, m: int, s: int) -> tuple[int, int]:
    """Calcula (W_x, W_y) na base topologica (m, -s), (0, n)."""
    sum_dx = 0
    sum_dy = 0

    for (uy, ux), (vy, vx) in zip(cycle, cycle[1:]):
        dx_apparent = (vx - ux) % m
        dx = decode_delta_x(dx_apparent, m)
        wrap_x = (ux + dx) // m

        dy_base = (vy - uy - wrap_x * s) % n
        dy = decode_delta_y(dy_base, n)

        sum_dx += dx
        sum_dy += dy

    if sum_dx % m != 0:
        raise RuntimeError(f"ciclo invalido: soma dx={sum_dx} nao e multiplo de m={m}")

    w_x = sum_dx // m
    adjusted_y = sum_dy + s * w_x
    if adjusted_y % n != 0:
        raise RuntimeError(
            "ciclo invalido: soma dy + s*W_x="
            f"{adjusted_y} nao e multiplo de n={n}"
        )

    return w_x, adjusted_y // n


def extract_fundamental_cycles(n: int, m: int, s: int) -> tuple[dict[str, int], list[FundamentalCycle]]:
    adjacency = build_sheared_graph(n, m, s)
    edges = graph_edges(adjacency)
    parent, _depth, tree_edges, back_edges = bfs_spanning_tree(adjacency)

    beta_1 = len(edges) - len(adjacency) + 1
    if len(back_edges) != beta_1:
        raise RuntimeError(
            f"contagem inconsistente: back_edges={len(back_edges)} beta_1={beta_1}"
        )

    cycles = []
    for u, v in back_edges:
        vertices = reconstruct_cycle(u, v, parent)
        winding = decode_winding(vertices, n, m, s)
        cycles.append(FundamentalCycle((u, v), vertices, winding))

    stats = {
        "vertices": len(adjacency),
        "edges": len(edges),
        "tree_edges": len(tree_edges),
        "beta_1": beta_1,
    }
    return stats, cycles


def format_vertex_sequence(vertices: list[Vertex]) -> str:
    return " -> ".join(f"({x},{y})" for y, x in vertices)


def format_edge(edge: Edge) -> str:
    return f"({format_vertex_sequence([edge[0]])}, {format_vertex_sequence([edge[1]])})"


def print_collapse_analysis(
    n: int,
    m: int,
    s: int,
    classes: Counter[tuple[int, int]],
) -> None:
    """Imprime a verificacao empirica da paridade vertical para s=n/2."""
    even_y_count = sum(count for (_wx, wy), count in classes.items() if wy % 2 == 0)
    odd_y_count = sum(count for (_wx, wy), count in classes.items() if wy % 2 != 0)
    even_y_classes = {winding: count for winding, count in sorted(classes.items()) if winding[1] % 2 == 0}
    odd_y_classes = {winding: count for winding, count in sorted(classes.items()) if winding[1] % 2 != 0}

    print("\nAnalise de colapso topologico:")
    if s * 2 == n:
        print(f"  Caso critico detectado: s={s}=n/2 em {n}x{m}.")
    else:
        print(f"  Caso nao critico: s={s}; a regra especial s=n/2 nao se aplica diretamente.")

    print(f"  Ciclos globais com W_y par   : {even_y_count}")
    print(f"  Ciclos globais com W_y impar : {odd_y_count}")

    if even_y_classes:
        print("  Classes sobreviventes candidatas (W_y par):")
        for winding, count in even_y_classes.items():
            print(f"    W={winding}: {count}")

    if odd_y_classes:
        print("  Classes que deveriam estar extintas (W_y impar), mas apareceram:")
        for winding, count in odd_y_classes.items():
            print(f"    W={winding}: {count}")
    else:
        print("  Classes extintas (W_y impar): exatamente 0 ciclos.")

    if s * 2 == n and odd_y_count == 0:
        print(
            "  Conclusao: CONFIRMADA na base de ciclos fundamentais; "
            "a paridade vertical impar colapsou."
        )
    elif s * 2 == n:
        print(
            "  Conclusao: NAO CONFIRMADA na base de ciclos fundamentais; "
            "existem ciclos fundamentais com W_y impar."
        )


def print_s1_comparison(n: int, m: int, s: int, classes: Counter[tuple[int, int]]) -> None:
    if not (n == 6 and m == 6 and s == 3):
        return

    s1_stats, s1_cycles = extract_fundamental_cycles(n, m, 1)
    s1_global = [cycle for cycle in s1_cycles if cycle.winding != (0, 0)]
    s1_classes = Counter(cycle.winding for cycle in s1_global)
    odd_y_count = sum(count for (_wx, wy), count in classes.items() if wy % 2 != 0)

    print("\nComparacao com o controle s=1:")
    print(
        f"  No controle s=1 foram observadas {len(s1_classes)} classes globais "
        f"em {len(s1_global)} ciclos globais (beta1={s1_stats['beta_1']})."
    )
    print(
        f"  Em s=3 foram observadas {len(classes)} classes globais, "
        f"mas {odd_y_count} ciclos globais possuem W_y impar."
    )
    if odd_y_count == 0:
        print(
            "  Declaracao: em relacao ao s=1, o colapso de paridade foi confirmado "
            "para esta base fundamental."
        )
    else:
        print(
            "  Declaracao: em relacao ao s=1, nao houve colapso de paridade na "
            "base de ciclos fundamentais; a restricao pode pertencer a uma "
            "subfamilia adicional, como tours Hamiltonianos, e nao ao H1 gerado "
            "pela BFS."
        )


def print_report(
    n: int,
    m: int,
    s: int,
    stats: dict[str, int],
    cycles: list[FundamentalCycle],
    show_vertices: bool,
    max_global: int | None,
) -> None:
    local_cycles = [cycle for cycle in cycles if cycle.winding == (0, 0)]
    global_cycles = [cycle for cycle in cycles if cycle.winding != (0, 0)]
    classes = Counter(cycle.winding for cycle in global_cycles)

    print("=" * 72)
    print(f"Extracao de ciclos fundamentais - toro cisalhado {n}x{m}, s={s}")
    print("=" * 72)
    print(f"Vertices (V)                         : {stats['vertices']}")
    print(f"Arestas (E)                          : {stats['edges']}")
    print(f"Arestas da arvore BFS                : {stats['tree_edges']}")
    print(f"Arestas de retorno encontradas beta1 : {stats['beta_1']}")
    print(f"Ciclos locais (W_x,W_y)=(0,0)        : {len(local_cycles)}")
    print(f"Ciclos globais                       : {len(global_cycles)}")

    if classes:
        print("\nClasses topologicas globais:")
        for winding, count in sorted(classes.items()):
            print(f"  W={winding}: {count}")

    print_collapse_analysis(n, m, s, classes)
    print_s1_comparison(n, m, s, classes)

    print("\nCiclos globais encontrados:")
    if not global_cycles:
        print("  nenhum")
        return

    listed = global_cycles if max_global is None else global_cycles[:max_global]
    for index, cycle in enumerate(listed, start=1):
        print(
            f"  #{index:03d} W={cycle.winding} tamanho={cycle.size} "
            f"back-edge={format_edge(cycle.back_edge)}"
        )
        if show_vertices:
            print(f"       {format_vertex_sequence(cycle.vertices)}")

    omitted = len(global_cycles) - len(listed)
    if omitted > 0:
        print(f"  ... {omitted} ciclos globais omitidos; use --max-global para ajustar.")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extrai e classifica ciclos fundamentais no toro cisalhado."
    )
    parser.add_argument("--n", type=int, default=6, help="altura / eixo Y modular")
    parser.add_argument("--m", type=int, default=6, help="largura / eixo X modular")
    parser.add_argument("--s", type=int, default=3, help="parametro de cisalhamento")
    parser.add_argument(
        "--show-vertices",
        action="store_true",
        help="imprime a sequencia fechada de vertices de cada ciclo global",
    )
    parser.add_argument(
        "--max-global",
        type=int,
        default=None,
        help="limita a quantidade de ciclos globais detalhados",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    stats, cycles = extract_fundamental_cycles(args.n, args.m, args.s)
    print_report(args.n, args.m, args.s, stats, cycles, args.show_vertices, args.max_global)


if __name__ == "__main__":
    main()
