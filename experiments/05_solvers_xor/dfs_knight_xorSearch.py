#!/usr/bin/env python3
import json
import sys
from collections import deque

# Configurações do Universo
BOARD = 8
START_NODE = (0, 0)

KNIGHT_DELTAS = [(1, 2), (2, 1), (-1, 2), (-2, 1),
                 (1, -2), (2, -1), (-1, -2), (-2, -1)]

def neighbors(x, y):
    out = []
    for dx, dy in KNIGHT_DELTAS:
        nx, ny = x + dx, y + dy
        if 0 <= nx < BOARD and 0 <= ny < BOARD:
            out.append((nx, ny))
    # Heurística de Warnsdorff embutida para guiar a DFS
    def degree(n):
        return sum(1 for dx, dy in KNIGHT_DELTAS
                   if 0 <= n[0]+dx < BOARD and 0 <= n[1]+dy < BOARD)
    return sorted(out, key=degree)

def build_graph():
    nodes = [(x, y) for x in range(BOARD) for y in range(BOARD)]
    adj = {n: neighbors(*n) for n in nodes}
    edges = set()
    for n, nbrs in adj.items():
        for m in nbrs:
            edges.add(tuple(sorted([n, m])))
    return nodes, adj, edges

def compute_spanning_tree(adj, start=START_NODE):
    """Gera a Árvore Cobra (DFS) para minimizar o emaranhado."""
    parent = {start: None}
    level = {start: 0}
    order = []
    stack = [start]
    tree_edges = set()
    visited = {start}

    while stack:
        u = stack.pop()
        order.append(u)
        for v in reversed(adj[u]):
            if v not in visited:
                visited.add(v)
                parent[v] = u
                level[v] = level[u] + 1
                tree_edges.add(tuple(sorted([u, v])))
                stack.append(v)
    return parent, level, order, tree_edges

def extract_loops(parent, tree_edges, edges):
    """Calcula a Base Fundamental H1 via XOR."""
    loops = []
    for e in sorted(edges):
        if e in tree_edges: continue
        u, v = e
        path_u = []
        curr = u
        while curr is not None:
            path_u.append(curr)
            curr = parent[curr]
        path_v = []
        curr = v
        while curr is not None:
            path_v.append(curr)
            curr = parent[curr]

        loops.append({
            "caminho1": [list(p) for p in path_u],
            "caminho2": [list(p) for p in path_v],
            "colisao": [list(u), list(v)],
            "tamanho_ciclo": len(path_u) + len(path_v),
        })
    return loops

def solve_with_z3(nodes, edges, loops):
    from z3 import Solver, Bool, Sum, If, sat, is_true
    solver = Solver()
    loop_vars = [Bool(f"L_{i}") for i in range(len(loops))]

    def get_loop_edge_set(loop):
        edges_set = set()
        def toggle(u, v):
            e = tuple(sorted([tuple(u), tuple(v)]))
            if e in edges_set: edges_set.remove(e)
            else: edges_set.add(e)
        c1, c2 = loop['caminho1'], loop['caminho2']
        for i in range(len(c1)-1): toggle(c1[i], c1[i+1])
        for i in range(len(c2)-1): toggle(c2[i], c2[i+1])
        toggle(loop['colisao'][0], loop['colisao'][1])
        return edges_set

    edge_to_loops = {e: [] for e in edges}
    for i, lp in enumerate(loops):
        for e in get_loop_edge_set(lp):
            edge_to_loops[e].append(i)

    edge_active = {}
    for e in edges:
        if not edge_to_loops[e]: edge_active[e] = False
        else: edge_active[e] = (Sum([If(loop_vars[i], 1, 0) for i in edge_to_loops[e]]) % 2 == 1)

    for n in nodes:
        incident = [e for e in edges if n in e]
        solver.add(Sum([If(edge_active[e], 1, 0) for e in incident]) == 2)

    while True:
        if solver.check() == sat:
            m = solver.model()
            active = [e for e in edges if edge_to_loops[e] and is_true(m.evaluate(edge_active[e]))]

            # Checar Sub-rotas (Ilhas)
            adj_res = {n: [] for n in nodes}
            for u, v in active:
                adj_res[u].append(v); adj_res[v].append(u)

            visited, components = set(), []
            for n in nodes:
                if n not in visited:
                    comp, q = [], [n]
                    visited.add(n)
                    while q:
                        curr = q.pop(0); comp.append(curr)
                        for nbr in adj_res[curr]:
                            if nbr not in visited: visited.add(nbr); q.append(nbr)
                    components.append(comp)

            if len(components) == 1: return active, [is_true(m.evaluate(v)) for v in loop_vars]

            for comp in components:
                comp_set = set(comp)
                in_edges = [e for e in edges if e[0] in comp_set and e[1] in comp_set]
                solver.add(Sum([If(edge_active[e], 1, 0) for e in in_edges]) <= len(comp) - 1)
        else: return None, None

def main():
    nodes, adj, edges = build_graph()
    parent, level, order, tree_edges = compute_spanning_tree(adj)
    loops = extract_loops(parent, tree_edges, edges)

    active_edges, loop_states = solve_with_z3(nodes, edges, loops)

    if active_edges:
        # Ordenação do caminho final
        adj_final = {n: [] for n in nodes}
        for u, v in active_edges:
            adj_final[u].append(v); adj_final[v].append(u)

        path = [START_NODE]
        while len(path) < len(nodes):
            curr = path[-1]
            next_node = [n for n in adj_final[curr] if n not in path][0]
            path.append(next_node)

        data = {
            "descricao": {"V": len(nodes), "E": len(edges), "H1": len(loops), "inicio": list(START_NODE)},
            "loops": loops,
            "solucao": {"caminho": [list(c) for c in path], "loop_activation": loop_states, "completo": True}
        }
        with open("cavalo_data.json", "w") as f: json.dump(data, f)
        print(f"Sucesso! {len(path)} lances exportados.")

if __name__ == "__main__":
    main()
