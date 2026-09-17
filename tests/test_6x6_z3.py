#!/usr/bin/env python3
"""
cavalo_engine_6x6_allsat.py
Extrai múltiplas soluções (All-SAT) do tabuleiro 6x6 para estudo de assinaturas.
"""
import json
import os
import sys

# MUDANÇA 1: Tabuleiro 6x6
BOARD = 6
START_NODE = (0, 0)
# Por padrão, extrai todas as soluções possíveis.
# Para limitar a execução, use a variável de ambiente MAX_SOLUTIONS (ex.: MAX_SOLUTIONS=1000).
MAX_SOLUTIONS = int(os.getenv("MAX_SOLUTIONS")) if os.getenv("MAX_SOLUTIONS") else None

# --- GRUPO DIÉDRICO D4 (Transformações Geométricas) ---
def rot90(x, y): return (y, BOARD - 1 - x)
def rot180(x, y): return (BOARD - 1 - x, BOARD - 1 - y)
def rot270(x, y): return (BOARD - 1 - y, x)
def refX(x, y): return (BOARD - 1 - x, y)
def refY(x, y): return (x, BOARD - 1 - y)
def refD1(x, y): return (y, x)
def refD2(x, y): return (BOARD - 1 - y, BOARD - 1 - x)

TRANSFORMACOES = [rot90, rot180, rot270, refX, refY, refD1, refD2]

def transpor_aresta(u, v, transformacao):
    """Aplica a rotação/espelhamento aos vértices de uma aresta."""
    nu = transformacao(*u)
    nv = transformacao(*v)
    return tuple(sorted([nu, nv]))

KNIGHT_DELTAS = [(1, 2), (2, 1), (-1, 2), (-2, 1),
                 (1, -2), (2, -1), (-1, -2), (-2, -1)]

def neighbors(x, y):
    out = []
    for dx, dy in KNIGHT_DELTAS:
        nx, ny = x + dx, y + dy
        if 0 <= nx < BOARD and 0 <= ny < BOARD:
            out.append((nx, ny))
    def degree(n):
        return sum(1 for dx, dy in KNIGHT_DELTAS
                   if 0 <= n[0]+dx < BOARD and 0 <= n[1]+dy < BOARD)
    return sorted(out, key=degree)

def build_graph():
    nodes = [(x, y) for x in range(BOARD) for y in range(BOARD)]
    adj = {n: neighbors(*n) for n in nodes}
    edges = set()
    for n, nbrs in adj.items():
        for m in nbrs: edges.add(tuple(sorted([n, m])))
    return nodes, adj, edges

def compute_spanning_tree(adj, start=START_NODE):
    parent = {start: None}; level = {start: 0}; order = []; stack = [start]
    tree_edges = set(); visited = {start}
    while stack:
        u = stack.pop()
        order.append(u)
        for v in reversed(adj[u]):
            if v not in visited:
                visited.add(v); parent[v] = u; level[v] = level[u] + 1
                tree_edges.add(tuple(sorted([u, v]))); stack.append(v)
    return parent, level, order, tree_edges

def extract_loops(parent, tree_edges, edges):
    loops = []
    for e in sorted(edges):
        if e in tree_edges: continue
        u, v = e
        path_u = []; curr = u
        while curr is not None: path_u.append(curr); curr = parent[curr]
        path_v = []; curr = v
        while curr is not None: path_v.append(curr); curr = parent[curr]
        loops.append({
            "caminho1": [list(p) for p in path_u],
            "caminho2": [list(p) for p in path_v],
            "colisao": [list(u), list(v)],
            "tamanho_ciclo": len(path_u) + len(path_v),
        })
    return loops

def solve_allsat_z3(nodes, edges, loops):
    from z3 import Solver, Bool, Sum, If, sat, is_true, Or
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
        for e in get_loop_edge_set(lp): edge_to_loops[e].append(i)

    edge_active = {}
    for e in edges:
        if not edge_to_loops[e]: edge_active[e] = False
        else: edge_active[e] = (Sum([If(loop_vars[i], 1, 0) for i in edge_to_loops[e]]) % 2 == 1)

    for n in nodes:
        incident = [e for e in edges if n in e]
        solver.add(Sum([If(edge_active[e], 1, 0) for e in incident]) == 2)

    # ---------------------------------------------------------
    # SYMMETRY BREAKING PREDICATES (SBP) - Quebra do Grupo D4
    # ---------------------------------------------------------
    print("Injetando Predicados de Quebra de Simetria (D4)...")

    arestas_ordenadas = sorted(list(edges))

    def lex_leq(z3_vars_A, z3_vars_B):
        condicao = True
        for a, b in reversed(list(zip(z3_vars_A, z3_vars_B))):
            condicao = If(a == b, condicao, b)
        return condicao

    vetor_original = [edge_active[e] for e in arestas_ordenadas]

    for transformacao in TRANSFORMACOES:
        vetor_transformado = []
        for e in arestas_ordenadas:
            aresta_girada = transpor_aresta(e[0], e[1], transformacao)
            vetor_transformado.append(edge_active[aresta_girada])
        solver.add(lex_leq(vetor_original, vetor_transformado))
    # ---------------------------------------------------------

    todas_assinaturas = []

    if MAX_SOLUTIONS is None:
        print("Buscando todas as soluções possíveis...")
    else:
        print(f"Buscando as primeiras {MAX_SOLUTIONS} soluções...")

    while True:
        if solver.check() == sat:
            m = solver.model()
            active = [e for e in edges if edge_to_loops[e] and is_true(m.evaluate(edge_active[e]))]

            adj_res = {n: [] for n in nodes}
            for u, v in active: adj_res[u].append(v); adj_res[v].append(u)

            visited, components = set(), []
            for n in nodes:
                if n not in visited:
                    comp, q = [], [n]; visited.add(n)
                    while q:
                        curr = q.pop(0); comp.append(curr)
                        for nbr in adj_res[curr]:
                            if nbr not in visited: visited.add(nbr); q.append(nbr)
                    components.append(comp)

            if len(components) == 1:
                # MUDANÇA 2: Salva a assinatura e aplica a Bomba de Bloqueio
                assinatura = [is_true(m.evaluate(v)) for v in loop_vars]
                todas_assinaturas.append(assinatura)

                print(f"-> Solução #{len(todas_assinaturas)} extraída!")

                if MAX_SOLUTIONS is not None and len(todas_assinaturas) >= MAX_SOLUTIONS:
                    print("\nLimite de amostras atingido. Finalizando extração.")
                    break

                # CLAÚSULA DE BLOQUEIO (Impede o Z3 de repetir a assinatura)
                clausula_bloqueio = []
                for i in range(len(loops)):
                    if assinatura[i]: clausula_bloqueio.append(loop_vars[i] == False)
                    else: clausula_bloqueio.append(loop_vars[i] == True)
                solver.add(Or(clausula_bloqueio))
                continue

            for comp in components:
                comp_set = set(comp)
                in_edges = [e for e in edges if e[0] in comp_set and e[1] in comp_set]
                solver.add(Sum([If(edge_active[e], 1, 0) for e in in_edges]) <= len(comp) - 1)
        else:
            print("\nEspaço de busca esgotado. Nenhuma outra solução existe.")
            break

    return todas_assinaturas

def main():
    nodes, adj, edges = build_graph()
    parent, level, order, tree_edges = compute_spanning_tree(adj)
    loops = extract_loops(parent, tree_edges, edges)

    assinaturas = solve_allsat_z3(nodes, edges, loops)

    if assinaturas:
        data = {
            "descricao": {"V": len(nodes), "E": len(edges), "H1": len(loops), "amostras": len(assinaturas)},
            "assinaturas_booleanas": assinaturas
        }
        with open("assinaturas_6x6.json", "w") as f: json.dump(data, f)
        print(f"\nSucesso! {len(assinaturas)} assinaturas exportadas para assinaturas_6x6.json.")

if __name__ == "__main__":
    main()
