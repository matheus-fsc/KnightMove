#!/usr/bin/env python3
"""
bishop_queen_6x6.py
===================
Teste de generalização: bispo e dama em 6×6.

Definição (convencional para "piece graph"):
  Aresta entre v1 e v2 sse a peça consegue ir de v1 a v2 em UM lance
  (i.e., toda a diagonal para bispo; linha + coluna + diagonal para dama).

Reporta:
  - Graus de cada vértice da primeira linha
  - Componentes conexos (bispo é bipartido por cor e desconexo)
  - Existência de ciclo Hamiltoniano (DFS com Warnsdorff)
  - Se existe, número de DOF combinatórios na borda (C(deg, 2))
"""

import sys
import time
from collections import defaultdict
from math import comb


BOARD = 6


def vid(r, c): return r * BOARD + c
def vrc(v): return divmod(v, BOARD)
def label(v):
    r, c = vrc(v)
    return chr(65 + c) + str(BOARD - r)


def build_bishop_edges():
    """Bispo: arestas entre quaisquer dois vértices na mesma diagonal."""
    edges = set()
    for v1 in range(BOARD * BOARD):
        r1, c1 = vrc(v1)
        for v2 in range(v1 + 1, BOARD * BOARD):
            r2, c2 = vrc(v2)
            if (r1 - r2) == (c1 - c2) or (r1 - r2) == -(c1 - c2):
                edges.add((v1, v2))
    return edges


def build_queen_edges():
    """Dama: bispo ∪ torre."""
    edges = set()
    for v1 in range(BOARD * BOARD):
        r1, c1 = vrc(v1)
        for v2 in range(v1 + 1, BOARD * BOARD):
            r2, c2 = vrc(v2)
            if (r1 == r2 or c1 == c2 or
                    abs(r1 - r2) == abs(c1 - c2)):
                edges.add((v1, v2))
    return edges


def adjacency(edges):
    adj = defaultdict(set)
    for u, v in edges:
        adj[u].add(v)
        adj[v].add(u)
    return adj


def connected_components(adj, n_verts):
    visited = [False] * n_verts
    comps = []
    for start in range(n_verts):
        if visited[start]:
            continue
        comp = []
        stack = [start]
        visited[start] = True
        while stack:
            u = stack.pop()
            comp.append(u)
            for v in adj[u]:
                if not visited[v]:
                    visited[v] = True
                    stack.append(v)
        comps.append(sorted(comp))
    return comps


def is_bipartite(adj, comp):
    """Testa se um componente é bipartido (BFS 2-color)."""
    color = {}
    start = comp[0]
    color[start] = 0
    stack = [start]
    while stack:
        u = stack.pop()
        for v in adj[u]:
            if v not in color:
                color[v] = 1 - color[u]
                stack.append(v)
            elif color[v] == color[u]:
                return False, None
    sides = ([v for v in comp if color.get(v) == 0],
             [v for v in comp if color.get(v) == 1])
    return True, sides


def find_ham_cycle(adj, vertices, time_limit=10.0):
    """DFS com Warnsdorff. Procura um ciclo Hamiltoniano dentro dos
    'vertices' dados (subconjunto). Retorna lista ou None."""
    V = set(vertices)
    n = len(V)
    start = vertices[0]
    path = [start]
    visited = {start}
    t0 = time.time()
    timeout = [False]

    def dfs():
        if timeout[0]:
            return None
        if time.time() - t0 > time_limit:
            timeout[0] = True
            return None
        if len(path) == n:
            if path[0] in adj[path[-1]]:
                return list(path)
            return None
        last = path[-1]
        candidates = adj[last] & V - visited
        nbrs = sorted(candidates,
                      key=lambda v: len((adj[v] & V) - visited))
        for nb in nbrs:
            path.append(nb)
            visited.add(nb)
            res = dfs()
            if res is not None:
                return res
            path.pop()
            visited.remove(nb)
        return None

    res = dfs()
    return res, timeout[0]


def report_piece(name, edges):
    print()
    print("=" * 72)
    print(f"GRAFO DA {name.upper()} {BOARD}×{BOARD}")
    print("=" * 72)
    n_verts = BOARD * BOARD
    adj = adjacency(edges)
    degs = [len(adj[v]) for v in range(n_verts)]
    print(f"\n  Vértices: {n_verts}    Arestas: {len(edges)}    "
          f"Grau min/máx/média: {min(degs)}/{max(degs)}/"
          f"{sum(degs)/n_verts:.2f}")

    # graus da primeira linha
    print(f"\n  Graus da primeira linha (top row):")
    for c in range(BOARD):
        v = vid(0, c)
        print(f"    {label(v)} = (0,{c}): grau {len(adj[v])}")

    # graus da segunda linha (B5..E5 são "segunda posição na borda" lateral)
    print(f"\n  Graus da segunda linha (row 1):")
    for c in range(BOARD):
        v = vid(1, c)
        print(f"    {label(v)} = (1,{c}): grau {len(adj[v])}")

    # componentes
    comps = connected_components(adj, n_verts)
    print(f"\n  Componentes conexos: {len(comps)}")
    for i, c in enumerate(comps):
        bip, sides = is_bipartite(adj, c) if c else (True, ([], []))
        bip_info = ""
        if bip and sides:
            bip_info = f"  bipartido |sides|={len(sides[0])},{len(sides[1])}"
        elif not bip:
            bip_info = f"  NÃO bipartido"
        print(f"    comp[{i}]: {len(c)} vértices{bip_info}")
        if len(c) <= 20:
            labs = [label(v) for v in c]
            print(f"      {labs}")

    # ciclo Hamiltoniano em cada componente
    print(f"\n  Busca por ciclo Hamiltoniano:")
    for i, c in enumerate(comps):
        if len(c) < 3:
            print(f"    comp[{i}] (n={len(c)}): trivial — não tem ciclo "
                  f"Hamiltoniano")
            continue
        cycle, to = find_ham_cycle(adj, c, time_limit=10.0)
        if cycle:
            print(f"    comp[{i}] (n={len(c)}): SIM — exemplo: "
                  f"{[label(v) for v in cycle[:8]]}...")
        else:
            ttag = " (timeout)" if to else ""
            print(f"    comp[{i}] (n={len(c)}): NÃO encontrado em 10s{ttag}")

    return adj, degs, comps


def main():
    # --- BISPO ---
    bishop_edges = build_bishop_edges()
    bishop_adj, bishop_degs, bishop_comps = report_piece("Bispo", bishop_edges)

    # --- DAMA ---
    queen_edges = build_queen_edges()
    queen_adj, queen_degs, queen_comps = report_piece("Dama", queen_edges)

    # --- comparação DOF na borda vs cavalo ---
    print()
    print("=" * 72)
    print("COMPARAÇÃO DE DOF NA BORDA — CAVALO vs BISPO vs DAMA (6×6)")
    print("=" * 72)
    print()
    print(f"  DOF combinatório = C(grau, 2) — pares possíveis de arestas")
    print(f"  incidentes que o ciclo Hamiltoniano pode escolher num vértice.")
    print()
    print(f"  {'vértice':<10s}  {'cavalo':>12s}  {'bispo':>12s}  {'dama':>12s}")
    print(f"  {'-'*10}  {'-'*12}  {'-'*12}  {'-'*12}")

    # graus do cavalo na primeira linha (hardcoded, conhecidos)
    cav_deg_top = {0: 2, 1: 3, 2: 4, 3: 4, 4: 3, 5: 2}
    for c in range(BOARD):
        v = vid(0, c)
        lab = label(v)
        dc = cav_deg_top[c]
        db = len(bishop_adj[v])
        dq = len(queen_adj[v])
        fc = f"deg={dc:>2}, C={comb(dc,2):>3}"
        fb = f"deg={db:>2}, C={comb(db,2):>3}"
        fq = f"deg={dq:>2}, C={comb(dq,2):>3}"
        print(f"  {lab:<10s}  {fc:>12s}  {fb:>12s}  {fq:>12s}")

    print()
    print(f"  Total DOF na borda (Σ C(deg,2) primeira linha):")
    sum_cav = sum(comb(d, 2) for d in cav_deg_top.values())
    sum_bishop = sum(comb(len(bishop_adj[vid(0,c)]), 2) for c in range(BOARD))
    sum_queen = sum(comb(len(queen_adj[vid(0,c)]), 2) for c in range(BOARD))
    print(f"    cavalo : {sum_cav}")
    print(f"    bispo  : {sum_bishop}  ({sum_bishop/sum_cav:.1f}× cavalo)")
    print(f"    dama   : {sum_queen}  ({sum_queen/sum_cav:.1f}× cavalo)")

    return 0


if __name__ == "__main__":
    sys.exit(main())
