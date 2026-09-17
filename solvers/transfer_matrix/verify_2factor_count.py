"""
Verificacao independente: conta 2-fatores do grafo do cavalo n x n
via backtracking baseado em ARESTAS (em vez de coluna),
com pruning de grau e de viabilidade.

Compara com a contagem de transfer_build.py.

n=4 deve rodar em segundos. n=6 em ~minuto se pruning for bom.
"""

import sys
import time

DRC = [(1, 2), (2, 1), (-1, 2), (-2, 1), (1, -2), (2, -1), (-1, -2), (-2, -1)]


def build_graph(n):
    V = [(r, c) for r in range(n) for c in range(n)]
    pos = {v: i for i, v in enumerate(V)}
    edges = []
    seen = set()
    inc = [[] for _ in range(len(V))]  # edges incidentes a cada vertice (idx)
    for v in V:
        r, c = v
        for dr, dc in DRC:
            u = (r + dr, c + dc)
            if 0 <= u[0] < n and 0 <= u[1] < n:
                e = tuple(sorted([v, u]))
                if e not in seen:
                    seen.add(e)
                    eid = len(edges)
                    edges.append(e)
                    inc[pos[v]].append(eid)
                    inc[pos[u]].append(eid)
    return V, edges, inc, pos


def count_2factors(n, verbose=False):
    V, edges, inc, pos = build_graph(n)
    nV = len(V)
    nE = len(edges)
    if verbose:
        print(f"  V={nV}, E={nE}")

    # Para cada vertice: lista de edges incidentes
    # Para podar: precisamos que cada vertice acabe com grau 2.
    # Iteramos pelos vertices em ordem; para cada um, escolhemos quais arestas
    # incidentes "comprometemos" a usar -- usando apenas arestas (u, v) com
    # u <= v na ordem (i.e., u tambem ja processado ou e' o atual).

    # Mas como uma aresta envolve dois vertices, e' mais simples:
    # iterar sobre arestas em ordem. Para cada aresta decidir 'in' ou 'out'.
    # Pruning: para cada vertice, ja_in + ainda_disponiveis >= 2 - ja_in. ok.

    deg = [0] * nV
    # remaining[v] = numero de arestas incidentes a v ainda nao decididas
    remaining = [len(x) for x in inc]
    count = [0]

    # ordenar arestas: por max-grau dos endpoints (vertices restritos primeiro)
    # ajuda pruning.
    edge_order = sorted(range(nE), key=lambda eid: -(len(inc[pos[edges[eid][0]]]) + len(inc[pos[edges[eid][1]]])))
    # ou simplesmente em ordem original
    edge_order = list(range(nE))

    def rec(i):
        if i == nE:
            if all(d == 2 for d in deg):
                count[0] += 1
            return

        eid = edge_order[i]
        u, v = edges[eid]
        ui = pos[u]
        vi = pos[v]

        # OPCAO 1: descartar aresta
        remaining[ui] -= 1
        remaining[vi] -= 1
        # pruning: cada vertice ainda precisa de (2 - deg) >= 0 e <= remaining
        if (2 - deg[ui] <= remaining[ui]) and (2 - deg[vi] <= remaining[vi]):
            rec(i + 1)
        # remaining[ui] += 1  # vamos restaurar depois
        # remaining[vi] += 1

        # OPCAO 2: incluir aresta
        if deg[ui] < 2 and deg[vi] < 2:
            deg[ui] += 1
            deg[vi] += 1
            # remaining ja foi decrementado acima
            if (2 - deg[ui] <= remaining[ui]) and (2 - deg[vi] <= remaining[vi]):
                rec(i + 1)
            deg[ui] -= 1
            deg[vi] -= 1

        # restaurar remaining
        remaining[ui] += 1
        remaining[vi] += 1

    rec(0)
    return count[0]


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 4

    print(f"n = {n}:")
    t0 = time.time()
    c = count_2factors(n, verbose=True)
    t1 = time.time()
    print(f"  2-fatores: {c}")
    print(f"  tempo: {t1 - t0:.2f}s")
