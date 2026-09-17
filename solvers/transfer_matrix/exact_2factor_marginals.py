"""
Marginais EXATAS de 2-fatores via DP finito (gold standard).

Conta, para cada aresta e, o numero de 2-fatores que contem e:
  N(e) = sum_F 1{e in F}
e o marginal:
  P(e) = N(e) / N(total)

Implementacao: DP coluna-por-coluna mantendo um v[state] = vetor de
counts, e em paralelo um e_count[edge] = numero de 2-fatores que passam
por aquela aresta.

Estrategia: para cada passo c, ao enumerar transicoes, separamos a
contribuicao com aresta e versus sem.
"""

import json
import sys
import time
from collections import defaultdict
from itertools import combinations
from pathlib import Path

import numpy as np

THIS = Path(__file__).parent
DATA = THIS / "data"

sys.path.insert(0, str(THIS))
from transfer_build import forward_options


def transitions_with_edges(state, n, forbid_cp2=False):
    """Mesma logica que compare_finf.transitions_with_edges -- duplicada para
    independencia."""
    deg_c, deg_cp1 = state
    row_choices = []
    for r in range(n):
        need = 2 - deg_c[r]
        if need < 0:
            return
        opts = forward_options(n, r)
        if forbid_cp2:
            opts = [o for o in opts if o[0] == 1]
        if need > len(opts):
            return
        row_choices.append((need, opts))

    def rec(r, cp1_inc, cp2_inc, partial_edges):
        if r == n:
            new_cp1 = tuple(deg_cp1[i] + cp1_inc[i] for i in range(n))
            new_cp2 = tuple(cp2_inc[i] for i in range(n))
            yield (new_cp1, new_cp2), tuple(partial_edges)
            return

        need, opts = row_choices[r]
        if need == 0:
            yield from rec(r + 1, cp1_inc, cp2_inc, partial_edges)
            return

        for combo in combinations(range(len(opts)), need):
            new_edges = []
            valid = True
            cp1_back, cp2_back = [], []
            for idx in combo:
                col_off, dst_r = opts[idx]
                if col_off == 1:
                    if deg_cp1[dst_r] + cp1_inc[dst_r] + 1 > 2:
                        valid = False
                        break
                    cp1_inc[dst_r] += 1
                    cp1_back.append(dst_r)
                else:
                    if cp2_inc[dst_r] + 1 > 2:
                        valid = False
                        break
                    cp2_inc[dst_r] += 1
                    cp2_back.append(dst_r)
                new_edges.append((r, col_off, dst_r))
            if valid:
                yield from rec(r + 1, cp1_inc, cp2_inc, partial_edges + new_edges)
            for dst_r in cp1_back:
                cp1_inc[dst_r] -= 1
            for dst_r in cp2_back:
                cp2_inc[dst_r] -= 1

    yield from rec(0, [0] * n, [0] * n, [])


def exact_marginals(n):
    """Retorna dict { (r1, c1, r2, c2) : count, ... }, total."""
    # Forward DP: v[state] = count of partial 2-factors leading to state
    # Backward DP: u[state] = count of completions from state
    # Total = sum_{state s} v[s] * u[s] * 1{s = transition_compatible} -- atraves
    # de transicoes (s, s').
    # Para cada aresta e que aparece em algum E_c na transicao (s, s'),
    # a contagem de 2-fatores contendo e = sum (sobre todas as ocorrencias)
    # de v[s] * u[s'].

    # Step 1: forward DP
    s0 = ((0,) * n, (0,) * n)
    v = {s0: 1}
    forward = [dict(v)]  # forward[c] = v antes do passo c

    for c in range(n):
        last = (c == n - 1)
        second_last = (c == n - 2)
        forbid_cp2 = second_last or last

        new_v = {}
        if last:
            for s, w in v.items():
                deg_c, deg_cp1 = s
                if all(d == 2 for d in deg_c):
                    s_next = (deg_cp1, (0,) * n)
                    new_v[s_next] = new_v.get(s_next, 0) + w
        else:
            for s, w in v.items():
                for s_next, _ in transitions_with_edges(s, n, forbid_cp2=forbid_cp2):
                    new_v[s_next] = new_v.get(s_next, 0) + w
        v = new_v
        forward.append(dict(v))

    total = v.get(((0,) * n, (0,) * n), 0)
    if total == 0:
        return {}, 0
    print(f"  total 2-fatores = {total}")

    # Step 2: backward DP -- u[state at step c] = # completions to s_final
    s_final = ((0,) * n, (0,) * n)
    u = {s_final: 1}
    backward = [None] * (n + 1)
    backward[n] = dict(u)

    for c in range(n - 1, -1, -1):
        last = (c == n - 1)
        second_last = (c == n - 2)
        forbid_cp2 = second_last or last

        new_u = defaultdict(int)
        # Para cada estado em forward[c], queremos saber quantas completioes
        # ele tem. Iteramos sobre forward[c] e enumeramos transicoes.
        for s, _ in forward[c].items():
            if last:
                deg_c, deg_cp1 = s
                if all(d == 2 for d in deg_c):
                    s_next = (deg_cp1, (0,) * n)
                    if s_next in u:
                        new_u[s] += u[s_next]
            else:
                for s_next, _ in transitions_with_edges(s, n, forbid_cp2=forbid_cp2):
                    if s_next in u:
                        new_u[s] += u[s_next]
        u = dict(new_u)
        backward[c] = u

    # Step 3: contar arestas
    edge_count = defaultdict(int)  # (r1, c1, r2, c2) sorted -> count
    for c in range(n - 1):  # passos que adicionam arestas; ultimo passo nao adiciona
        last = (c == n - 1)
        second_last = (c == n - 2)
        forbid_cp2 = second_last or last
        if last:
            continue

        v_c = forward[c]
        u_cp1 = backward[c + 1]

        for s, vw in v_c.items():
            for s_next, E_c in transitions_with_edges(s, n, forbid_cp2=forbid_cp2):
                uw = u_cp1.get(s_next, 0)
                if uw == 0:
                    continue
                contrib = vw * uw
                for (src_r, col_off, dst_r) in E_c:
                    edge = (min(src_r, dst_r), c, max(src_r, dst_r), c + col_off) \
                        if src_r != dst_r else (src_r, c, dst_r, c + col_off)
                    # canonical: vertex tuples
                    v1 = (src_r, c)
                    v2 = (dst_r, c + col_off)
                    key = tuple(sorted([v1, v2]))
                    edge_count[key] += contrib

    return dict(edge_count), total


# ---------------------------------------------------------------------------
# Helper: nivel L de uma aresta
# ---------------------------------------------------------------------------

def edge_level(n, v1, v2):
    r1, c1 = v1
    r2, c2 = v2
    L1 = min(r1, n - 1 - r1, c1, n - 1 - c1)
    L2 = min(r2, n - 1 - r2, c2, n - 1 - c2)
    return min(L1, L2)


def main(ns):
    for n in ns:
        print(f"\n=== n = {n} (marginais exatas de 2-fatores) ===")
        t0 = time.time()
        edge_count, total = exact_marginals(n)
        t1 = time.time()
        print(f"  contagem exata em {t1 - t0:.2f}s")
        print(f"  {len(edge_count)} arestas tem ocorrencia > 0")

        # Marginal por nivel L
        by_L = defaultdict(list)
        for (v1, v2), c in edge_count.items():
            p = c / total
            L = edge_level(n, v1, v2)
            by_L[L].append(p)

        print(f"  Marginal exata de 2-fatores por nivel L:")
        print(f"  {'L':>3} {'n_edges':>8} {'mean':>10} {'std':>10} {'min':>10} {'max':>10}")
        out = []
        for L in sorted(by_L):
            arr = np.array(by_L[L])
            print(f"  {L:>3} {len(arr):>8} {arr.mean():>10.4f} {arr.std():>10.4f}"
                  f" {arr.min():>10.4f} {arr.max():>10.4f}")
            out.append({"L": L, "n_edges": len(arr),
                        "mean": float(arr.mean()),
                        "std": float(arr.std()),
                        "min": float(arr.min()),
                        "max": float(arr.max())})

        # Salvar
        with open(DATA / f"exact_2factor_marginals_n{n}.json", "w") as f:
            json.dump({
                "n": n,
                "total_2factors": total,
                "by_level": out,
                "edges": [
                    {"src": list(v1), "dst": list(v2),
                     "count": c, "marginal": c / total,
                     "level_L": edge_level(n, v1, v2)}
                    for (v1, v2), c in edge_count.items()
                ],
            }, f, indent=2)
        print(f"  salvo em data/exact_2factor_marginals_n{n}.json")


if __name__ == "__main__":
    ns = [int(x) for x in sys.argv[1:]] if len(sys.argv) > 1 else [6]
    main(ns)
