#!/usr/bin/env python3
"""
hexagon_span.py
===============
Fase 3, tarefa extra: os hexagonos do bulk GERAM Z_bulk?

Motivacao. O Lema B garante que todo R em C_n^perp \\ X tem interseccao IMPAR
com algum ciclo do bulk -- mas nao diz o comprimento desse ciclo. A construcao
do switcher so' funciona com k IMPAR (restricao de cor: v_1 e v_{k+1} precisam
de cores opostas), isto e', |C| = 2k com k impar, ou seja |C| ≡ 2 (mod 4).
Se |C| ≡ 0 (mod 4), (S3) e' impossivel por paridade de cor.

Portanto a pergunta que fecha (S2.a) e':

    os ciclos do bulk de comprimento ≡ 2 (mod 4) geram Z_bulk?

Se os HEXAGONOS sozinhos ja geram, entao todo R fora de X tem interseccao
impar com algum hexagono do bulk -- exatamente o objeto que switcher_search.py
mostrou ser construtivel.

CLI:
  python hexagon_span.py --max-n 12
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Dict, List, Set, Tuple

from bulk_perp_dims import build_rect, gf2_rank_ints, rect_corners

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"


def cycles_of_length(adj, forbidden: Set[int], total: int,
                     length: int) -> List[Tuple[int, ...]]:
    """Ciclos simples de comprimento `length` evitando `forbidden`."""
    out = []
    for v1 in range(total):
        if v1 in forbidden:
            continue

        def rec(path: List[int]):
            cur = path[-1]
            if len(path) == length:
                if v1 in adj[cur] and path[1] < path[-1]:
                    out.append(tuple(path))
                return
            for u in adj[cur]:
                if u <= v1 or u in forbidden or u in path:
                    continue
                path.append(u)
                rec(path)
                path.pop()

        rec([v1])
    return out


def analyse(n: int, lengths: List[int]) -> dict:
    adj_list, edges = build_rect(n, n)
    adj = [sorted(set(a)) for a in adj_list]
    total = n * n
    eidx = {e: i for i, e in enumerate(edges)}
    cs = set(rect_corners(n, n))

    V_bulk = total - 4
    E_bulk = len(edges) - 8
    dim_Zbulk = E_bulk - V_bulk + 1        # bulk conexo (teorema)

    res = {"n": n, "dim_Z_bulk": dim_Zbulk, "by_length": {}}
    cumulative: List[int] = []
    print(f"\n G_{n}x{n}:  dim Z_bulk = {dim_Zbulk}")
    for L in lengths:
        cycs = cycles_of_length(adj, cs, total, L)
        vecs = []
        for cyc in cycs:
            x = 0
            for a, b in zip(cyc, cyc[1:] + cyc[:1]):
                x ^= 1 << eidx[(min(a, b), max(a, b))]
            vecs.append(x)
        r_alone = gf2_rank_ints(vecs)
        cumulative += vecs
        r_cum = gf2_rank_ints(cumulative)
        res["by_length"][str(L)] = {
            "n_cycles": len(cycs), "rank_alone": r_alone,
            "rank_cumulative": r_cum,
            "spans_Z_bulk": r_cum == dim_Zbulk,
        }
        mark = "  <-- GERA Z_bulk" if r_cum == dim_Zbulk else ""
        print(f"   |C|={L:>2}  ciclos={len(cycs):>6}  rank={r_alone:>4}  "
              f"rank acumulado={r_cum:>4}/{dim_Zbulk}{mark}")
    return res


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ns", type=str, default="6,8,10,12")
    ap.add_argument("--lengths", type=str, default="4,6,8,10")
    a = ap.parse_args()
    lengths = [int(x) for x in a.lengths.split(",")]

    print("=" * 68)
    print(" Os ciclos curtos do bulk geram Z_bulk?")
    print(" (so' |C| ≡ 2 mod 4 serve para o switcher: k impar)")
    print("=" * 68)
    t0 = time.time()
    out = []
    for n in [int(x) for x in a.ns.split(",")]:
        out.append(analyse(n, lengths))

    print("\n" + "=" * 68)
    print(" Resumo")
    print("=" * 68)
    for r in out:
        last = r["by_length"][str(lengths[-1])]
        status = ("GERA Z_bulk" if last["spans_Z_bulk"]
                  else f"NAO gera ({last['rank_cumulative']}/{r['dim_Z_bulk']})")
        print(f"   n={r['n']:>3}: comprimentos {lengths} -> {status}")
    if all(x % 4 == 2 for x in lengths):
        print("   (todos os comprimentos testados tem k IMPAR: utilizaveis)")

    DATA.mkdir(parents=True, exist_ok=True)
    p = DATA / "hexagon_span.json"
    p.write_text(json.dumps(out, indent=2))
    print(f"\n-> {p}  ({time.time() - t0:.1f}s)")


if __name__ == "__main__":
    main()
