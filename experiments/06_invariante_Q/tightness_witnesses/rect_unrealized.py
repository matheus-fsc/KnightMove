#!/usr/bin/env python3
"""
rect_unrealized.py
==================
Quais hexagonos do bulk NAO sao flip-realizaveis num retangulo R x C, e como
se organizam sob o grupo de simetria do tabuleiro?

Para R != C o grupo e' o de Klein de ordem 4 {id, flip-h, flip-v, rot180},
nao D4. A pergunta e' a mesma que em 6x6: os excepcionais formam orbitas
fechadas?

Reutiliza a enumeracao exaustiva de `flip_graph_rect.py`.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

from flip_graph import alternating_matchings, cycle_mask
from flip_graph_rect import build_rect, enumerate_tours, hexagons_bulk_rect

DATA = Path(__file__).resolve().parent / "data"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rows", type=int, default=6)
    ap.add_argument("--cols", type=int, default=7)
    args = ap.parse_args()
    R, C = args.rows, args.cols

    t0 = time.time()
    adj, edges, eidx = build_rect(R, C)
    hexes = hexagons_bulk_rect(R, C, adj)
    masks = enumerate_tours(R, C, adj, eidx)
    print(f"[{R}x{C}] tours={len(masks)} hex={len(hexes)} "
          f"({time.time()-t0:.0f}s)", flush=True)

    bad = []
    bad_match_only = []
    for cyc in hexes:
        cm = cycle_mask(cyc, eidx)
        m0, m1 = alternating_matchings(cyc, eidx)
        hit = [tm for tm in masks if (tm & cm) in (m0, m1)]
        if not hit:
            bad.append(cyc)
            continue
        # passa no emparelhamento: falha so' se todo tau XOR C for desconexo
        ok = False
        for tm in hit:
            nm = tm ^ cm
            nb = [[] for _ in range(R * C)]
            m = nm
            while m:
                b = m & -m
                u, v = edges[b.bit_length() - 1]
                nb[u].append(v)
                nb[v].append(u)
                m ^= b
            prev, cur, cnt = -1, 0, 0
            while True:
                a, b2 = nb[cur]
                nxt = a if a != prev else b2
                prev, cur = cur, nxt
                cnt += 1
                if cur == 0 or cnt > R * C:
                    break
            if cnt == R * C:
                ok = True
                break
        if not ok:
            bad_match_only.append(cyc)

    def cells(cyc):
        return frozenset(divmod(v, C) for v in cyc)

    G = [lambda r, c: (r, c),
         lambda r, c: (r, C - 1 - c),
         lambda r, c: (R - 1 - r, c),
         lambda r, c: (R - 1 - r, C - 1 - c)]
    allbad = {cells(c) for c in bad} | {cells(c) for c in bad_match_only}
    seen, orbits = set(), []
    for s in allbad:
        if s in seen:
            continue
        orb = {frozenset(g(*p) for p in s) for g in G}
        orbits.append((len(orb), len(orb & allbad)))
        seen |= orb

    res = {
        "R": R, "C": C, "n_tours": len(masks), "n_hexagons_bulk": len(hexes),
        "n_unrealized": len(bad) + len(bad_match_only),
        "n_no_matching": len(bad),
        "n_matching_only": len(bad_match_only),
        "orbits_klein": [{"orbit_size": a, "inside": b} for a, b in orbits],
        "all_orbits_closed": all(a == b for a, b in orbits),
        "unrealized_cells": [sorted(divmod(v, C) for v in c)
                             for c in bad + bad_match_only],
        "seconds": round(time.time() - t0, 1),
    }
    out = DATA / f"rect_unrealized_{R}x{C}.json"
    out.write_text(json.dumps(res, indent=2))
    print(json.dumps({k: v for k, v in res.items()
                      if k != "unrealized_cells"}, indent=2))
    for c in res["unrealized_cells"]:
        print("   ", c)
    print(f"-> {out}")


if __name__ == "__main__":
    main()
