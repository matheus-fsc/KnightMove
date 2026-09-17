#!/usr/bin/env python3
"""
n8_strong_statement.py
======================
Vale em n=8 o ENUNCIADO FORTE que falha em n=6?

    "todo hexagono do bulk e' flip-realizavel"
    (existe tour tau com tau∩C alternado e tau XOR C conexo)

Estrategia em dois estagios, para nao pagar busca onde a bola ja' responde:

  1. bola de BFS por flips -> certifica positivamente uma fatia grande
     (cada hexagono realizado vem com o par (tau, tau XOR C) verificado);
  2. busca hamiltoniana RESTRITA e EXAUSTIVA no residuo -- decisiva nos dois
     sentidos, com o orcamento reportado quando estoura.

O estagio 2 usa `flip_realizable.test_hexagon`, validado contra a verdade
exaustiva de n=6 (532/532 de acordo, 0 inconclusivos).
"""
from __future__ import annotations

import argparse
import json
import random
import time
from collections import deque
from pathlib import Path

from flip_graph import (
    alternating_matchings, build_board, cycle_mask, hexagons_bulk,
    is_hamiltonian, tour_mask,
)
from flip_graph_n8 import warnsdorff_tour
from flip_realizable import search_flip

DATA = Path(__file__).resolve().parent / "data"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=8)
    ap.add_argument("--ball", type=int, default=1_500_000)
    ap.add_argument("--budget", type=int, default=20_000_000)
    ap.add_argument("--seed", type=int, default=11)
    args = ap.parse_args()

    n = args.n
    t0 = time.time()
    adj, edges, eidx = build_board(n)
    hexes = hexagons_bulk(n, adj)
    hx = [(cycle_mask(c, eidx), *alternating_matchings(c, eidx)) for c in hexes]
    by_mask = {cycle_mask(c, eidx): c for c in hexes}

    rng = random.Random(args.seed)
    root = tour_mask(warnsdorff_tour(n, adj, rng), eidx)
    seen = {root}
    q = deque([root])
    realized = set()
    while q and len(seen) < args.ball:
        m = q.popleft()
        for cm, m0, m1 in hx:
            it = m & cm
            if it != m0 and it != m1:
                continue
            nm = m ^ cm
            if not is_hamiltonian(nm, n, edges):
                continue
            realized.add(cm)
            if nm not in seen:
                seen.add(nm)
                q.append(nm)
    print(f"[estagio 1] bola={len(seen)}  realizados={len(realized)}/{len(hexes)}"
          f"  ({time.time()-t0:.0f}s)", flush=True)

    residue = [by_mask[cm] for cm in by_mask if cm not in realized]
    print(f"[estagio 2] residuo: {len(residue)} hexagonos, busca restrita "
          f"exaustiva (orcamento {args.budget} nos por emparelhamento)",
          flush=True)

    ok_flip = []
    not_flip = []
    incon = []
    rng2 = random.Random(args.seed + 1)
    for k, cyc in enumerate(residue):
        # fase A: barata, aleatorizada -- so' pode dar positivo
        found, _ = search_flip(n, adj, edges, eidx, cyc,
                               args.budget // 20, rng2)
        ex = True
        if not found:
            # fase B: determinista e sem `want`, decisiva nos dois sentidos
            found, ex = search_flip(n, adj, edges, eidx, cyc, args.budget)
        if found:
            ok_flip.append(cyc)
        elif ex:
            not_flip.append(cyc)
        else:
            incon.append(cyc)
        if (k + 1) % 5 == 0:
            print(f"    {k+1}/{len(residue)}  flip={len(ok_flip)} "
                  f"NAO-flip={len(not_flip)} inconclusivo={len(incon)}  "
                  f"({time.time()-t0:.0f}s)", flush=True)

    res = {
        "n": n, "ball": len(seen), "budget": args.budget,
        "n_hexagons_bulk": len(hexes),
        "realized_in_ball": len(realized),
        "residue": len(residue),
        "residue_flip_realizable": len(ok_flip),
        "residue_not_flip_realizable": len(not_flip),
        "residue_inconclusive": len(incon),
        "total_flip_realizable": len(realized) + len(ok_flip),
        "strong_statement_holds": len(not_flip) == 0 and len(incon) == 0,
        "not_flip_cells": [[divmod(v, n) for v in c] for c in not_flip],
        "inconclusive_cells": [[divmod(v, n) for v in c] for c in incon],
        "seconds": round(time.time() - t0, 1),
    }
    out = DATA / f"n{n}_strong_statement.json"
    out.write_text(json.dumps(res, indent=2))
    print()
    print("=" * 66)
    print(f" ENUNCIADO FORTE EM n={n}")
    print("=" * 66)
    print(f" hexagonos do bulk           : {res['n_hexagons_bulk']}")
    print(f" flip-realizaveis (total)    : {res['total_flip_realizable']}")
    print(f"   via bola                  : {res['realized_in_ball']}")
    print(f"   via busca restrita        : {res['residue_flip_realizable']}")
    print(f" NAO flip-realizaveis (prova): {res['residue_not_flip_realizable']}")
    print(f" inconclusivos (orcamento)   : {res['residue_inconclusive']}")
    print(f" ENUNCIADO FORTE VALE?       : {res['strong_statement_holds']}")
    print(f" tempo {res['seconds']}s -> {out}")


if __name__ == "__main__":
    main()
