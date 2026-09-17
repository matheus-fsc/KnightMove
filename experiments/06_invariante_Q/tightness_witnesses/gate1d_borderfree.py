#!/usr/bin/env python3
"""
gate1d_borderfree.py
====================
Passo 1 da ordem revisada: existem MUITAS formas livres de borda como o T5?
E elas, sozinhas, geram Z_bulk?

Contexto. T5 = (0,1)(0,3)(1,2)(2,0)(2,2)(2,4) certifica em 100% das posicoes,
inclusive encostada na borda -- nao tem efeito de canto algum. Se as formas com
essa propriedade ja' geram Z_bulk, o conjunto excepcional de 18 posicoes
canto-relativas DESAPARECE: toda posicao de toda forma geradora certifica, e a
uniformidade vira um argumento so'.

Procedimento: certifica TODAS as posicoes de TODAS as 136 formas, filtra as de
taxa 100%, e testa o rank do span delas contra dim Z_bulk.

Vale mesmo se falhar: mede se "livre de borda" e' raro (so' T5, sugerindo
estrutura especial da caixa 3x5) ou comum.

CLI:
  python gate1d_borderfree.py --n 10
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Dict, List, Tuple

from gate1_templates import shape_translation
from gate2_certify_templates import certify, cycle_from_cells
from hexagon_locality import Basis, build_sub, cyc_vec, cycles_len
from switcher_search import build, corners

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"


def run(n: int, max_internal: int, max_w: int, budget: int) -> dict:
    t0 = time.time()
    adjk = build(n)
    cs = set(corners(n))
    cells = [(r, c) for r in range(n) for c in range(n)]
    idx, adj, edges = build_sub(cells)
    rev = {i: c for c, i in idx.items()}
    eidx = {e: i for i, e in enumerate(edges)}
    corner_cells = {(0, 0), (0, n - 1), (n - 1, 0), (n - 1, n - 1)}
    dim_Zbulk = (len(edges) - 8) - (n * n - 4) + 1

    # agrupa hexagonos sem canto por forma
    by_shape: Dict[Tuple, List[Tuple[List, int]]] = {}
    for h in cycles_len(adj, len(idx), 6):
        cl = [rev[v] for v in h]
        if any(x in corner_cells for x in cl):
            continue
        by_shape.setdefault(shape_translation(cl), []).append((cl, cyc_vec(h, eidx)))

    print(f"\n n={n}: {len(by_shape)} formas, dim Z_bulk={dim_Zbulk}")

    recs = []
    for si, (shape, items) in enumerate(sorted(by_shape.items())):
        H = max(r for r, _ in shape)
        W = max(c for _, c in shape)
        ok_vecs, fail_sigs = [], []
        for cl, vec in items:
            r0 = min(r for r, _ in cl)
            c0 = min(c for _, c in cl)
            hexa = cycle_from_cells(adjk, [(r0 + r) * n + (c0 + c)
                                           for r, c in shape])
            good = hexa is not None and certify(adjk, n, hexa, cs,
                                                max_internal, max_w, budget)
            if good:
                ok_vecs.append(vec)
            else:
                top, left = r0, c0
                bot, right = n - 1 - (r0 + H), n - 1 - (c0 + W)
                fail_sigs.append((min(top, bot), 0 if top <= bot else 1,
                                  min(left, right), 0 if left <= right else 1))
        rate = len(ok_vecs) / len(items)
        recs.append({"shape": [list(x) for x in shape],
                     "box": [H + 1, W + 1],
                     "positions": len(items), "certified": len(ok_vecs),
                     "rate": rate, "border_free": rate == 1.0,
                     "fail_sigs": sorted(set(fail_sigs))})
        if (si + 1) % 20 == 0:
            print(f"   ... {si + 1}/{len(by_shape)} formas "
                  f"({time.time() - t0:.0f}s)")

    bf = [r for r in recs if r["border_free"]]
    print(f"\n   formas livres de borda: {len(bf)}/{len(recs)}")

    # rank do span das formas livres de borda
    b = Basis()
    for r in bf:
        shape = tuple(tuple(x) for x in r["shape"])
        for cl, vec in by_shape[shape]:
            b.add(vec)
    spans = b.rank == dim_Zbulk
    print(f"   rank(span das livres de borda) = {b.rank}/{dim_Zbulk}  "
          f"{'GERA' if spans else 'FALTA %d' % (dim_Zbulk - b.rank)}")

    # distribuicao de caixas entre as livres de borda
    boxes: Dict[str, int] = {}
    for r in bf:
        boxes[str(tuple(r["box"]))] = boxes.get(str(tuple(r["box"])), 0) + 1
    print(f"   caixas das livres de borda: {boxes}")

    return {"n": n, "dim_Z_bulk": dim_Zbulk, "n_shapes": len(recs),
            "n_border_free": len(bf), "rank_border_free": b.rank,
            "border_free_spans": spans, "boxes_border_free": boxes,
            "shapes": recs, "elapsed_s": round(time.time() - t0, 1)}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=10)
    ap.add_argument("--max-internal", type=int, default=4)
    ap.add_argument("--max-w", type=int, default=20)
    ap.add_argument("--budget", type=int, default=400_000)
    a = ap.parse_args()

    print("=" * 72)
    print(" Passo 1 — formas livres de borda geram Z_bulk?")
    print("=" * 72)
    res = run(a.n, a.max_internal, a.max_w, a.budget)
    DATA.mkdir(parents=True, exist_ok=True)
    p = DATA / f"gate1d_borderfree_n{a.n}.json"
    p.write_text(json.dumps(res, indent=2))
    print(f"\n-> {p}   ({res['elapsed_s']}s)")


if __name__ == "__main__":
    main()
