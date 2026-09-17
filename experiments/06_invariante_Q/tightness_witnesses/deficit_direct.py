#!/usr/bin/env python3
"""
deficit_direct.py
=================
`deficit(R,C)` exato por enumeracao de tours, SEM o laco sobre hexagonos.

    deficit = beta_1 - rk(Ham) = beta_1 - (dim span{tau + tau_0} + 1)

O laco de `missing_directions.py` custa O(hexagonos x tours); este custa
O(tours) operacoes de base. Para saber apenas o deficit -- por exemplo, se
acrescentar COLUNAS conserta um tabuleiro de 5 linhas -- e' o caminho barato.

CLI:
  python deficit_direct.py --rows 5 --cols 10
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

from coset_tests import Basis as IntBasis, Board, Ctx, Q_paper, zbulk_basis
from flip_graph_rect import build_rect, enumerate_tours

DATA = Path(__file__).resolve().parent / "data"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rows", type=int, default=5)
    ap.add_argument("--cols", type=int, default=10)
    args = ap.parse_args()
    R, C = args.rows, args.cols

    t0 = time.time()
    b = Board(R, C)
    dim_zb = IntBasis(zbulk_basis(Ctx(R, C))).rank
    q = Q_paper(b)
    print(f"[{R}x{C}] V={b.V} E={b.E} beta1={b.beta1} dim_ZBulk={dim_zb} "
          f"Q={q} c_Punc={b.punc_components()}", flush=True)

    adj, edges, eidx = build_rect(R, C)
    masks = enumerate_tours(R, C, adj, eidx)
    t_enum = time.time() - t0
    print(f"[{R}x{C}] tours (EXAUSTIVO) = {len(masks)}  ({t_enum:.0f}s)",
          flush=True)

    t0m = masks[0]
    B = IntBasis(tm ^ t0m for tm in masks[1:])
    Zb = IntBasis(zbulk_basis(Ctx(R, C)))
    # invariante: diferencas de tours vivem em Z_bulk
    assert all(Zb.contains(v) for v in B.vectors()), \
        "tau+tau' fora de Z_bulk: impossivel (ambos tem y=1)"
    assert B.rank <= dim_zb

    rk_ham = B.rank + 1
    deficit = b.beta1 - rk_ham
    res = {
        "R": R, "C": C, "V": b.V, "E": b.E, "beta1": b.beta1,
        "dim_ZBulk": dim_zb, "Q": q, "n_tours": len(masks),
        "dim_span_diffs": B.rank, "rank_Ham": rk_ham,
        "deficit": deficit, "deficit_minus_Q": deficit - q,
        "codim_in_ZBulk": dim_zb - B.rank,
        "tight": deficit == q,
        "enum_seconds": round(t_enum, 1),
        "seconds": round(time.time() - t0, 1),
    }
    out = DATA / f"deficit_direct_{R}x{C}.json"
    out.write_text(json.dumps(res, indent=2))
    print("=" * 60)
    for k, v in res.items():
        print(f" {k:20s}: {v}")
    print(f" -> {out}")


if __name__ == "__main__":
    main()
