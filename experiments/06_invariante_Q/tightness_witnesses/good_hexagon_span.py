#!/usr/bin/env python3
"""
good_hexagon_span.py
====================
Fase 3, item 4.3: os hexagonos que ADMITEM switcher completo ainda geram
Z_bulk?

Por que importa. Numa prova por contradicao, R nao e' escolhido por nos: ele
"acende" (interseccao impar) algum hexagono, e nao controlamos qual. Se apenas
os hexagonos BONS (os que passam em (S3)) ja' geram Z_bulk, entao todo
R fora de X acende pelo menos um hexagono bom -- e a existencia de switcher
para os hexagonos ruins deixa de ser necessaria.

Procedimento: varre hexagonos do bulk em ordem aleatoria, testa (S2.b)+(S3)
em cada um, e acumula o rank do span dos que passam. Para assim que o rank
atinge dim Z_bulk (basta gerar, nao precisa testar todos).

CLI:
  python good_hexagon_span.py --n 8 --max-hex 800
"""
from __future__ import annotations

import argparse
import json
import random
import time
from pathlib import Path
from typing import Dict, List, Set

from bulk_perp_dims import build_rect, gf2_rank_ints, rect_corners
from switcher_search import (build, build_switchers, color, ham_path,
                             hexagons, label, verify_switcher)

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"


def run(n: int, max_hex: int, max_internal: int, max_w: int,
        budget: int, seed: int) -> dict:
    t0 = time.time()
    adj = build(n)
    _, edges = build_rect(n, n)
    eidx = {e: i for i, e in enumerate(edges)}
    total = n * n
    cs = set(rect_corners(n, n))
    rng = random.Random(seed)

    V_bulk, E_bulk = total - 4, len(edges) - 8
    dim_Zbulk = E_bulk - V_bulk + 1

    hexa = hexagons(adj, cs, total)
    rng.shuffle(hexa)
    hexa = hexa[:max_hex]

    print(f"\n G_{n}x{n}:  dim Z_bulk = {dim_Zbulk}, "
          f"hexagonos varridos = {len(hexa)}")

    good_vecs: List[int] = []
    n_good = n_bad = 0
    rank_now = 0
    spanned_after = None

    for h in hexa:
        cands: List[dict] = []
        for rot in range(3):
            cands += build_switchers(adj, h, rot, max_internal, max_w)
        cands.sort(key=lambda w: w["VW"])

        solved = False
        for w in cands:
            v1, v4 = w["poles"]
            removed = (set(w["cycle"]) | set(w["P2_internal"])
                       | set(w["P3_internal"])) - {v1, v4}
            allowed = set(range(total)) - removed
            ca = sum(1 for x in allowed if color(x, n) == 0)
            if not (color(v1, n) != color(v4, n) and ca * 2 == len(allowed)):
                continue
            if any(sum(1 for u in adj[c] if u in allowed) < 2 for c in cs):
                continue
            p, _exhausted, _nodes = ham_path(adj, allowed, v1, v4, budget)
            if p is not None and not verify_switcher(adj, w, p, allowed, total):
                solved = True
                break

        if solved:
            n_good += 1
            x = 0
            for a, b in zip(h, h[1:] + h[:1]):
                x ^= 1 << eidx[(min(a, b), max(a, b))]
            good_vecs.append(x)
            rank_now = gf2_rank_ints(good_vecs)
            if rank_now == dim_Zbulk and spanned_after is None:
                spanned_after = n_good
                print(f"   -> hexagonos BONS geram Z_bulk apos {n_good} deles "
                      f"({time.time() - t0:.0f}s)")
                break
        else:
            n_bad += 1

    spans = rank_now == dim_Zbulk
    print(f"   bons={n_good}  ruins={n_bad}  "
          f"rank(span dos bons) = {rank_now}/{dim_Zbulk}  "
          f"{'GERA' if spans else 'NAO GERA'}")

    return {
        "n": n, "dim_Z_bulk": dim_Zbulk,
        "hexagons_scanned": n_good + n_bad,
        "good": n_good, "bad": n_bad,
        "rank_good_span": rank_now,
        "good_hexagons_span_Zbulk": spans,
        "spanned_after_n_good": spanned_after,
        "params": {"max_hex": max_hex, "max_internal": max_internal,
                   "max_w": max_w, "budget": budget, "seed": seed},
        "elapsed_s": round(time.time() - t0, 1),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ns", type=str, default="8,10,12")
    ap.add_argument("--max-hex", type=int, default=800)
    ap.add_argument("--max-internal", type=int, default=4)
    ap.add_argument("--max-w", type=int, default=20)
    ap.add_argument("--budget", type=int, default=400_000)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()

    print("=" * 68)
    print(" Os hexagonos que passam em (S3) geram Z_bulk?")
    print("=" * 68)
    out = [run(n, a.max_hex, a.max_internal, a.max_w, a.budget, a.seed)
           for n in [int(x) for x in a.ns.split(",")]]

    print("\n" + "=" * 68)
    for r in out:
        print(f"   n={r['n']:>3}: {'SIM' if r['good_hexagons_span_Zbulk'] else 'NAO'}"
              f"  (rank {r['rank_good_span']}/{r['dim_Z_bulk']}, "
              f"bons={r['good']}, ruins={r['bad']})")

    DATA.mkdir(parents=True, exist_ok=True)
    p = DATA / "good_hexagon_span.json"
    p.write_text(json.dumps(out, indent=2))
    print(f"\n-> {p}")


if __name__ == "__main__":
    main()
