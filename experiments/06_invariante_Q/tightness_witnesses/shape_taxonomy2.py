#!/usr/bin/env python3
"""Complementos a shape_taxonomy.py: (i) onde vivem os 12 hexagonos nunca
realizados, por forma; (ii) cota inferior de cobertura para o minimo."""
from __future__ import annotations
from collections import defaultdict
from pathlib import Path
import json
import numpy as np
from flip_graph import (build_board, hexagons_bulk, alternating_matchings,
                        cycle_mask, tour_mask)

ROOT = Path(__file__).resolve().parent
N = 6

def normalize(cells):
    r0 = min(r for r, _ in cells); c0 = min(c for _, c in cells)
    return tuple(sorted((r - r0, c - c0) for r, c in cells))

adj, edges, eidx = build_board(N)
tours = np.load(ROOT.parent / "complex_orbit" / "data" / "tours_closed_6x6.npy")
masks = sorted({tour_mask(t, eidx) for t in tours})
idx = {m: i for i, m in enumerate(masks)}
T = len(masks)
hexa = hexagons_bulk(N, adj)

per_hex = {}
touch_by_shape = defaultdict(set)
for cyc in hexa:
    cm = cycle_mask(cyc, eidx); m1, m2 = alternating_matchings(cyc, eidx)
    sh = normalize([divmod(v, N) for v in cyc])
    nf = 0
    for m in masks:
        it = m & cm
        if it != m1 and it != m2:
            continue
        j = idx.get(m ^ cm)
        if j is not None:
            nf += 1
            touch_by_shape[sh].add(idx[m]); touch_by_shape[sh].add(j)
    per_hex[cyc] = (sh, nf)

dead = [(sh, cyc) for cyc, (sh, nf) in per_hex.items() if nf == 0]
by_shape_dead = defaultdict(list)
tot_by_shape = defaultdict(int)
for cyc, (sh, nf) in per_hex.items():
    tot_by_shape[sh] += 1
for sh, cyc in dead:
    by_shape_dead[sh].append(cyc)

print(f"hexagonos nunca realizados: {len(dead)}   formas afetadas: {len(by_shape_dead)}")
print(f"{'bbox':6s} {'mortos/total':>13s}   forma")
for sh in sorted(by_shape_dead, key=lambda s: (-len(by_shape_dead[s]), s)):
    r = max(x for x, _ in sh) + 1; c = max(y for _, y in sh) + 1
    print(f"{r}x{c:<4d} {len(by_shape_dead[sh]):5d}/{tot_by_shape[sh]:<7d}  {sh}")
print("\n-> nenhuma forma tem TODOS os translados mortos"
      if all(len(by_shape_dead[s]) < tot_by_shape[s] for s in by_shape_dead)
      else "\n-> ALGUMA forma esta' totalmente morta")

# cota inferior: cobrir os 9862 tours (grau >= 1) exige pelo menos k formas
shapes = sorted(touch_by_shape, key=lambda s: -len(touch_by_shape[s]))
cov, k = set(), 0
while len(cov) < T:
    best = max(shapes, key=lambda s: len(touch_by_shape[s] - cov))
    cov |= touch_by_shape[best]; shapes.remove(best); k += 1
print(f"\ncobertura gulosa de todos os {T} tours: {k} formas "
      f"(e' uma cobertura VALIDA, logo OPT_cobertura <= {k}; nao e' cota "
      f"inferior. Conexidade exige cobrir, entao OPT_conexo >= OPT_cobertura)")
maxcov = max(len(v) for v in touch_by_shape.values())
print(f"maior cobertura de uma forma so': {maxcov} tours "
      f"-> minimo >= ceil({T}/{maxcov}) = {-(-T//maxcov)}")
