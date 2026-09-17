#!/usr/bin/env python3
"""
shape_taxonomy.py
=================
Categoriza as 128 FORMAS de hexagono do bulk (6x6) quanto ao papel delas na
conectividade do grafo de flips sobre os 9862 tours.

Cuidado metodologico: o "101 formas" do single_shape_flip.py depende da ORDEM
(formas ordenadas por numero de hexagonos). As 27 restantes NAO sao
intrinsecamente especiais. Aqui calculamos noceos independentes de ordem:

  A. realizada      : a forma habilita >=1 flip valido
  B. essencial      : remover a forma do conjunto COMPLETO desconecta o grafo
  C. redundante     : remover a forma nao muda nada (ainda conexo)
  D. cobertura      : quantos flips / quantos tours a forma toca
  E. minimo guloso  : menor conjunto de formas que conecta (guloso por ganho)

  python shape_taxonomy.py
"""
from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

import numpy as np

from flip_graph import (build_board, hexagons_bulk, alternating_matchings,
                        cycle_mask, tour_mask)

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
N = 6


def normalize(cells):
    r0 = min(r for r, _ in cells)
    c0 = min(c for _, c in cells)
    return tuple(sorted((r - r0, c - c0) for r, c in cells))


def d4_min(cells, n):
    best = None
    cur = list(cells)
    for _ in range(4):
        cur = [(c, n - 1 - r) for (r, c) in cur]
        for img in (normalize(cur), normalize([(r, n - 1 - c) for r, c in cur])):
            if best is None or img < best:
                best = img
    return best


def ncomp(T, edge_lists):
    par = list(range(T))

    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x

    for el in edge_lists:
        for a, b in el:
            ra, rb = find(a), find(b)
            if ra != rb:
                par[ra] = rb
    sizes = defaultdict(int)
    for x in range(T):
        sizes[find(x)] += 1
    s = sorted(sizes.values(), reverse=True)
    return len(s), s[0]


def main():
    adj, edges, eidx = build_board(N)
    tours = np.load(ROOT.parent / "complex_orbit" / "data" / "tours_closed_6x6.npy")
    masks = sorted({tour_mask(t, eidx) for t in tours})
    idx = {m: i for i, m in enumerate(masks)}
    T = len(masks)

    hexa = hexagons_bulk(N, adj)
    shape_of, orbit_of = {}, {}
    for cyc in hexa:
        cells = [divmod(v, N) for v in cyc]
        shape_of[cyc] = normalize(cells)
        orbit_of[cyc] = d4_min(cells, N)

    # ---- calcula TODOS os flips validos uma unica vez, agrupados por forma
    edges_by_shape = defaultdict(list)
    hexcount = defaultdict(int)
    matching_ok = defaultdict(int)
    for cyc in hexa:
        sh = shape_of[cyc]
        hexcount[sh] += 1
        cm = cycle_mask(cyc, eidx)
        m1, m2 = alternating_matchings(cyc, eidx)
        for m in masks:
            inter = m & cm
            if inter != m1 and inter != m2:
                continue
            matching_ok[sh] += 1
            j = idx.get(m ^ cm)
            if j is not None and idx[m] < j:
                edges_by_shape[sh].append((idx[m], j))

    shapes = sorted(hexcount)
    print(f"tours={T}  hexagonos={len(hexa)}  formas={len(shapes)}  "
          f"orbitasD4={len({orbit_of[c] for c in hexa})}")

    base_nc, base_big = ncomp(T, edges_by_shape.values())
    print(f"conjunto COMPLETO: comps={base_nc} maior={base_big}")

    # ---- A/B/C/D por forma
    rows = []
    for sh in shapes:
        el = edges_by_shape[sh]
        others = [edges_by_shape[s] for s in shapes if s != sh]
        nc_wo, big_wo = ncomp(T, others)
        touched = len({v for e in el for v in e})
        r, c = max(x for x, _ in sh) + 1, max(y for _, y in sh) + 1
        rows.append({
            "shape": sh, "bbox": f"{r}x{c}", "n_hex": hexcount[sh],
            "orbit_size": len({orbit_of[cy] for cy in hexa
                               if shape_of[cy] == sh}),
            "matching_ok": matching_ok[sh],
            "n_flips": len(el), "tours_touched": touched,
            "realized": len(el) > 0,
            "comps_without": nc_wo, "biggest_without": big_wo,
            "essential": nc_wo > base_nc,
        })

    nao_real = [r for r in rows if not r["realized"]]
    essen = [r for r in rows if r["essential"]]
    redun = [r for r in rows if r["realized"] and not r["essential"]]

    print(f"\nA. NAO realizadas (0 flips)      : {len(nao_real)}")
    print(f"B. ESSENCIAIS (remocao desconecta): {len(essen)}")
    print(f"C. redundantes                    : {len(redun)}")

    if nao_real:
        print("\n  --- formas nao realizadas ---")
        for r in sorted(nao_real, key=lambda r: (-r["n_hex"], r["shape"])):
            print(f"   bbox={r['bbox']}  hex={r['n_hex']:2d}  "
                  f"matching_ok={r['matching_ok']:3d}  {r['shape']}")

    if essen:
        print("\n  --- formas essenciais ---")
        for r in sorted(essen, key=lambda r: -r["comps_without"]):
            print(f"   bbox={r['bbox']}  hex={r['n_hex']:2d}  flips={r['n_flips']:5d}  "
                  f"sem ela: comps={r['comps_without']:4d} maior={r['biggest_without']}")

    # ---- E. minimo guloso
    print("\nE. minimo guloso (maximiza reducao de componentes)")
    chosen, cur = [], []
    while True:
        best = None
        for sh in shapes:
            if sh in chosen:
                continue
            nc, big = ncomp(T, cur + [edges_by_shape[sh]])
            if best is None or (nc, -big) < (best[0], -best[1]):
                best = (nc, big, sh)
        nc, big, sh = best
        chosen.append(sh)
        cur.append(edges_by_shape[sh])
        if len(chosen) <= 5 or nc == 1 or len(chosen) % 5 == 0:
            print(f"   {len(chosen):3d} formas: comps={nc:5d} maior={big:5d}")
        if nc == 1:
            break
    print(f"   -> CONEXO com {len(chosen)} formas (vs 101 na ordem por tamanho)")

    out = {
        "n_tours": T, "n_shapes": len(shapes),
        "base_components": base_nc,
        "n_nao_realizadas": len(nao_real), "n_essenciais": len(essen),
        "n_redundantes": len(redun),
        "greedy_min_shapes": len(chosen),
        "shapes": [{k: (list(map(list, v)) if k == "shape" else v)
                    for k, v in r.items()} for r in rows],
    }
    DATA.mkdir(exist_ok=True)
    (DATA / "shape_taxonomy.json").write_text(json.dumps(out, indent=2))
    print(f"\n-> {DATA/'shape_taxonomy.json'}")


if __name__ == "__main__":
    main()
