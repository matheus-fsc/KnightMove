#!/usr/bin/env python3
"""
single_shape_flip.py
====================
Pergunta: um UNICO tipo de hexagono (uma unica FORMA, transladada por todo o
bulk) ja' gera todos os 9862 tours do 6x6 a partir de uma seed?

Isto e' a versao precisa da ideia "um unico losango como semente": a semente e'
um TOUR, e o movimento e' o XOR com hexagonos de UMA forma so'.

Para cada forma (classe de translacao) e cada orbita D4, restringe o conjunto
de movimentos e mede o numero de componentes conexas do grafo de flips.

  python single_shape_flip.py
"""
from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

import numpy as np

from flip_graph import (build_board, hexagons_bulk, alternating_matchings,
                        cycle_mask, tour_mask, is_hamiltonian)

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
N = 6


def d4_images(cells, n):
    """8 imagens dihedrais de um conjunto de celulas (r,c)."""
    out = []
    cur = list(cells)
    for _ in range(4):
        cur = [(c, n - 1 - r) for (r, c) in cur]          # rotacao 90
        out.append(normalize(cur))
        out.append(normalize([(r, n - 1 - c) for (r, c) in cur]))  # + espelho
    return out


def normalize(cells):
    r0 = min(r for r, _ in cells)
    c0 = min(c for _, c in cells)
    return tuple(sorted((r - r0, c - c0) for r, c in cells))


def main():
    adj, edges, eidx = build_board(N)
    tours = np.load(ROOT.parent / "complex_orbit" / "data" / "tours_closed_6x6.npy")

    masks = sorted({tour_mask(t, eidx) for t in tours})
    idx = {m: i for i, m in enumerate(masks)}
    T = len(masks)
    print(f"tours distintos: {T}")

    hexa = hexagons_bulk(N, adj)
    print(f"hexagonos do bulk: {len(hexa)}")

    # agrupa por forma (translacao) e por orbita D4
    by_shape = defaultdict(list)
    by_orbit = defaultdict(list)
    for cyc in hexa:
        cells = [divmod(v, N) for v in cyc]
        sh = normalize(cells)
        by_shape[sh].append(cyc)
        by_orbit[min(d4_images(cells, N))].append(cyc)
    print(f"formas (translacao): {len(by_shape)}   orbitas D4: {len(by_orbit)}")

    def components(cyc_list):
        """union-find sobre os T tours usando so' estes hexagonos."""
        par = list(range(T))

        def find(x):
            while par[x] != x:
                par[x] = par[par[x]]
                x = par[x]
            return x

        nflips = 0
        for cyc in cyc_list:
            cm = cycle_mask(cyc, eidx)
            m1, m2 = alternating_matchings(cyc, eidx)
            for m in masks:
                inter = m & cm
                if inter != m1 and inter != m2:
                    continue
                nm = m ^ cm
                j = idx.get(nm)
                if j is None:
                    if not is_hamiltonian(nm, N, edges):
                        continue
                    continue
                nflips += 1
                a, b = find(idx[m]), find(j)
                if a != b:
                    par[a] = b
        sizes = defaultdict(int)
        for x in range(T):
            sizes[find(x)] += 1
        s = sorted(sizes.values(), reverse=True)
        return len(s), s[0], nflips // 2

    rows = []
    for name, groups in (("forma", by_shape), ("orbitaD4", by_orbit)):
        print(f"\n=== restringindo a UMA {name} ===")
        res = []
        for key, cyc_list in groups.items():
            nc, big, nf = components(cyc_list)
            res.append((nc, big, nf, len(cyc_list), key))
        res.sort()
        for nc, big, nf, ncy, key in res[:8]:
            frac = 100.0 * big / T
            print(f"  {name}={str(key)[:38]:40s} hex={ncy:3d} flips={nf:6d} "
                  f"comps={nc:5d} maior={big:5d} ({frac:5.1f}%)")
        rows.append({"nivel": name, "melhor_comps": res[0][0],
                     "melhor_maior": res[0][1], "n_grupos": len(groups)})

    # cumulativo: quantas formas bastam?
    print("\n=== cumulativo (formas ordenadas por tamanho da componente) ===")
    order = sorted(by_shape.items(), key=lambda kv: -len(kv[1]))
    acc = []
    for i, (key, cl) in enumerate(order, 1):
        acc += cl
        nc, big, nf = components(acc)
        print(f"  {i:2d} formas ({len(acc):3d} hex): comps={nc:5d} maior={big:5d}")
        if nc == 1:
            print(f"  -> CONEXO com {i} formas")
            rows.append({"formas_para_conexo": i})
            break

    DATA.mkdir(exist_ok=True)
    (DATA / "single_shape_flip.json").write_text(json.dumps(rows, indent=2))
    print(f"\n-> {DATA/'single_shape_flip.json'}")


if __name__ == "__main__":
    main()
