#!/usr/bin/env python3
"""
missing_directions.py
=====================
Por que o 5x8 perde exatamente 6 direcoes?

A conta que liga as 6 ao deficit:

    rk(Ham)  = dim span{tau + tau'} + 1
    deficit  = beta_1 - rk(Ham) = 3 + (dim Z_bulk - dim span{tau + tau'})

e, quando o grafo de flips e' CONEXO, span{tau+tau'} = span(hexagonos
realizados). Logo, no 5x8, `deficit - 3 = 47 - 41 = 6`: as 6 direcoes sao a
codimensao dos hexagonos realizados dentro de Z_bulk.

Este script identifica QUAIS hexagonos geram essas direcoes e mede a caixa
(bbox) das suas formas, para testar a hipotese:

    os hexagonos que gerariam as direcoes faltantes exigem gadgets que nao
    cabem num tabuleiro de largura 5.

Distingue os dois modos de falha de realizacao (como `rect_unrealized.py`):
  - `no_match`   : nenhum tour intersecta o hexagono num emparelhamento
                   alternado -- a forma nao e' sequer *alcancavel*;
  - `disconnect` : ha' tour compativel, mas todo XOR quebra em varios ciclos.

CLI:
  python missing_directions.py --rows 5 --cols 8
  python missing_directions.py --rows 6 --cols 6      # controle (tight)
  python missing_directions.py --rows 6 --cols 7      # controle (tight)
"""
from __future__ import annotations

import argparse
import json
import time
from collections import Counter
from pathlib import Path
from typing import Dict, List, Tuple

from coset_tests import Board, Basis as IntBasis, zbulk_basis, Ctx
from flip_graph import alternating_matchings, cycle_mask
from flip_graph_rect import build_rect, enumerate_tours, hexagons_bulk_rect

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"


# ---------------------------------------------------------------- geometria

def cells(cyc: Tuple[int, ...], C: int) -> List[Tuple[int, int]]:
    return [divmod(v, C) for v in cyc]


def bbox(cl: List[Tuple[int, int]]) -> Tuple[int, int]:
    rs = [r for r, _ in cl]
    cs = [c for _, c in cl]
    return (max(rs) - min(rs) + 1, max(cs) - min(cs) + 1)


def shape_key(cl: List[Tuple[int, int]]) -> Tuple[Tuple[int, int], ...]:
    """Forma canonica: transladada para a origem, ordenada. Sem simetrias --
    queremos distinguir orientacao, que e' o ponto quando um lado e' 5."""
    r0 = min(r for r, _ in cl)
    c0 = min(c for _, c in cl)
    return tuple(sorted((r - r0, c - c0) for r, c in cl))


def touches(cl: List[Tuple[int, int]], R: int, C: int) -> Dict[str, bool]:
    return {
        "top": any(r == 0 for r, _ in cl),
        "bottom": any(r == R - 1 for r, _ in cl),
        "left": any(c == 0 for _, c in cl),
        "right": any(c == C - 1 for _, c in cl),
    }


# ---------------------------------------------------------------- principal

def run(R: int, C: int) -> dict:
    t0 = time.time()
    adj, edges, eidx = build_rect(R, C)
    V, E = R * C, len(edges)
    hexes = hexagons_bulk_rect(R, C, adj)
    masks = enumerate_tours(R, C, adj, eidx)
    print(f"[{R}x{C}] V={V} E={E} tours={len(masks)} hex={len(hexes)} "
          f"({time.time()-t0:.0f}s)", flush=True)

    # Z_bulk pela definicao (nucleo dos funcionais de canto), NAO por hexagonos
    ctx = Ctx(R, C)
    Zb = IntBasis(zbulk_basis(ctx))
    dim_zb = Zb.rank
    print(f"[{R}x{C}] dim Z_bulk = {dim_zb} (independente dos hexagonos)")

    def hamil(mask: int) -> bool:
        nb: List[List[int]] = [[] for _ in range(V)]
        m = mask
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
            if cur == 0:
                return cnt == V
            if cnt > V:
                return False

    # --- classifica cada hexagono ---------------------------------------
    records = []
    for cyc in hexes:
        cm = cycle_mask(cyc, eidx)
        m0, m1 = alternating_matchings(cyc, eidx)
        n_match = 0
        realized = False
        for tm in masks:
            it = tm & cm
            if it != m0 and it != m1:
                continue
            n_match += 1
            if hamil(tm ^ cm):
                realized = True
                break
        if realized:
            mode = "realized"
        elif n_match == 0:
            mode = "no_match"
        else:
            mode = "disconnect"
        cl = cells(cyc, C)
        bb = bbox(cl)
        records.append({
            "cells": cl, "mask": cm, "mode": mode,
            "bbox": bb, "bbox_max_side": max(bb),
            "shape": shape_key(cl), "touches": touches(cl, R, C),
        })

    n_real = sum(1 for r in records if r["mode"] == "realized")
    print(f"[{R}x{C}] realizados={n_real}  no_match="
          f"{sum(1 for r in records if r['mode']=='no_match')}  disconnect="
          f"{sum(1 for r in records if r['mode']=='disconnect')}")

    # --- span dos realizados, e a codimensao ----------------------------
    Breal = IntBasis(r["mask"] for r in records if r["mode"] == "realized")
    Ball = IntBasis(r["mask"] for r in records)
    # INVARIANTE BARATA. Todo hexagono do bulk *esta* em Z_bulk, entao o rank de
    # qualquer conjunto deles nao pode passar de dim Z_bulk. Foi exatamente esta
    # desigualdade que o bug do `zbulk_basis` violava (dava dim Z_bulk = 31 com
    # rank(realizados) = 41), e que teria delatado o erro sem nenhum valor de
    # referencia externo. Qualquer script desta familia deve checa-la.
    assert Ball.rank <= dim_zb, (
        f"rank(hexagonos)={Ball.rank} > dim Z_bulk={dim_zb}: impossivel, "
        "os hexagonos do bulk sao elementos de Z_bulk")
    assert all(Zb.contains(r["mask"]) for r in records), \
        "algum hexagono do bulk caiu fora de Z_bulk"

    codim = dim_zb - Breal.rank
    print(f"[{R}x{C}] rank(realizados)={Breal.rank}  rank(todos)={Ball.rank}"
          f"  dim Z_bulk={dim_zb}  ->  CODIMENSAO = {codim}")

    # confere que span{tau+tau'} = span(realizados) (vale se flip-grafo conexo)
    t0m = masks[0]
    Bdiff = IntBasis(tm ^ t0m for tm in masks[1:])
    diffs_eq_real = (Bdiff.rank == Breal.rank and
                     all(Bdiff.contains(v) for v in Breal.vectors()))

    # --- quais hexagonos sao NECESSARIOS (nao-nulos no quociente) --------
    probe = IntBasis(r["mask"] for r in records if r["mode"] == "realized")
    needed = []          # geradores minimais escolhidos gulosamente
    nonzero = []         # todos os que sao != 0 no quociente
    for r in records:
        if r["mode"] == "realized":
            continue
        if probe.reduce(r["mask"]) != 0:
            nonzero.append(r)
    # Ordem de escolha: INTERIORES PRIMEIRO. A escolha gulosa na ordem de
    # enumeracao dava 6 geradores todos encostados numa borda, o que fazia o
    # teste de borda parecer cobrir as 6 direcoes -- quando `interior_only`
    # mostra que 2 delas tem origem interior. Priorizar os interiores expoe
    # essas 2 explicitamente, em vez de esconde-las.
    def n_borders(r) -> int:
        return sum(r["touches"].values())

    probe2 = IntBasis(r["mask"] for r in records if r["mode"] == "realized")
    for r in sorted(nonzero, key=n_borders):
        if probe2.add(r["mask"]):
            r = dict(r)
            r["n_borders"] = n_borders(r)
            r["interior"] = n_borders(r) == 0
            needed.append(r)

    print(f"[{R}x{C}] nao-realizados nao-nulos no quociente: {len(nonzero)}"
          f"   geradores escolhidos: {len(needed)}")

    # --- as medidas que testam a hipotese -------------------------------
    def stats(rs):
        if not rs:
            return {}
        return {
            "n": len(rs),
            "bbox_hist": {f"{a}x{b}": k for (a, b), k in
                          sorted(Counter(r["bbox"] for r in rs).items(),
                                 key=lambda kv: -kv[1])},
            "max_side_hist": dict(sorted(Counter(r["bbox_max_side"]
                                                 for r in rs).items())),
            "n_shapes": len({r["shape"] for r in rs}),
            "touch_both_rows": sum(1 for r in rs
                                   if r["touches"]["top"] and
                                   r["touches"]["bottom"]),
            "touch_both_cols": sum(1 for r in rs
                                   if r["touches"]["left"] and
                                   r["touches"]["right"]),
        }

    # --- FILTRACAO: quanto das 6 direcoes cada classe de tamanho recupera? --
    # A escolha gulosa dos geradores acima e' arbitraria; o que discrimina a
    # hipotese "nao cabe" e' o MENOR tamanho capaz de gerar as 6 direcoes.
    def recovered(pred) -> int:
        B = IntBasis(r["mask"] for r in records if r["mode"] == "realized")
        base = B.rank
        for r in records:
            if r["mode"] != "realized" and pred(r):
                B.add(r["mask"])
        return B.rank - base

    filt = {
        "by_bbox_rows_le": {k: recovered(lambda r, k=k: r["bbox"][0] <= k)
                            for k in range(2, R + 1)},
        "by_bbox_cols_le": {k: recovered(lambda r, k=k: r["bbox"][1] <= k)
                            for k in range(2, C + 1)},
        "by_max_side_le": {k: recovered(lambda r, k=k: r["bbox_max_side"] <= k)
                           for k in range(2, max(R, C) + 1)},
        "interior_only": recovered(lambda r: not any(r["touches"].values())),
        "no_row_border": recovered(lambda r: not (r["touches"]["top"] or
                                                  r["touches"]["bottom"])),
        "no_col_border": recovered(lambda r: not (r["touches"]["left"] or
                                                  r["touches"]["right"])),
        "all_unrealized": recovered(lambda r: True),
    }

    real_rs = [r for r in records if r["mode"] == "realized"]
    res = {
        "R": R, "C": C, "V": V, "E": E, "beta1": E - V + 1,
        "n_tours": len(masks), "n_hexagons_bulk": len(hexes),
        "dim_ZBulk": dim_zb,
        "rank_realized": Breal.rank, "rank_all_hexagons": Ball.rank,
        "codim_missing": codim,
        "dim_span_diffs": Bdiff.rank,
        "span_diffs_eq_span_realized": diffs_eq_real,
        "rank_Ham": Bdiff.rank + 1,
        "deficit": (E - V + 1) - (Bdiff.rank + 1),
        "deficit_minus_3": (E - V + 1) - (Bdiff.rank + 1) - 3,
        "counts": {m: sum(1 for r in records if r["mode"] == m)
                   for m in ("realized", "no_match", "disconnect")},
        "n_nonzero_in_quotient": len(nonzero),
        "filtration": filt,
        "stats_realized": stats(real_rs),
        "stats_unrealized": stats([r for r in records
                                   if r["mode"] != "realized"]),
        "stats_nonzero_in_quotient": stats(nonzero),
        "generators": [
            {"cells": r["cells"], "bbox": list(r["bbox"]),
             "mode": r["mode"], "touches": r["touches"],
             "n_borders": r["n_borders"], "interior": r["interior"],
             "shape": [list(x) for x in r["shape"]]}
            for r in needed
        ],
        "seconds": round(time.time() - t0, 1),
    }
    return res


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rows", type=int, default=5)
    ap.add_argument("--cols", type=int, default=8)
    args = ap.parse_args()
    res = run(args.rows, args.cols)
    out = DATA / f"missing_directions_{args.rows}x{args.cols}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(res, indent=2, default=str))

    print()
    print("=" * 70)
    print(f" DIRECOES FALTANTES  {args.rows}x{args.cols}")
    print("=" * 70)
    for k in ("beta1", "dim_ZBulk", "rank_realized", "rank_all_hexagons",
              "dim_span_diffs", "span_diffs_eq_span_realized", "rank_Ham",
              "deficit", "deficit_minus_3", "codim_missing", "counts",
              "n_nonzero_in_quotient"):
        print(f" {k:30s}: {res[k]}")
    print()
    print(" --- FILTRACAO: direcoes recuperadas (de "
          f"{res['codim_missing']}) por classe de tamanho")
    for kk, vv in res["filtration"].items():
        print(f"     {kk:22s}: {vv}")
    print()
    for k in ("stats_realized", "stats_unrealized", "stats_nonzero_in_quotient"):
        print(f" --- {k}")
        for kk, vv in res[k].items():
            print(f"     {kk:22s}: {vv}")
    print()
    print(f" --- {len(res['generators'])} geradores das direcoes faltantes")
    for i, g in enumerate(res["generators"]):
        t = "".join(s[0].upper() for s in ("top", "bottom", "left", "right")
                    if g["touches"][s]) or "-"
        tag = "INTERIOR" if g["interior"] else "borda   "
        print(f"  [{i}] {tag} bbox={g['bbox']}  modo={g['mode']:10s} "
              f"bordas={t:4s} cells={g['cells']}")
    print(f" -> {out}")


if __name__ == "__main__":
    main()
