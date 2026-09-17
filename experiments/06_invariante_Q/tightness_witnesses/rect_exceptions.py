#!/usr/bin/env python3
"""
rect_exceptions.py
==================
Exceções ao enunciado forte num retângulo R x C qualquer, sem enumerar tours.

Estágio 1 -- bola de BFS por flips, várias sementes. Só produz POSITIVOS, e
cada um é um certificado exibido (o par tau, tau XOR C). Reduz o resíduo.

Estágio 2 -- `flip_paths.flip_realizable` no resíduo. Decisivo nos dois
sentidos; reporta estouro de orçamento quando não conclui.

⚠️ Direção do viés: a bola só afirma "realizável"; a busca por caminhos é a
única que pode afirmar "NÃO realizável", e essa é a afirmação frágil (poda
excessiva fabrica exceções). Por isso `flip_paths` foi validado contra duas
verdades exaustivas independentes -- n=6 (520/532, os dois sentidos) e o 6x7
(4 negativos genuínos). Ver FLIP_GRAPH_RESULTS.md.

CLI:
  python rect_exceptions.py --rows 7 --cols 8
"""
from __future__ import annotations

import argparse
import json
import random
import time
from collections import deque
from pathlib import Path
from typing import Dict, List

from board_symmetry import (board_group, close_certificates,
                            orbit_representatives)
from flip_graph import Basis, alternating_matchings, cycle_mask
from flip_graph_rect import build_rect, hexagons_bulk_rect, is_corner_rect
from flip_paths import flip_realizable
from width6_threshold import ham_rect, tour_mask_rect, warnsdorff_rect

DATA = Path(__file__).resolve().parent / "data"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rows", type=int, required=True)
    ap.add_argument("--cols", type=int, required=True)
    ap.add_argument("--ball", type=int, default=600000)
    ap.add_argument("--seeds", type=int, default=4)
    ap.add_argument("--budget", type=int, default=400_000_000)
    args = ap.parse_args()
    R, C = args.rows, args.cols

    t0 = time.time()
    adj, edges, eidx = build_rect(R, C)
    V = R * C
    hexes = hexagons_bulk_rect(R, C, adj)
    hx = [(cycle_mask(c, eidx), *alternating_matchings(c, eidx)) for c in hexes]
    by_mask = {cycle_mask(c, eidx): c for c in hexes}
    dimZb = (len(edges) - 8) - V + 5
    print(f"[{R}x{C}] V={V} E={len(edges)} beta1={len(edges)-V+1} "
          f"dim_ZBulk={dimZb} hex_bulk={len(hexes)}", flush=True)

    # (1) ANTES da busca: quocientar pelo grupo de automorfismos do tabuleiro.
    # flip-realizabilidade e' constante em cada orbita, entao basta um
    # representante. Fator ~4 em retangulo, ~8 em quadrado.
    reps, orbit_of = orbit_representatives(R, C, hexes, edges, eidx)
    G = board_group(R, C)
    print(f"[simetria] |G|={len(G)}  {len(hexes)} hexagonos -> {len(reps)} "
          f"representantes (fator {len(hexes)/len(reps):.2f})", flush=True)

    realized: Dict[int, tuple] = {}
    for s in range(args.seeds):
        rng = random.Random(500 + s)
        p = warnsdorff_rect(R, C, adj, rng)
        if p is None:
            continue
        root = tour_mask_rect(p, eidx)
        seen = {root}
        q = deque([root])
        while q and len(seen) < args.ball:
            m = q.popleft()
            for cm, m0, m1 in hx:
                it = m & cm
                if it != m0 and it != m1:
                    continue
                nm = m ^ cm
                if not ham_rect(nm, V, edges):
                    continue
                if cm not in realized:
                    realized[cm] = (m, nm)
                if nm not in seen:
                    seen.add(nm)
                    q.append(nm)
        print(f"  semente {s}: bola={len(seen)} realizados="
              f"{len(realized)}/{len(hexes)} ({time.time()-t0:.0f}s)", flush=True)
        if len(realized) == len(hexes):
            break

    # (2) DEPOIS da bola: fechar os certificados sob o grupo, ANTES de gastar
    # busca. A bola e' assimetrica, entao o residuo dela quase sempre contem
    # hexagonos cujo espelho ja' esta' certificado. No 7x8 isso resolveu em
    # milissegundos um caso que consumiu 4.781 s de busca e voltou
    # inconclusivo. Cada transporte e' VERIFICADO, nunca assumido.
    before = len(realized)
    hex_masks = set(by_mask)
    realized = close_certificates(R, C, realized, hex_masks, edges, eidx,
                                  lambda mk: ham_rect(mk, V, edges))
    print(f"[simetria] fecho do certificado: {before} -> {len(realized)} "
          f"(+{len(realized)-before} sem busca)", flush=True)

    B = Basis()
    for cm in realized:
        B.add(cm)
    # so' representantes de orbita ainda nao certificados vao para a busca
    residue = [by_mask[cm] for cm in by_mask
               if cm not in realized and orbit_of[cm] == cm]
    dropped = sum(1 for cm in by_mask
                  if cm not in realized and orbit_of[cm] != cm)
    if dropped:
        print(f"[simetria] {dropped} nao certificados sao imagens de "
              f"representantes -- herdam a resposta", flush=True)
    print(f"[estagio 2] residuo={len(residue)}  rank realizados="
          f"{B.rank}/{dimZb}", flush=True)

    rows = []
    for k, cyc in enumerate(residue):
        cells = sorted(divmod(v, C) for v in cyc)
        anch = (any(r in (0, R - 1) for r, _ in cells)
                and any(c in (0, C - 1) for _, c in cells))
        t1 = time.time()
        f, ex = flip_realizable(V, adj, edges, eidx, cyc, args.budget)
        rows.append({"cells": cells, "canto": anch, "flip": f,
                     "exaustivo": ex, "s": round(time.time() - t1, 1)})
        print(f"  {k+1}/{len(residue)} {cells} canto={anch} -> flip={f} "
              f"exaustivo={ex} ({time.time()-t1:.0f}s)", flush=True)
        json.dump(rows, open(DATA / f"rect_exceptions_{R}x{C}.json", "w"),
                  indent=2)

    # expandir de volta as orbitas: cada representante provado vale por toda
    # a sua orbita. Sem isso o sumario sub-reporta.
    from collections import Counter
    orbit_size = Counter(orbit_of.values())
    for x, cyc in zip(rows, residue):
        x["orbit_size"] = orbit_size[cycle_mask(cyc, eidx)]
    prov = [x for x in rows if not x["flip"] and x["exaustivo"]]
    prov_total = sum(x["orbit_size"] for x in prov)
    prov_corner = sum(x["orbit_size"] for x in prov if x["canto"])
    real_total = sum(x["orbit_size"] for x in rows if x["flip"])
    inc_total = sum(x["orbit_size"] for x in rows
                    if not x["flip"] and not x["exaustivo"])
    print()
    print("=" * 66)
    print(f" EXCECOES EM {R}x{C}   (min lado = {min(R,C)})")
    print("=" * 66)
    print(f" hexagonos do bulk           : {len(hexes)}")
    print(f" certificados pela bola      : {len(realized)}")
    print(f" representantes buscados     : {len(residue)}")
    print(f" NAO realizaveis PROVADOS    : {prov_total} hexagonos"
          f"  ({len(prov)} orbitas, em canto: {prov_corner})")
    print(f" realizaveis pela busca      : {real_total}"
          f"  ({sum(1 for x in rows if x['flip'])} orbitas)")
    print(f" inconclusivos (orcamento)   : {inc_total}"
          f"  ({sum(1 for x in rows if not x['flip'] and not x['exaustivo'])} orbitas)")
    print(f" rank realizados / dim Z_bulk: {B.rank}/{dimZb}")
    print(f" tempo {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
