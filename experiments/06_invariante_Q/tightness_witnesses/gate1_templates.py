#!/usr/bin/env python3
"""
gate1_templates.py
==================
Gate 1 do PLAN.md: qual a MENOR familia geradora de hexagonos fechada por
translacao?

Motivacao. Nao e' preciso certificar todo hexagono -- basta que ALGUMA familia
geradora seja certificada. Se essa familia for uniao de poucas CLASSES DE
TRANSLACAO ("templates"), entao (U1) colapsa de "teorema sobre todos os
hexagonos" para "verificacao finita de k templates + argumento de translacao".

Procedimento:
  1. enumerar hexagonos de m x m, classifica-los por forma (translacao)
  2. checar que o numero de formas ESTABILIZA em m (se nao estabiliza, nao ha
     template finito e o Gate 1 falha)
  3. busca gulosa/exaustiva pela menor colecao de formas cujas translacoes
     geram Z_1(G_m)
  4. reportar tambem a classificacao por translacao+D4

CLI:
  python gate1_templates.py --ms 8,10,12,14
"""
from __future__ import annotations

import argparse
import json
import time
from itertools import combinations
from pathlib import Path
from typing import Dict, List, Tuple

from hexagon_locality import Basis, build_sub, cyc_vec, cycles_len

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"


def shape_translation(cells: List[Tuple[int, int]]) -> Tuple:
    """Forma a menos de translacao: ordena e desloca para min = (0,0)."""
    s = sorted(cells)
    r0, c0 = s[0][0], min(c for _, c in s)
    return tuple((r - r0, c - c0) for r, c in s)


D4 = [lambda r, c: (r, c), lambda r, c: (c, -r), lambda r, c: (-r, -c),
      lambda r, c: (-c, r), lambda r, c: (r, -c), lambda r, c: (-r, c),
      lambda r, c: (c, r), lambda r, c: (-c, -r)]


def shape_d4(cells: List[Tuple[int, int]]) -> Tuple:
    """Forma a menos de translacao + D4: minimo lexicografico das 8 imagens."""
    return min(shape_translation([g(r, c) for r, c in cells]) for g in D4)


def collect(m: int):
    cells = [(r, c) for r in range(m) for c in range(m)]
    idx, adj, edges = build_sub(cells)
    rev = {i: c for c, i in idx.items()}
    eidx = {e: i for i, e in enumerate(edges)}
    V, E = len(idx), len(edges)
    dimZ = E - V + 1

    hexs = cycles_len(adj, V, 6)
    by_t: Dict[Tuple, List[int]] = {}
    by_d4: Dict[Tuple, List[Tuple]] = {}
    vecs: Dict[Tuple, List[int]] = {}
    for h in hexs:
        cs = [rev[v] for v in h]
        st = shape_translation(cs)
        by_t.setdefault(st, []).append(0)
        vecs.setdefault(st, []).append(cyc_vec(h, eidx))
        by_d4.setdefault(shape_d4(cs), []).append(st)
    return dimZ, vecs, {k: len(set(v)) for k, v in by_d4.items()}


def minimal_family(dimZ: int, vecs: Dict[Tuple, List[int]], max_k: int):
    """Menor colecao de formas cujas translacoes geram Z_1."""
    shapes = sorted(vecs, key=lambda s: -len(vecs[s]))

    # rank individual de cada forma
    solo = {}
    for s in shapes:
        b = Basis()
        for x in vecs[s]:
            b.add(x)
        solo[s] = b.rank

    # guloso: escolhe a forma que mais aumenta o rank
    b = Basis()
    chosen: List[Tuple] = []
    while b.rank < dimZ and len(chosen) < len(shapes):
        best, best_gain, best_basis = None, -1, None
        for s in shapes:
            if s in chosen:
                continue
            t = Basis()
            t.piv = dict(b.piv)
            for x in vecs[s]:
                t.add(x)
            if t.rank - b.rank > best_gain:
                best, best_gain, best_basis = s, t.rank - b.rank, t
        if best_gain <= 0:
            break
        chosen.append(best)
        b = best_basis
    greedy_k, greedy_rank = len(chosen), b.rank

    # exaustivo para k pequeno (confirma minimalidade do guloso)
    exact_k = None
    if greedy_rank == dimZ:
        for k in range(1, min(greedy_k, max_k) + 1):
            found = False
            for combo in combinations(shapes, k):
                t = Basis()
                for s in combo:
                    for x in vecs[s]:
                        t.add(x)
                if t.rank == dimZ:
                    exact_k = (k, combo)
                    found = True
                    break
            if found:
                break
    return solo, chosen, greedy_rank, exact_k


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ms", type=str, default="8,10,12")
    ap.add_argument("--max-k", type=int, default=3)
    a = ap.parse_args()

    print("=" * 72)
    print(" Gate 1 — templates de translacao de hexagonos")
    print("=" * 72)
    out = []
    for m in [int(x) for x in a.ms.split(",")]:
        t0 = time.time()
        dimZ, vecs, d4 = collect(m)
        n_hex = sum(len(v) for v in vecs.values())
        print(f"\n m={m}: hexagonos={n_hex}  dim Z_1={dimZ}  "
              f"formas(transl)={len(vecs)}  formas(transl+D4)={len(d4)}")

        solo, chosen, grank, exact = minimal_family(dimZ, vecs, a.max_k)
        top = sorted(solo.items(), key=lambda kv: -kv[1])[:6]
        print("   rank das translacoes de UMA forma (top 6): "
              + ", ".join(f"{len(vecs[s])}x->{r}" for s, r in top))
        print(f"   guloso: {len(chosen)} formas -> rank {grank}/{dimZ} "
              f"{'GERA' if grank == dimZ else 'NAO GERA'}")
        if exact:
            k, combo = exact
            print(f"   MINIMO EXATO: k={k} formas bastam")
            for s in combo:
                print(f"      {s}   ({len(vecs[s])} translacoes)")
        out.append({
            "m": m, "n_hex": n_hex, "dim_Z1": dimZ,
            "n_shapes_translation": len(vecs), "n_shapes_d4": len(d4),
            "greedy_k": len(chosen), "greedy_rank": grank,
            "greedy_spans": grank == dimZ,
            "minimal_k": exact[0] if exact else None,
            "minimal_shapes": [list(map(list, s)) for s in exact[1]] if exact else None,
            "shape_solo_rank": {str(s): r for s, r in sorted(
                solo.items(), key=lambda kv: -kv[1])},
            "elapsed_s": round(time.time() - t0, 1),
        })

    DATA.mkdir(parents=True, exist_ok=True)
    p = DATA / "gate1_templates.json"
    p.write_text(json.dumps(out, indent=2))
    print(f"\n-> {p}")


if __name__ == "__main__":
    main()
