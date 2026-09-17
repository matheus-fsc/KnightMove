#!/usr/bin/env python3
"""
hexagon_in_tourspace.py
=======================
Auditoria de ponta a ponta da afirmacao critica da Fase 3.

MOTIVACAO. A cadeia da Fase 3 depende de "os hexagonos bons geram Z_bulk".
Duas objecoes possiveis:

  (i)  amostragem -- a varredura viu so' parte dos hexagonos;
  (ii) falso positivo -- um "hexagono bom" que na verdade nao e' bom.

Sobre (i): a afirmacao e' MONOTONA na direcao segura. Cada hexagono bom entra
com certificado positivo, e rank so' cresce ao acrescentar vetores. Amostrar
pode apenas SUB-reportar. Ainda assim, aqui fazemos ENUMERACAO COMPLETA.

Sobre (ii): e' o risco real. Aqui a verificacao e' end-to-end e independente
do `verify_switcher` do switcher_search.py:

  H'   = caminho hamiltoniano de v1 a v4 em G \\ (V(W)\\{v1,v4})   [passo (S3)]
  H_A  = H' concatenado com o 1o caminho hamiltoniano de W
  H_B  = H' concatenado com o 2o caminho hamiltoniano de W

Verificamos do ZERO que H_A e H_B sao CICLOS HAMILTONIANOS de G_n (todo
vertice grau 2, conexo, cobre V, todas as arestas existem) e que

                       H_A  XOR  H_B  =  C.

Isso e' mais forte do que "o esquema do CNP e' instanciavel": exibe dois tours
do cavalo cuja diferenca simetrica e' exatamente o hexagono, logo prova
CONSTRUTIVAMENTE que C pertence a Span(Ham(n)).

CONSEQUENCIA. Se os hexagonos assim certificados geram Z_bulk, entao
Z_bulk esta contido em Span(Ham(n)), donde
      rank(Ham) >= dim Z_bulk = beta_1 - 4.
Note que a tightness pede beta_1 - 3: sobra exatamente UMA dimensao. Esse
gap de 1 e' a forma primal da codimensao 1 do Lema A.

CLI:
  python hexagon_in_tourspace.py --ns 8,10,12          # enumeracao completa
  python hexagon_in_tourspace.py --ns 8 --max-hex 500  # amostra (mais rapido)
"""
from __future__ import annotations

import argparse
import json
import random
import time
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from bulk_perp_dims import build_rect, gf2_rank_ints, rect_corners
from switcher_search import (build, build_switchers, color, ham_path,
                             hexagons, label)

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"


# ── auditoria independente ──────────────────────────────────────────

def edges_of_path(seq: List[int]) -> Set[frozenset]:
    return {frozenset((a, b)) for a, b in zip(seq, seq[1:])}


def audit(adj, n: int, h: Tuple[int, ...], w: dict, hp: List[int],
          corner_set: Set[int]) -> Tuple[Optional[Set[frozenset]], List[str]]:
    """
    Reconstroi H_A e H_B do zero e confere que sao ciclos hamiltonianos de G_n
    com H_A xor H_B = C. Retorna (arestas do hexagono, problemas).
    """
    bad: List[str] = []
    total = n * n
    v = w["cycle"]
    i2, i3 = w["P2_internal"], w["P3_internal"]

    # --- o hexagono e' mesmo um 6-ciclo do bulk? ---
    if len(set(h)) != 6:
        bad.append("hexagono com vertices repetidos")
    if set(h) & corner_set:
        bad.append("hexagono toca canto (nao esta no bulk)")
    for a, b in zip(h, h[1:] + h[:1]):
        if b not in adj[a]:
            bad.append(f"hexagono: {a}-{b} nao e' aresta de G_n")
    hex_edges = {frozenset((a, b)) for a, b in zip(h, h[1:] + h[:1])}
    if len(hex_edges) != 6:
        bad.append("hexagono nao tem 6 arestas distintas")
    # C como vetor deve estar no espaco de ciclos: todo vertice com grau par
    degc: Dict[int, int] = {}
    for e in hex_edges:
        for x in e:
            degc[x] = degc.get(x, 0) + 1
    if any(d != 2 for d in degc.values()):
        bad.append("hexagono nao e' 2-regular (nao e' ciclo simples)")
    # o conjunto de arestas de C deve coincidir com o de w["cycle"]
    if hex_edges != {frozenset((a, b)) for a, b in zip(v, v[1:] + v[:1])}:
        bad.append("w['cycle'] nao e' o mesmo hexagono")

    # --- os dois caminhos hamiltonianos de W ---
    pathA = [v[0], v[1]] + i2 + [v[5], v[4]] + i3[::-1] + [v[2], v[3]]
    pathB = [v[0], v[5]] + i2[::-1] + [v[1], v[2]] + i3 + [v[4], v[3]]
    VW = set(v) | set(i2) | set(i3)
    for name, p in (("W-hamA", pathA), ("W-hamB", pathB)):
        if len(p) != len(set(p)) or set(p) != VW:
            bad.append(f"{name}: nao cobre V(W) exatamente uma vez")
        for a, b in zip(p, p[1:]):
            if b not in adj[a]:
                bad.append(f"{name}: {a}-{b} nao e' aresta")

    # --- H' cobre o complemento e so' toca W em v1, v4 ---
    if hp[0] != v[0] or hp[-1] != v[3]:
        bad.append("(S3): extremos errados")
    if len(hp) != len(set(hp)):
        bad.append("(S3): vertice repetido")
    if set(hp) != (set(range(total)) - (VW - {v[0], v[3]})):
        bad.append("(S3): nao cobre exatamente o complemento de W")
    for a, b in zip(hp, hp[1:]):
        if b not in adj[a]:
            bad.append(f"(S3): {a}-{b} nao e' aresta")

    if bad:
        return None, bad

    # --- H_A e H_B como CICLOS HAMILTONIANOS de G_n ---
    eh = edges_of_path(hp)
    cycles = {}
    for name, p in (("H_A", pathA), ("H_B", pathB)):
        cyc = eh | edges_of_path(p)
        # todo vertice grau 2
        deg: Dict[int, int] = {}
        for e in cyc:
            for x in e:
                deg[x] = deg.get(x, 0) + 1
        if len(deg) != total or any(d != 2 for d in deg.values()):
            bad.append(f"{name}: nao e' 2-regular sobre todos os {total} vertices")
            continue
        if len(cyc) != total:
            bad.append(f"{name}: {len(cyc)} arestas, esperado {total}")
            continue
        # conexo (um unico ciclo, nao uniao de ciclos)
        nb: Dict[int, List[int]] = {x: [] for x in deg}
        for e in cyc:
            a, b = tuple(e)
            nb[a].append(b)
            nb[b].append(a)
        seen = {0}
        stack = [0]
        while stack:
            x = stack.pop()
            for u in nb[x]:
                if u not in seen:
                    seen.add(u)
                    stack.append(u)
        if len(seen) != total:
            bad.append(f"{name}: desconexo (uniao de sub-ciclos), "
                       f"componente tem {len(seen)}/{total}")
            continue
        for e in cyc:
            a, b = tuple(e)
            if b not in adj[a]:
                bad.append(f"{name}: aresta inexistente")
                break
        cycles[name] = cyc

    if len(cycles) == 2:
        sym = cycles["H_A"] ^ cycles["H_B"]
        if sym != hex_edges:
            bad.append(f"H_A xor H_B tem {len(sym)} arestas, "
                       f"nao e' o hexagono (6)")

    return (hex_edges if not bad else None), bad


# ── varredura ───────────────────────────────────────────────────────

def run(n: int, max_hex: Optional[int], max_internal: int, max_w: int,
        budget: int, seed: int) -> dict:
    t0 = time.time()
    adj = build(n)
    _, edges = build_rect(n, n)
    eidx = {e: i for i, e in enumerate(edges)}
    total = n * n
    cs = set(rect_corners(n, n))

    V_bulk, E_bulk = total - 4, len(edges) - 8
    dim_Zbulk = E_bulk - V_bulk + 1
    beta1 = len(edges) - total + 1

    hexa = hexagons(adj, cs, total)
    n_all = len(hexa)
    complete = max_hex is None or max_hex >= n_all
    if not complete:
        random.Random(seed).shuffle(hexa)
        hexa = hexa[:max_hex]

    print(f"\n G_{n}x{n}:  beta_1 = {beta1},  dim Z_bulk = {dim_Zbulk} "
          f"(= beta_1 - 4)")
    print(f"   hexagonos no bulk: {n_all}  -- varrendo "
          f"{'TODOS (enumeracao completa)' if complete else f'{len(hexa)} (amostra)'}")

    good_vecs: List[int] = []
    n_good = n_bad = 0
    audit_failures: List[str] = []

    for h in hexa:
        cands: List[dict] = []
        for rot in range(3):
            cands += build_switchers(adj, h, rot, max_internal, max_w)
        cands.sort(key=lambda w: w["VW"])

        certified = False
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
            p, _ex, _nd = ham_path(adj, allowed, v1, v4, budget)
            if p is None:
                continue
            hex_edges, problems = audit(adj, n, h, w, p, cs)
            if problems:
                audit_failures += problems[:2]
                continue
            certified = True
            x = 0
            for e in hex_edges:
                a, b = tuple(e)
                x ^= 1 << eidx[(min(a, b), max(a, b))]
            good_vecs.append(x)
            break

        if certified:
            n_good += 1
        else:
            n_bad += 1

    rank_good = gf2_rank_ints(good_vecs)
    spans = rank_good == dim_Zbulk

    print(f"   hexagonos CERTIFICADOS (2 tours com XOR = C): {n_good}")
    print(f"   hexagonos sem certificado                   : {n_bad}")
    print(f"   falhas de auditoria                         : "
          f"{len(audit_failures)}")
    if audit_failures:
        print(f"   !!! {audit_failures[:4]}")
    print(f"   rank(span dos certificados) = {rank_good}/{dim_Zbulk}  "
          f"{'GERA Z_bulk' if spans else 'NAO gera'}")
    if spans:
        print(f"   => Z_bulk contido em Span(Ham), logo rank(Ham) >= "
              f"{dim_Zbulk} = beta_1 - 4  (tightness pede beta_1 - 3: "
              f"falta 1 dimensao)")

    return {
        "n": n, "beta1": beta1, "dim_Z_bulk": dim_Zbulk,
        "hexagons_total_bulk": n_all,
        "hexagons_scanned": len(hexa),
        "complete_enumeration": complete,
        "certified": n_good, "uncertified": n_bad,
        "audit_failures": len(audit_failures),
        "rank_certified_span": rank_good,
        "spans_Z_bulk": spans,
        "implies_rank_Ham_at_least": dim_Zbulk if spans else None,
        "gap_to_tightness": 1 if spans else None,
        "params": {"max_hex": max_hex, "max_internal": max_internal,
                   "max_w": max_w, "budget": budget, "seed": seed},
        "elapsed_s": round(time.time() - t0, 1),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ns", type=str, default="8,10,12")
    ap.add_argument("--max-hex", type=int, default=None,
                    help="omitido = enumeracao COMPLETA")
    ap.add_argument("--max-internal", type=int, default=4)
    ap.add_argument("--max-w", type=int, default=20)
    ap.add_argument("--budget", type=int, default=400_000)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()

    print("=" * 72)
    print(" Auditoria end-to-end: cada hexagono bom = 2 tours com XOR = C")
    print("=" * 72)
    out = [run(n, a.max_hex, a.max_internal, a.max_w, a.budget, a.seed)
           for n in [int(x) for x in a.ns.split(",")]]

    print("\n" + "=" * 72)
    for r in out:
        print(f"   n={r['n']:>3}: {'COMPLETA' if r['complete_enumeration'] else 'amostra'}"
              f"  certificados={r['certified']}/{r['hexagons_scanned']}"
              f"  rank={r['rank_certified_span']}/{r['dim_Z_bulk']}"
              f"  {'GERA' if r['spans_Z_bulk'] else 'NAO GERA'}"
              f"  auditoria_falhas={r['audit_failures']}")

    DATA.mkdir(parents=True, exist_ok=True)
    p = DATA / "hexagon_in_tourspace.json"
    p.write_text(json.dumps(out, indent=2))
    print(f"\n-> {p}")


if __name__ == "__main__":
    main()
