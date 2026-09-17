#!/usr/bin/env python3
"""
bulk_perp_dims.py
=================
Fase 3, Tarefa 1: verificacao dimensional dos Lemas A e B.

Notacao (fixada em FASE3_PROMPT.md):
  C^perp      = espaco de cortes = row(d_1),           dim |V| - 1
  Mand        = as 8 arestas incidentes a canto
  X           = C^perp + Span(XOR_pairs),              dim |V| + 2
  Z_bulk      = espaco de ciclos de G menos os 4 cantos
  perp(Z_bulk)= {R em F_2^E : <z,R> = 0 para todo z em Z_bulk}

Lema A:  perp(Z_bulk) = X + <1_e> para qualquer e em Mand  (codim 1 sobre X).

Predicao falseavel:  dim perp(Z_bulk) - dim X = 1  SEMPRE.

Nota: nao e' preciso enumerar tours -- e' algebra linear pura sobre GF(2).
Como toda aresta de Mand incide num canto, e ciclos do bulk evitam cantos,
perp(Z_bulk) = <1_e : e em Mand> + (lift de row(d_1^bulk)).

CLI:
  python bulk_perp_dims.py
  python bulk_perp_dims.py --extra 18,20
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"

MOVES = [(2, 1), (2, -1), (-2, 1), (-2, -1),
         (1, 2), (1, -2), (-1, 2), (-1, -2)]


# ── grafo retangular (generaliza knight_gf2.build_graph) ────────────

def build_rect(rows: int, cols: int):
    """Grafo do cavalo rows x cols. Retorna (adj, edges) com edges ordenado."""
    total = rows * cols
    adj: List[List[int]] = [[] for _ in range(total)]
    eset = set()
    for r in range(rows):
        for c in range(cols):
            v = r * cols + c
            for dr, dc in MOVES:
                nr, nc = r + dr, c + dc
                if 0 <= nr < rows and 0 <= nc < cols:
                    u = nr * cols + nc
                    adj[v].append(u)
                    if u > v:
                        eset.add((v, u))
    return adj, sorted(eset)


def rect_corners(rows: int, cols: int) -> List[int]:
    return [0, cols - 1, (rows - 1) * cols, rows * cols - 1]


# ── rank GF(2) com inteiros (rapido o bastante ate n=16) ────────────

def gf2_rank_ints(vectors: Iterable[int]) -> int:
    """Rank de uma familia de vetores GF(2) codificados como inteiros."""
    pivots: Dict[int, int] = {}          # bit-pivo -> vetor reduzido
    rank = 0
    for v in vectors:
        while v:
            b = v.bit_length() - 1
            if b not in pivots:
                pivots[b] = v
                rank += 1
                break
            v ^= pivots[b]
    return rank


def bit(i: int) -> int:
    return 1 << i


# ── analise de um tabuleiro ─────────────────────────────────────────

def analyse_board(rows: int, cols: int) -> dict:
    adj, edges = build_rect(rows, cols)
    V, E = rows * cols, len(edges)
    eidx = {e: i for i, e in enumerate(edges)}
    cs = rect_corners(rows, cols)

    # arestas obrigatorias: as 2 de cada canto (exige grau 2 => min(n,m) >= 4)
    degs = {c: len(set(adj[c])) for c in cs}
    if any(d != 2 for d in degs.values()):
        return {"rows": rows, "cols": cols, "skipped": f"graus de canto {degs}"}
    mand = []
    for c in cs:
        for u in sorted(set(adj[c])):
            mand.append(eidx[(min(c, u), max(c, u))])
    assert len(set(mand)) == 8, "cantos com arestas compartilhadas"

    # C^perp = row(d_1): uma linha por vertice = 1_{delta(v)}
    cuts = []
    for v in range(V):
        x = 0
        for u in set(adj[v]):
            x ^= bit(eidx[(min(u, v), max(u, v))])
        cuts.append(x)
    dim_cut = gf2_rank_ints(cuts)

    # X = C^perp + Span(XOR_pairs)
    xors = [bit(mand[i]) ^ bit(mand[j])
            for i in range(8) for j in range(i + 1, 8)]
    dim_X = gf2_rank_ints(cuts + xors)
    Q = dim_X - dim_cut

    # perp(Z_bulk): 1_e para e em Mand, mais cortes do grafo-bulk (levantados).
    # Justificativa: Z_bulk vive em F_2^{E_bulk}, logo perp(Z_bulk) em F_2^E e'
    # {tudo em Mand} + {cortes do bulk}.
    corner_set = set(cs)
    cuts_bulk = []
    for v in range(V):
        if v in corner_set:
            continue
        x = 0
        for u in set(adj[v]):
            if u in corner_set:
                continue
            x ^= bit(eidx[(min(u, v), max(u, v))])
        cuts_bulk.append(x)
    perp_basis = cuts_bulk + [bit(m) for m in mand]
    dim_perp = gf2_rank_ints(perp_basis)

    # numero de componentes do bulk (deve ser 1 pelo teorema de conexidade)
    seen = set(cs)
    comps = 0
    for s in range(V):
        if s in seen:
            continue
        comps += 1
        stack = [s]
        seen.add(s)
        while stack:
            x = stack.pop()
            for u in adj[x]:
                if u not in seen:
                    seen.add(u)
                    stack.append(u)

    # Lema A: X esta contido em perp(Z_bulk)?  e a igualdade X + <1_e>?
    X_in_perp = gf2_rank_ints(cuts + xors + perp_basis) == dim_perp
    e0 = mand[0]
    dim_X_plus_e = gf2_rank_ints(cuts + xors + [bit(e0)])
    lemma_A = X_in_perp and dim_X_plus_e == dim_perp and dim_perp - dim_X == 1

    return {
        "rows": rows, "cols": cols, "V": V, "E": E,
        "beta1": E - V + 1,
        "bulk_components": comps,
        "dim_cut": dim_cut, "dim_X": dim_X, "dim_perp_Zbulk": dim_perp,
        "Q": Q,
        "codim_X_in_perp": dim_perp - dim_X,
        "X_subset_perp": X_in_perp,
        "dim_X_plus_1e": dim_X_plus_e,
        "lemma_A_holds": lemma_A,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--extra", type=str, default="",
                    help="tamanhos quadrados adicionais, ex: 18,20")
    a = ap.parse_args()

    squares = [6, 8, 10, 12, 14, 16]
    squares += [int(x) for x in a.extra.split(",") if x.strip()]
    rects = [(6, 8), (6, 10), (8, 10), (6, 7), (7, 9), (8, 12)]

    boards = [(n, n) for n in squares] + rects
    rows_out = []
    print(f"{'board':>8} {'V':>5} {'E':>5} {'comp':>4} {'dimC^p':>7} "
          f"{'dimX':>6} {'dimPerp':>8} {'Q':>3} {'codim':>6} {'LemaA':>6}")
    print("-" * 72)
    t0 = time.time()
    for (r, c) in boards:
        res = analyse_board(r, c)
        rows_out.append(res)
        if "skipped" in res:
            print(f"{r}x{c:<5} PULADO: {res['skipped']}")
            continue
        print(f"{r}x{c:<5} {res['V']:>5} {res['E']:>5} "
              f"{res['bulk_components']:>4} {res['dim_cut']:>7} "
              f"{res['dim_X']:>6} {res['dim_perp_Zbulk']:>8} {res['Q']:>3} "
              f"{res['codim_X_in_perp']:>6} "
              f"{'OK' if res['lemma_A_holds'] else 'FALHA':>6}")

    bad = [r for r in rows_out if not r.get("skipped")
           and not r["lemma_A_holds"]]
    print("-" * 72)
    if bad:
        print(f"!!! LEMA A FALHA em {len(bad)} tabuleiros -- PARE e reporte")
    else:
        print(f"Lema A verificado em {len(rows_out)} tabuleiros; "
              f"codim(X em perp(Z_bulk)) = 1 sempre.  ({time.time() - t0:.1f}s)")

    DATA.mkdir(parents=True, exist_ok=True)
    out = DATA / "bulk_perp_dims.json"
    out.write_text(json.dumps(rows_out, indent=2))
    print(f"-> {out}")


if __name__ == "__main__":
    main()
