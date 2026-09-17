#!/usr/bin/env python3
"""
b6_grouping_6x6.py
==================
Particiona as 9.862 soluções únicas do 6×6 (CICLO Hamiltoniano fechado) pelo
"estado do vértice B6" — qual par dos 3 saltos legais é usado em cada solução.

Convenção de label (de analise_6x6.py):
  label(v) = chr(65 + col) + str(BOARD - row), com BOARD = 6

Saltos de cavalo a partir de B6 = (0, 1):
  D5  = (1, 3)   interior
  C4  = (2, 2)   interior
  A4  = (2, 0)   borda esquerda

Estados possíveis (B6 é vértice interior do ciclo, com grau 2 no caminho):
  S0 = {D5, C4}  ambos interior
  S1 = {D5, A4}  um interior + uma borda
  S2 = {C4, A4}  um interior + uma borda

(O usuário escreveu "A5" no prompt — corrigido para A4, que é o salto real
de cavalo a partir de B6 nesta convenção.)

Fonte dos dados:
  destruction_catalogue.json contém 710_064 sequências direcionadas
  (= 9_862 ciclos × 36 rotações × 2 sentidos). Dedup por frozenset de arestas.

Saídas:
  1. Contagem por estado de B6
  2. Para cada grupo: top arestas com freq > 0.85 (quase obrigatórias) e
     < 0.15 (quase impossíveis)
  3. Arestas com maior |Δfreq| entre grupos
  4. Veredito sobre propagação da invariante de borda para o interior
"""

import json
import os
import sys
from collections import defaultdict, Counter

import numpy as np


BOARD = 6
TOTAL = BOARD * BOARD
MOVES = [(2, 1), (2, -1), (-2, 1), (-2, -1),
         (1, 2), (1, -2), (-1, 2), (-1, -2)]


def vid(r, c):
    return r * BOARD + c


def vrc(v):
    return divmod(v, BOARD)


def label(v):
    r, c = vrc(v)
    return chr(65 + c) + str(BOARD - r)


def build_edges():
    EDGES = set()
    for r in range(BOARD):
        for c in range(BOARD):
            v = vid(r, c)
            for dr, dc in MOVES:
                nr, nc = r + dr, c + dc
                if 0 <= nr < BOARD and 0 <= nc < BOARD:
                    u = vid(nr, nc)
                    e = (min(u, v), max(u, v))
                    EDGES.add(e)
    return sorted(EDGES)


def is_interior_vertex(v):
    """v=(r,c) é interior se 1≤r≤4 e 1≤c≤4 no 6×6."""
    r, c = vrc(v)
    return 1 <= r <= BOARD - 2 and 1 <= c <= BOARD - 2


def is_interior_edge(edge):
    a, b = edge
    return is_interior_vertex(a) and is_interior_vertex(b)


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    cat_path = os.path.join(here, "destruction_catalogue.json")

    print(f"Carregando destruction_catalogue.json …")
    with open(cat_path, "r", encoding="utf-8") as f:
        cat = json.load(f)
    print(f"  total de sequências direcionadas: {len(cat):,}")
    print(f"  esperado 9_862 × 36 × 2 = 710_064:  "
          f"{'OK' if len(cat) == 710_064 else 'INESPERADO'}")

    # ── edge index ────────────────────────────────────────────────────
    EDGES = build_edges()
    n_edges = len(EDGES)
    edge_to_idx = {e: i for i, e in enumerate(EDGES)}
    label_to_vid = {label(v): v for v in range(TOTAL)}
    label_to_edge = {f"{label(e[0])}-{label(e[1])}": e for e in EDGES}
    print(f"  arestas do grafo                 : {n_edges}")

    interior_edge_mask = np.array(
        [is_interior_edge(e) for e in EDGES], dtype=bool)
    print(f"  arestas com AMBOS endpoints interiores: "
          f"{int(interior_edge_mask.sum())}")

    # ── B6 e seus saltos ─────────────────────────────────────────────
    B6 = label_to_vid["B6"]
    D5 = label_to_vid["D5"]
    C4 = label_to_vid["C4"]
    A4 = label_to_vid["A4"]
    nbrs_B6 = [D5, C4, A4]
    nbr_labels = ["D5", "C4", "A4"]
    print(f"\n  B6 = {B6} (rc={vrc(B6)})")
    print(f"  Saltos de cavalo (vid, label):")
    for n, lab in zip(nbrs_B6, nbr_labels):
        e = (min(B6, n), max(B6, n))
        print(f"    {n:>3}  {lab}  edge={label(e[0])}-{label(e[1])}  "
              f"idx={edge_to_idx[e]}")

    e_B6_D5 = (min(B6, D5), max(B6, D5))
    e_B6_C4 = (min(B6, C4), max(B6, C4))
    e_B6_A4 = (min(B6, A4), max(B6, A4))
    i_D5 = edge_to_idx[e_B6_D5]
    i_C4 = edge_to_idx[e_B6_C4]
    i_A4 = edge_to_idx[e_B6_A4]

    # ── dedup ciclos (frozenset de arestas) ──────────────────────────
    print(f"\nDeduplicando ciclos por frozenset de arestas…")
    seen = {}
    for entry in cat:
        seq_labels = entry["sequence"]
        seq_v = [label_to_vid[s] for s in seq_labels]
        es = set()
        for i in range(len(seq_v)):
            a, b = seq_v[i], seq_v[(i + 1) % len(seq_v)]
            es.add((min(a, b), max(a, b)))
        fs = frozenset(es)
        if fs not in seen:
            seen[fs] = es
    n_unique = len(seen)
    print(f"  ciclos únicos: {n_unique:,}  (esperado 9_862: "
          f"{'OK' if n_unique == 9_862 else 'INESPERADO'})")

    # ── monta matriz de assinaturas + estado de B6 ───────────────────
    print(f"\nMontando matriz de assinaturas (n_unique × {n_edges})…")
    X = np.zeros((n_unique, n_edges), dtype=np.bool_)
    states = np.zeros(n_unique, dtype=np.int8)  # 0,1,2

    state_lookup = {
        (True, True, False):  0,    # {D5, C4}
        (True, False, True):  1,    # {D5, A4}
        (False, True, True):  2,    # {C4, A4}
    }

    bad = 0
    for k, es in enumerate(seen.values()):
        for e in es:
            X[k, edge_to_idx[e]] = True
        b = (X[k, i_D5], X[k, i_C4], X[k, i_A4])
        st = state_lookup.get(b)
        if st is None:
            bad += 1
            states[k] = -1
        else:
            states[k] = st
    if bad:
        print(f"  AVISO: {bad} ciclos com estado de B6 fora de {{S0,S1,S2}}")

    # ── grupos ──────────────────────────────────────────────────────
    groups = {0: np.where(states == 0)[0],
              1: np.where(states == 1)[0],
              2: np.where(states == 2)[0]}
    print(f"\n=== Grupos por estado de B6 ===")
    state_desc = {
        0: "S0 = {D5, C4}   (ambas interiores)",
        1: "S1 = {D5, A4}   (interior + borda)",
        2: "S2 = {C4, A4}   (interior + borda)",
    }
    for s in (0, 1, 2):
        idxs = groups[s]
        print(f"  {state_desc[s]}: {len(idxs):>5,} ciclos "
              f"({len(idxs)/n_unique:6.2%})")

    # ── frequência por grupo ─────────────────────────────────────────
    freqs = {}
    for s in (0, 1, 2):
        idxs = groups[s]
        freqs[s] = X[idxs].mean(axis=0) if len(idxs) > 0 \
            else np.zeros(n_edges)

    # ── arestas com freq extrema por grupo ──────────────────────────
    print(f"\n=== Arestas com frequência extrema dentro de cada grupo ===")
    print(f"    (limiares: > 0.85 = quase obrigatória, < 0.15 = quase ausente)")
    for s in (0, 1, 2):
        f = freqs[s]
        extremes = []
        for i, fi in enumerate(f):
            if fi > 0.85 or fi < 0.15:
                extremes.append((fi, i))
        extremes.sort(key=lambda x: x[0])
        print(f"\n  Grupo {s} ({state_desc[s].split('=')[1].strip()}):"
              f"  {len(extremes)}/{n_edges} arestas extremas")
        # imprime em duas seções: ≈0 e ≈1
        near0 = [x for x in extremes if x[0] < 0.15]
        near1 = [x for x in extremes if x[0] > 0.85]
        if near0:
            print(f"    quase ausentes ({len(near0)}):")
            for f, i in near0:
                lab = f"{label(EDGES[i][0])}-{label(EDGES[i][1])}"
                tag = "[int]" if interior_edge_mask[i] else "     "
                print(f"      {tag} {lab:>9s} idx={i:>3}  freq={f:.3f}")
        if near1:
            print(f"    quase obrigatórias ({len(near1)}):")
            for f, i in near1:
                lab = f"{label(EDGES[i][0])}-{label(EDGES[i][1])}"
                tag = "[int]" if interior_edge_mask[i] else "     "
                print(f"      {tag} {lab:>9s} idx={i:>3}  freq={f:.3f}")

    # ── arestas com maior delta entre grupos ─────────────────────────
    print(f"\n=== Arestas com maior variação de frequência entre grupos ===")
    f0, f1, f2 = freqs[0], freqs[1], freqs[2]
    fstack = np.stack([f0, f1, f2], axis=1)
    delta = fstack.max(axis=1) - fstack.min(axis=1)
    order = np.argsort(-delta)
    n_show = 25
    print(f"\n  top-{n_show} arestas por Δmax-min de freq entre os 3 grupos:")
    print(f"    {'edge':>9s}  {'int?':<5s}  "
          f"{'S0={D5,C4}':>10s}  {'S1={D5,A4}':>10s}  {'S2={C4,A4}':>10s}  "
          f"{'Δmax':>6s}")
    for k in order[:n_show]:
        if delta[k] < 0.05:
            break
        lab = f"{label(EDGES[k][0])}-{label(EDGES[k][1])}"
        tag = "[int]" if interior_edge_mask[k] else "     "
        print(f"    {lab:>9s}  {tag:<5s}  "
              f"{f0[k]:>10.3f}  {f1[k]:>10.3f}  {f2[k]:>10.3f}  "
              f"{delta[k]:>6.3f}")

    # ── arestas SOMENTE INTERIORES com maior delta ──────────────────
    delta_int = delta.copy()
    delta_int[~interior_edge_mask] = 0
    order_int = np.argsort(-delta_int)
    n_show_int = 15
    print(f"\n  top-{n_show_int} arestas INTERIORES (ambos endpoints interior) "
          f"por Δmax-min:")
    print(f"    {'edge':>9s}  "
          f"{'S0={D5,C4}':>10s}  {'S1={D5,A4}':>10s}  {'S2={C4,A4}':>10s}  "
          f"{'Δmax':>6s}")
    for k in order_int[:n_show_int]:
        if delta_int[k] < 0.05:
            break
        lab = f"{label(EDGES[k][0])}-{label(EDGES[k][1])}"
        print(f"    {lab:>9s}  "
              f"{f0[k]:>10.3f}  {f1[k]:>10.3f}  {f2[k]:>10.3f}  "
              f"{delta_int[k]:>6.3f}")

    # ── estatísticas agregadas: distância L1 entre distribuições ─────
    print(f"\n=== Distâncias entre distribuições de freq dos 3 grupos ===")
    pairs = [(0, 1), (0, 2), (1, 2)]
    print(f"    {'par':<10s}  {'Σ|Δ| total':>10s}  {'Σ|Δ| interior':>15s}  "
          f"{'#arest Δ>0.20':>15s}")
    for i, j in pairs:
        diff = np.abs(freqs[i] - freqs[j])
        s_total = diff.sum()
        s_int = diff[interior_edge_mask].sum()
        n_big = int((diff > 0.20).sum())
        print(f"    S{i}↔S{j}     {s_total:>10.3f}  {s_int:>15.3f}  "
              f"{n_big:>15d}")

    # ── veredito ─────────────────────────────────────────────────────
    print(f"\n=== VEREDITO ===")
    max_delta = float(delta.max())
    max_delta_int = float(delta_int.max())
    n_extreme_interior_edges = int(np.sum(
        (delta_int > 0.30) & interior_edge_mask))
    print(f"  Δmax sobre todas arestas        : {max_delta:.3f}")
    print(f"  Δmax sobre arestas interiores    : {max_delta_int:.3f}")
    print(f"  Arestas interiores com Δ > 0.30 : {n_extreme_interior_edges}")
    if max_delta_int < 0.10:
        print(f"  → INTERIOR INDEPENDENTE: o estado de B6 não organiza o "
              f"interior; invariante é local")
    elif n_extreme_interior_edges == 0 and max_delta_int < 0.30:
        print(f"  → ORGANIZAÇÃO FRACA: o estado de B6 influencia o interior "
              f"de forma fraca/difusa")
    else:
        print(f"  → PROPAGAÇÃO ESTRUTURAL: o estado de B6 organiza "
              f"significativamente o interior; "
              f"{n_extreme_interior_edges} arestas interiores trocam de "
              f"regime entre grupos")

    # ── salvar resultados ────────────────────────────────────────────
    out_path = os.path.join(here, "b6_grouping_6x6.json")
    out = {
        "board": BOARD,
        "n_unique_cycles": int(n_unique),
        "edges": [f"{label(e[0])}-{label(e[1])}" for e in EDGES],
        "is_interior_edge": interior_edge_mask.tolist(),
        "groups_count": {
            "S0_D5_C4": int(len(groups[0])),
            "S1_D5_A4": int(len(groups[1])),
            "S2_C4_A4": int(len(groups[2])),
        },
        "freq_by_state": {
            "S0_D5_C4": freqs[0].tolist(),
            "S1_D5_A4": freqs[1].tolist(),
            "S2_C4_A4": freqs[2].tolist(),
        },
        "delta_max_minus_min": delta.tolist(),
        "max_delta": max_delta,
        "max_delta_interior": max_delta_int,
        "n_interior_edges_with_delta_above_0.30": n_extreme_interior_edges,
    }
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    print(f"\nResultado salvo em: {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
