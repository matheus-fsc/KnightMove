#!/usr/bin/env python3
"""
edge_compatibility_6x6.py
=========================
Para cada uma das 9.862 soluções únicas (ciclo Hamiltoniano fechado no
cavalo 6×6), regista o estado dos 8 vértices de borda B6, C6, D6, E6,
B1, C1, D1, E1 — onde "estado" = qual par de arestas-cavalo está ativo
nesse vértice (cada vértice tem grau d, então C(d,2) estados possíveis).

Reporta:
  • Graus reais de cada vértice da borda (sanity check do enunciado)
  • Contingências (B6 × cada um dos outros 7) — quem acopla mais forte
  • Tabela (B6, E6) e (B6, B1) detalhadas
  • Combinações com ZERO soluções (proibidas pela topologia/paridade)
  • Distribuição: uniforme vs concentrada (chi² vs uniforme, Cramér's V)
  • Counts condicionais: |soluções dado B6| e |soluções dado B6+E6|
"""

import json
import os
import sys
from collections import defaultdict
from itertools import combinations
from math import comb

import numpy as np


BOARD = 6
TOTAL = BOARD * BOARD
MOVES = [(2, 1), (2, -1), (-2, 1), (-2, -1),
         (1, 2), (1, -2), (-1, 2), (-1, -2)]


def vid(r, c): return r * BOARD + c
def vrc(v): return divmod(v, BOARD)
def label(v):
    r, c = vrc(v)
    return chr(65 + c) + str(BOARD - r)


def knight_neighbors(v):
    r, c = vrc(v)
    nbrs = []
    for dr, dc in MOVES:
        nr, nc = r + dr, c + dc
        if 0 <= nr < BOARD and 0 <= nc < BOARD:
            nbrs.append(vid(nr, nc))
    return sorted(nbrs)


def build_edges():
    EDGES = set()
    for v in range(TOTAL):
        for u in knight_neighbors(v):
            EDGES.add((min(u, v), max(u, v)))
    return sorted(EDGES)


def load_unique_cycles():
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "destruction_catalogue.json")) as f:
        cat = json.load(f)
    label_to_vid = {label(v): v for v in range(TOTAL)}
    seen = {}
    for entry in cat:
        seq = [label_to_vid[s] for s in entry["sequence"]]
        es = set()
        for i in range(len(seq)):
            a, b = seq[i], seq[(i + 1) % len(seq)]
            es.add((min(a, b), max(a, b)))
        fs = frozenset(es)
        if fs not in seen:
            seen[fs] = es
    return seen


def vertex_state_encoder(v, edges, edge_to_idx):
    """Retorna (n_states, encoder_fn, decoder_dict) para o vértice v.

    Estado = bitmask das arestas incidentes ativas. Para vértice interior
    de um ciclo, sempre exatamente 2 bits ativos: state ∈ C(deg, 2) possíveis.
    Codificamos como inteiro = bitmask.
    """
    nbrs = knight_neighbors(v)
    incident_edges = [(min(u, v), max(u, v)) for u in nbrs]
    incident_idxs = [edge_to_idx[e] for e in incident_edges]
    deg = len(nbrs)

    def encode_state(sig_bits):
        s = 0
        for j, idx in enumerate(incident_idxs):
            if sig_bits[idx]:
                s |= (1 << j)
        return s

    # estados possíveis: subconjuntos de tamanho 2
    state_to_label = {}
    for i, j in combinations(range(deg), 2):
        st = (1 << i) | (1 << j)
        nbr_labels = sorted([label(nbrs[i]), label(nbrs[j])])
        state_to_label[st] = "{" + ",".join(nbr_labels) + "}"

    return deg, nbrs, encode_state, state_to_label


def main():
    EDGES = build_edges()
    n_edges = len(EDGES)
    edge_to_idx = {e: i for i, e in enumerate(EDGES)}
    label_to_vid = {label(v): v for v in range(TOTAL)}

    # ── vértices a estudar ─────────────────────────────────────────
    BORDER_TOP = ["B6", "C6", "D6", "E6"]
    BORDER_BOT = ["B1", "C1", "D1", "E1"]
    ALL = BORDER_TOP + BORDER_BOT

    print("=" * 72)
    print("VERIFICAÇÃO DE GRAUS (sanity check do enunciado)")
    print("=" * 72)
    print()
    print("  O enunciado diz 'cada um com 3 estados possíveis'.")
    print("  Estados = C(grau,2). Grau 3 → 3 estados. Grau 4 → 6 estados.\n")
    deg_info = {}
    for lab in ALL:
        v = label_to_vid[lab]
        nbrs = knight_neighbors(v)
        deg = len(nbrs)
        nbr_labels = ", ".join(sorted(label(u) for u in nbrs))
        n_states = comb(deg, 2)
        deg_info[lab] = (v, deg, n_states, nbrs)
        tag = "✓ deg 3" if deg == 3 else "✗ deg 4 (6 estados, não 3)"
        print(f"  {lab}: grau {deg} ({tag}) | vizinhos: {nbr_labels}")

    # ── carrega cycles ─────────────────────────────────────────────
    print(f"\n  Carregando 9.862 ciclos únicos…")
    seen = load_unique_cycles()
    n_unique = len(seen)
    print(f"  ciclos únicos: {n_unique:,}")

    # ── monta matriz de assinaturas ────────────────────────────────
    X = np.zeros((n_unique, n_edges), dtype=np.bool_)
    for k, es in enumerate(seen.values()):
        for e in es:
            X[k, edge_to_idx[e]] = True

    # ── codifica estados de cada vértice de borda ──────────────────
    state_data = {}  # lab -> dict
    for lab in ALL:
        v, deg, n_states, nbrs = deg_info[lab]
        deg2, _, encode_state, state_to_label = vertex_state_encoder(
            v, EDGES, edge_to_idx)
        states = np.array([encode_state(X[k]) for k in range(n_unique)],
                          dtype=np.int32)
        state_data[lab] = {
            "deg": deg, "n_states": n_states,
            "states": states,
            "state_to_label": state_to_label,
        }

    # ── distribuições marginais ────────────────────────────────────
    print(f"\n" + "=" * 72)
    print("DISTRIBUIÇÕES MARGINAIS (estado de cada vértice)")
    print("=" * 72)
    for lab in ALL:
        d = state_data[lab]
        states = d["states"]
        unique, counts = np.unique(states, return_counts=True)
        H = -np.sum(counts/n_unique * np.log2(counts/n_unique))
        Hmax = np.log2(d["n_states"])
        print(f"\n  {lab} (grau {d['deg']}, {d['n_states']} estados possíveis, "
              f"{len(unique)} ocorrem):")
        print(f"    H = {H:.3f} / Hmax={Hmax:.3f} bits  "
              f"(uniformidade = {H/Hmax:.1%})")
        for st, cnt in zip(unique, counts):
            print(f"      {d['state_to_label'][int(st)]:<14s}  "
                  f"{cnt:>5,} ({cnt/n_unique:6.2%})")

    # ── contingências (B6, X) para X ∈ outros 7 ────────────────────
    print(f"\n" + "=" * 72)
    print("CONTINGÊNCIAS B6 × X — quem acopla mais com B6?")
    print("=" * 72)

    def cramers_v(x, y):
        cats_x, _ = np.unique(x, return_counts=True)
        cats_y, _ = np.unique(y, return_counts=True)
        r, c = len(cats_x), len(cats_y)
        table = np.zeros((r, c), dtype=np.int64)
        for i, vx in enumerate(cats_x):
            mx = (x == vx)
            for j, vy in enumerate(cats_y):
                table[i, j] = int((mx & (y == vy)).sum())
        n = int(table.sum())
        rs = table.sum(axis=1, keepdims=True)
        cs = table.sum(axis=0, keepdims=True)
        exp = (rs * cs) / n
        with np.errstate(divide="ignore", invalid="ignore"):
            chi2 = float(np.where(exp > 0, (table - exp) ** 2 / exp, 0).sum())
        denom = n * max(min(r - 1, c - 1), 1)
        V = float(np.sqrt(chi2 / denom)) if denom > 0 else 0.0
        return V, chi2, table, cats_x, cats_y

    couplings = []
    b6_states = state_data["B6"]["states"]
    for lab in ALL:
        if lab == "B6":
            continue
        x_states = state_data[lab]["states"]
        V, chi2, table, cx, cy = cramers_v(b6_states, x_states)
        n_zero_cells = int(np.sum(table == 0))
        n_total_cells = table.size
        couplings.append((V, chi2, lab, table, cx, cy, n_zero_cells,
                          n_total_cells))

    couplings.sort(key=lambda x: -x[0])
    print(f"\n  ranking por Cramér's V (associação categórica, V ∈ [0,1]):\n")
    print(f"  {'X':<5s}  {'V (B6,X)':>9s}  {'χ²':>10s}  "
          f"{'#cells':>7s}  {'#zeros':>7s}")
    for V, chi2, lab, table, cx, cy, n_zero, n_total in couplings:
        print(f"  {lab:<5s}  {V:>9.4f}  {chi2:>10.1f}  "
              f"{n_total:>7d}  {n_zero:>7d}")

    # ── tabelas detalhadas: (B6, E6) e (B6, B1) ────────────────────
    print(f"\n" + "=" * 72)
    print("TABELAS DETALHADAS (pedidos explicitamente)")
    print("=" * 72)
    b6_decoder = state_data["B6"]["state_to_label"]
    for target in ["E6", "B1"]:
        print(f"\n  --- B6 × {target} ---")
        x_states = state_data[target]["states"]
        V, chi2, table, cx, cy = cramers_v(b6_states, x_states)
        x_decoder = state_data[target]["state_to_label"]
        col_labels = [x_decoder[int(c)] for c in cy]
        row_labels = [b6_decoder[int(c)] for c in cx]
        col_widths = max(14, max(len(l) for l in col_labels) + 2)
        hdr = " " * 18 + "  ".join(f"{l:>{col_widths}s}" for l in col_labels)
        print(f"  V = {V:.4f}, χ² = {chi2:.1f}")
        print()
        print(hdr)
        for i, rl in enumerate(row_labels):
            cells = "  ".join(f"{table[i,j]:>{col_widths},}"
                              for j in range(table.shape[1]))
            print(f"  {rl:<16s}{cells}")
        # combinações zero
        zeros = [(row_labels[i], col_labels[j])
                 for i in range(table.shape[0])
                 for j in range(table.shape[1]) if table[i, j] == 0]
        if zeros:
            print(f"\n  combinações ZERO ({len(zeros)}):")
            for rl, cl in zeros:
                print(f"    B6={rl}  &  {target}={cl}")
        else:
            print(f"\n  nenhuma combinação zero — toda combinatorial está povoada")

    # ── counts condicionais ─────────────────────────────────────────
    print(f"\n" + "=" * 72)
    print("COUNTS CONDICIONAIS (multiplicador de poda implícita)")
    print("=" * 72)
    print()
    b6_unique = np.unique(b6_states)
    print(f"  Total de soluções: {n_unique:,}\n")
    print(f"  Fixando B6 sozinho ({len(b6_unique)} valores):")
    for st in b6_unique:
        mask = (b6_states == st)
        cnt = int(mask.sum())
        print(f"    B6={b6_decoder[int(st)]:<14s} → {cnt:>5,} soluções  "
              f"({cnt/n_unique:6.2%})")
    avg_B6 = n_unique / len(b6_unique)
    print(f"    média se uniforme: {avg_B6:,.0f}")

    print(f"\n  Fixando (B6, E6) — {3*3}=9 combinações possíveis:")
    e6_states = state_data["E6"]["states"]
    e6_decoder = state_data["E6"]["state_to_label"]
    table_be = np.zeros((3, 3), dtype=np.int64)
    cells = []
    for i, b in enumerate(np.unique(b6_states)):
        for j, e in enumerate(np.unique(e6_states)):
            cnt = int(((b6_states == b) & (e6_states == e)).sum())
            table_be[i, j] = cnt
            cells.append((b, e, cnt))
    counts_only = [c for _, _, c in cells]
    nonzero = [c for c in counts_only if c > 0]
    print(f"    combinações observadas    : {len(nonzero)}/9")
    print(f"    combinações zero          : {9 - len(nonzero)}/9")
    print(f"    mediana das não-zero       : {int(np.median(nonzero))}")
    print(f"    média das não-zero         : {np.mean(nonzero):.0f}")
    print(f"    mín / máx não-zero         : {min(nonzero)} / {max(nonzero)}")

    print(f"\n  Fixando (B6, B1) — {3*3}=9 combinações possíveis:")
    b1_states = state_data["B1"]["states"]
    table_bb = np.zeros((3, 3), dtype=np.int64)
    cells_bb = []
    for i, b in enumerate(np.unique(b6_states)):
        for j, e in enumerate(np.unique(b1_states)):
            cnt = int(((b6_states == b) & (b1_states == e)).sum())
            table_bb[i, j] = cnt
            cells_bb.append((b, e, cnt))
    counts_only_bb = [c for _, _, c in cells_bb]
    nonzero_bb = [c for c in counts_only_bb if c > 0]
    print(f"    combinações observadas    : {len(nonzero_bb)}/9")
    print(f"    combinações zero          : {9 - len(nonzero_bb)}/9")
    print(f"    mediana das não-zero       : {int(np.median(nonzero_bb))}")
    print(f"    média das não-zero         : {np.mean(nonzero_bb):.0f}")
    print(f"    mín / máx não-zero         : {min(nonzero_bb)} / {max(nonzero_bb)}")

    # ── salvar ─────────────────────────────────────────────────────
    out = {
        "board": BOARD,
        "n_unique_cycles": n_unique,
        "degrees": {lab: state_data[lab]["deg"] for lab in ALL},
        "n_states_possible": {lab: state_data[lab]["n_states"] for lab in ALL},
        "couplings_with_B6": [
            {"X": lab, "cramers_v": V, "chi2": chi2,
             "zero_cells": n_zero, "total_cells": n_total}
            for V, chi2, lab, _, _, _, n_zero, n_total in couplings
        ],
        "table_B6_E6": table_be.tolist(),
        "table_B6_B1": table_bb.tolist(),
    }
    here = os.path.dirname(os.path.abspath(__file__))
    out_path = os.path.join(here, "edge_compatibility_6x6.json")
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    print(f"\nResultado salvo em: {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
