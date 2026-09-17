#!/usr/bin/env python3
"""
edge_compatibility_4vertices_6x6.py
====================================
Fixa o estado dos 4 vértices de grau 3 do cavalo 6×6:
  B6, E6  (borda superior)
  B1, E1  (borda inferior)

Cada um tem 3 estados → 3⁴ = 81 combinações.

Reporta:
  1. # combinações com 0 soluções
  2. Distribuição completa (counts por combinação)
  3. Combinação mais frequente / mais rara
  4. Redução média (9862 / mean) e spread (max/min)
  5. Teste de independência: contagem observada vs predita por
     produto das marginais (qui-quadrado + MI 4-way)
"""

import json
import os
import sys
from collections import Counter
from itertools import product

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
    l2v = {label(v): v for v in range(TOTAL)}
    seen = {}
    for entry in cat:
        seq = [l2v[s] for s in entry["sequence"]]
        es = set()
        for i in range(len(seq)):
            a, b = seq[i], seq[(i + 1) % len(seq)]
            es.add((min(a, b), max(a, b)))
        fs = frozenset(es)
        if fs not in seen:
            seen[fs] = es
    return seen


def encode_state(sig, v, edges, edge_to_idx):
    """Retorna inteiro = bitmask das arestas incidentes ativas em v."""
    nbrs = knight_neighbors(v)
    s = 0
    for j, u in enumerate(nbrs):
        e = (min(u, v), max(u, v))
        if sig[edge_to_idx[e]]:
            s |= (1 << j)
    return s, nbrs


def state_label(state_int, nbrs):
    parts = sorted(label(u) for j, u in enumerate(nbrs) if state_int & (1 << j))
    return "{" + ",".join(parts) + "}"


def main():
    EDGES = build_edges()
    edge_to_idx = {e: i for i, e in enumerate(EDGES)}
    l2v = {label(v): v for v in range(TOTAL)}

    VERTS = ["B6", "E6", "B1", "E1"]
    vids = {lab: l2v[lab] for lab in VERTS}
    nbrs_by_v = {lab: knight_neighbors(vids[lab]) for lab in VERTS}

    print("Carregando 9.862 ciclos…")
    seen = load_unique_cycles()
    n_unique = len(seen)
    n_edges = len(EDGES)

    # ── monta matriz e estados dos 4 vértices ──────────────────────
    X = np.zeros((n_unique, n_edges), dtype=np.bool_)
    states_by_v = {lab: np.zeros(n_unique, dtype=np.int32) for lab in VERTS}
    for k, es in enumerate(seen.values()):
        for e in es:
            X[k, edge_to_idx[e]] = True
        for lab in VERTS:
            st, _ = encode_state(X[k], vids[lab], EDGES, edge_to_idx)
            states_by_v[lab][k] = st

    # ── identifica os 3 estados de cada vértice (ordenados por freq) ──
    state_lists = {}
    state_marg = {}
    for lab in VERTS:
        uniq, cnt = np.unique(states_by_v[lab], return_counts=True)
        order = np.argsort(-cnt)
        state_lists[lab] = [int(uniq[i]) for i in order]
        state_marg[lab] = {int(uniq[i]): int(cnt[i]) for i in range(len(uniq))}

    # rótulos legíveis
    state_label_by_v = {
        lab: {st: state_label(st, nbrs_by_v[lab]) for st in state_lists[lab]}
        for lab in VERTS
    }

    # ── 4D contingência ────────────────────────────────────────────
    # build tuple (B6_state, E6_state, B1_state, E1_state) por ciclo
    tuples = list(zip(states_by_v["B6"], states_by_v["E6"],
                      states_by_v["B1"], states_by_v["E1"]))
    counter = Counter(tuples)
    # garante presença das 81 chaves (zero onde não ocorreu)
    all_keys = list(product(state_lists["B6"], state_lists["E6"],
                             state_lists["B1"], state_lists["E1"]))
    obs = {k: counter.get(k, 0) for k in all_keys}
    counts = np.array(list(obs.values()), dtype=np.int64)
    n_total = int(counts.sum())
    assert n_total == n_unique, f"{n_total} != {n_unique}"

    n_zero = int((counts == 0).sum())
    n_nonzero = int((counts > 0).sum())
    nonzero = counts[counts > 0]
    mean_all = float(counts.mean())                # incluindo zeros
    mean_nz = float(nonzero.mean())                # apenas não-zeros
    median_nz = int(np.median(nonzero))
    cmax = int(counts.max())
    cmin_nz = int(nonzero.min())
    spread = cmax / cmin_nz if cmin_nz else float("inf")
    reduction_mean = n_unique / mean_all if mean_all > 0 else float("inf")
    reduction_mean_nz = n_unique / mean_nz if mean_nz > 0 else float("inf")

    # ── predição sob independência ────────────────────────────────
    # P(combo) = P(B6) × P(E6) × P(B1) × P(E1)
    pred = {}
    for k in all_keys:
        p = 1.0
        for lab, st in zip(VERTS, k):
            p *= state_marg[lab][st] / n_unique
        pred[k] = p * n_unique
    pred_arr = np.array([pred[k] for k in all_keys], dtype=np.float64)
    # chi²
    with np.errstate(divide="ignore", invalid="ignore"):
        chi2_arr = np.where(pred_arr > 0,
                            (counts - pred_arr) ** 2 / pred_arr, 0)
    chi2 = float(chi2_arr.sum())
    # df = (3-1)^4 = 16 ? não: para teste de independência 4-way completo,
    # df = ∏cardinalidade - 1 - Σ(card_i - 1)
    #    = 81 - 1 - 4*(3-1) = 72
    df = 81 - 1 - 4 * 2
    # informação mútua 4-way (total correlation): Σ marginais - H joint
    H_marg = sum(
        -sum((state_marg[lab][st]/n_unique) * np.log2(state_marg[lab][st]/n_unique)
             for st in state_lists[lab])
        for lab in VERTS
    )
    H_joint = -sum((c/n_unique) * np.log2(c/n_unique)
                    for c in nonzero)
    total_corr = H_marg - H_joint  # información compartilhada entre os 4

    # ── relatório ─────────────────────────────────────────────────
    print()
    print("=" * 72)
    print(f"FIXAÇÃO DOS 4 VÉRTICES DE GRAU 3 (B6, E6, B1, E1)  ·  3⁴ = 81")
    print("=" * 72)
    print()
    print(f"  Total de ciclos                  : {n_unique:,}")
    print(f"  Combinações observadas (≥1)       : {n_nonzero}/81 "
          f"({n_nonzero/81:.1%})")
    print(f"  Combinações ZERO                  : {n_zero}/81 "
          f"({n_zero/81:.1%})")
    print(f"  Soluções/combo (incluindo zeros)  : média {mean_all:.1f}  "
          f"esperado uniforme 9862/81 = {9862/81:.1f}")
    print(f"  Soluções/combo (apenas não-zero)  : mediana {median_nz}, "
          f"média {mean_nz:.1f}")
    print(f"  Mín não-zero  → Máx              : {cmin_nz} → {cmax}  "
          f"(spread = {spread:.0f}×)")
    print()
    print(f"  Redução média (9862 / média)      : {reduction_mean:.2f}×  "
          f"(meta uniforme = 81)")
    print(f"  Redução média (apenas não-zero)   : {reduction_mean_nz:.2f}×")

    print()
    print("─── TESTE DE INDEPENDÊNCIA ───")
    print(f"  H(margins) — Σ entropias dos 4    : {H_marg:.4f} bits")
    print(f"  H(joint)   — entropia conjunta   : {H_joint:.4f} bits")
    print(f"  Total correlation (Σ H_i - H_J)  : {total_corr:.4f} bits "
          f"(0 = independente)")
    print(f"  χ² (vs produto das marginais)    : {chi2:.1f}  (df={df})")
    print(f"  desvio relativo médio Σ|O-E|/N    : "
          f"{float(np.abs(counts - pred_arr).sum()/n_unique):.4f}")

    # ── top-5 mais frequentes / mais raros ────────────────────────
    order = np.argsort(-counts)
    print()
    print("─── TOP-5 MAIS FREQUENTES ───")
    print(f"  {'B6':<10s} {'E6':<10s} {'B1':<10s} {'E1':<10s}  "
          f"{'obs':>5s}  {'pred':>5s}  {'(O-E)/√E':>9s}")
    for i in order[:5]:
        k = all_keys[i]
        labs = [state_label_by_v[lab][st] for lab, st in zip(VERTS, k)]
        o = counts[i]; p = pred_arr[i]
        z = (o - p) / np.sqrt(p) if p > 0 else float("nan")
        print(f"  {labs[0]:<10s} {labs[1]:<10s} {labs[2]:<10s} {labs[3]:<10s}  "
              f"{o:>5}  {p:>5.0f}  {z:>+9.2f}")

    print()
    print("─── TOP-5 MAIS RAROS (não-zero) ───")
    order_nz = np.argsort(counts)  # ascending
    shown = 0
    for i in order_nz:
        if counts[i] == 0:
            continue
        k = all_keys[i]
        labs = [state_label_by_v[lab][st] for lab, st in zip(VERTS, k)]
        o = counts[i]; p = pred_arr[i]
        z = (o - p) / np.sqrt(p) if p > 0 else float("nan")
        print(f"  {labs[0]:<10s} {labs[1]:<10s} {labs[2]:<10s} {labs[3]:<10s}  "
              f"{o:>5}  {p:>5.0f}  {z:>+9.2f}")
        shown += 1
        if shown >= 5: break

    if n_zero:
        print()
        print(f"─── COMBINAÇÕES ZERO ({n_zero}) ───")
        for i, k in enumerate(all_keys):
            if counts[i] == 0:
                labs = [state_label_by_v[lab][st] for lab, st in zip(VERTS, k)]
                print(f"  B6={labs[0]:<10s} E6={labs[1]:<10s} "
                      f"B1={labs[2]:<10s} E1={labs[3]:<10s}  "
                      f"(pred={pred_arr[i]:.1f})")

    # ── histograma ────────────────────────────────────────────────
    print()
    print("─── HISTOGRAMA DE FREQUÊNCIAS ───")
    bins = [0, 1, 10, 50, 100, 200, 500, 1000, 5000]
    bin_lbl = ["0", "1-9", "10-49", "50-99", "100-199",
               "200-499", "500-999", "1000+"]
    h = np.zeros(len(bin_lbl), dtype=int)
    for c in counts:
        if c == 0: h[0] += 1
        elif c < 10: h[1] += 1
        elif c < 50: h[2] += 1
        elif c < 100: h[3] += 1
        elif c < 200: h[4] += 1
        elif c < 500: h[5] += 1
        elif c < 1000: h[6] += 1
        else: h[7] += 1
    for lbl, cnt in zip(bin_lbl, h):
        bar = "█" * cnt
        print(f"  {lbl:>10s} | {cnt:>3d}  {bar}")

    # ── salva ─────────────────────────────────────────────────────
    out = {
        "board": BOARD,
        "n_unique_cycles": n_unique,
        "vertices_fixed": VERTS,
        "n_combinations_total": 81,
        "n_zero_combinations": n_zero,
        "n_nonzero_combinations": n_nonzero,
        "mean_count_all": mean_all,
        "mean_count_nonzero": mean_nz,
        "median_count_nonzero": median_nz,
        "max_count": cmax,
        "min_nonzero_count": cmin_nz,
        "spread_max_over_min": spread,
        "reduction_mean": reduction_mean,
        "H_marginals_bits": H_marg,
        "H_joint_bits": H_joint,
        "total_correlation_bits": total_corr,
        "chi2_vs_independence": chi2,
        "chi2_df": df,
        "all_counts": [
            {"B6": state_label_by_v["B6"][k[0]],
             "E6": state_label_by_v["E6"][k[1]],
             "B1": state_label_by_v["B1"][k[2]],
             "E1": state_label_by_v["E1"][k[3]],
             "observed": int(counts[i]),
             "predicted_indep": float(pred_arr[i]),
             "residual_z": float((counts[i] - pred_arr[i]) / np.sqrt(pred_arr[i])
                                  if pred_arr[i] > 0 else 0)}
            for i, k in enumerate(all_keys)
        ],
    }
    here = os.path.dirname(os.path.abspath(__file__))
    out_path = os.path.join(here, "edge_compatibility_4vertices_6x6.json")
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    print(f"\nResultado salvo em: {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
