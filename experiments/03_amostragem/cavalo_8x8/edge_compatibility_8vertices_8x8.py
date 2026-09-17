#!/usr/bin/env python3
"""
edge_compatibility_8vertices_8x8.py
====================================
Compatibilidade de borda no 8×8: fixa o estado dos 8 vértices de grau 3
(uma órbita D₄ completa: B8, G8, B1, G1, A7, A2, H7, H2).

Cada um tem 3 estados → 3⁸ = 6.561 combinações.

Usa Phase 1 (data_parallel/samples/) — pairs cujos endpoints NÃO estão
nos 8 vértices de grau 3 (par 0, 2, 3, 7). Subamostra 50k pre-D₄ por par
para uniformidade estatística (200k total).

Reporta + compara explicitamente com o 6×6:
  - % combinações zero  (6×6 foi 8.6%; previsão 10-15%)
  - Total correlation   (6×6 foi 0.221 bits = 3.9%; previsão 3-6%)
  - Spread max/min      (6×6 foi 131×; previsão >> 131×)
  - Top pairwise        (6×6 foi B6↔E6 same-row; previsão B8↔G8)
"""

import glob
import json
import os
import sys
from collections import Counter, defaultdict
from itertools import combinations, product

import numpy as np


BOARD = 8
TOTAL = BOARD * BOARD
MOVES = [(2, 1), (2, -1), (-2, 1), (-2, -1),
         (1, 2), (1, -2), (-1, 2), (-1, -2)]

# os 8 vértices de grau 3 (D₄-órbita)
DEG3_RC = [(0, 1), (0, 6), (7, 1), (7, 6),
           (1, 0), (6, 0), (1, 7), (6, 7)]
DEG3_LABELS = ["B8", "G8", "B1", "G1", "A7", "A2", "H7", "H2"]


def vid(r, c): return r * BOARD + c
def vrc(v): return divmod(v, BOARD)
def label(v):
    r, c = vrc(v)
    return chr(65 + c) + str(BOARD - r)


def knight_neighbors_rc(r, c):
    nbrs = []
    for dr, dc in MOVES:
        nr, nc = r + dr, c + dc
        if 0 <= nr < BOARD and 0 <= nc < BOARD:
            nbrs.append((nr, nc))
    return sorted(nbrs)


def vertex_incident_edge_idxs(vert_rc, edges_list, edge_to_idx):
    """Retorna lista (idx, neighbor_label) ordenada por neighbor_label."""
    r, c = vert_rc
    me = label(vid(r, c))
    items = []
    for nr, nc in knight_neighbors_rc(r, c):
        u_lab = label(vid(nr, nc))
        e1 = f"{me}-{u_lab}"
        e2 = f"{u_lab}-{me}"
        idx = edge_to_idx.get(e1, edge_to_idx.get(e2))
        if idx is None:
            raise KeyError(f"edge {e1}/{e2} não em edges_list")
        items.append((idx, u_lab))
    items.sort(key=lambda x: x[1])
    return items


def encode_states_batch(X, vertex_incident_idxs):
    """X: (n, n_edges) bool. vertex_incident_idxs: list of 3 ints.
    Retorna array (n,) com bitmask."""
    state = np.zeros(X.shape[0], dtype=np.int8)
    for j, idx in enumerate(vertex_incident_idxs):
        state |= (X[:, idx].astype(np.int8) << j)
    return state


def state_to_str(st, incident_items):
    parts = sorted(u for j, (_, u) in enumerate(incident_items)
                    if st & (1 << j))
    return "{" + ",".join(parts) + "}"


def main():
    # ── selecionar pares "OK" ──────────────────────────────────────
    PAIRS_OK = [0, 2, 3, 7]   # endpoints fora dos 8 deg-3 vertices
    MAX_PRE_D4 = 50_000        # por par

    # ── descobrir edges (primeiro batch) ───────────────────────────
    sample_dir = "data_parallel/samples"
    first_b = sorted(glob.glob(os.path.join(sample_dir, "batch_p00_*.json")))[0]
    with open(first_b) as f:
        b0 = json.load(f)
    edges_list = b0["edges"]
    edge_to_idx = {e: i for i, e in enumerate(edges_list)}
    n_edges = len(edges_list)
    print(f"edges: {n_edges}")

    # ── indices incidentes p/ cada um dos 8 vertices ───────────────
    incident_by_v = {}
    for rc, lab in zip(DEG3_RC, DEG3_LABELS):
        items = vertex_incident_edge_idxs(rc, edges_list, edge_to_idx)
        incident_by_v[lab] = items
        nbr_labels = [u for _, u in items]
        idxs = [i for i, _ in items]
        print(f"  {lab} ({rc}): vizinhos={nbr_labels}  idxs={idxs}")

    # ── carrega samples dos pares OK ───────────────────────────────
    print(f"\nCarregando ≤{MAX_PRE_D4:,} amostras/par dos pares {PAIRS_OK}…")
    all_sigs = []
    n_per_pair = {}
    for pidx in PAIRS_OK:
        batches = sorted(glob.glob(
            os.path.join(sample_dir, f"batch_p{pidx:02d}_*.json")))
        n_loaded = 0
        sigs_pair = []
        for bp in batches:
            if n_loaded >= MAX_PRE_D4:
                break
            with open(bp) as f:
                bj = json.load(f)
            if bj["edges"] != edges_list:
                raise ValueError(f"ordem inconsistente em {bp}")
            sigs = np.asarray(bj["signatures"], dtype=np.bool_)
            take = min(sigs.shape[0], MAX_PRE_D4 - n_loaded)
            sigs_pair.append(sigs[:take])
            n_loaded += take
        Xpair = np.concatenate(sigs_pair, axis=0)
        n_per_pair[pidx] = Xpair.shape[0]
        all_sigs.append(Xpair)
        print(f"  par {pidx}: {Xpair.shape[0]:,} amostras carregadas")
    X = np.concatenate(all_sigs, axis=0)
    n_unique = X.shape[0]
    print(f"\n  X.shape = {X.shape}  total={n_unique:,}")

    # ── estados dos 8 vértices ─────────────────────────────────────
    states_by_v = {}
    for lab in DEG3_LABELS:
        incident_idxs = [i for i, _ in incident_by_v[lab]]
        states_by_v[lab] = encode_states_batch(X, incident_idxs)
    # sanity: verifica que cada vertice tem 3 estados observados
    for lab in DEG3_LABELS:
        uniq, cnt = np.unique(states_by_v[lab], return_counts=True)
        rare = uniq[cnt == cnt.min()][0] if len(uniq) else None
        # filtrar estados inválidos: deve ter exatamente 2 bits ativos
        pop_invalid = [u for u in uniq if bin(int(u)).count("1") != 2]
        if pop_invalid:
            print(f"  AVISO {lab}: estados com !=2 bits ativos: {pop_invalid}")
        print(f"  {lab}: {len(uniq)} estados, distrib %: ",
              [f'{c/n_unique:.2%}' for c in cnt])

    # ── marginais para cada vertice ────────────────────────────────
    state_marg = {}
    state_lists = {}
    for lab in DEG3_LABELS:
        uniq, cnt = np.unique(states_by_v[lab], return_counts=True)
        order = np.argsort(-cnt)
        state_lists[lab] = [int(uniq[i]) for i in order]
        state_marg[lab] = {int(uniq[i]): int(cnt[i]) for i in range(len(uniq))}

    # ── contingência 8-way (3^8 = 6561) ────────────────────────────
    print(f"\nMontando contingência 8-way (3^8 = 6561 cells)…")
    tup_arr = np.stack([states_by_v[lab] for lab in DEG3_LABELS], axis=1)
    # convert each row to a single hash key (8 bytes packed)
    keys = (tup_arr[:, 0].astype(np.int64)
            | (tup_arr[:, 1].astype(np.int64) << 8)
            | (tup_arr[:, 2].astype(np.int64) << 16)
            | (tup_arr[:, 3].astype(np.int64) << 24)
            | (tup_arr[:, 4].astype(np.int64) << 32)
            | (tup_arr[:, 5].astype(np.int64) << 40)
            | (tup_arr[:, 6].astype(np.int64) << 48)
            | (tup_arr[:, 7].astype(np.int64) << 56))
    uniq_keys, obs_counts = np.unique(keys, return_counts=True)
    n_observed_combos = len(uniq_keys)
    n_total_combos = 1
    for lab in DEG3_LABELS:
        n_total_combos *= len(state_lists[lab])
    n_zero = n_total_combos - n_observed_combos
    print(f"  combinações possíveis: {n_total_combos}")
    print(f"  combinações observadas: {n_observed_combos}")
    print(f"  combinações zero      : {n_zero} "
          f"({n_zero/n_total_combos*100:.1f}%)")

    obs_arr = obs_counts.astype(np.int64)
    cmax = int(obs_arr.max())
    cmin_nz = int(obs_arr.min())
    mean_nz = float(obs_arr.mean())
    median_nz = float(np.median(obs_arr))
    spread = cmax / cmin_nz if cmin_nz else float("inf")
    mean_all = n_unique / n_total_combos
    reduction_mean = n_unique / mean_all  # = n_total_combos

    print(f"\n  N total amostras     : {n_unique:,}")
    print(f"  média se uniforme    : {mean_all:.1f}")
    print(f"  min não-zero → max   : {cmin_nz} → {cmax} (spread = {spread:.0f}×)")
    print(f"  mediana não-zero     : {median_nz:.0f}")
    print(f"  média não-zero        : {mean_nz:.1f}")

    # ── total correlation: H(margins) - H(joint) ───────────────────
    H_marg = 0.0
    for lab in DEG3_LABELS:
        for st, c in state_marg[lab].items():
            p = c / n_unique
            if p > 0:
                H_marg += -p * np.log2(p)
    H_joint = 0.0
    for c in obs_arr:
        p = c / n_unique
        if p > 0:
            H_joint += -p * np.log2(p)
    total_corr = H_marg - H_joint
    print(f"\n  H(marginais, 8 vértices) = {H_marg:.4f} bits")
    print(f"  H(joint)                  = {H_joint:.4f} bits")
    print(f"  Total correlation         = {total_corr:.4f} bits")
    print(f"  fraction shared          = {total_corr/H_marg*100:.2f}%")

    # ── pairwise Cramér's V (C(8,2)=28 pares) ──────────────────────
    def cramers_v(x, y):
        cx, _ = np.unique(x, return_counts=True)
        cy, _ = np.unique(y, return_counts=True)
        r, c = len(cx), len(cy)
        T = np.zeros((r, c), dtype=np.int64)
        for i, vx in enumerate(cx):
            mx = (x == vx)
            for j, vy in enumerate(cy):
                T[i, j] = int((mx & (y == vy)).sum())
        n = int(T.sum())
        rs = T.sum(axis=1, keepdims=True); cs = T.sum(axis=0, keepdims=True)
        exp = (rs * cs) / n
        with np.errstate(divide="ignore", invalid="ignore"):
            chi2 = float(np.where(exp > 0, (T - exp) ** 2 / exp, 0).sum())
        denom = n * max(min(r - 1, c - 1), 1)
        V = float(np.sqrt(chi2 / denom)) if denom > 0 else 0.0
        return V, chi2, T

    print(f"\n=== ACOPLAMENTO PAIRWISE Cramér's V (28 pares) ===")
    pair_results = []
    for i, j in combinations(range(8), 2):
        l1, l2 = DEG3_LABELS[i], DEG3_LABELS[j]
        V, chi2, T = cramers_v(states_by_v[l1], states_by_v[l2])
        # caracteriza relação geométrica
        rc1 = DEG3_RC[i]; rc2 = DEG3_RC[j]
        same_row = (rc1[0] == rc2[0])
        same_col = (rc1[1] == rc2[1])
        diag = abs(rc1[0] - rc2[0]) == abs(rc1[1] - rc2[1])
        rel = "same-row" if same_row else ("same-col" if same_col else
              ("diagonal" if diag else "other"))
        pair_results.append((V, chi2, l1, l2, rel, T))
    pair_results.sort(key=lambda x: -x[0])
    print(f"  {'par':<10s}  {'rel':<10s}  {'V':>8s}  {'χ²':>10s}")
    for V, chi2, l1, l2, rel, _ in pair_results:
        print(f"  {l1}↔{l2:<5s}  {rel:<10s}  {V:>8.4f}  {chi2:>10.1f}")

    # ── histograma ────────────────────────────────────────────────
    print(f"\n=== HISTOGRAMA DE FREQUÊNCIAS ===")
    bins_lo = [0, 1, 5, 10, 50, 100, 500, 1000, 5000]
    bin_lbl = ["0", "1-4", "5-9", "10-49", "50-99",
               "100-499", "500-999", "1000-4999", "5000+"]
    # contagem real (incluindo zeros)
    n_in_bin = [0] * 9
    n_in_bin[0] = n_zero
    for c in obs_arr:
        if c < 5: n_in_bin[1] += 1
        elif c < 10: n_in_bin[2] += 1
        elif c < 50: n_in_bin[3] += 1
        elif c < 100: n_in_bin[4] += 1
        elif c < 500: n_in_bin[5] += 1
        elif c < 1000: n_in_bin[6] += 1
        elif c < 5000: n_in_bin[7] += 1
        else: n_in_bin[8] += 1
    maxbar = max(n_in_bin)
    for lbl, cnt in zip(bin_lbl, n_in_bin):
        bar = "█" * int(40 * cnt / maxbar) if maxbar > 0 else ""
        print(f"  {lbl:>10s} | {cnt:>5d}  {bar}")

    # ── comparação explícita com 6×6 ──────────────────────────────
    print(f"\n" + "=" * 72)
    print(f"COMPARAÇÃO EXPLÍCITA — 6×6 vs 8×8")
    print("=" * 72)
    SIX = {
        "vertices": 4, "combos": 81,
        "zeros": 7, "zero_pct": 8.6,
        "spread": 131.0, "mean_nz": 133.3,
        "tc_bits": 0.2208, "H_marg_bits": 5.6319,
        "tc_pct": 3.92,
        "top_coupling": "B6↔E6 same-row V=0.2326",
    }
    top = pair_results[0]
    EIGHT = {
        "vertices": 8, "combos": n_total_combos,
        "zeros": int(n_zero), "zero_pct": float(n_zero/n_total_combos*100),
        "spread": float(spread), "mean_nz": mean_nz,
        "tc_bits": float(total_corr), "H_marg_bits": float(H_marg),
        "tc_pct": float(total_corr/H_marg*100),
        "top_coupling": f"{top[2]}↔{top[3]} {top[4]} V={top[0]:.4f}",
    }
    fields = [
        ("vértices fixados", "vertices", "{:d}"),
        ("combinações possíveis", "combos", "{:,}"),
        ("# zeros", "zeros", "{:,}"),
        ("% zeros", "zero_pct", "{:.1f}%"),
        ("spread max/min nonzero", "spread", "{:.0f}×"),
        ("média não-zero", "mean_nz", "{:.1f}"),
        ("H(marginais)", "H_marg_bits", "{:.4f} bits"),
        ("Total correlation", "tc_bits", "{:.4f} bits"),
        ("% info compartilhada", "tc_pct", "{:.2f}%"),
        ("top coupling", "top_coupling", "{}"),
    ]
    print(f"  {'métrica':<24s}  {'6×6':>24s}  {'8×8':>30s}")
    print(f"  {'-'*24}  {'-'*24}  {'-'*30}")
    for lbl, k, fmt in fields:
        v6 = fmt.format(SIX[k])
        v8 = fmt.format(EIGHT[k])
        print(f"  {lbl:<24s}  {v6:>24s}  {v8:>30s}")

    # ── salvar ────────────────────────────────────────────────────
    out = {
        "board": BOARD, "vertices_fixed": DEG3_LABELS,
        "n_pairs_used": len(PAIRS_OK), "pairs_used": PAIRS_OK,
        "n_samples_per_pair_cap": MAX_PRE_D4,
        "n_samples_total": n_unique,
        "n_combinations_total": n_total_combos,
        "n_zero": int(n_zero), "zero_pct": float(n_zero/n_total_combos*100),
        "spread": float(spread), "min_nonzero": cmin_nz, "max_count": cmax,
        "median_nonzero": float(median_nz), "mean_nonzero": float(mean_nz),
        "H_marginals_bits": float(H_marg),
        "H_joint_bits": float(H_joint),
        "total_correlation_bits": float(total_corr),
        "total_correlation_fraction": float(total_corr/H_marg),
        "pairwise_couplings": [
            {"v1": l1, "v2": l2, "relation": rel,
             "cramers_v": float(V), "chi2": float(chi2)}
            for V, chi2, l1, l2, rel, _ in pair_results
        ],
        "comparison_6x6": SIX,
    }
    out_path = "edge_compatibility_8vertices_8x8.json"
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    print(f"\nResultado salvo em: {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
