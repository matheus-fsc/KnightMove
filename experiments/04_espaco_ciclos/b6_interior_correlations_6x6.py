#!/usr/bin/env python3
"""
b6_interior_correlations_6x6.py
================================
Correlações de Pearson entre as 24 arestas estritamente interiores do 6×6,
calculadas DENTRO de cada grupo do estado de B6 (S0={D5,C4}, S1={D5,A4},
S2={C4,A4}).

Aresta estritamente interior = ambos endpoints com row ∈ {1..4} e col ∈ {1..4}.

Para cada grupo:
  - Matriz 24×24 de Pearson
  - Top-5 pares por |r|
  - Classificação dos pares: Tipo A (compartilham vértice) vs
    Tipo B (não compartilham vértice)

Cross-group:
  - União dos top-K pares — quantos aparecem em quais grupos
  - Sinais e magnitudes comparados entre grupos
  - Veredito: o interior tem estrutura própria condicional ao estado de B6?
"""

import json
import os
import sys
from collections import defaultdict

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
    r, c = vrc(v)
    return 1 <= r <= BOARD - 2 and 1 <= c <= BOARD - 2


def is_interior_edge(edge):
    a, b = edge
    return is_interior_vertex(a) and is_interior_vertex(b)


def edge_pair_type(e1, e2):
    """Tipo A se compartilham vértice, Tipo B se não."""
    vs1 = set(e1)
    vs2 = set(e2)
    shared = vs1 & vs2
    if shared:
        return "A", next(iter(shared))
    return "B", None


def pearson_matrix(X):
    """X: (n, k) bool/int. Retorna (k,k) Pearson, com NaN onde variância=0."""
    Xf = X.astype(np.float64)
    n = Xf.shape[0]
    mu = Xf.mean(axis=0)
    Xc = Xf - mu
    cov = (Xc.T @ Xc) / n
    var = np.diag(cov)
    sd = np.sqrt(var)
    constant = sd < 1e-12
    sd_safe = np.where(constant, 1.0, sd)
    corr = cov / (sd_safe[:, None] * sd_safe[None, :])
    for i, c in enumerate(constant):
        if c:
            corr[i, :] = np.nan
            corr[:, i] = np.nan
            corr[i, i] = np.nan
    return corr, mu, var


def top_pairs(corr, top_k=5):
    """Retorna lista [(|r|, r, i, j)] ordenada por |r| desc, j > i."""
    K = corr.shape[0]
    pairs = []
    for i in range(K):
        for j in range(i + 1, K):
            r = corr[i, j]
            if np.isnan(r):
                continue
            pairs.append((abs(r), float(r), i, j))
    pairs.sort(key=lambda x: -x[0])
    return pairs[:top_k] if top_k else pairs


def load_unique_cycles():
    """Carrega destruction_catalogue.json e deduplica ciclos."""
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


def main():
    EDGES = build_edges()
    n_edges = len(EDGES)
    edge_to_idx = {e: i for i, e in enumerate(EDGES)}
    label_to_vid = {label(v): v for v in range(TOTAL)}

    # ── identifica arestas interiores ─────────────────────────────────
    interior_edge_mask = np.array(
        [is_interior_edge(e) for e in EDGES], dtype=bool)
    interior_idxs = np.where(interior_edge_mask)[0]
    interior_edges = [EDGES[i] for i in interior_idxs]
    n_int = len(interior_edges)
    print(f"Arestas estritamente interiores: {n_int}")
    print("  vértices interiores: rows 1-4, cols 1-4 (16 vértices)")
    print("  Lista:")
    for k, e in enumerate(interior_edges):
        lab = f"{label(e[0])}-{label(e[1])}"
        print(f"    [{k:>2}]  {lab:>5}  vid=({e[0]:>2},{e[1]:>2})")

    # ── carrega e dedup ciclos ────────────────────────────────────────
    print(f"\nCarregando + dedup ciclos do destruction_catalogue.json …")
    seen = load_unique_cycles()
    n_unique = len(seen)
    print(f"  ciclos únicos: {n_unique:,}")

    # ── monta matriz de assinaturas + estado de B6 ────────────────────
    B6 = label_to_vid["B6"]
    D5 = label_to_vid["D5"]
    C4 = label_to_vid["C4"]
    A4 = label_to_vid["A4"]
    e_B6_D5 = (min(B6, D5), max(B6, D5))
    e_B6_C4 = (min(B6, C4), max(B6, C4))
    e_B6_A4 = (min(B6, A4), max(B6, A4))
    i_D5 = edge_to_idx[e_B6_D5]
    i_C4 = edge_to_idx[e_B6_C4]
    i_A4 = edge_to_idx[e_B6_A4]

    state_lookup = {
        (True, True, False):  0,    # {D5, C4}
        (True, False, True):  1,    # {D5, A4}
        (False, True, True):  2,    # {C4, A4}
    }

    X = np.zeros((n_unique, n_edges), dtype=np.bool_)
    states = np.zeros(n_unique, dtype=np.int8)
    for k, es in enumerate(seen.values()):
        for e in es:
            X[k, edge_to_idx[e]] = True
        b = (X[k, i_D5], X[k, i_C4], X[k, i_A4])
        states[k] = state_lookup[b]
    print(f"  matriz X: {X.shape}, dtype={X.dtype}")

    # ── particiona em 3 grupos + computa Pearson interior 24×24 ───────
    state_desc = {
        0: ("S0", "{D5,C4}", "ambas interiores"),
        1: ("S1", "{D5,A4}", "interior + borda"),
        2: ("S2", "{C4,A4}", "interior + borda"),
    }

    corrs = {}
    means = {}
    sizes = {}
    print(f"\n=== Matriz de Pearson 24×24 por grupo ===")
    for s in (0, 1, 2):
        idxs = np.where(states == s)[0]
        X_sub = X[idxs][:, interior_idxs]
        sizes[s] = int(len(idxs))
        corr, mu, var = pearson_matrix(X_sub)
        corrs[s] = corr
        means[s] = mu
        n_const = int(np.sum(var < 1e-12))
        print(f"\n  Grupo {state_desc[s][0]} {state_desc[s][1]}: "
              f"{sizes[s]:>5,} ciclos | {n_const} arestas com variância 0 "
              f"(NaN na matriz)")

    # ── matriz global (sem grupo) para comparação ─────────────────────
    print(f"\n=== Matriz de Pearson 24×24 — todos os 9.862 (sem grupo) ===")
    corr_all, mu_all, var_all = pearson_matrix(X[:, interior_idxs])
    n_const_all = int(np.sum(var_all < 1e-12))
    print(f"  {n_const_all} arestas constantes em todo o ensemble")

    # ── top-5 pares por grupo ─────────────────────────────────────────
    print(f"\n=== Top-5 pares interiores por |r| ===\n")
    top_by_group = {}
    for s in (0, 1, 2):
        top = top_pairs(corrs[s], top_k=5)
        top_by_group[s] = top
        sd = state_desc[s]
        print(f"  Grupo {sd[0]} {sd[1]} ({sizes[s]:,} ciclos):")
        for absr, r, i, j in top:
            e1 = interior_edges[i]
            e2 = interior_edges[j]
            l1 = f"{label(e1[0])}-{label(e1[1])}"
            l2 = f"{label(e2[0])}-{label(e2[1])}"
            tipo, shared = edge_pair_type(e1, e2)
            shared_lab = f"@{label(shared)}" if shared else "—"
            print(f"    {l1:>5} ↔ {l2:<5}  r={r:+.4f}  "
                  f"Tipo {tipo} (compart: {shared_lab})")
        print()

    # ── top-5 global (sem grupo) ──────────────────────────────────────
    print(f"  Global (sem grupo, 9.862 ciclos):")
    top_all = top_pairs(corr_all, top_k=5)
    for absr, r, i, j in top_all:
        e1 = interior_edges[i]
        e2 = interior_edges[j]
        l1 = f"{label(e1[0])}-{label(e1[1])}"
        l2 = f"{label(e2[0])}-{label(e2[1])}"
        tipo, shared = edge_pair_type(e1, e2)
        shared_lab = f"@{label(shared)}" if shared else "—"
        print(f"    {l1:>5} ↔ {l2:<5}  r={r:+.4f}  "
              f"Tipo {tipo} (compart: {shared_lab})")

    # ── união dos top-5: aparece em quais grupos? ────────────────────
    print(f"\n=== Pares no top-5 de algum grupo — cross-group ===\n")
    union_pairs = set()
    for s in (0, 1, 2):
        for absr, r, i, j in top_by_group[s]:
            union_pairs.add((i, j))
    print(f"  {len(union_pairs)} pares distintos no top-5 de algum grupo")
    print()
    hdr = (f"  {'i ↔ j (aresta)':<24s}  {'Tipo':<5s}  "
           f"{'r_S0':>8s}  {'r_S1':>8s}  {'r_S2':>8s}  {'r_global':>8s}  "
           f"{'in_top5':>10s}")
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))

    rows_for_summary = []
    for (i, j) in sorted(union_pairs):
        e1 = interior_edges[i]
        e2 = interior_edges[j]
        l1 = f"{label(e1[0])}-{label(e1[1])}"
        l2 = f"{label(e2[0])}-{label(e2[1])}"
        tipo, shared = edge_pair_type(e1, e2)
        r0 = corrs[0][i, j]
        r1 = corrs[1][i, j]
        r2 = corrs[2][i, j]
        rg = corr_all[i, j]
        in_top = []
        for s in (0, 1, 2):
            if (i, j) in {(a, b) for _, _, a, b in top_by_group[s]}:
                in_top.append(f"S{s}")
        # imprime
        def fmt(x):
            return f"{x:+.4f}" if not np.isnan(x) else "  nan  "
        pair_lab = f"{l1} ↔ {l2}"
        print(f"  {pair_lab:<24s}  {tipo:<5s}  "
              f"{fmt(r0):>8s}  {fmt(r1):>8s}  {fmt(r2):>8s}  {fmt(rg):>8s}  "
              f"{','.join(in_top):>10s}")
        rows_for_summary.append({
            "i": int(i), "j": int(j),
            "edge_i": l1, "edge_j": l2,
            "type": tipo, "shared_vertex": label(shared) if shared else None,
            "r_S0": float(r0) if not np.isnan(r0) else None,
            "r_S1": float(r1) if not np.isnan(r1) else None,
            "r_S2": float(r2) if not np.isnan(r2) else None,
            "r_global": float(rg) if not np.isnan(rg) else None,
            "in_top5": in_top,
        })

    # ── métricas agregadas: |r| médio e máximo por grupo ─────────────
    print(f"\n=== Resumo agregado por grupo ===\n")
    print(f"  {'grupo':<8s}  {'n':>6s}  {'max |r|':>9s}  {'média |r|':>10s}  "
          f"{'#pares |r|>0.30':>17s}  {'#pares |r|>0.20':>17s}")
    summary_per_group = {}
    for s in (0, 1, 2):
        C = corrs[s]
        K = C.shape[0]
        abs_vals = []
        for i in range(K):
            for j in range(i + 1, K):
                if not np.isnan(C[i, j]):
                    abs_vals.append(abs(C[i, j]))
        abs_arr = np.array(abs_vals)
        max_r = float(abs_arr.max()) if len(abs_arr) else 0.0
        mean_r = float(abs_arr.mean()) if len(abs_arr) else 0.0
        n30 = int((abs_arr > 0.30).sum())
        n20 = int((abs_arr > 0.20).sum())
        sd = state_desc[s]
        print(f"  {sd[0]+' '+sd[1]:<8s}  {sizes[s]:>6,}  {max_r:>9.4f}  "
              f"{mean_r:>10.4f}  {n30:>17d}  {n20:>17d}")
        summary_per_group[sd[0]] = {
            "size": sizes[s], "max_abs_r": max_r,
            "mean_abs_r": mean_r,
            "n_above_0.30": n30, "n_above_0.20": n20,
        }

    # global
    abs_global = []
    K = corr_all.shape[0]
    for i in range(K):
        for j in range(i + 1, K):
            if not np.isnan(corr_all[i, j]):
                abs_global.append(abs(corr_all[i, j]))
    abs_g = np.array(abs_global)
    print(f"  {'GLOBAL':<8s}  {n_unique:>6,}  {abs_g.max():>9.4f}  "
          f"{abs_g.mean():>10.4f}  {int((abs_g>0.30).sum()):>17d}  "
          f"{int((abs_g>0.20).sum()):>17d}")

    # ── classificação Tipo A vs B nos pares > 0.30 ────────────────────
    print(f"\n=== Tipos geométricos dos pares com |r| > 0.30 ===\n")
    for s in (0, 1, 2):
        C = corrs[s]
        K = C.shape[0]
        type_a = 0
        type_b = 0
        examples_a = []
        examples_b = []
        for i in range(K):
            for j in range(i + 1, K):
                if np.isnan(C[i, j]) or abs(C[i, j]) <= 0.30:
                    continue
                e1 = interior_edges[i]
                e2 = interior_edges[j]
                tipo, shared = edge_pair_type(e1, e2)
                if tipo == "A":
                    type_a += 1
                    if len(examples_a) < 3:
                        l1 = f"{label(e1[0])}-{label(e1[1])}"
                        l2 = f"{label(e2[0])}-{label(e2[1])}"
                        examples_a.append(
                            f"{l1}↔{l2} (@{label(shared)}, r={C[i,j]:+.3f})")
                else:
                    type_b += 1
                    if len(examples_b) < 3:
                        l1 = f"{label(e1[0])}-{label(e1[1])}"
                        l2 = f"{label(e2[0])}-{label(e2[1])}"
                        examples_b.append(f"{l1}↔{l2} (r={C[i,j]:+.3f})")
        sd = state_desc[s]
        print(f"  Grupo {sd[0]}: "
              f"{type_a} pares Tipo A, {type_b} pares Tipo B")
        if examples_a:
            print(f"    ex. Tipo A: {'; '.join(examples_a)}")
        if examples_b:
            print(f"    ex. Tipo B: {'; '.join(examples_b)}")

    # ── consistência cross-group: top-5 union recurrence ─────────────
    print(f"\n=== Recorrência dos pares entre grupos (top-5) ===\n")
    rec_count = defaultdict(list)
    for s in (0, 1, 2):
        for absr, r, i, j in top_by_group[s]:
            rec_count[(i, j)].append((s, r))
    n3, n2, n1 = 0, 0, 0
    for k, occ in rec_count.items():
        if len(occ) == 3: n3 += 1
        elif len(occ) == 2: n2 += 1
        else: n1 += 1
    print(f"  pares no top-5 dos 3 grupos: {n3}")
    print(f"  pares no top-5 de 2 grupos : {n2}")
    print(f"  pares no top-5 de 1 grupo  : {n1}")

    # ── veredito ─────────────────────────────────────────────────────
    print(f"\n=== VEREDITO ===\n")
    max_all_groups = max(summary_per_group[k]["max_abs_r"]
                          for k in summary_per_group)
    mean_all_groups = np.mean([summary_per_group[k]["mean_abs_r"]
                                for k in summary_per_group])
    if max_all_groups < 0.20:
        veredito = ("INTERIOR LIVRE: nenhum grupo apresenta correlação "
                    "interior |r| > 0.20")
    elif max_all_groups < 0.30:
        veredito = ("ESTRUTURA INTERIOR FRACA: max |r| ∈ [0.20, 0.30] em "
                    "algum grupo, mas abaixo do threshold de invariância")
    else:
        veredito = ("INTERIOR TEM ESTRUTURA PRÓPRIA: max |r| > 0.30 em "
                    "ao menos um grupo")
    print(f"  Max |r| interior em qualquer grupo : {max_all_groups:.4f}")
    print(f"  Média |r| interior média dos grupos: {mean_all_groups:.4f}")
    print(f"  Pares Tipo A vs Tipo B nos |r|>0.30:")
    for s in (0, 1, 2):
        cnt_a = sum(1 for absr, r, i, j in
                     top_pairs(corrs[s], top_k=None)
                     if abs(r) > 0.30
                     and edge_pair_type(interior_edges[i],
                                         interior_edges[j])[0] == "A")
        cnt_b = sum(1 for absr, r, i, j in
                     top_pairs(corrs[s], top_k=None)
                     if abs(r) > 0.30
                     and edge_pair_type(interior_edges[i],
                                         interior_edges[j])[0] == "B")
        print(f"    S{s}: {cnt_a} Tipo A, {cnt_b} Tipo B")
    print()
    print(f"  → {veredito}")

    # ── salvar ───────────────────────────────────────────────────────
    out = {
        "board": BOARD,
        "n_unique_cycles": n_unique,
        "interior_edges": [f"{label(e[0])}-{label(e[1])}"
                            for e in interior_edges],
        "group_sizes": {state_desc[s][0]: sizes[s] for s in (0, 1, 2)},
        "summary_per_group": summary_per_group,
        "cross_group_top5_union": rows_for_summary,
        "recurrence": {"n_in_3_groups": n3, "n_in_2": n2, "n_in_1": n1},
        "verdict": veredito,
    }
    here = os.path.dirname(os.path.abspath(__file__))
    out_path = os.path.join(here, "b6_interior_correlations_6x6.json")
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    print(f"\nResultado salvo em: {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
