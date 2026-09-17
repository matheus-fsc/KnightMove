#!/usr/bin/env python3
"""
bifurcation_independence.py
============================
Teste de independência entre as 13 bifurcações canônicas da Fase B.

Cada bifurcação é um par de arestas (A, B) fortemente anti-correlacionadas.
Para uma amostra Z3 ser "excitada" pela bifurcação, exatamente uma das duas
arestas deve estar ativa (XOR). A "escolha" da bifurcação é encodificada
como 0 (A ativa) ou 1 (B ativa).

A independência é avaliada pela matriz de correlação Pearson entre as 13
"escolhas" binárias, restrita às amostras onde TODAS as bifurcações estão
excitadas simultaneamente (máscara XOR global, conforme especificação).

Saídas:
  1. Matriz de correlação 13×13 completa
  2. Máximo |r| fora da diagonal
  3. Média |r| fora da diagonal
  4. Pares com |r| > 0.15
  5. Veredito: INDEPENDENTES / QUASE / DEPENDENTES

Notas:
  - Bifurcações que compartilham vértice podem colapsar (variância 0) na
    máscara global; reportamos isso explicitamente.
  - Como complemento, também rodamos o teste em modo "pairwise": para cada
    par (i, j), usa-se a máscara XOR restrita às bifurcações i e j.
"""

import argparse
import glob
import json
import os
import sys
from itertools import combinations

import numpy as np

import d4_orbits
from correlator import _build_d4_remaps


# Top-13 canônicos da Fase B (presentes em 18/18 pools).
# Ordenados por r_min global.
BIFURCATIONS = [
    ("B8-D7", "B8-C6"),   # B00: r=-0.7698
    ("C8-E7", "C8-D6"),   # B01: r=-0.5430
    ("D8-F7", "D8-E6"),   # B02: r=-0.4792
    ("B8-A6", "B8-C6"),   # B03: r=-0.4787
    ("C8-E7", "C8-B6"),   # B04: r=-0.4551
    ("C8-B6", "B6-A4"),   # B05: r=-0.4646
    ("D8-B7", "B7-C5"),   # B06: r=-0.4240
    ("D8-B7", "B7-D6"),   # B07: r=-0.4327
    ("C8-E7", "E7-D5"),   # B08: r=-0.3908
    ("D8-C6", "D8-E6"),   # B09: r=-0.3865
    ("D8-F7", "F7-E5"),   # B10: r=-0.3554
    ("D8-B7", "D8-C6"),   # B11: r=-0.3405
    ("A8-C7", "A8-B6"),   # B12: r=-0.3333
]


def load_batches_expanded(samples_dir, pair_idx, max_pre_d4=None):
    paths = sorted(glob.glob(
        os.path.join(samples_dir, f"batch_p{pair_idx:02d}_*.json")))
    if not paths:
        return None, None
    with open(paths[0]) as f:
        first = json.load(f)
    edges = first["edges"]
    remaps = _build_d4_remaps(edges)
    parts = []
    n_pre = 0
    for path in paths:
        if max_pre_d4 is not None and n_pre >= max_pre_d4:
            break
        with open(path) as f:
            b = json.load(f)
        if b["edges"] != edges:
            raise ValueError(f"ordem inconsistente em {path}")
        X = np.asarray(b["signatures"], dtype=np.bool_)
        if max_pre_d4 is not None:
            take = min(X.shape[0], max_pre_d4 - n_pre)
            X = X[:take]
        for t in range(8):
            Y = np.empty_like(X)
            Y[:, remaps[t]] = X
            parts.append(Y)
        n_pre += X.shape[0]
    Xall = np.concatenate(parts, axis=0)
    return edges, Xall


def excitation_stats(X, edges, bifurcations):
    edge_idx = {e: i for i, e in enumerate(edges)}
    N = X.shape[0]
    stats = []
    for k, (a_lbl, b_lbl) in enumerate(bifurcations):
        a = X[:, edge_idx[a_lbl]]
        b = X[:, edge_idx[b_lbl]]
        n_both = int((a & b).sum())
        n_none = int((~a & ~b).sum())
        n_xor = int((a ^ b).sum())
        n_b_only = int((~a & b).sum())
        stats.append({
            "idx": k, "a": a_lbl, "b": b_lbl,
            "p_a": float(a.mean()), "p_b": float(b.mean()),
            "n_xor": n_xor, "rate_xor": n_xor / N,
            "n_both": n_both, "n_neither": n_none,
            "p_b_given_xor": (n_b_only / n_xor) if n_xor else 0.0,
        })
    return stats


def build_bifurcation_matrix_global(X, edges, bifurcations, verbose=True):
    """Z[n, k] = 0/1 (escolha de bifurcação k) com máscara XOR GLOBAL.

    Conforme spec do usuário: descarta amostras onde QUALQUER bifurcação
    não está excitada (ambas ativas ou ambas inativas).
    """
    edge_idx = {e: i for i, e in enumerate(edges)}
    N = X.shape[0]
    K = len(bifurcations)
    cols = np.zeros((N, K), dtype=np.int8)
    valid = np.ones(N, dtype=bool)
    for k, (a_lbl, b_lbl) in enumerate(bifurcations):
        a = X[:, edge_idx[a_lbl]]
        b = X[:, edge_idx[b_lbl]]
        valid &= (a ^ b)
        cols[:, k] = b.astype(np.int8)
    Z = cols[valid]
    if verbose:
        print(f"  N original (pós-D₄)        : {N:,}")
        print(f"  N válido (XOR todas excit.) : {len(Z):,}  ({len(Z)/N:6.2%})")
    return Z


def correlation_matrix(Z, eps=1e-12):
    """Pearson, tolerando variância zero (retorna NaN onde indefinido)."""
    Zf = Z.astype(np.float64)
    n = Zf.shape[0]
    mu = Zf.mean(axis=0)
    sd = Zf.std(axis=0)
    constant = (sd < eps)
    sd_safe = np.where(constant, 1.0, sd)
    Zn = (Zf - mu) / sd_safe
    C = (Zn.T @ Zn) / n
    # NaN nas linhas/colunas constantes
    for i, c in enumerate(constant):
        if c:
            C[i, :] = np.nan
            C[:, i] = np.nan
            C[i, i] = np.nan
    return C, constant


def pairwise_correlation_matrix(X, edges, bifurcations):
    """Para cada par (i, j), restringe à máscara XOR conjunta {i, j} e
    computa Pearson. Diagonal = 1, células sem dados = NaN."""
    edge_idx = {e: i for i, e in enumerate(edges)}
    K = len(bifurcations)
    C = np.full((K, K), np.nan, dtype=np.float64)
    n_used = np.zeros((K, K), dtype=np.int64)

    bif_a = []
    bif_b = []
    bif_xor = []
    for a_lbl, b_lbl in bifurcations:
        a = X[:, edge_idx[a_lbl]]
        b = X[:, edge_idx[b_lbl]]
        bif_a.append(a)
        bif_b.append(b)
        bif_xor.append(a ^ b)

    for i in range(K):
        C[i, i] = 1.0
    for i in range(K):
        for j in range(i + 1, K):
            mask = bif_xor[i] & bif_xor[j]
            n = int(mask.sum())
            n_used[i, j] = n_used[j, i] = n
            if n < 30:
                continue
            zi = bif_b[i][mask].astype(np.float64)
            zj = bif_b[j][mask].astype(np.float64)
            si, sj = zi.std(), zj.std()
            if si < 1e-12 or sj < 1e-12:
                # ao menos um lado é constante na interseção → r=0
                # (independência operacional sobre o suporte conjunto)
                C[i, j] = C[j, i] = 0.0
                continue
            mi, mj = zi.mean(), zj.mean()
            r = ((zi - mi) * (zj - mj)).mean() / (si * sj)
            C[i, j] = C[j, i] = r
    return C, n_used


def print_matrix(C, bifurcations, title):
    print(f"\n{title}")
    print("        " + "  ".join(f"B{j:02d}" for j in range(len(bifurcations))))
    for i in range(len(bifurcations)):
        cells = []
        for j in range(len(bifurcations)):
            v = C[i, j]
            if np.isnan(v):
                cells.append(" nan")
            else:
                cells.append(f"{v:+.2f}")
        print(f"  B{i:02d}: " + "  ".join(cells))


def summarize(C, bifurcations, label):
    K = len(bifurcations)
    off_vals = []
    for i in range(K):
        for j in range(i + 1, K):
            if not np.isnan(C[i, j]):
                off_vals.append((float(C[i, j]), i, j))
    if not off_vals:
        print(f"\n{label}: nenhum par com correlação definida")
        return None
    abs_vals = [abs(v) for v, _, _ in off_vals]
    max_abs = max(abs_vals)
    mean_abs = float(np.mean(abs_vals))
    above_15 = sorted([(v, i, j) for v, i, j in off_vals if abs(v) > 0.15],
                      key=lambda x: -abs(x[0]))
    sorted_top = sorted(off_vals, key=lambda x: -abs(x[0]))[:15]
    print(f"\n--- {label} ---")
    print(f"  pares com correlação definida : {len(off_vals)} / {K*(K-1)//2}")
    print(f"  máximo |r| fora da diagonal   : {max_abs:.4f}")
    print(f"  média  |r| fora da diagonal   : {mean_abs:.4f}")
    print(f"  pares com |r| > 0.15          : {len(above_15)}")
    print(f"\n  top-15 |r| fora da diagonal:")
    for r, i, j in sorted_top:
        ai, bi = bifurcations[i]
        aj, bj = bifurcations[j]
        flag = "  <-- |r|>0.15" if abs(r) > 0.15 else ""
        print(f"    B{i:02d} ({ai:>5s}↔{bi:<5s})  ↔  "
              f"B{j:02d} ({aj:>5s}↔{bj:<5s})  : r={r:+.4f}{flag}")
    if above_15:
        print(f"\n  todos pares com |r| > 0.15 ({len(above_15)}):")
        for r, i, j in above_15:
            print(f"    B{i:02d} ↔ B{j:02d}  r={r:+.4f}  "
                  f"[{bifurcations[i][0]}↔{bifurcations[i][1]}  vs  "
                  f"{bifurcations[j][0]}↔{bifurcations[j][1]}]")
    if max_abs < 0.10:
        verdict = "INDEPENDENTES (espaço de soluções ≈ hipercubo 2^K)"
    elif max_abs < 0.25:
        verdict = "QUASE INDEPENDENTES (dependências fracas existem)"
    else:
        verdict = "DEPENDENTES (acoplamento significativo)"
    print(f"\n  veredito: {verdict}")
    return {
        "max_abs": max_abs, "mean_abs": mean_abs,
        "n_above_0.15": len(above_15),
        "pairs_above_0.15": [{"i": i, "j": j, "r": r,
                              "bif_i": list(bifurcations[i]),
                              "bif_j": list(bifurcations[j])}
                             for r, i, j in above_15],
        "verdict": verdict,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="data_8x8_phase_b")
    ap.add_argument("--pair-idx", type=int, default=2,
                    help="par do config (default 2 = (0,0)→(0,7))")
    ap.add_argument("--max-pre-d4", type=int, default=None,
                    help="limita amostras pré-D₄ (default = tudo)")
    ap.add_argument("--board", type=int, default=8)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    d4_orbits.set_board(args.board)

    samples_dir = os.path.join(args.data_dir, "samples")
    print(f"Carregando + expandindo D₄ de {samples_dir} pair_idx={args.pair_idx}")
    edges, X = load_batches_expanded(samples_dir, args.pair_idx, args.max_pre_d4)
    if X is None:
        print(f"Nenhum batch encontrado.")
        return 1
    print(f"  X.shape = {X.shape}  ({X.shape[0]:,} amostras × {X.shape[1]} arestas)")

    print("\nEstatísticas de excitação por bifurcação:")
    print(f"  {'idx':<4} {'A':>6s} {'B':>6s}  {'P(A)':>6s} {'P(B)':>6s}  "
          f"{'XOR rate':>9s}  {'P(B|XOR)':>9s}")
    stats = excitation_stats(X, edges, BIFURCATIONS)
    for s in stats:
        print(f"  B{s['idx']:02d}   {s['a']:>5s}  {s['b']:>5s}    "
              f"{s['p_a']:.3f}  {s['p_b']:.3f}     "
              f"{s['rate_xor']:6.2%}     {s['p_b_given_xor']:6.3f}")

    print("\n=== TESTE 1: máscara XOR GLOBAL (especificação do usuário) ===")
    Z = build_bifurcation_matrix_global(X, edges, BIFURCATIONS)
    n_valid = Z.shape[0]
    if n_valid < 30:
        print(f"  ATENÇÃO: apenas {n_valid} amostras válidas — teste global "
              f"inviável.")
        C_global = np.full((13, 13), np.nan)
        const_global = np.ones(13, dtype=bool)
        summary_global = None
    else:
        C_global, const_global = correlation_matrix(Z)
        print(f"  bifurcações com variância 0 (constantes na intersecção): "
              f"{list(np.where(const_global)[0])}")
        print_matrix(C_global, BIFURCATIONS,
                     "Matriz de correlação 13×13 (XOR global):")
        summary_global = summarize(C_global, BIFURCATIONS, "GLOBAL")

    print("\n\n=== TESTE 2: máscara XOR PAIRWISE (mais informativo) ===")
    print("  Para cada par (i, j) usa-se a máscara XOR_i ∩ XOR_j separadamente.")
    C_pair, n_used = pairwise_correlation_matrix(X, edges, BIFURCATIONS)
    print(f"  amostras usadas por par (min/median/max): "
          f"{int(np.nanmin(n_used[n_used>0]) if (n_used>0).any() else 0):,} / "
          f"{int(np.median(n_used[n_used>0]) if (n_used>0).any() else 0):,} / "
          f"{int(np.nanmax(n_used) if (n_used>0).any() else 0):,}")
    print_matrix(C_pair, BIFURCATIONS,
                 "Matriz de correlação 13×13 (XOR pairwise):")
    summary_pair = summarize(C_pair, BIFURCATIONS, "PAIRWISE")

    if args.out:
        out = {
            "data_dir": args.data_dir,
            "pair_idx": args.pair_idx,
            "n_samples_pos_d4": int(X.shape[0]),
            "bifurcations": [list(b) for b in BIFURCATIONS],
            "excitation_stats": stats,
            "global": {
                "n_valid": int(n_valid),
                "constants": [int(i) for i in np.where(const_global)[0]],
                "correlation_matrix": [[None if np.isnan(v) else float(v)
                                        for v in row] for row in C_global],
                "summary": summary_global,
            },
            "pairwise": {
                "n_used_matrix": n_used.tolist(),
                "correlation_matrix": [[None if np.isnan(v) else float(v)
                                        for v in row] for row in C_pair],
                "summary": summary_pair,
            },
        }
        with open(args.out, "w") as f:
            json.dump(out, f, indent=2, ensure_ascii=False)
        print(f"\nResultado JSON salvo em: {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
