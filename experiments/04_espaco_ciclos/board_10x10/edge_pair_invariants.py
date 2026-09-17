#!/usr/bin/env python3
"""
edge_pair_invariants.py
=======================
Análise direta aresta×aresta para encontrar bifurcações topológicas.

Replica o achado central do 6×6 (pares com r ≈ -0.77) sem depender da
escolha arbitrária de base GF(2). Para cada par de arestas vivas (e1, e2):
  - freq_a, freq_b
  - p_coexist = P(x_a=1 ∧ x_b=1)
  - corr_pearson = (E[XY] - E[X]E[Y]) / sqrt(var_a * var_b)
                 = (p_coexist - p_a * p_b) / sqrt(p_a(1-p_a) * p_b(1-p_b))

Filtra para candidatos NOT(A∧B):
  - corr < -0.5 (exclusão mútua)
  - ambas freq vivas (0.01 < freq < 0.99)

Saída:
  data/invariants/edge_pair_correlations.json
  data/invariants/candidate_clauses.json  (sobrescreve com base aresta×aresta)
"""

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
SAMPLES = DATA / "samples"
INV = DATA / "invariants"

sys.path.insert(0, str(ROOT))
from graph_10x10 import build_graph, label  # noqa: E402


def main():
    INV.mkdir(parents=True, exist_ok=True)

    _, edges = build_graph()
    E = len(edges)
    edge_labels = [f"{label(u)}-{label(v)}" for u, v in edges]

    files = sorted(SAMPLES.glob("tours_10x10_batch_*.npy"))
    T = np.concatenate([np.load(f) for f in files], axis=0).astype(np.uint8)
    N = T.shape[0]
    print(f"Amostras: {N}  E={E}")

    freq = T.mean(axis=0)
    live = np.where((freq > 0.01) & (freq < 0.99))[0]
    print(f"Arestas vivas: {len(live)}")

    # matriz de coexistência: P(A∧B) = (T.T @ T) / N para arestas vivas
    Tl = T[:, live].astype(np.float32)
    coex = (Tl.T @ Tl) / N  # (n_live, n_live)
    pl = freq[live].astype(np.float32)
    pa_pb = np.outer(pl, pl)
    var_a = pl * (1 - pl)
    denom = np.sqrt(np.outer(var_a, var_a))
    denom[denom == 0] = 1.0
    corr = (coex - pa_pb) / denom

    # filtra par (i<j) com corr < -0.5
    candidates = []
    for ii in range(len(live)):
        for jj in range(ii + 1, len(live)):
            r = float(corr[ii, jj])
            if r < -0.5:
                e1 = int(live[ii])
                e2 = int(live[jj])
                candidates.append({
                    "edge_a_idx": e1,
                    "edge_b_idx": e2,
                    "edge_a_label": edge_labels[e1],
                    "edge_b_label": edge_labels[e2],
                    "r": round(r, 4),
                    "p_coexist": round(float(coex[ii, jj]), 6),
                    "freq_a": round(float(pl[ii]), 4),
                    "freq_b": round(float(pl[jj]), 4),
                })

    candidates.sort(key=lambda c: c["r"])
    print(f"\nPares com r < -0.5: {len(candidates)}")
    if candidates:
        print("\nTop 20:")
        print(f"  {'edge_a':>9s} {'edge_b':>9s}  {'r':>8s}  {'P(AB)':>8s}  "
              f"{'freq_a':>7s} {'freq_b':>7s}")
        for c in candidates[:20]:
            print(f"  {c['edge_a_label']:>9s} {c['edge_b_label']:>9s}  "
                  f"{c['r']:+8.4f}  {c['p_coexist']:>8.4f}  "
                  f"{c['freq_a']:>7.4f} {c['freq_b']:>7.4f}")

    # distribuição de correlações
    upper = corr[np.triu_indices_from(corr, k=1)]
    bins = [-1, -0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 0.7, 1.0]
    hist, _ = np.histogram(upper, bins=bins)
    print("\nDistribuição de correlações aresta-aresta (triangle superior):")
    for i in range(len(hist)):
        print(f"  [{bins[i]:+.1f}, {bins[i+1]:+.1f}): {hist[i]:>6d}")

    with open(INV / "edge_pair_correlations.json", "w") as f:
        json.dump({
            "n_samples": int(N),
            "n_live_edges": int(len(live)),
            "n_negative_pairs": len(candidates),
            "candidates": candidates,
            "histogram_bins": bins,
            "histogram_counts": hist.tolist(),
        }, f, indent=2)
    print(f"\nSalvo: data/invariants/edge_pair_correlations.json")

    # sobrescreve candidate_clauses.json com a versão direta — mais útil p/ benchmark
    with open(INV / "candidate_clauses.json", "w") as f:
        json.dump({
            "n_samples": int(N),
            "n_candidates": len(candidates),
            "method": "edge_pair_direct",
            "thresholds": {"correlation_lt": -0.5},
            "candidates": candidates,
        }, f, indent=2)
    print(f"Salvo: data/invariants/candidate_clauses.json  (sobrescrito com par direto)")


if __name__ == "__main__":
    main()
