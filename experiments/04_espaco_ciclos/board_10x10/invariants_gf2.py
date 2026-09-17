#!/usr/bin/env python3
"""
invariants_gf2.py
=================
Descoberta de invariantes do espaço de ciclos a partir das amostras 10×10.

Pipeline:
  2.1  Constrói matriz T (N × E) das assinaturas amostradas
  2.2  Eliminação gaussiana GF(2) → base B (rank × E) do espaço gerado
  2.3  Coordenadas coords = T @ B.T mod 2  (N × rank)
  2.4  Identifica invariantes não-triviais (variância > limite)
  2.5  Correlação ponto-biserial coord_k × aresta_e → candidatos NOT(A∧B)

Entradas:
  data/samples/tours_10x10_batch_*.npy
  data/boundary_matrix_10x10.npy (não usado diretamente, mas confere E)

Saídas:
  data/invariants/h1_basis_sampled.npy
  data/invariants/h1_coordinates.npy
  data/invariants/edge_correlations.npy
  data/invariants/candidate_clauses.json
  data/invariants/summary.json
"""

import argparse
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


def load_all_samples():
    files = sorted(SAMPLES.glob("tours_10x10_batch_*.npy"))
    if not files:
        raise FileNotFoundError("Nenhum batch encontrado em data/samples/")
    arrs = [np.load(f) for f in files]
    T = np.concatenate(arrs, axis=0)
    return T, len(files)


def gf2_row_reduce(M):
    """Eliminação gaussiana GF(2). Retorna (rref, pivots) onde pivots é lista
    dos índices das colunas pivô e rref é a forma reduzida (rank, E)."""
    A = M.copy().astype(np.uint8)
    rows, cols = A.shape
    r = 0
    pivots = []
    for c in range(cols):
        if r >= rows:
            break
        pivot = None
        for rr in range(r, rows):
            if A[rr, c]:
                pivot = rr
                break
        if pivot is None:
            continue
        if pivot != r:
            A[[r, pivot]] = A[[pivot, r]]
        for rr in range(rows):
            if rr != r and A[rr, c]:
                A[rr] ^= A[r]
        pivots.append(c)
        r += 1
    return A[:r], pivots


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--variance-threshold", type=float, default=0.01,
                   help="frac da variância máxima para considerar invariante "
                        "como não-trivial (default 0.01)")
    p.add_argument("--correlation-threshold", type=float, default=0.6,
                   help="r mínimo para correlação aresta×invariante (default 0.6)")
    p.add_argument("--coexist-threshold", type=float, default=0.01,
                   help="P(A=1 ∧ B=1) máximo para candidato NOT(A∧B)")
    args = p.parse_args()

    INV.mkdir(parents=True, exist_ok=True)

    _, edges = build_graph()
    E = len(edges)
    edge_labels = [f"{label(u)}-{label(v)}" for u, v in edges]

    print("Carregando amostras ...")
    T, n_batches = load_all_samples()
    N = T.shape[0]
    print(f"  amostras: {N}  batches: {n_batches}  E={T.shape[1]}")
    assert T.shape[1] == E, "E inconsistente com graph_10x10"

    # ── 2.2 Base do espaço amostrado via row reduction ──────────────────
    print("Eliminação gaussiana GF(2) ...")
    B, pivots = gf2_row_reduce(T)
    rank = B.shape[0]
    print(f"  rank(T)         : {rank}")
    print(f"  β₁ teórico       : {E - 100 + 1}")
    print(f"  pivots          : {len(pivots)} colunas (cols pivô)")

    np.save(INV / "h1_basis_sampled.npy", B)
    print(f"  Salvo: data/invariants/h1_basis_sampled.npy  shape={B.shape}")

    # ── 2.3 Coordenadas: cada amostra → vetor de rank dimensões ─────────
    print("Calculando coordenadas ...")
    # coords[i, k] = T[i] @ B[k].T mod 2 → produto interno GF(2)
    coords = (T.astype(np.uint8) @ B.T.astype(np.uint8)) % 2
    coords = coords.astype(np.uint8)
    np.save(INV / "h1_coordinates.npy", coords)
    print(f"  Salvo: data/invariants/h1_coordinates.npy  shape={coords.shape}")

    # ── 2.4 Invariantes não-triviais ────────────────────────────────────
    # variância de uma variável binária é p(1-p)
    p_inv = coords.mean(axis=0)
    var_inv = p_inv * (1 - p_inv)
    var_max = var_inv.max()
    threshold = args.variance_threshold * var_max

    nontrivial = [k for k in range(rank) if var_inv[k] > threshold]
    print(f"\nInvariantes não-triviais ({len(nontrivial)}/{rank}):")
    if len(nontrivial) <= 20:
        for k in nontrivial:
            n1 = int(coords[:, k].sum())
            n0 = N - n1
            print(f"  k={k:>3d}  p(1)={p_inv[k]:.4f}  n0={n0:>5d}  n1={n1:>5d}")
    else:
        print(f"  (muitos para listar; mostrando primeiros 20)")
        for k in nontrivial[:20]:
            n1 = int(coords[:, k].sum())
            n0 = N - n1
            print(f"  k={k:>3d}  p(1)={p_inv[k]:.4f}  n0={n0:>5d}  n1={n1:>5d}")

    # ── 2.5 Correlação ponto-biserial coord_k × aresta_e ────────────────
    # Para vars binárias, ponto-biserial = Pearson. Usamos np.corrcoef.
    print("\nCorrelando arestas × invariantes ...")
    edge_freq = T.mean(axis=0)
    live_edges = [e for e in range(E)
                  if 0.01 < edge_freq[e] < 0.99]
    print(f"  arestas vivas (0.01 < freq < 0.99): {len(live_edges)}")
    print(f"  arestas obrigatórias candidatas (freq > 0.99): "
          f"{int(np.sum(edge_freq > 0.99))}")

    # corr[k, e] = corr(coords[:, k], T[:, e])
    # cálculo vetorizado, usando float32 para economizar memória
    Cf = coords.astype(np.float32)
    Tf = T.astype(np.float32)
    Cf -= Cf.mean(axis=0, keepdims=True)
    Tf -= Tf.mean(axis=0, keepdims=True)
    std_c = np.linalg.norm(Cf, axis=0)
    std_t = np.linalg.norm(Tf, axis=0)
    std_c[std_c == 0] = 1.0
    std_t[std_t == 0] = 1.0
    corr = (Cf.T @ Tf) / (std_c[:, None] * std_t[None, :])
    np.save(INV / "edge_correlations.npy", corr.astype(np.float32))
    print(f"  Salvo: data/invariants/edge_correlations.npy  shape={corr.shape}")

    # ── candidatos NOT(A∧B): pares de arestas com correlação forte ao mesmo invariante ──
    print(f"\nProcurando candidatos NOT(A∧B) "
          f"(r>{args.correlation_threshold}, P(A∧B)<{args.coexist_threshold}):")
    candidates = []
    for k in nontrivial:
        idx_high = [e for e in live_edges if abs(corr[k, e]) > args.correlation_threshold]
        if len(idx_high) < 2:
            continue
        # para cada par, verifica coexistência empírica
        for i in range(len(idx_high)):
            for j in range(i + 1, len(idx_high)):
                e1, e2 = idx_high[i], idx_high[j]
                coexist = float(((T[:, e1] == 1) & (T[:, e2] == 1)).mean())
                if coexist < args.coexist_threshold:
                    candidates.append({
                        "invariant": int(k),
                        "edge_a_idx": int(e1),
                        "edge_a_label": edge_labels[e1],
                        "edge_b_idx": int(e2),
                        "edge_b_label": edge_labels[e2],
                        "r_a": round(float(corr[k, e1]), 4),
                        "r_b": round(float(corr[k, e2]), 4),
                        "p_coexist": round(coexist, 6),
                        "freq_a": round(float(edge_freq[e1]), 4),
                        "freq_b": round(float(edge_freq[e2]), 4),
                    })

    # dedup por par (a, b) — pode ser detectado por múltiplos invariantes
    seen = set()
    unique_cands = []
    for c in sorted(candidates, key=lambda x: x["p_coexist"]):
        key = tuple(sorted((c["edge_a_idx"], c["edge_b_idx"])))
        if key not in seen:
            seen.add(key)
            unique_cands.append(c)

    print(f"  pares candidatos: {len(unique_cands)} (após dedup)")
    if unique_cands[:15]:
        print("  Top 15 (menor coexistência):")
        for c in unique_cands[:15]:
            print(f"    {c['edge_a_label']:>7s} ↔ {c['edge_b_label']:>7s}  "
                  f"P(∧)={c['p_coexist']:.4f}  "
                  f"freq_a={c['freq_a']:.3f}  freq_b={c['freq_b']:.3f}  "
                  f"r(inv{c['invariant']})=({c['r_a']:+.3f},{c['r_b']:+.3f})")

    with open(INV / "candidate_clauses.json", "w") as f:
        json.dump({
            "n_samples": int(N),
            "n_candidates": len(unique_cands),
            "thresholds": {
                "variance_frac": args.variance_threshold,
                "correlation": args.correlation_threshold,
                "coexist": args.coexist_threshold,
            },
            "candidates": unique_cands,
        }, f, indent=2)
    print(f"\n  Salvo: data/invariants/candidate_clauses.json")

    # ── sumário ─────────────────────────────────────────────────────────
    summary = {
        "n_samples": int(N),
        "n_batches": int(n_batches),
        "E": int(E),
        "beta_1_theoretical": int(E - 100 + 1),
        "rank_sampled": int(rank),
        "rank_coverage_pct": round(rank / (E - 100 + 1) * 100, 2),
        "n_invariants_nontrivial": len(nontrivial),
        "n_candidates_not_and": len(unique_cands),
        "n_edges_mandatory_candidate": int((edge_freq > 0.99).sum()),
        "n_edges_impossible_candidate": int((edge_freq < 0.01).sum()),
        "edge_freq_min": round(float(edge_freq.min()), 4),
        "edge_freq_max": round(float(edge_freq.max()), 4),
    }
    with open(INV / "summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print("\n─" * 30)
    print("RESUMO INVARIANTES")
    print("─" * 30)
    print(f"  N amostras           : {N}")
    print(f"  E arestas            : {E}")
    print(f"  rank(T) amostrado    : {rank}")
    print(f"  β₁ teórico           : {E - 100 + 1}")
    print(f"  cobertura            : {rank / (E - 100 + 1) * 100:.1f}%")
    print(f"  invariantes não-trivs: {len(nontrivial)}")
    print(f"  candidatos NOT(A∧B)  : {len(unique_cands)}")
    print(f"  candidatas mandatory : {(edge_freq > 0.99).sum()}")


if __name__ == "__main__":
    main()
