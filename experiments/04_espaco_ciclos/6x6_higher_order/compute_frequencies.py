#!/usr/bin/env python3
"""
compute_frequencies.py
======================
Constrói a matriz de incidência T ∈ {0,1}^{9862×E} a partir do ground
truth exaustivo (tours_closed_6x6.npy) e calcula:

  - freq[e]       — marginal de cada aresta
  - freq[e1,e2]   — coexistência aresta×aresta (triangular superior)
  - exclusões de pares: P(A∧B)=0 com freq individual > 0

Saídas:
  data/freq_singles.npy
  data/freq_pairs.npy   (matriz simétrica float32 (E,E))
  data/edges_6x6.json   (ordenação canônica das arestas)
  data/incidence_matrix_6x6.npy  (T, uint8 (9862, E))
  data/exclusions_pairs.json
"""

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"

# importa o grafo canônico do projeto 6×6
sys.path.insert(0, str(ROOT.parent))
from cavalo_loop_destruicao_6x6 import (  # noqa: E402
    ADJ, EDGES_LIST, TOTAL, BOARD, vid, label
)

TOURS_PATH = ROOT.parent / "complex_orbit" / "data" / "tours_closed_6x6.npy"


def build_incidence_matrix(tours, edges):
    """T[i, e] = 1 sse a aresta e está no tour i (ciclo fechado)."""
    edge_idx = {e: i for i, e in enumerate(edges)}
    E = len(edges)
    N, n_vertices = tours.shape
    T = np.zeros((N, E), dtype=np.uint8)
    for i, tour in enumerate(tours):
        for k in range(n_vertices):
            u, v = int(tour[k]), int(tour[(k + 1) % n_vertices])
            e = (min(u, v), max(u, v))
            j = edge_idx.get(e)
            if j is not None:
                T[i, j] = 1
    return T


def main():
    DATA.mkdir(parents=True, exist_ok=True)

    print("Carregando ground truth ...")
    tours = np.load(TOURS_PATH).astype(np.int8)
    N = tours.shape[0]
    print(f"  tours: {N}  comprimento: {tours.shape[1]}")
    print(f"  esperado: 9862 tours × 36 vértices")
    assert N == 9862 and tours.shape[1] == 36

    edges = list(EDGES_LIST)
    E = len(edges)
    edge_labels = [f"{label(u)}-{label(v)}" for u, v in edges]
    print(f"\n|E| = {E}  (canônico, lex order)")

    # ── matriz de incidência ─────────────────────────────────────────
    print("\nConstruindo matriz de incidência T ...")
    T = build_incidence_matrix(tours, edges)
    # sanity: cada tour usa exatamente 36 arestas (Hamiltoniano fechado)
    sums = T.sum(axis=1)
    assert (sums == TOTAL).all(), f"Tour com !=36 arestas? sums={sums[:5]}"
    print(f"  shape={T.shape}  cada linha soma {TOTAL} (Hamilton).")

    # ── frequências ──────────────────────────────────────────────────
    freq = T.mean(axis=0).astype(np.float64)
    n_zero = int((freq == 0).sum())
    n_one = int((freq == 1).sum())
    n_live = int(((freq > 0) & (freq < 1)).sum())
    print(f"\nfreq[e]:")
    print(f"  min={freq.min():.4f}  max={freq.max():.4f}")
    print(f"  arestas nunca usadas (freq=0): {n_zero}")
    print(f"  arestas sempre usadas (freq=1, obrigatórias): {n_one}")
    print(f"  arestas vivas (0 < freq < 1): {n_live}")

    # histograma
    print("\nHistograma de freq[e] (10 bins):")
    hist, edges_h = np.histogram(freq, bins=10, range=(0, 1))
    for i in range(10):
        bar = "█" * (hist[i] * 30 // max(hist.max(), 1))
        print(f"  [{edges_h[i]:.1f}, {edges_h[i+1]:.1f}): {hist[i]:>3d}  {bar}")

    np.save(DATA / "freq_singles.npy", freq.astype(np.float32))

    # ── coexistência (E × E) ─────────────────────────────────────────
    print("\nCalculando freq[e1, e2] ...")
    # T.T @ T conta coexistências; dividir por N = probabilidade conjunta
    coex_int = (T.astype(np.int32).T @ T.astype(np.int32))
    coex = coex_int.astype(np.float64) / N
    np.save(DATA / "freq_pairs.npy", coex.astype(np.float32))
    print(f"  shape={coex.shape}  bytes={coex.astype(np.float32).nbytes:,}")

    # ── exclusões de pares ───────────────────────────────────────────
    print("\nIdentificando pares P(A∧B)=0 ...")
    excl_pairs = []
    for i in range(E):
        if freq[i] == 0:
            continue
        for j in range(i + 1, E):
            if freq[j] == 0:
                continue
            if coex[i, j] == 0:
                excl_pairs.append({
                    "edges": [i, j],
                    "edge_names": [edge_labels[i], edge_labels[j]],
                    "freq_individual": [round(float(freq[i]), 6),
                                        round(float(freq[j]), 6)],
                    "p_coexist": 0.0,
                })

    print(f"  pares excluídos (P=0, freq_i > 0, freq_j > 0): {len(excl_pairs)}")

    if excl_pairs[:15]:
        print("\nTop 15 (primeiros por ordem canônica):")
        for p in excl_pairs[:15]:
            print(f"  e{p['edges'][0]:>3d}—e{p['edges'][1]:>3d}  "
                  f"{p['edge_names'][0]:>7s} ↔ {p['edge_names'][1]:>7s}  "
                  f"freq=({p['freq_individual'][0]:.3f}, "
                  f"{p['freq_individual'][1]:.3f})")

    # ── salvar artefatos ─────────────────────────────────────────────
    with open(DATA / "edges_6x6.json", "w") as f:
        json.dump({
            "n_edges": E,
            "ordering": "lex (u,v) com u<v",
            "edges_uv": [[int(u), int(v)] for u, v in edges],
            "edge_labels": edge_labels,
        }, f, indent=2)

    np.save(DATA / "incidence_matrix_6x6.npy", T)

    with open(DATA / "exclusions_pairs.json", "w") as f:
        json.dump({
            "n_tours": N,
            "n_edges": E,
            "n_exclusions": len(excl_pairs),
            "exclusions": excl_pairs,
        }, f, indent=2)

    # ── relatório ────────────────────────────────────────────────────
    print(f"\n{'─' * 65}")
    print("RESUMO TAREFA 0")
    print("─" * 65)
    print(f"  N tours        : {N}")
    print(f"  |E|            : {E}")
    print(f"  obrigatórias   : {n_one}")
    print(f"  vivas          : {n_live}")
    print(f"  nunca usadas   : {n_zero}")
    print(f"  pares excluídos: {len(excl_pairs)}")


if __name__ == "__main__":
    main()
