#!/usr/bin/env python3
"""
vertex_state_independence.py
============================
Reduz cada um dos 4 clusters identificados no teste de bifurcações (A8, B8,
C8, D8 — vértices da linha superior do tabuleiro 8×8) a uma única variável
categórica: o subconjunto de arestas-cavalo incidentes ativas no caminho
(estado do vértice).

Para cada par de vértices (v_i, v_j), calcula Cramér's V (estatística de
associação categórica baseada em chi-quadrado de Pearson, ∈ [0, 1]; V = 0
sob independência, V = 1 sob dependência funcional).

Também reporta entropia, informação mútua normalizada (NMI) como métrica
complementar.

Threshold de veredito (análogo a |r| < 0.15 do teste de bifurcações):
  V < 0.10  → INDEPENDENTES
  V < 0.20  → QUASE INDEPENDENTES
  V ≥ 0.20  → DEPENDENTES (mapear acoplamento)
"""

import argparse
import glob
import json
import os
import sys

import numpy as np

import d4_orbits
from correlator import _build_d4_remaps


VERTICES = ["A8", "B8", "C8", "D8"]


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


def incident_edges(vertex_label, edges):
    """Lista as arestas (e índices, e label-do-outro-vértice) incidentes."""
    out = []
    for i, e in enumerate(edges):
        a, b = e.split("-")
        if a == vertex_label:
            out.append((i, e, b))
        elif b == vertex_label:
            out.append((i, e, a))
    out.sort(key=lambda x: x[2])  # ordem estável por vizinho
    return out


def encode_vertex_state(X, edges, vertex_label):
    """Codifica cada amostra como inteiro = bitmask das arestas ativas
    incidentes ao vértice (k bits onde k = grau do vértice)."""
    inc = incident_edges(vertex_label, edges)
    bits = X[:, [i for i, _, _ in inc]].astype(np.int32)
    state = np.zeros(X.shape[0], dtype=np.int32)
    for j in range(bits.shape[1]):
        state |= (bits[:, j] << j)
    return state, inc


def decode_state_label(state_int, inc):
    """Retorna string descritiva: 'AB-CD,EF-GH' das arestas ativas."""
    parts = []
    for j, (_, edge, neigh) in enumerate(inc):
        if state_int & (1 << j):
            parts.append(neigh)
    if not parts:
        return "∅"
    return "{" + ",".join(parts) + "}"


def cramers_v(x, y):
    """Cramér's V entre categóricos. Retorna (V, chi2, df, n, n_rows, n_cols).

    Implementação manual (sem scipy): constrói contingência, calcula chi²
    de Pearson, e V = sqrt(chi² / (n * min(r-1, c-1))).
    """
    cats_x = np.unique(x)
    cats_y = np.unique(y)
    r, c = len(cats_x), len(cats_y)
    table = np.zeros((r, c), dtype=np.int64)
    for i, vx in enumerate(cats_x):
        mask_x = (x == vx)
        for j, vy in enumerate(cats_y):
            table[i, j] = int((mask_x & (y == vy)).sum())
    n = int(table.sum())
    row_sum = table.sum(axis=1, keepdims=True)
    col_sum = table.sum(axis=0, keepdims=True)
    expected = (row_sum * col_sum) / n
    with np.errstate(divide="ignore", invalid="ignore"):
        contrib = np.where(expected > 0,
                           (table - expected) ** 2 / expected, 0.0)
    chi2 = float(contrib.sum())
    denom = n * min(r - 1, c - 1) if min(r - 1, c - 1) > 0 else 1
    V = float(np.sqrt(chi2 / denom)) if denom > 0 else 0.0
    df = (r - 1) * (c - 1)
    return V, chi2, df, n, r, c, table, cats_x, cats_y


def normalized_mutual_info(x, y):
    """NMI(x, y) = 2 I(x; y) / (H(x) + H(y)), ∈ [0, 1]."""
    cats_x, cnt_x = np.unique(x, return_counts=True)
    cats_y, cnt_y = np.unique(y, return_counts=True)
    n = len(x)
    px = cnt_x / n
    py = cnt_y / n
    # joint
    r, c = len(cats_x), len(cats_y)
    table = np.zeros((r, c), dtype=np.int64)
    for i, vx in enumerate(cats_x):
        mask_x = (x == vx)
        for j, vy in enumerate(cats_y):
            table[i, j] = int((mask_x & (y == vy)).sum())
    pxy = table / n
    Hx = -float(np.sum(px * np.log2(np.maximum(px, 1e-300))))
    Hy = -float(np.sum(py * np.log2(np.maximum(py, 1e-300))))
    pxy_safe = np.maximum(pxy, 1e-300)
    px_col = px[:, None]
    py_row = py[None, :]
    log_term = np.log2(pxy_safe / np.maximum(px_col * py_row, 1e-300))
    I = float(np.sum(pxy * log_term))
    nmi = (2 * I / (Hx + Hy)) if (Hx + Hy) > 0 else 0.0
    return nmi, I, Hx, Hy


def report_vertex_marginals(states_by_v, inc_by_v):
    print("\n=== Distribuições marginais dos estados de vértice ===")
    for v in VERTICES:
        s = states_by_v[v]
        inc = inc_by_v[v]
        cats, cnt = np.unique(s, return_counts=True)
        N = len(s)
        print(f"\n  {v} (grau={len(inc)}, vizinhos={[x[2] for x in inc]}):")
        # entropia
        p = cnt / N
        H = -float(np.sum(p * np.log2(np.maximum(p, 1e-300))))
        print(f"    {len(cats)} estados distintos, H = {H:.3f} bits "
              f"(máx possível = {np.log2(len(cats)):.3f} se uniforme)")
        for cat, n in zip(cats, cnt):
            label = decode_state_label(int(cat), inc)
            print(f"      {label:>20s} : {n:>9,} ({n/N:6.2%})")


def report_pairwise(states_by_v, inc_by_v):
    print("\n=== Matriz Cramér's V (vértice × vértice) ===")
    K = len(VERTICES)
    V_mat = np.zeros((K, K), dtype=np.float64)
    NMI_mat = np.zeros((K, K), dtype=np.float64)
    details = {}
    for i in range(K):
        V_mat[i, i] = 1.0
        NMI_mat[i, i] = 1.0
    for i in range(K):
        for j in range(i + 1, K):
            vi, vj = VERTICES[i], VERTICES[j]
            V, chi2, df, n, r, c, table, cats_x, cats_y = cramers_v(
                states_by_v[vi], states_by_v[vj])
            nmi, mi, Hx, Hy = normalized_mutual_info(
                states_by_v[vi], states_by_v[vj])
            V_mat[i, j] = V_mat[j, i] = V
            NMI_mat[i, j] = NMI_mat[j, i] = nmi
            details[(vi, vj)] = {
                "cramers_v": V, "nmi": nmi,
                "mutual_info_bits": mi, "Hx": Hx, "Hy": Hy,
                "chi2": chi2, "df": df, "n": n,
                "table_shape": [int(r), int(c)],
                "table": table.tolist(),
                "cats_x": [int(x) for x in cats_x],
                "cats_y": [int(y) for y in cats_y],
            }

    print(f"\n  Cramér's V (V = 0 ⇒ independentes, V = 1 ⇒ funcional):\n")
    print(f"          " + "    ".join(f"{v:>4s}" for v in VERTICES))
    for i, vi in enumerate(VERTICES):
        cells = "    ".join(f"{V_mat[i,j]:>4.2f}" for j in range(K))
        print(f"  {vi:>4s}:   " + cells)

    print(f"\n  NMI (informação mútua normalizada, ∈ [0,1]):\n")
    print(f"          " + "    ".join(f"{v:>4s}" for v in VERTICES))
    for i, vi in enumerate(VERTICES):
        cells = "    ".join(f"{NMI_mat[i,j]:>4.2f}" for j in range(K))
        print(f"  {vi:>4s}:   " + cells)

    print(f"\n  Detalhe de cada par:")
    for (vi, vj), d in details.items():
        print(f"\n    {vi} × {vj}:  V={d['cramers_v']:.4f}  "
              f"NMI={d['nmi']:.4f}  MI={d['mutual_info_bits']:.4f} bits  "
              f"χ²={d['chi2']:.1f} (df={d['df']}, n={d['n']:,})")
        # tabela de contingência decodificada
        table = np.array(d["table"])
        cats_x = d["cats_x"]
        cats_y = d["cats_y"]
        inc_x = inc_by_v[vi]
        inc_y = inc_by_v[vj]
        labels_x = [decode_state_label(c, inc_x) for c in cats_x]
        labels_y = [decode_state_label(c, inc_y) for c in cats_y]
        # imprime tabela
        col_hdr = "  ".join(f"{l:>14s}" for l in labels_y)
        print(f"        {' ':<20s}  {col_hdr}")
        for r_i, row in enumerate(table):
            cells = "  ".join(f"{c:>14,}" for c in row)
            print(f"        {labels_x[r_i]:<20s}  {cells}")

    return V_mat, NMI_mat, details


def verdict(V_mat):
    K = V_mat.shape[0]
    off = np.array([V_mat[i, j] for i in range(K) for j in range(i + 1, K)])
    max_v = float(off.max())
    mean_v = float(off.mean())
    print(f"\n=== VEREDITO ===")
    print(f"  Máximo V fora da diagonal : {max_v:.4f}")
    print(f"  Média  V fora da diagonal : {mean_v:.4f}")
    if max_v < 0.10:
        verd = "INDEPENDENTES — os 4 vértices não compartilham informação"
    elif max_v < 0.20:
        verd = "QUASE INDEPENDENTES — acoplamentos fracos entre vértices"
    else:
        verd = "DEPENDENTES — vértices têm acoplamento significativo"
    print(f"  Veredito: {verd}")
    return max_v, mean_v, verd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="data_8x8_phase_b")
    ap.add_argument("--pair-idx", type=int, default=2)
    ap.add_argument("--max-pre-d4", type=int, default=None)
    ap.add_argument("--board", type=int, default=8)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    d4_orbits.set_board(args.board)

    samples_dir = os.path.join(args.data_dir, "samples")
    print(f"Carregando + D₄ de {samples_dir} pair_idx={args.pair_idx}")
    edges, X = load_batches_expanded(samples_dir, args.pair_idx, args.max_pre_d4)
    if X is None:
        print(f"Nenhum batch encontrado.")
        return 1
    print(f"  X.shape = {X.shape}")

    states_by_v = {}
    inc_by_v = {}
    for v in VERTICES:
        s, inc = encode_vertex_state(X, edges, v)
        states_by_v[v] = s
        inc_by_v[v] = inc

    report_vertex_marginals(states_by_v, inc_by_v)
    V_mat, NMI_mat, details = report_pairwise(states_by_v, inc_by_v)
    max_v, mean_v, verd = verdict(V_mat)

    if args.out:
        out = {
            "data_dir": args.data_dir,
            "pair_idx": args.pair_idx,
            "n_samples_pos_d4": int(X.shape[0]),
            "vertices": VERTICES,
            "vertex_neighbors": {v: [n for _, _, n in inc_by_v[v]]
                                  for v in VERTICES},
            "marginals": {
                v: {
                    "n_states": int(len(np.unique(states_by_v[v]))),
                    "entropy_bits": float(
                        -np.sum(np.unique(states_by_v[v], return_counts=True)[1]
                                / len(states_by_v[v])
                                * np.log2(np.maximum(
                                    np.unique(states_by_v[v],
                                              return_counts=True)[1]
                                    / len(states_by_v[v]), 1e-300)))),
                } for v in VERTICES
            },
            "cramers_v_matrix": V_mat.tolist(),
            "nmi_matrix": NMI_mat.tolist(),
            "max_v_off_diag": max_v,
            "mean_v_off_diag": mean_v,
            "verdict": verd,
            "pairwise_details": {
                f"{vi}_{vj}": d for (vi, vj), d in details.items()
            },
        }
        with open(args.out, "w") as f:
            json.dump(out, f, indent=2, ensure_ascii=False)
        print(f"\nResultado salvo em: {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
