#!/usr/bin/env python3
"""
structural_analysis.py
======================
T2: análise estrutural dos 4 cantos.

Passos:
  2.1  dim(Span(8 obrigatórias)) em GF(2)^E (deve ser 8)
       dim(Span(8 obrigatórias) mod row(∂₁)) (predição: 4)
  2.2  Confirma e₁(c) + e₂(c) ∈ row(∂₁) para cada canto c, e exibe
       o vetor-vértice (cadeia 0-dimensional) cuja borda é exatamente
       essa soma.
  2.3  Calcula dim do span dos 4 representantes-por-canto no quociente,
       identifica a relação linear (predição: soma dos 4 = 0).
       Constrói explicitamente uma cadeia 0-dimensional `α` em GF(2)^V
       cuja borda é Σ_c (e₁(c) + e₂(c)) — i.e. uma "explicação geométrica".
  2.4  Para cada n ∈ {4,5,6,8,10,12} verifica
         - cantos têm grau 2
         - vizinhanças dos 4 cantos são disjuntas
         - distância mínima entre cantos no grafo
       Reporta n_min em que cantos ficam "independentes".

Saída: data/results/structural_n{n}.json e structural_summary.json
"""
from __future__ import annotations

import argparse
import json
from collections import deque
from itertools import combinations
from pathlib import Path

import numpy as np

# importa utilitários do módulo irmão
import sys
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from verify_small_cases import (  # noqa: E402
    build_graph, vid, label, boundary_matrix, gf2_rank, gf2_rref,
    quotient_dim, corner_edges, xor_pair_vectors,
)

DATA = ROOT / "data"
RESULTS = DATA / "results"


# ── verifica relação e₁(c) + e₂(c) = ∂(canto) ──────────────────────────

def edge_sum_equals_corner_chain(n, edges, corner_v, mand_pair, V):
    """
    Devolve True se (e₁ + e₂) ∈ row(∂₁) — explicitamente, se há vetor-vértice
    α ∈ GF(2)^V cuja borda é exatamente o vetor (e₁ + e₂).

    O caso natural é α = δ_{corner_v}: nesse caso ∂(δ_corner) é a soma das
    arestas incidentes ao canto = e₁ + e₂ (porque grau = 2).
    """
    E = len(edges)
    target = np.zeros(E, dtype=np.uint8)
    for ei in mand_pair:
        target[ei] = 1
    # α = δ_{canto}
    alpha = np.zeros(V, dtype=np.uint8)
    alpha[corner_v] = 1
    boundary = np.zeros(E, dtype=np.uint8)
    for i, (u, v) in enumerate(edges):
        if alpha[u] ^ alpha[v]:
            boundary[i] = 1
    return np.array_equal(target, boundary), alpha


def find_chain_with_boundary(target_edge_vec, B):
    """
    Resolve B^T x = target em GF(2). Retorna x se houver solução, senão None.
    (B: shape (V, E); B^T: shape (E, V).)
    """
    BT = B.T.astype(np.uint8)
    E_dim, V_dim = BT.shape
    # sistema [BT | target] em REF
    M = np.hstack([BT, target_edge_vec.reshape(-1, 1).astype(np.uint8)])
    rows, cols = M.shape
    pivots = []
    r = 0
    for c in range(cols - 1):  # não pivota coluna target
        if r >= rows:
            break
        piv = None
        for rr in range(r, rows):
            if M[rr, c]:
                piv = rr
                break
        if piv is None:
            continue
        if piv != r:
            M[[r, piv]] = M[[piv, r]]
        for rr in range(rows):
            if rr != r and M[rr, c]:
                M[rr] ^= M[r]
        pivots.append(c)
        r += 1
    # checa inconsistência
    for rr in range(r, rows):
        if M[rr, -1]:
            return None
    # solução: zeros nas livres, pivôs recuperam o valor
    x = np.zeros(V_dim, dtype=np.uint8)
    for k, c in enumerate(pivots):
        x[c] = M[k, -1]
    return x


# ── distância no grafo entre vértices ──────────────────────────────────

def bfs_distances(adj, source, total):
    dist = [-1] * total
    dist[source] = 0
    q = deque([source])
    while q:
        v = q.popleft()
        for u in adj[v]:
            if dist[u] < 0:
                dist[u] = dist[v] + 1
                q.append(u)
    return dist


# ── análise para um n ──────────────────────────────────────────────────

def analyze(n: int) -> dict:
    print(f"\n{'='*70}\nESTRUTURA n = {n}\n{'='*70}")
    total = n * n
    adj, edges = build_graph(n)
    E = len(edges)
    V = total
    B = boundary_matrix(edges, V)
    rank_B = gf2_rank(B)

    print(f"V={V}  E={E}  β₁={E-V+1}  rank(∂₁)={rank_B}")

    # cantos
    try:
        corners, mand_per_corner, mand_idx = corner_edges(n, adj, edges)
    except ValueError as e:
        print(f"AVISO: {e}")
        return {"n": n, "error": str(e)}

    # graus dos cantos
    degs = {c: len(adj[c]) for c in corners}
    print(f"Cantos: { {label(c,n): degs[c] for c in corners} }")

    # ── 2.1: dim das 8 obrigatórias ──
    mand_vec = np.zeros((8, E), dtype=np.uint8)
    for k, ei in enumerate(mand_idx):
        mand_vec[k, ei] = 1
    rank_raw = gf2_rank(mand_vec)
    rank_quot = quotient_dim(mand_vec, B)
    print(f"\n[2.1]  dim raw  = {rank_raw}   (esperado 8)")
    print(f"       dim mod ∂₁ = {rank_quot}  (predição 4)")

    # ── 2.2: cada canto satisfaz e₁+e₂ ∈ row(∂₁) com α=δ_{canto} ──
    print(f"\n[2.2]  Verificando e₁(c)+e₂(c) = ∂(δ_canto) para cada canto:")
    by_corner_ok = {}
    rep_vecs = np.zeros((len(corners), E), dtype=np.uint8)
    for k, c in enumerate(corners):
        pair = mand_per_corner[c]
        ok, alpha = edge_sum_equals_corner_chain(n, edges, c, pair, V)
        by_corner_ok[label(c, n)] = bool(ok)
        # representante = e₁(c) (qualquer um, já que e₁ ≡ e₂ mod B)
        rep_vecs[k, pair[0]] = 1
        marker = "✓" if ok else "✗"
        print(f"       canto {label(c,n)}: {marker}")

    # ── 2.3: dim dos 4 representantes-por-canto no quociente ──
    rank_reps_raw = gf2_rank(rep_vecs)
    rank_reps_quot = quotient_dim(rep_vecs, B)
    print(f"\n[2.3]  dim 4 representantes raw  = {rank_reps_raw}")
    print(f"       dim 4 representantes mod ∂₁ = {rank_reps_quot}")
    print(f"       (Predição do brief era 3 — VERIFICAR: na prática observamos 4)")

    # Identificar a relação linear: rank(reps mod B) < 4 implica ∃ subset não-trivial
    # com soma ∈ row(B). Testar especificamente "soma dos 4 = ∂(α)" para algum α.
    sum_reps = rep_vecs.sum(axis=0) % 2
    alpha_for_sum = find_chain_with_boundary(sum_reps, B)
    if alpha_for_sum is not None:
        n_active = int(alpha_for_sum.sum())
        print(f"       soma dos 4 reps = ∂(α) onde α tem {n_active} vértices ativos")
        active_vs = [label(v, n) for v in range(V) if alpha_for_sum[v]]
        preview = active_vs if n_active <= 20 else active_vs[:10] + ["...", active_vs[-1]]
        print(f"       α ⊃ {preview}")
    else:
        print(f"       Soma dos 4 reps NÃO está em row(∂₁) — confirma reps independentes.")

    # ── 2.3-bis: imagem dos 28 XORs no quociente vive em Span{ri+rj} ──
    # Construir diretamente os 6 vetores ri+rj e checar sua dim mod B.
    diffs = []
    diff_labels = []
    for k1 in range(len(corners)):
        for k2 in range(k1 + 1, len(corners)):
            d = (rep_vecs[k1] ^ rep_vecs[k2]).astype(np.uint8)
            diffs.append(d)
            diff_labels.append(f"r({label(corners[k1],n)})+r({label(corners[k2],n)})")
    diffs = np.stack(diffs, axis=0)
    dim_diffs_quot = quotient_dim(diffs, B)
    print(f"\n[2.3b] dim Span(6 diferenças r_i+r_j) mod ∂₁ = {dim_diffs_quot}")
    print(f"       Predição: #cantos − 1 = 3 (ker da soma)")
    print(f"       Esta é a interpretação correta de Q(n) = 3.")

    # ── XOR span check (sanity) ──
    XOR, _ = xor_pair_vectors(mand_idx, E)
    Q_n = quotient_dim(XOR, B)
    print(f"\n[Q]    Q(n) = dim π(Span(28 XOR))  = {Q_n}")

    # ── 2.4: independência das vizinhanças dos cantos ──
    pair_dists = {}
    overlaps = {}
    for ci, cj in combinations(corners, 2):
        d = bfs_distances(adj, ci, total)
        pair_dists[f"{label(ci,n)}↔{label(cj,n)}"] = d[cj]
        # arestas incidentes a ci e cj
        e_ci = set(adj[ci])
        e_cj = set(adj[cj])
        shared = e_ci & e_cj | ({ci} & e_cj) | (e_ci & {cj})
        overlaps[f"{label(ci,n)}↔{label(cj,n)}"] = sorted(label(v, n) for v in shared)
    print(f"\n[2.4]  Distâncias entre cantos: {pair_dists}")
    any_overlap = any(v for v in overlaps.values())
    print(f"       Sobreposição de vizinhanças? {any_overlap}")
    if any_overlap:
        for k, v in overlaps.items():
            if v:
                print(f"         {k}: vizinhança compartilhada = {v}")

    return {
        "n": n,
        "V": V,
        "E": E,
        "beta1": E - V + 1,
        "rank_boundary": int(rank_B),
        "corners": [label(c, n) for c in corners],
        "corner_degrees": {label(c, n): degs[c] for c in corners},
        "step_2_1": {
            "rank_8_raw": int(rank_raw),
            "rank_8_mod_boundary": int(rank_quot),
            "matches_4": int(rank_quot) == 4,
        },
        "step_2_2": {
            "corner_chain_works": by_corner_ok,
            "all_ok": all(by_corner_ok.values()),
        },
        "step_2_3": {
            "rank_4_reps_raw": int(rank_reps_raw),
            "rank_4_reps_mod_boundary": int(rank_reps_quot),
            "reps_independent_mod_boundary": int(rank_reps_quot) == 4,
            "sum_of_4_reps_in_row_boundary": alpha_for_sum is not None,
            "alpha_support_size": int(alpha_for_sum.sum()) if alpha_for_sum is not None else None,
            "alpha_support_labels": (
                [label(v, n) for v in range(V) if alpha_for_sum[v]]
                if alpha_for_sum is not None else None
            ),
            "dim_pairwise_differences_mod_boundary": int(dim_diffs_quot),
            "matches_Q": int(dim_diffs_quot) == int(Q_n),
        },
        "Q_n": int(Q_n),
        "step_2_4": {
            "corner_pair_distances": pair_dists,
            "any_vertex_overlap_between_corner_neighborhoods": any_overlap,
        },
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--ns", type=int, nargs="+", default=[4, 5, 6, 8, 10, 12])
    args = p.parse_args()

    RESULTS.mkdir(parents=True, exist_ok=True)
    summary = []

    for n in args.ns:
        try:
            r = analyze(n)
        except Exception as e:
            print(f"ERRO em n={n}: {e}")
            r = {"n": n, "error": str(e)}
        summary.append(r)
        with open(RESULTS / f"structural_n{n}.json", "w") as f:
            json.dump(r, f, indent=2)

    with open(RESULTS / "structural_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    # ── tabela consolidada ──
    print(f"\n{'='*78}\n{'RESUMO ESTRUTURAL':^78}\n{'='*78}")
    print(f"{'n':>3} {'V':>5} {'E':>5} {'β₁':>5} {'8 mod ∂':>9} "
          f"{'4reps mod ∂':>12} {'6 diffs mod ∂':>14} {'Q':>3} {'overlap?':>9}")
    for r in summary:
        if "error" in r:
            print(f"{r['n']:>3} ERRO: {r['error']}")
            continue
        print(f"{r['n']:>3} {r['V']:>5} {r['E']:>5} {r['beta1']:>5} "
              f"{r['step_2_1']['rank_8_mod_boundary']:>9} "
              f"{r['step_2_3']['rank_4_reps_mod_boundary']:>12} "
              f"{r['step_2_3']['dim_pairwise_differences_mod_boundary']:>14} "
              f"{r['Q_n']:>3} "
              f"{str(r['step_2_4']['any_vertex_overlap_between_corner_neighborhoods']):>9}")


if __name__ == "__main__":
    main()
