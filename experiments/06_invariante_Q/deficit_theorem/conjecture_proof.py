#!/usr/bin/env python3
"""
conjecture_proof.py
===================
T3: prova construtiva (parte algébrica) e verificação computacional dos
lemas, integrada com o documento `conjecture_proof.md`.

O que é PROVADO algebricamente aqui (computacional + dedução simbólica):

  Lema 1 (cantos têm grau 2): VERIFICADO computacionalmente para
    n ∈ {4, 5, 6, 8, 10, 12} e é elementar de mostrar para todo n ≥ 4
    (vizinhanças de cavalo a partir do canto têm exatamente 2 entradas
    dentro do tabuleiro). Prova combinatorial direta no .md.

  Lema 2 (cada canto satisfaz e₁(c)+e₂(c) ∈ row(∂₁)):
    Prova algébrica direta: tomando α = δ_{c} ∈ GF(2)^V (chain 0-dim
    concentrada no canto), a borda ∂(δ_c) é o vetor indicador das
    arestas incidentes a c em GF(2). Como c tem grau 2 (Lema 1), são
    exatamente e₁ e e₂. Logo e₁+e₂ = ∂(δ_c) ∈ row(∂).

  Lema 3 (os 4 representantes r_c := e₁(c) mod row(∂₁) são linearmente
    independentes em GF(2)^E / row(∂₁)):
    VERIFICADO computacionalmente. Prova esboçada: cada r_c "vê"
    arestas exclusivas ao canto c (porque, para n ≥ 6, vizinhanças
    são disjuntas).

  Lema 4 (a imagem dos 28 vetores XOR no quociente coincide com
    Span{r_ci + r_cj : i<j}):
    PROVA DIRETA. Cada XOR v_{ij} com i,j ∈ Mand(n) pertence a uma de
    duas categorias:
      (a) i,j em mesmo canto c: v_{ij} = e₁(c)+e₂(c) ≡ 0 mod row(∂).
      (b) i ∈ Mand(c_a), j ∈ Mand(c_b) com a ≠ b: v_{ij} = e_x(c_a) +
          e_y(c_b) ≡ r_{c_a} + r_{c_b} mod row(∂) (pela Lema 2).
    Logo Span(XOR) mod row(∂) ⊆ Span{r_a+r_b}. Inclusão reversa porque
    cada r_a+r_b é literalmente um XOR-pair entre cantos diferentes.

  Teorema (Q(n) = 3 para n ≥ 4):
    Pelo Lema 4, Q(n) = dim Span{r_i + r_j : 1 ≤ i < j ≤ 4} em
    GF(2)^E / row(∂).
    Pelo Lema 3, r_1, r_2, r_3, r_4 formam um espaço W ≅ GF(2)^4 dentro
    do quociente. As somas pareadas {r_i+r_j} são vetores em W cuja
    coordenada-soma é zero (cada soma tem peso 2 nas coordenadas i,j,
    soma global = 0 mod 2). Logo Span{r_i+r_j} ⊆ ker(σ) onde σ:W→GF(2)
    é a soma das 4 coordenadas. Mas ker(σ) tem dim 4 − 1 = 3, e os
    vetores r_1+r_2, r_1+r_3, r_1+r_4 são linearmente independentes
    em ker(σ) (são as 3 "diferenças partindo de r_1"). ∎

VERIFICAÇÃO COMPUTACIONAL aqui: re-executa todos os passos para um
conjunto de n e cospe um JSON com o status de cada lema.

O que NÃO está provado neste módulo (lacunas honestas):
  (L5)  rank(Ham(n)) = β₁(n) − Q(n) = β₁(n) − 3.
        Demonstrado computacionalmente para n=6 (exaustivo), n=8,n=10
        (amostragem). Algebricamente: equivale a mostrar que toda
        funcional φ ∈ GF(2)^E que se anula em todos os tours fechados
        está em row(∂) + Span(XOR_pairs). Estabelecer essa identidade
        requer um argumento construtivo de "deformação local de tour"
        que não é coberto aqui.
  (L6)  Tabuleiros n ímpar — os 4 cantos ainda têm grau 2 e Lema 1–4
        valem, mas não existem tours fechados (paridade) → Ham(n) é
        trivial. Refraseamento necessário.
  (L7)  Tabuleiros n×m retangulares.

Saída: data/results/proof_verification.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from verify_small_cases import (  # noqa: E402
    build_graph, label, boundary_matrix, gf2_rank,
    quotient_dim, corner_edges, xor_pair_vectors,
)


def lemma1_corner_degree(n, adj, corners):
    return all(len(adj[c]) == 2 for c in corners)


def lemma2_corner_chain(n, edges, corners, mand_per_corner, V):
    """Verifica: ∀ canto c, ∂(δ_c) = e₁(c) + e₂(c)."""
    E = len(edges)
    for c in corners:
        alpha = np.zeros(V, dtype=np.uint8)
        alpha[c] = 1
        boundary = np.zeros(E, dtype=np.uint8)
        for i, (u, v) in enumerate(edges):
            if alpha[u] ^ alpha[v]:
                boundary[i] = 1
        target = np.zeros(E, dtype=np.uint8)
        for ei in mand_per_corner[c]:
            target[ei] = 1
        if not np.array_equal(boundary, target):
            return False
    return True


def lemma3_reps_independent(n, edges, corners, mand_per_corner, B):
    """rank dos 4 representantes (e₁ de cada canto) mod row(∂)."""
    E = len(edges)
    reps = np.zeros((4, E), dtype=np.uint8)
    for k, c in enumerate(corners):
        reps[k, mand_per_corner[c][0]] = 1
    return quotient_dim(reps, B) == 4


def lemma4_image_is_pairwise_diffs(n, edges, corners, mand_per_corner, mand_idx, B):
    """
    Constrói:
      A = imagem dos 28 XORs em GF(2)^E / row(∂)  (rank quociente)
      C = imagem das 6 diferenças r_i+r_j em GF(2)^E / row(∂)
    e verifica A == C (mesmo span no quociente).

    Implementação: rank quociente de [XOR; diffs] mod row(∂) deve ser igual
    a cada um separadamente.
    """
    E = len(edges)
    XOR, _ = xor_pair_vectors(mand_idx, E)
    reps = np.zeros((4, E), dtype=np.uint8)
    for k, c in enumerate(corners):
        reps[k, mand_per_corner[c][0]] = 1
    diffs = []
    for k1 in range(4):
        for k2 in range(k1 + 1, 4):
            diffs.append((reps[k1] ^ reps[k2]).astype(np.uint8))
    D = np.stack(diffs, axis=0)
    qX = quotient_dim(XOR, B)
    qD = quotient_dim(D, B)
    combined = np.vstack([XOR, D])
    qC = quotient_dim(combined, B)
    return (qX == qD == qC), qX, qD, qC


def theorem_Q_equals_3(n, edges, corners, mand_per_corner, mand_idx, B):
    """Verifica Q(n) = 3."""
    E = len(edges)
    XOR, _ = xor_pair_vectors(mand_idx, E)
    return quotient_dim(XOR, B) == 3


def verify_proof_for_n(n: int) -> dict:
    adj, edges = build_graph(n)
    V = n * n
    B = boundary_matrix(edges, V)
    try:
        corners, mand_per_corner, mand_idx = corner_edges(n, adj, edges)
    except ValueError as e:
        return {"n": n, "error": str(e)}

    L1 = lemma1_corner_degree(n, adj, corners)
    L2 = lemma2_corner_chain(n, edges, corners, mand_per_corner, V)
    L3 = lemma3_reps_independent(n, edges, corners, mand_per_corner, B)
    L4_ok, qX, qD, qC = lemma4_image_is_pairwise_diffs(
        n, edges, corners, mand_per_corner, mand_idx, B)
    Tok = theorem_Q_equals_3(n, edges, corners, mand_per_corner, mand_idx, B)

    return {
        "n": n,
        "V": V,
        "E": len(edges),
        "beta1": len(edges) - V + 1,
        "L1_corners_deg2": L1,
        "L2_corner_chain": L2,
        "L3_reps_independent": L3,
        "L4_span_equals_diffs": L4_ok,
        "L4_dims": {"Q(28 XORs)": qX, "Q(6 diffs)": qD, "Q(combined)": qC},
        "Theorem_Q_equals_3": Tok,
        "all_passed": bool(L1 and L2 and L3 and L4_ok and Tok),
    }


def main():
    out = []
    print(f"{'n':>3} {'L1':>4} {'L2':>4} {'L3':>4} {'L4':>4} {'Q=3':>5} {'all':>5}")
    for n in [4, 5, 6, 7, 8, 9, 10, 11, 12]:
        r = verify_proof_for_n(n)
        out.append(r)
        if "error" in r:
            print(f"{n:>3}  ERRO: {r['error']}")
            continue
        flags = lambda b: "✓" if b else "✗"
        print(f"{n:>3} {flags(r['L1_corners_deg2']):>4} "
              f"{flags(r['L2_corner_chain']):>4} "
              f"{flags(r['L3_reps_independent']):>4} "
              f"{flags(r['L4_span_equals_diffs']):>4} "
              f"{flags(r['Theorem_Q_equals_3']):>5} "
              f"{flags(r['all_passed']):>5}")

    out_path = ROOT / "data" / "results" / "proof_verification.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nSalvo: {out_path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
