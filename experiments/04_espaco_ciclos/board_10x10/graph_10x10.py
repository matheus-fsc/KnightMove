#!/usr/bin/env python3
"""
graph_10x10.py
==============
Construção do grafo do cavalo 10×10 e matriz de bordas ∂₁ sobre GF(2).

Saídas:
  data/boundary_matrix_10x10.npy   — ∂₁ shape (100, E), dtype uint8
  data/graph_meta_10x10.json       — metadados (V, E, β₁, graus)
"""

import json
from pathlib import Path

import numpy as np

BOARD = 10
TOTAL = BOARD * BOARD
MOVES = [(2, 1), (2, -1), (-2, 1), (-2, -1),
         (1, 2), (1, -2), (-1, 2), (-1, -2)]

DATA_DIR = Path(__file__).resolve().parent / "data"


def vid(r, c):
    return r * BOARD + c


def vrc(v):
    return divmod(v, BOARD)


def label(v):
    r, c = vrc(v)
    return chr(ord('A') + c) + str(BOARD - r)


def build_graph():
    """Retorna (ADJ, EDGES). EDGES é lista ordenada lexicograficamente de tuplas (u, v) com u<v."""
    adj = [[] for _ in range(TOTAL)]
    edge_set = set()
    for r in range(BOARD):
        for c in range(BOARD):
            v = vid(r, c)
            for dr, dc in MOVES:
                nr, nc = r + dr, c + dc
                if 0 <= nr < BOARD and 0 <= nc < BOARD:
                    u = vid(nr, nc)
                    adj[v].append(u)
                    if u > v:
                        edge_set.add((v, u))
    edges = sorted(edge_set)
    return adj, edges


def boundary_matrix(edges, n_vertices=TOTAL):
    """∂₁ sobre GF(2): coluna e tem 1s em u e v se e={u,v}."""
    E = len(edges)
    B = np.zeros((n_vertices, E), dtype=np.uint8)
    for i, (u, v) in enumerate(edges):
        B[u, i] = 1
        B[v, i] = 1
    return B


def gf2_rank(M):
    """Rank de M sobre GF(2) via eliminação gaussiana."""
    A = M.copy().astype(np.uint8)
    rows, cols = A.shape
    rank = 0
    r = 0
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
        r += 1
        rank += 1
    return rank


def is_connected(adj):
    """BFS a partir do vértice 0."""
    visited = {0}
    stack = [0]
    while stack:
        v = stack.pop()
        for u in adj[v]:
            if u not in visited:
                visited.add(u)
                stack.append(u)
    return len(visited) == TOTAL


def vertex_class(r, c):
    """Classifica vértice por posição: corner / edge / near-corner / interior."""
    on_h_edge = r in (0, BOARD - 1)
    on_v_edge = c in (0, BOARD - 1)
    if on_h_edge and on_v_edge:
        return "corner"
    if on_h_edge or on_v_edge:
        return "edge"
    if r in (1, BOARD - 2) or c in (1, BOARD - 2):
        return "near_edge"
    return "interior"


def main():
    print("=" * 65)
    print("GRAFO DO CAVALO 10×10")
    print("=" * 65)

    adj, edges = build_graph()
    V = TOTAL
    E = len(edges)

    degrees = [len(a) for a in adj]
    deg_min = min(degrees)
    deg_max = max(degrees)
    deg_mean = sum(degrees) / V

    deg_counter = {d: degrees.count(d) for d in range(deg_min, deg_max + 1)}

    # classes posicionais
    class_counts = {}
    class_degrees = {}
    for v in range(V):
        r, c = vrc(v)
        cls = vertex_class(r, c)
        class_counts[cls] = class_counts.get(cls, 0) + 1
        class_degrees.setdefault(cls, []).append(degrees[v])

    print(f"\n|V|                 : {V}")
    print(f"|E|                 : {E}")
    print(f"β₁ esperado (E-V+1) : {E - V + 1}")

    print(f"\nGraus: min={deg_min} max={deg_max} mean={deg_mean:.3f}")
    print("Distribuição de graus:")
    for d, n in sorted(deg_counter.items()):
        print(f"  grau {d:2d}: {n:>3d} vértices")

    print("\nVértices por classe posicional:")
    for cls, n in class_counts.items():
        degs = class_degrees[cls]
        print(f"  {cls:12s}: n={n:>3d}  graus={sorted(set(degs))}")

    connected = is_connected(adj)
    print(f"\nConectividade       : {'CONEXO' if connected else 'NÃO CONEXO'}")

    # ── matriz de bordas ∂₁ sobre GF(2) ──────────────────────────────
    print("\nConstruindo matriz ∂₁ ...")
    B = boundary_matrix(edges, V)
    print(f"  shape           : {B.shape}")
    print(f"  dtype           : {B.dtype}")
    print(f"  células ativas  : {int(B.sum())} (esperado {2*E})")

    # rank GF(2) — para confirmar β₁ via teorema do posto
    rank_B = gf2_rank(B)
    beta1_via_rank = E - rank_B
    print(f"  rank GF(2)(∂₁)  : {rank_B}  (esperado V-1 = {V-1})")
    print(f"  β₁ via posto    : E - rank = {beta1_via_rank}")
    if beta1_via_rank != E - V + 1:
        print(f"  AVISO: β₁ via posto difere do esperado")

    # ── salva artefatos ──────────────────────────────────────────────
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    out_npy = DATA_DIR / "boundary_matrix_10x10.npy"
    np.save(out_npy, B)
    print(f"\nSalvo: {out_npy.relative_to(Path(__file__).resolve().parent)}")

    edges_path = DATA_DIR / "edges_10x10.json"
    with open(edges_path, "w") as f:
        json.dump({
            "edges_uv": [[int(u), int(v)] for u, v in edges],
            "edge_labels": [f"{label(u)}-{label(v)}" for u, v in edges],
        }, f)
    print(f"Salvo: {edges_path.relative_to(Path(__file__).resolve().parent)}")

    meta_path = DATA_DIR / "graph_meta_10x10.json"
    meta = {
        "board": BOARD,
        "V": V,
        "E": E,
        "beta_1": E - V + 1,
        "rank_boundary_gf2": int(rank_B),
        "connected": bool(connected),
        "degree_distribution": {str(d): n for d, n in sorted(deg_counter.items())},
        "degree_min": deg_min,
        "degree_max": deg_max,
        "degree_mean": round(deg_mean, 4),
        "class_counts": class_counts,
        "class_degree_ranges": {
            cls: {"min": min(d), "max": max(d)} for cls, d in class_degrees.items()
        },
    }
    with open(meta_path, "w") as f:
        json.dump(meta, f, indent=2)
    print(f"Salvo: {meta_path.relative_to(Path(__file__).resolve().parent)}")

    print()
    print("─" * 65)
    print(f"RESUMO  →  V={V}  E={E}  β₁={E-V+1}")
    print("─" * 65)


if __name__ == "__main__":
    main()
