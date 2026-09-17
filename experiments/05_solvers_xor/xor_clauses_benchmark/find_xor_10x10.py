#!/usr/bin/env python3
"""
find_xor_10x10.py
=================
T2: identifica as cláusulas XOR de paridade do 10×10.

Mesma análise que `forbidden_cycles_6x6/obstruction_theorem.py` faz no 6×6:

  1. Construir os C(8,2) = 28 vetores XOR entre as 8 arestas obrigatórias
     dos cantos do 10×10.
  2. Calcular dim(Span(28 XORs) mod row_space(∂)) sobre GF(2).
     - row_space(∂) é o anulador de H₁ — funcionais que zeram em todo ciclo.
     - dim do quociente == deficit detectado por XORs entre obrigatórias.
  3. Validar contra as 5000 amostras T existentes: T @ v_{ij}^T = 0?
     Conta quantos dos 28 pares são empiricamente válidos.
  4. Selecionar 3 pares linearmente independentes mod row_space(∂)
     como base mínima.
  5. (Opcional) verificação Z3: forçar x_a=1, x_b=0 para cada par escolhido
     e testar UNSAT.

Salva: data/xor_10x10_clauses.json   (3 pares escolhidos + diagnostic)
"""

from __future__ import annotations

import argparse
import itertools
import json
import sys
import time
from pathlib import Path

import numpy as np
from z3 import Bool, Or, PbEq, Solver, sat, unsat, is_true

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
REPO = ROOT.parent

sys.path.insert(0, str(REPO / "board_10x10"))
from graph_10x10 import build_graph, TOTAL, label  # noqa: E402


# ── GF(2) linear algebra ──────────────────────────────────────────────

def gf2_rref(M):
    """Row-reduced echelon form em GF(2). Retorna (R, pivots, rank)."""
    A = M.copy().astype(np.uint8)
    rows, cols = A.shape
    pivots = []
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
        pivots.append(c)
        r += 1
    return A, pivots, r


def gf2_rank(M):
    return gf2_rref(M)[2]


# ── carrega 10×10 ─────────────────────────────────────────────────────

def load_10x10_data():
    adj, edges = build_graph()
    E = len(edges)

    # ∂₁: linhas = vértices, colunas = arestas
    inc = np.zeros((TOTAL, E), dtype=np.uint8)
    for i, (u, v) in enumerate(edges):
        inc[u, i] = 1
        inc[v, i] = 1

    # tours amostrados
    samples_dir = REPO / "board_10x10" / "data" / "samples"
    batches = sorted(samples_dir.glob("tours_10x10_batch_*.npy"))
    if not batches:
        raise FileNotFoundError("Nenhum batch de tours encontrado.")
    Ts = [np.load(b) for b in batches]
    T = np.concatenate(Ts, axis=0).astype(np.uint8)

    # mandatórias
    mand_p = REPO / "board_10x10" / "data" / "invariants" / "mandatory_edges.json"
    with open(mand_p) as f:
        md = json.load(f)
    mandatory_idx = sorted([m["idx"] for m in md["mandatory"]])
    mandatory_labels = {m["idx"]: m["label"] for m in md["mandatory"]}

    return adj, edges, inc, T, mandatory_idx, mandatory_labels


# ── T2.2: dim quociente ────────────────────────────────────────────────

def project_quotient(vectors_E, row_space_E):
    """
    Retorna a dimensão de Span(vectors) mod Span(row_space) e a
    sub-base de `vectors` que forma uma base independente do quociente.

    Implementação: empilha `row_space` (rank R) acima dos `vectors` e
    faz RREF. As linhas que viram pivôs nas posições dos `vectors`
    (índices >= R) compõem a base do quociente; o número delas é
    dim(quociente).
    """
    R = row_space_E.shape[0]
    M = np.vstack([row_space_E, vectors_E]).astype(np.uint8)
    A = M.copy()
    rows, cols = A.shape
    pivot_rows = []  # indices em A
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
        pivot_rows.append((c, r))
        r += 1
    total_rank = r
    quotient_dim = total_rank - gf2_rank(row_space_E)

    # Selecionar quais dos `vectors` (linhas com índice >= R no M original)
    # são linearmente independentes mod row_space.
    # Repete a redução adicionando vetores um a um:
    chosen = []
    used_basis = row_space_E.copy()
    for i, v in enumerate(vectors_E):
        # se v não está em Span(used_basis), adiciona
        cand = np.vstack([used_basis, v[None, :]])
        if gf2_rank(cand) > gf2_rank(used_basis):
            chosen.append(i)
            used_basis = cand
        if len(chosen) >= quotient_dim:
            break
    return quotient_dim, chosen


def main_find(args):
    print("─" * 70)
    print("T2: identificar cláusulas XOR do 10×10")
    print("─" * 70)

    adj, edges, inc, T, mandatory_idx, mandatory_labels = load_10x10_data()
    V = TOTAL
    E = len(edges)
    print(f"  V={V}  E={E}  β₁={E - V + 1}")
    print(f"  Tours amostrados (T): {T.shape}")
    print(f"  Mandatórias: {len(mandatory_idx)} arestas")
    for m in mandatory_idx:
        print(f"    idx {m:3d}  {mandatory_labels[m]}")

    # rank de T amostrado
    rank_T = gf2_rank(T)
    beta1 = E - V + 1
    print(f"  rank(T amostrado) = {rank_T}  β₁ = {beta1}  "
          f"deficit empírico = {beta1 - rank_T}")

    # row_space(∂) sobre GF(2): rank = V-1 (∂ tem nullity 1 — V conexo)
    rs, _, rank_inc = gf2_rref(inc)
    rs_reduced = rs[:rank_inc]   # shape (V-1, E)
    print(f"  rank(∂) = {rank_inc}  (esperado V-1 = {V-1})")

    # ── construir os 28 vetores XOR pares ────────────────────────────
    pairs = list(itertools.combinations(mandatory_idx, 2))
    XOR = np.zeros((len(pairs), E), dtype=np.uint8)
    for k, (i, j) in enumerate(pairs):
        XOR[k, i] = 1
        XOR[k, j] = 1
    print(f"  C({len(mandatory_idx)},2) = {len(pairs)} vetores XOR pares")

    # ── validação empírica: cada XOR zera em todos os tours? ─────────
    prod = (T @ XOR.T) % 2  # shape (n_samples, 28)
    pair_valid = (prod.sum(axis=0) == 0)
    n_valid = int(pair_valid.sum())
    print(f"  XOR pares válidos em todas as {T.shape[0]} amostras: "
          f"{n_valid}/{len(pairs)}")

    # ── dim quociente XOR mod row_space(∂) ───────────────────────────
    qd, chosen_idx = project_quotient(XOR, rs_reduced)
    print(f"  dim(Span(XOR_pares) mod row_space(∂)) = {qd}")

    # ── base mínima: 3 pares linearmente independentes ───────────────
    if qd == 0:
        print("  AVISO: dim quociente = 0; XOR entre obrigatórias é redundante.")
    chosen_pairs = [pairs[i] for i in chosen_idx]
    chosen_labels = [(mandatory_labels[a], mandatory_labels[b])
                     for (a, b) in chosen_pairs]
    print(f"  Base mínima: {len(chosen_idx)} pares")
    for (a, b), (la, lb) in zip(chosen_pairs, chosen_labels):
        print(f"    XOR: x_{{{la}}} ⊕ x_{{{lb}}} = 0   (idx {a}, {b})")

    # ── extras: dimensões para o quadro do README ────────────────────
    span_xor_dim = gf2_rank(XOR)
    rs_plus_xor = np.vstack([rs_reduced, XOR])
    rs_plus_xor_dim = gf2_rank(rs_plus_xor)

    print(f"\n  Quadro de dimensões (GF(2)^{E}):")
    print(f"    Span(XOR pares)                 dim {span_xor_dim}")
    print(f"    row_space(∂)                    dim {rank_inc}")
    print(f"    Span(XOR pares) + row_space(∂)  dim {rs_plus_xor_dim}")
    print(f"    quociente (deficit detectado)   dim {qd}")

    # ── verificação Z3: forçar x_a=1, x_b=0 → esperar UNSAT ──────────
    z3_results = []
    if args.verify_z3 and chosen_pairs:
        print("\n  Verificação Z3 — UNSAT esperada ao violar cada XOR:")
        for (a, b) in chosen_pairs:
            la, lb = mandatory_labels[a], mandatory_labels[b]
            t0 = time.perf_counter()
            res = verify_xor_pair(edges, a, b, timeout_ms=args.z3_timeout_ms,
                                  max_iters=args.z3_max_iters)
            dt = time.perf_counter() - t0
            print(f"    XOR({la}, {lb}): {res['status']}  "
                  f"iters={res['iters']}  t={dt:.1f}s")
            z3_results.append({
                "pair": [int(a), int(b)],
                "labels": [la, lb],
                "status": res["status"],
                "iters": res["iters"],
                "elapsed_s": round(dt, 3),
            })

    out = {
        "board": 10,
        "V": V,
        "E": E,
        "beta_1": beta1,
        "rank_T_sampled": int(rank_T),
        "rank_boundary": int(rank_inc),
        "n_samples": int(T.shape[0]),
        "n_mandatory": len(mandatory_idx),
        "mandatory_idx": list(map(int, mandatory_idx)),
        "mandatory_labels": {str(k): v for k, v in mandatory_labels.items()},
        "n_pair_xors": len(pairs),
        "n_pair_xors_empirically_valid": n_valid,
        "span_xor_dim": int(span_xor_dim),
        "rs_plus_xor_dim": int(rs_plus_xor_dim),
        "quotient_dim": int(qd),
        "chosen_pairs": [list(map(int, p)) for p in chosen_pairs],
        "chosen_labels": [list(l) for l in chosen_labels],
        "z3_verification": z3_results,
    }
    DATA.mkdir(parents=True, exist_ok=True)
    outp = DATA / "xor_10x10_clauses.json"
    with open(outp, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nSalvo: {outp.relative_to(REPO)}")
    return out


def verify_xor_pair(edges, a, b, *, timeout_ms=30_000, max_iters=80):
    """
    Tenta provar UNSAT para o tour 10×10 com x_a=1, x_b=0
    (viola o XOR x_a⊕x_b=0).
    """
    E = len(edges)
    xvars = [Bool(f"x_{i}") for i in range(E)]
    s = Solver()
    s.set("timeout", timeout_ms)
    inc = {v: [] for v in range(TOTAL)}
    for i, (u, v) in enumerate(edges):
        inc[u].append(i)
        inc[v].append(i)
    for v, lst in inc.items():
        s.add(PbEq([(xvars[i], 1) for i in lst], 2))
    s.add(xvars[a])
    s.add(Not(xvars[b]))

    for it in range(max_iters):
        res = s.check()
        if res == unsat:
            return {"status": "unsat", "iters": it + 1}
        if res != sat:
            return {"status": "inconclusive", "iters": it + 1}
        m = s.model()
        sig = np.zeros(E, dtype=np.uint8)
        active = []
        for i, (u, v) in enumerate(edges):
            if is_true(m.evaluate(xvars[i])):
                sig[i] = 1
                active.append((u, v))
        ok, _ = is_single_tour_local(active)
        if ok:
            return {"status": "SAT_counterexample", "iters": it + 1, "sig": sig}
        s.add(Or([xvars[i] != bool(sig[i]) for i in range(E)]))
    return {"status": "inconclusive_max_iters", "iters": max_iters}


def is_single_tour_local(active_edges, n_vertices=TOTAL):
    adj = [[] for _ in range(n_vertices)]
    for u, v in active_edges:
        adj[u].append(v)
        adj[v].append(u)
    if any(len(a) != 2 for a in adj):
        return False, []
    visited = [False] * n_vertices
    visited[0] = True
    cur, prev = 0, -1
    cnt = 1
    for _ in range(n_vertices - 1):
        nxt = adj[cur][0] if adj[cur][0] != prev else adj[cur][1]
        if visited[nxt]:
            return False, []
        visited[nxt] = True
        cnt += 1
        prev, cur = cur, nxt
    return (adj[cur][0] == 0 or adj[cur][1] == 0), []


# z3 helpers (import after definition for clarity)
from z3 import Not  # noqa: E402


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--verify-z3", action="store_true",
                   help="Verificar UNSAT no Z3 ao violar cada XOR escolhido")
    p.add_argument("--z3-timeout-ms", type=int, default=60_000)
    p.add_argument("--z3-max-iters", type=int, default=80)
    args = p.parse_args()
    main_find(args)


if __name__ == "__main__":
    main()
