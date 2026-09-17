#!/usr/bin/env python3
"""
verify_small_cases.py
=====================
T1 do Teorema do Deficit 3.

Para cada n ∈ {6, 8, 10}:
  1. Constrói G_n (V, E, β₁)
  2. Identifica Mand(n) (8 arestas dos 4 cantos, grau-2)
  3. Amostra K tours hamiltonianos fechados via Z3
       (grau-2 + quebra de simetria aleatória + sub-tour elimination)
  4. Calcula rank_{GF(2)}(T)
  5. Constrói os 28 vetores XOR e calcula
       Q(n) = dim π(Span(XOR)) onde π : GF(2)^E → GF(2)^E / R(∂₁)
  6. Amostra K' 2-fatores (relaxa conectividade) e calcula rank_{GF(2)}(F)
       — deve coincidir com rank(T).

Saída: data/results/verify_n{n}.json e impressão tabulada.

CLI:
  python verify_small_cases.py --n 8 --tours 3000 --twofactors 500
  python verify_small_cases.py --all-default   # roda n=6,8,10 com defaults
"""
from __future__ import annotations

import argparse
import json
import random
import sys
import time
from pathlib import Path
from typing import List, Tuple

import numpy as np
from z3 import Bool, Or, PbEq, Solver, sat, is_true

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
RESULTS = DATA / "results"

MOVES = [(2, 1), (2, -1), (-2, 1), (-2, -1),
         (1, 2), (1, -2), (-1, 2), (-1, -2)]


# ── grafo parametrizado ─────────────────────────────────────────────

def vid(r, c, n):
    return r * n + c


def label(v, n):
    r, c = divmod(v, n)
    return f"{chr(ord('A') + c)}{n - r}"


def build_graph(n: int):
    """Retorna (adj, edges). edges é lista ordenada de (u, v) com u < v."""
    total = n * n
    adj = [[] for _ in range(total)]
    edge_set = set()
    for r in range(n):
        for c in range(n):
            v = vid(r, c, n)
            for dr, dc in MOVES:
                nr, nc = r + dr, c + dc
                if 0 <= nr < n and 0 <= nc < n:
                    u = vid(nr, nc, n)
                    adj[v].append(u)
                    if u > v:
                        edge_set.add((v, u))
    return adj, sorted(edge_set)


def boundary_matrix(edges, total):
    """∂₁ ∈ GF(2)^{V × E} sobre GF(2)."""
    E = len(edges)
    B = np.zeros((total, E), dtype=np.uint8)
    for i, (u, v) in enumerate(edges):
        B[u, i] = 1
        B[v, i] = 1
    return B


# ── álgebra GF(2) ───────────────────────────────────────────────────

def gf2_rref(M):
    """Row-reduced echelon form em GF(2). Retorna (rref, pivots, rank)."""
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
    return A[:r], pivots, r


def gf2_rank(M):
    return gf2_rref(M)[2]


def quotient_dim(V_rows, R_rows):
    """
    Calcula dim(π(Span(V_rows))) onde π : GF(2)^E → GF(2)^E / Span(R_rows).
    Estratégia: rank([V; R]) − rank(R) = dim do span "novo" mod R.
    """
    if V_rows.size == 0:
        return 0
    if R_rows.size == 0:
        return gf2_rank(V_rows)
    M = np.vstack([V_rows, R_rows]).astype(np.uint8)
    return gf2_rank(M) - gf2_rank(R_rows)


# ── identificação dos cantos e Mand(n) ──────────────────────────────

def corner_edges(n: int, adj, edges):
    """
    Retorna:
      corners = [v_NW, v_NE, v_SW, v_SE]    (ids dos cantos)
      mand_per_corner = { corner_v : [idx_e1, idx_e2] }
      mand_idx_sorted = lista ordenada dos 8 índices de arestas obrigatórias
    """
    edge_to_idx = {e: i for i, e in enumerate(edges)}
    corners = [
        vid(0, 0, n),
        vid(0, n - 1, n),
        vid(n - 1, 0, n),
        vid(n - 1, n - 1, n),
    ]
    mand_per_corner = {}
    for c in corners:
        if len(adj[c]) != 2:
            raise ValueError(f"canto {label(c,n)} tem grau {len(adj[c])} != 2 (n={n})")
        es = []
        for u in adj[c]:
            a, b = (c, u) if c < u else (u, c)
            es.append(edge_to_idx[(a, b)])
        mand_per_corner[c] = sorted(es)
    all_mand = sorted({i for lst in mand_per_corner.values() for i in lst})
    if len(all_mand) != 8:
        raise ValueError(f"esperado 8 arestas obrigatórias, obtive {len(all_mand)} (n={n})")
    return corners, mand_per_corner, all_mand


def xor_pair_vectors(mand_idx: List[int], E: int):
    """Os 28 vetores v_{ij} para i<j em Mand."""
    rows = []
    pairs = []
    for k in range(len(mand_idx)):
        for l in range(k + 1, len(mand_idx)):
            i, j = mand_idx[k], mand_idx[l]
            v = np.zeros(E, dtype=np.uint8)
            v[i] = 1
            v[j] = 1
            rows.append(v)
            pairs.append((i, j))
    return np.stack(rows, axis=0), pairs


# ── conectividade ───────────────────────────────────────────────────

def is_single_tour(active_edges, total):
    adj = [[] for _ in range(total)]
    for u, v in active_edges:
        adj[u].append(v)
        adj[v].append(u)
    if any(len(a) != 2 for a in adj):
        return False, []
    visited = [False] * total
    visited[0] = True
    tour = [0]
    cur, prev = 0, -1
    for _ in range(total - 1):
        nxt = adj[cur][0] if adj[cur][0] != prev else adj[cur][1]
        if visited[nxt]:
            return False, []
        visited[nxt] = True
        tour.append(nxt)
        prev, cur = cur, nxt
    return (adj[cur][0] == 0 or adj[cur][1] == 0), tour


# ── amostragem ──────────────────────────────────────────────────────

def build_base_solver(edges, total):
    E = len(edges)
    xvars = [Bool(f"x_{i}") for i in range(E)]
    inc = {v: [] for v in range(total)}
    for i, (u, v) in enumerate(edges):
        inc[u].append(i)
        inc[v].append(i)
    s = Solver()
    for v, lst in inc.items():
        s.add(PbEq([(xvars[i], 1) for i in lst], 2))
    return s, xvars, inc


def sample_solutions(n: int, K: int, *, connected: bool, seed: int,
                     timeout_per_sample_s: int = 30,
                     max_sym_attempts: int = 4,
                     verbose: bool = True,
                     log_every: int = 25):
    """
    Amostra K assinaturas distintas. Se connected=True, exige tour único
    (Hamiltoniano fechado). Caso contrário, aceita qualquer 2-fator.
    """
    total = n * n
    adj, edges = build_graph(n)
    E = len(edges)
    rng = random.Random(seed)

    solver, xvars, _ = build_base_solver(edges, total)

    sigs: List[np.ndarray] = []
    n_attempts = 0
    n_rej_subtours = 0
    n_unsat = 0
    t0 = time.perf_counter()

    while len(sigs) < K:
        n_attempts += 1
        solver.push()

        local_unsat = True
        for _ in range(max_sym_attempts):
            edge_pick = rng.randrange(E)
            val_pick = rng.randrange(2)
            solver.push()
            solver.add(xvars[edge_pick] == (val_pick == 1))
            solver.set("timeout", int(timeout_per_sample_s * 1000))
            res = solver.check()
            if res == sat:
                local_unsat = False
                break
            solver.pop()

        if local_unsat:
            n_unsat += 1
            solver.pop()
            continue

        m = solver.model()
        sig = np.zeros(E, dtype=np.uint8)
        active = []
        for i, (u, v) in enumerate(edges):
            if is_true(m.evaluate(xvars[i])):
                sig[i] = 1
                active.append((u, v))

        if connected:
            ok, _ = is_single_tour(active, total)
        else:
            ok = True

        solver.pop()  # tira sym
        solver.pop()  # tira outer

        # corte permanente
        solver.add(Or([xvars[i] != bool(sig[i]) for i in range(E)]))

        if ok:
            sigs.append(sig)
            if verbose and len(sigs) % log_every == 0:
                elapsed = time.perf_counter() - t0
                rate = len(sigs) / max(elapsed, 1e-6)
                print(f"    [{len(sigs)}/{K}] attempts={n_attempts} "
                      f"rej_subtour={n_rej_subtours} unsat={n_unsat} "
                      f"elapsed={elapsed:.1f}s rate={rate:.2f}/s",
                      flush=True)
        else:
            n_rej_subtours += 1

    elapsed = time.perf_counter() - t0
    stats = {
        "K_requested": K,
        "K_obtained": len(sigs),
        "n_attempts": n_attempts,
        "n_rejected_subtours": n_rej_subtours,
        "n_unsat_local": n_unsat,
        "elapsed_s": round(elapsed, 3),
        "connected": connected,
    }
    return np.stack(sigs, axis=0), stats


# ── verificação completa para um n ─────────────────────────────────

def verify(n: int, *, K_tours: int, K_2fact: int, seed: int,
           timeout_per_sample_s: int) -> dict:
    print(f"\n{'='*70}\nVERIFICAÇÃO n = {n}\n{'='*70}")
    total = n * n
    adj, edges = build_graph(n)
    E = len(edges)
    V = total
    beta1 = E - V + 1
    B = boundary_matrix(edges, V)
    rank_B = gf2_rank(B)
    beta1_via_rank = E - rank_B

    print(f"\nV = {V}    E = {E}    β₁ = E − V + 1 = {beta1}")
    print(f"rank(∂₁) = {rank_B}    β₁ via posto = {beta1_via_rank}")
    assert beta1 == beta1_via_rank, "inconsistência β₁ vs rank(∂₁)"

    corners, mand_per_corner, mand_idx = corner_edges(n, adj, edges)
    print(f"\nCantos ({len(corners)}): {[label(c, n) for c in corners]}")
    for c, lst in mand_per_corner.items():
        labs = []
        for ei in lst:
            u, v = edges[ei]
            labs.append(f"{label(u,n)}-{label(v,n)} (idx {ei})")
        print(f"  {label(c,n)}: {labs}")
    print(f"Mand(n) índices (ordenados): {mand_idx}")

    # ── Q(n) via 28 XORs vs row(∂₁) ──
    XOR, pairs = xor_pair_vectors(mand_idx, E)
    rank_xor = gf2_rank(XOR)
    Q_n = quotient_dim(XOR, B)
    print(f"\nXOR_pairs(n): 28 vetores construídos")
    print(f"  rank(Span(28 XOR))            = {rank_xor}")
    print(f"  rank(B)                       = {rank_B}")
    print(f"  rank([XOR; B])                = {gf2_rank(np.vstack([XOR, B]))}")
    print(f"  Q(n) = rank([XOR;B]) − rank(B) = {Q_n}")
    print(f"  Predição (conjectura): 3")

    # ── amostragem de tours hamiltonianos fechados ──
    print(f"\nAmostrando {K_tours} tours (connected=True) ...")
    T, stats_T = sample_solutions(
        n, K_tours, connected=True, seed=seed,
        timeout_per_sample_s=timeout_per_sample_s,
        verbose=True)
    rank_T = gf2_rank(T)
    deficit_T = beta1 - rank_T
    print(f"\n  K = {T.shape[0]} tours obtidos em {stats_T['elapsed_s']}s")
    print(f"  rank_{{GF(2)}}(T) = {rank_T}")
    print(f"  deficit(n) = β₁ − rank(T) = {beta1} − {rank_T} = {deficit_T}")
    print(f"  Predição: deficit = 3 → rank esperado = {beta1 - 3}")

    # ── amostragem de 2-fatores (sem conectividade) ──
    rank_F = None
    deficit_F = None
    stats_F = None
    F_shape = None
    if K_2fact > 0:
        print(f"\nAmostrando {K_2fact} 2-fatores (connected=False) ...")
        F, stats_F = sample_solutions(
            n, K_2fact, connected=False, seed=seed + 7919,
            timeout_per_sample_s=timeout_per_sample_s,
            verbose=True)
        rank_F = gf2_rank(F)
        deficit_F = beta1 - rank_F
        F_shape = list(F.shape)
        print(f"  K' = {F.shape[0]} 2-fatores em {stats_F['elapsed_s']}s")
        print(f"  rank_{{GF(2)}}(F) = {rank_F}")
        print(f"  deficit(F) = β₁ − rank(F) = {deficit_F}")
        print(f"  rank(F) == rank(T)? {rank_F == rank_T}")

    # ── projeção: rank das 8 obrigatórias mod B (sanity) ──
    mand_vec = np.zeros((8, E), dtype=np.uint8)
    for k, ei in enumerate(mand_idx):
        mand_vec[k, ei] = 1
    rank_mand_raw = gf2_rank(mand_vec)
    rank_mand_quot = quotient_dim(mand_vec, B)
    print(f"\nSanity: 8 vetores indicadores de obrigatórias")
    print(f"  rank raw     = {rank_mand_raw}  (esperado 8)")
    print(f"  rank mod ∂₁  = {rank_mand_quot}  (predição 4 = #cantos)")

    out = {
        "n": n,
        "V": V,
        "E": E,
        "beta1": beta1,
        "rank_boundary_gf2": int(rank_B),
        "corners_labels": [label(c, n) for c in corners],
        "mand_indices": mand_idx,
        "mand_per_corner": {label(c, n): lst for c, lst in mand_per_corner.items()},
        "rank_xor_raw": int(rank_xor),
        "Q_n": int(Q_n),
        "rank_mand_raw": int(rank_mand_raw),
        "rank_mand_mod_boundary": int(rank_mand_quot),
        "tours": {
            "K_obtained": int(T.shape[0]),
            "rank_T": int(rank_T),
            "deficit": int(deficit_T),
            "stats": stats_T,
        },
        "two_factors": {
            "K_obtained": int(F_shape[0]) if F_shape else 0,
            "rank_F": int(rank_F) if rank_F is not None else None,
            "deficit": int(deficit_F) if deficit_F is not None else None,
            "matches_tour_rank": (rank_F == rank_T) if rank_F is not None else None,
            "stats": stats_F,
        },
        "predictions": {
            "Q_n_predicted": 3,
            "deficit_predicted": 3,
        },
        "verdict": {
            "Q_matches": int(Q_n) == 3,
            "deficit_matches": int(deficit_T) == 3,
        },
    }
    return out


# ── runner ──────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--n", type=int, default=None,
                   help="tamanho do tabuleiro (ímpar permitido — emite aviso)")
    p.add_argument("--tours", type=int, default=3000,
                   help="K = quantos tours hamiltonianos amostrar")
    p.add_argument("--twofactors", type=int, default=500,
                   help="K' = quantos 2-fatores amostrar (0 para pular)")
    p.add_argument("--seed", type=int, default=2026)
    p.add_argument("--timeout", type=int, default=60,
                   help="timeout por sample em segundos")
    p.add_argument("--all-default", action="store_true",
                   help="roda n=6,8,10 com defaults sensatos")
    args = p.parse_args()

    RESULTS.mkdir(parents=True, exist_ok=True)

    if args.all_default:
        plan = [
            dict(n=6, tours=2000, twofactors=500),
            dict(n=8, tours=3000, twofactors=500),
            dict(n=10, tours=1500, twofactors=300),
        ]
    elif args.n is not None:
        plan = [dict(n=args.n, tours=args.tours, twofactors=args.twofactors)]
    else:
        p.error("forneça --n ou --all-default")

    all_results = []
    for cfg in plan:
        n = cfg["n"]
        if n % 2 != 0:
            print(f"AVISO: n={n} é ímpar. Conjectura só foi formulada para n par.")
        res = verify(n,
                     K_tours=cfg["tours"],
                     K_2fact=cfg["twofactors"],
                     seed=args.seed,
                     timeout_per_sample_s=args.timeout)
        all_results.append(res)

        out_path = RESULTS / f"verify_n{n}.json"
        with open(out_path, "w") as f:
            json.dump(res, f, indent=2)
        print(f"\nSalvo: {out_path.relative_to(ROOT)}")

    # tabela consolidada
    print(f"\n{'='*70}")
    print(f"{'CONSOLIDADO':^70}")
    print(f"{'='*70}")
    print(f"{'n':>4} {'V':>5} {'E':>5} {'β₁':>5} {'rank(T)':>8} {'def':>4} {'Q':>3}  {'OK':>4}")
    for r in all_results:
        ok = "✓" if r["verdict"]["Q_matches"] and r["verdict"]["deficit_matches"] else "✗"
        print(f"{r['n']:>4} {r['V']:>5} {r['E']:>5} {r['beta1']:>5} "
              f"{r['tours']['rank_T']:>8} {r['tours']['deficit']:>4} {r['Q_n']:>3}  {ok:>4}")

    consolidated_path = RESULTS / "consolidated.json"
    with open(consolidated_path, "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"\nSalvo: {consolidated_path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
