#!/usr/bin/env python3
"""tightness_intermediate.py — Verifica rank(Ham(G(H))) = β₁(G(H)) − Q(H)
em pontos intermediários da família híbrida toro→plano no 6×6.

Configs: 10 pontos com k(H) ∈ {0,1,2,3,4} mais cantos parciais.
"""

# --- bootstrap: torna core/ importavel a partir de qualquer subpasta ---
import sys as _sys, pathlib as _pathlib
for _anc in _pathlib.Path(__file__).resolve().parents:
    if (_anc / "core").is_dir() and (_anc / "pyproject.toml").exists():
        _sys.path.insert(0, str(_anc / "core"))
        break
# --- fim do bootstrap ---

import sys
import time
import numpy as np

sys.path.insert(0, "/home/math/Dev/knight_tour")
from knight_tours_torus import (
    build_graph, State,
    FREE, ACTIVE, INACTIVE, OK, CONTRADICTION, SUBTOUR, COMPLETE_TOUR,
    _process_queue, _propagate_initial, _backtrack,
    _is_complete_tour, _extract_tour,
)
from Q_locality_theorem import (
    build_plane_edges, build_torus_edges, classify_edges, compute_Q,
)
from rank_ham_torus import IncrementalRank


# ── motor de tours em grafo híbrido ───────────────────────────────────────────

def knight_tours_hybrid(n, K, forbidden_edges, seed=None):
    """Tours fechados no toro n×n com forbidden_edges removidas."""
    ctx = build_graph(n, n)
    V, E = ctx["V"], ctx["E"]
    edge_endpoints = ctx["edge_endpoints"]

    edge_to_idx = {(int(edge_endpoints[ei, 0]), int(edge_endpoints[ei, 1])): ei
                   for ei in range(E)}

    forbidden_ei = []
    for u, v in forbidden_edges:
        key = (min(u, v), max(u, v))
        if key in edge_to_idx:
            forbidden_ei.append(edge_to_idx[key])

    rng = np.random.default_rng(seed)
    state = State(V, E)

    # Marca arestas proibidas como INACTIVE via propagação
    initial_queue = [(ei, 0) for ei in forbidden_ei]
    status = _process_queue(state, ctx, initial_queue)
    if status in (CONTRADICTION, SUBTOUR):
        return []
    if status == COMPLETE_TOUR:
        if _is_complete_tour(state, ctx):
            return [_extract_tour(state, ctx)]
        return []

    # Propagação de cantos / forçados
    status = _propagate_initial(state, ctx)
    if status in (CONTRADICTION, SUBTOUR):
        return []

    tours = []
    if status == COMPLETE_TOUR:
        if _is_complete_tour(state, ctx):
            tours.append(_extract_tour(state, ctx))
    else:
        _backtrack(state, ctx, rng, tours, K)

    return tours[:K]


# ── indexação unificada de arestas (sempre sobre E_torus) ─────────────────────

def edge_indexing(n):
    T = build_torus_edges(n)
    edge_list = sorted(T)
    edge_to_idx = {e: i for i, e in enumerate(edge_list)}
    return edge_list, edge_to_idx, len(edge_list)


def tour_signature(tour, n, edge_to_idx, E_torus):
    V = n * n
    sig = np.zeros(E_torus, dtype=np.uint8)
    for i in range(V):
        u, v = int(tour[i]), int(tour[(i + 1) % V])
        sig[edge_to_idx[(min(u, v), max(u, v))]] = 1
    return sig


def collect_tours_hybrid(n, K_per_seed, n_seeds, forbidden_set, log=False):
    tours_all = []
    seen = set()
    t0 = time.time()
    for seed in range(n_seeds):
        batch = knight_tours_hybrid(n, K_per_seed, forbidden_set, seed=seed)
        n_new = 0
        for t in batch:
            key = tuple(int(x) for x in t)
            if key not in seen:
                seen.add(key)
                tours_all.append(t)
                n_new += 1
        if log and (seed % 5 == 0 or seed == n_seeds - 1):
            print(f"      seed={seed:>3}: +{n_new:>4} novos, "
                  f"total={len(tours_all):>5}, t={time.time()-t0:.1f}s")
    return tours_all


# ── teste de um config ────────────────────────────────────────────────────────

def test_config(name, forbidden_set, P, W, n, K_per_seed, n_seeds,
                edge_list, edge_to_idx, E_torus, log=False):
    print(f"\n{'─'*72}")
    print(f"Config {name}")
    print(f"{'─'*72}")
    t_start = time.time()

    V = n * n
    # H = wraps que PERMANECEM no grafo híbrido
    H = W - forbidden_set
    edge_set_hybrid = P | H
    E_hybrid = len(edge_set_hybrid)
    beta1 = E_hybrid - V + 1

    Q, deg2 = compute_Q(n, edge_set_hybrid)
    corners_n = [0, n-1, n*(n-1), n*n-1]
    k_H = sum(1 for c in corners_n if c in deg2)

    print(f"  |forbidden|={len(forbidden_set)}, E_hybrid={E_hybrid}, "
          f"β₁={beta1}, k(H)={k_H}, Q(H)={Q}")

    # Coleta tours
    tours = collect_tours_hybrid(n, K_per_seed, n_seeds, forbidden_set, log=log)
    print(f"  tours únicos: {len(tours)}")

    if len(tours) == 0:
        print(f"  ⚠ Sem tours.")
        return {"name": name, "k_H": k_H, "Q": Q, "beta1": beta1,
                "rank": None, "deficit": None, "match": None,
                "n_tours": 0, "time": time.time()-t_start}

    # Rank incremental
    rank_eng = IncrementalRank(E_torus)
    rank_log = []
    log_every = max(20, len(tours) // 40)
    for i, t in enumerate(tours):
        sig = tour_signature(t, n, edge_to_idx, E_torus)
        rank_eng.add(sig)
        if (i + 1) % log_every == 0:
            rank_log.append((i + 1, rank_eng.rank))

    rank = rank_eng.rank
    deficit = beta1 - rank
    match = (deficit == Q)

    # Convergência: últimos 5 pontos iguais?
    last_vals = [r for _, r in rank_log[-5:]]
    converged = len(set(last_vals)) == 1 if len(last_vals) >= 2 else False

    print(f"  rank={rank}, deficit={deficit}, Q={Q}, match={match}")
    print(f"  últimos 5: {last_vals}  convergiu={converged}")
    print(f"  tempo: {time.time()-t_start:.1f}s")

    return {
        "name": name, "k_H": k_H, "Q": Q, "beta1": beta1,
        "rank": rank, "deficit": deficit, "match": match,
        "n_tours": len(tours), "converged": converged,
        "last_vals": last_vals, "time": time.time()-t_start,
    }


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    n = 6
    print(f"tightness_intermediate.py — n={n}, 2026-05-22")
    print("=" * 72)

    P, T, W_set, cw, W_corners, corners = classify_edges(n)
    c_TL, c_TR, c_BL, c_BR = corners
    print(f"V={n*n}, E_torus={len(T)}, E_plane={len(P)}, "
          f"|W|={len(W_set)}, |W_corners|={len(W_corners)}")
    print(f"corners: TL={c_TL}, TR={c_TR}, BL={c_BL}, BR={c_BR}")
    print(f"per corner |W(c)|: TL={len(cw[c_TL])}, TR={len(cw[c_TR])}, "
          f"BL={len(cw[c_BL])}, BR={len(cw[c_BR])}")

    edge_list, edge_to_idx, E_torus = edge_indexing(n)

    # Configurações
    cw_TL = cw[c_TL]
    cw_TR = cw[c_TR]
    cw_BL = cw[c_BL]
    cw_BR = cw[c_BR]
    cw_TL_list = sorted(cw_TL)

    configs = [
        ("TORO",  set()),
        # Grupo A: k=1
        ("A1",    cw_TL),
        ("A2",    cw_BR),
        # Grupo B: k=2
        ("B1",    cw_TL | cw_TR),     # adjacentes
        ("B2",    cw_TL | cw_BR),     # opostos
        # Grupo C: k=3
        ("C1",    cw_TL | cw_TR | cw_BL),
        ("C2",    cw_TL | cw_TR | cw_BR),
        # Grupo D: k=4 (plano)
        ("D1",    W_set),
        # Grupo E: parciais (k=0)
        ("E1",    {cw_TL_list[0]}),
        ("E2",    set(cw_TL_list[:3])),
    ]

    K_per_seed = 1000
    n_seeds = 20

    results = []
    t_total0 = time.time()
    for name, forbidden in configs:
        r = test_config(name, forbidden, P, W_set, n,
                        K_per_seed, n_seeds,
                        edge_list, edge_to_idx, E_torus, log=False)
        results.append(r)

    # ── Tabela final ────────────────────────────────────────────────────────
    print(f"\n{'═'*72}")
    print(f"TABELA FINAL — Tightness na família híbrida G(H) (n={n})")
    print(f"{'═'*72}")
    print(f"{'Config':<7} {'k(H)':>4} {'Q':>3} {'β₁':>4} {'rank':>5} "
          f"{'def':>4} {'match':>6} {'tours':>5} {'conv':>5} {'time':>6}")
    print("─" * 72)
    for r in results:
        rank_s = str(r["rank"]) if r["rank"] is not None else "—"
        def_s = str(r["deficit"]) if r["deficit"] is not None else "—"
        match_s = ("✓" if r["match"] else "✗") if r["match"] is not None else "—"
        conv_s = ("Y" if r.get("converged") else "N") if r["rank"] is not None else "—"
        print(f"{r['name']:<7} {r['k_H']:>4} {r['Q']:>3} {r['beta1']:>4} "
              f"{rank_s:>5} {def_s:>4} {match_s:>6} {r['n_tours']:>5} "
              f"{conv_s:>5} {r['time']:>5.1f}s")
    print("─" * 72)

    all_match = all(r["match"] for r in results if r["match"] is not None)
    print(f"\nTodos os pontos batem deficit = Q: {all_match}")
    print(f"Tempo total: {time.time()-t_total0:.1f}s")

    return results


if __name__ == "__main__":
    main()
