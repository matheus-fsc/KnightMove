#!/usr/bin/env python3
"""benchmark_finf.py — Benchmark BT_base vs BT_theory (heurística f∞(L)).

Hipótese: ganho do f∞ ∝ fração de vértices com L ≥ 4 = max(0, (n-8)²)/n²
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
import statistics

import numpy as np

sys.path.insert(0, "/home/math/Dev/knight_tour")
sys.setrecursionlimit(20000)

from knight_tours import (
    build_graph, State, edge_level, edge_priority,
    FREE, ACTIVE, INACTIVE, OK, CONTRADICTION, SUBTOUR, COMPLETE_TOUR,
    _process_queue, _propagate_initial, _is_complete_tour, _extract_tour,
    verify_tour,
)


# ── seleção SEM f∞ (BT_base) ──────────────────────────────────────────────────

def _choose_next_edge_base(state, ctx, rng):
    """Pressão de vértice mas sem heurística f∞ — primeira aresta livre,
    valor inicial fixo (0 = INACTIVE primeiro, idem ao default global)."""
    V = ctx["V"]
    total_inc = ctx["total_incident"]
    adj_edges = ctx["adj_edges"]

    best_score = -1
    candidates = []
    for v in range(V):
        free_v = int(total_inc[v]) - int(state.degree[v]) - int(state.n_inactive[v])
        if free_v == 0:
            continue
        score = int(state.degree[v]) * 100 + \
                int(total_inc[v] - state.n_inactive[v] - 2)
        if score > best_score:
            best_score = score
            candidates = [v]
        elif score == best_score:
            candidates.append(v)

    if not candidates:
        return -1, 0

    if len(candidates) == 1:
        v_star = candidates[0]
    else:
        v_star = int(candidates[int(rng.integers(len(candidates)))])

    # primeira aresta livre incidente a v_star — SEM heurística f∞
    for ei in adj_edges[v_star]:
        ei_i = int(ei)
        if state.fixed[ei_i] == FREE:
            return ei_i, 0  # val=0 default (igual a F_INF_DEFAULT < 0.5)

    return -1, 0


# ── seleção COM f∞ (BT_theory) ────────────────────────────────────────────────

def _choose_next_edge_theory(state, ctx, rng):
    """Idêntico ao de knight_tours.py (pressão de vértice + f∞(L))."""
    V = ctx["V"]
    n = ctx["n"]
    total_inc = ctx["total_incident"]
    adj_edges = ctx["adj_edges"]
    edge_endpoints = ctx["edge_endpoints"]

    best_score = -1
    candidates = []
    for v in range(V):
        free_v = int(total_inc[v]) - int(state.degree[v]) - int(state.n_inactive[v])
        if free_v == 0:
            continue
        score = int(state.degree[v]) * 100 + \
                int(total_inc[v] - state.n_inactive[v] - 2)
        if score > best_score:
            best_score = score
            candidates = [v]
        elif score == best_score:
            candidates.append(v)

    if not candidates:
        return -1, 0

    if len(candidates) == 1:
        v_star = candidates[0]
    else:
        v_star = int(candidates[int(rng.integers(len(candidates)))])

    best_pri = -1.0
    best_edge = -1
    best_first = 0
    for ei in adj_edges[v_star]:
        ei_i = int(ei)
        if state.fixed[ei_i] != FREE:
            continue
        a = int(edge_endpoints[ei_i, 0])
        b = int(edge_endpoints[ei_i, 1])
        L = edge_level(a, b, n)
        pri, fv = edge_priority(L)
        if pri > best_pri:
            best_pri = pri
            best_edge = ei_i
            best_first = fv

    return best_edge, best_first


# ── backtracking com contador ────────────────────────────────────────────────

def _backtrack_with_count(state, ctx, rng, tours, K, choose_fn, counters):
    if len(tours) >= K:
        return True

    counters["nodes"] += 1
    e, first_val = choose_fn(state, ctx, rng)

    if e == -1:
        if _is_complete_tour(state, ctx):
            tours.append(_extract_tour(state, ctx))
        return len(tours) >= K

    for val in (first_val, 1 - first_val):
        snap = state.snapshot()
        # fix_and_propagate inline (evita import circular)
        status = _process_queue(state, ctx, [(e, val)])
        if status == OK:
            if _backtrack_with_count(state, ctx, rng, tours, K, choose_fn, counters):
                state.restore(snap)
                return True
        elif status == COMPLETE_TOUR:
            if _is_complete_tour(state, ctx):
                tours.append(_extract_tour(state, ctx))
            if len(tours) >= K:
                state.restore(snap)
                return True
        state.restore(snap)

    return len(tours) >= K


def knight_tours_variant(n, K, seed, use_theory):
    """Runs backtracking with or without f∞. Returns (tours, n_nodes, t_total)."""
    ctx = build_graph(n)
    rng = np.random.default_rng(seed)
    state = State(ctx["V"], ctx["E"])

    t0 = time.perf_counter()
    status = _propagate_initial(state, ctx)
    if status in (CONTRADICTION, SUBTOUR):
        return [], 0, time.perf_counter() - t0

    tours = []
    counters = {"nodes": 0}

    if status == COMPLETE_TOUR:
        if _is_complete_tour(state, ctx):
            tours.append(_extract_tour(state, ctx))
    else:
        choose_fn = _choose_next_edge_theory if use_theory else _choose_next_edge_base
        _backtrack_with_count(state, ctx, rng, tours, K, choose_fn, counters)

    t_total = time.perf_counter() - t0
    return tours[:K], counters["nodes"], t_total


# ── execução do benchmark ────────────────────────────────────────────────────

def fracao_L_ge_4(n):
    """Fração de vértices com L ≥ 4 num tabuleiro n×n. Vértice (r,c) tem
    L ≥ 4 ⟺ r,c ∈ [4, n-5] ⟺ n ≥ 9 com (n-8)² vértices interiores."""
    return max(0, n - 8) ** 2 / (n * n)


def run_benchmark(ns=(6, 8, 10, 12, 14, 16), K=200, n_runs=3, base_seed=2026):
    print("=" * 78)
    print(f"BENCHMARK BT_base vs BT_theory (f∞)  —  K={K}, runs={n_runs}")
    print("=" * 78)

    results = {}

    for n in ns:
        frac = fracao_L_ge_4(n)
        speedup_pred = 1.0 + frac  # hipótese linear
        print(f"\n── n={n}  (fração L≥4 = {100*frac:.1f}%, "
              f"speedup_pred ~ {speedup_pred:.2f}) ──")
        t_base_list, t_theory_list = [], []
        nodes_base_list, nodes_theory_list = [], []

        for run in range(n_runs):
            seed = base_seed + 1000 * run + n
            # rodar BASE
            tours_b, nodes_b, tb = knight_tours_variant(n, K, seed, use_theory=False)
            # rodar THEORY (mesma seed → mesmas escolhas de v_star quando empate)
            tours_t, nodes_t, tt = knight_tours_variant(n, K, seed, use_theory=True)
            assert len(tours_b) == K, f"BT_base não atingiu K em n={n}, run={run}"
            assert len(tours_t) == K, f"BT_theory não atingiu K em n={n}, run={run}"
            # verificação cruzada (sample)
            assert verify_tour(tours_b[0], n), f"BT_base produziu tour inválido em n={n}"
            assert verify_tour(tours_t[0], n), f"BT_theory produziu tour inválido em n={n}"

            t_base_list.append(tb)
            t_theory_list.append(tt)
            nodes_base_list.append(nodes_b)
            nodes_theory_list.append(nodes_t)
            print(f"  run {run+1}: t_base={tb:.3f}s ({nodes_b/K:.2f} nós/tour)  "
                  f"t_theory={tt:.3f}s ({nodes_t/K:.2f} nós/tour)  "
                  f"speedup={tb/tt:.3f}")

        tb_mean = statistics.mean(t_base_list)
        tb_std = statistics.stdev(t_base_list) if n_runs > 1 else 0
        tt_mean = statistics.mean(t_theory_list)
        tt_std = statistics.stdev(t_theory_list) if n_runs > 1 else 0
        nb_mean = statistics.mean(nodes_base_list) / K
        nt_mean = statistics.mean(nodes_theory_list) / K
        speedup = tb_mean / tt_mean

        all_speedups = [tb / tt for tb, tt in zip(t_base_list, t_theory_list)]
        speedup_median = statistics.median(all_speedups)

        results[n] = {
            "frac": frac,
            "speedup_pred": speedup_pred,
            "t_base": (tb_mean, tb_std),
            "t_theory": (tt_mean, tt_std),
            "nodes_base": nb_mean,
            "nodes_theory": nt_mean,
            "speedup_mean": speedup,
            "speedup_median": speedup_median,
            "all_speedups": all_speedups,
        }
        print(f"  → t_base   = {tb_mean:.3f} ± {tb_std:.3f} s")
        print(f"  → t_theory = {tt_mean:.3f} ± {tt_std:.3f} s")
        print(f"  → speedup  mean={speedup:.3f}  median={speedup_median:.3f}  "
              f"(predito ~{speedup_pred:.3f})")

    # ── Tabela final ────────────────────────────────────────────────────────
    print()
    print("=" * 90)
    print("TABELA FINAL")
    print("=" * 90)
    hdr = (f"{'n':>3} | {'frac':>6} | {'t_base':>14} | {'t_theory':>14} | "
           f"{'nó/tour':>14} | {'su_med':>7} | {'su_avg':>7} | {'pred':>6} | {'match':>5}")
    print(hdr)
    print("─" * 100)
    for n in ns:
        r = results[n]
        tb_str = f"{r['t_base'][0]:.3f}±{r['t_base'][1]:.3f}"
        tt_str = f"{r['t_theory'][0]:.3f}±{r['t_theory'][1]:.3f}"
        nodes_str = f"{r['nodes_base']:.2f}/{r['nodes_theory']:.2f}"
        # match: speedup mediano dentro de ±20% do predito
        ratio = r["speedup_median"] / r["speedup_pred"]
        match = "≈" if 0.8 <= ratio <= 1.2 else ("↑" if ratio > 1.2 else "↓")
        print(f"{n:>3} | {100*r['frac']:>5.1f}% | {tb_str:>14} | {tt_str:>14} | "
              f"{nodes_str:>14} | {r['speedup_median']:>7.3f} | "
              f"{r['speedup_mean']:>7.3f} | {r['speedup_pred']:>6.3f} | "
              f"{match:>5}")
    print("─" * 100)

    return results


if __name__ == "__main__":
    # K=200 + 5 runs (outliers ainda existem; reportar mediana)
    results = run_benchmark(ns=(6, 8, 10, 12, 14, 16), K=200, n_runs=5)
