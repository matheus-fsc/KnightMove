#!/usr/bin/env python3
"""
xor_10x10.py
============
T3: benchmark Z3 do 10×10 com 5 configurações:

  A: Z3 puro (grau-2 + sub-tour elimination iterativo)
  B: A + 8 obrigatórias dos cantos
  C: B + 3 cláusulas XOR (base mínima do quociente, calculada em
     find_xor_10x10.py: data/xor_10x10_clauses.json)
  D: B + 28 XOR (todos os pares de obrigatórias)
  E: B + 3 XOR + NOT(A∧B) para top-K pares mais correlacionados

K amostras por config; cada sample tem timeout configurável.
"""

from __future__ import annotations

import argparse
import itertools
import json
import random
import sys
import time
from pathlib import Path

import numpy as np
from z3 import Bool, If, Not, Or, And, Sum, PbEq, Solver, sat, is_true

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
REPO = ROOT.parent

sys.path.insert(0, str(REPO / "board_10x10"))
from graph_10x10 import build_graph, TOTAL  # noqa: E402


# ── helpers reaproveitados ────────────────────────────────────────────

def is_single_tour(active_edges, n_vertices=TOTAL):
    adj = [[] for _ in range(n_vertices)]
    for u, v in active_edges:
        adj[u].append(v)
        adj[v].append(u)
    if any(len(a) != 2 for a in adj):
        return False
    visited = [False] * n_vertices
    visited[0] = True
    cur, prev = 0, -1
    for _ in range(n_vertices - 1):
        nxt = adj[cur][0] if adj[cur][0] != prev else adj[cur][1]
        if visited[nxt]:
            return False
        visited[nxt] = True
        prev, cur = cur, nxt
    return adj[cur][0] == 0 or adj[cur][1] == 0


def build_base_solver(edges):
    E = len(edges)
    xvars = [Bool(f"x_{i}") for i in range(E)]
    inc = {v: [] for v in range(TOTAL)}
    for i, (u, v) in enumerate(edges):
        inc[u].append(i)
        inc[v].append(i)
    s = Solver()
    for v, lst in inc.items():
        s.add(PbEq([(xvars[i], 1) for i in lst], 2))
    return s, xvars


def add_xor_clause(solver, xvars, support):
    if len(support) == 2:
        a, b = support
        solver.add(xvars[a] == xvars[b])
        return
    terms = [If(xvars[e], 1, 0) for e in support]
    solver.add(Sum(terms) % 2 == 0)


def build_solver_for_config(cfg, edges, mand, xor_min, xor_all,
                            top_excl_pairs):
    s, xvars = build_base_solver(edges)
    if cfg != "A":
        for e in mand:
            s.add(xvars[e])
    if cfg == "C":
        for supp in xor_min:
            add_xor_clause(s, xvars, supp)
    if cfg == "D":
        for supp in xor_all:
            add_xor_clause(s, xvars, supp)
    if cfg == "E":
        for supp in xor_min:
            add_xor_clause(s, xvars, supp)
        for (i, j) in top_excl_pairs:
            s.add(Not(And(xvars[i], xvars[j])))
    return s, xvars


def sample_k_tours(solver, xvars, edges, K, *, seed, timeout_ms=60_000,
                   max_sym_tries=4):
    rng = random.Random(seed)
    E = len(edges)
    sigs = []
    stats = {"K": K, "n_attempts": 0, "n_subtour_rejects": 0,
             "n_local_unsat": 0, "t_first": None, "t_total": 0.0,
             "times_per_valid": []}
    t0 = time.perf_counter()
    while len(sigs) < K:
        stats["n_attempts"] += 1
        t_attempt = time.perf_counter()
        solver.push()
        local_unsat = True
        for _ in range(max_sym_tries):
            e_idx = rng.randrange(E)
            val = rng.randrange(2)
            solver.push()
            solver.add(xvars[e_idx] == (val == 1))
            solver.set("timeout", timeout_ms)
            res = solver.check()
            if res == sat:
                local_unsat = False
                break
            solver.pop()
        if local_unsat:
            stats["n_local_unsat"] += 1
            solver.pop()
            continue
        m = solver.model()
        sig = np.zeros(E, dtype=np.uint8)
        active = []
        for i, (u, v) in enumerate(edges):
            if is_true(m.evaluate(xvars[i])):
                sig[i] = 1
                active.append((u, v))
        ok = is_single_tour(active)
        solver.pop()
        solver.pop()
        if ok:
            sigs.append(sig)
            dt = time.perf_counter() - t_attempt
            stats["times_per_valid"].append(round(dt, 4))
            if stats["t_first"] is None:
                stats["t_first"] = round(time.perf_counter() - t0, 4)
            solver.add(Or([xvars[i] != bool(sig[i]) for i in range(E)]))
        else:
            stats["n_subtour_rejects"] += 1
            solver.add(Or([xvars[i] != bool(sig[i]) for i in range(E)]))
    stats["t_total"] = round(time.perf_counter() - t0, 3)
    stats["sol_per_s"] = round(K / stats["t_total"], 3) if stats["t_total"] > 0 else None
    stats["attempts_per_sol"] = round(stats["n_attempts"] / K, 3)
    stats["mean_time_per_valid"] = round(
        sum(stats["times_per_valid"]) / len(stats["times_per_valid"]), 4)
    sigs = np.stack(sigs, axis=0)
    return sigs, stats


def top_corr_pairs_10x10(K_pairs=8):
    """Top-K pares STRICTAMENTE excluídos (p_coexist=0).

    No 10×10 atual NÃO existem pares com p_coexist=0 (verificado em
    edge_pair_correlations.json: o mais negativo tem p_coexist=0.10).
    Retorna [] e a Config E acaba sendo equivalente a Config C.
    """
    p = REPO / "board_10x10" / "data" / "invariants" / "edge_pair_correlations.json"
    if not p.exists():
        return []
    with open(p) as f:
        d = json.load(f)
    items = d.get("candidates", [])
    strict = []
    for item in items:
        if float(item.get("p_coexist", 1.0)) == 0.0:
            a = item.get("edge_a_idx")
            b = item.get("edge_b_idx")
            if a is not None and b is not None:
                strict.append((int(a), int(b)))
        if len(strict) >= K_pairs:
            break
    return strict[:K_pairs]


def freq_per_edge(sigs):
    return sigs.mean(axis=0)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--K", type=int, default=200)
    p.add_argument("--seed", type=int, default=2026)
    p.add_argument("--configs", type=str, default="A,B,C,D,E")
    p.add_argument("--timeout-ms", type=int, default=60_000)
    args = p.parse_args()
    configs = args.configs.split(",")

    adj, edges = build_graph()
    E = len(edges)

    # carrega XOR identificadas
    xclauses = json.loads((DATA / "xor_10x10_clauses.json").read_text())
    mand = list(map(int, xclauses["mandatory_idx"]))
    xor_min = [list(map(int, p)) for p in xclauses["chosen_pairs"]]
    xor_all = list(itertools.combinations(mand, 2))

    # top-K NOT(A∧B) — pares mais informativos
    top_excl = top_corr_pairs_10x10(K_pairs=8)
    if not top_excl:
        print("  (nenhum par de exclusão NOT(A∧B) disponível para Config E — "
              "vai usar 0)")

    print(f"V={TOTAL}  E={E}  mandatórias={len(mand)}  "
          f"XOR mínimas={len(xor_min)}  XOR todas={len(xor_all)}  "
          f"NOT pares={len(top_excl)}")
    print(f"K={args.K}  seed={args.seed}  timeout_per_sample={args.timeout_ms}ms")

    # baseline empírico para KL
    samples_dir = REPO / "board_10x10" / "data" / "samples"
    batches = sorted(samples_dir.glob("tours_10x10_batch_*.npy"))
    T = np.concatenate([np.load(b) for b in batches], axis=0)
    p_truth = T.mean(axis=0)

    results = {}
    for cfg in configs:
        print(f"\n— Config {cfg} —")
        s, xvars = build_solver_for_config(cfg, edges, mand, xor_min, xor_all, top_excl)
        t0 = time.perf_counter()
        sigs, stats = sample_k_tours(s, xvars, edges, args.K,
                                     seed=args.seed,
                                     timeout_ms=args.timeout_ms)
        dt = time.perf_counter() - t0
        p_smp = freq_per_edge(sigs)
        # KL Bernoulli por aresta
        eps = 1e-6
        p_ = np.clip(p_smp, eps, 1 - eps)
        q_ = np.clip(p_truth, eps, 1 - eps)
        kl = float(np.sum(p_ * np.log(p_ / q_) + (1 - p_) * np.log((1 - p_) / (1 - q_))))
        stats["kl_vs_truth"] = round(kl, 4)
        stats["wall_total"] = round(dt, 3)
        results[cfg] = stats
        print(f"  t_total={stats['t_total']}s  sol/s={stats['sol_per_s']}  "
              f"attempts/sol={stats['attempts_per_sol']}  KL={stats['kl_vs_truth']:.3f}")

    # speedups
    if "A" in results:
        baseA = results["A"]["t_total"]
        for cfg, st in results.items():
            st["speedup_vs_A"] = round(baseA / st["t_total"], 3)
    if "B" in results:
        baseB = results["B"]["t_total"]
        for cfg, st in results.items():
            st["speedup_vs_B"] = round(baseB / st["t_total"], 3)

    out = {
        "board": 10,
        "K": args.K,
        "seed": args.seed,
        "configs": list(configs),
        "n_mandatory": len(mand),
        "mandatory_idx": mand,
        "xor_min_pairs": xor_min,
        "n_xor_min": len(xor_min),
        "n_xor_all_pairs": len(xor_all),
        "top_excl_pairs": [list(map(int, p)) for p in top_excl],
        "results": results,
    }
    DATA.mkdir(parents=True, exist_ok=True)
    outp = DATA / "benchmark_10x10_xor.json"
    with open(outp, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nSalvo: {outp.relative_to(REPO)}")


if __name__ == "__main__":
    main()
