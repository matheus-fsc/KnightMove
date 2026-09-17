#!/usr/bin/env python3
"""
benchmark_triples.py
====================
Compara 5 configurações de Z3 para amostrar tours hamiltonianos 10×10:

  A — Z3 puro (grau-2 + BFS)
  B — A + 8 mandatory edges
  C — B + NOT(A∧B) dos top-8 pares (do benchmark anterior)
  D — B + NOT(A∧B∧C) das triplas PROVEN
  E — B + NOT(A∧B) top-8 + NOT(A∧B∧C) proven (combinado)

Métricas (K=200 por config):
  t_first, t_total, sol/s, attempts/sol, rejection_rate

Saída:
  data/benchmark_triples.json
  data/plots/benchmark_comparison.png
"""

import argparse
import json
import random
import sys
import time
from pathlib import Path

import numpy as np
from z3 import Bool, Or, PbEq, Solver, sat, is_true, Not, And

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
PLOTS = DATA / "plots"
PARENT = ROOT.parent

sys.path.insert(0, str(PARENT))
from graph_10x10 import build_graph, TOTAL  # noqa: E402
from sample_tours import is_single_tour  # noqa: E402

ADJ, EDGES = build_graph()


def build_solver(*, mandatory_idx=None, forbid_pairs=None, forbid_triples=None):
    E = len(EDGES)
    xvars = [Bool(f"x_{i}") for i in range(E)]
    s = Solver()

    inc = {v: [] for v in range(TOTAL)}
    for i, (u, v) in enumerate(EDGES):
        inc[u].append(i)
        inc[v].append(i)
    for v, lst in inc.items():
        s.add(PbEq([(xvars[i], 1) for i in lst], 2))

    if mandatory_idx:
        for i in mandatory_idx:
            s.add(xvars[i] == True)

    if forbid_pairs:
        for i, j in forbid_pairs:
            s.add(Not(And(xvars[i], xvars[j])))

    if forbid_triples:
        for i, j, k in forbid_triples:
            s.add(Not(And(xvars[i], xvars[j], xvars[k])))

    return s, xvars


def run_config(name, *, target_k=200, seed=2026,
               mandatory_idx=None, forbid_pairs=None, forbid_triples=None,
               timeout_per_check_s=30, verbose=False):
    rng = random.Random(seed)
    E = len(EDGES)

    s, xvars = build_solver(mandatory_idx=mandatory_idx,
                            forbid_pairs=forbid_pairs,
                            forbid_triples=forbid_triples)

    n_valid = 0
    n_attempts = 0
    n_rejected = 0
    time_first = None
    solution_times = []
    t0 = time.perf_counter()

    while n_valid < target_k:
        n_attempts += 1
        s.push()
        edge_pick = rng.randrange(E)
        val_pick = rng.randrange(2)
        s.add(xvars[edge_pick] == (val_pick == 1))
        s.set("timeout", int(timeout_per_check_s * 1000))
        res = s.check()
        if res != sat:
            s.pop()
            continue

        m = s.model()
        active = []
        bits = np.zeros(E, dtype=np.uint8)
        for i, (u, v) in enumerate(EDGES):
            if is_true(m.evaluate(xvars[i])):
                active.append((u, v))
                bits[i] = 1
        ok, _ = is_single_tour(active)
        s.pop()

        if ok:
            now = time.perf_counter() - t0
            if time_first is None:
                time_first = now
            solution_times.append(now)
            n_valid += 1
            s.add(Or([xvars[i] != bool(bits[i]) for i in range(E)]))
            if verbose and n_valid % 25 == 0:
                print(f"    [{name}] {n_valid}/{target_k}  "
                      f"attempts={n_attempts}  elapsed={now:.1f}s",
                      flush=True)
        else:
            n_rejected += 1
            s.add(Or([xvars[i] != bool(bits[i]) for i in range(E)]))

    elapsed = time.perf_counter() - t0
    return {
        "config": name,
        "target_k": target_k,
        "n_valid": n_valid,
        "n_attempts": n_attempts,
        "n_rejected": n_rejected,
        "rejection_rate": round(n_rejected / max(n_attempts, 1), 4),
        "elapsed_s": round(elapsed, 3),
        "time_to_first_s": round(time_first, 4) if time_first else None,
        "mean_s_per_valid": round(elapsed / n_valid, 4) if n_valid else None,
        "attempts_per_valid": round(n_attempts / n_valid, 4) if n_valid else None,
        "solutions_per_second": round(n_valid / elapsed, 4) if elapsed else None,
        "n_mandatory": len(mandatory_idx) if mandatory_idx else 0,
        "n_forbid_pairs": len(forbid_pairs) if forbid_pairs else 0,
        "n_forbid_triples": len(forbid_triples) if forbid_triples else 0,
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--target-k", type=int, default=200)
    p.add_argument("--seed", type=int, default=2026)
    p.add_argument("--max-pairs", type=int, default=8)
    p.add_argument("--max-triples", type=int, default=None,
                   help="Limita triplas proven injetadas (None=todas)")
    args = p.parse_args()

    PLOTS.mkdir(parents=True, exist_ok=True)

    # ── carrega dados ────────────────────────────────────────────────
    mand = json.loads((PARENT / "data" / "invariants" /
                       "mandatory_edges.json").read_text())
    mandatory_idx = [e["idx"] for e in mand["mandatory"]]
    print(f"mandatory: {len(mandatory_idx)}")

    cand = json.loads((PARENT / "data" / "invariants" /
                       "candidate_clauses.json").read_text())
    pair_clauses = cand["candidates"][: args.max_pairs]
    forbid_pairs = [(c["edge_a_idx"], c["edge_b_idx"]) for c in pair_clauses]
    print(f"forbid_pairs (top {args.max_pairs}): {len(forbid_pairs)}")

    ver = json.loads((DATA / "triples_verified.json").read_text())
    proven = [r for r in ver["results"] if r["z3_status"] == "proven"]
    if args.max_triples:
        proven = proven[: args.max_triples]
    forbid_triples = [tuple(r["edges"]) for r in proven]
    print(f"forbid_triples (proven): {len(forbid_triples)}")

    if not forbid_triples:
        print("\nERRO: zero triplas proven — não há nada para D/E. "
              "Rode z3_verify_triples.py primeiro.")
        return 1

    # ── 5 configurações ─────────────────────────────────────────────
    results = []

    print("\n[A] Z3 puro ...")
    rA = run_config("A_puro", target_k=args.target_k, seed=args.seed,
                    verbose=True)
    results.append(rA)
    print(f"  t_total={rA['elapsed_s']}s  t_first={rA['time_to_first_s']}s  "
          f"sol/s={rA['solutions_per_second']}  reject={rA['rejection_rate']:.1%}")

    print("\n[B] + mandatory ...")
    rB = run_config("B_mandatory", target_k=args.target_k, seed=args.seed,
                    mandatory_idx=mandatory_idx, verbose=True)
    results.append(rB)
    print(f"  t_total={rB['elapsed_s']}s  t_first={rB['time_to_first_s']}s  "
          f"sol/s={rB['solutions_per_second']}  reject={rB['rejection_rate']:.1%}")

    print(f"\n[C] + mandatory + NOT(A∧B) top-{args.max_pairs} ...")
    rC = run_config("C_pairs", target_k=args.target_k, seed=args.seed,
                    mandatory_idx=mandatory_idx,
                    forbid_pairs=forbid_pairs, verbose=True)
    results.append(rC)
    print(f"  t_total={rC['elapsed_s']}s  t_first={rC['time_to_first_s']}s  "
          f"sol/s={rC['solutions_per_second']}  reject={rC['rejection_rate']:.1%}")

    print(f"\n[D] + mandatory + NOT(A∧B∧C) ({len(forbid_triples)} triples proven) ...")
    rD = run_config("D_triples", target_k=args.target_k, seed=args.seed,
                    mandatory_idx=mandatory_idx,
                    forbid_triples=forbid_triples, verbose=True)
    results.append(rD)
    print(f"  t_total={rD['elapsed_s']}s  t_first={rD['time_to_first_s']}s  "
          f"sol/s={rD['solutions_per_second']}  reject={rD['rejection_rate']:.1%}")

    print(f"\n[E] + mandatory + pairs + triples (combinado) ...")
    rE = run_config("E_combined", target_k=args.target_k, seed=args.seed,
                    mandatory_idx=mandatory_idx,
                    forbid_pairs=forbid_pairs,
                    forbid_triples=forbid_triples, verbose=True)
    results.append(rE)
    print(f"  t_total={rE['elapsed_s']}s  t_first={rE['time_to_first_s']}s  "
          f"sol/s={rE['solutions_per_second']}  reject={rE['rejection_rate']:.1%}")

    # ── speedups ─────────────────────────────────────────────────────
    speedups = {}
    for r in results[1:]:
        speedups[r["config"]] = (round(rA["elapsed_s"] / r["elapsed_s"], 3)
                                  if r["elapsed_s"] else None)

    with open(DATA / "benchmark_triples.json", "w") as f:
        json.dump({
            "args": vars(args),
            "results": results,
            "speedups_vs_A": speedups,
        }, f, indent=2)
    print(f"\nSalvo: data/benchmark_triples.json")

    # ── plot ─────────────────────────────────────────────────────────
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        names = [r["config"] for r in results]
        t_total = [r["elapsed_s"] for r in results]
        t_first = [r["time_to_first_s"] or 0 for r in results]

        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        x = np.arange(len(names))

        axes[0].bar(x, t_total, color=["#888", "#5b9", "#5b5", "#e85", "#a37"])
        axes[0].set_xticks(x)
        axes[0].set_xticklabels(names, rotation=15)
        axes[0].set_ylabel("t_total (s)")
        axes[0].set_title(f"K={args.target_k} samples — tempo total")
        axes[0].axhline(t_total[0], color="red", linestyle="--", alpha=0.5,
                        label="baseline A")
        axes[0].legend()

        axes[1].bar(x, t_first, color=["#888", "#5b9", "#5b5", "#e85", "#a37"])
        axes[1].set_xticks(x)
        axes[1].set_xticklabels(names, rotation=15)
        axes[1].set_ylabel("t_first (s)")
        axes[1].set_title("Tempo para a primeira solução")
        axes[1].axhline(t_first[0], color="red", linestyle="--", alpha=0.5,
                        label="baseline A")
        axes[1].legend()

        fig.suptitle(f"Benchmark 10×10  —  pairs={len(forbid_pairs)} "
                     f"triples={len(forbid_triples)}")
        fig.tight_layout()
        plt.savefig(PLOTS / "benchmark_comparison.png", dpi=120)
        plt.close(fig)
        print(f"Salvo: data/plots/benchmark_comparison.png")
    except Exception as ex:
        print(f"AVISO: erro ao plotar: {ex}")

    # ── sumário ─────────────────────────────────────────────────────
    print(f"\n{'─' * 60}")
    print("RESUMO")
    print("─" * 60)
    print(f"  {'config':18s} {'t_total':>10s} {'t_first':>10s} "
          f"{'sol/s':>10s} {'reject':>8s}  speedup")
    for r in results:
        sp = speedups.get(r["config"])
        print(f"  {r['config']:18s} {r['elapsed_s']:>10.2f} "
              f"{r['time_to_first_s'] or 0:>10.4f} "
              f"{r['solutions_per_second'] or 0:>10.3f} "
              f"{r['rejection_rate']:>8.1%}  "
              f"{sp if sp else '—':>6}")


if __name__ == "__main__":
    sys.exit(main())
