#!/usr/bin/env python3
"""
benchmark_10x10.py
==================
Compara três configurações de Z3 para amostrar tours hamiltonianos 10×10:

  A) Z3 puro            — só grau-2 + pós-filtro BFS
  B) + mandatory         — A com arestas obrigatórias fixadas em 1
  C) + mandatory + ¬(A∧B) — B com cláusulas NOT(A∧B) dos candidatos

Métricas (K=200 amostras por config):
  - time-to-first
  - tempo médio por solução
  - attempts/solution
  - taxa de rejeição

Saídas:
  data/benchmark_10x10.json
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
INV = DATA / "invariants"

sys.path.insert(0, str(ROOT))
from graph_10x10 import build_graph, TOTAL  # noqa: E402
from sample_tours import is_single_tour  # noqa: E402


def build_solver(edges, *, mandatory_idx=None, forbid_pairs=None):
    """
    edges: lista de arestas (u, v)
    mandatory_idx: índices de arestas a fixar em 1
    forbid_pairs: lista de pares (i, j) de índices com restrição ¬(x_i ∧ x_j)
    """
    E = len(edges)
    xvars = [Bool(f"x_{i}") for i in range(E)]
    s = Solver()

    inc = {v: [] for v in range(TOTAL)}
    for i, (u, v) in enumerate(edges):
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

    return s, xvars


def run_config(edges, name, *, target_k=200, seed=2026,
               mandatory_idx=None, forbid_pairs=None,
               timeout_per_check_s=30, verbose=False):
    """Executa um benchmark e retorna métricas."""
    rng = random.Random(seed)
    E = len(edges)

    s, xvars = build_solver(edges,
                            mandatory_idx=mandatory_idx,
                            forbid_pairs=forbid_pairs)

    n_valid = 0
    n_attempts = 0
    n_rejected = 0
    time_first = None
    solution_times = []

    t0 = time.perf_counter()
    while n_valid < target_k:
        n_attempts += 1
        s.push()
        # quebra de simetria aleatória
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
        for i, (u, v) in enumerate(edges):
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
        "n_rejected_subtours": n_rejected,
        "rejection_rate": round(n_rejected / max(n_attempts, 1), 4),
        "elapsed_s": round(elapsed, 3),
        "time_to_first_s": round(time_first, 4) if time_first else None,
        "mean_s_per_valid": round(elapsed / n_valid, 4) if n_valid else None,
        "attempts_per_valid": round(n_attempts / n_valid, 4) if n_valid else None,
        "solutions_per_second": round(n_valid / elapsed, 4) if elapsed > 0 else None,
        "mandatory_applied": len(mandatory_idx) if mandatory_idx else 0,
        "forbid_pairs_applied": len(forbid_pairs) if forbid_pairs else 0,
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--target-k", type=int, default=200)
    p.add_argument("--seed", type=int, default=2026)
    p.add_argument("--max-pairs", type=int, default=20,
                   help="Limita pares NOT(A∧B) injetados (top por menor coexist)")
    args = p.parse_args()

    _, edges = build_graph()

    # carrega mandatory e candidates dos passos anteriores
    mand_path = INV / "mandatory_edges.json"
    cand_path = INV / "candidate_clauses.json"

    mandatory_idx = []
    if mand_path.exists():
        d = json.loads(mand_path.read_text())
        mandatory_idx = [e["idx"] for e in d["mandatory"]]
        print(f"mandatory carregadas: {len(mandatory_idx)}")

    forbid_pairs = []
    if cand_path.exists():
        d = json.loads(cand_path.read_text())
        cands = d["candidates"][: args.max_pairs]
        forbid_pairs = [(c["edge_a_idx"], c["edge_b_idx"]) for c in cands]
        print(f"forbid_pairs carregados: {len(forbid_pairs)} "
              f"(top {args.max_pairs} por coexist)")

    results = []

    print("\n[A] Z3 puro ...")
    rA = run_config(edges, "A_puro", target_k=args.target_k,
                    seed=args.seed, verbose=True)
    results.append(rA)
    print(f"  elapsed={rA['elapsed_s']}s  t_first={rA['time_to_first_s']}s  "
          f"sol/s={rA['solutions_per_second']}  "
          f"reject={rA['rejection_rate']:.1%}")

    print("\n[B] Z3 + mandatory ...")
    rB = run_config(edges, "B_mandatory", target_k=args.target_k,
                    seed=args.seed, mandatory_idx=mandatory_idx, verbose=True)
    results.append(rB)
    print(f"  elapsed={rB['elapsed_s']}s  t_first={rB['time_to_first_s']}s  "
          f"sol/s={rB['solutions_per_second']}  "
          f"reject={rB['rejection_rate']:.1%}")

    print("\n[C] Z3 + mandatory + NOT(A∧B) ...")
    rC = run_config(edges, "C_full_invariants", target_k=args.target_k,
                    seed=args.seed,
                    mandatory_idx=mandatory_idx,
                    forbid_pairs=forbid_pairs, verbose=True)
    results.append(rC)
    print(f"  elapsed={rC['elapsed_s']}s  t_first={rC['time_to_first_s']}s  "
          f"sol/s={rC['solutions_per_second']}  "
          f"reject={rC['rejection_rate']:.1%}")

    # speedups
    speedups = {}
    for r in results[1:]:
        speedups[r["config"]] = round(
            rA["elapsed_s"] / r["elapsed_s"], 3) if r["elapsed_s"] else None

    out_path = DATA / "benchmark_10x10.json"
    with open(out_path, "w") as f:
        json.dump({
            "args": vars(args),
            "results": results,
            "speedups_vs_A": speedups,
        }, f, indent=2)
    print(f"\nSalvo: {out_path.relative_to(ROOT)}")

    print("\n─" * 30)
    print("RESUMO")
    print("─" * 30)
    print(f"  {'config':22s} {'t_total':>10s} {'t_first':>10s} "
          f"{'sol/s':>10s} {'att/sol':>10s} {'reject':>10s}  speedup")
    for r in results:
        sp = speedups.get(r["config"])
        print(f"  {r['config']:22s} {r['elapsed_s']:>10.2f} "
              f"{r['time_to_first_s'] or 0:>10.4f} "
              f"{r['solutions_per_second'] or 0:>10.3f} "
              f"{r['attempts_per_valid'] or 0:>10.3f} "
              f"{r['rejection_rate']:>10.2%}  "
              f"{sp if sp else '—':>6}")


if __name__ == "__main__":
    main()
