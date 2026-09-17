"""
Benchmark reprodutível da Tabela 9 do paper — 10×10, K=200 tours.
3 repetições (mediana) para métodos BT; Z3 single-shot com timeout.
"""
from __future__ import annotations
import json, sys, time, random
from pathlib import Path
from statistics import median

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "incremental_subtour"))
sys.path.insert(0, str(ROOT / "residual_search_10x10"))
sys.path.insert(0, str(ROOT / "local_phase_heuristic"))

N = 10
K = 200
SEED = 42
TIMEOUT_Z3 = 300.0  # 5 min per Z3 method
TIMEOUT_BT = 300.0  # 5 min per BT method
REPS = 3


def run_z3_puro(K, timeout):
    from scaling_minimal_v2 import build_knight_graph
    from z3 import Bool, Not, Or, PbEq, Solver, Sum, sat
    from collections import deque

    edges_uv, V, v2e = build_knight_graph(N)
    E = len(edges_uv)

    x = [Bool(f"x_{i}") for i in range(E)]
    s = Solver()
    for vtx, es in v2e.items():
        s.add(PbEq([(x[e], 1) for e in es], 2))

    t0 = time.perf_counter()
    tours = []
    tentativas = 0
    while len(tours) < K:
        if time.perf_counter() - t0 > timeout:
            break
        tentativas += 1
        if s.check() != sat:
            break
        m = s.model()
        ativos = [e for e in range(E) if m[x[e]] is not None and bool(m[x[e]])]
        adj = [[] for _ in range(V)]
        for e in ativos:
            u, v = edges_uv[e]
            adj[u].append(v)
            adj[v].append(u)
        bfs = [0]
        seen = {0}
        qi = 0
        while qi < len(bfs):
            node = bfs[qi]; qi += 1
            for nb in adj[node]:
                if nb not in seen:
                    seen.add(nb)
                    bfs.append(nb)
        if len(seen) == V:
            tours.append(ativos)
            s.add(Sum([x[e] if e in ativos else Not(x[e]) for e in range(E)]) < E)
        else:
            ac = [e for e in range(E) if edges_uv[e][0] in seen and edges_uv[e][1] in seen
                  and m[x[e]] is not None and bool(m[x[e]])]
            if ac:
                s.add(Or([Not(x[e]) for e in ac]))
    t_total = time.perf_counter() - t0
    return {"tours": len(tours), "t_total": t_total, "tentativas": tentativas}


def run_z3_mandatory(K, timeout):
    from scaling_minimal_v2 import build_knight_graph, mandatory_corner_edges
    from z3 import Bool, Not, Or, PbEq, Solver, Sum, sat
    from collections import deque

    edges_uv, V, v2e = build_knight_graph(N)
    E = len(edges_uv)
    mand = mandatory_corner_edges(N, v2e)

    x = [Bool(f"x_{i}") for i in range(E)]
    s = Solver()
    for vtx, es in v2e.items():
        s.add(PbEq([(x[e], 1) for e in es], 2))
    for e in mand:
        s.add(x[e])

    t0 = time.perf_counter()
    tours = []
    tentativas = 0
    while len(tours) < K:
        if time.perf_counter() - t0 > timeout:
            break
        tentativas += 1
        if s.check() != sat:
            break
        m = s.model()
        ativos = [e for e in range(E) if m[x[e]] is not None and bool(m[x[e]])]
        adj = [[] for _ in range(V)]
        for e in ativos:
            u, v = edges_uv[e]
            adj[u].append(v)
            adj[v].append(u)
        bfs = [0]
        seen = {0}
        qi = 0
        while qi < len(bfs):
            node = bfs[qi]; qi += 1
            for nb in adj[node]:
                if nb not in seen:
                    seen.add(nb)
                    bfs.append(nb)
        if len(seen) == V:
            tours.append(ativos)
            s.add(Sum([x[e] if e in ativos else Not(x[e]) for e in range(E)]) < E)
        else:
            ac = [e for e in range(E) if edges_uv[e][0] in seen and edges_uv[e][1] in seen
                  and m[x[e]] is not None and bool(m[x[e]])]
            if ac:
                s.add(Or([Not(x[e]) for e in ac]))
    t_total = time.perf_counter() - t0
    return {"tours": len(tours), "t_total": t_total, "tentativas": tentativas}


def run_bt_v1(K, timeout):
    """BT v1: propagação R2+R3+R6, pressão de vértice, SEM Union-Find."""
    from propagation_engine_10x10 import carregar_dados
    from backtracking_10x10 import Backtracker10
    dados = carregar_dados(strict_pair_mode="todos")
    bt = Backtracker10(dados, alvo_tours=K, timeout_s=timeout)
    r = bt.executar()
    return {"tours": r["n_tours"], "t_total": r["tempo_s"], "nos": r["n_nos"],
            "nos_por_tour": r["nos_por_tour"]}


def run_bt_v2(K, timeout):
    """BT v2: R2 + Union-Find incremental (MinimalV2 de scaling_minimal_v2)."""
    from scaling_minimal_v2 import MinimalV2
    bt = MinimalV2(n=N, alvo=K, timeout=timeout)
    r = bt.executar()
    return {"tours": r["n_tours"], "t_total": r["tempo_s"], "nos": r["n_nos"],
            "nos_por_tour": r["nos_por_tour"]}


def run_bt_theory(K, timeout):
    """BT theory: v2 + heurística f∞ teórica."""
    from theory_heuristic import TheoryHeuristicV2
    bt = TheoryHeuristicV2(n=N, alvo=K, timeout=timeout)
    t0 = time.perf_counter()
    r = bt.executar()
    dt = time.perf_counter() - t0
    return {"tours": r["n_tours"], "t_total": dt, "nos": r["n_nos"],
            "nos_por_tour": r["nos_por_tour"]}


def run_with_reps(func, K, timeout, reps=REPS, label=""):
    times = []
    last_result = None
    for i in range(reps):
        np.random.seed(SEED + i)
        random.seed(SEED + i)
        print(f"    rep {i+1}/{reps}...", end=" ", flush=True)
        r = func(K, timeout)
        print(f"t={r['t_total']:.3f}s  tours={r['tours']}")
        times.append(r['t_total'])
        last_result = r
        if r['t_total'] > timeout * 0.99:
            print(f"    TIMEOUT — interrompendo reps")
            break
    med = median(times)
    last_result['t_total_median'] = med
    last_result['all_times'] = times
    return last_result


def main():
    print(f"=== Benchmark Tabela 9: 10×10, K={K}, seed={SEED}, reps={REPS} ===\n")
    results = {}

    # BT theory (fastest first for quick feedback)
    print("[5] BT theory (f∞ + Union-Find):")
    results['BT_theory'] = run_with_reps(run_bt_theory, K, TIMEOUT_BT, label="BT theory")

    print(f"\n[4] BT v2 (Union-Find incremental):")
    results['BT_v2'] = run_with_reps(run_bt_v2, K, TIMEOUT_BT, label="BT v2")

    print(f"\n[3] BT v1 (pressão de vértice + R2/R3/R6):")
    results['BT_v1'] = run_with_reps(run_bt_v1, K, TIMEOUT_BT, label="BT v1")

    print(f"\n[1] Z3 puro:")
    results['Z3_puro'] = run_with_reps(run_z3_puro, K, TIMEOUT_Z3, reps=1, label="Z3 puro")

    print(f"\n[2] Z3 + 8 obrigatórias:")
    results['Z3_mandatory'] = run_with_reps(run_z3_mandatory, K, TIMEOUT_Z3, reps=1, label="Z3 mandatory")

    # Summary
    print("\n" + "="*70)
    print(f"{'Método':<25} {'t_total (s)':>12} {'speedup':>10} {'tours':>6}")
    print("-"*70)

    z3_t = results['Z3_puro']['t_total_median']
    paper_values = {
        'Z3_puro': 33.11,
        'Z3_mandatory': 54.42,
        'BT_v1': 7.18,
        'BT_v2': 1.31,
        'BT_theory': 0.67,
    }
    paper_speedups = {
        'Z3_puro': 1.0,
        'Z3_mandatory': 0.61,
        'BT_v1': 4.6,
        'BT_v2': 25.0,
        'BT_theory': 48.0,
    }

    for key, label in [('Z3_puro', 'Z3 puro'), ('Z3_mandatory', 'Z3 + 8 obrigatórias'),
                       ('BT_v1', 'BT v1 (pressão+R2)'), ('BT_v2', 'BT v2 (+Union-Find)'),
                       ('BT_theory', 'BT theory (+f∞)')]:
        r = results[key]
        t = r['t_total_median']
        sp = z3_t / t if t > 0 else float('inf')
        tours = r['tours']

        # compare with paper
        paper_t = paper_values[key]
        ratio = t / paper_t
        if ratio < 0.5 or ratio > 2.0:
            flag = "❌"
        elif abs(ratio - 1.0) > 0.2:
            flag = "⚠️"
        else:
            flag = "✅"

        print(f"{label:<25} {t:>10.3f}s  {sp:>8.1f}×  {tours:>5d}  "
              f"| paper={paper_t:.2f}s  ratio={ratio:.2f}  {flag}")

    # Save
    out = ROOT / "data" / "benchmark_tabela9_rerun.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(results, indent=2, default=str))
    print(f"\nSalvo em {out}")


if __name__ == "__main__":
    main()
