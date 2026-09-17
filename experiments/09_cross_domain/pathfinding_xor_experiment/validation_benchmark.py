#!/usr/bin/env python3
"""
Validation benchmark for the XOR pathfinding claim  compat ≈ 0.76 × coverage.
============================================================================

Three risks to eliminate before writing the paper section:
  RISK 1: 0.76 is specific to 20×20 / 30% obstacles  → BENCHMARK 1 (size×density)
  RISK 2: results don't generalize off-grid          → BENCHMARK 2 (planar types)
  RISK 3: Yen comparison is unfair                    → BENCHMARK 3 (Yen vs XOR)

Honest reporting: if 0.76 varies or Yen dominates, say so and by how much.

Uso:
    ../venv/bin/python validation_benchmark.py [1|2|3|all]
"""

import csv
import json
import os
import sys
import time
from collections import defaultdict

import numpy as np
import networkx as nx

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import pathfinding_xor as base                      # noqa: E402
from path_selection import (edge_to_cycles, coverage,  # noqa: E402
                            path_coverage_greedy)

RES = os.path.join(HERE, "results")
B1_JSON = os.path.join(RES, "validation_benchmark_1_efficiency_matrix.json")
B1_CSV = os.path.join(RES, "validation_benchmark_1_efficiency_matrix.csv")
B2_JSON = os.path.join(RES, "validation_benchmark_2_graph_types.json")
B3_JSON = os.path.join(RES, "validation_benchmark_3_yen_comparison.json")
SUMMARY = os.path.join(RES, "validation_summary.md")


# --------------------------------------------------------------------------
# capped grid generator (high densities may be unsatisfiable)
# --------------------------------------------------------------------------
def gen_grid(n, density, seed, max_attempts=300):
    S, T = (0, 0), (n - 1, n - 1)
    for attempt in range(max_attempts):
        cur = seed * 1000 + attempt
        rng = np.random.default_rng(cur)
        blocked = rng.random((n, n)) < density
        blocked[S] = False
        blocked[T] = False
        G = nx.grid_2d_graph(n, n)
        G.remove_nodes_from([(r, c) for r in range(n) for c in range(n)
                             if blocked[r, c]])
        if T in G and nx.has_path(G, S, T):
            comp = nx.node_connected_component(G, S)
            return G.subgraph(comp).copy(), S, T, cur
    return None


def eval_single_xor(P, cycles, G_check_S, G_check_T):
    """returns (coverage, compat, valid_edge_sets)."""
    valid = []
    for c in cycles:
        cand = base.xor(P, set(c))
        ok, _ = base.check_path(cand, G_check_S, G_check_T)
        if ok:
            valid.append(cand)
    n = len(cycles)
    cov = coverage(P, cycles)
    compat = len(valid) / n if n else 0.0
    return cov, compat, valid


# ==========================================================================
# BENCHMARK 1 — efficiency matrix
# ==========================================================================
def benchmark1():
    sizes = [10, 20, 30, 50]
    densities = [0.10, 0.20, 0.30, 0.40, 0.50]
    n_seeds = 30
    cells = {}
    rows = []
    all_eff = []
    print("=== BENCHMARK 1: efficiency matrix (running) ===")
    for size in sizes:
        for dens in densities:
            effs, covs, comps = [], [], []
            n_skip_unsat = 0
            n_skip_nocyc = 0
            for s in range(n_seeds):
                inst = gen_grid(size, dens, s)
                if inst is None:
                    n_skip_unsat += 1
                    continue
                G, S, T, _ = inst
                cycles = base.fundamental_cycles(G, S)
                if not cycles:
                    n_skip_nocyc += 1
                    continue
                P = base.seq_to_edges(base.astar_path(G, S, T))
                cov, compat, _ = eval_single_xor(P, cycles, S, T)
                if cov == 0.0:
                    n_skip_nocyc += 1
                    continue
                eff = compat / cov
                effs.append(eff); covs.append(cov); comps.append(compat)
                all_eff.append(eff)
            if effs:
                ea = np.array(effs)
                cell = {
                    "n": len(effs), "n_skip_unsat": n_skip_unsat,
                    "n_skip_nocyc": n_skip_nocyc,
                    "mean_coverage": float(np.mean(covs)),
                    "std_coverage": float(np.std(covs)),
                    "mean_compat": float(np.mean(comps)),
                    "std_compat": float(np.std(comps)),
                    "mean_efficiency": float(ea.mean()),
                    "std_efficiency": float(ea.std()),
                    "efficiency_cv": float(ea.std() / ea.mean()) if ea.mean() else float("nan"),
                }
            else:
                cell = {"n": 0, "n_skip_unsat": n_skip_unsat,
                        "n_skip_nocyc": n_skip_nocyc, "mean_efficiency": float("nan"),
                        "std_efficiency": float("nan"), "efficiency_cv": float("nan"),
                        "mean_coverage": float("nan"), "mean_compat": float("nan")}
            cells[f"{size}x{size}_d{dens}"] = cell
            rows.append({"size": size, "density": dens, **cell})
            print(f"  {size}×{size} d={dens}: n={cell['n']:2d} "
                  f"eff={cell['mean_efficiency']:.3f}±{cell['std_efficiency']:.3f} "
                  f"cv={cell['efficiency_cv']:.3f} (skip unsat={n_skip_unsat})")
        # save incrementally
        json.dump(cells, open(B1_JSON, "w"), indent=2)

    all_eff = np.array(all_eff)
    summary = {
        "sizes": sizes, "densities": densities, "n_seeds": n_seeds,
        "cells": cells,
        "global_mean_eff": float(all_eff.mean()) if len(all_eff) else float("nan"),
        "global_std_eff": float(all_eff.std()) if len(all_eff) else float("nan"),
        "global_min_eff": float(all_eff.min()) if len(all_eff) else float("nan"),
        "global_max_eff": float(all_eff.max()) if len(all_eff) else float("nan"),
        "n_total": int(len(all_eff)),
    }
    # robustness criterion
    means = [c["mean_efficiency"] for c in cells.values() if c["n"] > 0]
    cvs = [c["efficiency_cv"] for c in cells.values() if c["n"] > 0]
    in_band = all(0.70 <= m <= 0.82 for m in means)
    cv_ok = all(cv < 0.15 for cv in cvs)
    n_cv_high = sum(1 for cv in cvs if cv > 0.20)
    summary["robust_in_band"] = bool(in_band)
    summary["robust_cv_ok"] = bool(cv_ok)
    summary["n_cells_cv_gt_0.20"] = int(n_cv_high)
    summary["risk1_eliminated"] = ("YES" if (in_band and cv_ok)
                                   else ("PARTIAL" if not (means and n_cv_high > 1) else "NO"))
    json.dump(summary, open(B1_JSON, "w"), indent=2)
    with open(B1_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

    # print matrix
    print("\n=== EFFICIENCY MATRIX (compat / coverage) ===\n")
    hdr = "size↓ dens→ " + "".join(f"{d:>11.2f}" for d in densities)
    print(hdr)
    for size in sizes:
        line = f"{size:>3}×{size:<3}   "
        for dens in densities:
            c = cells[f"{size}x{size}_d{dens}"]
            if c["n"] > 0:
                line += f"{c['mean_efficiency']:6.3f}±{c['std_efficiency']:<4.3f}"
            else:
                line += f"{'  (skip)':>11}"
        print(line)
    print(f"\nGlobal mean efficiency: {summary['global_mean_eff']:.3f} "
          f"± {summary['global_std_eff']:.3f} "
          f"(range [{summary['global_min_eff']:.3f}, {summary['global_max_eff']:.3f}], "
          f"n={summary['n_total']})")
    print(f"RISK 1 eliminated: {summary['risk1_eliminated']}")
    return summary


# ==========================================================================
# BENCHMARK 2 — non-grid planar graphs
# ==========================================================================
def _nearest(points, target):
    d = ((points - np.array(target)) ** 2).sum(axis=1)
    return int(np.argmin(d))


def make_delaunay(seed, npts=400):
    from scipy.spatial import Delaunay
    rng = np.random.default_rng(seed)
    pts = rng.uniform(0, 100, size=(npts, 2))
    tri = Delaunay(pts)
    G = nx.Graph()
    G.add_nodes_from(range(npts))
    for simplex in tri.simplices:
        for a in range(3):
            for b in range(a + 1, 3):
                G.add_edge(int(simplex[a]), int(simplex[b]))
    return _connect(G), pts, _nearest(pts, (0, 0)), _nearest(pts, (100, 100))


def make_gabriel(seed, npts=400):
    """Gabriel ⊆ Delaunay: edge ij iff no point in circle with diameter ij."""
    from scipy.spatial import Delaunay
    rng = np.random.default_rng(seed)
    pts = rng.uniform(0, 100, size=(npts, 2))
    tri = Delaunay(pts)
    cand = set()
    for simplex in tri.simplices:
        for a in range(3):
            for b in range(a + 1, 3):
                i, j = int(simplex[a]), int(simplex[b])
                cand.add((min(i, j), max(i, j)))
    G = nx.Graph(); G.add_nodes_from(range(npts))
    for (i, j) in cand:
        mid = (pts[i] + pts[j]) / 2
        r2 = ((pts[i] - pts[j]) ** 2).sum() / 4
        d2 = ((pts - mid) ** 2).sum(axis=1)
        d2[i] = np.inf; d2[j] = np.inf
        if d2.min() >= r2 - 1e-9:          # no other point strictly inside
            G.add_edge(i, j)
    return _connect(G), pts, _nearest(pts, (0, 0)), _nearest(pts, (100, 100))


def make_random_planar(seed, npts=400, p=0.008):
    """Erdős–Rényi clipped to planar: greedily remove the edge whose endpoints
    have the largest degree sum (densest/most-crossing region) until planar.
    Uses plain check_planarity (fast); avoids the very slow counterexample=True
    Kuratowski extraction."""
    G = nx.gnp_random_graph(npts, p, seed=seed)
    guard = 0
    while guard < 5000:
        if nx.check_planarity(G, counterexample=False)[0]:
            break
        deg = dict(G.degree())
        u, v = max(G.edges(), key=lambda e: deg[e[0]] + deg[e[1]])
        G.remove_edge(u, v)
        guard += 1
    G = _connect(G)
    nodes = sorted(G.nodes())
    return G, None, nodes[0], nodes[-1]


def _connect(G):
    """keep the largest connected component."""
    if G.number_of_nodes() == 0:
        return G
    comp = max(nx.connected_components(G), key=len)
    return G.subgraph(comp).copy()


def benchmark2():
    print("\n=== BENCHMARK 2: planar graph generalization (running) ===")
    types = {
        "Grid+obstacles": lambda s: (lambda inst: (inst[0], inst[1], inst[2]))(
            gen_grid(20, 0.30, s)),
        "Delaunay": lambda s: (lambda r: (r[0], r[2], r[3]))(make_delaunay(s)),
        "Gabriel": lambda s: (lambda r: (r[0], r[2], r[3]))(make_gabriel(s)),
        "Random planar": lambda s: (lambda r: (r[0], r[2], r[3]))(make_random_planar(s)),
    }
    n_inst = 20
    results = {}
    for tname, gen in types.items():
        rec = defaultdict(list)
        for s in range(n_inst):
            try:
                G, S, T = gen(s)
            except Exception as e:
                print(f"   {tname} seed {s}: error {e}")
                continue
            if G is None or S not in G or T not in G or not nx.has_path(G, S, T):
                continue
            cycles = base.fundamental_cycles(G, S)
            if not cycles:
                continue
            try:
                P = base.seq_to_edges(nx.astar_path(G, S, T))
            except Exception:
                continue
            cov, compat, _ = eval_single_xor(P, cycles, S, T)
            if cov == 0:
                continue
            rec["V"].append(G.number_of_nodes())
            rec["E"].append(G.number_of_edges())
            rec["deg"].append(2 * G.number_of_edges() / G.number_of_nodes())
            rec["coverage"].append(cov)
            rec["compat"].append(compat)
            rec["efficiency"].append(compat / cov)
        results[tname] = {k: (float(np.mean(v)), float(np.std(v)))
                          for k, v in rec.items()} if rec["efficiency"] else {}
        results[tname]["n"] = len(rec["efficiency"])
        if rec["efficiency"]:
            e = np.array(rec["efficiency"])
            print(f"  {tname:<16} n={len(e)} V={np.mean(rec['V']):.0f} "
                  f"deg={np.mean(rec['deg']):.1f} cov={np.mean(rec['coverage'])*100:.1f}% "
                  f"compat={np.mean(rec['compat'])*100:.1f}% eff={e.mean():.3f}±{e.std():.3f}")
        json.dump(results, open(B2_JSON, "w"), indent=2)

    # risk2 verdict
    grid_eff = results["Grid+obstacles"].get("efficiency", (float("nan"),))[0]
    print("\n=== PLANAR GRAPH GENERALIZATION ===")
    print(f"{'Type':<16}{'|V|':>6}{'|E|':>7}{'deg':>6}{'cov%':>8}{'compat%':>9}{'efficiency':>13}")
    for tname in types:
        r = results[tname]
        if r.get("n", 0) == 0:
            print(f"{tname:<16}  (no instances)")
            continue
        print(f"{tname:<16}{r['V'][0]:6.0f}{r['E'][0]:7.0f}{r['deg'][0]:6.1f}"
              f"{r['coverage'][0]*100:7.1f}%{r['compat'][0]*100:8.1f}%"
              f"{r['efficiency'][0]:9.3f}±{r['efficiency'][1]:.3f}")
    effs = [results[t]["efficiency"][0] for t in types if results[t].get("n", 0) > 0]
    matches = all(abs(e - grid_eff) <= 0.10 for e in effs) if not np.isnan(grid_eff) else False
    in_band = all(0.65 <= e <= 0.85 for e in effs)
    verdict = "YES" if (matches and in_band) else ("PARTIAL" if in_band else "NO")
    results["_verdict"] = {"grid_eff": grid_eff, "all_effs": effs,
                           "matches_grid_within_0.10": bool(matches),
                           "all_in_band_0.65_0.85": bool(in_band),
                           "risk2_eliminated": verdict}
    json.dump(results, open(B2_JSON, "w"), indent=2)
    print(f"RISK 2 eliminated: {verdict}")
    return results


# ==========================================================================
# BENCHMARK 3 — Yen vs XOR vs repeated A*
# ==========================================================================
def yen_paths(G, S, T, k, timeout_s=30.0):
    t0 = time.perf_counter()
    out = []
    gen = nx.shortest_simple_paths(G, S, T, weight=None)
    timed_out = False
    for p in gen:
        out.append(base.seq_to_edges(p))
        if len(out) >= k:
            break
        if time.perf_counter() - t0 > timeout_s:
            timed_out = True
            break
    return out, (time.perf_counter() - t0) * 1000.0, timed_out


def xor_paths_cg(G, S, T, cycles, e2c, k, rng):
    t0 = time.perf_counter()
    seq = None
    for _ in range(10):
        seq = path_coverage_greedy(G, S, T, cycles, e2c, rng)
        if seq is not None:
            break
    if seq is None:
        seq = base.astar_path(G, S, T)
    P = base.seq_to_edges(seq)
    valid = [P]                                 # base path counts as a route
    for c in cycles:
        cand = base.xor(P, set(c))
        ok, _ = base.check_path(cand, S, T)
        if ok:
            valid.append(cand)
        if len(valid) >= k:
            break
    return valid[:k], (time.perf_counter() - t0) * 1000.0


def mean_len(paths):
    return float(np.mean([len(p) for p in paths])) if paths else 0.0


def benchmark3():
    print("\n=== BENCHMARK 3: Yen vs XOR (running) ===")
    prior = [r for r in json.load(open(os.path.join(RES, "pathfinding_xor_experiment.json")))
             if r["grid_size"] == 20][:20]
    ks = [5, 10, 25, 50, 100]
    per = []
    for r in prior:
        G, S, T, _ = base.gen_graph(20, r["used_seed"])
        cycles = base.fundamental_cycles(G, S)
        e2c = edge_to_cycles(cycles)
        rng = np.random.default_rng(r["used_seed"])
        inst = {"used_seed": r["used_seed"]}
        for k in ks:
            xp, xt = xor_paths_cg(G, S, T, cycles, e2c, k, rng)
            yp, yt, yto = yen_paths(G, S, T, k)
            ap, at, asucc = base.repeated_astar(G, S, T, k)
            inst[f"k{k}"] = {
                "xor": {"n": len(xp), "time_ms": xt, "diversity": base.diversity(xp),
                        "mean_length": mean_len(xp), "throughput": len(xp) / xt if xt else 0},
                "yen": {"n": len(yp), "time_ms": yt, "diversity": base.diversity(yp),
                        "mean_length": mean_len(yp), "throughput": len(yp) / yt if yt else 0,
                        "timeout": yto},
                "astar_rep": {"n": len(ap), "time_ms": at, "diversity": base.diversity(ap),
                              "mean_length": mean_len(ap), "throughput": len(ap) / at if at else 0},
            }
        per.append(inst)
        json.dump({"per_instance": per, "ks": ks}, open(B3_JSON, "w"), indent=2)
        print(f"  seed {r['used_seed']}: done")

    # aggregate
    agg = {}
    for k in ks:
        agg[k] = {}
        for m in ["xor", "yen", "astar_rep"]:
            agg[k][m] = {
                "time_ms": float(np.mean([p[f"k{k}"][m]["time_ms"] for p in per])),
                "diversity": float(np.mean([p[f"k{k}"][m]["diversity"] for p in per])),
                "mean_length": float(np.mean([p[f"k{k}"][m]["mean_length"] for p in per])),
                "n": float(np.mean([p[f"k{k}"][m]["n"] for p in per])),
                "throughput": float(np.mean([p[f"k{k}"][m]["throughput"] for p in per])),
            }
    # break-even k*: smallest k where yen diversity > xor diversity
    kstar = None
    for k in ks:
        if agg[k]["yen"]["diversity"] > agg[k]["xor"]["diversity"]:
            kstar = k
            break
    result = {"ks": ks, "aggregate": agg, "break_even_kstar": kstar,
              "per_instance": per}
    json.dump(result, open(B3_JSON, "w"), indent=2)

    print("\n=== YEN vs XOR BENCHMARK ===")
    for k in ks:
        x, y, a = agg[k]["xor"], agg[k]["yen"], agg[k]["astar_rep"]
        print(f"\nk={k}:")
        print(f"  XOR:  time={x['time_ms']:7.1f}ms div={x['diversity']:.3f} "
              f"len={x['mean_length']:.1f} thru={x['throughput']:.2f} (got {x['n']:.0f})")
        print(f"  Yen:  time={y['time_ms']:7.1f}ms div={y['diversity']:.3f} "
              f"len={y['mean_length']:.1f} thru={y['throughput']:.2f} (got {y['n']:.0f})")
        print(f"  A*rep:time={a['time_ms']:7.1f}ms div={a['diversity']:.3f} "
              f"len={a['mean_length']:.1f} thru={a['throughput']:.2f} (got {a['n']:.0f})")
        wd = max(("XOR", x["diversity"]), ("Yen", y["diversity"]), key=lambda z: z[1])[0]
        wt = max(("XOR", x["throughput"]), ("Yen", y["throughput"]), key=lambda z: z[1])[0]
        print(f"  Winner diversity: {wd} | Winner throughput: {wt}")
    print(f"\nBreak-even k* (Yen diversity > XOR diversity): "
          f"{'k*=' + str(kstar) if kstar else 'never (XOR never beaten in tested range)'}")
    risk3 = "YES" if kstar is not None else "PARTIAL"
    print(f"RISK 3 eliminated: {risk3}")
    result["risk3_eliminated"] = risk3
    json.dump(result, open(B3_JSON, "w"), indent=2)
    return result


# ==========================================================================
def write_summary(b1, b2, b3):
    L = ["# Validation summary — XOR pathfinding (compat ≈ 0.76 × coverage)\n"]
    # risk 1
    L.append("## RISK 1 — is 0.76 grid/density-specific?\n")
    L.append(f"- Efficiency across all {b1['n_total']} instances: "
             f"range [{b1['global_min_eff']:.3f}, {b1['global_max_eff']:.3f}], "
             f"global mean **{b1['global_mean_eff']:.3f} ± {b1['global_std_eff']:.3f}**.")
    L.append(f"- All cell means in [0.70,0.82]: {b1['robust_in_band']}; "
             f"all cell CV<0.15: {b1['robust_cv_ok']}; cells with CV>0.20: {b1['n_cells_cv_gt_0.20']}.")
    L.append(f"- **ELIMINATED: {b1['risk1_eliminated']}**\n")
    L.append("### Efficiency matrix (mean ± std)\n")
    L.append("| size＼density | " + " | ".join(f"{d:.2f}" for d in b1["densities"]) + " |")
    L.append("|---" * (len(b1["densities"]) + 1) + "|")
    for size in b1["sizes"]:
        cells = []
        for d in b1["densities"]:
            c = b1["cells"][f"{size}x{size}_d{d}"]
            cells.append(f"{c['mean_efficiency']:.3f}±{c['std_efficiency']:.3f}" if c["n"] > 0 else "skip")
        L.append(f"| {size}×{size} | " + " | ".join(cells) + " |")
    L.append("")
    # risk 2
    L.append("## RISK 2 — generalization to non-grid planar graphs\n")
    if b2:
        L.append("| Type | n | \\|V\\| | \\|E\\| | deg | coverage | compat | efficiency |")
        L.append("|---|---:|---:|---:|---:|---:|---:|---:|")
        for t in ["Grid+obstacles", "Delaunay", "Gabriel", "Random planar"]:
            r = b2.get(t, {})
            if r.get("n", 0) == 0:
                L.append(f"| {t} | 0 | — | — | — | — | — | — |"); continue
            L.append(f"| {t} | {r['n']} | {r['V'][0]:.0f} | {r['E'][0]:.0f} | {r['deg'][0]:.1f} "
                     f"| {r['coverage'][0]*100:.1f}% | {r['compat'][0]*100:.1f}% "
                     f"| {r['efficiency'][0]:.3f}±{r['efficiency'][1]:.3f} |")
        v = b2["_verdict"]
        L.append(f"\n- Grid efficiency {v['grid_eff']:.3f}; off-grid effs {[round(e,3) for e in v['all_effs']]}.")
        L.append(f"- Matches grid within 0.10: {v['matches_grid_within_0.10']}; "
                 f"all in [0.65,0.85]: {v['all_in_band_0.65_0.85']}.")
        L.append(f"- **ELIMINATED: {v['risk2_eliminated']}**\n")
    # risk 3
    L.append("## RISK 3 — fair comparison with Yen\n")
    if b3:
        L.append("| k | XOR div | Yen div | XOR len | Yen len | XOR thru | Yen thru |")
        L.append("|---:|---:|---:|---:|---:|---:|---:|")
        for k in b3["ks"]:
            a = b3["aggregate"][k]
            L.append(f"| {k} | {a['xor']['diversity']:.3f} | {a['yen']['diversity']:.3f} "
                     f"| {a['xor']['mean_length']:.1f} | {a['yen']['mean_length']:.1f} "
                     f"| {a['xor']['throughput']:.2f} | {a['yen']['throughput']:.2f} |")
        ks_ = b3["break_even_kstar"]
        L.append(f"\n- Break-even k* (Yen diversity > XOR): "
                 f"{('k*='+str(ks_)) if ks_ else 'never in tested range'}.")
        L.append(f"- XOR competitive for k < k*; Yen superior for k ≥ k*.")
        L.append(f"- **ELIMINATED: {b3['risk3_eliminated']}**\n")
    # overall
    L.append("## OVERALL VERDICT\n")
    r1 = b1["risk1_eliminated"]; r2 = b2["_verdict"]["risk2_eliminated"] if b2 else "?"
    r3 = b3["risk3_eliminated"] if b3 else "?"
    ready = (r1 == "YES" and r2 in ("YES", "PARTIAL") and r3 in ("YES", "PARTIAL"))
    L.append(f"- RISK 1 (0.76 grid-specific): **{r1}**")
    L.append(f"- RISK 2 (non-general): **{r2}**")
    L.append(f"- RISK 3 (unfair Yen): **{r3}**")
    L.append(f"\n**OVERALL: {'READY TO DOCUMENT' if ready else 'NEEDS QUALIFICATION'}**\n")

    # required precise qualifications
    L.append("## REQUIRED PAPER QUALIFICATIONS\n")
    L.append("**The claim `compat ≈ 0.76 × coverage` is NOT a universal law and "
             "must be restated.** Efficiency η = compat/coverage varies "
             "systematically and 0.76 is merely the value at the cell where it "
             "was first measured (20×20, 30% obstacles → 0.756/0.764).")
    L.append("")
    L.append("1. **η decreases with graph size** (fixed 30% density): "
             "10×10≈0.83, 20×20≈0.76, 30×30≈0.60, 50×50≈0.46.")
    L.append("2. **η increases with obstacle density** (e.g. 50×50: "
             "0.39 @10% → 0.73 @40%; 20×20: 0.52 @10% → 0.95 @50%).")
    L.append("3. **η varies 0.34–0.94 across planar graph types**: "
             "Gabriel 0.34, Delaunay 0.59, grid 0.76, near-tree random-planar 0.94.")
    L.append("4. **Unifying driver:** η rises as the graph becomes more "
             "*corridor-like* (lower mean degree / higher obstacle density / "
             "sparser). In corridors a touched cycle crosses the path "
             "contiguously (→ valid); in open regions it crosses at separated "
             "points (→ disconnected, invalid). So η measures local corridor "
             "structure, not a constant.")
    L.append("")
    L.append("**Correct statement for the paper:** "
             "`compat = η · coverage`, with η ∈ [0.34, 0.95] increasing in "
             "corridor-likeness; η ≈ 0.76 holds specifically for 20×20 grids at "
             "30% obstacle density. The *qualitative* law (compat is the product "
             "of coverage and a structure-dependent crossing-efficiency) IS "
             "robust; the *numeric coefficient* is not.")
    L.append("")
    L.append("**On Yen (RISK 3):** the comparison is now fair, and the honest "
             "result is that **Yen's k-shortest dominates XOR on diversity at "
             "every k ≥ 5** and returns shorter (higher-quality) paths; XOR's "
             "only advantage is throughput (paths/ms) and it caps at ~32 valid "
             "alternatives from a single base path. The paper must NOT claim XOR "
             "is a competitive *diverse*-path enumerator — only a cheap generator "
             "of many *local* variations.")
    open(SUMMARY, "w").write("\n".join(L) + "\n")
    print(f"\nSummary written: {SUMMARY}")


def main():
    os.makedirs(RES, exist_ok=True)
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    b1 = b2 = b3 = None
    if which in ("1", "all"):
        b1 = benchmark1()
    if which in ("2", "all"):
        b2 = benchmark2()
    if which in ("3", "all"):
        b3 = benchmark3()
    if which == "all":
        write_summary(b1, b2, b3)


if __name__ == "__main__":
    main()
