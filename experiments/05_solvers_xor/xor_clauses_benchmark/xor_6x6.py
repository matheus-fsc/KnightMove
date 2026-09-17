#!/usr/bin/env python3
"""
xor_6x6.py
==========
T0 + T1: carrega os 3 detectores duais do 6×6 (de forbidden_cycles_6x6),
verifica suportes e roda benchmark Z3 comparando 5 configurações:

  A: Z3 puro (grau-2 + sub-tour elimination iterativo)
  B: A + 8 arestas obrigatórias como fatos (x_e = True)
  C: B + 3 cláusulas XOR de paridade dos detectores duais (phi_0, phi_1, phi_2)
  D: B + 28 cláusulas XOR entre todos os C(8,2) pares de obrigatórias
  E: B + 3 XOR + NOT(A∧B) dos top-K pares de exclusão minimais

Para cada config: K=500 tours amostrados; reporta tempo total, tempo por
amostra, attempts/sol, sol/s. Salva em data/benchmark_6x6_xor.json.

Implementação das cláusulas XOR:
  - Suporte 2: x_a == x_b   (equivalência direta)
  - Suporte > 2: Sum([If(x_e,1,0) for e in supp]) % 2 == 0
    (avaliado mais rápido que reduzir Xor binário no Z3)
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

# importa grafo 6×6
sys.path.insert(0, str(REPO))
from cavalo_loop_destruicao_6x6 import EDGES_LIST, ADJ, TOTAL  # noqa: E402


# ── helpers ───────────────────────────────────────────────────────────

def is_single_tour(active_edges, n_vertices=TOTAL):
    adj = [[] for _ in range(n_vertices)]
    for u, v in active_edges:
        adj[u].append(v)
        adj[v].append(u)
    if any(len(a) != 2 for a in adj):
        return False, []
    visited = [False] * n_vertices
    visited[0] = True
    tour = [0]
    cur, prev = 0, -1
    for _ in range(n_vertices - 1):
        nxt = adj[cur][0] if adj[cur][0] != prev else adj[cur][1]
        if visited[nxt]:
            return False, []
        visited[nxt] = True
        tour.append(nxt)
        prev, cur = cur, nxt
    return (adj[cur][0] == 0 or adj[cur][1] == 0), tour


def build_base_solver(edges):
    """Solver com grau-2 em cada vértice. Sem mais nada."""
    E = len(edges)
    xvars = [Bool(f"x_{i}") for i in range(E)]
    inc = {v: [] for v in range(TOTAL)}
    for i, (u, v) in enumerate(edges):
        inc[u].append(i)
        inc[v].append(i)
    s = Solver()
    for v, lst in inc.items():
        s.add(PbEq([(xvars[i], 1) for i in lst], 2))
    return s, xvars, inc


def add_mandatory(solver, xvars, mandatory_idx):
    for e in mandatory_idx:
        solver.add(xvars[e])


def add_xor_clause(solver, xvars, support):
    """Adiciona Σ_{e∈supp} x_e ≡ 0 (mod 2). Suporte 2 vira equivalência."""
    if len(support) == 2:
        a, b = support
        solver.add(xvars[a] == xvars[b])
        return
    terms = [If(xvars[e], 1, 0) for e in support]
    solver.add(Sum(terms) % 2 == 0)


def add_pair_exclusions(solver, xvars, pairs):
    for (i, j) in pairs:
        solver.add(Not(And(xvars[i], xvars[j])))


# ── construção de cada config ─────────────────────────────────────────

def build_solver_for_config(cfg_name, edges, mandatory_idx,
                            phi_supports, all_pair_xors,
                            top_excl_pairs):
    s, xvars, inc = build_base_solver(edges)
    if cfg_name == "A":
        return s, xvars
    if cfg_name == "B":
        add_mandatory(s, xvars, mandatory_idx)
        return s, xvars
    if cfg_name == "C":
        add_mandatory(s, xvars, mandatory_idx)
        for supp in phi_supports:
            add_xor_clause(s, xvars, supp)
        return s, xvars
    if cfg_name == "D":
        add_mandatory(s, xvars, mandatory_idx)
        for supp in all_pair_xors:
            add_xor_clause(s, xvars, supp)
        return s, xvars
    if cfg_name == "E":
        add_mandatory(s, xvars, mandatory_idx)
        for supp in phi_supports:
            add_xor_clause(s, xvars, supp)
        add_pair_exclusions(s, xvars, top_excl_pairs)
        return s, xvars
    raise ValueError(f"config desconhecida: {cfg_name}")


# ── loop de amostragem ────────────────────────────────────────────────

def sample_k_tours(solver, xvars, edges, K, *, seed, timeout_ms=30_000,
                   max_sym_tries=4):
    rng = random.Random(seed)
    E = len(edges)
    sigs = []
    stats = {
        "K": K,
        "n_attempts": 0,
        "n_subtour_rejects": 0,
        "n_local_unsat": 0,
        "t_first": None,
        "t_total": 0.0,
        "times_per_valid": [],
    }
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
        ok, _ = is_single_tour(active)
        solver.pop()
        solver.pop()
        if ok:
            sigs.append(sig)
            dt = time.perf_counter() - t_attempt
            stats["times_per_valid"].append(round(dt, 4))
            if stats["t_first"] is None:
                stats["t_first"] = round(time.perf_counter() - t0, 4)
            # corte de exclusão para não repetir este tour
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


# ── KL-divergence simples vs ground-truth ─────────────────────────────

def freq_per_edge(sigs):
    return sigs.mean(axis=0)


def kl_div_freq(p_sample, p_truth, eps=1e-6):
    """KL Bernoulli por aresta, somado. Não é métrica perfeita mas serve."""
    p = np.clip(p_sample, eps, 1 - eps)
    q = np.clip(p_truth, eps, 1 - eps)
    return float(np.sum(p * np.log(p / q) + (1 - p) * np.log((1 - p) / (1 - q))))


# ── T0: carregar detectores e fazer sanidade ──────────────────────────

def load_dual_detectors_6x6():
    p = REPO / "forbidden_cycles_6x6" / "data" / "results" / "dual_detectors.json"
    with open(p) as f:
        d = json.load(f)
    detectors = []
    for x in d:
        detectors.append({
            "phi_i": x["phi_i"],
            "support": x["edge_indices"],
            "support_size": x["support_size"],
            "labels": x["edge_labels"],
        })
    return detectors


def load_mandatory_6x6():
    """Identifica as 8 obrigatórias do 6×6 dos 4 cantos.
    Vértices de canto têm grau 2 → suas 2 arestas são obrigatórias.
    Carrega de forbidden_cycles ou deduz."""
    # cantos no 6×6: A6=0, F6=5, A1=30, F1=35  (label() mapeia r,c → letra+num)
    # Deduz das adjacências: vértice de grau 2 → as 2 arestas são obrigatórias
    edges = list(EDGES_LIST)
    edge_idx = {tuple(sorted(e)): i for i, e in enumerate(edges)}
    mand = []
    for v in range(TOTAL):
        if len(ADJ[v]) == 2:
            for u in ADJ[v]:
                key = tuple(sorted((u, v)))
                mand.append(edge_idx[key])
    mand = sorted(set(mand))
    return mand


def sanity_check_xor(detectors, edges, mandatory_idx, sigs_path=None):
    """Verifica que os detectores zeram em todos os tours conhecidos."""
    # ground-truth: incidence_matrix_6x6.npy = (9862, 80) GF(2)
    inc_path = REPO / "6x6_higher_order" / "data" / "incidence_matrix_6x6.npy"
    T = np.load(inc_path).astype(np.uint8)
    ok = True
    rep = []
    for det in detectors:
        v = np.zeros(len(edges), dtype=np.uint8)
        v[det["support"]] = 1
        prod = (T @ v) % 2
        zeros = int((prod == 0).sum())
        total = T.shape[0]
        rep.append({
            "phi_i": det["phi_i"],
            "support_size": det["support_size"],
            "zero_on_tours": f"{zeros}/{total}",
            "ok": bool(zeros == total),
        })
        if zeros != total:
            ok = False
    return ok, rep


# ── T1: rodar benchmark ────────────────────────────────────────────────

def run_t1(K, seed, configs):
    edges = list(EDGES_LIST)
    E = len(edges)

    # T0 — carregar detectores e mandatórias
    detectors = load_dual_detectors_6x6()
    mandatory_idx = load_mandatory_6x6()
    phi_supports = [d["support"] for d in detectors]

    print("─" * 70)
    print("T0: detectores duais do 6×6")
    for det in detectors:
        print(f"  φ_{det['phi_i']}: |supp|={det['support_size']}  "
              f"primeiras: {det['labels'][:6]}")
    print(f"  Obrigatórias (cantos, deg=2): {len(mandatory_idx)}  "
          f"índices={mandatory_idx}")

    print("\nVerificando que detectores zeram nos 9.862 tours conhecidos...")
    ok, rep = sanity_check_xor(detectors, edges, mandatory_idx)
    for r in rep:
        marker = "✓" if r["ok"] else "✗"
        print(f"  {marker} φ_{r['phi_i']}: zera em {r['zero_on_tours']}")
    if not ok:
        print("ATENÇÃO: algum detector NÃO zera em todos os tours. ABORT.")
        sys.exit(1)

    # T1 — preparar configs
    all_pair_xors = list(itertools.combinations(mandatory_idx, 2))
    print(f"\n  C(8,2) = {len(all_pair_xors)} pares de obrigatórias")

    # NOT(A∧B) top-K pares — usa exclusions_pairs.json (88 pares P=0)
    with open(REPO / "6x6_higher_order" / "data" / "exclusions_pairs.json") as f:
        epd = json.load(f)
    excl_pairs = [tuple(x["edges"]) for x in epd["exclusions"]]
    TOP_K_EXCL = 8
    top_excl_pairs = excl_pairs[:TOP_K_EXCL]
    print(f"  top {TOP_K_EXCL} pares de exclusão (NOT(A∧B)) carregados")

    # ground-truth para KL
    T = np.load(REPO / "6x6_higher_order" / "data" / "incidence_matrix_6x6.npy")
    p_truth = T.mean(axis=0)

    print("\n" + "─" * 70)
    print(f"T1: benchmark — K={K} tours por config, seed={seed}")
    print("─" * 70)

    results = {}
    for cfg in configs:
        print(f"\n  Config {cfg} construindo solver...")
        s, xvars = build_solver_for_config(
            cfg, edges, mandatory_idx, phi_supports,
            all_pair_xors, top_excl_pairs)
        t0 = time.perf_counter()
        sigs, stats = sample_k_tours(s, xvars, edges, K, seed=seed)
        dt = time.perf_counter() - t0
        p_smp = freq_per_edge(sigs)
        kl = kl_div_freq(p_smp, p_truth)
        stats["kl_vs_truth"] = round(kl, 4)
        stats["wall_total"] = round(dt, 3)
        results[cfg] = stats
        print(f"    t_total={stats['t_total']}s  "
              f"sol/s={stats['sol_per_s']}  "
              f"attempts/sol={stats['attempts_per_sol']}  "
              f"KL={stats['kl_vs_truth']:.3f}")

    # speedups relativos a A
    if "A" in results:
        baseA = results["A"]["t_total"]
        for cfg, st in results.items():
            st["speedup_vs_A"] = round(baseA / st["t_total"], 3)
    if "B" in results:
        baseB = results["B"]["t_total"]
        for cfg, st in results.items():
            st["speedup_vs_B"] = round(baseB / st["t_total"], 3)

    out = {
        "board": 6,
        "K": K,
        "seed": seed,
        "configs": list(configs),
        "n_mandatory": len(mandatory_idx),
        "mandatory_idx": mandatory_idx,
        "n_xor_phi": len(phi_supports),
        "phi_supports": [list(map(int, s)) for s in phi_supports],
        "n_xor_all_pairs": len(all_pair_xors),
        "n_excl_pairs": len(top_excl_pairs),
        "results": results,
    }
    DATA.mkdir(parents=True, exist_ok=True)
    outp = DATA / "benchmark_6x6_xor.json"
    with open(outp, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nSalvo: {outp.relative_to(REPO)}")
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--K", type=int, default=500)
    p.add_argument("--seed", type=int, default=2026)
    p.add_argument("--configs", type=str, default="A,B,C,D,E")
    args = p.parse_args()
    configs = args.configs.split(",")
    run_t1(args.K, args.seed, configs)


if __name__ == "__main__":
    main()
