#!/usr/bin/env python3
"""
z3_verify.py
============
Verifica via Z3 que as exclusões minimais encontradas (pares/triplas/quadras
com P(∧)=0 nas 9.862 soluções exaustivas) são *formalmente* UNSAT no
problema de tour hamiltoniano fechado do cavalo 6×6.

Formulação:
  - x_e ∈ {0,1} por aresta
  - grau-2 em cada vértice
  - subtour elimination iterativo (corte para cada componente desconexa
    encontrada)
  - Para a exclusão (e_1,...,e_k): adicionar x_{e_i} = 1 ∀ i e checar
    UNSAT após exaurir sub-tours

Como temos ground truth exaustivo, P=0 ⇒ não existe tour Hamiltoniano
usando essas arestas simultaneamente. Z3 *deve* confirmar UNSAT.
Se algum for SAT, indica bug na construção de T.

Por padrão, amostra n_sample por ordem (sample uniforme); use
--full pra checar todas (lento para quadras: 17k × ~poucos seg).

Saída:
  data/z3_verification.json
"""

import argparse
import json
import random
import sys
import time
from pathlib import Path

import numpy as np
from z3 import Bool, Or, PbEq, Solver, sat, unsat, is_true

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"

sys.path.insert(0, str(ROOT.parent))
from cavalo_loop_destruicao_6x6 import EDGES_LIST, ADJ, TOTAL  # noqa: E402


def is_single_tour_local(active_edges, n_vertices=TOTAL):
    """BFS para detectar tour único vs união de sub-tours."""
    adj = [[] for _ in range(n_vertices)]
    for u, v in active_edges:
        adj[u].append(v)
        adj[v].append(u)
    if any(len(a) != 2 for a in adj):
        return False, set()
    visited = [False] * n_vertices
    visited[0] = True
    cur, prev = 0, -1
    comp = {0}
    for _ in range(n_vertices - 1):
        nxt = adj[cur][0] if adj[cur][0] != prev else adj[cur][1]
        if visited[nxt]:
            return False, comp
        visited[nxt] = True
        comp.add(nxt)
        prev, cur = cur, nxt
    return (0 in adj[cur]), comp


def verify_unsat(forced_edges, *, timeout_ms=30000, max_iters=50):
    """
    Tenta provar que forçar todas as forced_edges = 1 é UNSAT para
    um tour Hamiltoniano fechado (grau-2 + conectividade).

    Retorna: ("mandatory_unsat", iters)  se UNSAT confirmado
              ("sat_counterexample", iters) se SAT (problema!)
              ("inconclusive", iters)  se timeout/max_iters
    """
    edges = list(EDGES_LIST)
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

    for e in forced_edges:
        s.add(xvars[e] == True)

    for it in range(max_iters):
        res = s.check()
        if res == unsat:
            return "unsat", it + 1
        if res != sat:
            return "inconclusive", it + 1
        m = s.model()
        active = []
        bits = np.zeros(E, dtype=np.uint8)
        for i, (u, v) in enumerate(edges):
            if is_true(m.evaluate(xvars[i])):
                active.append((u, v))
                bits[i] = 1
        ok, _ = is_single_tour_local(active)
        if ok:
            return "sat_counterexample", it + 1
        # corte: essa configuração não pode reaparecer
        s.add(Or([xvars[i] != bool(bits[i]) for i in range(E)]))

    return "inconclusive", max_iters


def sample_indices(n, k, rng):
    if n <= k:
        return list(range(n))
    return sorted(rng.sample(range(n), k))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--n-sample-pairs", type=int, default=None,
                   help="default: TODAS (88 pares)")
    p.add_argument("--n-sample-triples", type=int, default=50)
    p.add_argument("--n-sample-quads", type=int, default=50)
    p.add_argument("--seed", type=int, default=2026)
    p.add_argument("--timeout-ms", type=int, default=30000)
    p.add_argument("--full", action="store_true",
                   help="verifica TODAS as exclusões (caro)")
    args = p.parse_args()

    rng = random.Random(args.seed)
    results = {"pairs": [], "triples": [], "quads": []}
    totals = {"pairs": {"unsat": 0, "sat": 0, "inconclusive": 0},
              "triples": {"unsat": 0, "sat": 0, "inconclusive": 0},
              "quads": {"unsat": 0, "sat": 0, "inconclusive": 0}}

    # ── pares ────────────────────────────────────────────────────────
    pairs_data = json.loads((DATA / "exclusions_pairs.json").read_text())
    pairs = pairs_data["exclusions"]
    n_pairs = len(pairs)
    pair_indices = (list(range(n_pairs)) if args.full or args.n_sample_pairs is None
                    else sample_indices(n_pairs, args.n_sample_pairs, rng))
    print(f"Verificando pares: {len(pair_indices)} / {n_pairs}")
    t0 = time.perf_counter()
    for k, idx in enumerate(pair_indices):
        forced = pairs[idx]["edges"]
        verdict, iters = verify_unsat(forced, timeout_ms=args.timeout_ms)
        entry = {"index": idx, "edges": forced,
                 "names": pairs[idx]["edge_names"],
                 "verdict": verdict, "iters": iters}
        results["pairs"].append(entry)
        if verdict == "unsat":
            totals["pairs"]["unsat"] += 1
        elif verdict == "sat_counterexample":
            totals["pairs"]["sat"] += 1
        else:
            totals["pairs"]["inconclusive"] += 1
        if (k + 1) % 10 == 0 or k + 1 == len(pair_indices):
            elapsed = time.perf_counter() - t0
            print(f"  [{k+1}/{len(pair_indices)}] {verdict:>20s} "
                  f"({iters} iters) elapsed={elapsed:.1f}s")

    # ── triplas ──────────────────────────────────────────────────────
    triples_data = json.loads((DATA / "exclusions_triples.json").read_text())
    triples = triples_data["minimal_triples"]
    n_tr = len(triples)
    tr_indices = (list(range(n_tr)) if args.full
                  else sample_indices(n_tr, args.n_sample_triples, rng))
    print(f"\nVerificando triplas: {len(tr_indices)} / {n_tr}")
    t0 = time.perf_counter()
    for k, idx in enumerate(tr_indices):
        forced = triples[idx]["edges"]
        verdict, iters = verify_unsat(forced, timeout_ms=args.timeout_ms)
        entry = {"index": idx, "edges": forced,
                 "names": triples[idx]["edge_names"],
                 "verdict": verdict, "iters": iters}
        results["triples"].append(entry)
        if verdict == "unsat":
            totals["triples"]["unsat"] += 1
        elif verdict == "sat_counterexample":
            totals["triples"]["sat"] += 1
        else:
            totals["triples"]["inconclusive"] += 1
        if (k + 1) % 10 == 0 or k + 1 == len(tr_indices):
            elapsed = time.perf_counter() - t0
            print(f"  [{k+1}/{len(tr_indices)}] {verdict:>20s} "
                  f"({iters} iters) elapsed={elapsed:.1f}s")

    # ── quadras ──────────────────────────────────────────────────────
    quads_data = json.loads((DATA / "exclusions_quads.json").read_text())
    quads = quads_data["minimal_quads"]
    n_q = len(quads)
    q_indices = (list(range(n_q)) if args.full
                 else sample_indices(n_q, args.n_sample_quads, rng))
    print(f"\nVerificando quadras: {len(q_indices)} / {n_q}")
    t0 = time.perf_counter()
    for k, idx in enumerate(q_indices):
        forced = quads[idx]["edges"]
        verdict, iters = verify_unsat(forced, timeout_ms=args.timeout_ms)
        entry = {"index": idx, "edges": forced,
                 "names": quads[idx]["edge_names"],
                 "verdict": verdict, "iters": iters}
        results["quads"].append(entry)
        if verdict == "unsat":
            totals["quads"]["unsat"] += 1
        elif verdict == "sat_counterexample":
            totals["quads"]["sat"] += 1
        else:
            totals["quads"]["inconclusive"] += 1
        if (k + 1) % 10 == 0 or k + 1 == len(q_indices):
            elapsed = time.perf_counter() - t0
            print(f"  [{k+1}/{len(q_indices)}] {verdict:>20s} "
                  f"({iters} iters) elapsed={elapsed:.1f}s")

    # ── sumário ──────────────────────────────────────────────────────
    print(f"\n{'─' * 60}")
    print("VEREDICTOS")
    print("─" * 60)
    for order, t in totals.items():
        total = sum(t.values())
        print(f"  {order:>8s}: total={total}  unsat={t['unsat']}  "
              f"sat={t['sat']}  inconclusive={t['inconclusive']}")

    output = {
        "args": vars(args),
        "totals": totals,
        "results": results,
    }
    with open(DATA / "z3_verification.json", "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nSalvo: data/z3_verification.json")

    # alerta se SAT
    any_sat = any(t["sat"] > 0 for t in totals.values())
    if any_sat:
        print("\n*** ALERTA: ao menos uma exclusão retornou SAT — "
              "investigar construção de T ***")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
