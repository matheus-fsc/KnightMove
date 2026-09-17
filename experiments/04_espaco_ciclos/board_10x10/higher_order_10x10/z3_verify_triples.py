#!/usr/bin/env python3
"""
z3_verify_triples.py
====================
Prova formal das triplas candidatas via Z3.

Para cada candidata, força as 3 arestas a 1 e tenta provar UNSAT via
subtour elimination iterativo (igual ao usado em mandatory_edges.py).

Classifica:
  proven      → UNSAT confirmado
  false_pos   → SAT (encontrou tour com as 3 arestas)
  unverified  → timeout / max_iters não convergiram

Saída:
  data/triples_verified.json
"""

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
from z3 import Bool, Or, PbEq, Solver, sat, unsat, is_true

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
PARENT = ROOT.parent

sys.path.insert(0, str(PARENT))
from graph_10x10 import build_graph, label, TOTAL  # noqa: E402
from sample_tours import is_single_tour  # noqa: E402

ADJ, EDGES = build_graph()


def verify_unsat(forced_edges, *, timeout_ms=30000, max_iters=100):
    E = len(EDGES)
    xvars = [Bool(f"x_{i}") for i in range(E)]
    s = Solver()
    s.set("timeout", timeout_ms)

    inc = {v: [] for v in range(TOTAL)}
    for i, (u, v) in enumerate(EDGES):
        inc[u].append(i)
        inc[v].append(i)
    for v, lst in inc.items():
        s.add(PbEq([(xvars[i], 1) for i in lst], 2))

    for e in forced_edges:
        s.add(xvars[e] == True)

    for it in range(max_iters):
        res = s.check()
        if res == unsat:
            return "proven", it + 1
        if res != sat:
            return "unverified", it + 1
        m = s.model()
        active = []
        bits = np.zeros(E, dtype=np.uint8)
        for i, (u, v) in enumerate(EDGES):
            if is_true(m.evaluate(xvars[i])):
                active.append((u, v))
                bits[i] = 1
        ok, _ = is_single_tour(active)
        if ok:
            return "false_pos", it + 1
        s.add(Or([xvars[i] != bool(bits[i]) for i in range(E)]))

    return "unverified", max_iters


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--max-p-zero", type=int, default=999_999,
                   help="máx triples com P=0 a verificar")
    p.add_argument("--sample-nonzero", type=int, default=50,
                   help="amostra triples com P>0 a verificar")
    p.add_argument("--timeout-ms", type=int, default=30000)
    p.add_argument("--max-iters", type=int, default=100)
    p.add_argument("--seed", type=int, default=2026)
    args = p.parse_args()

    DATA.mkdir(parents=True, exist_ok=True)

    data = json.loads((DATA / "triples_candidates.json").read_text())
    cands = data["candidates"]
    print(f"Total candidatas: {len(cands)}")

    zero_p = [c for c in cands if c["P_triple"] == 0]
    nonzero = [c for c in cands if c["P_triple"] > 0]
    print(f"  com P=0       : {len(zero_p)}")
    print(f"  com P>0       : {len(nonzero)}")

    to_verify = zero_p[: args.max_p_zero]
    if nonzero:
        import random
        rng = random.Random(args.seed)
        n_sample = min(args.sample_nonzero, len(nonzero))
        to_verify += rng.sample(nonzero, n_sample)

    print(f"\nVerificando {len(to_verify)} candidatas com Z3 ...")
    verdicts = {"proven": 0, "false_pos": 0, "unverified": 0}
    results = []
    t0 = time.perf_counter()
    for k, c in enumerate(to_verify):
        verdict, iters = verify_unsat(c["edges"],
                                       timeout_ms=args.timeout_ms,
                                       max_iters=args.max_iters)
        verdicts[verdict] += 1
        c_out = dict(c)
        c_out["z3_status"] = verdict
        c_out["z3_iters"] = iters
        results.append(c_out)
        if (k + 1) % 20 == 0 or k + 1 == len(to_verify):
            elapsed = time.perf_counter() - t0
            print(f"  [{k+1}/{len(to_verify)}] last={verdict:>12s} "
                  f"({iters} iters) "
                  f"proven={verdicts['proven']} "
                  f"false_pos={verdicts['false_pos']} "
                  f"unverified={verdicts['unverified']}  "
                  f"elapsed={elapsed:.1f}s")

    elapsed = time.perf_counter() - t0
    print(f"\nTempo total: {elapsed:.1f}s  ({elapsed/len(to_verify):.2f}s/cand)")
    print(f"\nVeredictos finais:")
    print(f"  proven      : {verdicts['proven']}")
    print(f"  false_pos   : {verdicts['false_pos']}")
    print(f"  unverified  : {verdicts['unverified']}")

    with open(DATA / "triples_verified.json", "w") as f:
        json.dump({
            "args": vars(args),
            "totals": verdicts,
            "elapsed_s": round(elapsed, 2),
            "results": results,
        }, f, indent=2)
    print(f"\nSalvo: data/triples_verified.json")

    proven = [r for r in results if r["z3_status"] == "proven"]
    if proven:
        print(f"\nTriplas PROVEN (top 10):")
        for r in proven[:10]:
            print(f"  {r['edge_names']}  P={r['P_triple']:.5f}  "
                  f"src={r.get('source', '—')}")


if __name__ == "__main__":
    main()
