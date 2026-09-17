#!/usr/bin/env python3
"""
benchmark_z3.py
===============
Benchmark de eficiência do Z3 GF(2) no 6×6.

Duas fontes de dados:

  A) Agregador dos batches já existentes em data_6x6_sampled/
     - 400 batches × 200 amostras × 8 pares
     - Cada batch grava n_attempts, n_found, elapsed_s
     - Permite extrair attempts/sol e sols/s sem re-executar

  B) Runs frescos via engine_6x6_sampler.sample_signatures
     - Para K controlado (10, 100, 1000) com medição precisa
     - Mede time-to-first via instrumentação fina do solver
     - Roda também com invariantes (H4) — ver --invariants

Saída:
  benchmark/results/benchmark_z3_aggregated.json
  benchmark/results/benchmark_z3_fresh.json
"""

import argparse
import json
import sys
import time
import tracemalloc
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "cavalo_8x8"))
sys.path.insert(0, str(ROOT))

from engine_6x6_sampler import (  # noqa: E402
    BOARD, build_graph, get_ciclos_de_4, spanning_tree,
    extract_loops, get_loop_edges, edge_label, node_label,
)
from z3 import Solver, Bool, Sum, If, Or, sat, is_true  # noqa: E402


SAMPLED_DIR = ROOT / "cavalo_8x8" / "data_6x6_sampled" / "samples"
OUT_AGG = ROOT / "benchmark" / "results" / "benchmark_z3_aggregated.json"
OUT_FRESH = ROOT / "benchmark" / "results" / "benchmark_z3_fresh.json"

# Pares canônicos usados na validação 6×6 (data_6x6_sampled/config.json).
CANONICAL_PAIRS = [
    [[0, 0], [0, 5]],
    [[0, 0], [5, 4]],
    [[0, 0], [3, 2]],
    [[0, 0], [2, 3]],
    [[0, 1], [4, 4]],
    [[0, 1], [5, 3]],
    [[0, 1], [2, 2]],
    [[1, 2], [5, 3]],
]

# Invariantes 6×6 do summary.json — 8 pares com r ≈ -0.771.
# Esses são os "bifurcations" com correlação negativa muito forte.
INVARIANT_EDGES_LABELS = [
    ("B6-D5", "B6-C4"),
    ("E6-C5", "E6-D4"),
    ("A5-C4", "A5-B3"),
    ("F5-D4", "F5-E3"),
    ("B4-A2", "C3-A2"),
    ("E4-F2", "D3-F2"),
    ("C3-B1", "D2-B1"),
    ("D3-E1", "C2-E1"),
]


# ── A) Agregador dos batches existentes ──────────────────────────────────

def aggregate_existing_batches():
    """Lê todos os batches em data_6x6_sampled/samples e agrega métricas."""
    if not SAMPLED_DIR.exists():
        return None

    files = sorted(SAMPLED_DIR.glob("batch_p*_b*.json"))
    if not files:
        return None

    by_pair = {}
    total = {
        "n_files": 0,
        "n_samples": 0,
        "n_attempts": 0,
        "elapsed_s": 0.0,
        "failed_attempts": 0,
    }

    for f in files:
        try:
            d = json.loads(f.read_text())
        except Exception:
            continue
        pair_key = (tuple(d["start"]), tuple(d["end"]))
        entry = by_pair.setdefault(str(pair_key), {
            "start": d["start"], "end": d["end"],
            "n_batches": 0, "n_samples": 0, "n_attempts": 0,
            "elapsed_s": 0.0, "failed_attempts": 0,
        })
        entry["n_batches"] += 1
        entry["n_samples"] += d["n_found"]
        entry["n_attempts"] += d["n_attempts"]
        entry["elapsed_s"] += d["elapsed_s"]
        entry["failed_attempts"] += (d["n_attempts"] - d["n_found"])

        total["n_files"] += 1
        total["n_samples"] += d["n_found"]
        total["n_attempts"] += d["n_attempts"]
        total["elapsed_s"] += d["elapsed_s"]
        total["failed_attempts"] += (d["n_attempts"] - d["n_found"])

    # taxas derivadas por par
    for entry in by_pair.values():
        n = entry["n_samples"] or 1
        entry["attempts_per_solution"] = round(entry["n_attempts"] / n, 4)
        entry["failed_per_solution"] = round(entry["failed_attempts"] / n, 4)
        entry["s_per_solution"] = round(entry["elapsed_s"] / n, 4)
        entry["solutions_per_second"] = round(n / entry["elapsed_s"], 4)

    n = total["n_samples"] or 1
    total["attempts_per_solution"] = round(total["n_attempts"] / n, 4)
    total["failed_per_solution"] = round(total["failed_attempts"] / n, 4)
    total["s_per_solution"] = round(total["elapsed_s"] / n, 4)
    total["solutions_per_second"] = round(n / total["elapsed_s"], 4)

    return {
        "source": str(SAMPLED_DIR.relative_to(ROOT)),
        "method": "z3_gf2_sym_false",
        "board": BOARD,
        "totals": total,
        "by_pair": by_pair,
    }


# ── B) Runs frescos com medição fina ─────────────────────────────────────

def build_solver(start_t, end_t, *, break_symmetry=False, extra_constraints=None):
    """
    Constrói o solver Z3 espelhando engine_6x6_sampler.sample_signatures.
    Retorna (solver, edges_orig, edges_aug, lvars, ea, edges_s, virtual).
    """
    start = tuple(start_t)
    end = tuple(end_t)

    nodes, adj, edges_orig = build_graph()
    virtual = tuple(sorted([start, end]))
    edges_aug = set(edges_orig) | {virtual}
    adj_aug = {n: list(adj[n]) for n in nodes}
    adj_aug[start].append(end)
    adj_aug[end].append(start)

    par, te = spanning_tree(adj_aug, nodes[0])
    loops = extract_loops(par, te, edges_aug)

    solver = Solver()
    lvars = [Bool(f"L{i}") for i in range(len(loops))]
    e2l = {e: [] for e in edges_aug}
    for i, lp in enumerate(loops):
        for e in get_loop_edges(lp):
            e2l[e].append(i)

    def ea(e):
        idxs = e2l.get(e, [])
        if not idxs:
            return False
        return (Sum([If(lvars[i], 1, 0) for i in idxs]) % 2 == 1)

    for n in nodes:
        inc = [e for e in edges_aug if n in e]
        solver.add(Sum([If(ea(e), 1, 0) for e in inc]) == 2)

    solver.add(ea(virtual) == True)

    if break_symmetry and start == (0, 0):
        edge_sym = tuple(sorted([(0, 0), (1, 2)]))
        if edge_sym in edges_aug:
            solver.add(ea(edge_sym) == True)

    c4 = get_ciclos_de_4(nodes, adj)
    for ciclo in c4:
        es4 = [
            tuple(sorted([ciclo[0], ciclo[1]])),
            tuple(sorted([ciclo[1], ciclo[2]])),
            tuple(sorted([ciclo[2], ciclo[3]])),
            tuple(sorted([ciclo[3], ciclo[0]])),
        ]
        if all(e in edges_aug for e in es4):
            solver.add(Sum([If(ea(e), 1, 0) for e in es4]) <= 3)

    # Constraints extras (ex.: invariantes para H4)
    extra_applied = 0
    if extra_constraints:
        for c in extra_constraints(ea, edges_aug):
            if c is not None:
                solver.add(c)
                extra_applied += 1

    edges_s = sorted(edges_orig)
    return solver, nodes, edges_orig, edges_aug, lvars, ea, edges_s, virtual, extra_applied


def invariant_constraints_factory(invariant_pairs_labels):
    """
    Fábrica de extra_constraints que adiciona NOT(A AND B) para cada par
    invariante com correlação fortemente negativa.

    invariant_pairs_labels: lista de tuplas (label_a, label_b) tipo ("B6-D5","B6-C4")
    """
    def factory(ea, edges_aug):
        # mapa label -> aresta (tuple)
        lbl_to_edge = {}
        for e in edges_aug:
            lbl_to_edge[edge_label(e)] = e
        out = []
        for la, lb in invariant_pairs_labels:
            ea_lbl = la if la in lbl_to_edge else f"{la.split('-')[1]}-{la.split('-')[0]}"
            eb_lbl = lb if lb in lbl_to_edge else f"{lb.split('-')[1]}-{lb.split('-')[0]}"
            if ea_lbl not in lbl_to_edge or eb_lbl not in lbl_to_edge:
                # tenta normalizar pelas pontas
                continue
            from z3 import Not, And
            a, b = lbl_to_edge[ea_lbl], lbl_to_edge[eb_lbl]
            out.append(Not(And(ea(a) == True, ea(b) == True)))
        return out
    return factory


def z3_timed(start, end, target_k, *,
             max_attempts_factor=10,
             timeout_s=None,
             break_symmetry=False,
             use_invariants=False,
             record_solution_times=False,
             verbose=False):
    """
    Roda Z3 amostrando ciclos hamiltonianos start→end até target_k soluções.
    Mede tempo, attempts, time-to-first, time-to-K e memória pico.
    """
    extra = invariant_constraints_factory(INVARIANT_EDGES_LABELS) if use_invariants else None
    tracemalloc.start()
    t_setup0 = time.perf_counter()
    (solver, nodes, edges_orig, edges_aug, lvars, ea, edges_s,
     virtual, n_extra) = build_solver(start, end,
                                      break_symmetry=break_symmetry,
                                      extra_constraints=extra)
    t_setup = time.perf_counter() - t_setup0

    max_attempts = max_attempts_factor * target_k

    sigs = []
    attempts = 0
    failed_attempts = 0
    exhausted = False
    timed_out = False
    time_to_first = None
    solution_times = [] if record_solution_times else None

    t0 = time.perf_counter()
    while len(sigs) < target_k:
        if attempts >= max_attempts:
            break
        if timeout_s is not None and (time.perf_counter() - t0) > timeout_s:
            timed_out = True
            break

        attempts += 1
        if solver.check() != sat:
            exhausted = True
            break

        m = solver.model()
        active = [e for e in edges_orig if is_true(m.evaluate(ea(e)))]
        adj_r = {n: [] for n in nodes}
        for u, v in active:
            adj_r[u].append(v)
            adj_r[v].append(u)

        start_t = tuple(start)
        vis = {start_t}
        q = [start_t]
        while q:
            c = q.pop(0)
            for nb in adj_r[c]:
                if nb not in vis:
                    vis.add(nb)
                    q.append(nb)

        if len(vis) == len(nodes):
            sigs.append([is_true(m.evaluate(ea(e))) for e in edges_s])
            now = time.perf_counter() - t0
            if time_to_first is None:
                time_to_first = now
            if record_solution_times:
                solution_times.append(now)
            if verbose and len(sigs) % 100 == 0:
                print(f"  -> {len(sigs)}/{target_k}  "
                      f"attempts={attempts}  elapsed={now:.2f}s",
                      flush=True)
        else:
            failed_attempts += 1

        solver.add(Or([lvars[i] != is_true(m.evaluate(lvars[i]))
                       for i in range(len(lvars))]))

        if len(vis) < len(nodes):
            seen = set(vis)
            for seed in [n for n in nodes if n not in seen]:
                if seed in seen:
                    continue
                comp = {seed}
                q2 = [seed]
                while q2:
                    c = q2.pop(0)
                    for nb in adj_r[c]:
                        if nb not in comp:
                            comp.add(nb)
                            q2.append(nb)
                in_e = [e for e in edges_aug if e[0] in comp and e[1] in comp]
                solver.add(Sum([If(ea(e), 1, 0) for e in in_e]) <= len(comp) - 1)
                seen |= comp

    elapsed = time.perf_counter() - t0
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    n = len(sigs)
    result = {
        "method": "z3_gf2" + ("_invariants" if use_invariants else ""),
        "board": BOARD,
        "start": list(start), "end": list(end),
        "target_k": target_k,
        "found_k": n,
        "setup_s": round(t_setup, 4),
        "elapsed_s": round(elapsed, 4),
        "n_attempts": attempts,
        "failed_attempts": failed_attempts,
        "failed_per_solution": round(failed_attempts / n, 4) if n else None,
        "attempts_per_solution": round(attempts / n, 4) if n else None,
        "peak_memory_mb": round(peak / 1024 / 1024, 3),
        "solutions_per_second": round(n / elapsed, 4) if elapsed > 0 and n else None,
        "time_to_first_s": (round(time_to_first, 6)
                            if time_to_first is not None else None),
        "exhausted": exhausted,
        "timed_out": timed_out,
        "break_symmetry": break_symmetry,
        "use_invariants": use_invariants,
        "extra_constraints_applied": n_extra,
    }
    if record_solution_times:
        result["solution_times_s"] = [round(t, 6) for t in solution_times]
        # curva time-to-K
        out_curve = {}
        for k in (1, 10, 100, 1000):
            if k <= len(solution_times):
                out_curve[str(k)] = round(solution_times[k - 1], 6)
            else:
                out_curve[str(k)] = None
        result["time_to_k_s"] = out_curve

    # também guarda assinaturas (úteis p/ KL-divergência)
    result["signatures"] = sigs
    result["edges"] = [edge_label(e) for e in edges_s]
    return result


# ── runner ────────────────────────────────────────────────────────────────

def run_fresh_benchmark(pairs, target_ks, *, use_invariants=False,
                        save_signatures=True, verbose=True):
    out = {
        "method": "z3_gf2" + ("_invariants" if use_invariants else ""),
        "board": BOARD,
        "target_ks": list(target_ks),
        "pairs": [],
    }
    for pair in pairs:
        start, end = tuple(pair[0]), tuple(pair[1])
        pair_results = {"start": list(start), "end": list(end), "runs": []}
        for k in target_ks:
            if verbose:
                tag = "+inv" if use_invariants else "    "
                print(f"  [{tag}] par {start}→{end}  K={k} ...", flush=True)
            r = z3_timed(start, end, k, record_solution_times=True,
                         use_invariants=use_invariants, verbose=False)
            if not save_signatures:
                r.pop("signatures", None)
            pair_results["runs"].append(r)
            if verbose:
                print(f"        found={r['found_k']}/{k}  "
                      f"elapsed={r['elapsed_s']:.2f}s  "
                      f"attempts/sol={r['attempts_per_solution']}  "
                      f"failed/sol={r['failed_per_solution']}  "
                      f"t_first={r['time_to_first_s']}s",
                      flush=True)
        out["pairs"].append(pair_results)
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--aggregate-only", action="store_true",
                   help="Só agregar batches existentes (rápido)")
    p.add_argument("--fresh-only", action="store_true",
                   help="Só rodar Z3 fresco (pula agregação)")
    p.add_argument("--invariants", action="store_true",
                   help="Rodar variante H4 com invariantes injetados")
    p.add_argument("--ks", type=int, nargs="+", default=[10, 100, 1000],
                   help="Valores de K para o run fresco")
    p.add_argument("--n-pairs", type=int, default=3,
                   help="Quantos pares canônicos usar (até 8)")
    p.add_argument("--no-save-signatures", action="store_true",
                   help="Não gravar assinaturas (economiza espaço)")
    args = p.parse_args()

    do_agg = not args.fresh_only
    do_fresh = not args.aggregate_only

    if do_agg:
        print("=" * 65)
        print("Z3 6×6 — agregação dos batches existentes")
        print("=" * 65)
        agg = aggregate_existing_batches()
        if agg is None:
            print("  (nenhum batch encontrado em data_6x6_sampled/samples)")
        else:
            OUT_AGG.parent.mkdir(parents=True, exist_ok=True)
            with open(OUT_AGG, "w") as f:
                json.dump(agg, f, indent=2)
            t = agg["totals"]
            print(f"  arquivos      : {t['n_files']}")
            print(f"  amostras      : {t['n_samples']:,}")
            print(f"  attempts      : {t['n_attempts']:,}")
            print(f"  tempo total   : {t['elapsed_s']:.1f}s")
            print(f"  attempts/sol  : {t['attempts_per_solution']}")
            print(f"  failed/sol    : {t['failed_per_solution']}")
            print(f"  sols/s        : {t['solutions_per_second']}")
            print(f"  Salvo: {OUT_AGG.relative_to(ROOT)}")
        print()

    if do_fresh:
        print("=" * 65)
        print(f"Z3 6×6 — run fresco (invariants={args.invariants})")
        print("=" * 65)
        pairs = CANONICAL_PAIRS[: args.n_pairs]
        out = run_fresh_benchmark(
            pairs, args.ks,
            use_invariants=args.invariants,
            save_signatures=not args.no_save_signatures,
        )
        target_path = OUT_FRESH
        if args.invariants:
            target_path = OUT_FRESH.with_name("benchmark_z3_fresh_invariants.json")
        target_path.parent.mkdir(parents=True, exist_ok=True)
        with open(target_path, "w") as f:
            json.dump(out, f, indent=2)
        print(f"  Salvo: {target_path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
