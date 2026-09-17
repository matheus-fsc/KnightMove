#!/usr/bin/env python3
"""
baseline_bt_6x6.py
==================
Baseline de eficiência do backtracking 6×6.

Re-roda o backtracking de cavalo_loop_destruicao_6x6.py adicionando:
  - contagem precisa de dead-ends (becos sem saída)
  - contagem de nós visitados (chamadas recursivas)
  - timestamp de cada solução encontrada (para curva time-to-K)
  - pico de memória via tracemalloc

Roda em três variantes:
  - bt_warnsdorff      : ordenação por grau ascendente (a do código original)
  - bt_random_order    : ordem arbitrária (sem heurística) — controle
  - bt_warnsdorff_v0   : só busca a partir de v=0 (19.724 dirigidos),
                         o resto é gerado por rotação como no script original

Saída:
  benchmark/results/baseline_bt_6x6.json
"""

import json
import os
import resource
import sys
import time
from pathlib import Path

# importa o grafo do módulo original
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from cavalo_loop_destruicao_6x6 import (  # noqa: E402
    ADJ, EDGES_LIST, TOTAL, BOARD, label
)

OUT_PATH = ROOT / "benchmark" / "results" / "baseline_bt_6x6.json"


def bt_run(start, *, use_warnsdorff=True, target_k=None,
           record_solution_times=False, progress_every=2_000_000,
           timeout_s=None):
    """
    Backtracking instrumentado a partir de `start`.

    Retorna:
      {
        "method": str,
        "start": int,
        "target_k": int|None,
        "n_solutions": int,
        "elapsed_s": float,
        "dead_ends": int,
        "nodes_visited": int,
        "peak_memory_mb": float,
        "time_to_first_s": float|None,
        "solution_times": list[float] (se record_solution_times=True),
      }
    """
    solutions_count = 0
    dead_ends = 0
    nodes_visited = 0
    solution_times = [] if record_solution_times else None
    time_to_first = None
    stopped_at_k = False

    path = [start]
    alive = bytearray([1] * TOTAL)
    alive[start] = 0
    deg = [sum(1 for u in ADJ[v] if alive[u]) for v in range(TOTAL)]

    rss_before = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    t0 = time.perf_counter()

    # bandeira mutável para parar quando atingir target_k
    stop = [False]

    def backtrack():
        nonlocal solutions_count, dead_ends, nodes_visited, time_to_first
        if stop[0]:
            return
        nodes_visited += 1
        if nodes_visited % progress_every == 0:
            elapsed = time.perf_counter() - t0
            rate = nodes_visited / elapsed if elapsed > 0 else 0
            print(f"    [progresso] nós={nodes_visited:>12,}  "
                  f"sols={solutions_count:>6,}  "
                  f"dead_ends={dead_ends:>10,}  "
                  f"taxa={rate:>9,.0f}/s  elapsed={elapsed:5.0f}s",
                  flush=True)
            if timeout_s is not None and elapsed > timeout_s:
                stop[0] = True
                return

        n = len(path)
        current = path[-1]

        if n == TOTAL:
            if start in ADJ[current]:
                solutions_count += 1
                if record_solution_times:
                    solution_times.append(time.perf_counter() - t0)
                if time_to_first is None:
                    time_to_first = time.perf_counter() - t0
                if target_k is not None and solutions_count >= target_k:
                    stop[0] = True
            else:
                dead_ends += 1
            return

        candidates = [u for u in ADJ[current] if alive[u]]
        if not candidates:
            dead_ends += 1
            return

        if use_warnsdorff:
            candidates.sort(key=lambda u: deg[u])
            # poda: candidato com deg=0 só serve no passo de fechamento
            if deg[candidates[0]] == 0 and n < TOTAL - 1:
                dead_ends += 1
                return

        for nxt in candidates:
            if stop[0]:
                return
            alive[nxt] = 0
            for w in ADJ[nxt]:
                deg[w] -= 1
            path.append(nxt)
            backtrack()
            path.pop()
            alive[nxt] = 1
            for w in ADJ[nxt]:
                deg[w] += 1

    backtrack()
    elapsed = time.perf_counter() - t0
    rss_after = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    # ru_maxrss em KB no Linux; pegamos o pico do processo desde o início
    peak_kb = rss_after

    method = "bt_warnsdorff" if use_warnsdorff else "bt_no_heuristic"
    aborted_timeout = stop[0] and (target_k is None or solutions_count < target_k)
    result = {
        "method": method,
        "board": BOARD,
        "start_vertex": start,
        "start_label": label(start),
        "target_k": target_k,
        "timeout_s": timeout_s,
        "aborted_timeout": aborted_timeout,
        "n_solutions": solutions_count,
        "elapsed_s": round(elapsed, 4),
        "dead_ends": dead_ends,
        "nodes_visited": nodes_visited,
        "peak_memory_mb": round(peak_kb / 1024, 3),
        "time_to_first_s": (round(time_to_first, 6)
                            if time_to_first is not None else None),
        "dead_ends_per_solution": (round(dead_ends / solutions_count, 4)
                                   if solutions_count > 0 else None),
        "nodes_per_solution": (round(nodes_visited / solutions_count, 4)
                               if solutions_count > 0 else None),
        "solutions_per_second": (round(solutions_count / elapsed, 4)
                                 if elapsed > 0 else None),
        "stopped_at_target": (target_k is not None
                              and solutions_count >= target_k),
    }
    if record_solution_times:
        result["solution_times_s"] = [round(t, 6) for t in solution_times]
    return result


def time_to_k_curve(solution_times, ks=(1, 10, 100, 1000, 10_000)):
    """Extrai pontos da curva time-to-K para os valores em `ks`."""
    out = {}
    for k in ks:
        if k <= len(solution_times):
            out[str(k)] = round(solution_times[k - 1], 6)
        else:
            out[str(k)] = None
    return out


def main():
    print("=" * 65)
    print("BASELINE BACKTRACKING 6×6")
    print("=" * 65)
    print(f"V={TOTAL}  E={len(EDGES_LIST)}")
    print()

    runs = []

    # ── 1. BT + Warnsdorff target_k=1000 (curva rápida) ─────────────────
    print("[1/3] BT + Warnsdorff a partir de v=0, target_k=1000...")
    r1 = bt_run(0, use_warnsdorff=True, target_k=1000,
                record_solution_times=True)
    r1["time_to_k_s"] = time_to_k_curve(r1["solution_times_s"])
    sol_times_1 = r1.pop("solution_times_s")
    runs.append(r1)
    print(f"      → {r1['n_solutions']} sols em {r1['elapsed_s']}s")
    print(f"      → dead-ends={r1['dead_ends']:,}  "
          f"de/sol={r1['dead_ends_per_solution']}  "
          f"sols/s={r1['solutions_per_second']}")
    print()

    # ── 2. BT puro (sem Warnsdorff) — controle, com timeout ─────────────
    # Limite: 120s — só para ver quanto "tempo perdido" custa não usar
    # heurística. Resultado fica parcial se não atingir K.
    print("[2/3] BT sem heurística a partir de v=0, target_k=1000, timeout=120s...")
    r2 = bt_run(0, use_warnsdorff=False, target_k=1000,
                record_solution_times=True, timeout_s=120)
    r2["time_to_k_s"] = time_to_k_curve(r2["solution_times_s"])
    sol_times_2 = r2.pop("solution_times_s")
    runs.append(r2)
    print(f"      → {r2['n_solutions']} sols em {r2['elapsed_s']}s "
          f"{'(TIMEOUT)' if r2['aborted_timeout'] else ''}")
    print(f"      → dead-ends={r2['dead_ends']:,}  "
          f"de/sol={r2['dead_ends_per_solution']}  "
          f"sols/s={r2['solutions_per_second']}")
    print()

    # ── 3. BT + Warnsdorff exaustivo (target_k=None) ────────────────────
    # Roda até esgotar a busca a partir de v=0. Este é o run que dá o
    # baseline "completo" comparável aos 19.724 ciclos dirigidos do script
    # original. Timeout 1800s = 30min como margem.
    print("[3/3] BT + Warnsdorff a partir de v=0, exaustivo (timeout 1800s)...")
    r3 = bt_run(0, use_warnsdorff=True, record_solution_times=True,
                timeout_s=1800)
    r3["time_to_k_s"] = time_to_k_curve(
        r3["solution_times_s"], ks=(1, 10, 100, 1000, 5000, 10000, 19000))
    sol_times_3 = r3.pop("solution_times_s")
    runs.append(r3)
    print(f"      → {r3['n_solutions']} ciclos dirigidos em {r3['elapsed_s']}s "
          f"{'(TIMEOUT)' if r3['aborted_timeout'] else ''}")
    print(f"      → dead-ends={r3['dead_ends']:,}  "
          f"de/sol={r3['dead_ends_per_solution']}  "
          f"sols/s={r3['solutions_per_second']}")
    print()

    # ── salva resultado consolidado ───────────────────────────────────────
    output = {
        "board": BOARD,
        "total_directed_cycles_at_v0": r3["n_solutions"],
        "total_solutions_with_rotation": r3["n_solutions"] * TOTAL,
        "h1_initial": 45,
        "notes": (
            "r3 é o run exaustivo a partir de v=0 (com Warnsdorff). "
            "Para comparar com o sumário original (710.064 soluções "
            "com rotação) multiplique r3.n_solutions por 36."
        ),
        "runs": runs,
    }

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_PATH, "w") as f:
        json.dump(output, f, indent=2)

    # salva os timestamps detalhados separadamente (podem ser grandes)
    detail_path = OUT_PATH.parent / "baseline_bt_6x6_solution_times.json"
    with open(detail_path, "w") as f:
        json.dump({
            "warnsdorff_k1000": sol_times_1,
            "no_heuristic_k1000_timeout120": sol_times_2,
            "warnsdorff_exhaustive": sol_times_3,
        }, f)

    print(f"Salvo: {OUT_PATH.relative_to(ROOT)}")
    print(f"Salvo: {detail_path.relative_to(ROOT)}")
    print()
    print("─" * 65)
    print("RESUMO")
    print("─" * 65)
    for r in runs:
        de = r['dead_ends_per_solution']
        de_str = f"{de:>10.3f}" if isinstance(de, (int, float)) else f"{'—':>10s}"
        print(f"  {r['method']:20s}  target={str(r['target_k']):>6s}  "
              f"sols={r['n_solutions']:>6,}  "
              f"tempo={r['elapsed_s']:>8.2f}s  "
              f"dead-ends/sol={de_str}")
    print()


if __name__ == "__main__":
    main()
