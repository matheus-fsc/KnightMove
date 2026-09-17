"""Benchmark f∞ = 0.25 fixo vs. tabela empírica.

Compara 3 métodos:
  A) BT_original: tabela empírica {0:0.528, 1:0.193, 2:0.198, 3:0.294, 4:0.247, 5:0.261}
  B) BT_torus:    L≥4 forçado a 0.25 exato
  C) BT_uniform:  f∞ = 0.25 para todos os L
"""

# --- bootstrap: torna core/ importavel a partir de qualquer subpasta ---
import sys as _sys, pathlib as _pathlib
for _anc in _pathlib.Path(__file__).resolve().parents:
    if (_anc / "core").is_dir() and (_anc / "pyproject.toml").exists():
        _sys.path.insert(0, str(_anc / "core"))
        break
# --- fim do bootstrap ---

import time
import json
import sys
import os

sys.path.insert(0, os.path.dirname(__file__))
import knight_tours as kt
import numpy as np

CONFIGS = {
    "A": {0: 0.528, 1: 0.193, 2: 0.198, 3: 0.294, 4: 0.247, 5: 0.261},
    "B": {0: 0.528, 1: 0.193, 2: 0.198, 3: 0.294, 4: 0.25,  5: 0.25},
    "C": {0: 0.25,  1: 0.25,  2: 0.25,  3: 0.25,  4: 0.25,  5: 0.25},
}
DEFAULTS = {"A": 0.25, "B": 0.25, "C": 0.25}

K = 100
RUNS = 3
BOARDS = [12, 14, 16, 18]

def frac_L_ge4(n):
    count = 0
    for r in range(n):
        for c in range(n):
            if min(r, c, n-1-r, n-1-c) >= 4:
                count += 1
    return count / (n * n)

results = {}

for n in BOARDS:
    results[n] = {"frac_L4": round(frac_L_ge4(n) * 100)}
    for label, finf in CONFIGS.items():
        times = []
        for run in range(RUNS):
            kt.F_INF = dict(finf)
            kt.F_INF_DEFAULT = DEFAULTS[label]
            seed = 42 + run * 1000
            t0 = time.perf_counter()
            tours = kt.knight_tours(n, K, seed=seed)
            t1 = time.perf_counter()
            if len(tours) < K:
                times.append(float('inf'))
            else:
                times.append((t1 - t0) * 1000 / K)
            print(f"  n={n} method={label} run={run} K={len(tours)} "
                  f"t={t1-t0:.2f}s ms/tour={times[-1]:.2f}", flush=True)
        med = sorted(times)[len(times) // 2]
        results[n][label] = round(med, 2)

kt.F_INF = {0: 0.528, 1: 0.193, 2: 0.198, 3: 0.294, 4: 0.247, 5: 0.261}
kt.F_INF_DEFAULT = 0.25

print("\n" + "=" * 70)
print(f"{'n':>4} | {'%L≥4':>5} | {'t_A ms':>8} | {'t_B ms':>8} | {'t_C ms':>8} | {'B/A':>5} | {'C/A':>5}")
print("-" * 70)
for n in BOARDS:
    r = results[n]
    tA, tB, tC = r["A"], r["B"], r["C"]
    ba = tB / tA if tA > 0 else float('inf')
    ca = tC / tA if tA > 0 else float('inf')
    print(f"{n:>4} | {r['frac_L4']:>4}% | {tA:>8.2f} | {tB:>8.2f} | {tC:>8.2f} | {ba:>5.2f} | {ca:>5.2f}")

print("=" * 70)

with open("benchmark_f025_results.json", "w") as f:
    json.dump(results, f, indent=2)
print("Resultados salvos em benchmark_f025_results.json")
