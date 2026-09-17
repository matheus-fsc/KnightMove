"""
Micro-benchmark: rollback vs cópia.

Cenário: 10×10, replicar padrão de fix+undo do backtracking real
(~1908 nós, ~10k fixes/undos no agregado v1). Mede tempo total e
operações por segundo.
"""
from __future__ import annotations

import json
import random
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from subtour_detector import (
    SubtourDetectorRollback, SubtourDetectorCopy, OK,
)


def bench(DetCls, edges_uv, V, scenarios, label):
    d = DetCls(edges_uv, n_vertices=V)
    t0 = time.perf_counter()
    n_ops = 0
    for plan in scenarios:
        # plan = lista de arestas (idx) a fixar em sequência, depois rollback
        cp = d.checkpoint()
        for e_idx in plan:
            d.fix(e_idx)
            n_ops += 1
        d.rollback(cp)
    dt = time.perf_counter() - t0
    print(f"  {label:28s}  ops={n_ops:6d}  t={dt*1000:7.2f}ms  ops/ms={n_ops/(dt*1000):6.1f}")
    return dt


def main():
    edges_json = json.loads(
        (ROOT.parent / "board_10x10" / "data" / "edges_10x10.json").read_text()
    )
    edges_uv = [tuple(uv) for uv in edges_json["edges_uv"]]
    V = 100
    samples = np.load(
        ROOT.parent / "board_10x10" / "data" / "samples"
        / "tours_10x10_batch_000.npy"
    )
    rng = random.Random(0)

    # cenário: simular ~2000 nós (≈ v1) cada um fixando ~50 arestas e revertendo
    scenarios = []
    for _ in range(2000):
        ti = rng.randrange(samples.shape[0])
        edges_in_tour = [i for i in range(len(edges_uv)) if samples[ti][i] == 1]
        rng.shuffle(edges_in_tour)
        scenarios.append(edges_in_tour[:60])

    print("Bench: 2000 nós × ~60 arestas/nó × fix+rollback")
    bench(SubtourDetectorRollback, edges_uv, V, scenarios, "Rollback")
    bench(SubtourDetectorCopy, edges_uv, V, scenarios, "Copy")


if __name__ == "__main__":
    main()
