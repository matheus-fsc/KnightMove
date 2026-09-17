"""Roda n=14 com K=500 e log de progresso a cada 50 tours.

Reaproveita MinimalV2 e instrumenta para registrar (nodes, tours, t)
quando um tour é encontrado.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from scaling_minimal_v2 import MinimalV2, FREE  # noqa: E402


class MinimalV2Logged(MinimalV2):
    def __init__(self, n: int, alvo: int, timeout: float, log_step: int = 50):
        super().__init__(n=n, alvo=alvo, timeout=timeout)
        self.log_step = log_step
        self.node_log: list[tuple[int, int, float]] = []
        self.t_first: float | None = None

    def branch(self):
        if self.parar:
            return
        self.nos += 1
        if (self.nos & 0xFFFF) == 0:
            if time.perf_counter() - self.t0 > self.timeout:
                self.parar = True
                return

        snap_est = self.estado.copy()
        snap_seen = self.last_seen.copy()
        cp = self.detector.checkpoint()

        try:
            self.propagar()
        except ValueError:
            self.podas["R2"] += 1
            self.estado = snap_est
            self.last_seen = snap_seen
            self.detector.rollback(cp)
            return

        from subtour_detector import OK  # type: ignore
        if self.sync_det() != OK:
            self.estado = snap_est
            self.last_seen = snap_seen
            self.detector.rollback(cp)
            return

        e = self.escolher_var()
        if e == -1:
            self.folhas += 1
            if (self.estado == 1).sum() == self.V:
                self.fatores2 += 1
                if self.conexo():
                    self.tours += 1
                    if self.t_first is None:
                        self.t_first = time.perf_counter() - self.t0
                    if (
                        self.tours % self.log_step == 0
                        or self.tours == 1
                        or self.tours == self.alvo
                    ):
                        dt = time.perf_counter() - self.t0
                        self.node_log.append((self.nos, self.tours, dt))
                        print(
                            f"  [t={dt:7.2f}s] tours={self.tours:4d}  nós={self.nos:8d}  "
                            f"nós/tour={self.nos/max(1,self.tours):8.2f}",
                            flush=True,
                        )
                    if self.tours >= self.alvo:
                        self.parar = True
            self.estado = snap_est
            self.last_seen = snap_seen
            self.detector.rollback(cp)
            return

        for val in (1, 0):
            if self.parar:
                break
            inner_est = self.estado.copy()
            inner_seen = self.last_seen.copy()
            inner_cp = self.detector.checkpoint()
            self.estado[e] = val
            self.branch()
            self.estado = inner_est
            self.last_seen = inner_seen
            self.detector.rollback(inner_cp)

        self.estado = snap_est
        self.last_seen = snap_seen
        self.detector.rollback(cp)


def main():
    n = 14
    alvo = 500
    timeout = 600.0

    print(f"=== n={n}  alvo={alvo}  timeout={timeout}s ===", flush=True)
    v2 = MinimalV2Logged(n=n, alvo=alvo, timeout=timeout, log_step=50)
    print(
        f"V={v2.V}  E={v2.E}  β₁={v2.E-v2.V+1}  "
        f"mand={len(v2.mand_idx)}  n_free_inicial={int((v2.estado==FREE).sum())}",
        flush=True,
    )
    r = v2.executar()

    out = {
        "n": n,
        "V": v2.V,
        "E": v2.E,
        "beta1": v2.E - v2.V + 1,
        "n_free": int(r.get("n_free_inicial", 0)),
        "Q": 3,
        "alvo": alvo,
        "K_found": r["n_tours"],
        "t_total": r["tempo_s"],
        "t_first": v2.t_first,
        "nodes_total": r["n_nos"],
        "nodes_per_tour": r["nos_por_tour"],
        "ratio_2fat_tours": r["razao_2fat_tour"],
        "n_2fatores": r["n_2fatores"],
        "n_folhas": r["n_folhas"],
        "podas": r["podas"],
        "parou_por_timeout": r["parou_por_timeout"],
        "node_log": [
            {"n_nodes": int(a), "n_tours": int(b), "t": float(c)}
            for (a, b, c) in v2.node_log
        ],
    }
    out_path = ROOT / "data" / "scaling_n14.json"
    out_path.write_text(json.dumps(out, indent=2))
    print(f"\nSalvo em {out_path}", flush=True)
    print(
        f"\nRESUMO: K_found={out['K_found']}  t_total={out['t_total']:.2f}s  "
        f"nodes/tour={out['nodes_per_tour']:.2f}  ratio={out['ratio_2fat_tours']:.3f}",
        flush=True,
    )


if __name__ == "__main__":
    main()
