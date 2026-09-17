"""
T1 — logged_backtracking: v2 mínimo com event log por nó.

Adaptado de `incremental_subtour/scaling_minimal_v2.py` com a adição de um
event_log que registra, para CADA chamada de branch(), o estado e o desfecho:

  {
    "event_id": int,         # ordem DFS de visita
    "n": int,
    "depth": float,          # n_fixadas / E (proporção de arestas fixadas)
    "n_components": int,     # nº de componentes no UF (arestas == 1)
    "max_degree_free": int,  # max grau livre sobre vértices ainda livres
    "L_current": int,        # nível do vértice escolhido neste nó (-1 se podou)
    "L_parent": int,         # nível do vértice escolhido no nó pai (-1 na raiz)
    "score_pressure": float, # pressão do vértice alvo; NaN se podou antes
    "outcome": str           # TOUR | SUBTOUR_EARLY | R2 | CONTRADICTION_DET | RECURSE_INTERNAL
  }

Saída: data/event_log_n{n}.json — lista de dicts.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from collections import deque
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
INCREMENTAL = ROOT.parent / "incremental_subtour"
sys.path.insert(0, str(INCREMENTAL))
from subtour_detector import (
    SubtourDetectorCopy, OK, SUBTOUR, COMPLETE_TOUR, CONTRADICTION,
)

FREE = -1


def build_knight_graph(n: int):
    V = n * n
    pos = [(r, c) for r in range(n) for c in range(n)]
    moves = [(1, 2), (2, 1), (-1, 2), (-2, 1),
             (1, -2), (2, -1), (-1, -2), (-2, -1)]
    edges = set()
    for i, (r, c) in enumerate(pos):
        for dr, dc in moves:
            nr, nc = r + dr, c + dc
            if 0 <= nr < n and 0 <= nc < n:
                j = nr * n + nc
                edges.add((min(i, j), max(i, j)))
    edges_uv = sorted(edges)
    v2e: dict[int, list[int]] = {v: [] for v in range(V)}
    for ei, (u, v) in enumerate(edges_uv):
        v2e[u].append(ei)
        v2e[v].append(ei)
    return edges_uv, V, v2e


def mandatory_corner_edges(n: int, v2e: dict) -> list[int]:
    corners = [0, n - 1, n * (n - 1), n * n - 1]
    mand = set()
    for c in corners:
        if len(v2e[c]) == 2:
            mand.update(v2e[c])
    return sorted(mand)


def vertex_level(v: int, n: int) -> int:
    """Anel do vértice no tabuleiro n×n: L = min(r, c, n-1-r, n-1-c)."""
    r, c = divmod(v, n)
    return int(min(r, c, n - 1 - r, n - 1 - c))


class LoggedMinimalV2:
    def __init__(self, n: int, alvo: int, timeout: float):
        self.n = n
        self.V = n * n
        self.alvo = alvo
        self.timeout = timeout

        self.edges_uv, V, self.v2e = build_knight_graph(n)
        self.E = len(self.edges_uv)
        max_grau = max(len(es) for es in self.v2e.values())
        self.v2e_mat = np.full((self.V, max_grau), -1, dtype=np.int32)
        self.v2e_len = np.zeros(self.V, dtype=np.int32)
        for v, es in self.v2e.items():
            self.v2e_mat[v, : len(es)] = es
            self.v2e_len[v] = len(es)
        self.grau_arr = np.array(
            [len(self.v2e[v]) for v in range(self.V)], dtype=np.int32
        )
        self.level_arr = np.array(
            [vertex_level(v, n) for v in range(self.V)], dtype=np.int32
        )

        self.estado = np.full(self.E, FREE, dtype=np.int8)
        self.mand_idx = mandatory_corner_edges(n, self.v2e)
        for e in self.mand_idx:
            self.estado[e] = 1

        self.detector = SubtourDetectorCopy(self.edges_uv, n_vertices=self.V)
        for e in self.mand_idx:
            self.detector.fix(int(e))
        self.last_seen = self.estado.copy()

        self.nos = 0
        self.tours = 0
        self.fatores2 = 0
        self.folhas = 0
        self.podas = {"R2": 0, "SUBTOUR_EARLY": 0, "CONTRADICTION_DET": 0}
        self.parar = False
        self.t0 = 0.0

        self.event_log: list[dict] = []

    # ------------------------------------------------------------------
    # propagação R2 (idêntica a scaling_minimal_v2)
    # ------------------------------------------------------------------
    def propagar_R2(self) -> int:
        s = self.estado
        mat = self.v2e_mat
        vmask = mat >= 0
        idx = np.where(vmask, mat, 0)
        vals = s[idx]
        vals = np.where(vmask, vals, 99)
        n_um = (vals == 1).sum(axis=1)
        n_zero = (vals == 0).sum(axis=1)
        if (n_um > 2).any():
            raise ValueError("R2: n_um>2")
        if (self.grau_arr - n_zero < 2).any():
            raise ValueError("R2: grau efetivo <2")
        n = 0
        for v in np.where(n_um == 2)[0]:
            row = mat[v, : self.v2e_len[v]]
            for e in row:
                if self.estado[e] == FREE:
                    self.estado[e] = 0
                    n += 1
        for v in np.where(self.grau_arr - n_zero == 2)[0]:
            row = mat[v, : self.v2e_len[v]]
            for e in row:
                if self.estado[e] == FREE:
                    self.estado[e] = 1
                    n += 1
        return n

    def propagar(self) -> None:
        while True:
            if self.propagar_R2() == 0:
                return

    def sync_det(self) -> str:
        diff = np.where((self.estado == 1) & (self.last_seen != 1))[0]
        for e in diff:
            r = self.detector.fix(int(e))
            self.last_seen[e] = 1
            if r == SUBTOUR:
                self.podas["SUBTOUR_EARLY"] += 1
                return SUBTOUR
            if r == CONTRADICTION:
                self.podas["CONTRADICTION_DET"] += 1
                return CONTRADICTION
        self.last_seen = self.estado.copy()
        return OK

    def escolher_var(self) -> tuple[int, int, float]:
        """Retorna (edge_idx, v_alvo, pressao_v). edge_idx=-1 se não há livre."""
        s = self.estado
        mat = self.v2e_mat
        vmask = mat >= 0
        idx = np.where(vmask, mat, 0)
        vals = s[idx]
        vals = np.where(vmask, vals, 99)
        n_um = (vals == 1).sum(axis=1)
        n_zero = (vals == 0).sum(axis=1)
        n_free_v = (vals == -1).sum(axis=1)
        pressao = n_um.astype(np.int64) * 100 + (self.grau_arr - n_zero - 2)
        pressao = np.where(n_free_v > 0, pressao, -10**9)
        if pressao.max() <= -10**8:
            return -1, -1, float("nan")
        v_alvo = int(np.argmax(pressao))
        row = mat[v_alvo, : self.v2e_len[v_alvo]]
        for e in row:
            if self.estado[e] == FREE:
                return int(e), v_alvo, float(pressao[v_alvo])
        return -1, -1, float("nan")

    def conexo(self) -> bool:
        ativos = np.where(self.estado == 1)[0]
        if len(ativos) != self.V:
            return False
        adj = [[] for _ in range(self.V)]
        for e in ativos:
            u, v = self.edges_uv[e]
            adj[u].append(v)
            adj[v].append(u)
        bfs = deque([0])
        seen = {0}
        while bfs:
            x = bfs.popleft()
            for y in adj[x]:
                if y not in seen:
                    seen.add(y)
                    bfs.append(y)
        return len(seen) == self.V

    # ------------------------------------------------------------------
    # métricas de estado para o log
    # ------------------------------------------------------------------
    def _n_components_uf(self) -> int:
        """Conta componentes não-triviais do UF (cobrindo vértices com degree>=1).
        Vértices isolados (degree=0) contam como uma componente cada.
        Retorna o nº total = isolados + componentes não-triviais."""
        deg = self.detector.degree
        non_iso_roots = set()
        n_iso = 0
        for v in range(self.V):
            if deg[v] == 0:
                n_iso += 1
            else:
                non_iso_roots.add(self.detector.find(v))
        return n_iso + len(non_iso_roots)

    def _max_degree_free(self) -> int:
        s = self.estado
        mat = self.v2e_mat
        vmask = mat >= 0
        idx = np.where(vmask, mat, 0)
        vals = s[idx]
        vals = np.where(vmask, vals, 99)
        n_zero = (vals == 0).sum(axis=1)
        n_free_v = (vals == -1).sum(axis=1)
        livre_eff = self.grau_arr - n_zero  # arestas não-zero remanescentes
        livre_eff = np.where(n_free_v > 0, livre_eff, -1)
        return int(livre_eff.max())

    def _log_event(self, outcome: str, L_current: int,
                   L_parent: int, score_pressure: float) -> None:
        n_fixadas = int((self.estado != FREE).sum())
        self.event_log.append({
            "event_id": len(self.event_log),
            "n": self.n,
            "depth": n_fixadas / self.E,
            "n_components": self._n_components_uf(),
            "max_degree_free": self._max_degree_free(),
            "L_current": L_current,
            "L_parent": L_parent,
            "score_pressure": (
                float(score_pressure)
                if not (isinstance(score_pressure, float)
                        and math.isnan(score_pressure))
                else None
            ),
            "outcome": outcome,
        })

    # ------------------------------------------------------------------
    # branch com logging
    # ------------------------------------------------------------------
    def branch(self, L_parent: int = -1):
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

        # propagação R2
        try:
            self.propagar()
        except ValueError:
            self.podas["R2"] += 1
            self._log_event("R2", L_current=-1, L_parent=L_parent,
                            score_pressure=float("nan"))
            self.estado = snap_est
            self.last_seen = snap_seen
            self.detector.rollback(cp)
            return

        # detector incremental
        det_status = self.sync_det()
        if det_status != OK:
            outcome = "SUBTOUR_EARLY" if det_status == SUBTOUR else "CONTRADICTION_DET"
            self._log_event(outcome, L_current=-1, L_parent=L_parent,
                            score_pressure=float("nan"))
            self.estado = snap_est
            self.last_seen = snap_seen
            self.detector.rollback(cp)
            return

        # escolha de variável
        e, v_alvo, pressao_v = self.escolher_var()
        if e == -1:
            # folha
            self.folhas += 1
            is_tour = False
            if (self.estado == 1).sum() == self.V:
                self.fatores2 += 1
                if self.conexo():
                    self.tours += 1
                    is_tour = True
                    if self.tours >= self.alvo:
                        self.parar = True
            outcome = "TOUR" if is_tour else "LEAF_NON_TOUR"
            self._log_event(outcome, L_current=-1, L_parent=L_parent,
                            score_pressure=float("nan"))
            self.estado = snap_est
            self.last_seen = snap_seen
            self.detector.rollback(cp)
            return

        # nó interno — registra antes de recursar
        L_cur = int(self.level_arr[v_alvo])
        self._log_event("RECURSE_INTERNAL", L_current=L_cur,
                        L_parent=L_parent, score_pressure=pressao_v)

        for val in (1, 0):
            if self.parar:
                break
            inner_est = self.estado.copy()
            inner_seen = self.last_seen.copy()
            inner_cp = self.detector.checkpoint()
            self.estado[e] = val
            self.branch(L_parent=L_cur)
            self.estado = inner_est
            self.last_seen = inner_seen
            self.detector.rollback(inner_cp)

        self.estado = snap_est
        self.last_seen = snap_seen
        self.detector.rollback(cp)

    def executar(self) -> dict:
        self.t0 = time.perf_counter()
        try:
            self.propagar()
        except ValueError as exc:
            return {"erro": f"propagação inicial: {exc}"}
        self.sync_det()
        self.branch()
        dt = time.perf_counter() - self.t0
        n_evt = len(self.event_log)
        outc_count: dict[str, int] = {}
        for ev in self.event_log:
            outc_count[ev["outcome"]] = outc_count.get(ev["outcome"], 0) + 1
        mu = outc_count.get("TOUR", 0) / max(1, n_evt)
        return {
            "n": self.n, "V": self.V, "E": self.E,
            "beta1": self.E - self.V + 1,
            "n_mand": len(self.mand_idx),
            "alvo": self.alvo,
            "tempo_s": dt,
            "n_nos": self.nos, "n_tours": self.tours,
            "n_2fatores": self.fatores2, "n_folhas": self.folhas,
            "nos_por_tour": self.nos / max(1, self.tours),
            "podas": dict(self.podas),
            "n_eventos": n_evt,
            "outcome_count": outc_count,
            "mu": mu,
            "parou_por_timeout": dt > self.timeout * 0.99,
        }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ns", type=str, default="6,8,10,12,14")
    parser.add_argument("--alvo", type=int, default=500)
    parser.add_argument("--timeout", type=float, default=600.0)
    parser.add_argument("--out-dir", type=str, default="data")
    args = parser.parse_args()

    out_dir = ROOT / args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    resumo = []
    for n in [int(x) for x in args.ns.split(",")]:
        print(f"\n=== n={n} (alvo={args.alvo}, timeout={args.timeout}s) ===",
              flush=True)
        bt = LoggedMinimalV2(n=n, alvo=args.alvo, timeout=args.timeout)
        print(
            f"  E={bt.E}  V={bt.V}  β₁={bt.E-bt.V+1}  mand={len(bt.mand_idx)}",
            flush=True,
        )
        r = bt.executar()
        if "erro" in r:
            print(f"  ERRO: {r['erro']}")
            resumo.append(r)
            continue
        print(
            f"  → eventos={r['n_eventos']:6d}  tours={r['n_tours']:3d}  "
            f"μ(n)={r['mu']:.4f}  nós/tour={r['nos_por_tour']:6.2f}  "
            f"t={r['tempo_s']:6.2f}s",
            flush=True,
        )
        print(f"  outcome_count: {r['outcome_count']}", flush=True)

        log_path = out_dir / f"event_log_n{n}.json"
        log_path.write_text(json.dumps(bt.event_log))
        size_mb = log_path.stat().st_size / 1024 / 1024
        print(f"  log salvo: {log_path.name} ({size_mb:.2f} MB)", flush=True)
        resumo.append(r)

    resumo_path = out_dir / "summary_T1.json"
    resumo_path.write_text(json.dumps(resumo, indent=2))
    print(f"\nResumo salvo em {resumo_path}")

    print("\n=== TABELA μ(n) ===")
    print(f"{'n':>3} {'eventos':>9} {'tours':>6} {'μ(n)':>8} {'nós/tour':>9}")
    for r in resumo:
        if "erro" in r:
            continue
        print(f"{r['n']:>3} {r['n_eventos']:>9} {r['n_tours']:>6} "
              f"{r['mu']:>8.4f} {r['nos_por_tour']:>9.2f}")


if __name__ == "__main__":
    main()
