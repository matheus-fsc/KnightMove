"""
T5 — Escalonamento via "v2 mínimo": grafo do cavalo + apenas R2 + detector.

Diferença para `backtracking_v2.py`: NÃO depende de amostras Z3/pairs/xor.
Constrói o grafo n×n do zero, fixa as 8 mandatórias dos cantos (R1
estruturalmente trivial), e roda v2 (R2 + detector incremental).

Como R3/R6 são removidos, espera-se nós/tour MAIORES que no benchmark
completo do 10×10 — mas a razão será mais limpa, pois capta o efeito
do DETECTOR isoladamente.

Tamanhos testados: n ∈ {6, 8, 10, 12}.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from collections import deque
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from subtour_detector import (
    SubtourDetectorCopy, OK, SUBTOUR, COMPLETE_TOUR, CONTRADICTION,
)


FREE = -1


def build_knight_graph(n: int):
    V = n * n
    pos = [(r, c) for r in range(n) for c in range(n)]
    moves = [(1, 2), (2, 1), (-1, 2), (-2, 1), (1, -2), (2, -1), (-1, -2), (-2, -1)]
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


class MinimalV2:
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
            try:
                if self.propagar_R2() == 0:
                    return
            except ValueError as exc:
                raise

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

    def escolher_var(self) -> int:
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
            return -1
        v_alvo = int(np.argmax(pressao))
        row = mat[v_alvo, : self.v2e_len[v_alvo]]
        for e in row:
            if self.estado[e] == FREE:
                return int(e)
        return -1

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

    def executar(self) -> dict:
        self.t0 = time.perf_counter()
        try:
            self.propagar()
        except ValueError as exc:
            return {"erro": f"propagação inicial: {exc}"}
        self.sync_det()
        self.branch()
        dt = time.perf_counter() - self.t0
        return {
            "n": self.n, "V": self.V, "E": self.E,
            "beta1": self.E - self.V + 1,
            "n_mand": len(self.mand_idx),
            "n_free_inicial": int((self.estado == FREE).sum()),
            "alvo": self.alvo,
            "tempo_s": dt,
            "n_nos": self.nos, "n_tours": self.tours,
            "n_2fatores": self.fatores2, "n_folhas": self.folhas,
            "nos_por_tour": self.nos / max(1, self.tours),
            "razao_2fat_tour": self.fatores2 / max(1, self.tours),
            "podas": dict(self.podas),
            "parou_por_timeout": dt > self.timeout * 0.99,
        }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ns", type=str, default="6,8,10,12")
    parser.add_argument("--alvo", type=int, default=50)
    parser.add_argument("--timeout", type=float, default=900.0)
    args = parser.parse_args()

    out = []
    for n in [int(x) for x in args.ns.split(",")]:
        print(f"\n=== n={n} (alvo={args.alvo}, timeout={args.timeout}s) ===")
        # contar n_free inicial após propagação para reportar
        v2 = MinimalV2(n=n, alvo=args.alvo, timeout=args.timeout)
        # snapshot n_free antes do branch
        v2._n_free_inicial_pre = int((v2.estado == FREE).sum())
        print(
            f"  E={v2.E}  V={v2.V}  β₁={v2.E-v2.V+1}  "
            f"mand={len(v2.mand_idx)}  n_free_inicial={v2._n_free_inicial_pre}"
        )
        r = v2.executar()
        if "erro" in r:
            print(f"  ERRO: {r['erro']}")
            out.append(r)
            continue
        print(
            f"  → tours={r['n_tours']:3d}  2fat={r['n_2fatores']:5d}  "
            f"nós={r['n_nos']:7d}  nós/tour={r['nos_por_tour']:8.1f}  "
            f"razão={r['razao_2fat_tour']:.2f}×  t={r['tempo_s']:7.2f}s"
        )
        print(f"  podas: {r['podas']}")
        out.append(r)

    out_path = ROOT / "data" / "scaling_minimal_v2.json"
    out_path.write_text(json.dumps(out, indent=2))
    print(f"\nSalvo em {out_path}")


if __name__ == "__main__":
    main()
