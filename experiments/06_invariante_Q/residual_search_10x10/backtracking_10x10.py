"""
Backtracking com propagação R1..R6 para encontrar tours fechados 10×10.

Objetivo: encontrar K tours (default 50) e comparar com Z3.
Métricas: nós explorados, folhas, 2-fatores válidos, tours conexos,
tempo, razão nós/tour.
"""

from __future__ import annotations

import argparse
import json
import time
from collections import deque
from pathlib import Path

import numpy as np

from propagation_engine_10x10 import (
    FREE,
    PropagationEngine10,
    ROOT,
    carregar_dados,
)


class Backtracker10:
    def __init__(self, dados, alvo_tours: int = 50, timeout_s: float = 3600.0):
        self.dados = dados
        self.alvo = alvo_tours
        self.timeout = timeout_s
        self.eng = PropagationEngine10(
            edges_uv=dados["edges_uv"],
            freq=dados["freq"],
            pares=dados["pares"],
            triplas=dados["triplas"],
            quadras=dados["quadras"],
            xor=dados["xor"],
            xor_paridade=dados["xor_paridade_alvo"],
        )
        self.eng.run_to_fixpoint([1, 2, 3, 6])
        self.eng.gravar_log = False
        self.eng.log = []

        freq = self.eng.freq
        score = np.minimum(freq, 1 - freq)
        self.prioridade = np.argsort(score)  # fallback global

        self.nos = 0
        self.podas_R = {"R2": 0, "R3": 0, "R6": 0, "outras": 0}
        self.folhas = 0
        self.fatores2 = 0
        self.tours = 0
        self.tours_encontrados: list[list[int]] = []
        self.t0 = 0.0
        self.log_progresso: list[tuple[float, int, int, int]] = []
        self.parar = False
        # heartbeat para diagnóstico
        self.heartbeat_passo = 10000
        self.proximo_heartbeat = self.heartbeat_passo

    def viavel(self) -> bool:
        s = self.eng.estado
        mat = self.eng.v2e_mat
        vmask = mat >= 0
        idx = np.where(vmask, mat, 0)
        vals = s[idx]
        vals = np.where(vmask, vals, 99)
        n_um = (vals == 1).sum(axis=1)
        n_zero = (vals == 0).sum(axis=1)
        if (n_um > 2).any():
            return False
        if (self.eng.grau_arr - n_zero < 2).any():
            return False
        return True

    def escolher_variavel(self) -> int:
        """Escolha dinâmica: aresta livre incidente ao vértice de maior
        pressão (max(n_um, deg-n_zero-2) — quanto mais saturado, mais
        forte é a propagação ao branchar lá)."""
        s = self.eng.estado
        mat = self.eng.v2e_mat
        vmask = mat >= 0
        idx = np.where(vmask, mat, 0)
        vals = s[idx]
        vals = np.where(vmask, vals, 99)
        n_um = (vals == 1).sum(axis=1)
        n_zero = (vals == 0).sum(axis=1)

        # pressão = n_um (já em 1, prestes a saturar grau-2)
        # vértices com alguma livre incidente
        n_free_v = (vals == -1).sum(axis=1)
        # candidatos: vértices com livres E com pressão
        pressao = n_um.astype(np.int64) * 100 + (
            self.eng.grau_arr - n_zero - 2
        )  # menor (deg-n_zero-2) = mais perto de forçar; somar n_um*100 prioriza n_um
        pressao = np.where(n_free_v > 0, pressao, -10**9)
        if pressao.max() <= -10**8:
            return -1
        v_alvo = int(np.argmax(pressao))
        # aresta livre em v_alvo com freq mais decidida
        row = mat[v_alvo, : self.eng.v2e_len[v_alvo]]
        livres = [int(e) for e in row if s[e] == FREE]
        if not livres:
            return -1
        # entre as livres, escolher a com freq mais distante de 0.5
        score = [(abs(float(self.eng.freq[e]) - 0.5), e) for e in livres]
        score.sort(reverse=True)
        return score[0][1]

    def conexo(self) -> bool:
        s = self.eng.estado
        ativos = np.where(s == 1)[0]
        if len(ativos) != 100:  # |V| = 100 → tour deve ter 100 arestas
            return False
        adj = [[] for _ in range(self.eng.n_vertices)]
        for e in ativos:
            u, v = self.eng.edges_uv[e]
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
        return len(seen) == 100

    def registrar_tour(self):
        estado = self.eng.estado.copy()
        self.tours_encontrados.append([int(i) for i in np.where(estado == 1)[0]])
        self.tours += 1
        if self.tours % 5 == 0 or self.tours == 1:
            dt = time.perf_counter() - self.t0
            self.log_progresso.append(
                (dt, self.nos, self.tours, self.fatores2)
            )
            print(
                f"  [{dt:7.1f}s] tours={self.tours:3d}  2fatores={self.fatores2:4d}  "
                f"nós={self.nos:7d}  nós/tour={self.nos/max(1,self.tours):8.1f}"
            )
        if self.tours >= self.alvo:
            self.parar = True

    def branch(self) -> None:
        if self.parar:
            return
        self.nos += 1
        if self.nos >= self.proximo_heartbeat:
            dt = time.perf_counter() - self.t0
            n_fixed = int((self.eng.estado != FREE).sum())
            print(
                f"  [heartbeat {dt:6.1f}s] nós={self.nos:8d}  tours={self.tours:3d}  "
                f"folhas={self.folhas:5d}  2fatores={self.fatores2:5d}  "
                f"fixed_atual={n_fixed:3d}/288",
                flush=True,
            )
            self.proximo_heartbeat += self.heartbeat_passo
            if dt > self.timeout:
                self.parar = True
                return

        snapshot_estado = self.eng.estado.copy()
        try:
            self.eng.run_to_fixpoint([2, 3, 6])
        except ValueError as exc:
            tag = "R2" if "R2" in str(exc) else "R3" if "R3" in str(exc) else "R6" if "R6" in str(exc) else "outras"
            self.podas_R[tag] += 1
            self.eng.estado = snapshot_estado
            return

        if not self.viavel():
            self.podas_R["R2"] += 1
            self.eng.estado = snapshot_estado
            return

        e = self.escolher_variavel()
        if e == -1:
            self.folhas += 1
            if (self.eng.estado == 1).sum() == 100 and self.viavel():
                self.fatores2 += 1
                if self.conexo():
                    self.registrar_tour()
            self.eng.estado = snapshot_estado
            return

        f = float(self.eng.freq[e])
        ordem_val = (1, 0) if f >= 0.5 else (0, 1)
        for val in ordem_val:
            if self.parar:
                break
            inner = self.eng.estado.copy()
            self.eng.estado[e] = val
            self.branch()
            self.eng.estado = inner

        self.eng.estado = snapshot_estado

    def executar(self) -> dict:
        self.t0 = time.perf_counter()
        self.branch()
        dt = time.perf_counter() - self.t0
        return {
            "tempo_s": dt,
            "alvo_tours": self.alvo,
            "n_tours": self.tours,
            "n_2fatores": self.fatores2,
            "n_folhas": self.folhas,
            "n_nos": self.nos,
            "nos_por_tour": self.nos / max(1, self.tours),
            "nos_por_2fator": self.nos / max(1, self.fatores2),
            "podas": dict(self.podas_R),
            "log_progresso": self.log_progresso,
            "tours": self.tours_encontrados,
            "parou_por_alvo": self.tours >= self.alvo,
            "parou_por_timeout": dt >= self.timeout * 0.99,
        }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--alvo", type=int, default=50)
    parser.add_argument("--timeout", type=float, default=3600.0)
    parser.add_argument("--saida", type=str, default="data/backtracking_results.json")
    parser.add_argument("--mode_pares", choices=["todos", "candidatos"], default="todos")
    args = parser.parse_args()

    dados = carregar_dados(strict_pair_mode=args.mode_pares)
    print(
        f"Pares estritos: {len(dados['pares'])}, XORs: {len(dados['xor'])}, "
        f"mandatory: {len(dados['mandatory_idx'])}"
    )
    bt = Backtracker10(dados, alvo_tours=args.alvo, timeout_s=args.timeout)
    print(
        f"Iniciando: alvo={args.alvo} tours, timeout={args.timeout}s"
    )
    res = bt.executar()
    print("\n=== RESULTADO ===")
    print(json.dumps({k: v for k, v in res.items() if k != "tours"}, indent=2)[:1500])

    out = ROOT / args.saida
    # Não salvar tours inteiros se forem muitos para evitar JSON gigante
    res_save = dict(res)
    out.write_text(json.dumps(res_save, indent=2))
    print(f"\nSalvo em {out}")


if __name__ == "__main__":
    main()
