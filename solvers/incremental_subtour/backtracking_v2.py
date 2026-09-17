"""
Backtracking v2 — backtracking_10x10 + detector incremental de sub-ciclos.

Diferença sobre v1: após cada fix+propagação (R2/R3/R6), sincroniza o
detector com as arestas que viraram 1 e podada IMEDIATAMENTE se o
detector reportar SUBTOUR. v1 só detectava 2-fatores desconexos na FOLHA.

Métricas adicionais (por categoria):
  podas["R2"]            — propagação detectou n_um>2 ou grau efetivo <2
  podas["R3"]            — par estrito violado
  podas["R6"]            — XOR violado
  podas["SUBTOUR_EARLY"] — detector incremental disparou SUBTOUR
  podas["SUBTOUR_LEAF"]  — chegou na folha como 2-fator desconexo (residual)

Por padrão usa SubtourDetectorCopy (mais rápida no micro-bench).
Flag --rollback alterna p/ SubtourDetectorRollback (validação cruzada).
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
sys.path.insert(0, str(ROOT.parent / "residual_search_10x10"))

from propagation_engine_10x10 import (
    FREE,
    PropagationEngine10,
    ROOT as ROOT_RS,
    carregar_dados,
)
from subtour_detector import (
    SubtourDetectorCopy,
    SubtourDetectorRollback,
    OK, SUBTOUR, COMPLETE_TOUR, CONTRADICTION,
)


class BacktrackerV2:
    def __init__(
        self,
        dados,
        alvo_tours: int = 50,
        timeout_s: float = 3600.0,
        usar_rollback: bool = False,
    ):
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

        DetCls = SubtourDetectorRollback if usar_rollback else SubtourDetectorCopy
        self.detector = DetCls(
            edges_uv=[tuple(uv) for uv in dados["edges_uv"]],
            n_vertices=self.eng.n_vertices,
        )
        # sincronizar detector com arestas já fixadas em 1 pela propagação inicial
        for e in np.where(self.eng.estado == 1)[0]:
            r = self.detector.fix(int(e))
            if r in (SUBTOUR, CONTRADICTION):
                raise RuntimeError(
                    f"detector inconsistente na inicialização: e={e} retornou {r}"
                )
        self.last_seen = self.eng.estado.copy()

        # estatísticas
        self.nos = 0
        self.podas = {
            "R2": 0, "R3": 0, "R6": 0,
            "SUBTOUR_EARLY": 0,
            "CONTRADICTION_DET": 0,
            "outras": 0,
        }
        # registro por sub-ciclo detectado (para análise T3)
        self.subtour_log: list[dict] = []
        self.folhas = 0
        self.fatores2 = 0  # 2-fatores que chegaram à folha
        self.subtour_leaf = 0
        self.tours = 0
        self.tours_encontrados: list[list[int]] = []
        self.log_progresso: list[tuple[float, int, int, int]] = []
        self.t0 = 0.0
        self.parar = False
        self.heartbeat_passo = 10000
        self.proximo_heartbeat = self.heartbeat_passo

    # ------------------------------------------------------------------
    # Sincronização detector ↔ estado de arestas
    # ------------------------------------------------------------------
    def sync_detector(self) -> str:
        """Empurra para o detector toda aresta nova com estado==1 desde a
        última sincronização. Retorna OK / SUBTOUR / CONTRADICTION.
        Atualiza self.last_seen ao longo do caminho."""
        est = self.eng.estado
        diff = np.where((est == 1) & (self.last_seen != 1))[0]
        for e in diff:
            r = self.detector.fix(int(e))
            self.last_seen[e] = 1
            if r == SUBTOUR:
                # registrar tamanho do sub-ciclo para T3
                ue, ve = self.detector.edges[int(e)]
                root = self.detector.find(int(ue))
                self.subtour_log.append({
                    "no": self.nos,
                    "depth": int((est != FREE).sum()),
                    "edge_idx": int(e),
                    "subcycle_size": int(self.detector.size[root]),
                })
                self.podas["SUBTOUR_EARLY"] += 1
                return SUBTOUR
            if r == CONTRADICTION:
                self.podas["CONTRADICTION_DET"] += 1
                return CONTRADICTION
        # marcar também as 0s (para coerência do diff futuro)
        self.last_seen = est.copy()
        return OK

    # ------------------------------------------------------------------
    # viabilidade local (cópia da v1)
    # ------------------------------------------------------------------
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
        s = self.eng.estado
        mat = self.eng.v2e_mat
        vmask = mat >= 0
        idx = np.where(vmask, mat, 0)
        vals = s[idx]
        vals = np.where(vmask, vals, 99)
        n_um = (vals == 1).sum(axis=1)
        n_zero = (vals == 0).sum(axis=1)

        n_free_v = (vals == -1).sum(axis=1)
        pressao = n_um.astype(np.int64) * 100 + (
            self.eng.grau_arr - n_zero - 2
        )
        pressao = np.where(n_free_v > 0, pressao, -10**9)
        if pressao.max() <= -10**8:
            return -1
        v_alvo = int(np.argmax(pressao))
        row = mat[v_alvo, : self.eng.v2e_len[v_alvo]]
        livres = [int(e) for e in row if s[e] == FREE]
        if not livres:
            return -1
        score = [(abs(float(self.eng.freq[e]) - 0.5), e) for e in livres]
        score.sort(reverse=True)
        return score[0][1]

    def conexo(self) -> bool:
        s = self.eng.estado
        ativos = np.where(s == 1)[0]
        if len(ativos) != 100:
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
                f"  [{dt:7.1f}s] tours={self.tours:3d}  2fat={self.fatores2:4d}  "
                f"sub_early={self.podas['SUBTOUR_EARLY']:5d}  "
                f"nós={self.nos:6d}  nós/tour={self.nos/max(1,self.tours):7.1f}",
                flush=True,
            )
        if self.tours >= self.alvo:
            self.parar = True

    # ------------------------------------------------------------------
    # branch
    # ------------------------------------------------------------------
    def branch(self) -> None:
        if self.parar:
            return
        self.nos += 1
        if self.nos >= self.proximo_heartbeat:
            dt = time.perf_counter() - self.t0
            n_fixed = int((self.eng.estado != FREE).sum())
            print(
                f"  [hb {dt:6.1f}s] nós={self.nos:7d} tours={self.tours:3d} "
                f"folhas={self.folhas:5d} 2fat={self.fatores2:5d} "
                f"sub_early={self.podas['SUBTOUR_EARLY']:5d} "
                f"fixed={n_fixed:3d}/288",
                flush=True,
            )
            self.proximo_heartbeat += self.heartbeat_passo
            if dt > self.timeout:
                self.parar = True
                return

        # ---- snapshot ----
        snap_estado = self.eng.estado.copy()
        snap_last_seen = self.last_seen.copy()
        det_cp = self.detector.checkpoint()

        # ---- propagação ----
        try:
            self.eng.run_to_fixpoint([2, 3, 6])
        except ValueError as exc:
            tag = (
                "R2" if "R2" in str(exc)
                else "R3" if "R3" in str(exc)
                else "R6" if "R6" in str(exc)
                else "outras"
            )
            self.podas[tag] += 1
            self.eng.estado = snap_estado
            self.last_seen = snap_last_seen
            self.detector.rollback(det_cp)
            return

        # ---- detector incremental ----
        det_status = self.sync_detector()
        if det_status != OK:
            self.eng.estado = snap_estado
            self.last_seen = snap_last_seen
            self.detector.rollback(det_cp)
            return

        if not self.viavel():
            self.podas["R2"] += 1
            self.eng.estado = snap_estado
            self.last_seen = snap_last_seen
            self.detector.rollback(det_cp)
            return

        # ---- escolher próxima ----
        e = self.escolher_variavel()
        if e == -1:
            self.folhas += 1
            if (self.eng.estado == 1).sum() == 100 and self.viavel():
                self.fatores2 += 1
                if self.conexo():
                    self.registrar_tour()
                else:
                    # 2-fator desconexo que escapou ao detector — registrar
                    self.subtour_leaf += 1
            self.eng.estado = snap_estado
            self.last_seen = snap_last_seen
            self.detector.rollback(det_cp)
            return

        f = float(self.eng.freq[e])
        ordem_val = (1, 0) if f >= 0.5 else (0, 1)
        for val in ordem_val:
            if self.parar:
                break
            inner_estado = self.eng.estado.copy()
            inner_last_seen = self.last_seen.copy()
            inner_cp = self.detector.checkpoint()
            self.eng.estado[e] = val
            self.branch()
            self.eng.estado = inner_estado
            self.last_seen = inner_last_seen
            self.detector.rollback(inner_cp)

        self.eng.estado = snap_estado
        self.last_seen = snap_last_seen
        self.detector.rollback(det_cp)

    def executar(self) -> dict:
        self.t0 = time.perf_counter()
        self.branch()
        dt = time.perf_counter() - self.t0
        return {
            "tempo_s": dt,
            "alvo_tours": self.alvo,
            "n_tours": self.tours,
            "n_2fatores": self.fatores2,
            "n_subtour_leaf": self.subtour_leaf,
            "n_folhas": self.folhas,
            "n_nos": self.nos,
            "nos_por_tour": self.nos / max(1, self.tours),
            "nos_por_2fator": self.nos / max(1, self.fatores2),
            "razao_residual_2fat_tour": self.fatores2 / max(1, self.tours),
            "podas": dict(self.podas),
            "log_progresso": self.log_progresso,
            "n_subtours_logados": len(self.subtour_log),
            "subtour_log_amostra": self.subtour_log[:200],  # truncar
            "parou_por_alvo": self.tours >= self.alvo,
            "parou_por_timeout": dt >= self.timeout * 0.99,
        }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--alvo", type=int, default=50)
    parser.add_argument("--timeout", type=float, default=1800.0)
    parser.add_argument(
        "--saida", type=str, default="data/backtracking_v2_results.json"
    )
    parser.add_argument("--rollback", action="store_true",
                        help="usa SubtourDetectorRollback ao invés de Copy")
    parser.add_argument("--mode_pares", choices=["todos", "candidatos"],
                        default="todos")
    args = parser.parse_args()

    dados = carregar_dados(strict_pair_mode=args.mode_pares)
    print(
        f"Pares estritos: {len(dados['pares'])}, XORs: {len(dados['xor'])}, "
        f"mandatory: {len(dados['mandatory_idx'])}"
    )
    bt = BacktrackerV2(
        dados, alvo_tours=args.alvo, timeout_s=args.timeout,
        usar_rollback=args.rollback,
    )
    print(
        f"Iniciando v2: alvo={args.alvo} tours  timeout={args.timeout}s  "
        f"detector={'Rollback' if args.rollback else 'Copy'}"
    )
    res = bt.executar()
    print("\n=== RESULTADO v2 ===")
    print(json.dumps({
        k: v for k, v in res.items()
        if k not in ("subtour_log_amostra", "log_progresso")
    }, indent=2)[:2000])

    out = ROOT / args.saida
    out.parent.mkdir(parents=True, exist_ok=True)
    # também salvar a lista completa de subtours em arquivo separado
    sub_log = {
        "subtour_log_full": bt.subtour_log,
        "n_total": len(bt.subtour_log),
    }
    (out.parent / "subtour_log.json").write_text(json.dumps(sub_log, indent=2))
    out.write_text(json.dumps(res, indent=2))
    print(f"\nSalvo em {out}")


if __name__ == "__main__":
    main()
