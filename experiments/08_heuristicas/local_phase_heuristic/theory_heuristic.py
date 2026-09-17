"""
Heurística teórica baseada em fase local f∞(L).

Substitui freq[e] amostrado (v2) por uma tabela teórica de 6 valores:

    f∞(L) = {0: 0.528, 1: 0.193, 2: 0.198, 3: 0.294, 4: 0.247, 5: 0.261}

onde L(e) = min(level(u), level(v)) e level(v) é a distância de v à borda
do tabuleiro n×n.

Reutiliza a estrutura de backtracking de scaling_minimal_v2.MinimalV2:
  - propagação R2 vetorizada (NumPy)
  - Union-Find incremental (SubtourDetectorCopy)
  - escolha de aresta por pressão de vértice + |f - 0.5|

ÚNICA mudança: o vetor freq passa a ser f∞(L(e)) construído do grafo,
sem amostragem prévia. Sem dados externos.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
INC = ROOT.parent / "incremental_subtour"
sys.path.insert(0, str(INC))

from scaling_minimal_v2 import MinimalV2, FREE, build_knight_graph  # noqa: E402
from subtour_detector import OK  # noqa: E402


F_INF_DEFAULT = {
    0: 0.528,
    1: 0.193,
    2: 0.198,
    3: 0.294,
    4: 0.247,
    5: 0.261,
}
F_INF_DEFAULT_FALLBACK = 0.25  # para L > 5 (n ≥ 16)


def level_vertice(i: int, n: int) -> int:
    r, c = i // n, i % n
    return min(r, c, n - 1 - r, n - 1 - c)


def level_aresta(u: int, v: int, n: int) -> int:
    return min(level_vertice(u, n), level_vertice(v, n))


def construir_freq_teorica(
    n: int,
    edges_uv,
    f_inf: dict[int, float] | None = None,
    fallback: float = F_INF_DEFAULT_FALLBACK,
) -> tuple[np.ndarray, np.ndarray]:
    """Retorna (freq, levels) onde freq[e] = f_inf[L(e)]."""
    if f_inf is None:
        f_inf = F_INF_DEFAULT
    E = len(edges_uv)
    levels = np.zeros(E, dtype=np.int32)
    freq = np.zeros(E, dtype=np.float64)
    for ei, (u, v) in enumerate(edges_uv):
        L = level_aresta(u, v, n)
        levels[ei] = L
        freq[ei] = f_inf.get(L, fallback)
    return freq, levels


class TheoryHeuristicV2(MinimalV2):
    """v2 mínimo com escolha de aresta guiada por f∞(L(e))."""

    def __init__(
        self,
        n: int,
        alvo: int,
        timeout: float,
        f_inf: dict[int, float] | None = None,
        fallback: float = F_INF_DEFAULT_FALLBACK,
    ):
        super().__init__(n=n, alvo=alvo, timeout=timeout)
        self.freq, self.levels = construir_freq_teorica(
            n, self.edges_uv, f_inf=f_inf, fallback=fallback
        )
        self._f_inf_used = dict(f_inf) if f_inf is not None else dict(F_INF_DEFAULT)
        self._fallback_used = fallback

    def escolher_var(self) -> int:
        """Idêntico ao v2, mas usando self.freq (teórico) ao desempatar."""
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
        livres = [int(e) for e in row if self.estado[e] == FREE]
        if not livres:
            return -1
        score = [(abs(float(self.freq[e]) - 0.5), e) for e in livres]
        score.sort(reverse=True)
        return score[0][1]

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

        # ordem de tentativa pela fase local teórica
        f = float(self.freq[e])
        ordem = (1, 0) if f >= 0.5 else (0, 1)
        for val in ordem:
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


def rodar(n: int, alvo: int, timeout: float,
          f_inf: dict[int, float] | None = None,
          fallback: float = F_INF_DEFAULT_FALLBACK,
          verbose: bool = True) -> dict:
    bt = TheoryHeuristicV2(
        n=n, alvo=alvo, timeout=timeout, f_inf=f_inf, fallback=fallback
    )
    t_first = None
    # patch para capturar t_first
    parent_tours = [0]

    def hook_branch(orig_branch):
        def _wrapped():
            orig_branch()
            if parent_tours[0] == 0 and bt.tours >= 1:
                parent_tours[0] = 1
        return _wrapped

    if verbose:
        print(
            f"  n={n}  V={bt.V}  E={bt.E}  β₁={bt.E-bt.V+1}  "
            f"mand={len(bt.mand_idx)}", flush=True,
        )
        # mostrar quantas arestas em cada nível
        Lcounts = {}
        for L in bt.levels:
            Lcounts[int(L)] = Lcounts.get(int(L), 0) + 1
        Lstr = ", ".join(f"L={L}:{c}" for L, c in sorted(Lcounts.items()))
        print(f"    arestas por nível → {Lstr}", flush=True)
    t0 = time.perf_counter()
    r = bt.executar()
    dt = time.perf_counter() - t0
    if "erro" in r:
        return r
    r["t_first"] = None  # captura precisa fica para um benchmark dedicado
    r["t_setup"] = 0.0   # sem amostragem
    r["f_inf_used"] = bt._f_inf_used
    return r


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--ns", type=str, default="6,8,10,12,14")
    p.add_argument("--alvo", type=int, default=500)
    p.add_argument("--timeout", type=float, default=600.0)
    p.add_argument("--saida", type=str, default="data/benchmark_theory.json")
    p.add_argument("--validar-6x6", action="store_true",
                   help="rodar n=6 com K=500 e parar se nós/tour > 8")
    args = p.parse_args()

    ns = [int(x) for x in args.ns.split(",")]

    if args.validar_6x6 and 6 in ns:
        print("\n=== Validação n=6 (K=500) ===")
        r6 = rodar(6, alvo=500, timeout=120.0)
        if "erro" in r6:
            print(f"  ERRO: {r6['erro']}")
            return
        print(
            f"  → tours={r6['n_tours']}  2fat={r6['n_2fatores']}  "
            f"nós/tour={r6['nos_por_tour']:.2f}  razão={r6['razao_2fat_tour']:.2f}× "
            f" t={r6['tempo_s']:.2f}s"
        )
        if r6["nos_por_tour"] > 8.0:
            print(f"\nABORT: nós/tour={r6['nos_por_tour']:.2f} > 8.0 no 6×6")
            print("       a heurística teórica ou o mapeamento está incorreto.")
            (ROOT / "data" / "benchmark_theory_aborted.json").write_text(
                json.dumps([r6], indent=2)
            )
            return
        if r6["razao_2fat_tour"] > 1.05:
            print(f"\nABORT: razão 2-fat/tour={r6['razao_2fat_tour']:.3f} > 1.05")
            print("       Union-Find não está funcionando.")
            return
        print(
            f"  ✓ validação 6×6 OK (nós/tour≤8, razão≈1.00×) — prosseguindo"
        )
        ns_rest = [n for n in ns if n != 6]
        resultados = [r6]
    else:
        resultados = []
        ns_rest = ns

    for n in ns_rest:
        print(f"\n=== n={n}  (K_alvo={args.alvo}  timeout={args.timeout}s) ===")
        r = rodar(n, alvo=args.alvo, timeout=args.timeout)
        if "erro" in r:
            print(f"  ERRO: {r['erro']}")
            resultados.append(r)
            continue
        if r["razao_2fat_tour"] > 1.05:
            print(
                f"  ⚠️  razão 2-fat/tour={r['razao_2fat_tour']:.3f} > 1.05 — "
                "bug no Union-Find?"
            )
        print(
            f"  → tours={r['n_tours']:4d}  2fat={r['n_2fatores']:6d}  "
            f"nós={r['n_nos']:8d}  nós/tour={r['nos_por_tour']:7.2f}  "
            f"razão={r['razao_2fat_tour']:.3f}×  t={r['tempo_s']:7.2f}s"
        )
        print(f"  podas: {r['podas']}")
        resultados.append(r)

    out = ROOT / args.saida
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(resultados, indent=2))
    print(f"\nSalvo em {out}")


if __name__ == "__main__":
    main()
