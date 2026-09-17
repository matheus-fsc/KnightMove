"""
Benchmark completo Tabela 9 expandida — 10×10, K=2000.
Inclui backtracking ingênuo (sem otimizações) como baseline absoluto.
"""
from __future__ import annotations
import json, sys, time, random
from pathlib import Path
from statistics import median
from collections import deque

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "incremental_subtour"))
sys.path.insert(0, str(ROOT / "residual_search_10x10"))
sys.path.insert(0, str(ROOT / "local_phase_heuristic"))

N = 10
K = 999_999_999  # sem limite de tours — controlado por tempo
SEED = 42
T_PER_METHOD = 300  # 5 min por método (total ~35 min com 7 métodos)


# ====================================================================
# MÉTODO 0: Backtracking ingênuo — sem NENHUMA otimização
# ====================================================================
class NaiveBacktracker:
    """
    Backtracking puro sobre arestas do grafo do cavalo.
    - Sem propagação R2 (não força arestas quando grau residual = 2)
    - Sem Union-Find (não detecta sub-ciclos prematuros)
    - Sem pressão de vértice (escolhe arestas em ordem fixa de índice)
    - Sem heurística f∞
    - Único corte: grau > 2 ou grau efetivo < 2 (viabilidade mínima)
    """
    def __init__(self, n, alvo, timeout):
        self.n = n
        self.V = n * n
        self.alvo = alvo
        self.timeout = timeout
        self.edges_uv = self._build_edges(n)
        self.E = len(self.edges_uv)
        self.v2e = [[] for _ in range(self.V)]
        for ei, (u, v) in enumerate(self.edges_uv):
            self.v2e[u].append(ei)
            self.v2e[v].append(ei)
        self.grau = [len(self.v2e[v]) for v in range(self.V)]
        self.estado = [-1] * self.E  # -1=livre, 0=excluída, 1=ativa
        self.deg_ativo = [0] * self.V
        self.deg_excl = [0] * self.V
        self.nos = 0
        self.tours = 0
        self.folhas = 0
        self.parar = False
        self.t0 = 0.0

    def _build_edges(self, n):
        moves = [(1,2),(2,1),(-1,2),(-2,1),(1,-2),(2,-1),(-1,-2),(-2,-1)]
        edges = set()
        for i in range(n*n):
            r, c = i // n, i % n
            for dr, dc in moves:
                nr, nc = r+dr, c+dc
                if 0 <= nr < n and 0 <= nc < n:
                    j = nr*n + nc
                    edges.add((min(i,j), max(i,j)))
        return sorted(edges)

    def _viavel_rapido(self):
        for v in range(self.V):
            if self.deg_ativo[v] > 2:
                return False
            if self.grau[v] - self.deg_excl[v] < 2:
                return False
        return True

    def _conexo(self):
        ativos = [ei for ei in range(self.E) if self.estado[ei] == 1]
        if len(ativos) != self.V:
            return False
        adj = [[] for _ in range(self.V)]
        for ei in ativos:
            u, v = self.edges_uv[ei]
            adj[u].append(v)
            adj[v].append(u)
        seen = {0}
        q = deque([0])
        while q:
            x = q.popleft()
            for y in adj[x]:
                if y not in seen:
                    seen.add(y)
                    q.append(y)
        return len(seen) == self.V

    def _escolher_var(self):
        for ei in range(self.E):
            if self.estado[ei] == -1:
                return ei
        return -1

    def branch(self):
        if self.parar:
            return
        self.nos += 1
        if (self.nos & 0x3FFF) == 0:
            if time.perf_counter() - self.t0 > self.timeout:
                self.parar = True
                return

        if not self._viavel_rapido():
            return

        e = self._escolher_var()
        if e == -1:
            self.folhas += 1
            n_ativas = sum(1 for x in self.estado if x == 1)
            if n_ativas == self.V:
                ok = all(self.deg_ativo[v] == 2 for v in range(self.V))
                if ok and self._conexo():
                    self.tours += 1
                    if self.tours >= self.alvo:
                        self.parar = True
            return

        u, v = self.edges_uv[e]
        for val in (1, 0):
            if self.parar:
                break
            self.estado[e] = val
            if val == 1:
                self.deg_ativo[u] += 1
                self.deg_ativo[v] += 1
            else:
                self.deg_excl[u] += 1
                self.deg_excl[v] += 1
            self.branch()
            if val == 1:
                self.deg_ativo[u] -= 1
                self.deg_ativo[v] -= 1
            else:
                self.deg_excl[u] -= 1
                self.deg_excl[v] -= 1
            self.estado[e] = -1

    def executar(self):
        self.t0 = time.perf_counter()
        self.branch()
        dt = time.perf_counter() - self.t0
        return {
            "tempo_s": dt, "n_tours": self.tours, "n_nos": self.nos,
            "nos_por_tour": self.nos / max(1, self.tours),
            "n_folhas": self.folhas,
            "parou_timeout": self.parar and self.tours < self.alvo,
        }


def run_naive(K, timeout):
    bt = NaiveBacktracker(n=N, alvo=K, timeout=timeout)
    return bt.executar()


def run_bt_v1(K, timeout):
    from propagation_engine_10x10 import carregar_dados
    from backtracking_10x10 import Backtracker10
    dados = carregar_dados(strict_pair_mode="todos")
    bt = Backtracker10(dados, alvo_tours=K, timeout_s=timeout)
    r = bt.executar()
    return {"tempo_s": r["tempo_s"], "n_tours": r["n_tours"], "n_nos": r["n_nos"],
            "nos_por_tour": r["nos_por_tour"]}


def run_bt_v2(K, timeout):
    from propagation_engine_10x10 import carregar_dados
    from backtracking_v2 import BacktrackerV2
    dados = carregar_dados(strict_pair_mode="todos")
    bt = BacktrackerV2(dados, alvo_tours=K, timeout_s=timeout)
    r = bt.executar()
    return {"tempo_s": r["tempo_s"], "n_tours": r["n_tours"], "n_nos": r["n_nos"],
            "nos_por_tour": r["nos_por_tour"]}


def run_bt_uniforme(K, timeout):
    from theory_heuristic import TheoryHeuristicV2
    f_unif = {L: 0.25 for L in range(10)}
    bt = TheoryHeuristicV2(n=N, alvo=K, timeout=timeout, f_inf=f_unif)
    t0 = time.perf_counter()
    r = bt.executar()
    dt = time.perf_counter() - t0
    return {"tempo_s": dt, "n_tours": r["n_tours"], "n_nos": r["n_nos"],
            "nos_por_tour": r["nos_por_tour"]}


def run_bt_theory(K, timeout):
    from theory_heuristic import TheoryHeuristicV2
    bt = TheoryHeuristicV2(n=N, alvo=K, timeout=timeout)
    t0 = time.perf_counter()
    r = bt.executar()
    dt = time.perf_counter() - t0
    return {"tempo_s": dt, "n_tours": r["n_tours"], "n_nos": r["n_nos"],
            "nos_por_tour": r["nos_por_tour"]}


def run_z3(K, timeout, com_mandatory=False):
    from scaling_minimal_v2 import build_knight_graph, mandatory_corner_edges
    from z3 import Bool, Not, Or, PbEq, Solver, Sum, sat

    edges_uv, V, v2e = build_knight_graph(N)
    E = len(edges_uv)
    if com_mandatory:
        mand = mandatory_corner_edges(N, v2e)

    x = [Bool(f"x_{i}") for i in range(E)]
    s = Solver()
    for vtx, es in v2e.items():
        s.add(PbEq([(x[e], 1) for e in es], 2))
    if com_mandatory:
        for e in mand:
            s.add(x[e])

    t0 = time.perf_counter()
    tours = []
    tentativas = 0
    while len(tours) < K:
        if time.perf_counter() - t0 > timeout:
            break
        tentativas += 1
        if s.check() != sat:
            break
        m = s.model()
        ativos = [e for e in range(E) if m[x[e]] is not None and bool(m[x[e]])]
        adj = [[] for _ in range(V)]
        for e in ativos:
            u, v = edges_uv[e]
            adj[u].append(v)
            adj[v].append(u)
        bfs = [0]; seen = {0}; qi = 0
        while qi < len(bfs):
            node = bfs[qi]; qi += 1
            for nb in adj[node]:
                if nb not in seen:
                    seen.add(nb)
                    bfs.append(nb)
        if len(seen) == V:
            tours.append(ativos)
            s.add(Sum([x[e] if e in ativos else Not(x[e]) for e in range(E)]) < E)
        else:
            ac = [e for e in range(E) if edges_uv[e][0] in seen and edges_uv[e][1] in seen
                  and m[x[e]] is not None and bool(m[x[e]])]
            if ac:
                s.add(Or([Not(x[e]) for e in ac]))
    dt = time.perf_counter() - t0
    return {"tempo_s": dt, "n_tours": len(tours), "tentativas": tentativas,
            "nos_por_tour": None}


def run_reps(func, K, timeout, reps, label, **kwargs):
    times = []
    last = None
    for i in range(reps):
        np.random.seed(SEED + i)
        random.seed(SEED + i)
        print(f"    rep {i+1}/{reps}...", end=" ", flush=True)
        r = func(K, timeout, **kwargs)
        t = r["tempo_s"]
        tours = r["n_tours"]
        npt = r.get("nos_por_tour")
        npt_s = f"nos/tour={npt:.2f}" if npt else ""
        timeout_flag = " [TIMEOUT]" if r.get("parou_timeout") else ""
        print(f"t={t:.3f}s  tours={tours}{timeout_flag}  {npt_s}")
        times.append(t)
        last = r
        if r.get("parou_timeout"):
            break
    med = median(times)
    last["t_median"] = med
    last["all_times"] = times
    return last


def main():
    t_global_start = time.perf_counter()
    print(f"=== Benchmark Tabela 9 Completo: 10×10, {T_PER_METHOD}s/método, seed={SEED} ===")
    print(f"=== Cada método roda até {T_PER_METHOD}s — métrica: tours/s ===\n")
    results = {}

    T = T_PER_METHOD

    print(f"[1/7] BT theory (f∞ por nível) — {T}s:")
    results["BT_theory"] = run_reps(run_bt_theory, K, T, 1, "BT theory")

    print(f"\n[2/7] BT uniforme (f=0.25 ∀L) — {T}s:")
    results["BT_uniforme"] = run_reps(run_bt_uniforme, K, T, 1, "BT uniforme")

    print(f"\n[3/7] BT v2 (freq amostrada + UF) — {T}s:")
    results["BT_v2"] = run_reps(run_bt_v2, K, T, 1, "BT v2")

    print(f"\n[4/7] BT v1 (pressão+R2/R3/R6, sem UF) — {T}s:")
    results["BT_v1"] = run_reps(run_bt_v1, K, T, 1, "BT v1")

    print(f"\n[5/7] Z3 puro — {T}s:")
    results["Z3_puro"] = run_reps(run_z3, K, T, 1, "Z3 puro")

    print(f"\n[6/7] Z3 + 8 mandatory — {T}s:")
    results["Z3_mandatory"] = run_reps(
        lambda k, t: run_z3(k, t, com_mandatory=True), K, T, 1, "Z3+mand")

    print(f"\n[7/7] Backtracking ingênuo (sem otimizações) — {T}s:")
    results["BT_naive"] = run_reps(run_naive, K, T, 1, "BT naive")

    dt_global = time.perf_counter() - t_global_start

    # === TABELA FINAL ===
    print("\n" + "=" * 95)
    print(f"{'Método':<35} {'t(s)':>8} {'tours':>8} {'tours/s':>10} {'speedup':>10} {'nos/tour':>10}")
    print("-" * 95)

    z3_rate = results["Z3_puro"]["n_tours"] / results["Z3_puro"]["t_median"]

    order = [
        ("BT_naive", "BT ingênuo (sem otimizações)"),
        ("Z3_mandatory", "Z3 + 8 mandatory"),
        ("Z3_puro", "Z3 puro"),
        ("BT_v1", "BT v1 (pressão+R2/R3/R6)"),
        ("BT_v2", "BT v2 (+UF, freq amostrada)"),
        ("BT_uniforme", "BT uniforme (f=0.25 ∀L)"),
        ("BT_theory", "BT theory (f∞ por nível)"),
    ]

    for key, label in order:
        r = results[key]
        t = r["t_median"]
        tours = r["n_tours"]
        npt = r.get("nos_por_tour")

        rate = tours / t if t > 0 else 0
        sp = rate / z3_rate if z3_rate > 0 else 0

        npt_s = f"{npt:.1f}" if isinstance(npt, (int, float)) and npt else "-"
        print(f"{label:<35} {t:>7.1f}s {tours:>8} {rate:>9.1f}/s {sp:>9.1f}× {npt_s:>10}")

    print(f"\nTempo total do benchmark: {dt_global:.1f}s ({dt_global/60:.1f} min)")

    out = ROOT / "data" / "benchmark_tabela9_K2000.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(results, indent=2, default=str))
    print(f"Salvo em {out}")


if __name__ == "__main__":
    main()
