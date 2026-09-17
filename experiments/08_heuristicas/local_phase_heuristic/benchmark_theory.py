"""
Benchmark comparativo no 10×10 (K=200):

  A) Z3 puro                              — baseline SAT
  B) BT v2 amostrado (freq pré-calculada) — melhor anterior
  C) BT theory       (f∞ teórico)         — proposta
  D) BT H_uniforme   (f=0.25 ∀L)          — só pressão de vértice

Métricas: t_first, t_total, nós/tour, razão 2-fat/tour, t_setup.

t_setup mede o custo PRÉVIO necessário para a abordagem funcionar em
uma instância NOVA do problema (n×n nunca visto):
  - Z3:        0 (não precisa de dados externos)
  - v2:        custo de amostragem Z3 para estimar freq[e]
               (referência documentada: ~417s para 5000 amostras no 10×10)
  - theory:    0 (tabela f∞ de 6 valores é fixa)
  - H_uniforme: 0 (só pressão de vértice, sem fase local)

t_total_justo = t_setup + t_busca.
"""
from __future__ import annotations

import json
import sys
import time
from collections import deque
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent
INC = ROOT.parent / "incremental_subtour"
RS10 = ROOT.parent / "residual_search_10x10"
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(INC))
sys.path.insert(0, str(RS10))

from theory_heuristic import (  # noqa: E402
    TheoryHeuristicV2, F_INF_DEFAULT,
)

N = 10
K = 200
TIMEOUT = 900.0

T_SETUP_V2_REF = 417.0  # ~5000 amostras Z3 no 10×10 (documentado em memória do projeto)


def bfs_conexo(arestas_em_um, edges_uv, n_vertices):
    adj = [[] for _ in range(n_vertices)]
    for e in arestas_em_um:
        u, v = edges_uv[e]
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
    return len(seen) == n_vertices, seen


def metodo_A_z3(K: int, timeout: float) -> dict:
    """Z3 puro com grau-2 + sub-tour elim por BFS, sobre o grafo n×n."""
    try:
        from z3 import Bool, Not, Or, PbEq, Solver, Sum, sat
    except Exception as e:
        return {"erro": f"z3 não disponível: {e}"}

    from scaling_minimal_v2 import build_knight_graph
    edges_uv, V, _ = build_knight_graph(N)
    E = len(edges_uv)
    v2e: dict = {}
    for ei, (u, v) in enumerate(edges_uv):
        v2e.setdefault(u, []).append(ei)
        v2e.setdefault(v, []).append(ei)

    x = [Bool(f"x_{i}") for i in range(E)]
    s = Solver()
    for vtx, es in v2e.items():
        s.add(PbEq([(x[e], 1) for e in es], 2))

    t0 = time.perf_counter()
    tours: list[list[int]] = []
    t_first = None
    tentativas = 0
    tempo_z3 = 0.0
    while len(tours) < K:
        if time.perf_counter() - t0 > timeout:
            break
        tentativas += 1
        ts = time.perf_counter()
        if s.check() != sat:
            tempo_z3 += time.perf_counter() - ts
            break
        tempo_z3 += time.perf_counter() - ts
        m = s.model()
        ativos = [e for e in range(E) if m[x[e]] is not None and bool(m[x[e]])]
        conexo, seen = bfs_conexo(ativos, edges_uv, V)
        if conexo:
            tours.append(ativos)
            if t_first is None:
                t_first = time.perf_counter() - t0
            s.add(Sum([x[e] if e in ativos else Not(x[e]) for e in range(E)])
                  < E)
        else:
            arestas_comp = [
                e for e in range(E)
                if edges_uv[e][0] in seen and edges_uv[e][1] in seen
                and m[x[e]] is not None and bool(m[x[e]])
            ]
            if arestas_comp:
                s.add(Or([Not(x[e]) for e in arestas_comp]))
    t_total = time.perf_counter() - t0
    return {
        "metodo": "A_z3_puro",
        "K_alvo": K, "n_tours": len(tours),
        "t_first": t_first, "t_total": t_total,
        "tempo_z3_check_total_s": tempo_z3,
        "tentativas": tentativas,
        "t_setup": 0.0,
        "t_total_justo": t_total,
        "nos_por_tour": None, "razao_2fat_tour": None,
    }


def metodo_BT(label: str, f_inf, custom_freq=None) -> dict:
    """Backtracking teórico (theory ou uniforme)."""
    bt = TheoryHeuristicV2(n=N, alvo=K, timeout=TIMEOUT, f_inf=f_inf)
    if custom_freq is not None:
        bt.freq = custom_freq.astype(np.float64)
    t0 = time.perf_counter()
    r = bt.executar()
    dt = time.perf_counter() - t0
    return {
        "metodo": label,
        "K_alvo": K, "n_tours": r["n_tours"],
        "t_first": None, "t_total": dt,
        "n_nos": r["n_nos"], "nos_por_tour": r["nos_por_tour"],
        "razao_2fat_tour": r["razao_2fat_tour"],
        "podas": r["podas"],
        "t_setup": 0.0,
        "t_total_justo": dt,
    }


def metodo_B_v2_amostrado() -> dict:
    """v2 amostrado: carrega resultado de incremental_subtour/data/v2_K200.json
    e adiciona t_setup ≈ 417s (custo de amostragem Z3 das 5000 amostras
    usadas para estimar freq[e])."""
    p = INC / "data" / "v2_K200.json"
    r = json.loads(p.read_text())
    t_setup = T_SETUP_V2_REF
    return {
        "metodo": "B_v2_amostrado",
        "K_alvo": r["alvo_tours"], "n_tours": r["n_tours"],
        "t_first": (r["log_progresso"][0][0] if r["log_progresso"] else None),
        "t_total": r["tempo_s"],
        "n_nos": r["n_nos"], "nos_por_tour": r["nos_por_tour"],
        "razao_2fat_tour": r["razao_residual_2fat_tour"],
        "podas": r["podas"],
        "t_setup": t_setup,
        "t_setup_obs": "5000 amostras Z3 no 10×10 (~7 min, referência documentada)",
        "t_total_justo": t_setup + r["tempo_s"],
    }


def plot_comparacao(metodos: list[dict], out_path: Path):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    labels_short = {
        "A_z3_puro": "Z3 puro",
        "B_v2_amostrado": "BT v2\n(amostrado)",
        "C_BT_theory": "BT theory\n(f∞ teórico)",
        "D_BT_uniforme": "BT H_uniforme\n(f=0.25)",
    }
    nomes = [labels_short[m["metodo"]] for m in metodos]
    nos = [m.get("nos_por_tour") for m in metodos]
    cores = ["#7a7a7a", "#e07a3a", "#1f77b4", "#a07ac8"]

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # painel 1: nós/tour
    ax = axes[0]
    xs = list(range(len(nomes)))
    bars = ax.bar(xs, [v if v is not None else 0 for v in nos], color=cores)
    ax.set_xticks(xs)
    ax.set_xticklabels(nomes)
    ax.set_ylabel("nós/tour")
    ax.set_title(f"Nós/tour no 10×10 (K={K})")
    for x, v in zip(xs, nos):
        if v is not None:
            ax.text(x, v + 0.1, f"{v:.2f}", ha="center", fontsize=10)
        else:
            ax.text(x, 0.5, "n/a", ha="center", fontsize=10)
    ax.grid(axis="y", linestyle=":", alpha=0.4)

    # painel 2: t_total_justo (com setup)
    ax2 = axes[1]
    tt = [m["t_total_justo"] for m in metodos]
    bars2 = ax2.bar(xs, tt, color=cores)
    ax2.set_xticks(xs)
    ax2.set_xticklabels(nomes)
    ax2.set_ylabel("t_total_justo = t_setup + t_busca (s)")
    ax2.set_title(f"Tempo total honesto (com t_setup)")
    ax2.set_yscale("log")
    for x, v in zip(xs, tt):
        ax2.text(x, v * 1.1, f"{v:.2f}s", ha="center", fontsize=9)
    ax2.grid(axis="y", linestyle=":", alpha=0.4, which="both")

    # rodapé com nota
    fig.text(0.5, -0.02,
             f"Nota: t_setup(B) inclui ~417s de amostragem Z3 prévia "
             f"(5000 tours no 10×10).  Outros métodos têm t_setup=0.",
             ha="center", fontsize=9, style="italic")
    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"plot salvo: {out_path}")


def main():
    print(f"=== Benchmark comparativo no 10×10  K={K} ===\n")
    resultados = []

    print("--- A) Z3 puro ---")
    rA = metodo_A_z3(K=K, timeout=TIMEOUT)
    print(json.dumps({k: v for k, v in rA.items() if k != "podas"}, indent=2)[:400])
    resultados.append(rA)

    print("\n--- B) BT v2 amostrado (carregado) ---")
    rB = metodo_B_v2_amostrado()
    print(json.dumps({k: v for k, v in rB.items() if k != "podas"}, indent=2)[:400])
    resultados.append(rB)

    print("\n--- C) BT theory ---")
    rC = metodo_BT("C_BT_theory", f_inf=None)
    print(json.dumps({k: v for k, v in rC.items() if k != "podas"}, indent=2)[:400])
    resultados.append(rC)

    print("\n--- D) BT H_uniforme ---")
    rD = metodo_BT("D_BT_uniforme", f_inf={L: 0.25 for L in range(0, 10)})
    print(json.dumps({k: v for k, v in rD.items() if k != "podas"}, indent=2)[:400])
    resultados.append(rD)

    print("\n=== Resumo ===")
    print(f"{'método':18s}  {'nós/tour':>10s}  {'t_busca':>10s}  "
          f"{'t_setup':>10s}  {'t_total_justo':>14s}")
    for r in resultados:
        if "erro" in r:
            print(f"  {r.get('metodo', '?'):18s} ERRO: {r['erro']}")
            continue
        n_str = (f"{r.get('nos_por_tour'):8.2f}"
                 if r.get('nos_por_tour') is not None else "      n/a")
        print(f"{r['metodo']:18s}  {n_str}   "
              f"{r['t_total']:9.2f}s  {r['t_setup']:9.2f}s  "
              f"{r['t_total_justo']:13.2f}s")

    out = ROOT / "data" / "benchmark_comparison.json"
    out.write_text(json.dumps(resultados, indent=2))
    print(f"\nSalvo em {out}")

    plot_comparacao(resultados, ROOT / "data" / "plots" / "nodes_per_tour_comparison.png")


if __name__ == "__main__":
    main()
