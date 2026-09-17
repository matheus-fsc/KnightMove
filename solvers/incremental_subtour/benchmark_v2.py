"""
T4 — Benchmark v2 vs v1 vs Z3.

Métodos:
  A) Z3 puro                — grau-2 + sub-tour elimination via BFS na folha
  B) Z3 + 8 obrigatórias    — A + pré-fixação R1
  C) Backtracking v1         — propagação R1..R6, sem detector incremental
  D) Backtracking v2 (Copy)  — v1 + detector incremental (variante cópia)
  E) Backtracking v2 (Roll)  — v1 + detector incremental (variante rollback)

Carrega resultados pré-calculados para C/D/E (mais rápidos);
mede A/B aqui. Salva o resumo e gera plots.
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

from propagation_engine_10x10 import carregar_dados

try:
    from z3 import Bool, Not, Or, PbEq, Solver, Sum, And, sat
    HAS_Z3 = True
except Exception:
    HAS_Z3 = False


def bfs_tour(ativos, edges_uv, V):
    adj = [[] for _ in range(V)]
    for e in ativos:
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
    return len(seen) == V, seen


def z3_busca(dados, alvo: int, com_mandatory: bool, timeout_s: float) -> dict:
    if not HAS_Z3:
        return {"erro": "z3 não disponível"}
    edges_uv = dados["edges_uv"]
    E = len(edges_uv)
    V = max(max(u, v) for u, v in edges_uv) + 1
    v2e = {}
    for ei, (u, v) in enumerate(edges_uv):
        v2e.setdefault(u, []).append(ei)
        v2e.setdefault(v, []).append(ei)

    x = [Bool(f"x_{i}") for i in range(E)]
    solver = Solver()
    for vtx, es in v2e.items():
        solver.add(PbEq([(x[e], 1) for e in es], 2))
    if com_mandatory:
        for e in dados["mandatory_idx"]:
            solver.add(x[e])

    t0 = time.perf_counter()
    tours = []
    t_first = None
    tentativas = 0
    while len(tours) < alvo:
        if time.perf_counter() - t0 > timeout_s:
            break
        tentativas += 1
        if solver.check() != sat:
            break
        m = solver.model()
        ativos = [e for e in range(E) if m[x[e]] is not None and bool(m[x[e]])]
        is_conex, seen = bfs_tour(ativos, edges_uv, V)
        if is_conex:
            tours.append(ativos)
            if t_first is None:
                t_first = time.perf_counter() - t0
            solver.add(Sum(
                [x[e] if e in ativos else Not(x[e]) for e in range(E)]
            ) < E)
        else:
            ac = [
                e for e in range(E)
                if edges_uv[e][0] in seen and edges_uv[e][1] in seen
                and m[x[e]] is not None and bool(m[x[e]])
            ]
            if ac:
                solver.add(Or([Not(x[e]) for e in ac]))
            else:
                solver.add(Not(And([x[e] for e in ativos])))
    t_total = time.perf_counter() - t0
    return {
        "alvo": alvo, "tours": len(tours),
        "t_first": t_first, "t_total": t_total,
        "tentativas": tentativas, "com_mandatory": com_mandatory,
    }


def carregar_bt(json_path: Path) -> dict:
    if not json_path.exists():
        return {"erro": f"{json_path} ausente"}
    r = json.loads(json_path.read_text())
    t_first = None
    for entry in r.get("log_progresso", []):
        if entry[2] >= 1:
            t_first = entry[0]
            break
    out = {
        "tours": r["n_tours"], "t_first": t_first, "t_total": r["tempo_s"],
        "nos": r["n_nos"], "nos_por_tour": r["nos_por_tour"],
        "podas": r["podas"],
    }
    for k in ("n_2fatores", "n_subtour_leaf", "razao_residual_2fat_tour"):
        if k in r:
            out[k] = r[k]
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--alvo", type=int, default=200)
    parser.add_argument("--timeout", type=float, default=1800.0)
    parser.add_argument("--pular_z3", action="store_true")
    args = parser.parse_args()

    dados = carregar_dados(strict_pair_mode="todos")
    out = {"alvo": args.alvo, "timeout": args.timeout}

    if not args.pular_z3:
        print(f"\n=== A) Z3 puro (K={args.alvo}) ===")
        out["A_z3_puro"] = z3_busca(dados, args.alvo, False, args.timeout)
        print(json.dumps(out["A_z3_puro"], indent=2))

        print(f"\n=== B) Z3 + mandatory (K={args.alvo}) ===")
        out["B_z3_mandatory"] = z3_busca(dados, args.alvo, True, args.timeout)
        print(json.dumps(out["B_z3_mandatory"], indent=2))

    print("\n=== C) Backtracking v1 ===")
    out["C_v1"] = carregar_bt(
        ROOT.parent / "residual_search_10x10" / "data"
        / ("v1_K200.json" if args.alvo >= 200 else "backtracking_results.json")
    )
    print(json.dumps({k: v for k, v in out["C_v1"].items() if k != "tours"}, indent=2))

    print("\n=== D) Backtracking v2 (Copy) ===")
    out["D_v2_copy"] = carregar_bt(ROOT / "data" / "v2_K200.json")
    print(json.dumps({k: v for k, v in out["D_v2_copy"].items() if k != "tours"}, indent=2))

    print("\n=== E) Backtracking v2 (Rollback) ===")
    out["E_v2_rollback"] = carregar_bt(ROOT / "data" / "v2_K200_rollback.json")
    print(json.dumps({k: v for k, v in out["E_v2_rollback"].items() if k != "tours"}, indent=2))

    # speedups
    if "A_z3_puro" in out and out["A_z3_puro"].get("t_total"):
        z3_t = out["A_z3_puro"]["t_total"]
        z3m_t = out["B_z3_mandatory"]["t_total"]
        v1_t = out["C_v1"]["t_total"]
        v2_t = out["D_v2_copy"]["t_total"]
        out["speedups_vs_z3_puro"] = {
            "B_mandatory": z3_t / z3m_t if z3m_t else None,
            "C_v1": z3_t / v1_t,
            "D_v2_copy": z3_t / v2_t,
            "E_v2_rollback": z3_t / out["E_v2_rollback"]["t_total"],
        }
        out["speedups_v1_vs_v2"] = v1_t / v2_t

    (ROOT / "data" / "benchmark_results.json").write_text(json.dumps(out, indent=2))
    print(f"\nSalvo em {ROOT / 'data' / 'benchmark_results.json'}")
    if "speedups_vs_z3_puro" in out:
        print("\n=== SPEEDUPS ===")
        print(json.dumps(out["speedups_vs_z3_puro"], indent=2))
        print(f"v1 → v2: {out['speedups_v1_vs_v2']:.2f}×")

    _make_plots(out)


def _make_plots(out):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib indisponível — pulando plots")
        return

    plots = ROOT / "data" / "plots"
    plots.mkdir(parents=True, exist_ok=True)

    methods = []
    t_totals = []
    t_firsts = []
    for key, label in [
        ("A_z3_puro", "Z3 puro"),
        ("B_z3_mandatory", "Z3 + mand"),
        ("C_v1", "BT v1"),
        ("D_v2_copy", "BT v2 (Copy)"),
        ("E_v2_rollback", "BT v2 (Rb)"),
    ]:
        d = out.get(key, {})
        if d.get("t_total") is not None and "erro" not in d:
            methods.append(label)
            t_totals.append(float(d["t_total"]))
            t_firsts.append(float(d.get("t_first") or 0))

    fig, ax = plt.subplots(figsize=(9, 5))
    xs = np.arange(len(methods))
    w = 0.36
    ax.bar(xs - w/2, t_totals, w, label=f"t_total (K={out['alvo']})",
           color="steelblue")
    ax.bar(xs + w/2, t_firsts, w, label="t_first", color="indianred")
    ax.set_yscale("log")
    ax.set_xticks(xs)
    ax.set_xticklabels(methods, rotation=15)
    ax.set_ylabel("tempo (s) — escala log")
    ax.set_title(f"Speedup comparativo no 10×10  (K={out['alvo']} tours)")
    ax.grid(alpha=0.3, axis="y")
    ax.legend()
    fig.tight_layout()
    fig.savefig(plots / "speedup_v2.png", dpi=120)
    plt.close(fig)

    # poda breakdown v1 vs v2
    v1 = out.get("C_v1", {})
    v2 = out.get("D_v2_copy", {})
    if v1.get("podas") and v2.get("podas"):
        cats_v1 = ["R2", "R3", "R6", "leaf_2fat_disconn", "tour_valido"]
        n_v1 = [
            v1["podas"].get("R2", 0),
            v1["podas"].get("R3", 0),
            v1["podas"].get("R6", 0),
            v1.get("n_2fatores", 0) - v1["tours"] if isinstance(v1.get("n_2fatores"), int) else 0,
            v1["tours"],
        ]
        cats_v2 = ["R2", "R3", "R6", "SUBTOUR_EARLY", "leaf_2fat_disconn", "tour_valido"]
        n_v2 = [
            v2["podas"].get("R2", 0),
            v2["podas"].get("R3", 0),
            v2["podas"].get("R6", 0),
            v2["podas"].get("SUBTOUR_EARLY", 0),
            v2.get("n_subtour_leaf", 0),
            v2["tours"],
        ]

        fig, axes = plt.subplots(1, 2, figsize=(12, 5))

        # barras horizontais (mais legíveis quando há categoria 0)
        axes[0].barh(cats_v1, n_v1, color="indianred")
        axes[0].set_title(f"v1 (sem detector)  total={sum(n_v1)}")
        axes[0].grid(alpha=0.3, axis="x")
        axes[0].set_xscale("symlog")
        for i, v_ in enumerate(n_v1):
            axes[0].text(max(v_, 0.5), i, f" {v_}", va="center")

        axes[1].barh(cats_v2, n_v2, color="steelblue")
        axes[1].set_title(f"v2 (com detector)  total={sum(n_v2)}")
        axes[1].grid(alpha=0.3, axis="x")
        axes[1].set_xscale("symlog")
        for i, v_ in enumerate(n_v2):
            axes[1].text(max(v_, 0.5), i, f" {v_}", va="center")

        fig.suptitle(f"Breakdown por tipo de evento  (K={out['alvo']}, 10×10)")
        fig.tight_layout()
        fig.savefig(plots / "poda_breakdown.png", dpi=120)
        plt.close(fig)

    print(f"Plots salvos em {plots}")


if __name__ == "__main__":
    main()
