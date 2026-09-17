"""Gera data/plots/scaling_theory.png — comparação theory vs v2 amostrado."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
INC = ROOT.parent / "incremental_subtour"


def carregar_theory():
    p = ROOT / "data" / "benchmark_theory.json"
    res = json.loads(p.read_text())
    pontos = {}
    for r in res:
        if "erro" in r:
            continue
        pontos[r["n"]] = {
            "nos_por_tour": r["nos_por_tour"],
            "t_total": r["tempo_s"],
        }
    return pontos


def carregar_v2():
    """Reúne n∈{6..12} de scaling_minimal_v2.json + n=14 de scaling_n14.json."""
    p1 = INC / "data" / "scaling_minimal_v2.json"
    res1 = json.loads(p1.read_text())
    pontos = {}
    for r in res1:
        if "erro" in r:
            continue
        pontos[r["n"]] = {
            "nos_por_tour": r["nos_por_tour"],
            "t_total": r["tempo_s"],
        }
    p2 = INC / "data" / "scaling_n14.json"
    if p2.exists():
        r14 = json.loads(p2.read_text())
        pontos[r14["n"]] = {
            "nos_por_tour": r14["nodes_per_tour"],
            "t_total": r14["t_total"],
        }
    return pontos


def main():
    theory = carregar_theory()
    v2 = carregar_v2()
    ns = sorted(set(theory) | set(v2))

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # painel 1: nós/tour
    ax = axes[0]
    xs_t = [n for n in ns if n in theory]
    ys_t = [theory[n]["nos_por_tour"] for n in xs_t]
    xs_v = [n for n in ns if n in v2]
    ys_v = [v2[n]["nos_por_tour"] for n in xs_v]
    ax.plot(xs_v, ys_v, "o-", lw=2, color="#e07a3a", label="v2 amostrado",
            markersize=8)
    ax.plot(xs_t, ys_t, "s-", lw=2, color="#1f77b4", label="theory (f∞)",
            markersize=8)
    for x, y in zip(xs_t, ys_t):
        ax.annotate(f"{y:.2f}", (x, y), xytext=(0, -16), textcoords="offset points",
                    ha="center", fontsize=9, color="#1f4e79")
    for x, y in zip(xs_v, ys_v):
        ax.annotate(f"{y:.2f}", (x, y), xytext=(0, 8), textcoords="offset points",
                    ha="center", fontsize=9, color="#7a3e15")
    ax.set_xticks(ns)
    ax.set_xlabel("n (lado do tabuleiro)")
    ax.set_ylabel("nós explorados por tour")
    ax.set_title("Eficiência: nós/tour vs n  (K=500)")
    ax.set_ylim(2.5, 6.5)
    ax.grid(linestyle=":", alpha=0.4)
    ax.legend(loc="upper right", fontsize=10)

    # painel 2: t_total
    ax2 = axes[1]
    xs_t = [n for n in ns if n in theory]
    yt_t = [theory[n]["t_total"] for n in xs_t]
    xs_v = [n for n in ns if n in v2]
    yt_v = [v2[n]["t_total"] for n in xs_v]
    ax2.plot(xs_v, yt_v, "o-", lw=2, color="#e07a3a", label="v2 amostrado",
             markersize=8)
    ax2.plot(xs_t, yt_t, "s-", lw=2, color="#1f77b4", label="theory (f∞)",
             markersize=8)
    ax2.set_xticks(ns)
    ax2.set_xlabel("n")
    ax2.set_ylabel("t_total (s)")
    ax2.set_title("Tempo de busca vs n  (K=500)")
    ax2.grid(linestyle=":", alpha=0.4)
    ax2.legend(loc="upper left", fontsize=10)

    fig.suptitle(
        "Heurística de fase local f∞(L) vs v2 amostrado  "
        "[sem amostragem prévia]",
        fontsize=12, y=1.02,
    )
    plt.tight_layout()
    out = ROOT / "data" / "plots" / "scaling_theory.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"plot salvo: {out}")

    print("\n=== Tabela final ===")
    print(f"{'n':>3s}  {'theory n/t':>11s}  {'v2 n/t':>9s}  {'razão':>8s}  "
          f"{'t_theory':>9s}  {'t_v2':>9s}")
    for n in ns:
        nt = theory.get(n, {}).get("nos_por_tour")
        nv = v2.get(n, {}).get("nos_por_tour")
        razao = (nt / nv) if (nt and nv) else None
        tt = theory.get(n, {}).get("t_total")
        tv = v2.get(n, {}).get("t_total")
        razao_s = f"{razao:.2f}×" if razao else "  n/a"
        nt_s = f"{nt:8.2f}" if nt else "    n/a"
        nv_s = f"{nv:6.2f}" if nv else "  n/a"
        tt_s = f"{tt:7.2f}s" if tt else "    n/a"
        tv_s = f"{tv:7.2f}s" if tv else "    n/a"
        print(f"{n:3d}  {nt_s}   {nv_s}   {razao_s:>8s}  {tt_s}  {tv_s}")


if __name__ == "__main__":
    main()
