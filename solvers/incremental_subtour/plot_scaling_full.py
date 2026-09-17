"""Plot escalonamento completo n ∈ {6, 8, 10, 12, 14} — 3 painéis."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent


def coletar():
    base = json.loads((ROOT / "data" / "scaling_minimal_v2.json").read_text())
    n14 = json.loads((ROOT / "data" / "scaling_n14.json").read_text())

    linhas = []
    for r in base:
        linhas.append(
            {
                "n": r["n"],
                "V": r["V"],
                "E": r["E"],
                "beta1": r["beta1"],
                "n_free": r["n_free_inicial"],
                "K": r["n_tours"],
                "t_total": r["tempo_s"],
                "nodes_per_tour": r["nos_por_tour"],
                "ratio": r["razao_2fat_tour"],
            }
        )
    linhas.append(
        {
            "n": n14["n"],
            "V": n14["V"],
            "E": n14["E"],
            "beta1": n14["beta1"],
            "n_free": n14["n_free"],
            "K": n14["K_found"],
            "t_total": n14["t_total"],
            "nodes_per_tour": n14["nodes_per_tour"],
            "ratio": n14["ratio_2fat_tours"],
        }
    )
    linhas.sort(key=lambda r: r["n"])
    return linhas


def main():
    linhas = coletar()
    ns = np.array([r["n"] for r in linhas])
    nodes_per = np.array([r["nodes_per_tour"] for r in linhas])
    t_total = np.array([r["t_total"] for r in linhas])
    ratio = np.array([r["ratio"] for r in linhas])

    fig, axes = plt.subplots(1, 3, figsize=(15, 6.5))

    # ── Painel 1: nós/tour vs n ──
    ax1 = axes[0]
    cores = [
        "#4682b4" if 3.0 <= v <= 7.0 else "#b22222" for v in nodes_per
    ]
    ax1.bar(ns, nodes_per, color=cores, width=1.0, edgecolor="black")
    ax1.axhline(5.0, color="gray", linestyle="--", alpha=0.5, label="banda ~4-5")
    ax1.axhline(3.0, color="lightgray", linestyle=":", alpha=0.6)
    ax1.axhline(7.0, color="lightgray", linestyle=":", alpha=0.6)
    for x, y in zip(ns, nodes_per):
        ax1.text(x, y + 0.15, f"{y:.2f}", ha="center", fontsize=9)
    ax1.set_xticks(ns)
    ax1.set_xlabel("n (lado do tabuleiro)")
    ax1.set_ylabel("nós explorados / tour")
    ax1.set_title("Nós/tour vs n — escala constante")
    ax1.set_ylim(0, max(nodes_per.max() + 1.5, 8))
    ax1.grid(axis="y", linestyle=":", alpha=0.4)
    ax1.legend(loc="upper right")

    # ── Painel 2: t_total vs n ──
    ax2 = axes[1]
    ax2.plot(ns, t_total, "o-", color="#1f4e79", lw=2, markersize=8)
    for x, y in zip(ns, t_total):
        ax2.text(x, y + 0.1, f"{y:.2f}s", ha="center", fontsize=9)
    ax2.set_xticks(ns)
    ax2.set_xlabel("n")
    ax2.set_ylabel("t_total (s) p/ K=500 tours")
    ax2.set_title("Tempo total vs n (alvo K=500)")
    ax2.set_ylim(0, t_total.max() + 1.5)
    ax2.grid(linestyle=":", alpha=0.4)
    # ajuste polinomial leve só para o leitor
    if len(ns) >= 3:
        p = np.polyfit(ns, t_total, 2)
        xfine = np.linspace(ns.min(), ns.max(), 50)
        ax2.plot(xfine, np.polyval(p, xfine), "--", color="gray",
                 alpha=0.5, label=f"fit O(n²): {p[0]:+.3f}n²{p[1]:+.3f}n{p[2]:+.3f}")
        ax2.legend(loc="upper left", fontsize=9)

    # ── Painel 3: razão 2-fat / tours ──
    ax3 = axes[2]
    cores3 = ["#4682b4" if abs(r - 1.0) < 0.05 else "#b22222" for r in ratio]
    ax3.bar(ns, ratio, color=cores3, width=1.0, edgecolor="black")
    ax3.axhline(1.0, color="black", linestyle="--", alpha=0.5,
                label="alvo: 1.00× (detector ótimo)")
    for x, y in zip(ns, ratio):
        ax3.text(x, y + 0.02, f"{y:.2f}×", ha="center", fontsize=9)
    ax3.set_xticks(ns)
    ax3.set_xlabel("n")
    ax3.set_ylabel("razão 2-fatores / tours")
    ax3.set_title("Razão 2-fatores/tours — detector cobre todos sub-ciclos")
    ax3.set_ylim(0, max(ratio.max() * 1.15, 1.2))
    ax3.grid(axis="y", linestyle=":", alpha=0.4)
    ax3.legend(loc="upper right")

    fig.suptitle(
        "Escalonamento incremental_subtour v2 — n ∈ {6, 8, 10, 12, 14}, K=500",
        fontsize=12,
        y=1.02,
    )
    plt.tight_layout()
    out = ROOT / "data" / "plots" / "scaling_full.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"salvo: {out}")

    # imprimir tabela
    print("\n| n  | V   | E   | β₁  | n_free | nós/tour | razão | t(s) |")
    print("|----|-----|-----|-----|--------|----------|-------|------|")
    for r in linhas:
        print(
            f"| {r['n']:>2} | {r['V']:>3} | {r['E']:>3} | {r['beta1']:>3} | "
            f"{r['n_free']:>6} | {r['nodes_per_tour']:>8.2f} | "
            f"{r['ratio']:>5.2f} | {r['t_total']:>4.2f} |"
        )


if __name__ == "__main__":
    main()
