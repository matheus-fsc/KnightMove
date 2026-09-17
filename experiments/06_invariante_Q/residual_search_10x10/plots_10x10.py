"""
Gera os três plots do experimento 10×10:
  data/plots/propagation_cascade_10x10.png
  data/plots/residual_structure_10x10.png
  data/plots/benchmark_10x10.png
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

from propagation_engine_10x10 import ROOT, carregar_dados


def plot_cascata():
    rel = json.loads((ROOT / "data" / "propagation_log.json").read_text())
    niveis = rel["niveis"]
    n_total = rel["n_arestas"]

    rotulos = [f"N{n['nivel']}\n{n['descricao'].split(':',1)[0]}" for n in niveis]
    f1 = [n["n_fixadas_1"] for n in niveis]
    f0 = [n["n_fixadas_0"] for n in niveis]
    fr = [n["n_free"] for n in niveis]

    fig, ax = plt.subplots(figsize=(10, 5))
    x = np.arange(len(niveis))
    ax.bar(x, f1, label="fixadas = 1", color="#b22222")
    ax.bar(x, f0, bottom=f1, label="fixadas = 0", color="#444444")
    ax.bar(x, fr, bottom=np.array(f1) + np.array(f0), label="FREE", color="#4682b4")

    ax.set_xticks(x)
    ax.set_xticklabels(rotulos, fontsize=9)
    ax.set_ylabel(f"nº de arestas (de {n_total})")
    ax.set_title("Cascata de propagação 10×10 — R1..R6 cumulativos")
    ax.set_ylim(0, n_total + 12)
    for i, (a, b, c) in enumerate(zip(f1, f0, fr)):
        ax.text(i, a + b + c + 3, f"FREE={c}", ha="center", fontsize=9)
    ax.legend(loc="lower right")
    ax.grid(axis="y", linestyle=":", alpha=0.4)
    plt.tight_layout()
    out = ROOT / "data" / "plots" / "propagation_cascade_10x10.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"salvo: {out}")


def plot_residual():
    dados = carregar_dados(strict_pair_mode="todos")
    edges = dados["edges_uv"]
    freq = dados["freq"]

    rel = json.loads((ROOT / "data" / "residual_variables.json").read_text())
    fixadas_1 = set(rel["fixadas_1"])
    n_free = rel["n_free"]
    rho_max = rel["correlacoes_stats"]["corr_abs_max"]

    fig, ax = plt.subplots(figsize=(9, 9))
    ax.set_xlim(-0.6, 9.6)
    ax.set_ylim(-0.6, 9.6)
    ax.set_aspect("equal")

    for r in range(10):
        for c in range(10):
            cor = "#f0d9b5" if (r + c) % 2 == 0 else "#b58863"
            ax.add_patch(plt.Rectangle((c - 0.5, r - 0.5), 1, 1, color=cor, zorder=0))

    def pos(i):
        r, c = i // 10, i % 10
        return (c, 9 - r)

    for ei, (u, v) in enumerate(edges):
        xu, yu = pos(u)
        xv, yv = pos(v)
        if ei in fixadas_1:
            ax.plot([xu, xv], [yu, yv], color="#b22222", lw=3.0, alpha=0.95, zorder=2)
        else:
            f = float(freq[ei])
            lw = 0.4 + 3.2 * f
            alpha = 0.20 + 0.55 * f
            ax.plot([xu, xv], [yu, yv], color="#1f4e79", lw=lw, alpha=alpha, zorder=1)

    for i in range(100):
        x, y = pos(i)
        ax.scatter([x], [y], color="white", edgecolor="black", s=80, zorder=3)
        cols = "ABCDEFGHIJ"
        r, c = i // 10, i % 10
        ax.annotate(
            f"{cols[c]}{10-r}",
            (x, y),
            ha="center",
            va="center",
            fontsize=5.5,
            zorder=4,
        )

    handles = [
        mpatches.Patch(color="#b22222", label="fixadas = 1 (8 cantos)"),
        mpatches.Patch(color="#1f4e79", label=f"livres ({n_free}) — espessura ∝ freq"),
    ]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.02), ncol=2)
    ax.set_title(
        f"Estrutura residual 10×10 — 1 componente, {n_free} livres, |ρ|max={rho_max:.2f}",
        fontsize=11,
    )
    ax.set_xticks([])
    ax.set_yticks([])
    plt.tight_layout()
    out = ROOT / "data" / "plots" / "residual_structure_10x10.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"salvo: {out}")


def plot_benchmark():
    rel = json.loads((ROOT / "data" / "benchmark_comparison.json").read_text())
    A = rel["A_z3_puro"]
    B = rel["B_z3_mandatory"]
    C = rel["C_backtracking"]

    metodos = ["Z3 puro", "Z3 + mandatory\n(R1 apenas)", "Backtracking\n+ R1..R6"]
    t_first = [A["t_first"], B["t_first"], C["t_first"]]
    t_total = [A["t_total"], B["t_total"], C["t_total"]]

    x = np.arange(len(metodos))
    w = 0.35

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    ax1.bar(x - w / 2, t_first, w, label="t_first (1º tour)", color="#4682b4")
    ax1.bar(x + w / 2, t_total, w, label=f"t_total (K={rel['alvo']} tours)", color="#b22222")
    ax1.set_xticks(x)
    ax1.set_xticklabels(metodos, fontsize=9)
    ax1.set_ylabel("tempo (s)")
    ax1.set_title(f"Tempo p/ encontrar {rel['alvo']} tours fechados — 10×10")
    for i, (a, b) in enumerate(zip(t_first, t_total)):
        ax1.text(i - w / 2, a + 0.05, f"{a:.2f}", ha="center", fontsize=8)
        ax1.text(i + w / 2, b + 0.05, f"{b:.2f}", ha="center", fontsize=8)
    ax1.legend(loc="upper right")
    ax1.grid(axis="y", linestyle=":", alpha=0.4)

    speedup_first = [A["t_first"] / t for t in t_first]
    speedup_total = [A["t_total"] / t for t in t_total]
    ax2.bar(x - w / 2, speedup_first, w, label="speedup t_first", color="#4682b4")
    ax2.bar(x + w / 2, speedup_total, w, label="speedup t_total", color="#b22222")
    ax2.set_xticks(x)
    ax2.set_xticklabels(metodos, fontsize=9)
    ax2.set_ylabel("speedup vs Z3 puro (×)")
    ax2.set_title("Aceleração relativa")
    ax2.axhline(1.0, color="black", linestyle="--", alpha=0.5)
    for i, (a, b) in enumerate(zip(speedup_first, speedup_total)):
        ax2.text(i - w / 2, a + 0.05, f"{a:.2f}×", ha="center", fontsize=8)
        ax2.text(i + w / 2, b + 0.05, f"{b:.2f}×", ha="center", fontsize=8)
    ax2.legend(loc="upper left")
    ax2.grid(axis="y", linestyle=":", alpha=0.4)

    plt.tight_layout()
    out = ROOT / "data" / "plots" / "benchmark_10x10.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"salvo: {out}")


if __name__ == "__main__":
    (ROOT / "data" / "plots").mkdir(parents=True, exist_ok=True)
    plot_cascata()
    plot_residual()
    plot_benchmark()
