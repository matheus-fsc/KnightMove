"""
Gera os dois plots do experimento:
  data/plots/propagation_cascade.png
  data/plots/residual_structure.png
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

from propagation_engine import ROOT, carregar_dados, rotulo_aresta


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
    ax.set_ylabel("nº de arestas (de 80)")
    ax.set_title("Cascata de propagação 6×6 — R1..R6 cumulativos")
    ax.set_ylim(0, n_total + 4)
    for i, (a, b, c) in enumerate(zip(f1, f0, fr)):
        ax.text(i, a + b + c + 0.8, f"FREE={c}", ha="center", fontsize=9)
    ax.legend(loc="lower right")
    ax.grid(axis="y", linestyle=":", alpha=0.4)
    plt.tight_layout()
    out = ROOT / "data" / "plots" / "propagation_cascade.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"salvo: {out}")


def plot_residual():
    dados = carregar_dados()
    edges = dados["edges_uv"]
    freq = dados["freq"]
    n_arestas = len(edges)

    rel = json.loads((ROOT / "data" / "residual_variables.json").read_text())
    fixadas_1 = set(rel["fixadas_1"])

    fig, ax = plt.subplots(figsize=(7.5, 7.5))
    ax.set_xlim(-0.6, 5.6)
    ax.set_ylim(-0.6, 5.6)
    ax.set_aspect("equal")

    # tabuleiro
    for r in range(6):
        for c in range(6):
            cor = "#f0d9b5" if (r + c) % 2 == 0 else "#b58863"
            ax.add_patch(plt.Rectangle((c - 0.5, r - 0.5), 1, 1, color=cor, zorder=0))

    # vértice idx -> (linha,coluna) em coords plot (col=x, linha=y de cima p/ baixo)
    def pos(i):
        r, c = i // 6, i % 6
        # invertido: linha 0 em cima → y=5
        return (c, 5 - r)

    # Desenha as 80 arestas: livres em azul (espessura ∝ freq), fixadas_1 em vermelho
    for ei, (u, v) in enumerate(edges):
        xu, yu = pos(u)
        xv, yv = pos(v)
        if ei in fixadas_1:
            ax.plot([xu, xv], [yu, yv], color="#b22222", lw=3.2, alpha=0.95, zorder=2)
        else:
            f = float(freq[ei])
            lw = 0.6 + 4.0 * f
            alpha = 0.25 + 0.55 * f
            ax.plot([xu, xv], [yu, yv], color="#1f4e79", lw=lw, alpha=alpha, zorder=1)

    # vértices
    for i in range(36):
        x, y = pos(i)
        ax.scatter([x], [y], color="white", edgecolor="black", s=120, zorder=3)
        cols = "ABCDEF"
        r, c = i // 6, i % 6
        ax.annotate(
            f"{cols[c]}{6-r}", (x, y), ha="center", va="center", fontsize=7, zorder=4
        )

    handles = [
        mpatches.Patch(color="#b22222", label="fixadas = 1 (8 cantos)"),
        mpatches.Patch(color="#1f4e79", label="livres (72) — espessura ∝ freq"),
    ]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.04), ncol=2)
    ax.set_title(
        "Estrutura residual 6×6 — 1 componente, 72 livres, ρmax=0.77",
        fontsize=11,
    )
    ax.set_xticks([])
    ax.set_yticks([])
    plt.tight_layout()
    out = ROOT / "data" / "plots" / "residual_structure.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"salvo: {out}")


if __name__ == "__main__":
    plot_cascata()
    plot_residual()
