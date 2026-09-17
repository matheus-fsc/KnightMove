#!/usr/bin/env python3
"""
Visualizações para resultados_8x8_path.json.

Gera (por par com amostras suficientes):
  - heatmap de correlação entre loops livres
  - curva de distribuição de pesos Hamming

E sempre gera:
  - painel de resumo (amostras, hamming média/std, loops obrigatórios/impossíveis)
"""
import json
import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path

INPUT = "resultados_8x8_path.json"
OUT_DIR = Path("plots_8x8_path")


def load(path=INPUT):
    with open(path) as f:
        return json.load(f)


def plot_heatmap(label, assinaturas, out_dir):
    if len(assinaturas) < 2:
        print(f"  [{label}] sem amostras suficientes para heatmap.")
        return

    n_loops = len(assinaturas[0])
    cols = [f"L{i}" for i in range(n_loops)]
    df = pd.DataFrame(assinaturas, columns=cols)

    freq = df.mean()
    livres = freq[(freq > 0.0) & (freq < 1.0)].index.tolist()
    if len(livres) < 2:
        print(f"  [{label}] menos de 2 loops livres — heatmap ignorado.")
        return

    corr = df[livres].corr()
    ordem = corr.abs().mean().sort_values(ascending=False).index
    corr = corr.loc[ordem, ordem]

    mask = np.triu(np.ones_like(corr, dtype=bool), k=1)
    labels_ann = corr.apply(
        lambda col: col.map(lambda v: f"{v:.2f}" if abs(v) >= 0.35 else "")
    )

    n = len(corr)
    fig, ax = plt.subplots(figsize=(max(9, n * 0.8), max(8, n * 0.75)))
    sns.heatmap(
        corr, mask=mask, annot=labels_ann, fmt="",
        cmap="RdBu_r", center=0, vmin=-1, vmax=1,
        square=True, linewidths=0.35, linecolor="#e8e8e8",
        cbar_kws={"shrink": 0.82, "label": "Correlação de Pearson"},
        annot_kws={"fontsize": 9}, ax=ax,
    )
    ax.set_title(f"Correlação loops livres — {label}", fontsize=14, pad=12)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
    fig.tight_layout()

    out = out_dir / f"heatmap_{label}.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  [{label}] heatmap → {out}")


def plot_curva(label, assinaturas, out_dir):
    if len(assinaturas) < 2:
        return

    n_loops = len(assinaturas[0])
    cols = [f"L{i}" for i in range(n_loops)]
    df = pd.DataFrame(assinaturas, columns=cols)
    pesos = df.sum(axis=1)
    curva = pesos.value_counts().sort_index()

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(curva.index, curva.values, marker="o", linewidth=2)
    ax.set_xlabel("Peso Hamming (# loops ativos)")
    ax.set_ylabel("Nº de amostras")
    ax.set_title(f"Distribuição pesos — {label}")
    ax.grid(True, alpha=0.4)
    fig.tight_layout()

    out = out_dir / f"curva_{label}.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  [{label}] curva   → {out}")


def plot_resumo(resultados, out_dir):
    """Painel com n_amostras, hamming média/std, nº loops obrig/impossíveis por par."""
    rows = []
    for r in resultados:
        rows.append({
            "par":        r.get("label", "?"),
            "amostras":   r.get("n_amostras", 0),
            "hamming_μ":  r.get("hamming_media", float("nan")),
            "hamming_σ":  r.get("hamming_std",   float("nan")),
            "obrig":      len(r.get("obrigatorios", [])),
            "impos":      len(r.get("impossiveis",  [])),
        })
    df = pd.DataFrame(rows).set_index("par")

    metricas = ["amostras", "hamming_μ", "hamming_σ", "obrig", "impos"]
    n = len(metricas)
    fig, axes = plt.subplots(1, n, figsize=(4 * n, 5), sharey=False)

    for ax, col in zip(axes, metricas):
        vals = df[col]
        colors = ["#4C72B0" if not np.isnan(v) and v > 0 else "#cccccc" for v in vals]
        ax.barh(df.index, vals.fillna(0), color=colors)
        ax.set_title(col, fontsize=11)
        ax.set_xlabel("")
        ax.invert_yaxis()
        ax.grid(axis="x", alpha=0.3)
        for i, v in enumerate(vals):
            if not np.isnan(v):
                ax.text(v, i, f" {v:.1f}" if isinstance(v, float) else f" {int(v)}",
                        va="center", fontsize=8)

    fig.suptitle("Resumo por par canônico — 8×8 knight path", fontsize=13, y=1.02)
    fig.tight_layout()

    out = out_dir / "resumo_pares.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  resumo → {out}")


def main():
    src = INPUT if len(sys.argv) < 2 else sys.argv[1]
    print(f"Carregando {src}...")
    data = load(src)

    OUT_DIR.mkdir(exist_ok=True)
    sns.set_theme(style="white", context="talk")

    resultados = data.get("resultados", [])
    assinaturas_por_par = data.get("assinaturas_por_par", {})

    print(f"\n--- Heatmaps e curvas por par ---")
    for label, assinaturas in assinaturas_por_par.items():
        plot_heatmap(label, assinaturas, OUT_DIR)
        plot_curva(label, assinaturas, OUT_DIR)

    print(f"\n--- Painel de resumo ---")
    plot_resumo(resultados, OUT_DIR)

    print(f"\nPronto. Gráficos em: {OUT_DIR}/")


if __name__ == "__main__":
    main()
