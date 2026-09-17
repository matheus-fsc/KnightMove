"""
T3 — Análise dos sub-ciclos detectados pelo detector incremental.

Lê o `subtour_log.json` produzido por `backtracking_v2.py` e gera:

  1. Histograma de tamanhos de sub-ciclos (|C|)
  2. Distribuição de profundidade na árvore de busca (n_fixed quando disparou)
  3. Estatísticas: média, mediana, p25/p50/p75/p95 de cada métrica
  4. Razão (sub_early / 2-fatores residuais na folha) — verifica se o
     detector está cobrindo a maior parte dos sub-ciclos
  5. Plot data/plots/subtour_size_hist.png + depth_hist.png

Saída: data/subtour_analysis.json + 2 PNGs.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent


def main():
    log = json.loads((ROOT / "data" / "subtour_log.json").read_text())
    full = log["subtour_log_full"]
    if not full:
        print("Sem sub-ciclos no log — abortando.")
        return

    sizes = np.array([s["subcycle_size"] for s in full])
    depths = np.array([s["depth"] for s in full])
    nodes = np.array([s["no"] for s in full])

    def stats(arr):
        return {
            "n": int(len(arr)),
            "min": int(arr.min()),
            "max": int(arr.max()),
            "mean": float(arr.mean()),
            "median": float(np.median(arr)),
            "p25": float(np.percentile(arr, 25)),
            "p75": float(np.percentile(arr, 75)),
            "p95": float(np.percentile(arr, 95)),
        }

    # buckets de tamanho
    bins_size = [0, 5, 10, 20, 40, 60, 80, 100]
    hist_size, _ = np.histogram(sizes, bins=bins_size)
    hist_size_d = {
        f"{bins_size[i]}-{bins_size[i+1]}": int(hist_size[i])
        for i in range(len(bins_size) - 1)
    }

    bins_depth = [0, 50, 100, 150, 200, 250, 288]
    hist_depth, _ = np.histogram(depths, bins=bins_depth)
    hist_depth_d = {
        f"{bins_depth[i]}-{bins_depth[i+1]}": int(hist_depth[i])
        for i in range(len(bins_depth) - 1)
    }

    res_path = ROOT / "data" / "v2_K200.json"
    if res_path.exists():
        res = json.loads(res_path.read_text())
        razao_residual = res["razao_residual_2fat_tour"]
        n_tours = res["n_tours"]
        n_nos = res["n_nos"]
        n_2fat = res["n_2fatores"]
        n_sub_leaf = res["n_subtour_leaf"]
    else:
        razao_residual = n_tours = n_nos = n_2fat = n_sub_leaf = None

    cobertura = {
        "n_subtour_early": int(len(full)),
        "n_subtour_leaf": int(n_sub_leaf) if n_sub_leaf is not None else None,
        "razao_early_leaf": (
            float(len(full)) / max(1, n_sub_leaf)
            if n_sub_leaf is not None else None
        ),
        "razao_residual_pos_v2": razao_residual,
    }

    # profundidade relativa: (depth ao disparar) / 288 (total)
    rel_depth = depths.astype(float) / 288
    rel_stats = {
        "frac_evitada_media": float(1.0 - rel_depth.mean()),
        "frac_evitada_min": float(1.0 - rel_depth.max()),
        "frac_evitada_max": float(1.0 - rel_depth.min()),
    }

    relatorio = {
        "n_subtours": int(len(full)),
        "tamanho_subciclo": stats(sizes),
        "profundidade_disparo_n_fixed": stats(depths),
        "no_busca_quando_disparou": stats(nodes),
        "histograma_tamanho": hist_size_d,
        "histograma_profundidade_n_fixed": hist_depth_d,
        "fracao_arvore_evitada": rel_stats,
        "cobertura": cobertura,
        "interpretacao": _interpret(sizes, depths, rel_depth, cobertura),
    }

    (ROOT / "data" / "subtour_analysis.json").write_text(
        json.dumps(relatorio, indent=2)
    )
    print(json.dumps(relatorio, indent=2))

    _make_plots(sizes, depths, rel_depth)


def _interpret(sizes, depths, rel_depth, cobertura):
    lines = []
    if sizes.mean() < 20:
        lines.append(
            f"Maioria dos sub-ciclos é PEQUENA (média {sizes.mean():.1f} "
            f"vértices) → detector ataca o problema na raiz."
        )
    elif sizes.mean() > 50:
        lines.append(
            f"Sub-ciclos GRANDES (média {sizes.mean():.1f}) → detector "
            f"está disparando tarde; investigar look-ahead."
        )
    else:
        lines.append(
            f"Sub-ciclos de tamanho médio ({sizes.mean():.1f} vértices)."
        )

    if rel_depth.mean() < 0.5:
        lines.append(
            f"Profundidade média {rel_depth.mean()*100:.0f}% da árvore → "
            f"evitamos {(1-rel_depth.mean())*100:.0f}% do trabalho residual."
        )
    else:
        lines.append(
            f"Profundidade média {rel_depth.mean()*100:.0f}% — detecção "
            f"acontece tarde mas ainda assim antes da folha (288)."
        )

    if cobertura["n_subtour_leaf"] == 0:
        lines.append(
            "COBERTURA 100%: zero 2-fatores desconexos chegando à folha. "
            "Razão residual = 1.0 (vs 17.4 em v1)."
        )
    elif cobertura["n_subtour_leaf"] is not None:
        lines.append(
            f"Cobertura imperfeita: {cobertura['n_subtour_leaf']} 2-fatores "
            f"desconexos ainda escapam — candidato a look-ahead."
        )

    return lines


def _make_plots(sizes, depths, rel_depth):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib não disponível — pulando plots")
        return

    plots = ROOT / "data" / "plots"
    plots.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.hist(sizes, bins=np.arange(0, 110, 5), edgecolor="black", color="steelblue")
    ax.set_xlabel("|C| — tamanho do sub-ciclo detectado (vértices)")
    ax.set_ylabel("Frequência")
    ax.set_title(
        f"Distribuição de tamanhos de sub-ciclos detectados (N={len(sizes)})\n"
        f"média={sizes.mean():.1f}  mediana={np.median(sizes):.0f}  "
        f"max={sizes.max()}"
    )
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(plots / "subtour_size_hist.png", dpi=120)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.hist(depths, bins=np.arange(0, 300, 15), edgecolor="black", color="indianred")
    ax.set_xlabel("Profundidade (n_fixed) ao disparar SUBTOUR_EARLY")
    ax.set_ylabel("Frequência")
    ax.set_title(
        f"Profundidade do disparo na árvore de busca (N={len(depths)})\n"
        f"média={depths.mean():.0f}/288 ({rel_depth.mean()*100:.0f}%) — "
        f"evitamos {(1-rel_depth.mean())*100:.0f}% da árvore"
    )
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(plots / "subtour_depth_hist.png", dpi=120)
    plt.close(fig)
    print(f"Plots salvos em {plots}")


if __name__ == "__main__":
    main()
