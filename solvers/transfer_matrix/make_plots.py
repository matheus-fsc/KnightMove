"""
Gera os 3 plots:
  data/plots/spectral_gap.png
  data/plots/eigenvector_vs_finf.png
  data/plots/tour_count_scaling.png
"""

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

THIS = Path(__file__).parent
DATA = THIS / "data"
PLOTS = DATA / "plots"
PLOTS.mkdir(exist_ok=True)


def plot_spectral_gap(ns):
    lam1s, lam2s, gaps = [], [], []
    for n in ns:
        path = DATA / f"eigenvalues_n{n}.json"
        if not path.exists():
            continue
        with open(path) as f:
            e = json.load(f)
        lam1s.append((n, e["eigenvalues_abs"][0]))
        lam2s.append((n, e["eigenvalues_abs"][1]))
        gaps.append((n, e["gap_abs"]))

    if not lam1s:
        print("Nenhum eigenvalue file ainda.")
        return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))

    xs1 = [r[0] for r in lam1s]
    ys1 = [r[1] for r in lam1s]
    ys2 = [r[1] for r in lam2s]
    ax1.plot(xs1, ys1, "o-", label="$|\\lambda_1|$", color="C0")
    ax1.plot(xs1, ys2, "s-", label="$|\\lambda_2|$", color="C1")
    ax1.set_yscale("log")
    ax1.set_xlabel("n")
    ax1.set_ylabel("autovalor (|.|)")
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_title("Autovalores dominantes de T_bulk(n)")

    xs2 = [r[0] for r in gaps]
    ys_gap = [r[1] for r in gaps]
    ax2.plot(xs2, ys_gap, "o-", color="C2")
    ax2.set_xlabel("n")
    ax2.set_ylabel("$\\Delta = |\\lambda_1| - |\\lambda_2|$")
    ax2.grid(True, alpha=0.3)
    ax2.set_title("Gap espectral")

    plt.tight_layout()
    out = PLOTS / "spectral_gap.png"
    plt.savefig(out, dpi=130)
    plt.close()
    print(f"  salvo: {out}")


def plot_eigenvector_vs_finf(ns):
    # Carrega marginais por edge shape + nivel
    fig, ax = plt.subplots(1, 1, figsize=(7, 6))

    finf_path = THIS.parent / "incremental_subtour" / "data" / "edge_freq_scaling.json"
    with open(finf_path) as f:
        finf_data = json.load(f)
    finf_by_n = {r["n"]: {int(L): st["freq_mean"] for L, st in r["por_nivel"].items()}
                 for r in finf_data["resultados"]}

    colors = {0: "C0", 1: "C1", 2: "C2", 3: "C3", 4: "C4", 5: "C5"}
    markers = {6: "o", 8: "s", 10: "D"}

    for n in ns:
        path = DATA / f"eigenvector_marginals_n{n}.json"
        if not path.exists():
            continue
        with open(path) as f:
            d = json.load(f)
        finf_n = finf_by_n.get(n, {})
        for L_data in d["by_level_all"]:
            L = L_data["L"]
            eig_mean = L_data["eig_mean"]
            finf = finf_n.get(L)
            if finf is None:
                continue
            ax.scatter(
                finf, eig_mean,
                color=colors.get(L, "gray"),
                marker=markers.get(n, "x"),
                s=80,
                label=f"n={n}, L={L}"
            )

    lim = [0, 0.7]
    ax.plot(lim, lim, "k--", alpha=0.5, label="y=x")
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel("$f_\\infty(L)$ empirico (tours, incremental_subtour)")
    ax.set_ylabel("marginal autovetor v_1 (2-fatores, bulk)")
    ax.set_title("Marginal de aresta: autovetor vs $f_\\infty$")
    ax.legend(fontsize=8, loc="lower right")
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    out = PLOTS / "eigenvector_vs_finf.png"
    plt.savefig(out, dpi=130)
    plt.close()
    print(f"  salvo: {out}")


def plot_tour_count_scaling():
    path = DATA / "tour_count_estimates.json"
    if not path.exists():
        return
    with open(path) as f:
        d = json.load(f)
    rows = d.get("rows", [])
    if not rows:
        return

    fig, ax = plt.subplots(1, 1, figsize=(7, 5))

    ns_known = [r["n"] for r in rows if r["count"] is not None]
    counts_known = [r["count"] for r in rows if r["count"] is not None]
    if ns_known:
        ax.plot(ns_known, counts_known, "o", color="C0", markersize=10,
                label="contagem exata (DP)")
    ax.set_yscale("log")
    ax.set_xlabel("n")
    ax.set_ylabel("# 2-fatores")
    ax.set_title("Scaling do numero de 2-fatores")
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    out = PLOTS / "tour_count_scaling.png"
    plt.savefig(out, dpi=130)
    plt.close()
    print(f"  salvo: {out}")


if __name__ == "__main__":
    ns = [int(x) for x in sys.argv[1:]] if len(sys.argv) > 1 else [6]
    print("Gerando plots...")
    plot_spectral_gap(ns)
    plot_eigenvector_vs_finf(ns)
    plot_tour_count_scaling()
