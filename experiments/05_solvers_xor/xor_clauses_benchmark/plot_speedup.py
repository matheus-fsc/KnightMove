#!/usr/bin/env python3
"""
plot_speedup.py
===============
T4: figura comparativa (6×6 e 10×10) do speedup relativo das 5
configurações do benchmark Z3 com cláusulas XOR de paridade.

Lê:
  data/benchmark_6x6_xor.json
  data/benchmark_10x10_xor.json

Salva:
  data/plots/speedup_comparison.png
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
PLOTS = DATA / "plots"

LABELS = {
    "A": "A: Z3 puro",
    "B": "B: + 8 mandatory",
    "C": "C: + 3 XOR mín.",
    "D": "D: + 28 XOR pares",
    "E": "E: + 3 XOR + NOT(A∧B)",
}

COLORS = {"A": "#bdbdbd", "B": "#6baed6", "C": "#fd8d3c",
          "D": "#74c476", "E": "#9e9ac8"}


def plot_panel(ax, benchmark, title):
    results = benchmark["results"]
    cfgs = benchmark["configs"]
    speed = [results[c].get("speedup_vs_A", 1.0) for c in cfgs]
    times = [results[c]["t_total"] for c in cfgs]

    xs = np.arange(len(cfgs))
    bars = ax.bar(xs, speed, color=[COLORS[c] for c in cfgs],
                  edgecolor="black", linewidth=0.5)

    # destaque Config C
    if "C" in cfgs:
        ic = cfgs.index("C")
        bars[ic].set_edgecolor("red")
        bars[ic].set_linewidth(2.0)

    ax.axhline(1.0, color="gray", linestyle="--", linewidth=1)
    for i, (s, t) in enumerate(zip(speed, times)):
        ax.text(i, s + 0.02, f"{s:.2f}×\n({t:.1f}s)", ha="center",
                va="bottom", fontsize=8)

    ax.set_xticks(xs)
    ax.set_xticklabels([LABELS[c].replace(": ", ":\n", 1) for c in cfgs],
                       fontsize=8)
    ax.set_ylabel("speedup vs Config A")
    ax.set_title(title)
    ax.set_ylim(0, max(speed) * 1.18 + 0.1)
    ax.grid(axis="y", alpha=0.3)


def main():
    p6 = json.loads((DATA / "benchmark_6x6_xor.json").read_text())
    p10 = json.loads((DATA / "benchmark_10x10_xor.json").read_text())

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    plot_panel(axes[0], p6,
               f"6×6  (K={p6['K']}, β₁=45, deficit=3)")
    plot_panel(axes[1], p10,
               f"10×10  (K={p10['K']}, β₁=189, deficit=3)")
    fig.suptitle("Impacto das cláusulas XOR de paridade derivadas de H₁",
                 fontsize=12, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.96])

    PLOTS.mkdir(parents=True, exist_ok=True)
    out = PLOTS / "speedup_comparison.png"
    fig.savefig(out, dpi=140, bbox_inches="tight")
    print(f"Salvo: {out.relative_to(ROOT.parent)}")


if __name__ == "__main__":
    main()
