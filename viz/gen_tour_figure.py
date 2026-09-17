#!/usr/bin/env python3
"""Gera figura de um tour hamiltoniano FECHADO do cavalo no 8×8.

Usa knight_tours.py para gerar um tour verificado (verify_tour),
depois desenha com setas direcionais e numeração clara.
"""

# --- bootstrap: torna core/ importavel a partir de qualquer subpasta ---
import sys as _sys, pathlib as _pathlib
for _anc in _pathlib.Path(__file__).resolve().parents:
    if (_anc / "core").is_dir() and (_anc / "pyproject.toml").exists():
        _sys.path.insert(0, str(_anc / "core"))
        break
# --- fim do bootstrap ---
import sys
sys.path.insert(0, "/home/math/Dev/knight_tour")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import knight_tours as kt

N = 8

tour = kt.knight_tours(N, 1, seed=42)[0]
assert kt.verify_tour(tour, N), "Tour inválido!"

path = [divmod(int(v), N) for v in tour]

fig, ax = plt.subplots(1, 1, figsize=(7, 7))

for r in range(N):
    for c in range(N):
        color = "#F0D9B5" if (r + c) % 2 == 0 else "#B58863"
        ax.add_patch(patches.Rectangle((c, N-1-r), 1, 1,
                     facecolor=color, edgecolor="#444", linewidth=0.4))

xs = [c + 0.5 for r, c in path]
ys = [N - 1 - r + 0.5 for r, c in path]

for i in range(len(path)):
    x0, y0 = xs[i], ys[i]
    x1, y1 = xs[(i+1) % len(path)], ys[(i+1) % len(path)]
    color = "#922B21" if i == len(path)-1 else "#1a5276"
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                arrowprops=dict(arrowstyle="->,head_width=0.15,head_length=0.12",
                                color=color, lw=1.1),
                zorder=2)

for idx, (r, c) in enumerate(path):
    x, y = c + 0.5, N - 1 - r + 0.5
    circle = plt.Circle((x, y), 0.28, color="#2c3e50", zorder=3)
    ax.add_patch(circle)
    ax.text(x, y, str(idx + 1), ha="center", va="center",
            fontsize=6.2, fontweight="bold", color="white", zorder=4)

cols = "abcdefgh"
for c in range(N):
    ax.text(c + 0.5, -0.3, cols[c], ha="center", va="center", fontsize=11,
            fontfamily="serif")
for r in range(N):
    ax.text(-0.3, N - 1 - r + 0.5, str(r + 1), ha="center", va="center",
            fontsize=11, fontfamily="serif")

ax.set_xlim(-0.5, N + 0.3)
ax.set_ylim(-0.5, N + 0.3)
ax.set_aspect("equal")
ax.axis("off")

fig.tight_layout()
fig.savefig("knight_tour_8x8.pdf", bbox_inches="tight")
fig.savefig("knight_tour_8x8.png", bbox_inches="tight", dpi=300)
print(f"Tour verificado: {len(path)} casas, fecha ciclo ✓")
print("Figuras salvas: knight_tour_8x8.pdf / .png")
