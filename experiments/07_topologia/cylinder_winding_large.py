"""cylinder_winding_large.py — Amostra grande para confirmar a paridade
do winding no cilindro 6×6 (e probe 8×8)."""

# --- bootstrap: torna core/ importavel a partir de qualquer subpasta ---
import sys as _sys, pathlib as _pathlib
for _anc in _pathlib.Path(__file__).resolve().parents:
    if (_anc / "core").is_dir() and (_anc / "pyproject.toml").exists():
        _sys.path.insert(0, str(_anc / "core"))
        break
# --- fim do bootstrap ---

import time
from collections import Counter

import numpy as np

from knight_tours_cylinder import knight_tours, verify_tour
from cylinder_winding import _build_dx_decoder, tour_winding, dedup_tours


def amostra_grande(n, m, K_por_seed, n_seeds, label=""):
    decoder = _build_dx_decoder(m)
    chaves = set()
    paridades_unicos = Counter()
    windings_unicos = Counter()

    t0 = time.time()
    total_brutos = 0
    for s in range(n_seeds):
        tours = knight_tours(n, m, K_por_seed, seed=s)
        total_brutos += len(tours)
        for t in tours:
            V = n * m
            ar = []
            for i in range(V):
                u, w = int(t[i]), int(t[(i + 1) % V])
                ar.append(tuple(sorted([u, w])))
            ar.sort()
            k = tuple(ar)
            if k in chaves:
                continue
            chaves.add(k)
            w, _ = tour_winding(t, n, m, decoder)
            windings_unicos[w] += 1
            paridades_unicos[abs(w) % 2] += 1
    t1 = time.time()

    total_unicos = sum(paridades_unicos.values())
    print(f"\n{'='*60}")
    print(f"  {label}  [{n}x{m}]   {n_seeds} seeds × K={K_por_seed}")
    print(f"{'='*60}")
    print(f"  tempo:           {t1-t0:.2f}s")
    print(f"  brutos:          {total_brutos}")
    print(f"  únicos:          {total_unicos}")
    print(f"  distribuição w:")
    for w in sorted(windings_unicos):
        print(f"    w = {w:+d}  : {windings_unicos[w]:6d}")
    print(f"  paridade par   : {paridades_unicos[0]}")
    print(f"  paridade ímpar : {paridades_unicos[1]}")
    if total_unicos > 0:
        pct0 = 100 * paridades_unicos[0] / total_unicos
        pct1 = 100 * paridades_unicos[1] / total_unicos
        print(f"  %: par={pct0:.2f}%  ímpar={pct1:.2f}%")
        diff = paridades_unicos[0] - paridades_unicos[1]
        # erro padrão se ~ Binomial(n, 0.5)
        sigma = (0.25 * total_unicos) ** 0.5
        print(f"  Δ = {diff}  (σ esperado sob H0=50/50: ±{sigma:.1f}, "
              f"z={diff/sigma:.2f})")


if __name__ == '__main__':
    amostra_grande(6, 6, K_por_seed=1000, n_seeds=20,
                   label="6x6 grande")
    amostra_grande(8, 8, K_por_seed=500, n_seeds=10,
                   label="8x8 sondagem")
    amostra_grande(6, 8, K_por_seed=500, n_seeds=10,
                   label="6x8 sondagem")
    amostra_grande(8, 6, K_por_seed=500, n_seeds=10,
                   label="8x6 sondagem")
