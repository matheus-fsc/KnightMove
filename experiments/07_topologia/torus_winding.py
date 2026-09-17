"""torus_winding.py — Mede os dois winding numbers (w_x, w_y) e a
distribuição conjunta de paridades nos tours fechados do cavalo no toro.

H₁(T²; ℤ/2) = (ℤ/2)² → há 4 classes (par,par)(par,ímp)(ímp,par)(ímp,ímp).
Se Q_toro = 0, as 4 classes devem ser realizadas (idealmente ≈ 25% cada).
"""

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

from knight_tours_torus import knight_tours, verify_tour, _KNIGHT_MOVES


def _build_torus_decoder(n, m):
    """(diff_y, diff_x) ∈ (ℤ/n × ℤ/m) → (dy, dx) ∈ knight_moves.
    Para N,M ≥ 5 sem colisões."""
    table = {}
    collisions = set()
    for dy, dx in _KNIGHT_MOVES:
        key = (dy % n, dx % m)
        if key in table and table[key] != (dy, dx):
            collisions.add(key)
        else:
            table[key] = (dy, dx)
    for k in collisions:
        table[k] = None
    return table


def tour_windings(tour, n, m, decoder):
    """Retorna (w_y, w_x). Lança se houver passo ambíguo ou sanity check falhar."""
    V = n * m
    sum_dy = 0
    sum_dx = 0
    for i in range(V):
        u = int(tour[i])
        w = int(tour[(i + 1) % V])
        uy, ux = divmod(u, m)
        wy, wx = divmod(w, m)
        diff_y = (wy - uy) % n
        diff_x = (wx - ux) % m
        decoded = decoder.get((diff_y, diff_x))
        if decoded is None:
            raise ValueError(f"passo ambíguo (dy={diff_y}, dx={diff_x})")
        dy, dx = decoded
        sum_dy += dy
        sum_dx += dx
    if sum_dy % n != 0 or sum_dx % m != 0:
        raise RuntimeError(
            f"sanity falhou: Σdy={sum_dy} mod {n}; Σdx={sum_dx} mod {m}")
    return sum_dy // n, sum_dx // m


def dedup_tour_key(t, n, m):
    V = n * m
    ar = sorted(tuple(sorted([int(t[i]), int(t[(i + 1) % V])]))
                for i in range(V))
    return tuple(ar)


def medir_torus(n, m, K_por_seed, n_seeds, label=""):
    decoder = _build_torus_decoder(n, m)
    chaves = set()
    par_conj = Counter()
    w_dist = Counter()
    t0 = time.time()
    bruto = 0
    for s in range(n_seeds):
        tours = knight_tours(n, m, K_por_seed, seed=s)
        bruto += len(tours)
        for t in tours:
            assert verify_tour(t, n, m), "tour toroidal inválido"
            k = dedup_tour_key(t, n, m)
            if k in chaves:
                continue
            chaves.add(k)
            w_y, w_x = tour_windings(t, n, m, decoder)
            par_conj[(abs(w_y) % 2, abs(w_x) % 2)] += 1
            w_dist[(w_y, w_x)] += 1
    t1 = time.time()

    total = sum(par_conj.values())
    print(f"\n{'='*64}")
    print(f"  {label}  [toro {n}x{m}]   {n_seeds} seeds × K={K_por_seed}")
    print(f"{'='*64}")
    print(f"  tempo: {t1-t0:.2f}s  brutos={bruto}  únicos={total}")
    print(f"  distribuição conjunta (par_y, par_x):")
    classes = [(0, 0), (0, 1), (1, 0), (1, 1)]
    for c in classes:
        cnt = par_conj.get(c, 0)
        pct = 100 * cnt / total if total else 0
        print(f"    (y={c[0]}, x={c[1]}) : {cnt:6d}   {pct:5.2f}%")
    sigma = (total * 0.25 * 0.75) ** 0.5  # se H0: 25% uniforme
    print(f"  σ esperado por classe sob H0=25%: ±{sigma:.1f}")
    for c in classes:
        cnt = par_conj.get(c, 0)
        z = (cnt - total / 4) / sigma if sigma > 0 else 0
        print(f"    z(y={c[0]},x={c[1]}) = {z:+.2f}")

    # Marginal: w_y mod 2 e w_x mod 2 separados
    marg_y = Counter()
    marg_x = Counter()
    for (py, px), cnt in par_conj.items():
        marg_y[py] += cnt
        marg_x[px] += cnt
    print(f"  marginal w_y mod 2: par={marg_y[0]} ímpar={marg_y[1]}")
    print(f"  marginal w_x mod 2: par={marg_x[0]} ímpar={marg_x[1]}")

    # Teste de independência: par(y) ⊥ par(x) ?
    if total > 0:
        esperado = {(py, px): marg_y[py] * marg_x[px] / total
                    for py in (0, 1) for px in (0, 1)}
        chi2 = sum((par_conj.get(c, 0) - esperado[c]) ** 2 / esperado[c]
                   for c in classes if esperado[c] > 0)
        print(f"  χ² (indep. das paridades, df=1): {chi2:.3f}  "
              f"(crítico α=0.05: 3.84)")

    # Suporte completo de (w_y, w_x)? Quantas classes (w_y, w_x) distintas
    print(f"  classes (w_y, w_x) distintas observadas: {len(w_dist)}")
    return par_conj, w_dist


if __name__ == '__main__':
    medir_torus(6, 6, K_por_seed=1000, n_seeds=20, label="6×6 grande")
    medir_torus(6, 8, K_por_seed=500, n_seeds=10, label="6×8 sondagem")
    medir_torus(8, 8, K_por_seed=500, n_seeds=10, label="8×8 sondagem")
    medir_torus(8, 6, K_por_seed=500, n_seeds=10, label="8×6 sondagem")
    medir_torus(10, 10, K_por_seed=200, n_seeds=10, label="10×10 sondagem")
