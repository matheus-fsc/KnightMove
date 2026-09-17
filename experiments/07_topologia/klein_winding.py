"""klein_winding.py — Mede w_x (Z winding) e w_y (paridade de wraps com
flip) nos tours fechados na "garrafa de Klein" (Möbius com Y rígido).

Predição teórica (argumento bipartido):
  c(y,x) = (y+x) mod 2. Aresta wrap muda cor em (n+m) mod 2.
  Ciclo Ham fechado de comprimento V com W wraps:
    no_wraps flipam cor; wraps flipam cor sse n+m ímpar.
    Retornando à cor inicial:
       n+m ímpar : V par (sempre OK, sem restrição em W)
       n+m par   : V - W par  ⇒  W ≡ V (mod 2)

  No recobrimento (Möbius), g = (x+m, n-1-y), g² = (x+2m, y).
  Lift termina em g^k(0,0); Σ dx_signed = k·m;  k mod 2 = W mod 2.
  ⇒ par(w_x) = par(w_y).

  Para n,m ambos pares: V par ⇒ W par ⇒ w_x mod 2 = w_y mod 2 = 0.
  Apenas a classe (0,0) é realizada.  Q_klein = 1.
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

from knight_tours_klein import (
    knight_tours, verify_tour, build_graph, klein_step, _KNIGHT_MOVES,
)


def decode_step(u, v, n, m):
    """Determina (dy, dx, wrap) do movimento de cavalo que vai de u a v na Klein.
    Retorna lista de realizações (pode ser >1 em casos degenerados; usar
    o primeiro). Cada elemento: (dy, dx, wrap)."""
    y_u, x_u = divmod(u, m)
    y_v, x_v = divmod(v, m)
    realizacoes = []
    for dy, dx in _KNIGHT_MOVES:
        step = klein_step(y_u, x_u, dy, dx, n, m)
        if step is None:
            continue
        ny, nx, wrap = step
        if (ny, nx) == (y_v, x_v):
            realizacoes.append((dy, dx, wrap))
    return realizacoes


def tour_winding(tour, n, m):
    """Retorna (w_x, w_y, wraps_total, sum_dx, ambig_count).
       w_x = sum_dx // m   (deve ser inteiro)
       w_y = wraps_total mod 2"""
    V = n * m
    sum_dx = 0
    wraps = 0
    ambig = 0
    for i in range(V):
        u = int(tour[i])
        w = int(tour[(i + 1) % V])
        reals = decode_step(u, w, n, m)
        if not reals:
            raise RuntimeError(f"sem realização para passo {u}->{w}")
        if len(reals) > 1:
            ambig += 1
        dy, dx, wrap_flag = reals[0]
        sum_dx += dx
        wraps += wrap_flag
    if sum_dx % m != 0:
        raise RuntimeError(
            f"sanity falhou: Σdx={sum_dx} não múltiplo de m={m}")
    return sum_dx // m, wraps % 2, wraps, sum_dx, ambig


def dedup_key(t, n, m):
    V = n * m
    return tuple(sorted(tuple(sorted([int(t[i]), int(t[(i + 1) % V])]))
                        for i in range(V)))


def medir(n, m, K_por_seed, n_seeds, label=""):
    chaves = set()
    par_conj = Counter()
    w_dist = Counter()
    wraps_dist = Counter()
    sanidade_ok = 0
    sanidade_total = 0
    ambig_total = 0
    t0 = time.time()
    bruto = 0
    for s in range(n_seeds):
        tours = knight_tours(n, m, K_por_seed, seed=s)
        bruto += len(tours)
        for t in tours:
            assert verify_tour(t, n, m)
            k = dedup_key(t, n, m)
            if k in chaves:
                continue
            chaves.add(k)
            try:
                w_x, w_y_par, wraps_total, sum_dx, ambig = tour_winding(t, n, m)
                sanidade_total += 1
                if sum_dx % m == 0:
                    sanidade_ok += 1
                ambig_total += ambig
                par_conj[(w_y_par, abs(w_x) % 2)] += 1
                w_dist[w_x] += 1
                wraps_dist[wraps_total] += 1
            except RuntimeError as e:
                print(f"  WARN: {e}")
    t1 = time.time()

    total = sum(par_conj.values())
    print(f"\n{'=' * 70}")
    print(f"  {label}  [klein {n}x{m}]   {n_seeds} seeds × K={K_por_seed}")
    print(f"  paridades:  n%2={n%2}  m%2={m%2}  (n+m)%2={(n+m)%2}")
    print(f"{'=' * 70}")
    print(f"  tempo:                {t1-t0:.2f}s")
    print(f"  brutos:               {bruto}")
    print(f"  únicos:               {total}")
    print(f"  sanity Σdx≡0 mod m:   {sanidade_ok}/{sanidade_total}"
          f"  ({100*sanidade_ok/max(sanidade_total,1):.1f}%)")
    print(f"  passos ambíguos:      {ambig_total}")

    print(f"\n  distribuição conjunta (par_y, par_x):")
    classes = [(0, 0), (0, 1), (1, 0), (1, 1)]
    sigma = (total * 0.25 * 0.75) ** 0.5 if total else 0
    for c in classes:
        cnt = par_conj.get(c, 0)
        pct = 100 * cnt / total if total else 0
        z = (cnt - total / 4) / sigma if sigma else 0
        print(f"    (y={c[0]}, x={c[1]}) : {cnt:6d}   {pct:6.2f}%   z={z:+.2f}")

    # w_x ímpar
    impares = sum(cnt for c, cnt in par_conj.items() if c[1] == 1)
    print(f"\n  tours com w_x ÍMPAR:  {impares}/{total}  "
          f"({100*impares/max(total,1):.2f}%)")

    # marginais
    marg_y = Counter()
    marg_x = Counter()
    for (py, px), cnt in par_conj.items():
        marg_y[py] += cnt
        marg_x[px] += cnt
    print(f"  marginal par(w_y):    par={marg_y[0]}  ímpar={marg_y[1]}")
    print(f"  marginal par(w_x):    par={marg_x[0]}  ímpar={marg_x[1]}")

    # distribuição de w_x
    print(f"\n  distribuição w_x  (suporte: {len(w_dist)} valores)")
    for w in sorted(w_dist):
        print(f"    w_x = {w:+3d}  : {w_dist[w]:5d}")

    # distribuição de wraps
    print(f"\n  distribuição #wraps (paridade=w_y, suporte: {len(wraps_dist)})")
    for w in sorted(wraps_dist):
        print(f"    wraps = {w:3d}  : {wraps_dist[w]:5d}")

    return par_conj


if __name__ == '__main__':
    print("=" * 70)
    print(" HIPÓTESE: para n,m ambos pares, par(w_x) = par(w_y) = 0  (Q_klein=1)")
    print("           apenas a classe (0,0) deve ser realizada.")
    print("=" * 70)

    medir(6, 6, K_por_seed=2000, n_seeds=10, label="6×6 grande")
    medir(6, 8, K_por_seed=500, n_seeds=10, label="6×8 sondagem")
    medir(8, 6, K_por_seed=500, n_seeds=10, label="8×6 sondagem")
    medir(8, 8, K_por_seed=500, n_seeds=10, label="8×8 sondagem")
