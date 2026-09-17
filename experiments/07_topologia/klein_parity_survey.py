"""klein_parity_survey.py — Varre as 4 paridades (n%2, m%2) em Klein
para validar a previsão de Q_klein:

  (n%2, m%2)   n+m   V    bipartido?   classe forçada    Q_klein
  ------------------------------------------------------------------
  (par, par)    0    par  intra-cor    (0,0) só          1
  (par, ímp)    1    par  preservado   todas 4           0
  (ímp, par)    1    par  preservado   todas 4           0
  (ímp, ímp)    0    ímp  intra-cor    (1,1) só          1
                                                          (W ≡ V mod 2)
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

from knight_tours_klein import knight_tours, verify_tour
from klein_winding import tour_winding, dedup_key


def varre(n, m, K_por_seed, n_seeds):
    chaves = set()
    par_conj = Counter()
    wraps_dist = Counter()
    wx_dist = Counter()
    bruto = 0
    sanity_ok = 0
    t0 = time.time()
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
                w_x, w_y_par, wraps_total, sum_dx, _ambig = \
                    tour_winding(t, n, m)
                par_conj[(w_y_par, abs(w_x) % 2)] += 1
                wraps_dist[wraps_total] += 1
                wx_dist[w_x] += 1
                if sum_dx % m == 0:
                    sanity_ok += 1
            except RuntimeError:
                pass
    t1 = time.time()
    return {
        'n': n, 'm': m,
        'tempo': t1 - t0,
        'bruto': bruto,
        'unicos': sum(par_conj.values()),
        'sanity_ok': sanity_ok,
        'par_conj': par_conj,
        'wraps_dist': wraps_dist,
        'wx_dist': wx_dist,
    }


def imprime(r):
    n, m = r['n'], r['m']
    total = r['unicos']
    print(f"\n{'=' * 72}")
    print(f"  [klein {n}×{m}]   n%2={n%2}  m%2={m%2}  (n+m)%2={(n+m)%2}  "
          f"V={n*m}  V%2={(n*m)%2}")
    print(f"{'=' * 72}")
    print(f"  tempo={r['tempo']:.2f}s  brutos={r['bruto']}  únicos={total}")
    print(f"  sanity Σdx≡0 mod m: {r['sanity_ok']}/{total}")
    if total == 0:
        print("  >>> NENHUM TOUR ENCONTRADO <<<")
        return
    classes = [(0, 0), (0, 1), (1, 0), (1, 1)]
    sigma = (total * 0.25 * 0.75) ** 0.5
    for c in classes:
        cnt = r['par_conj'].get(c, 0)
        pct = 100 * cnt / total
        z = (cnt - total / 4) / sigma if sigma else 0
        marker = "  <<<" if cnt == 0 else ""
        print(f"    (y={c[0]}, x={c[1]}) : {cnt:6d}   {pct:6.2f}%   "
              f"z={z:+7.2f}{marker}")
    classes_real = sum(1 for c in classes if r['par_conj'].get(c, 0) > 0)
    print(f"  classes realizadas: {classes_real}/4")
    wraps_set = sorted(r['wraps_dist'].keys())
    parity_wraps = {x % 2 for x in wraps_set}
    print(f"  #wraps observado : {wraps_set}  paridades={parity_wraps}")
    wx_set = sorted(r['wx_dist'].keys())
    parity_wx = {abs(x) % 2 for x in wx_set}
    print(f"  w_x observado    : {wx_set}  paridades={parity_wx}")


def predicao(n, m):
    """Retorna predição teórica (par_y, par_x) forçados ou None se livres."""
    nm_par = (n + m) % 2
    V = n * m
    if nm_par == 0:
        # bipartido quebrado pelo wrap intra-cor
        # W ≡ V mod 2
        forced = V % 2
        return (forced, forced)
    else:
        # bipartido preservado, sem restrição em W
        return None


if __name__ == '__main__':
    # Casos suficientemente grandes para Knight ser viável
    casos = [
        (6, 6),   # par, par   → Q=1, classe (0,0)
        (6, 8),   # par, par   → Q=1, classe (0,0)
        (5, 6),   # ímp, par   → Q=0, todas 4
        (6, 5),   # par, ímp   → Q=0, todas 4
        (5, 8),   # ímp, par   → Q=0
        (8, 5),   # par, ímp   → Q=0
        (7, 6),   # ímp, par   → Q=0
        (6, 7),   # par, ímp   → Q=0
        (5, 5),   # ímp, ímp   → V ímpar, Q=1 classe (1,1)
        (5, 7),   # ímp, ímp
        (7, 7),   # ímp, ímp
    ]
    resultados = []
    for n, m in casos:
        K_seed = 1000 if n * m <= 36 else 500
        n_seeds = 8
        try:
            r = varre(n, m, K_seed, n_seeds)
            imprime(r)
            pred = predicao(n, m)
            print(f"  PREDIÇÃO: classe forçada = {pred}  "
                  f"({'Q=0' if pred is None else 'Q=1'})")
            resultados.append((n, m, r, pred))
        except Exception as e:
            print(f"\n[{n}×{m}] ERRO: {e}")

    # Tabela resumo
    print()
    print("=" * 80)
    print(" TABELA RESUMO")
    print("=" * 80)
    print(f"  {'n':>2} {'m':>2} {'V%2':>3} {'(n+m)%2':>7}  "
          f"{'únicos':>7}  {'(0,0)':>6} {'(0,1)':>6} {'(1,0)':>6} {'(1,1)':>6}  "
          f"{'Q_emp':>5}  pred")
    for n, m, r, pred in resultados:
        total = r['unicos']
        cnts = [r['par_conj'].get(c, 0) for c in [(0, 0), (0, 1), (1, 0), (1, 1)]]
        classes_real = sum(1 for x in cnts if x > 0)
        q_emp = 4 - classes_real  # número de classes proibidas
        # Mostra % se houver tours
        if total > 0:
            cnt_s = [f"{100*c/total:5.1f}%" for c in cnts]
        else:
            cnt_s = ["—" for _ in cnts]
        pred_s = "Q=0" if pred is None else f"Q=1 cl{pred}"
        print(f"  {n:>2} {m:>2} {(n*m)%2:>3} {(n+m)%2:>7}  "
              f"{total:>7}  {cnt_s[0]:>6} {cnt_s[1]:>6} {cnt_s[2]:>6} {cnt_s[3]:>6}  "
              f"{q_emp:>5}  {pred_s}")
