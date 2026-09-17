"""cylinder_winding.py — Mede o winding number (e sua paridade) dos
tours fechados do cavalo no cilindro N×M.

Para cada tour fechado: w = Σ dx_signed / m  (inteiro, pois o tour fecha
em x). Reportamos a distribuição de w e de w mod 2.

Hipóteses-teste:
  - Q_cyl = 0  ⇒  ambas paridades realizadas (idealmente ≈ 50/50).
  - Q_cyl ≥ 1 ⇒  uma paridade é forçada (split degenera para 100/0
                  ou para padrão fixo).
"""

# --- bootstrap: torna core/ importavel a partir de qualquer subpasta ---
import sys as _sys, pathlib as _pathlib
for _anc in _pathlib.Path(__file__).resolve().parents:
    if (_anc / "core").is_dir() and (_anc / "pyproject.toml").exists():
        _sys.path.insert(0, str(_anc / "core"))
        break
# --- fim do bootstrap ---

import sys
import time
from collections import Counter

import numpy as np

from knight_tours_cylinder import knight_tours, verify_tour, _KNIGHT_MOVES


# Mapeamento (dy, dx_mod_m) -> dx_signed
def _build_dx_decoder(m):
    """Para cada (dy, diff = dx mod m), devolve o dx signed do movimento
    de cavalo que atinge aquele destino. Se houver colisão (M ≤ 4),
    devolve None nessa entrada (multigrafo)."""
    table = {}
    collisions = set()
    for dy, dx in _KNIGHT_MOVES:
        key = (dy, dx % m)
        if key in table and table[key] != dx:
            collisions.add(key)
        else:
            table[key] = dx
    for k in collisions:
        table[k] = None  # ambíguo
    return table


def tour_winding(tour, n, m, decoder):
    """Soma os dx signed ao longo do tour fechado. Retorna w = soma/m
    e a soma bruta para sanity check. Lança se houver passo ambíguo."""
    V = n * m
    soma = 0
    for i in range(V):
        u = int(tour[i])
        w_v = int(tour[(i + 1) % V])
        uy, ux = divmod(u, m)
        wy, wx = divmod(w_v, m)
        dy = wy - uy
        diff = (wx - ux) % m
        dx_signed = decoder.get((dy, diff))
        if dx_signed is None:
            raise ValueError(f"passo ambíguo em (dy={dy}, diff={diff}) — "
                             f"M={m} é pequeno demais para decodificar")
        soma += dx_signed
    if soma % m != 0:
        raise RuntimeError(
            f"sanity check falhou: Σdx={soma} não é múltiplo de m={m}")
    return soma // m, soma


def medir(n, m, K, seed=0, verbose=True):
    decoder = _build_dx_decoder(m)
    t0 = time.time()
    tours = knight_tours(n, m, K, seed=seed)
    t1 = time.time()
    if verbose:
        print(f"[{n}x{m}] {len(tours)} tours em {t1-t0:.2f}s")

    windings = []
    paridades = []
    for tour in tours:
        assert verify_tour(tour, n, m), "tour cilíndrico inválido"
        w, _ = tour_winding(tour, n, m, decoder)
        windings.append(w)
        paridades.append(abs(w) % 2)

    dist_w = Counter(windings)
    dist_par = Counter(paridades)

    if verbose:
        # Ordena windings para impressão
        for w in sorted(dist_w):
            print(f"  w = {w:+d}  : {dist_w[w]:5d}")
        print(f"  paridade 0 (par)   : {dist_par[0]:5d}")
        print(f"  paridade 1 (ímpar) : {dist_par[1]:5d}")
        total = dist_par[0] + dist_par[1]
        if total > 0:
            pct0 = 100 * dist_par[0] / total
            pct1 = 100 * dist_par[1] / total
            print(f"  % par={pct0:.1f}%  % ímpar={pct1:.1f}%")
    return dist_w, dist_par, tours


def teste_seeds_multiplas(n, m, K_por_seed, seeds, verbose=True):
    """Roda múltiplas seeds para reduzir viés do amostrador."""
    dist_w_total = Counter()
    dist_par_total = Counter()
    todos_tours = []
    for s in seeds:
        dw, dp, ts = medir(n, m, K_por_seed, seed=s, verbose=False)
        dist_w_total.update(dw)
        dist_par_total.update(dp)
        todos_tours.extend(ts)

    if verbose:
        print(f"\n[{n}x{m}] acumulado em {len(seeds)} seeds × K={K_por_seed} "
              f"= {sum(dist_par_total.values())} tours (com repetição entre seeds)")
        for w in sorted(dist_w_total):
            print(f"  w = {w:+d}  : {dist_w_total[w]:6d}")
        print(f"  par  : {dist_par_total[0]}")
        print(f"  ímpar: {dist_par_total[1]}")

    return dist_w_total, dist_par_total, todos_tours


def dedup_tours(tours, n, m):
    """Deduplica tours (mesmo conjunto de arestas)."""
    chaves = set()
    unicos = []
    for t in tours:
        V = n * m
        ar = []
        for i in range(V):
            u, w = int(t[i]), int(t[(i + 1) % V])
            ar.append(tuple(sorted([u, w])))
        ar.sort()
        k = tuple(ar)
        if k not in chaves:
            chaves.add(k)
            unicos.append(t)
    return unicos


if __name__ == '__main__':
    print("=" * 60)
    print(" Cilindro 6×6 — distribuição de winding")
    print("=" * 60)
    medir(6, 6, 200, seed=0)

    print()
    print("=" * 60)
    print(" Cilindro 6×6 — 8 seeds × 100 tours (depois dedupado)")
    print("=" * 60)
    _, _, todos = teste_seeds_multiplas(
        6, 6, 100, seeds=list(range(8)), verbose=True)
    unicos = dedup_tours(todos, 6, 6)
    print(f"\n  total bruto: {len(todos)}  únicos: {len(unicos)}")

    decoder = _build_dx_decoder(6)
    w_unicos = [abs(tour_winding(t, 6, 6, decoder)[0]) % 2 for t in unicos]
    dist = Counter(w_unicos)
    print(f"  paridade nos {len(unicos)} únicos: par={dist[0]} ímpar={dist[1]}")
    if dist[0] + dist[1] > 0:
        print(f"  % par={100*dist[0]/(dist[0]+dist[1]):.1f}%  "
              f"% ímpar={100*dist[1]/(dist[0]+dist[1]):.1f}%")

    print()
    print("=" * 60)
    print(" Cilindro 8×6 e 6×8 e 8×8 — sondas adicionais")
    print("=" * 60)
    for (n, m) in [(8, 6), (6, 8), (8, 8)]:
        try:
            medir(n, m, 200, seed=0)
            print()
        except Exception as e:
            print(f"[{n}x{m}] erro: {e}")
            print()
