"""Q2 — Por que 100% de sucesso?

Hipótese a verificar:
  "Todo par de perfis compatíveis na fronteira forma um tour válido."

Testes:
  1. cap=1 (apenas 1 path/par): taxa de sucesso mantém?
  2. Paths aleatórios do catálogo (não primeiro): mantém?
  3. Separar falhas: construção (sem compatível) vs verify_tour (tour inválido)
"""

# --- bootstrap: torna core/ importavel a partir de qualquer subpasta ---
import sys as _sys, pathlib as _pathlib
for _anc in _pathlib.Path(__file__).resolve().parents:
    if (_anc / "core").is_dir() and (_anc / "pyproject.toml").exists():
        _sys.path.insert(0, str(_anc / "core"))
        break
# --- fim do bootstrap ---

import os
import time
import pickle
from collections import defaultdict

import numpy as np

import knight_tours as kt
import knight_tours_dnc as dnc

N = dnc.N_GLOBAL


def random_matching_detailed(catalog, cross, K, rng_seed=0, max_attempts=None):
    """Como dnc_q1, mas separa falhas: construction vs verify."""
    rng = np.random.default_rng(rng_seed)
    fb = {b: list(dnc.get_frontier_vertices(b)) for b in dnc.BLOCKS}
    fb_BL_set = set(fb['BL']); fb_BR_set = set(fb['BR'])
    fb_TR_set = set(fb['TR']); fb_TL_set = set(fb['TL'])

    idx_BL_BR = defaultdict(list)
    for (a, b) in cross[('BL', 'BR')]:
        idx_BL_BR[a].append(b)
    idx_BR_TR = defaultdict(list)
    for (a, b) in cross[('BR', 'TR')]:
        idx_BR_TR[a].append(b)
    idx_TR_TL = defaultdict(list)
    for (a, b) in cross[('TR', 'TL')]:
        idx_TR_TL[a].append(b)
    idx_TL_BL = defaultdict(list)
    for (a, b) in cross[('TL', 'BL')]:
        idx_TL_BL[a].append(b)

    def pairs_for_block(fb_set):
        return [(s, e) for (s, e) in catalog.keys()
                if s in fb_set and e in fb_set and catalog[(s, e)]]

    pairs_BL = [(s, e) for (s, e) in pairs_for_block(fb_BL_set) if e in idx_BL_BR]
    pairs_BR = [(s, e) for (s, e) in pairs_for_block(fb_BR_set) if e in idx_BR_TR]
    pairs_TR = [(s, e) for (s, e) in pairs_for_block(fb_TR_set) if e in idx_TR_TL]
    pairs_TL = [(s, e) for (s, e) in pairs_for_block(fb_TL_set) if e in idx_TL_BL]

    by_start = defaultdict(list)
    for (s, e) in pairs_BR:
        by_start[('BR', s)].append(e)
    for (s, e) in pairs_TR:
        by_start[('TR', s)].append(e)
    for (s, e) in pairs_TL:
        by_start[('TL', s)].append(e)

    tours = []
    n_attempts = 0
    n_construct_fail = 0
    n_verify_fail = 0
    if max_attempts is None:
        max_attempts = K * 200

    while len(tours) < K and n_attempts < max_attempts:
        n_attempts += 1
        s_BL, e_BL = pairs_BL[rng.integers(len(pairs_BL))]
        BR_starts = [v for v in idx_BL_BR[e_BL] if v in fb_BR_set]
        if not BR_starts:
            n_construct_fail += 1
            continue
        s_BR = BR_starts[rng.integers(len(BR_starts))]
        BR_ends = by_start.get(('BR', s_BR))
        if not BR_ends:
            n_construct_fail += 1
            continue
        e_BR = BR_ends[rng.integers(len(BR_ends))]
        TR_starts = [v for v in idx_BR_TR[e_BR] if v in fb_TR_set]
        if not TR_starts:
            n_construct_fail += 1
            continue
        s_TR = TR_starts[rng.integers(len(TR_starts))]
        TR_ends = by_start.get(('TR', s_TR))
        if not TR_ends:
            n_construct_fail += 1
            continue
        e_TR = TR_ends[rng.integers(len(TR_ends))]
        TL_starts = [v for v in idx_TR_TL[e_TR] if v in fb_TL_set]
        if not TL_starts:
            n_construct_fail += 1
            continue
        s_TL = TL_starts[rng.integers(len(TL_starts))]
        TL_candidates = [e for e in by_start.get(('TL', s_TL), [])
                         if s_BL in idx_TL_BL.get(e, [])]
        if not TL_candidates:
            n_construct_fail += 1
            continue
        e_TL = TL_candidates[rng.integers(len(TL_candidates))]

        p_BL = catalog[(s_BL, e_BL)][rng.integers(len(catalog[(s_BL, e_BL)]))]
        p_BR = catalog[(s_BR, e_BR)][rng.integers(len(catalog[(s_BR, e_BR)]))]
        p_TR = catalog[(s_TR, e_TR)][rng.integers(len(catalog[(s_TR, e_TR)]))]
        p_TL = catalog[(s_TL, e_TL)][rng.integers(len(catalog[(s_TL, e_TL)]))]

        tour = dnc._assemble_tour(p_BL, p_BR, p_TR, p_TL)
        if kt.verify_tour(tour, N):
            tours.append(tour)
        else:
            n_verify_fail += 1

    return {
        'tours': tours,
        'attempts': n_attempts,
        'construct_fail': n_construct_fail,
        'verify_fail': n_verify_fail,
    }


if __name__ == '__main__':
    if not os.path.exists(dnc.CATALOG_FILE):
        raise SystemExit("Catálogo não existe. Rode knight_tours_dnc.py primeiro.")
    with open(dnc.CATALOG_FILE, 'rb') as f:
        obj = pickle.load(f)
    catalog_full = obj['catalog']
    cross = dnc.precompute_cross_edge_index()

    print("=" * 70)
    print(" Q2 — Teste estrutural: compatibilidade local ⇒ tour válido?")
    print("=" * 70)

    print("\n[Teste 1] Catálogo full (cap=50), paths aleatórios, K=500")
    t0 = time.time()
    res1 = random_matching_detailed(catalog_full, cross, K=500, rng_seed=42)
    t1 = time.time()
    print(f"  tours={len(res1['tours'])}, "
          f"tentativas={res1['attempts']}, "
          f"falhas_construção={res1['construct_fail']}, "
          f"falhas_verify={res1['verify_fail']} | {t1 - t0:.2f}s")

    # Teste 2: cap=1 (só 1 path/par)
    catalog_cap1 = {k: v[:1] for k, v in catalog_full.items()}
    print("\n[Teste 2] Catálogo reduzido cap=1, K=500")
    t0 = time.time()
    res2 = random_matching_detailed(catalog_cap1, cross, K=500, rng_seed=42)
    t1 = time.time()
    print(f"  tours={len(res2['tours'])}, "
          f"tentativas={res2['attempts']}, "
          f"falhas_construção={res2['construct_fail']}, "
          f"falhas_verify={res2['verify_fail']} | {t1 - t0:.2f}s")

    # Teste 3: forçar falha — mexer em paths para violar verify
    # Esse teste é só sanidade do detector: dropar 1 vértice de um path interno
    print("\n[Teste 3] Sanidade — path interno deliberadamente quebrado")
    bad_catalog = {k: list(v) for k, v in catalog_full.items()}
    # Corromper um path: trocar 2 vértices adjacentes
    k0 = next(iter(bad_catalog))
    p = list(bad_catalog[k0][0])
    p[10], p[20] = p[20], p[10]  # quebra knight move
    bad_catalog[k0] = [tuple(p)] + bad_catalog[k0][1:]
    # Esse teste é manual; só mostra que verify_tour pegaria a quebra
    rebuilt = dnc._assemble_tour(bad_catalog[k0][0], catalog_full[(7, 25)][0],
                                    catalog_full[(0, 6)][0], catalog_full[(4, 1)][0])
    print(f"  verify_tour(broken) = {kt.verify_tour(rebuilt, N)}")
    # Em condições normais o test deve dar False

    # Conclusão
    print("\n" + "=" * 70)
    print(" Conclusão estrutural")
    print("=" * 70)
    total_verify_fail = res1['verify_fail'] + res2['verify_fail']
    if total_verify_fail == 0:
        print("  Hipótese CONFIRMADA: 0 falhas em verify_tour em "
              f"{res1['attempts'] + res2['attempts']} tentativas.")
        print("  ⇒ Compatibilidade local (knight edges + endpoints válidos)")
        print("    IMPLICA tour hamiltoniano global no 12×12.")
        print("  As 'falhas' do random_matching são exclusivamente de")
        print("  CONSTRUÇÃO (par de perfis sem cross-edge), não de validade.")
    else:
        print(f"  Hipótese REFUTADA: {total_verify_fail} falhas em verify_tour.")
