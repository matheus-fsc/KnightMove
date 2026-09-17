"""Q1 — Diversidade dos tours D&C com matching aleatório.

Gera K tours com:
  - escolha aleatória de endpoints (s_BL, e_BL, s_BR, ...)
  - escolha aleatória de path interno por bloco (não sempre o primeiro)

Reporta:
  - Tours únicos por frozenset de arestas
  - Distribuição de multiplicidade
  - Órbitas D4 cobertas
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
from collections import Counter, defaultdict

import numpy as np

import knight_tours as kt
import knight_tours_dnc as dnc

N = dnc.N_GLOBAL


def edge_set_of_tour(tour):
    V = len(tour)
    return frozenset(
        frozenset((int(tour[i]), int(tour[(i + 1) % V])))
        for i in range(V)
    )


# ---------------------------------------------------------------------------
# D4 simetrias no tabuleiro 12×12
# ---------------------------------------------------------------------------

def _coord(v, n=N):
    return divmod(v, n)


def _to_idx(r, c, n=N):
    return r * n + c


def apply_d4(tour, op, n=N):
    """Aplica simetria do quadrado em todos os vértices.
    op ∈ {0..7}: 0=id, 1=rot90, 2=rot180, 3=rot270,
                  4=flip_h, 5=flip_v, 6=transp, 7=antitransp."""
    out = np.empty_like(tour)
    for i, v in enumerate(tour):
        r, c = _coord(int(v), n)
        if op == 0:
            rr, cc = r, c
        elif op == 1:
            rr, cc = c, n - 1 - r
        elif op == 2:
            rr, cc = n - 1 - r, n - 1 - c
        elif op == 3:
            rr, cc = n - 1 - c, r
        elif op == 4:
            rr, cc = r, n - 1 - c
        elif op == 5:
            rr, cc = n - 1 - r, c
        elif op == 6:
            rr, cc = c, r
        elif op == 7:
            rr, cc = n - 1 - c, n - 1 - r
        out[i] = _to_idx(rr, cc, n)
    return out


def orbit_canonical(edge_set, n=N):
    """Retorna o frozenset canônico (mínimo lexicográfico) da órbita D4."""
    # Edge set on vertices; aplicar todas 8 simetrias e pegar mínimo
    edges_arr = list(edge_set)
    # Para cada op, transforma cada aresta
    candidates = []
    for op in range(8):
        # Mapear cada vértice
        mapping = {}
        new_edges = set()
        for e in edges_arr:
            u, v = tuple(e)
            for x in (u, v):
                if x not in mapping:
                    r, c = _coord(x, n)
                    if op == 0:
                        rr, cc = r, c
                    elif op == 1:
                        rr, cc = c, n - 1 - r
                    elif op == 2:
                        rr, cc = n - 1 - r, n - 1 - c
                    elif op == 3:
                        rr, cc = n - 1 - c, r
                    elif op == 4:
                        rr, cc = r, n - 1 - c
                    elif op == 5:
                        rr, cc = n - 1 - r, c
                    elif op == 6:
                        rr, cc = c, r
                    elif op == 7:
                        rr, cc = n - 1 - c, n - 1 - r
                    mapping[x] = _to_idx(rr, cc, n)
            new_edges.add(frozenset((mapping[u], mapping[v])))
        candidates.append(frozenset(new_edges))
    # Pegar o canônico
    return min(candidates, key=lambda s: tuple(sorted((min(e), max(e)) for e in s)))


# ---------------------------------------------------------------------------
# Matching aleatório
# ---------------------------------------------------------------------------

def random_matching(catalog, cross, K, rng_seed=0, max_attempts=None):
    """Tenta K tentativas aleatórias. Retorna (tours, n_attempts, n_skel_found).
    Cada tentativa amostra (s_BL, e_BL) e tenta encontrar cross compatível
    em (BR, TR, TL). Se falha, conta como skel rejeitado.
    """
    rng = np.random.default_rng(rng_seed)
    fb = {b: list(dnc.get_frontier_vertices(b)) for b in dnc.BLOCKS}
    paths_by_se = catalog  # (s, e) -> [paths]

    bl_to_br = cross[('BL', 'BR')]
    br_to_tr = cross[('BR', 'TR')]
    tr_to_tl = cross[('TR', 'TL')]
    tl_to_bl = cross[('TL', 'BL')]

    idx_BL_BR = defaultdict(list)
    for (a, b) in bl_to_br:
        idx_BL_BR[a].append(b)
    idx_BR_TR = defaultdict(list)
    for (a, b) in br_to_tr:
        idx_BR_TR[a].append(b)
    idx_TR_TL = defaultdict(list)
    for (a, b) in tr_to_tl:
        idx_TR_TL[a].append(b)
    idx_TL_BL = defaultdict(list)
    for (a, b) in tl_to_bl:
        idx_TL_BL[a].append(b)

    fb_BL = fb['BL']; fb_BR = fb['BR']; fb_TR = fb['TR']; fb_TL = fb['TL']
    fb_BL_set = set(fb_BL); fb_BR_set = set(fb_BR)
    fb_TR_set = set(fb_TR); fb_TL_set = set(fb_TL)

    tours = []
    attempts = 0
    skels_failed = 0
    if max_attempts is None:
        max_attempts = K * 200

    # Pré-computar pares (s, e) por bloco (filtrando catálogo por frontier do bloco)
    def pairs_for_block(fb_set):
        return [(s, e) for (s, e) in paths_by_se.keys()
                if s in fb_set and e in fb_set and paths_by_se[(s, e)]]

    pairs_BL = pairs_for_block(fb_BL_set)
    pairs_BR = pairs_for_block(fb_BR_set)
    pairs_TR = pairs_for_block(fb_TR_set)
    pairs_TL = pairs_for_block(fb_TL_set)

    pairs_BL = [(s, e) for (s, e) in pairs_BL if e in idx_BL_BR]
    pairs_BR = [(s, e) for (s, e) in pairs_BR if e in idx_BR_TR]
    pairs_TR = [(s, e) for (s, e) in pairs_TR if e in idx_TR_TL]
    pairs_TL = [(s, e) for (s, e) in pairs_TL if e in idx_TL_BL]

    while len(tours) < K and attempts < max_attempts:
        attempts += 1
        # Amostrar (s_BL, e_BL)
        s_BL, e_BL = pairs_BL[rng.integers(len(pairs_BL))]
        candidates_BR_s = [v for v in idx_BL_BR[e_BL] if v in fb_BR_set]
        if not candidates_BR_s:
            skels_failed += 1
            continue
        s_BR = candidates_BR_s[rng.integers(len(candidates_BR_s))]
        # paths BR começando em s_BR
        cand_BR = [(s, e) for (s, e) in pairs_BR if s == s_BR]
        if not cand_BR:
            skels_failed += 1
            continue
        _, e_BR = cand_BR[rng.integers(len(cand_BR))]
        candidates_TR_s = [v for v in idx_BR_TR[e_BR] if v in fb_TR_set]
        if not candidates_TR_s:
            skels_failed += 1
            continue
        s_TR = candidates_TR_s[rng.integers(len(candidates_TR_s))]
        cand_TR = [(s, e) for (s, e) in pairs_TR if s == s_TR]
        if not cand_TR:
            skels_failed += 1
            continue
        _, e_TR = cand_TR[rng.integers(len(cand_TR))]
        candidates_TL_s = [v for v in idx_TR_TL[e_TR] if v in fb_TL_set]
        if not candidates_TL_s:
            skels_failed += 1
            continue
        s_TL = candidates_TL_s[rng.integers(len(candidates_TL_s))]
        # e_TL deve estar no catálogo com s_TL, e e_TL deve conectar a s_BL via TL→BL
        cand_TL = [(s, e) for (s, e) in pairs_TL
                   if s == s_TL and s_BL in idx_TL_BL.get(e, [])]
        if not cand_TL:
            skels_failed += 1
            continue
        _, e_TL = cand_TL[rng.integers(len(cand_TL))]

        # Sorteia paths internos
        p_BL = catalog[(s_BL, e_BL)][rng.integers(len(catalog[(s_BL, e_BL)]))]
        p_BR = catalog[(s_BR, e_BR)][rng.integers(len(catalog[(s_BR, e_BR)]))]
        p_TR = catalog[(s_TR, e_TR)][rng.integers(len(catalog[(s_TR, e_TR)]))]
        p_TL = catalog[(s_TL, e_TL)][rng.integers(len(catalog[(s_TL, e_TL)]))]

        tour = dnc._assemble_tour(p_BL, p_BR, p_TR, p_TL)
        if kt.verify_tour(tour, N):
            tours.append(tour)
        else:
            skels_failed += 1

    return tours, attempts, skels_failed


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    print("=" * 70)
    print(" Q1 — DIVERSIDADE DE 500 TOURS D&C (matching aleatório)")
    print("=" * 70)

    if not os.path.exists(dnc.CATALOG_FILE):
        raise SystemExit("Catálogo não existe. Rode knight_tours_dnc.py primeiro.")
    with open(dnc.CATALOG_FILE, 'rb') as f:
        obj = pickle.load(f)
    catalog = obj['catalog']
    cross = dnc.precompute_cross_edge_index()
    print(f"  Catálogo: {len(catalog)} pares, "
          f"{sum(len(v) for v in catalog.values())} paths")

    K = 500
    t0 = time.time()
    tours, attempts, fails = random_matching(catalog, cross, K, rng_seed=0)
    t1 = time.time()
    print(f"\n  Gerados: {len(tours)}/{K} tours em {t1 - t0:.2f}s")
    print(f"  Tentativas: {attempts}, falhas: {fails} (taxa rej={fails/attempts:.1%})")

    # Diversidade
    edge_hashes = [hash(edge_set_of_tour(t)) for t in tours]
    counts = Counter(edge_hashes)
    n_unique = len(counts)
    multiplicidade = Counter(counts.values())
    print(f"\n  Tours únicos: {n_unique}/{len(tours)} ({n_unique/len(tours):.1%})")
    print("  Distribuição (multiplicidade → #tours com essa mult):")
    for mult in sorted(multiplicidade.keys()):
        print(f"    aparece {mult}x: {multiplicidade[mult]} tours")

    # Órbitas D4 (sobre os primeiros 200 únicos para ser rápido)
    print("\n  Computando órbitas D4 (sobre primeiros 100 únicos)...")
    seen_canon = set()
    edge_sets_seen = {}
    t0 = time.time()
    distinct_count = 0
    for t in tours:
        h = hash(edge_set_of_tour(t))
        if h not in edge_sets_seen:
            edge_sets_seen[h] = edge_set_of_tour(t)
            if distinct_count < 100:
                canon = orbit_canonical(edge_sets_seen[h])
                seen_canon.add(canon)
            distinct_count += 1
    t1 = time.time()
    n_unique_for_orbit = min(100, n_unique)
    print(f"  Órbitas D4 distintas em {n_unique_for_orbit} amostras: "
          f"{len(seen_canon)} ({t1 - t0:.1f}s)")

    if n_unique / len(tours) > 0.9:
        veredito = "matching é EXPLORATÓRIO genuíno"
    elif n_unique / len(tours) < 0.5:
        veredito = "cap=50 é GARGALO — diversidade baixa"
    else:
        veredito = "diversidade intermediária"
    print(f"\n  Veredito: {veredito}")
