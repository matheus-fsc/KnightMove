"""Contagem exata N(12×12) via D&C.

FASE 1: count_catalog[(s,e)] = #HPs no 6×6 isolado, para todos pares
         (s,e) ∈ frontier × frontier com cor(s) ≠ cor(e).
         Usa simetria D4 do 6×6 para reduzir 512 → ~64 representantes.

FASE 2: enumera configurações de 4 cross-edges no 12×12 formando
         meta-Hamilton no 4-bloco (K_4): 3 topologias possíveis
         (1 perimetral + 2 diagonais).

FASE 3: N(12×12) = Σ_configs Π_blocos count[(s_i, e_i)].
"""

# --- bootstrap: torna core/ importavel a partir de qualquer subpasta ---
import sys as _sys, pathlib as _pathlib
for _anc in _pathlib.Path(__file__).resolve().parents:
    if (_anc / "core").is_dir() and (_anc / "pyproject.toml").exists():
        _sys.path.insert(0, str(_anc / "core"))
        break
# --- fim do bootstrap ---

import os
import sys
import time
import pickle
from collections import defaultdict, Counter

import numpy as np

import knight_tours as kt
import knight_tours_dnc as dnc


COUNT_CACHE_FILE = 'count_catalog_6x6_exact.pkl'
CONFIG_CACHE_FILE = 'count_n12x12_result.pkl'


# ===========================================================================
# Contador exato de HPs no 6×6 (zero storage de paths)
# ===========================================================================

def _count_backtrack(state, ctx, rng, counter):
    edge, first_val = kt._choose_next_edge(state, ctx, rng)
    if edge == -1:
        if kt._is_complete_tour(state, ctx):
            counter[0] += 1
        return
    for val in (first_val, 1 - first_val):
        snap = state.snapshot()
        status = kt.fix_and_propagate(state, ctx, edge, val)
        if status == kt.OK:
            _count_backtrack(state, ctx, rng, counter)
        elif status == kt.COMPLETE_TOUR:
            if kt._is_complete_tour(state, ctx):
                counter[0] += 1
        state.restore(snap)


def count_paths(s, e):
    """Conta HPs no 6×6 isolado de s a e via aresta virtual + R2 + UF."""
    if dnc._color(s) == dnc._color(e):
        return 0
    ctx, virtual = dnc._build_modified_ctx(s, e, n=dnc.N_BLOCK)
    V = ctx['V']
    state = kt.State(V, ctx['E'])
    queue = [(virtual, 1)]
    total_inc = ctx['total_incident']
    for w in range(V):
        if int(total_inc[w]) == 2:
            for ei in ctx['adj_edges'][w]:
                if int(ei) != virtual:
                    queue.append((int(ei), 1))
    status = kt._process_queue(state, ctx, queue)
    if status in (kt.CONTRADICTION, kt.SUBTOUR):
        return 0
    if status == kt.COMPLETE_TOUR:
        return 1 if kt._is_complete_tour(state, ctx) else 0
    counter = [0]
    rng = np.random.default_rng(0)
    _count_backtrack(state, ctx, rng, counter)
    return counter[0]


# ===========================================================================
# D4 sobre o 6×6
# ===========================================================================

def _d4_apply(v, op, n=dnc.N_BLOCK):
    r, c = divmod(v, n)
    if op == 0:
        rr, cc = r, c
    elif op == 1:  # rot 90 CW
        rr, cc = c, n - 1 - r
    elif op == 2:  # rot 180
        rr, cc = n - 1 - r, n - 1 - c
    elif op == 3:  # rot 270
        rr, cc = n - 1 - c, r
    elif op == 4:  # flip horiz (espelha colunas)
        rr, cc = r, n - 1 - c
    elif op == 5:  # flip vert
        rr, cc = n - 1 - r, c
    elif op == 6:  # transposta
        rr, cc = c, r
    elif op == 7:  # antitransposta
        rr, cc = n - 1 - c, n - 1 - r
    return rr * n + cc


def d4_pair_orbits(pairs):
    """Retorna lista de (representante, conjunto_orbit) cobrindo `pairs`."""
    pairs = set(pairs)
    visited = set()
    orbits = []
    for p in sorted(pairs):
        if p in visited:
            continue
        s, e = p
        orb = set()
        for op in range(8):
            sp = _d4_apply(s, op)
            ep = _d4_apply(e, op)
            orb.add((sp, ep))
        orb &= pairs
        # Representante = min
        rep = min(orb)
        if rep in visited:
            continue
        orbits.append((rep, frozenset(orb)))
        visited |= orb
    return orbits


# ===========================================================================
# FASE 1 — Catálogo exato
# ===========================================================================

def build_count_catalog(force=False, verbose=True):
    if not force and os.path.exists(COUNT_CACHE_FILE):
        with open(COUNT_CACHE_FILE, 'rb') as f:
            data = pickle.load(f)
        if verbose:
            print(f"  carregado de {COUNT_CACHE_FILE}: "
                  f"{len(data['catalog'])} pares")
        return data['catalog'], data['stats']

    # Pares (s,e) com s,e ∈ frontier_union e cores opostas
    frontier = sorted(set().union(
        *[dnc.get_frontier_vertices(b) for b in dnc.BLOCKS]))
    pairs = [(s, e) for s in frontier for e in frontier
             if s != e and dnc._color(s) != dnc._color(e)]
    if verbose:
        print(f"  {len(pairs)} pares candidatos (cor oposta, frontier∪)")

    orbits = d4_pair_orbits(pairs)
    if verbose:
        print(f"  D4 reduz para {len(orbits)} órbitas distintas")

    catalog = {}
    t_start = time.time()
    for i, (rep, orb) in enumerate(orbits):
        s, e = rep
        t0 = time.time()
        c = count_paths(s, e)
        t1 = time.time()
        # Propagar para toda órbita
        for p in orb:
            catalog[p] = c
        if verbose:
            elapsed = time.time() - t_start
            print(f"  [{i+1}/{len(orbits)}] orbit rep ({s},{e}): "
                  f"count={c}, orbit_size={len(orb)}, "
                  f"t={t1-t0:.2f}s, total={elapsed:.1f}s")
        # Checkpoint a cada 10 órbitas
        if (i + 1) % 10 == 0:
            with open(COUNT_CACHE_FILE, 'wb') as f:
                pickle.dump({'catalog': catalog,
                             'stats': {'partial': True, 'done': i + 1}}, f)

    counts = list(catalog.values())
    stats = {
        'n_pairs': len(catalog),
        'n_orbits': len(orbits),
        'total_paths_summed': sum(counts),
        'count_min': min(counts) if counts else 0,
        'count_max': max(counts) if counts else 0,
        'count_mean': sum(counts) / len(counts) if counts else 0,
        'count_median': sorted(counts)[len(counts) // 2] if counts else 0,
        'n_zero': sum(1 for c in counts if c == 0),
        'build_time_s': time.time() - t_start,
    }
    with open(COUNT_CACHE_FILE, 'wb') as f:
        pickle.dump({'catalog': catalog, 'stats': stats}, f)
    return catalog, stats


# ===========================================================================
# FASE 2 — Configurações cross-edge no 12×12
# ===========================================================================

# Topologias de meta-Hamilton em K_4 (BL, BR, TR, TL como nós):
# Existem 3 ciclos Hamilton em K_4 (não-orientados):
#   T1: BL-BR-TR-TL-BL (perimetral)
#   T2: BL-BR-TL-TR-BL (cruzado-a)
#   T3: BL-TR-BR-TL-BL (cruzado-b)
META_CYCLES = [
    ['BL', 'BR', 'TR', 'TL'],   # T1 perimetral
    ['BL', 'BR', 'TL', 'TR'],   # T2 cruzado-a
    ['BL', 'TR', 'BR', 'TL'],   # T3 cruzado-b
]


def enumerate_configs_for_cycle(catalog, cycle, verbose=False):
    """Enumera todas configs de cross-edges para a topologia `cycle`.

    Retorna (n_configs_validas, soma_de_contribuicoes).
    Cada cross-edge é (a→b) = (vertice em a, vertice em b) global.
    Para cada config:
      - Determina (s_i, e_i) local em cada bloco
      - Verifica s_i ≠ e_i e count(s_i, e_i) > 0
      - Contribuição = Π count(s_i, e_i)
    """
    n_blocks = len(cycle)
    # Cross-edges direcionais por par (a, b)
    cross = {}
    for i in range(n_blocks):
        a = cycle[i]
        b = cycle[(i + 1) % n_blocks]
        cross[(a, b)] = dnc.get_cross_edges(a, b)

    # Estratégia: iterar sobre as 4 cross-edges. Cada uma define
    # (vertice fora, vertice dentro) para 2 blocos.

    # Especificamente, para o ciclo a0 → a1 → ... → a_{n-1} → a0:
    # cross-edge i (entre a_i e a_{i+1}): (g_out_i, g_in_{i+1})
    # Para bloco a_i: e_i = g_out_i (local), s_i = g_in_i (local)

    contrib_total = 0
    n_valid = 0
    # Para n=4, iteramos exaustivamente
    cross_lists = [cross[(cycle[i], cycle[(i + 1) % n_blocks])]
                   for i in range(n_blocks)]
    if verbose:
        sizes = [len(cl) for cl in cross_lists]
        print(f"  cycle {cycle}: cross sizes {sizes} → "
              f"produto bruto {np.prod(sizes)}")

    # Pré-converter para local
    cross_local = []
    for i in range(n_blocks):
        a = cycle[i]
        b = cycle[(i + 1) % n_blocks]
        cl = [(dnc.global_to_local(ga, a), dnc.global_to_local(gb, b))
              for (ga, gb) in cross_lists[i]]
        cross_local.append(cl)

    # Iterar
    for ce0 in cross_local[0]:        # cycle[0] → cycle[1]
        e0_out, e1_in = ce0
        for ce1 in cross_local[1]:    # cycle[1] → cycle[2]
            e1_out, e2_in = ce1
            if e1_out == e1_in:
                continue
            for ce2 in cross_local[2]:  # cycle[2] → cycle[3]
                e2_out, e3_in = ce2
                if e2_out == e2_in:
                    continue
                for ce3 in cross_local[3]:  # cycle[3] → cycle[0]
                    e3_out, e0_in = ce3
                    if e3_out == e3_in or e0_out == e0_in:
                        continue

                    # (s_i, e_i) por bloco
                    se = [(e0_in, e0_out), (e1_in, e1_out),
                          (e2_in, e2_out), (e3_in, e3_out)]
                    # Lookup
                    c0 = catalog.get(se[0], 0)
                    if c0 == 0:
                        continue
                    c1 = catalog.get(se[1], 0)
                    if c1 == 0:
                        continue
                    c2 = catalog.get(se[2], 0)
                    if c2 == 0:
                        continue
                    c3 = catalog.get(se[3], 0)
                    if c3 == 0:
                        continue
                    n_valid += 1
                    contrib_total += c0 * c1 * c2 * c3
    return n_valid, contrib_total


def count_N_12x12(catalog, verbose=True):
    total = 0
    n_valid_total = 0
    per_cycle = []
    for ti, cycle in enumerate(META_CYCLES):
        t0 = time.time()
        n_valid, contrib = enumerate_configs_for_cycle(
            catalog, cycle, verbose=verbose)
        t1 = time.time()
        per_cycle.append((cycle, n_valid, contrib, t1 - t0))
        total += contrib
        n_valid_total += n_valid
        if verbose:
            print(f"  T{ti+1} {cycle}: configs={n_valid:,}, "
                  f"contrib={contrib:,}, t={t1-t0:.2f}s")
    return total, n_valid_total, per_cycle


# ===========================================================================
# Main
# ===========================================================================

if __name__ == '__main__':
    print("=" * 70)
    print(" Contagem EXATA N(12×12) via D&C")
    print("=" * 70)

    # Fase 1
    print("\n[FASE 1] count_catalog exato no 6×6 (com D4)")
    t0 = time.time()
    catalog, stats = build_count_catalog(force=False, verbose=True)
    t_phase1 = time.time() - t0
    print(f"\n  Estatísticas do catálogo:")
    for k, v in stats.items():
        if isinstance(v, float):
            print(f"    {k}: {v:.2f}")
        else:
            print(f"    {k}: {v:,}" if isinstance(v, int) else f"    {k}: {v}")
    print(f"  Tempo Fase 1: {t_phase1:.1f}s")

    # Fase 2 + 3
    print("\n[FASE 2+3] Enumeração de configs cross-edge + soma final")
    t0 = time.time()
    N, n_configs, per_cycle = count_N_12x12(catalog, verbose=True)
    t_phase2 = time.time() - t0

    print(f"\n  Configurações totais (3 topologias): {n_configs:,}")
    print(f"  Tempo Fase 2+3: {t_phase2:.2f}s")

    print("\n" + "=" * 70)
    print(" RESULTADO")
    print("=" * 70)
    print(f"  N_DnC(12×12) = {N:,}")
    print(f"  log10(N) = {np.log10(float(N)) if N > 0 else float('nan'):.2f}")
    print(f"  N / 8 (D4-canônicos) = {N // 8:,}")
    print(f"  N mod 8 = {N % 8}")

    # Top configs (por cycle)
    print("\n  Contribuições por topologia:")
    for (cycle, nv, ct, t) in per_cycle:
        frac = ct / N if N > 0 else 0
        print(f"    {cycle}: {nv:,} configs, "
              f"contrib={ct:,} ({frac:.4%}), {t:.2f}s")

    with open(CONFIG_CACHE_FILE, 'wb') as f:
        pickle.dump({
            'N': int(N),
            'n_configs': n_configs,
            'per_cycle': per_cycle,
            't_phase1': t_phase1,
            't_phase2': t_phase2,
        }, f)
    print(f"\n  Resultado salvo em {CONFIG_CACHE_FILE}")
