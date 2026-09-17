"""Fase 4 — Verificação parcial e análise do resultado N(12×12)."""

# --- bootstrap: torna core/ importavel a partir de qualquer subpasta ---
import sys as _sys, pathlib as _pathlib
for _anc in _pathlib.Path(__file__).resolve().parents:
    if (_anc / "core").is_dir() and (_anc / "pyproject.toml").exists():
        _sys.path.insert(0, str(_anc / "core"))
        break
# --- fim do bootstrap ---

import os
import pickle
from collections import defaultdict, Counter
import time

import numpy as np

import knight_tours as kt
import knight_tours_dnc as dnc
import dnc_count_exact as dce


def verify_random_tours(K=100, rng_seed=42):
    """Gera K tours D&C com matching aleatório e verifica que todos
    passam em verify_tour."""
    with open(dnc.CATALOG_FILE, 'rb') as f:
        obj = pickle.load(f)
    catalog_paths = obj['catalog']
    cross = dnc.precompute_cross_edge_index()
    import dnc_q1_diversity
    tours, atts, fails = dnc_q1_diversity.random_matching(
        catalog_paths, cross, K, rng_seed=rng_seed)
    valid = sum(1 for t in tours if kt.verify_tour(t, dnc.N_GLOBAL))
    return {
        'K': K,
        'gerados': len(tours),
        'validos': valid,
        'tentativas': atts,
        'falhas_construct': fails,
    }


def top_configs_by_contribution(count_catalog, topk=5):
    """Re-enumera configs por topologia e identifica os top-K e bottom-K."""
    results = []
    for cycle in dce.META_CYCLES:
        n_blocks = len(cycle)
        cross_lists = [dnc.get_cross_edges(cycle[i], cycle[(i+1) % n_blocks])
                       for i in range(n_blocks)]
        cross_local = []
        for i in range(n_blocks):
            a = cycle[i]
            b = cycle[(i+1) % n_blocks]
            cl = [(dnc.global_to_local(ga, a), dnc.global_to_local(gb, b))
                  for (ga, gb) in cross_lists[i]]
            cross_local.append(cl)
        configs_top = []
        for ce0 in cross_local[0]:
            e0_out, e1_in = ce0
            for ce1 in cross_local[1]:
                e1_out, e2_in = ce1
                if e1_out == e1_in: continue
                for ce2 in cross_local[2]:
                    e2_out, e3_in = ce2
                    if e2_out == e2_in: continue
                    for ce3 in cross_local[3]:
                        e3_out, e0_in = ce3
                        if e3_out == e3_in or e0_out == e0_in: continue
                        se = [(e0_in, e0_out), (e1_in, e1_out),
                              (e2_in, e2_out), (e3_in, e3_out)]
                        c0 = count_catalog.get(se[0], 0)
                        if c0 == 0: continue
                        c1 = count_catalog.get(se[1], 0)
                        if c1 == 0: continue
                        c2 = count_catalog.get(se[2], 0)
                        if c2 == 0: continue
                        c3 = count_catalog.get(se[3], 0)
                        if c3 == 0: continue
                        contrib = c0 * c1 * c2 * c3
                        configs_top.append(
                            (contrib, cycle, se, (c0, c1, c2, c3)))
        configs_top.sort(reverse=True)
        results.append({
            'cycle': cycle,
            'n_valid': len(configs_top),
            'top5': configs_top[:topk],
            'bottom5': configs_top[-topk:] if len(configs_top) >= topk else configs_top,
            'sum': sum(c[0] for c in configs_top),
        })
    return results


if __name__ == '__main__':
    print("=" * 70)
    print(" FASE 4 — Verificação e análise")
    print("=" * 70)

    # Verificação parcial
    print("\n[V1] K=100 tours D&C aleatórios — verificar todos válidos")
    res = verify_random_tours(K=100, rng_seed=42)
    print(f"  {res}")

    # Carregar count_catalog
    if not os.path.exists(dce.COUNT_CACHE_FILE):
        print("Catálogo exato não existe; rode dnc_count_exact.py primeiro.")
        raise SystemExit(1)
    with open(dce.COUNT_CACHE_FILE, 'rb') as f:
        data = pickle.load(f)
    catalog = data['catalog']
    print(f"\n  count_catalog: {len(catalog)} pares")

    # Top configs
    print("\n[V2] Top configurações por contribuição")
    t0 = time.time()
    per_cycle = top_configs_by_contribution(catalog, topk=5)
    print(f"  computado em {time.time()-t0:.2f}s")
    for r in per_cycle:
        print(f"\n  Topologia {r['cycle']}: {r['n_valid']:,} configs válidas, "
              f"soma={r['sum']:,}")
        print(f"    Top 5:")
        for (contrib, cyc, se, counts) in r['top5']:
            print(f"      contrib={contrib:.3e}  se={se}  counts={counts}")
        print(f"    Bottom 5:")
        for (contrib, cyc, se, counts) in r['bottom5']:
            print(f"      contrib={contrib:.3e}  se={se}  counts={counts}")

    # Carregar resultado final
    if os.path.exists(dce.CONFIG_CACHE_FILE):
        with open(dce.CONFIG_CACHE_FILE, 'rb') as f:
            result = pickle.load(f)
        N = result['N']
        print("\n" + "=" * 70)
        print(" Comparações")
        print("=" * 70)
        print(f"  N_DnC(12×12) = {N:,}")
        print(f"  log10(N) = {np.log10(float(N)):.2f}")
        print(f"  N / 8 (D4-canônicos) = {N // 8:,}")
        print(f"  N % 8 = {N % 8}")
        # Knuth estimate ~1.3e33 para n=12 (do contexto da memória)
        knuth_est = 1.3e33
        print(f"\n  Estimativa Knuth (memória): {knuth_est:.2e}")
        ratio = float(N) / knuth_est
        print(f"  Razão N_DnC / Knuth: {ratio:.2e}")
        if ratio < 1:
            print(f"  ⇒ N_DnC < N_total (esperado: D&C conta apenas tours")
            print(f"    com decomposição 1-path/bloco; tours multi-segmento")
            print(f"    são exclusivamente do BT direto).")
