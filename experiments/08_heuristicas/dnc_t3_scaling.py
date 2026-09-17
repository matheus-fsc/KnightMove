"""T3 — Mapa de scaling N ∈ {12, 18, 24} para D&C vs BT.

Métricas:
  - tempo/tour (D&C amortizado, BT)
  - taxa de sucesso, fragilidade ao seed (BT)
  - razão entre N² (proporção esperada D&C: k² blocos, k²
    cross-edges, lookup O(1))
"""

# --- bootstrap: torna core/ importavel a partir de qualquer subpasta ---
import sys as _sys, pathlib as _pathlib
for _anc in _pathlib.Path(__file__).resolve().parents:
    if (_anc / "core").is_dir() and (_anc / "pyproject.toml").exists():
        _sys.path.insert(0, str(_anc / "core"))
        break
# --- fim do bootstrap ---

import time
import signal
import pickle
from collections import defaultdict

import numpy as np

import knight_tours as kt
import knight_tours_dnc as dnc
import dnc_q3_18x18 as q3
import dnc_t2_24x24 as t2


def bt_with_seeds(n, K, seeds, timeout_s):
    """Roda BT direto para vários seeds; retorna lista de (seed, tempo, sucesso)."""
    out = []

    class TO(Exception):
        pass

    def _h(s, f):
        raise TO()

    signal.signal(signal.SIGALRM, _h)

    for seed in seeds:
        signal.alarm(timeout_s)
        t0 = time.time()
        try:
            tours = kt.knight_tours(n, K, seed=seed)
            t1 = time.time()
            signal.alarm(0)
            out.append((seed, t1 - t0, len(tours)))
        except TO:
            signal.alarm(0)
            out.append((seed, float('inf'), 0))
    return out


def matcher_for_n(n):
    """Devolve função de matching D&C apropriada para n ∈ {12,18,24}."""
    if n == 12:
        import dnc_q1_diversity
        # Q1 usa o matching genérico de 4 blocos
        return lambda catalog, K, seed: dnc_q1_diversity.random_matching(
            catalog, dnc.precompute_cross_edge_index(), K, rng_seed=seed)[0]
    if n == 18:
        return lambda catalog, K, seed: q3.random_matching_18(
            catalog, K=K, rng_seed=seed)['tours']
    if n == 24:
        return lambda catalog, K, seed: t2.random_matching_24(
            catalog, K=K, rng_seed=seed)['tours']
    raise ValueError(n)


if __name__ == '__main__':
    print("=" * 70)
    print(" T3 — Mapa de scaling")
    print("=" * 70)

    with open(dnc.CATALOG_FILE, 'rb') as f:
        obj = pickle.load(f)
    catalog = obj['catalog']

    Ns = [12, 18, 24]
    K = 20
    seeds = [0, 1, 2, 7, 42]
    timeout_per_seed = 15  # s

    results = {}
    for n in Ns:
        print(f"\n[n={n}]")
        # D&C
        matcher = matcher_for_n(n)
        # Catalog é carregado em pickle (~0.01s — desprezar)
        dnc_times = []
        for seed in seeds:
            t0 = time.time()
            tours = matcher(catalog, K, seed)
            t1 = time.time()
            valid = sum(1 for t in tours if kt.verify_tour(t, n))
            dnc_times.append((seed, t1 - t0, len(tours), valid))
        med_dnc = np.median([t for (_, t, _, _) in dnc_times])
        print(f"  D&C: median {med_dnc:.4f}s, "
              f"per-tour {med_dnc/K*1000:.3f} ms")
        for (seed, t, ntours, valid) in dnc_times:
            print(f"    seed={seed}: {t:.4f}s, "
                  f"{ntours} tours, {valid} válidos")

        # BT direto
        print(f"  BT direto (timeout {timeout_per_seed}s/seed)...")
        bt = bt_with_seeds(n, K, seeds, timeout_per_seed)
        bt_success_times = [t for (_, t, ntours) in bt
                            if t != float('inf') and ntours == K]
        med_bt = (np.median(bt_success_times)
                  if bt_success_times else float('inf'))
        n_success = sum(1 for (_, t, n_t) in bt
                        if t != float('inf') and n_t == K)
        print(f"  BT: {n_success}/{len(seeds)} sucessos, "
              f"median (sucessos)={med_bt if med_bt != float('inf') else 'N/A'}")
        for (seed, t, ntours) in bt:
            ts = f"{t:.4f}s" if t != float('inf') else "TIMEOUT"
            print(f"    seed={seed}: {ts}, tours={ntours}")

        results[n] = {
            'dnc_med': med_dnc,
            'bt_med': med_bt if med_bt != float('inf') else None,
            'bt_success_frac': n_success / len(seeds),
        }

    print("\n" + "=" * 70)
    print(" RESUMO — Tabela de scaling")
    print("=" * 70)
    print(f"{'N':>4} {'blocos':>7} {'cross-edges':>12} {'D&C ms/tour':>13} "
          f"{'BT ms/tour':>13} {'BT robustez':>13} {'speedup D&C':>13}")
    print("-" * 95)
    for n in Ns:
        kb = (n // 6) ** 2
        r = results[n]
        dnc_ms = r['dnc_med'] / K * 1000
        if r['bt_med'] is not None:
            bt_ms = r['bt_med'] / K * 1000
            speedup = f"{bt_ms/dnc_ms:.1f}×"
        else:
            bt_ms = float('inf')
            speedup = "∞"
        bt_robust = f"{int(r['bt_success_frac']*100)}% seeds"
        bt_str = f"{bt_ms:.3f}" if bt_ms != float('inf') else "N/A"
        print(f"{n:>4} {kb:>7} {kb:>12} {dnc_ms:>13.4f} "
              f"{bt_str:>13} {bt_robust:>13} {speedup:>13}")

    print("\n  Scaling D&C esperado: tempo/tour ∝ k² (k=N/6)")
    print(f"    k=2: 1×;  k=3: 2.25×;  k=4: 4×")
    # Normalizar à n=12
    if results[12]['dnc_med'] > 0:
        base = results[12]['dnc_med']
        print(f"\n  Normalizado a n=12:")
        for n in Ns:
            k = n // 6
            obs = results[n]['dnc_med'] / base
            pred = (k / 2) ** 2
            print(f"    n={n}: observado {obs:.2f}×, "
                  f"previsto k²/4 = {pred:.2f}×")
