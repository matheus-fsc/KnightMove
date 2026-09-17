"""Pipeline de execução do experimento de razão.

Faz checkpoint após cada n e usa unbuffered stdout (chamado com python -u).
"""
import json
import os
import sys
import time
import numpy as np

from ratio_analysis import sample_2factor_restart, sample_2factors_and_tours

CHECKPOINT = 'data/ratio_main_results.json'
os.makedirs('data', exist_ok=True)


def save(results):
    with open(CHECKPOINT, 'w') as f:
        json.dump(results, f, indent=2, default=str)


def load():
    if os.path.exists(CHECKPOINT):
        with open(CHECKPOINT) as f:
            return json.load(f)
    return {}


def run_exact(n, results):
    if str(n) in results and results[str(n)].get('mode') == 'exact':
        print(f"[skip] n={n} já executado (exact)")
        return
    print(f"\n=== n={n}: ENUMERAÇÃO EXATA ===", flush=True)
    t0 = time.time()
    r = sample_2factors_and_tours(n, K_target=10**9, seed=42,
                                    time_limit=120, use_heuristic=False,
                                    randomize_branch=False)
    dt = time.time() - t0
    print(f"  n_2fat={r['n_2fatores']}, tours={r['n_tours']}, "
          f"ratio={r['ratio']:.4f}, dt={dt:.1f}s", flush=True)
    results[str(n)] = {
        'mode': 'exact',
        'n_2fatores': int(r['n_2fatores']),
        'n_tours': int(r['n_tours']),
        'ratio': float(r['ratio']),
        'component_distribution': {str(k): int(v) for k, v in
                                    r['component_distribution'].items()},
        'elapsed_s': float(r['elapsed_s']),
    }
    save(results)


def run_restart(n, K, seeds, time_limit_per_seed, results):
    key = str(n)
    if key in results and results[key].get('mode') == 'restart' \
            and results[key].get('K_per_seed') == K \
            and results[key].get('seeds') == list(seeds):
        print(f"[skip] n={n} já executado (restart K={K} seeds={seeds})")
        return
    print(f"\n=== n={n}: RESTART K={K}, {len(seeds)} seeds ===", flush=True)
    seed_results = []
    for seed in seeds:
        t0 = time.time()
        r = sample_2factor_restart(n, K_target=K, seed=seed,
                                     time_limit=time_limit_per_seed,
                                     use_heuristic=False)
        dt = time.time() - t0
        seed_results.append(r)
        print(f"  seed={seed:3d}: K_obtido={r['n_2fatores']:5d}, "
              f"tours={r['n_tours']:5d}, ratio={r['ratio']:.4f}, "
              f"succ={r['success_rate']:.3f}, dt={dt:.1f}s, "
              f"stopped_by_time={r['stopped_by_time']}", flush=True)
    rs = [s['ratio'] for s in seed_results]
    print(f"  → ratio = {np.mean(rs):.4f} ± {np.std(rs):.4f} "
          f"(min={min(rs):.4f}, max={max(rs):.4f})", flush=True)
    # consolidar comp_dist da seed central
    comp_dist = seed_results[0]['component_distribution']
    comp_dist = {str(k): int(v) for k, v in comp_dist.items()}
    results[key] = {
        'mode': 'restart',
        'K_per_seed': K,
        'seeds': list(seeds),
        'ratios': [float(s['ratio']) for s in seed_results],
        'ratio_mean': float(np.mean(rs)),
        'ratio_std': float(np.std(rs)),
        'ratio_min': float(min(rs)),
        'ratio_max': float(max(rs)),
        'n_2fatores_per_seed': [int(s['n_2fatores']) for s in seed_results],
        'n_tours_per_seed':    [int(s['n_tours'])    for s in seed_results],
        'success_rates':       [float(s['success_rate']) for s in seed_results],
        'comp_dist': comp_dist,
    }
    save(results)


def main():
    results = load()
    print(f"Checkpoint inicial: {list(results.keys())}", flush=True)

    # n=6: enumeração exata
    run_exact(6, results)

    # n=6 RESTART (validação do estimador)
    if 'restart_n6_validation' not in results:
        print("\n=== n=6 VALIDAÇÃO RESTART (K=20000, 5 seeds) ===", flush=True)
        vals = []
        for seed in (1, 2, 3, 42, 100):
            r = sample_2factor_restart(6, K_target=20000, seed=seed,
                                         time_limit=120, use_heuristic=False)
            print(f"  seed={seed:3d}: ratio={r['ratio']:.4f}, "
                  f"succ={r['success_rate']:.3f}", flush=True)
            vals.append(r['ratio'])
        print(f"  → restart={np.mean(vals):.4f} ± {np.std(vals):.4f} "
              f"(true=0.2722, viés={(np.mean(vals)-0.2722)/0.2722*100:+.1f}%)",
              flush=True)
        results['restart_n6_validation'] = {
            'ratios': [float(v) for v in vals],
            'mean': float(np.mean(vals)),
            'std': float(np.std(vals)),
            'true_ratio': 0.2722,
            'bias_rel': float((np.mean(vals) - 0.2722) / 0.2722),
        }
        save(results)

    # n=8 — K=10000, 5 seeds, 180s per seed (largo)
    run_restart(8, K=10000, seeds=(1, 2, 3, 42, 100),
                 time_limit_per_seed=180, results=results)

    # n=10 — K=5000 (reduzido), 5 seeds, 180s
    run_restart(10, K=5000, seeds=(1, 2, 3, 42, 100),
                 time_limit_per_seed=180, results=results)

    # n=12 — K=2000 (reduzido), 5 seeds, 240s
    run_restart(12, K=2000, seeds=(1, 2, 3, 42, 100),
                 time_limit_per_seed=240, results=results)

    print("\n=== RESUMO FINAL ===", flush=True)
    for n in ('6', '8', '10', '12'):
        if n not in results:
            continue
        r = results[n]
        if r['mode'] == 'exact':
            print(f"  n={n}: ratio = {r['ratio']:.4f} (EXATO; "
                  f"N(2-fat)={r['n_2fatores']}, N(tours)={r['n_tours']})",
                  flush=True)
        else:
            print(f"  n={n}: ratio = {r['ratio_mean']:.4f} ± "
                  f"{r['ratio_std']:.4f} (RESTART K={r['K_per_seed']}×"
                  f"{len(r['seeds'])})", flush=True)
    print(f"\nSalvo: {CHECKPOINT}", flush=True)


if __name__ == '__main__':
    main()
