"""parallel_prefix.py — Paralelização por partição da árvore de busca.

Em vez de seeds aleatórias (que produzem subárvores sobrepostas e exigem
dedup), particiona o espaço por prefixos de atribuição. Cada worker
recebe uma subárvore disjunta → zero duplicatas, eficiência ≈100%.
"""

import os
import numpy as np
from multiprocessing import Pool

from .core import (
    State, OK, CONTRADICTION, SUBTOUR, COMPLETE,
    build_graph, fix_and_propagate, propagate_initial, choose_next_edge,
)


# ---------------------------------------------------------------------------
# T5 — Geração de prefixos
# ---------------------------------------------------------------------------

def generate_prefixes(ctx, depth=3, seed=42):
    """Gera prefixos de profundidade ≤ depth cobrindo TODO o espaço de busca.

    Cada prefixo é lista [(edge_idx, value), ...]. As subárvores definidas
    pelos prefixos são disjuntas e cobrem toda a árvore (depois da
    propagação inicial).

    Prefixos de comprimento < depth aparecem quando a propagação
    completa o tour (COMPLETE) ou descobre que o ramo é morto antes
    de atingir a profundidade — esses casos são tratados separadamente.
    """
    state = State(ctx['V'], ctx['E'])
    status = propagate_initial(state, ctx, mode='cycle')
    if status in (CONTRADICTION, SUBTOUR):
        return [], []  # nada a explorar
    if status == COMPLETE:
        return [[]], []   # tour completo na propagação inicial

    rng = np.random.default_rng(seed)
    prefixes = []
    complete_tours = []  # tours encontrados durante a expansão

    def expand(prefix, remaining_depth):
        if remaining_depth == 0:
            prefixes.append(prefix[:])
            return
        e, _ = choose_next_edge(state, ctx, rng)
        if e == -1:
            # Sem aresta para ramificar — guarda o prefixo (worker checa)
            prefixes.append(prefix[:])
            return
        for val in (0, 1):
            snap = state.snapshot()
            status = fix_and_propagate(state, ctx, e, val, mode='cycle')
            prefix.append((e, val))
            if status == OK:
                expand(prefix, remaining_depth - 1)
            elif status == COMPLETE:
                from .tours import _is_complete_tour, _extract_tour
                if _is_complete_tour(state, ctx):
                    complete_tours.append(_extract_tour(state, ctx))
            # CONTRADICTION/SUBTOUR: ramo morto, ignora
            prefix.pop()
            state.restore(snap)

    expand([], depth)
    return prefixes, complete_tours


# ---------------------------------------------------------------------------
# T6 — Worker e knight_tours_prefix_parallel
# ---------------------------------------------------------------------------

def _worker_prefix(args):
    n, prefix, K_target, seed, use_numba = args

    if use_numba:
        from .core_numba import (
            build_graph_numba, StateNumba,
            fix_and_propagate_nb, propagate_initial_nb,
            _is_complete_tour_nb, _extract_tour_nb, _backtrack_nb,
        )
        ctx = build_graph_numba(n)
        state = StateNumba(ctx['V'], ctx['E'])
        status = propagate_initial_nb(state, ctx, mode='cycle')
        if status in (CONTRADICTION, SUBTOUR):
            return []
        if status == COMPLETE:
            return [_extract_tour_nb(state, ctx)] if _is_complete_tour_nb(state, ctx) else []
        for (e, val) in prefix:
            status = fix_and_propagate_nb(state, ctx, e, val, mode='cycle')
            if status in (CONTRADICTION, SUBTOUR):
                return []
            if status == COMPLETE:
                if _is_complete_tour_nb(state, ctx):
                    return [_extract_tour_nb(state, ctx)]
                return []
        tours = []
        rng = np.random.default_rng(seed)
        _backtrack_nb(state, ctx, rng, tours, K_target)
        return tours

    # Fallback NumPy puro
    from .tours import _backtrack, _is_complete_tour, _extract_tour
    ctx = build_graph(n)
    state = State(ctx['V'], ctx['E'])
    status = propagate_initial(state, ctx, mode='cycle')
    if status in (CONTRADICTION, SUBTOUR):
        return []
    if status == COMPLETE:
        return [_extract_tour(state, ctx)] if _is_complete_tour(state, ctx) else []
    for (e, val) in prefix:
        status = fix_and_propagate(state, ctx, e, val, mode='cycle')
        if status in (CONTRADICTION, SUBTOUR):
            return []
        if status == COMPLETE:
            if _is_complete_tour(state, ctx):
                return [_extract_tour(state, ctx)]
            return []
    tours = []
    rng = np.random.default_rng(seed)
    _backtrack(state, ctx, rng, tours, K_target)
    return tours


def knight_tours_prefix_parallel(n, K, depth=3, n_workers=None, seed=None,
                                  use_numba=False, prefix_seed=42, verbose=False,
                                  exhaust=False):
    """Gera K tours via partição por prefixo.

    Args:
      depth: profundidade de partição (2^depth subárvores no máximo)
      n_workers: número de processos
      use_numba: usar core_numba.py nos workers
      prefix_seed: seed da geração de prefixos (determinístico)
      exhaust: se True, cada worker exaure sua subárvore (ignora K como limite
               por worker). Total pode ser muito > K; truncamos no fim.
               Use para enumeração exaustiva. K ainda atua como TETO global.
    """
    if K <= 0:
        return []
    if n_workers is None:
        n_workers = os.cpu_count()

    ctx = build_graph(n)
    prefixes, eager_tours = generate_prefixes(ctx, depth=depth, seed=prefix_seed)
    if verbose:
        print(f"  Prefixos: {len(prefixes)} subárvores "
              f"+ {len(eager_tours)} tours completos na geração")
    if not prefixes and not eager_tours:
        return []

    K_remaining = max(0, K - len(eager_tours))
    if exhaust:
        K_per = K  # teto suficientemente alto: cada worker exaure sua subárvore
    else:
        # margem 2× para compensar subárvores pequenas
        K_per = max(1, (2 * K_remaining + len(prefixes) - 1) // len(prefixes)) if prefixes else 0

    rng = np.random.default_rng(seed)
    seeds = rng.integers(0, 2**31 - 1, size=len(prefixes)).tolist() if prefixes else []

    out = list(eager_tours)

    if prefixes and K_remaining > 0:
        args = [(n, prefixes[i], K_per, int(seeds[i]), use_numba)
                for i in range(len(prefixes))]
        workers = min(n_workers, len(prefixes))
        with Pool(workers) as pool:
            results = pool.map(_worker_prefix, args)
        for batch in results:
            out.extend(batch)

    return out[:K]


# ---------------------------------------------------------------------------
# T7 — Benchmark final combinado
# ---------------------------------------------------------------------------

def full_benchmark_optimized(n=10, K=5000, n_workers=8):
    """Benchmark final cobrindo TODOS os modos: seq numpy/numba, paralelo seeds,
    paralelo prefixo (cap e exaustivo), com e sem Numba.

    O resultado reflete uma realidade importante:
      - par_seeds: rápido em amostragem, mas REDUNDANTE
        (workers exploram regiões sobrepostas)
      - par_prefix: cobertura única (subárvores disjuntas);
        beneficia exaustão / contagem
    """
    import time
    from .tours import knight_tours
    from .core_numba import knight_tours_numba, NUMBA_AVAILABLE
    from .parallel import knight_tours_parallel

    print(f"=== Benchmark completo n={n} K={K} workers={n_workers} ===\n")
    results = {}

    def report(label, t, n_tours, n_unique):
        cov = n_unique / max(n_tours, 1) * 100
        sp = results['numpy_seq'] / max(t, 1e-9)
        rate = n_tours / max(t, 1e-9)
        rate_u = n_unique / max(t, 1e-9)
        print(f"  {label:30s} t={t:6.3f}s  "
              f"tours={n_tours:>6} únicos={n_unique:>6} ({cov:5.1f}%)  "
              f"speedup={sp:.2f}×  ({rate_u:.0f} únicos/s)")
        return t

    def uniq(tours):
        return len(set(tuple(int(x) for x in t) for t in tours))

    # 1. NumPy seq
    t0 = time.perf_counter()
    tours = knight_tours(n, K, seed=0)
    results['numpy_seq'] = time.perf_counter() - t0
    print(f"  {'1. NumPy seq':30s} t={results['numpy_seq']:6.3f}s  "
          f"tours={len(tours):>6} únicos={uniq(tours):>6}  baseline 1.00×  "
          f"({len(tours)/results['numpy_seq']:.0f} tours/s)")

    # 2. Numba seq
    if NUMBA_AVAILABLE:
        knight_tours_numba(n, 5, seed=0)  # warmup
        t0 = time.perf_counter()
        tours = knight_tours_numba(n, K, seed=0)
        results['numba_seq'] = report('2. Numba seq', time.perf_counter() - t0,
                                       len(tours), uniq(tours))

    # 3. Par seeds (atual)
    t0 = time.perf_counter()
    tours = knight_tours_parallel(n, K, n_workers=n_workers, seed=0)
    results['par_seeds'] = report(f'3. Par seeds ({n_workers}w)',
                                   time.perf_counter() - t0,
                                   len(tours), uniq(tours))

    # 4. Par prefix (cap)
    t0 = time.perf_counter()
    tours = knight_tours_prefix_parallel(n, K, depth=4, n_workers=n_workers,
                                          seed=0, use_numba=False)
    results['par_prefix'] = report(f'4. Par prefix ({n_workers}w)',
                                    time.perf_counter() - t0,
                                    len(tours), uniq(tours))

    # 5. Par prefix + Numba (cap)
    if NUMBA_AVAILABLE:
        t0 = time.perf_counter()
        tours = knight_tours_prefix_parallel(n, K, depth=4, n_workers=n_workers,
                                              seed=0, use_numba=True)
        results['par_prefix_numba'] = report('5. Par prefix + Numba',
                                              time.perf_counter() - t0,
                                              len(tours), uniq(tours))

    # 6. Par prefix exaustivo + Numba — caso de uso de enumeração
    if NUMBA_AVAILABLE:
        t0 = time.perf_counter()
        tours = knight_tours_prefix_parallel(n, K, depth=4, n_workers=n_workers,
                                              seed=0, use_numba=True,
                                              exhaust=True)
        results['par_prefix_numba_exh'] = report('6. Par prefix+Numba exhaust',
                                                  time.perf_counter() - t0,
                                                  len(tours), uniq(tours))

    return results


if __name__ == '__main__':
    full_benchmark_optimized()
