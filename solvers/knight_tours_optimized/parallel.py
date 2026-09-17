"""parallel.py — Wrappers de multiprocessing para tours e caminhos.

Cada worker roda knight_tours/knight_path independentemente com seed própria.
Combina os resultados — não deduplica (tours possivelmente repetidos entre seeds).
"""

import numpy as np
from multiprocessing import Pool, cpu_count


def _worker_tours(args):
    n, K, seed = args
    from .tours import knight_tours
    return knight_tours(n, K, seed=seed)


def knight_tours_parallel(n, K, n_workers=None, seed=None):
    if K == 0:
        return []
    if n_workers is None:
        n_workers = cpu_count()
    n_workers = max(1, min(n_workers, K))
    rng = np.random.default_rng(seed)
    chunk = K // n_workers
    leftover = K - chunk * n_workers
    chunks = [chunk + (1 if i < leftover else 0) for i in range(n_workers)]
    seeds = rng.integers(0, 2**31 - 1, size=n_workers).tolist()
    args = [(n, chunks[i], int(seeds[i])) for i in range(n_workers)]
    with Pool(n_workers) as pool:
        results = pool.map(_worker_tours, args)
    out = []
    for r in results:
        out.extend(r)
    return out[:K]


def _worker_paths(args):
    n, start, end, K, seed = args
    from .paths import knight_path
    return knight_path(n, start, end, K=K, seed=seed)


def knight_path_parallel(n, start, end, K, n_workers=None, seed=None):
    if K == 0:
        return []
    if n_workers is None:
        n_workers = cpu_count()
    n_workers = max(1, min(n_workers, K))
    rng = np.random.default_rng(seed)
    chunk = K // n_workers
    leftover = K - chunk * n_workers
    chunks = [chunk + (1 if i < leftover else 0) for i in range(n_workers)]
    seeds = rng.integers(0, 2**31 - 1, size=n_workers).tolist()
    args = [(n, start, end, chunks[i], int(seeds[i])) for i in range(n_workers)]
    with Pool(n_workers) as pool:
        results = pool.map(_worker_paths, args)
    out = []
    for r in results:
        out.extend(r)
    return out[:K]
