"""benchmark.py — Benchmark completo."""

import time

from .tours import knight_tours
from .paths import knight_path
from .symmetry import knight_tours_canonical, expand_d4
from .parallel import knight_tours_parallel


def benchmark_parallel(n, K):
    t0 = time.perf_counter()
    knight_tours(n, K, seed=0)
    t_seq = time.perf_counter() - t0
    print(f"  Sequencial:        {t_seq:.3f}s")
    for nw in (2, 4, 8):
        t0 = time.perf_counter()
        knight_tours_parallel(n, K, n_workers=nw, seed=0)
        t_par = time.perf_counter() - t0
        sp = t_seq / max(t_par, 1e-9)
        print(f"  Paralelo ({nw} wk): {t_par:.3f}s  speedup={sp:.2f}×")


def full_benchmark():
    print("=== knight_tours_optimized — benchmark completo ===\n")

    print("1. Paralelização (n=10, K=500):")
    benchmark_parallel(n=10, K=500)

    print("\n2. Simetria D₄ (n=6):")
    t0 = time.perf_counter()
    # 1245 = contagem correta de órbitas D₄ (1232 do enunciado é inconsistente com 9862)
    canon = knight_tours_canonical(6, 1245, seed=0)
    t_canon = time.perf_counter() - t0

    t0 = time.perf_counter()
    expanded = expand_d4(canon, 6)
    t_expand = time.perf_counter() - t0

    total = t_canon + t_expand
    print(f"  Gerar {len(canon)} canônicos: {t_canon:.3f}s")
    print(f"  Expandir → {len(expanded)} tours: {t_expand:.3f}s")
    print(f"  Total:                       {total:.3f}s")
    print(f"  vs original (5.8s):           {5.8 / max(total, 1e-9):.2f}×")

    print("\n3. Caminhos abertos:")
    # Pares originais do enunciado (canto→canto): todos com cores iguais → aviso/✗
    # Pares válidos (cor oposta para V par): start=0 → end de cor 1
    cases = [
        ("enunciado canto-canto (mesma cor)", [(6, 0, 35), (8, 0, 63), (10, 0, 99)]),
        ("válidos (cor oposta)",              [(6, 0, 17), (8, 0, 33), (10, 0, 98)]),
    ]
    for label, pairs in cases:
        print(f"  -- {label}:")
        for n, s, e in pairs:
            t0 = time.perf_counter()
            paths = knight_path(n, s, e, K=1, seed=0)
            t1 = time.perf_counter()
            status = "✓" if paths else "✗"
            print(f"     n={n} ({s}→{e}): {status} em {t1 - t0:.3f}s")


if __name__ == '__main__':
    full_benchmark()
