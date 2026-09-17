"""tests.py — Testes para o pacote knight_tours_optimized."""

import time
import numpy as np

from .tours import knight_tours, verify_tour
from .paths import knight_path, verify_path
from .symmetry import (
    d4_symmetries, canonical_form, is_canonical,
    expand_d4, verify_d4_decomposition,
)
from .parallel import knight_tours_parallel


def test_d4():
    print("=== Teste D₄ ===")
    # NOTA: o enunciado pedia 1232 canônicos, mas a contagem correta para
    # 9862 tours fechados no 6×6 modulo D₄ é 1245 (1232·8 = 9856 ≠ 9862;
    # 1245·8 - 98 = 9862 explicado por órbitas com simetria interna).
    canonical, unique = verify_d4_decomposition(
        6, expected_total=9862, expected_canonical=1245)
    syms = d4_symmetries(6)
    tours = knight_tours(6, 100, seed=0)
    for t in tours:
        cf = canonical_form(t, syms)
        assert is_canonical(cf, syms), \
            f"canonical_form retornou algo não-canônico"
    print("✓ D₄ correto")
    return canonical


def test_paths():
    print("\n=== Teste caminhos abertos ===")
    for start, end in [(0, 35), (0, 17), (1, 34)]:
        paths = knight_path(6, start, end, K=1, seed=0)
        if paths:
            assert verify_path(paths[0], 6, start, end), \
                f"verify_path falhou para {start}→{end}"
            print(f"  path {start}→{end}: ✓")
        else:
            print(f"  path {start}→{end}: nenhum encontrado")

    # n=10: par com cor compatível (V=100 par → cores opostas)
    # 0 cor 0; pedir end com cor 1: v=99 tem (9+9)%2=0 (mesma cor!),
    # então usamos (0, 98): 98 → (9,8) → (9+8)%2=1 (cor oposta).
    start, end = 0, 98
    t0 = time.perf_counter()
    paths = knight_path(10, start, end, K=1, seed=0)
    t1 = time.perf_counter()
    if paths:
        assert verify_path(paths[0], 10, start, end), \
            f"verify_path falhou para 10×10 {start}→{end}"
        print(f"  path 10×10 ({start}→{end}): ✓ em {t1 - t0:.3f}s")
    else:
        print(f"  path 10×10 ({start}→{end}): não encontrado em {t1 - t0:.3f}s")


def test_parallel():
    print("\n=== Teste paralelização ===")
    tours = knight_tours_parallel(6, 50, n_workers=2, seed=0)
    assert len(tours) == 50, f"esperava 50 tours, obteve {len(tours)}"
    assert all(verify_tour(t, 6) for t in tours)

    t0 = time.perf_counter()
    knight_tours(10, 200, seed=0)
    t_seq = time.perf_counter() - t0

    t0 = time.perf_counter()
    knight_tours_parallel(10, 200, n_workers=4, seed=0)
    t_par = time.perf_counter() - t0

    print(f"  Sequencial:       {t_seq:.3f}s")
    print(f"  Paralelo (4 wk):  {t_par:.3f}s")
    print(f"  Speedup:          {t_seq / max(t_par, 1e-9):.2f}×")
    print("✓ Paralelização correta")


if __name__ == '__main__':
    test_d4()
    test_paths()
    test_parallel()
