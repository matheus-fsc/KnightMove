"""knight_tours_optimized — pacote de geração otimizada de tours do cavalo.

Adições sobre knight_tours.py:
  - Decomposição por simetria D₄ (canonical_form, expand_d4)
  - Paralelização via multiprocessing
  - Caminhos Hamiltonianos abertos (knight_path)
"""

from .tours import knight_tours, verify_tour
from .paths import knight_path, verify_path
from .symmetry import (
    d4_symmetries, canonical_form, is_canonical,
    expand_d4, verify_d4_decomposition, knight_tours_canonical,
)
from .parallel import knight_tours_parallel, knight_path_parallel

try:
    from .core_numba import knight_tours_numba, NUMBA_AVAILABLE
    from .parallel_prefix import (
        generate_prefixes, knight_tours_prefix_parallel,
        full_benchmark_optimized,
    )
    _OPTIMIZED = True
except ImportError:
    _OPTIMIZED = False
    NUMBA_AVAILABLE = False

__all__ = [
    'knight_tours', 'verify_tour',
    'knight_path', 'verify_path',
    'd4_symmetries', 'canonical_form', 'is_canonical',
    'expand_d4', 'verify_d4_decomposition', 'knight_tours_canonical',
    'knight_tours_parallel', 'knight_path_parallel',
]
if _OPTIMIZED:
    __all__ += [
        'knight_tours_numba', 'NUMBA_AVAILABLE',
        'generate_prefixes', 'knight_tours_prefix_parallel',
        'full_benchmark_optimized',
    ]
