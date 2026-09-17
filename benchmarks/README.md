# benchmarks

Comparações entre solvers, todas com a mesma semântica de contagem (tours
fechados não-direcionados) para que os números sejam comparáveis.

| arquivo | o que mede |
|---|---|
| `benchmark_exhaustive_6x6.py` | **teste de regressão canônico**: enumeração completa do 6×6 deve dar exatamente 9862 tours |
| `benchmark_8x8.py`, `benchmark_8x8_parallel.py` | comparação cruzada no 8×8: flagship, minimal_v2, naive, Warnsdorff, Z3 puro, Z3 com arestas obrigatórias |
| `benchmark_hybrid.py`, `benchmark_heavy.py` | cargas maiores |
| `benchmark_tabela9.py`, `benchmark_tabela9_completo.py` | reprodução da tabela comparativa do relatório |
| `cross_solver/` | baseline de backtracking, benchmark de Z3, divergência KL entre distribuições de tours, e os resultados em `cross_solver/results/` |

Aviso de interpretação: vários speedups reportados em iterações anteriores
deste projeto vinham de comparações com baselines fracos. Os números aqui usam
`knight_tours.py` como baseline.
