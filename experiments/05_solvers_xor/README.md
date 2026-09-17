# 05 — Solvers SAT/SMT com cláusulas XOR

Com o espaço de ciclos definido, a tentativa natural: modelar no Z3, adicionar
as restrições como cláusulas XOR e deixar o solver enumerar. O ganho esperado
era eliminar os dead ends que dominam o custo do backtracking.

| diretório / arquivo | conteúdo |
|---|---|
| `knight_8x8_allsat_async.py` | AllSAT assíncrono no 8×8 |
| `dfs_knight_xorSearch.py` | busca DFS guiada por XOR |
| `cavalo_engine_z3` | engine Z3 inicial |
| `xor_clauses_benchmark/` | mede o speedup de adicionar cláusulas XOR ao Z3 |
| `benchmark_yen_vs_xor.py` + `benchmark_yen_vs_xor_results/` | algoritmo de Yen (k-shortest paths) contra a abordagem XOR |

**Veredito:** speedup ≈ 1,00x em 6×6 e em 10×10. As cláusulas XOR não ajudam o
Z3 — DPLL(XOR) já é padrão nos solvers modernos, então não havia ganho a extrair.
O benchmark teve, porém, um subproduto valioso: confirmou que
`deficit(10×10) = 3`, a mesma estrutura do 6×6.

O backtracking com union-find (`solvers/incremental_subtour/`) acabou sendo
~25x mais rápido que o Z3.
