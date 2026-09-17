# solvers — implementações de produção

| diretório | descrição |
|---|---|
| `incremental_subtour/` | **O resultado algorítmico principal.** Backtracking sobre o espaço de ciclos com detecção incremental de sub-ciclos por union-find. 25,2x sobre Z3 e 5,5x sobre a versão anterior no 10×10. A razão 2-fatores/tour cai de 18,1x para 1,00x, o que torna as podas R3/R6 redundantes. Custo por tour praticamente constante em `n ∈ {6,8,10,12,14}`. |
| `knight_tours_optimized/` | Empacotamento do solver com numba e paralelismo; é o que os benchmarks cruzados chamam de *flagship*. |
| `transfer_matrix/` | DP de perfil quebrado / matriz de transferência. `n=6` validado (9862 tours, 0.5s, `lambda_1 = 70.48`). `n=8` explode em número de estados (>187M) e estoura a memória — inviável em Python puro, ver `knight_transfer.cpp`. |

```bash
./venv/bin/python -m solvers.knight_tours_optimized.benchmark
```
