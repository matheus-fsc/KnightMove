# Estimativas de contagem

Contar tours exatamente é viável até `n = 8`. Acima disso, estima-se.

## Estimador de Knuth

O estimador de Knuth para o tamanho de uma árvore de backtracking: percorre um
caminho aleatório da raiz até uma folha, multiplicando o número de opções em
cada nó, e usa a média dessas estimativas.

| `n` | estimativa | conferência |
|---|---|---|
| 6 | 9862 | exato, bate |
| 8 | bate a literatura | dentro do intervalo de confiança |
| 10 | **≈ 2,4 × 10^22** | `log10` em `[20,75 ; 22,83]` |
| 12 | **≈ 1,3 × 10^33** | — |

O coeficiente de variação é ~150 em `n = 10`, o que explica a largura do
intervalo: a árvore é muito desbalanceada.

## Contagem exata por matriz de transferência

DP de perfil quebrado: `n = 6` validado em 0,5s, com `lambda_1 = 70,48` e
`rho = 0,555`. O autovetor dominante diverge ~30% de `f∞`, efeito de borda —
`n = 6` é estreito demais para o regime assintótico.

`n = 8` explode: >187 milhões de estados no vértice 35, memória insuficiente em
Python puro. A cota de estados é ~1,63M mas a contagem real chega a ~2,5M; o
gargalo é velocidade de Python, não memória. Viável numa máquina só com
implementação em C ou numba.

## Contagem por divide-and-conquer

`dnc_count_exact.py` faz contagem exata de `N(12×12)` compondo blocos 6×6.
Ver [[Algoritmo]].

## Razão tours / 2-fatores

Relacionado, e importante para entender o limite do método:

| `n` | `r(n)` |
|---|---|
| 6 | 0,272 |
| 12 | 0,10 |

Ajuste `r ≈ 0.84 · e^(-0.17n)`; `log2(1/r)` cresce de 1,88 para 3,31. O
decaimento é exponencial, o que refuta a leitura de que "o deficit 3 explica a
conectividade" — 3 é constante, o gap não é.

**Código:** `experiments/08_heuristicas/tour_count_estimator.py`,
`solvers/transfer_matrix/`
