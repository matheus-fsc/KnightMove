# 08 — Heurísticas, estimativas e divide-and-conquer

## Fase local f∞ — resultado negativo

A frequência de uso de cada aresta estabiliza por anel de distância à borda
quando `n` cresce (`L=1 ≈ 0.19`, `L=2 ≈ 0.20`, `L=3 ≈ 0.30`, `L=4 ≈ 0.25`).
A tabela teórica de 6 valores fixos substitui a amostragem e é 623x mais rápida
de montar.

**Mas:** o benchmark honesto (`benchmark_finf.py`) refuta a hipótese de que ela
acelera a busca. O speedup é constante ~1,0 em `n ∈ {6..16}`, com outliers
catastróficos para `n >= 12`. O ganho real vinha da pressão de vértice, não de f∞.

## Martingale

`H_MART` forte (custo i.i.d. constante por nó) **refutada**: spread de 28%,
chi-quadrado rejeita homogeneidade. A versão fraca se sustenta:
`mu_min >= 0.173` (Wilson 95%), logo `E[nós/tour] <= 5.79`. Variância linear.

## Razão tours / 2-fatores

`r(n)` decai exponencialmente: 0,272 em `n=6` para 0,10 em `n=12`, com ajuste
`r ≈ 0.84 · e^(-0.17n)`. Isso refuta a leitura de que "o deficit 3 explica a
conectividade" — `log2(1/r)` cresce de 1,88 para 3,31.

## Divide-and-conquer

Catálogo de 25,6k caminhos em blocos 6×6, universal (serve para 12×12 e 18×18).
Teorema de suficiência local (compatibilidade local implica hamiltonicidade
global) confirmado em 1639 ensaios. 18×18 em 0,01s contra timeout de 60s do
backtracking direto. Filtro direcional remove 42% das rejeições.

## Estimador de contagem

Estimador de Knuth (`tour_count_estimator.py`): `n=6` e `n=8` batem o ground
truth e a literatura dentro do intervalo de confiança. `N(10) ≈ 2.4e22`
(`log10` em `[20.75, 22.83]`), `N(12) ≈ 1.3e33`. Coeficiente de variação ~150 em `n=10`.

## Patch LUT

`investigate_patch.py`: a corretude do lema XOR pré-computado está OK
(9862 de 9862), mas a taxa de resgate de sub-tours é **0%** em 16k casos —
o union-find incremental funde componentes cedo demais para o patch disparar.
