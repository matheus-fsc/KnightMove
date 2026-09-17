# Benchmarks

## Aviso de metodologia

Iterações anteriores deste projeto reportaram speedups que vinham de
comparações com baselines fracos. Os números abaixo usam
`core/knight_tours.py` como baseline, e o teste de regressão canônico (9862
tours no 6×6) roda antes de qualquer medição.

```bash
./venv/bin/python benchmarks/benchmark_exhaustive_6x6.py   # deve dar 9862
```

## Solvers comparados

| solver | descrição |
|---|---|
| `naive` | backtracking sem poda |
| `warnsdorff` | guloso por menor grau |
| `uniforme` | amostragem uniforme |
| `z3_puro` | modelo SMT direto |
| `z3_mandatory` | Z3 com as arestas obrigatórias dos cantos pré-fixadas |
| `minimal_v2` | backtracking com podas de grau |
| `flagship` | espaço de ciclos + union-find incremental |
| `theory` | flagship com a tabela `f∞` |

## Resultados

| comparação | resultado |
|---|---|
| flagship vs Z3, 10×10 | **25,2x** mais rápido |
| flagship vs versão anterior, 10×10 | 5,5x |
| cláusulas XOR no Z3 | speedup **1,00x** (sem efeito) |
| `f∞` (theory) vs flagship | ~1,0x, com outliers catastróficos em `n >= 12` |
| D&C vs backtracking, 18×18 | 0,01s vs timeout de 60s |

## Escalabilidade do flagship

| `n` | nós por tour | tempo total |
|---|---|---|
| 6 | 3,85 | 0,76s |
| 8 | ~4,3 | n/d |
| 10 | ~4,9 | n/d |
| 12 | ~5,1 | n/d |
| 14 | 5,35 | 3,34s |

Custo por tour praticamente constante: o crescimento vem do número de tours,
não da dificuldade de achar cada um.

## Limite teórico

Análise de martingale: `H_MART` forte (custo i.i.d. constante por nó) é
**refutada** (spread 28%, chi-quadrado rejeita homogeneidade). A versão fraca
se sustenta, com `mu_min >= 0.173` (Wilson 95%), logo `E[nós por tour] <=
5.79`, consistente com a tabela acima.

**Código:** `benchmarks/`, `experiments/08_heuristicas/martingale_analysis/`
