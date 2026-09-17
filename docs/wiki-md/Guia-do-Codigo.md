# Guia do código

Mapa do repositório [KnightMove](https://github.com/matheus-fsc/KnightMove).

## Estrutura

```
KnightMove/
├── core/                  engines de enumeracao (plano, toro, cilindro, Klein, cisalhado, D&C)
├── solvers/               union-find incremental, otimizado, matriz de transferencia
├── experiments/01..09/    um diretorio por etapa cronologica
├── benchmarks/            comparacoes entre solvers
├── formalization/lean/    Lean 4
├── viz/                   geracao de figuras
├── data/                  resultados + archives (Git LFS)
├── docs/                  historia, prior art, wiki, relatorios
└── tests/                 regressao
```

Os `experiments/` são numerados na ordem em que as ideias apareceram, não por
tema. A numeração é a mesma do fluxograma no README e de [[Historia]].

| # | diretório | etapa |
|---|---|---|
| 01 | `01_backtracking/` | backtracking e podas |
| 02 | `02_simetria_d4/` | simetria diedral |
| 03 | `03_amostragem/` | fixar início/fim, amostragem em massa |
| 04 | `04_espaco_ciclos/` | descoberta dos loops, GF(2) |
| 05 | `05_solvers_xor/` | Z3 e cláusulas XOR |
| 06 | `06_invariante_Q/` | o deficit `Q = 3` |
| 07 | `07_topologia/` | cilindro, toro, Klein |
| 08 | `08_heuristicas/` | `f∞`, martingale, razão, D&C, estimador |
| 09 | `09_cross_domain/` | TSP, pathfinding, proteínas |

## Rodando

```bash
python -m venv venv && ./venv/bin/pip install -e ".[fast,analysis]"

./venv/bin/python benchmarks/benchmark_exhaustive_6x6.py     # regressao: 9862
./venv/bin/python -m solvers.knight_tours_optimized.benchmark
for t in tests/*.py; do ./venv/bin/python "$t"; done
```

## Convenção de imports

Scripts em `experiments/` e `benchmarks/` que usam as engines de `core/`
carregam um bootstrap de 5 linhas no topo do arquivo:

```python
import sys as _sys, pathlib as _pathlib
for _anc in _pathlib.Path(__file__).resolve().parents:
    if (_anc / "core").is_dir() and (_anc / "pyproject.toml").exists():
        _sys.path.insert(0, str(_anc / "core"))
        break
```

Ele localiza a raiz do repositório subindo a árvore de diretórios, o que
permite rodar qualquer script de qualquer lugar sem `PYTHONPATH`.

## Dados

`data/archives/` guarda as amostras brutas das Fases A/B/C via **Git LFS**:
2,4 GB de JSON comprimidos para 66 MB.

`data/raw/` não é versionado. Inclui o catálogo de destruição de 352 MB,
regenerável por `experiments/01_backtracking/cavalo_loop_destruicao_6x6.py`.

## O que não está aqui

Os artigos em LaTeX e os PDFs ficam fora do repositório, por decisão
explícita: são rascunhos em revisão.

O frontend interativo está em
[knight-tour-visualizer](https://github.com/matheus-fsc/knight-tour-visualizer).
