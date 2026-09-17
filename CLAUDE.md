# KnightMove

Pesquisa sobre o passeio do cavalo: espaço de ciclos sobre GF(2) e o
invariante de deficit `Q(n) = 3`. Ver `README.md` para o mapa completo e
`docs/historia.md` para como o projeto chegou aqui.

## Idioma

Português (pt-BR), inclusive nomes de arquivo, comentários e relatórios.
É o padrão do repositório; código novo deve seguir.

## Ambiente

```bash
./venv/bin/python                      # tem numpy, networkx, numba, z3
./venv/bin/python benchmarks/benchmark_exhaustive_6x6.py   # regressao: 9862 tours
```

`9862` é o número canônico de tours fechados do 6×6. Qualquer mudança nos
solvers ou em `core/` deve manter esse valor.

## Convenção de imports

Scripts em `experiments/` e `benchmarks/` que usam as engines de `core/`
carregam um bootstrap de 5 linhas no topo, que sobe a árvore de diretórios
até achar a raiz e insere `core/` no `sys.path`. Ao criar um script novo que
importe do `core/`, copie esse bloco de um script vizinho.

## Organização

`experiments/` é numerado por **ordem cronológica das ideias**, não por
tema. Um experimento novo entra na pasta da etapa a que pertence, com
README próprio que termina em um **veredito**.

## Regras de honestidade

Este repositório vale pelo que documenta de negativo. Ao escrever qualquer
README, relatório ou commit:

- **Resultado negativo é resultado.** Documente com o mesmo cuidado dos
  positivos, e diga por que falhou. Ver `docs/wiki-md/Resultados-Negativos.md`.
- **Não reivindique novidade sem checar `docs/prior-art.md`.** O algoritmo,
  o XOR de ciclos e o enquadramento (Hamilton space) já são prior art.
- **Baseline honesto.** Speedups se medem contra `core/knight_tours.py`, não
  contra uma implementação ingênua.
- **Amostragem engana.** A verificação de tightness precisa de `K >= 250k`;
  rank no toro precisa de muitas sementes. Com pouco, aparecem falsos
  contraexemplos, e isso já aconteceu.

## Dados

`data/archives/` está em Git LFS. `data/raw/` não é versionado: contém
saídas regeneráveis, incluindo o catálogo de destruição de 352 MB.

## O que não está aqui

Os artigos em LaTeX ficam em `/home/math/Dev/knight_tour_papers/`, fora da
árvore de código, com seu próprio `CLAUDE.md` e agents de escrita. O estado
do paper está espelhado em `docs/paper-status.md`.
