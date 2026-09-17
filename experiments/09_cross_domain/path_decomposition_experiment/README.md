# Path-Finding via DFS Loop Decomposition (XOR)

Experimento de busca de **caminhos hamiltonianos `s → t`** no grafo do passeio
do cavalo, comparando uma abordagem de **decomposição em ciclos + XOR** contra
heurísticas clássicas (backtracking exaustivo e Warnsdorff repetido).

## Para que serve este diretório

Testar, de forma **honesta e mensurável**, a seguinte hipótese:

> Dado um grafo `G` e extremos `s, t`, o conjunto de todos os caminhos
> hamiltonianos `s → t` pode ser enumerado de forma mais eficiente:
> 1. achando **um** caminho-base `s → t`;
> 2. identificando todos os **ciclos fundamentais** de `G` (via DFS / base de ciclos);
> 3. gerando novos caminhos via **XOR** (diferença simétrica de conjuntos de arestas)
>    do caminho-base com combinações de ciclos;
> 4. usando backtracking apenas dentro de subconjuntos de ciclos relevantes.

O grafo do cavalo é ideal por ter estrutura **esparsa e fixa**, com *ground truth*
parcialmente conhecido (6×6 é enumerável exaustivamente).

## Conteúdo

| Arquivo | Descrição |
|---|---|
| `path_decomposition.py` | Script principal do experimento (Passos 1–6). |
| `RESULTS.md` | Análise dos resultados + respostas às 3 perguntas-chave. |
| `results/path_decomposition_experiment.json` | Resultados quantitativos (gerado pelo script). |

## TL;DR dos resultados

O XOR de ciclos fundamentais **falha como enumerador**: cobre **0,0086 %** dos
46 666 caminhos `(0,0)→(0,1)` do 6×6, satura em ~10 caminhos no 8×8, e os
encontra **menos diversos** que Warnsdorff e que a população real. 73–84 % das
combinações XOR sequer têm `V−1` arestas (ciclos da base não *alternam* com o
caminho-base). Detalhes em [`RESULTS.md`](RESULTS.md).

## Como rodar

```bash
../venv/bin/python path_decomposition.py
```

(Requer `networkx`; já instalado no `venv/` do projeto.)

## Nota de paridade (importante)

Num tabuleiro `n×n` **par**, `(0,0)` e `(n-1,n-1)` têm a **mesma cor**, mas um
caminho hamiltoniano tem `n²-1` arestas (número **ímpar**), exigindo extremos de
**cores opostas**. Logo `s=(0,0) → t=(n-1,n-1)` é **impossível** em 6×6 e 8×8.
O script detecta isso e usa o par viável `s=(0,0) → t=(0,1)` (open tour),
documentando a troca — conforme previsto no enunciado.

## Perguntas-chave (respondidas honestamente no output)

1. O XOR encontra caminhos válidos, ou a maioria das combinações produz
   subgrafos desconexos / com grau errado?
2. A diversidade dos caminhos achados via XOR é maior que a do Warnsdorff
   com partidas aleatórias repetidas?
3. Que fração dos ciclos é "compatível" com o caminho-base (gera caminho
   válido ao ser XOR-ado)?

As respostas são reportadas independentemente do resultado ser favorável
ou não à hipótese.
