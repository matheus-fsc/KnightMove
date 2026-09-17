# TSP Cycle-Space vs 2-opt

Experimento de pesquisa que conecta a **estrutura do espaço de ciclos GF(2)**,
o mesmo arcabouço `H₁(G; F₂)` usado no projeto do passeio do cavalo, à
**busca local para o TSP Euclidiano**.

## Hipótese

A decomposição do espaço de ciclos `H₁(G; F₂)` pode guiar a busca local do TSP
de forma mais eficiente que o k-opt cego, identificando **quais ciclos
fundamentais têm o maior potencial de redução de custo**.

Cada aresta fora da árvore geradora mínima (MST) define um **ciclo
fundamental**. Um tour Hamiltoniano é um vetor em `F₂^|E|`. Somar
(XOR / diferença simétrica) o tour com um ciclo fundamental é um movimento de
busca local: arestas comuns são removidas, arestas exclusivas do ciclo são
adicionadas. Se o resultado continuar sendo um ciclo Hamiltoniano válido
(grau 2 em todo vértice + conexo) e tiver custo menor, aplicamos.

## Pergunta-chave

> Ordenar os ciclos fundamentais por **delta estimado** (estratégia *d*) torna
> a busca guiada por ciclos competitiva com o 2-opt?

O relatório é **honesto**: se o 2-opt vencer, dizemos por quanto.

## Pipeline (`tsp_cycle_experiment.py`)

1. **Instâncias**: pontos uniformes em `[0,100]²`, grafo completo com pesos
   Euclidianos. Seed mestre = 42.
2. **Tour inicial**: vizinho mais próximo (nearest neighbor) a partir do
   vértice 0.
3. **Base do espaço de ciclos**: MST (Kruskal/networkx); cada aresta não-árvore
   gera um ciclo fundamental `C_e = caminho(u→v na MST) + (u,v)`, representado
   como vetor binário em `F₂^|E|`. Total = `|E| - |V| + 1`.
4. **Melhoria guiada por ciclos**: para cada `C_e`, calcula
   `delta = custo(tour ⊕ C_e) - custo(tour)`; valida (grau 2 + conexo); aplica
   se `delta < 0`. Quatro ordenações testadas:
   - (a) aleatória
   - (b) por comprimento do ciclo (mais curto primeiro)
   - (c) por menor peso de aresta no ciclo (mais barato primeiro)
   - (d) por delta estimado (guloso)
   Repete até a convergência.
5. **Baseline 2-opt**: implementação padrão, mesma condição de parada.
6. **Benchmark**: `product([15,20,25], range(10))`, ambos os métodos a partir
   do **mesmo** tour NN. Registra custos, gap, iterações e tempos.

## Como rodar

```bash
./venv/bin/python tsp_cycle_space/tsp_cycle_experiment.py
```

Saída: tabela agregada no terminal + `results/tsp_cycle_experiment.json`
(resultados completos, salvos incrementalmente a cada 5 instâncias).

## Dependências

`numpy`, `networkx` (já presentes em `venv/`).
