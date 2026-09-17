# Pathfinding: GF(2) cycle-space XOR vs A* repetido

Testa se a decomposição em **espaço de ciclos GF(2)** enumera caminhos
alternativos `s→t` em mapas-grade de forma eficiente, comparada a rodar A*
várias vezes.

## Hipótese

O experimento anterior no TSP (grafo completo **denso**) mostrou que o XOR de
ciclos fundamentais falha: só ~3% dos ciclos preservavam a hamiltonicidade.
Aqui a aposta é o **oposto**: em grafos **planares esparsos** (grades com
obstáculos) os ciclos fundamentais seriam "corredores" geometricamente locais,
e o XOR do caminho-base com um ciclo quase sempre daria um caminho válido.

## Como rodar

```bash
../venv/bin/python pathfinding_xor.py
```

(Requer `networkx`; já instalado no `venv/` do projeto. Python 3.14.)

## Setup

- Tamanhos: 20×20, 30×30, 50×50; 20 instâncias cada (seeds 0..19).
- 30% de células bloqueadas; `S=(0,0)`, `T=(n-1,n-1)` sempre livres.
- Garante caminho `S→T`; se não houver, tenta `seed+1` (registra `used_seed`).
- Grafo: grade 4-conectada, peso 1; restrito à componente conexa de `S`.
- Caminho-base: A* com heurística Manhattan.
- Ciclos fundamentais: árvore DFS a partir de `S` + back edges.
- Candidato = `base XOR ciclo` (diferença simétrica de arestas).
- Validade: `S,T` grau 1; demais grau 2; conexo; sem célula bloqueada.
- Baseline: A* repetido com remoção da aresta-do-meio (k-shortest aprox.).

## Saídas

| Arquivo | Descrição |
|---|---|
| `pathfinding_xor.py` | Script principal (Passos 1-9). |
| `results/pathfinding_xor_experiment.json` | Métricas completas por instância. |
| `results/pathfinding_xor_per_instance.csv` | Mesmas métricas em CSV. |
| `RESULTADOS.md` | Relatório e respostas às perguntas-chave. |

## Nota sobre o tempo (honestidade)

`time_xor_ms` cobre **DFS + todas as verificações XOR** (single + 500 pares),
enumerando dezenas a >100 caminhos, não apenas 10. `time_astar_repeated_ms`
cobre achar **até 10** caminhos. Logo o "speedup" reportado **não** é
caminho-a-caminho; veja `RESULTADOS.md` para a leitura correta.
