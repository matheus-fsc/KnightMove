# O algoritmo

O solver mais rápido do projeto: **backtracking + espaço de ciclos + union-find
incremental**.

## O caminho até ele

1. **Backtracking puro** — dead ends dominam o custo. Inviável acima do 6×6.
2. **Z3 com cláusulas XOR** — elimina dead ends por construção, mas speedup
   ≈ 1,00x. DPLL(XOR) já é padrão nos solvers; não havia ganho a extrair.
3. **Backtracking sobre o espaço de ciclos com union-find** — o que funcionou.

## A ideia

Ao construir um tour incrementalmente, o erro caro não é escolher a aresta
errada — é descobrir tarde demais que as escolhas fecharam um **sub-ciclo**
antes de cobrir todos os vértices. Um 2-fator com dois ciclos disjuntos é um
beco sem saída que o backtracking ingênuo só detecta no final.

A solução é manter uma estrutura **union-find** sobre os segmentos já
construídos e testar, a cada aresta adicionada, se ela fecharia um ciclo
prematuro. O teste é O(α(n)) e poda o ramo imediatamente.

## Resultados

| métrica | valor |
|---|---|
| speedup sobre Z3 no 10×10 | **25,2x** |
| speedup sobre a versão anterior (v1) | 5,5x |
| razão 2-fatores explorados / tour encontrado | cai de 18,1x para **1,00x** |
| nós por tour, `n ∈ {6,8,10,12,14}` | 3,85 a 5,35 (praticamente constante) |
| tempo total na varredura de `n` | 0,76s a 3,34s |

A razão 1,00x é o ponto importante: o solver praticamente não gera mais
2-fatores desconexos. Isso torna as podas anteriores (R3, R6) redundantes.

## Podas que o compõem

| poda | efeito |
|---|---|
| grau de vértice | vértices de grau 2 forçam suas arestas (os 4 cantos, sempre) |
| pressão de vértice | prioriza vértices com menos opções restantes; é a fonte real do ganho de velocidade |
| union-find incremental | detecta sub-ciclos no momento em que se formariam |
| paridade de cor | descarta caminhos abertos entre casas de mesma cor |

Nota honesta: a heurística de fase local `f∞`, que parecia promissora, **não**
contribui — ver [[Resultados-Negativos]]. O ganho atribuído a ela vinha da
pressão de vértice.

## Divide-and-conquer

Linha paralela: catalogar todos os caminhos internos de blocos 6×6 (25,6k
entradas) e compor tabuleiros grandes por compatibilidade nas bordas.

- O catálogo é **universal**: serve para 12×12, 18×18, 24×24.
- Teorema de suficiência local (compatibilidade local implica hamiltonicidade
  global) confirmado em 1639 ensaios.
- 18×18 em **0,01s**, contra timeout de 60s do backtracking direto.
- Um filtro direcional nas bordas remove 42% das rejeições.

## Matriz de transferência

DP de perfil quebrado. `n=6` validado (9862 tours em 0,5s, `lambda_1 = 70,48`).
`n=8` explode em número de estados (>187M) e estoura a memória em Python puro.
Viável em C/numba, não implementado.

**Código:** `solvers/`, `experiments/08_heuristicas/dnc_*.py`
