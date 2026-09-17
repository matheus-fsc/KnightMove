# Teoria: GF(2) e o espaço de ciclos

De onde sai o corpo `F_2` neste projeto. A história começa com **loops** —
caminhos que saem de uma casa e voltam a ela — e termina com tours sendo
vetores binários.

## 1. O que é um loop

No grafo do cavalo, um **loop** é um caminho fechado: sai de uma casa, segue
uma sequência de movimentos válidos, e retorna à casa inicial sem repetir
aresta. Visualmente, um círculo no grafo.

Os loops não foram procurados na teoria — foram notados desenhando a árvore de
recursão do backtracking à mão e percebendo que ramos diferentes chegam à mesma
configuração. Ver [[Historia]].

O modo "Loops" do
[visualizador](https://github.com/matheus-fsc/knight-tour-visualizer) anima
cada ciclo fundamental do 8×8 — são **105** no total.

## 2. Como o algoritmo acha os loops

O tabuleiro `n×n` tem `beta_1(n) = |E| - |V| + 1` loops independentes. Achá-los
é clássico: construir uma **árvore geradora por BFS** a partir de uma casa raiz.

> Cada aresta **fora** da árvore geradora fecha **exatamente um** ciclo
> fundamental — e esses ciclos são todas as `beta_1` dimensões do espaço.

```python
def compute_bfs(adj, start=(0, 0)):
    """BFS a partir de start; devolve a arvore geradora."""
    parent, level, order = {start: None}, {start: 0}, [start]
    queue, tree_edges = deque([start]), set()
    while queue:
        u = queue.popleft()
        for v in adj[u]:
            if v not in parent:              # aresta de arvore
                parent[v] = u
                level[v] = level[u] + 1
                order.append(v)
                queue.append(v)
                tree_edges.add(tuple(sorted([u, v])))
    return parent, level, order, tree_edges
```

A BFS visita os `V` vértices e seleciona `V-1` arestas. Sobram
`|E| - (V-1) = beta_1` arestas extras; cada uma, combinada com o caminho de
árvore entre seus extremos, fecha um ciclo. Esses `beta_1` ciclos formam uma
**base** do espaço de ciclos.

Cada loop aparece no visualizador em dois pedaços — um caminho laranja e um
azul — que se encontram na aresta amarela que fecha o ciclo. Essa forma em dois
caminhos é exatamente a estrutura "caminho de árvore + aresta extra".

## 3. Números concretos

| tabuleiro | `V` | `E` | `beta_1 = E - V + 1` |
|---|---|---|---|
| 6×6 | 36 | 80 | 45 |
| 8×8 | 64 | 168 | 105 |
| 10×10 | 100 | 288 | 189 |

## 4. Por que GF(2)

Represente um subconjunto de arestas como um vetor em `F_2^E`: coordenada 1 se
a aresta está presente, 0 se não. Então:

- a **soma** de dois subconjuntos é o XOR, que é a diferença simétrica;
- o conjunto dos subconjuntos em que **todo vértice tem grau par** é um
  subespaço vetorial — o **espaço de ciclos** `Z_1 = H_1(G; F_2)`;
- sua dimensão é exatamente `beta_1`.

Um tour fechado é um subconjunto de arestas em que todo vértice tem grau
exatamente 2. Grau 2 implica grau par, então **todo tour está no espaço de
ciclos**. Mais: a diferença simétrica de dois tours é sempre um elemento do
espaço de ciclos, e o conjunto dos tours vive num coset.

Isso reformula o problema: em vez de construir tours movimento a movimento,
combiná-los por XOR de ciclos.

## 5. Onde a reformulação quebra

Grau par é uma condição **linear** sobre `F_2^E`. Conectividade **não é**.

O espaço de ciclos contém todos os 2-fatores (uniões de ciclos disjuntos
cobrindo todos os vértices), não só os tours. E a proporção de 2-fatores que
são conexos cai exponencialmente:

| `n` | razão tours / 2-fatores |
|---|---|
| 6 | 0,272 |
| 10 | ~0,057 |
| 12 | ~0,10 (multi-seed) |

Ajuste: `r ≈ 0.84 · e^(-0.17n)`.

Essa é a limitação estrutural do método, e é a mesma que aparece em todos os
[[Resultados-Negativos]].

## 6. O que sobra: o deficit

Se todo tour está no espaço de ciclos, uma pergunta natural é: os tours
**geram** o espaço de ciclos inteiro?

Não. O rank do espaço gerado pelos tours fica sistematicamente **3 abaixo** de
`beta_1`. Esse é o assunto de [[Invariante-Q]].
