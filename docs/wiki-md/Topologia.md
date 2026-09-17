# Topologia: como Q muda com a superfície

Se o deficit `Q = 3` vem dos cantos do tabuleiro (ver [[Invariante-Q]]), então
mudar a topologia deve mudar `Q`. Muda, e é isso que dá conteúdo geométrico ao
invariante.

## A hierarquia

| superfície | identificações | `Q` |
|---|---|---|
| **plano** `n×n` | nenhuma | **3** |
| **cilindro** | X modular, Y rígido | **0** |
| **toro** | ambas modulares | **0** |
| **Klein** | uma com inversão | **1** se `n+m` par, **0** se ímpar |

Interpretação: `Q` mede a obstrução introduzida pela **borda**. Numa
superfície sem borda, ela desaparece.

## Por que o toro tem deficit zero

A explicação não é a topologia: é a **simetria**.

O grafo do cavalo toroidal é um **grafo de Cayley** sobre `Z_n × Z_n`. Para
grafos de Cayley, deficit nulo é corolário de um teorema de **Alspach, Locke e
Witte (1990)**.

O que elimina o deficit é a **transitividade por vértices**: toda aresta é
equivalente a toda outra, então nenhuma pode ser forçada por um argumento
local. No plano, a fronteira quebra a transitividade e cria exatamente os
quatro vértices de grau 2 que geram `Q = 3`.

> Plano versus toro é, em essência, **com fronteira** versus **sem fronteira**.

Isso reposiciona a hierarquia abaixo: ela não é uma lista de medições soltas,
é uma consequência de quanto de transitividade cada superfície preserva.

## Cilindro

`beta_1 = 1`. Grafo 6×6: `V=36`, `E=108`, sem vértices de grau 2. A ancoragem
usa `C_m` em `x = 0`.

`Q_cilindro = 0`: a paridade do winding number se distribui 50/50 em 18k tours
únicos. Custo ~5 nós por tour.

## Toro

`beta_1 = 2`. Grau uniforme 8 em todo vértice: não há canto, não há borda.

`Q_toro = 0` confirmado: as 4 classes `(paridade_y, paridade_x)` aparecem com
~25% cada em 20k tours no 6×6. O rank é **cheio** (`rank = beta_1`) em `n ∈
{4, 6, 8, 10}`.

Isso refutou uma conjectura levantada durante o projeto, de que `rank = número
de órbitas`. Nota metodológica: uma medição anterior deu `rank ≈ 21` por usar
uma única semente; a medição correta exige muitas sementes.

## Garrafa de Klein

O resultado mais interessante da seção. Varrendo 11 combinações de paridade
(76k tours):

- a paridade off-diagonal `(0,1)` e `(1,0)` é **sempre proibida**;
- logo `H_1 = Z`, não `Z²`;
- `n+m` par implica `Q = 1`; `n+m` ímpar implica `Q = 0`.

## Toro cisalhado

Toro com parâmetro de cisalhamento `s`. Serve como família interpolante e para
testar se os resultados dependem da geometria específica ou só da topologia.

## Família híbrida toro→plano

Interpolando entre toro e plano em 10 pontos, `rank = beta_1 - Q` vale em
todos. Além disso, a geometria da identificação **não importa**, só o número
`k` de cantos com grau 2. Configurações geometricamente distintas com o mesmo
`k` dão o mesmo resultado.

**Código:** `experiments/07_topologia/`, `core/knight_tours_{torus,cylinder,klein,sheared}.py`
