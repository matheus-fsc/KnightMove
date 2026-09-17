# O invariante de deficit Q

O resultado central do projeto.

## A anomalia

Todo tour fechado está no espaço de ciclos `Z_1 = H_1(G; F_2)`, de dimensão
`beta_1 = |E| - |V| + 1` (ver [[Teoria-GF2]]). A pergunta: os tours **geram**
esse espaço?

Não. Medindo o rank do espaço gerado por todos os tours de um tabuleiro:

| `n` | `beta_1` | rank(Ham) | deficit |
|---|---|---|---|
| 6 | 45 | 42 | **3** |
| 8 | 105 | 102 | **3** |
| 10 | 189 | 186 | **3** |
| 12 | 297 | 294 | **3** |

Sempre 3. Esse deficit é o invariante `Q`.

## Ressalva de prior art (leia antes)

Que o grafo do cavalo **não** seja *Hamilton-generated*, isto é, que `deficit
> 0`, **é conhecido**. Heinig (2013) observou que `delta(G) >= 3` é necessário
para os ciclos hamiltonianos gerarem o espaço de ciclos, e o cavalo tem quatro
cantos de grau 2. O argumento é elementar: se `deg(v) = 2` com arestas `e_1,
e_2`, todo tour contém ambas, logo todo tour satisfaz `x_{e1} = x_{e2}`, um
hiperplano próprio.

O conteúdo desta página não é a *existência* da obstrução. É:

1. o **valor exato** 3, via o colapso `4 -> 3` pelo núcleo do funcional soma;
2. a **constância em `n`**;
3. a **localidade** nos cantos;
4. a **tightness**: não há outras obstruções.

O enquadramento geral tem nome na literatura: **Hamilton space**. Ver
[[Prior-Art]].

## Q(n) = 3

Provado estruturalmente para `n >= 6`, via 4 representantes independentes mais
o kernel da soma, e verificado computacionalmente em `n = 8` (rank 102) e `n =
12` (rank 294).

A prova geral depende da **conexidade do bulk** (o tabuleiro menos os cantos),
demonstrada por indução `n -> n+2` com bases em `n = 6` e `n = 7`, e
verificada computacionalmente até `n = 30`.

## Q é local, e mora nos cantos

O resultado mais informativo não é o valor 3: é *onde* os 3 graus de
liberdade se perdem.

> **Teorema de localidade.** Para um conjunto `S` de vértices removidos,
> `Q(G_T \ S) = 3` se e somente se `W_corners ⊆ S`.
> Em geral, `Q = max(0, k_deg2 - 1)`, onde `k_deg2` é o número de cantos
> com grau 2.

Verificado em `n ∈ {6, 8, 10}`. Existem 24 configurações de canto fixas para
`n >= 6`.

Em palavras: os 4 cantos do tabuleiro têm grau 2 no grafo do cavalo, e
**arestas incidentes a vértices de grau 2 são obrigatórias** em qualquer tour.
Cada canto assim força uma restrição; quatro cantos forçam quatro, mas uma
delas é dependente das outras, daí `4 - 1 = 3`.

## Identidade fechada

A reformulação por coset dá

```
Q = 8 - c(Punc)
```

onde `c(Punc)` é o número de componentes do grafo perfurado. Essa forma evita
o quociente `F_2^E / row(∂_1)` e é a que foi levada para o Lean.

## Tightness

`Q` diz quanto se perde. **Tightness** é a afirmação de que não se perde mais
nada: `rank(Ham) = beta_1 - Q`.

- Reduzida a: tightness `<=>` `Z_bulk ⊆ Span(Ham)`, com `dim Z_bulk = beta_1 - 4`.
- **Provada por certificados** em `n ∈ {6, 8, 10, 12, 14}` e nos retângulos
  `6×7` e `6×8`, exibindo pares de tours cujo XOR é um hexágono do bulk, cada
  certificado auditado individualmente (grau 2 em todo vértice, existência de
  cada aresta, componente única) e **sem amostragem**.
- Rank cheio com folga: 101/101 em `n=8`, 185/185 em `n=10`, 293/293 em `n=12`.
  Os 1,7-13,2% de hexágonos sem certificado não são necessários.
- Duas rotas independentes (certificação dirigida e flips hexagonais)
  concordam, e erram de formas diferentes.
- O caso `6×8` é o mais forte: 48 casas, enumeração inviável, e a tightness foi
  estabelecida **sem conhecer o conjunto de tours**, bastou certificar
  hexágonos até o rank atingir `dim Z_bulk = 65`.
- Verificada em todos os 10 pontos da família híbrida toro→plano.

Cuidado metodológico registrado: a verificação empírica da tightness precisa
de `K >= 250k` amostras. Com 20k aparecem falsos contraexemplos.

## Os cantos dão só a cota superior

Este é o ponto que separa o que está entendido do que não está.

Os tabuleiros **5×6** e **5×8** têm exatamente a mesma estrutura de cantos,
quatro vértices de grau 2, bulk conexo, `Q = 3`, e deficit **27** e **9**.

Portanto bulk conexo e `Q = 3` **não implicam** tightness. Os cantos
determinam a cota **superior** (`dim Ham <= beta_1 - 3`); o que separa os
casos
*tight* é a cota **inferior**, `Z_bulk ⊆ Ham`.

Qual parâmetro governa essa cota permanece **em aberto**. As duas hipóteses
naturais (bounding box e comportamento de borda) foram ambas refutadas: as
mesmas formas se realizam em 5×10, na mesma posição. A obstrução do `5×8`
**desaparece ao alongar o tabuleiro**, o que é incompatível com qualquer
explicação por gadget local.

## O que Q não é

A obstrução `rank < beta_1` **não é hamiltonicidade**. O rank do espaço gerado
pelos 2-fatores do 6×6 é o mesmo 42. O que falta aos 2-fatores para serem
tours é conectividade, e conectividade é invisível para o rank.

Ver também [[Topologia]], onde `Q` muda quando a borda desaparece, e
[[Formalizacao-Lean]].

**Código:** `experiments/06_invariante_Q/`
