# Como este projeto chegou onde chegou

Este documento existe porque a ordem em que as ideias apareceram explica o
repositório melhor do que qualquer organização temática. Nada aqui foi
planejado: cada etapa é a consequência de a etapa anterior ter travado.

## 1. A inconformidade inicial

O problema do passeio do cavalo apareceu primeiro numa disciplina de Projeto e
Análise de Algoritmos, resolvido por backtracking ou por guloso (Warnsdorff).
A pergunta que ficou foi: um problema clássico, da mesma família do caixeiro
viajante, realmente se resolve só assim?

O objetivo declarado era enumerar todas as soluções, não achar uma.

## 2. Backtracking, podas, e uma parede

Meses de backtracking com poda, poda inversa e ordenação por grau. O único
conceito estrutural que saiu daí foi a **paridade das cores do tabuleiro**:
o cavalo alterna cores a cada movimento, então um caminho aberto entre duas
casas da mesma cor é impossível num tabuleiro com número par de casas.

Isso corta casos inteiros de graça, mas não muda a classe de complexidade. O
8×8 continuou inalcançável por enumeração. O projeto parou aqui por um tempo.

Ver `experiments/01_backtracking/`.

## 3. Simetria D4

Ao retomar, a primeira melhoria real: o grafo do cavalo no tabuleiro quadrado é
invariante sob o grupo diedral D4. Trabalhar com um representante por órbita
reduz o espaço por um fator de até 8.

Fator constante, mas foi o que tornou a amostragem em tabuleiros grandes viável
na prática. Ver `experiments/02_simetria_d4/`.

## 4. Fixar início e fim

A ideia que o autor chamou de "paridade" (nome próprio, provavelmente já existe
na literatura com outro nome): em vez de enumerar o tabuleiro inteiro, fixar um
par (casa inicial, casa final) e amostrar. Isso compacta o espaço de busca o
suficiente para procurar padrões, e permite construir podas específicas contra
dead ends, já que o destino é conhecido.

Foi o que tornou possíveis as campanhas de amostragem das Fases A, B e C
(8×8 e 10×10). Ver `experiments/03_amostragem/`.

## 5. Os loops

A descoberta que mudou o projeto veio de desenhar a árvore de recursão à mão:
partindo da casa (0,0), abrindo os ramos, e percebendo que **existem caminhos
que colidem** — que chegam à mesma configuração por rotas diferentes.

Esses pontos de colisão foram apelidados de "loops".

## 6. Os loops são o espaço de ciclos

Contando quantos loops existem em um tabuleiro, o número bateu exatamente com
uma quantidade conhecida da teoria dos grafos: a dimensão do espaço de ciclos,

```
beta_1 = |E| - |V| + 1
```

Ou seja: a estrutura descoberta empiricamente era `H_1(G; F_2)`, o primeiro
grupo de homologia do grafo sobre GF(2). Um tour é um vetor em `F_2^E`; o
conjunto dos tours vive num coset do espaço de ciclos; cada ciclo fundamental é
um grau de liberdade.

Isso reformula o problema: em vez de construir tours movimento a movimento,
combinar ciclos por XOR. Ver `experiments/04_espaco_ciclos/`.

## 7. Z3 e cláusulas XOR

Com o espaço de ciclos definido, a rota óbvia: modelar no Z3 com as restrições
como cláusulas XOR, e deixar o solver enumerar. A promessa era eliminar os dead
ends que dominavam o custo do backtracking.

Funcionou, mas o speedup foi ≈ 1,00x. Motivo, descoberto depois: DPLL(XOR) já é
técnica padrão nos solvers modernos — não havia ganho a extrair.
Ver `experiments/05_solvers_xor/`.

## 8. Backtracking + espaço de ciclos + union-find

A combinação que acabou sendo a mais rápida: manter o backtracking, mas operar
sobre o espaço de ciclos e detectar sub-ciclos incrementalmente com union-find.
Isso poda o ramo no instante em que ele se fecha prematuramente, em vez de
descobrir no fim.

Resultado: 25,2x sobre o Z3 no 10×10, custo por tour praticamente constante em
`n ∈ {6,8,10,12,14}`, e a razão 2-fatores/tour caindo de 18,1x para 1,00x.
Ver `solvers/incremental_subtour/`.

## 9. A revisão de literatura

E aqui o projeto encontrou sua parede real. A revisão mostrou que:

- a enumeração por composição de segmentos com union-find é o *multi-path
  method* de Rubin (1974), Christofides (1975) e Kocay (1992) — com
  correspondência termo a termo com o código escrito aqui;
- o XOR de ciclos é o *circuit vector space method* de Welch (1966) e
  Mateti-Deo (1976), e **o resultado negativo já era teorema em 1976**: a
  razão tours/2-fatores tende a zero, logo enumerar por uniões de circuitos é
  desperdício assintótico;
- ambos já tinham sido aplicados **ao grafo do cavalo especificamente** — o
  6x6 com 9862 tours é linha de tabela publicada;
- flips por face sobre tours são a *Z-transformation* de Sheffield (2000);
- DPLL(XOR) é padrão nos solvers.

Uma segunda rodada de auditoria, já na fase de escrita, encontrou mais:

- o enquadramento inteiro tem nome — **Hamilton space** `C_n(G)`, e o
  "deficit" é sua codimensão. É área ativa (Heinig 2013, Hou-Yin 2025,
  Hefetz-Krivelevich 2025, Christoph-Nenadov-Petrova 2024);
- **`deficit > 0` também é prior art**: Heinig observou que `delta(G) >= 3` é
  necessário para Hamilton-generation, e o cavalo tem cantos de grau 2;
- tours em superfícies já foram classificados por classe de homotopia
  (Watkins 2004, Forrest-Teehan 2015, Forrest-Lague 2024).

Tudo que tinha sido construído já existia. Ver `docs/prior-art.md`.

O que sobreviveu foi mais estreito e mais preciso: não a *existência* da
obstrução, mas seu **valor exato**, sua **constância em `n`**, sua
**localidade** nos cantos, e — nas superfícies — **equidistribuição** em vez
de existência.

## 10. O pivô

Sem avanço algorítmico possível, a pergunta mudou. Em vez de "como enumerar mais
rápido", passou a ser **"o que essa estrutura diz sobre o problema"**.

Durante o trabalho com GF(2) tinha aparecido uma anomalia que nunca foi
explicada: o rank do espaço gerado pelos tours ficava sistematicamente **3
abaixo** de `beta_1`. Sempre 3, em todo tabuleiro testado.

Esse deficit virou o objeto central.

## 11. Q(n) = 3

O invariante foi batizado `Q`. O que se estabeleceu:

- `Q(n) = 3` para `n >= 6`, provado estruturalmente e verificado em `n=8` e `n=12`;
- os 3 graus de liberdade perdidos são **locais**, e vivem nos 4 cantos:
  `Q(G_T \ S) = 3 <=> W_corners ⊆ S`, com `Q = max(0, k_deg2 - 1)`;
- identidade fechada `Q = 8 - c(Punc)`;
- tightness (`rank = beta_1 - Q`) provada por construção em `n ∈ {8,10,12}`.

Ver `experiments/06_invariante_Q/`.

## 12. Topologia

Se `Q` vem dos cantos, mudar a topologia deve mudar `Q`. Mudou:

| superfície | Q |
|---|---|
| plano | 3 |
| cilindro | 0 |
| toro | 0 |
| Klein | 1 se `n+m` par, 0 se ímpar |

`Q` mede a obstrução introduzida pela borda. Sem borda, ela some. Isso deu ao
invariante um conteúdo geométrico, não só numérico. Ver `experiments/07_topologia/`.

## 13. Lean 4

Os resultados que sobreviveram à verificação foram formalizados: `Q(n)=3` e
`Q(n,m)=3`, ambos para `n,m >= 6`, sem `sorry`. Ver `formalization/lean/`.

## 14. Transferência, e o limite do método

A última pergunta foi se o método transfere. Foram testados TSP, pathfinding em
grades, decomposição de caminhos, e busca de sítios ativos em proteínas.

Todos negativos, todos pelo mesmo motivo: **GF(2) captura a condição de grau par
com perfeição, e não captura conectividade**. Onde a solução precisa ser conexa,
o XOR produz candidatos de grau correto que se quebram em componentes.

Essa é a mesma obstrução que limita o método no passeio do cavalo, agora vista
de fora. Ver `experiments/09_cross_domain/`.

## Estado atual

O código e os resultados estão estáveis. O que continua em aberto é onde a
contribuição se posiciona em relação à literatura — a revisão foi feita sem
orientação formal, e é exatamente aí que ajuda externa vale mais.
