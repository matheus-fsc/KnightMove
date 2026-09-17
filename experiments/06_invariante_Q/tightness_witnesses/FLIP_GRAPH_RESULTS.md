# Grafo de flips hexagonais — resultados

Scripts: `flip_graph.py` (Fase 1, n=6 exaustivo), `flip_graph_n8.py` (Fase 2,
bola de BFS, n ∈ {8,10,12,14}). Dados em `data/flip_graph_*.json`.

## O critério

Um flip `τ → τ ⊕ C` (C hexágono do bulk) é válido se e somente se

1. `τ ∩ C` é um dos **dois emparelhamentos alternados** de C — equivalente a
   "τ contém exatamente uma das duas arestas de C em cada vértice de C", que é
   exatamente a condição de `τ ⊕ C` continuar 2-regular; e
2. `τ ⊕ C` é conexo.

O passo 1 é `AND` + duas comparações de inteiro; o passo 2 é um passeio de |V|
passos. **Nenhuma busca hamiltoniana, nenhum gadget W, nenhum corte.**

Comparação de custo com a rota `switcher_search.py`: lá, uma posição que
*falha* custa ≈450 mil nós de busca; aqui, um par (τ,C) que falha custa duas
comparações. É a diferença entre 9.765 s (n=12, certificação por templates) e
1,77 s (n=6 inteiro, exaustivo) / 71 s (n=12, bola de 600 mil tours).

## Fase 1 — n=6, exaustivo (9.862 tours, sem amostragem)

| medida | valor |
|---|---|
| tours distintos como conjunto de arestas | 9.862 (0 duplicados) |
| hexágonos do bulk | 532 |
| pares (τ,C) com emparelhamento ok | 208.664 |
| — rejeitados por desconexão | 103.992 (49,8 %) |
| — que caem fora dos 9.862 | **0** |
| grau de flip min / médio / máx | **2** / 10,61 / 20 |
| tours isolados (grau 0) | **0** |
| **componentes conexas** | **1** |
| diâmetro do grafo de flips | **10** (excentricidades 9–10, double sweep = 10) |

`pairs_external = 0` é um teste de consistência: todo flip válido pousa dentro
do conjunto dos 9.862, o que confirma que a lista é completa.

**Leitura:** o caso 1 da tabela de interpretação — grafo conexo. Não há
invariante escondido separando componentes em n=6. Um sampler ergódico por
flips hexagonais é viável; e como o diâmetro é 10 com grau médio 10,6 sobre
9.862 vértices, o grafo é praticamente um expansor.

### Os 12 hexágonos nunca realizados — resolvidos

520 dos 532 hexágonos aparecem em algum flip. Os 12 restantes **falham no
passo 1, não no passo 2**: nenhum dos 9.862 tours intersecta qualquer um deles
num emparelhamento alternado (`matching_ok = 0` para todos os 12). A obstrução
é de grau/paridade, não de conectividade. Não é aresta morta: nenhuma das 80
arestas está em zero tours.

**Decomposição D₄.** Os 12 se organizam em exatamente **duas órbitas fechadas,
4 + 8** — nenhuma órbita parcial. A de tamanho 4 é o anel compacto (bbox 3×4,
com simetria de ordem 2, logo 8/2 = 4); a de tamanho 8 é uma forma alongada de
bbox 4×5. Entre os 12 há 10 formas de translação distintas.

```
 órbita 4 (anel, bbox 3x4)      órbita 8 (bbox 4x5, uma das formas)
        .##.                            .#...
        #..#                            ...#.
        .##.                            #.#..
                                        ..#.#
```

**São efeito de tamanho finito.** Colocando as 10 formas no bulk de tabuleiros
maiores e procurando um flip válido:

| n | translados no bulk | **flip-realizáveis** | formas cobertas |
|---|---|---|---|
| 8  | 212 | 175 | **10/10** |
| 10 | 440 | 195 | **10/10** |
| 12 | 748 | 218 | **10/10** |

Cada positivo é um certificado exato — o par (τ, τ⊕C) com ambos hamiltonianos.
Os translados não cobertos são resíduo do truncamento da bola de busca, não
falhas: a bola realizou só uma fração dos hexágonos totais nesses tamanhos.

**Isto inverte a leitura.** A conclusão anterior — "a indução tem de ser pelo
span, porque o enunciado forte falha" — vale **só para n = 6**. Se o enunciado
forte valer para n ≥ 8, a rota fácil de indução reabre. Ver a seção seguinte.

## Fase 2 — n ∈ {8,10,12,14}, bola de BFS

BFS por flips a partir de um tour semente, sem base de tours. A bola é
truncada por um teto explícito; o truncamento **não enfraquece nada**, porque
cada certificado produzido é verificado exatamente — a bola é só o dispositivo
de busca.

| n | hex. bulk | teto da bola | camadas | grau min/méd/máx (bola) | grau min/méd/máx (Warnsdorff) | hex. realizados | rank | dim Z_bulk | tempo |
|---|---|---|---|---|---|---|---|---|---|
| 8  | 2.264  | 200.000 | 1, 53, 1.606, 37.941, 160.401 | 30 / 50,9 / 75 | 27 / 47,9 / 68 | 1.870 | **101** | 101 | 18 s |
| 10 | 5.088  | 250.000 | 1, 99, 5.633, 244.311 | 75 / 96,5 / 121 | 74 / 103,8 / 136 | 3.220 | **185** | 185 | 25 s |
| 12 | 9.000  | 600.000 | 1, 186, 19.389, 580.435 | 159 / 185,8 / 220 | 143 / 176,8 / 209 | 3.987 | **293** | 293 | 71 s |
| 14 | 14.000 | 900.000 | 1, 263, 38.015, 861.738 | 226 / 258,8 / 298 | 223 / 266,2 / 300 | 4.956 | **425** | 425 | 145 s |

- **Crescimento exponencial** em todos os n: a camada 2 é ~30× a camada 1, a
  camada 3 é ~20× a camada 2. Nenhum platô. Evidência a favor de
  conectividade, não prova dela.
- **Grau 0 não aparece em nenhuma amostra**, nem na bola nem na amostra
  independente de Warnsdorff. Achar grau 0 seria refutação conclusiva; não
  achar não é confirmação.
- As duas amostras têm vieses **diferentes** (a bola é enviesada para grau
  alto por construção; Warnsdorff é enviesado de outro jeito) e mesmo assim as
  distribuições ficam a ≈8 % uma da outra. Nenhuma amostra de Warnsdorff caiu
  dentro da bola em n ≥ 10 (0/500, 0/300, 0/200), o que mostra que a bola é
  uma fração desprezível do espaço — como esperado.

## O enunciado forte em n=8 — VALE

`n8_strong_statement.py`. A pergunta que os 12 de n=6 deixaram aberta:

> todo hexágono do bulk é flip-realizável?

Falsa em n=6 (520 de 532). **Em n=8 é verdadeira: 2.264 de 2.264**, todos com
certificado exato (o par τ, τ⊕C com ambos hamiltonianos).

| via | quantidade |
|---|---|
| bola de BFS, 2M tours | 2.198 |
| busca restrita exaustiva no resíduo | 62 |
| bola de BFS com **outra semente**, 1,56M tours | 4 |
| **não flip-realizáveis provados** | **0** |

### Uma lição de método

Os 4 últimos resistiram a buscas restritas exaustivas com orçamento de **60
milhões de nós cada** — e caíram em **110 s** numa bola de BFS com semente
diferente. Não é acaso: os 4 são todos interiores (nenhum toca canto), e é
justamente aí que a busca restrita perde o poder de poda, porque não há
vértice de grau baixo perto das 3 arestas forçadas. **Para certificados
positivos, a bola domina a busca dirigida; a busca restrita só vale a pena
quando se precisa de um negativo conclusivo.**

### O que isso reabre

A rota fácil de indução volta a ser candidata para n ≥ 8: o enunciado
uniforme pode ser o forte ("todo hexágono do bulk é flip-realizável"), que é
muito mais fácil de induzir que o de span. n = 6 é a única exceção conhecida,
e a Fase 1 mostrou por quê — os excepcionais de lá são efeito de tamanho
finito, e a espécie sem explicação geométrica (a órbita do anel, no fundo do
bulk) já evapora em 6×7.

**Dois pontos, não uma lei.** Falso em n=6, verdadeiro em n=8. Falta n=10 em
diante.

## O limiar é na largura, não no tamanho — 6×8 decide

Os dados de exceções (6×6 → 12, 6×7 → 4, 8×8 → 0) eram compatíveis com duas
leituras: limiar em `n` (o tabuleiro cresce) ou em `min(n,m)` (o lado estreito
cresce). Elas só se separam num tabuleiro com um lado 6 e o outro ≥ 8.

**6×8: as 4 exceções persistem, provadas.** Uma por canto, todas com busca
exaustiva (`flip=False, exaustivo=True`), ~17 min cada:

```
(1,4)(2,6)(3,4)(3,5)(4,7)(5,5)   (0,5)(1,7)(2,4)(2,5)(3,6)(4,4)
(0,2)(1,0)(2,2)(2,3)(3,1)(4,3)   (1,3)(2,1)(3,2)(3,3)(4,0)(5,2)
```

| tabuleiro | min(n,m) | exceções | como |
|---|---|---|---|
| 6×6 | 6 | 12 | exaustivo (9.862 tours) |
| 6×7 | 6 | **4** | exaustivo (1.067.638 tours) **e** busca por caminhos |
| **6×8** | **6** | **4** | **busca por caminhos, exaustiva** |
| 8×8 | 8 | **0** | bola + busca, todos certificados |

Um lado de comprimento 8 **não** basta: 6×8 tem tantas exceções quanto 6×7.
O que as elimina é o lado estreito passar de 6. Portanto **o limiar é em
`min(n,m)`, não em `n`** — alargar o tabuleiro na direção longa não ajuda.

### 7×8 fecha o valor: o limiar é `min(n,m) ≥ 7`

7×8 (min = 7) tem **zero exceções**, provado. Dos 1.694 hexágonos do bulk,
1.693 foram certificados pela bola e o último **por transporte de simetria**:
ele é invariante por reflexão vertical mas não horizontal, e sua imagem
horizontal já estava certificada — como reflexão é automorfismo do grafo, o
par (τ, τ⊕C) transporta. Verificado, não assumido.

As 4 exceções do 6×8, transportadas ao 7×8 mantendo a posição relativa ao
canto, são **todas realizáveis**. A família de canto morre em min = 7.

| min(n,m) | 6 | **7** | 8 |
|---|---|---|---|
| exceções | 12 (6×6), 4 (6×7), 4 (6×8) | **0** | 0 (8×8) |

**A largura 6 é o único valor ruim conhecido.**

⚠️ Contraste de custo que vale registrar: nesse mesmo hexágono, a busca
exaustiva consumiu **4.781 s e voltou inconclusiva** (estourou 400M nós),
enquanto o transporte por simetria levou milissegundos. Ver `board_symmetry.py`
e §4-bis do `PLAN.md`.

### A reformulação que tornou isso possível

`flip_paths.py`. Se τ e τ⊕C são ambos tours, `F = τ∖C` tem os 6 vértices de C
com grau 1 e todos os outros com grau 2: é uma **união de 3 caminhos disjuntos
cobrindo o tabuleiro**, com extremos exatamente em C, e nenhum caminho
atravessa vértice de C. Sendo σ o emparelhamento induzido pelos caminhos,

> `F ∪ m` é um único ciclo ⟺ `σ ∪ m` é um 6-ciclo (idem para `m′`),

e **só 4 dos 15** emparelhamentos de 6 pontos servem às duas condições (8
servem só a `m`). Ganhos: as 6 arestas de C ficam proibidas em vez de 3;
tocar um vértice de C fecha o caminho na hora; e o teste do flip sai da folha
e vira escolha de σ na raiz.

O ganho não foi velocidade por nó (54k → 92k nós/s, 1,7×) e sim **tamanho da
árvore**: 6×7 esgota com menos de 2M nós, onde a formulação anterior não
terminava.

### Direção do viés do bug que a validação pegou

A primeira versão de `viable` podava demais: ao inundar a partir da cabeça do
caminho, os vértices de C entravam em `seen` sem serem expandidos, e um
extremo de caminho futuro que caísse ali tinha sua própria inundação
**pulada** — de onde vinham "vértices inalcançáveis" que na verdade eram
alcançáveis. Resultado: 43 divergências, **todas** `got False, want True`.

Poda excessiva faz a busca **perder flips que existem**. Logo o bug
**fabrica exceções**, não as apaga:

| afirmação | sobrevive ao bug? |
|---|---|
| "C é flip-realizável" (positiva) | **sim** — o certificado é exibido e verificado |
| "C **não** é flip-realizável" (negativa) | **não** — pode ser artefato da poda |

E a conclusão do 6×8 é justamente da segunda forma. Por isso a validação aqui
é carga estrutural, e o que a sustenta é ela cobrir **os dois sentidos**:

1. n=6 — reproduz 520/532 exato contra os 9.862 tours, ou seja acerta tanto os
   positivos quanto os 12 negativos;
2. 6×7 — confirma um **negativo genuíno** (4 hexágonos provados irrealizáveis
   pela enumeração completa dos 1.067.638 tours) em 17 s contra 742 s.

Sem um negativo verdadeiro no conjunto de validação, o "4 exceções provadas em
6×8" não teria proteção nenhuma. Este é o mesmo teste de direção que
desmontou o `max_w = 20` em §G1.12 — lá o corte tornava a certificação mais
difícil e preservava só as afirmações positivas; aqui vale o mesmo, e a
afirmação de interesse mudou de lado.

## Fase 3 — retângulos, segundo ponto exaustivo

`flip_graph_rect.py` enumera **todos** os tours fechados por backtracking com
poda de conectividade e monta o grafo de flips completo. Validação: reproduz
9.862 em 6×6 (todos os números da Fase 1 batem) e **44.202 em 5×8**, que é o
valor de literatura para esse tabuleiro.

| tabuleiro | tours | hex. bulk | grau min/méd/máx | isolados | **componentes** | rank realizados | dim Z_bulk | déficit |
|---|---|---|---|---|---|---|---|---|
| 6×6 | 9.862 | 532 | 2 / 10,6 / 20 | 0 | **1** | 41 | 41 | **3** |
| 6×7 | 1.067.638 | 828 | 2 / 17,2 / 38 | 0 | **1** | 53 | 53 | **3** |
| 5×8 | 44.202 | 588 | 3 / 14,6 / 33 | 0 | **1** | 41 | 47 | **9** |
| 5×6 | 8 | 264 | 0 / 0 / 0 | **8** | **8** | 0 | 29 | — |

Os três valores de contagem de tours — 9.862 (6×6), 1.067.638 (6×7) e 44.202
(5×8) — batem com os valores de literatura. Déficit exato medido diretamente
do span dos tours em 6×6 e 5×8; em 6×7 ele vem da redução (`Z_bulk ⊆ Span(Ham)`
mais `Q(6,7)=3`, provado para n,m ≥ 6).

### O que os retângulos mostram

**Três pontos exaustivos, todos conexos.** 6×6 (quadrado), 6×7 (lado ímpar) e
5×8 (estreito) dão uma única componente cada, com grau mínimo 2, 2 e 3. Em
6×7 o grafo de flips tem mais de um milhão de vértices e ainda assim nenhum
tour isolado.

**A linha divisória é min(n,m) ≥ 6, exatamente a hipótese do teorema.** Em
6×7 os hexágonos realizados geram Z_bulk (rank 53 = 53) e a tightness com
déficit 3 se certifica. Em 5×8 o span falha: rank 41 contra dim Z_bulk = 47.

### 5×8 — contraexemplo exato à tightness, com Q = 3 e bulk conexo

Este é o achado mais forte do conjunto, e merece enunciado próprio.

Calculando diretamente (mesma rotina de `Q_locality_theorem.compute_Q`,
adaptada a retângulos — `Q` é a dimensão do span dos pares XOR das arestas
obrigatórias módulo `row(∂₁)`):

| | 6×6 | 6×7 | **5×8** | 5×6 |
|---|---|---|---|---|
| β₁ | 45 | 57 | **51** | 33 |
| vértices de grau 2 / arestas obrigatórias | 4 / 8 | 4 / 8 | **4 / 8** | 4 / 8 |
| **Q** | 3 | 3 | **3** | 3 |
| componentes do bulk | 1 | 1 | **1** | 1 |
| dim Z_bulk | 41 | 53 | **47** = β₁−4 | 29 |
| **déficit medido** | 3 | 3 | **9** | — |

Em 5×8 o bulk é conexo, `Q = 3` exatamente, `dim Z_bulk = β₁ − 4` — todas as
hipóteses estruturais que o teorema `Q(n,m)=3` usa estão satisfeitas — e mesmo
assim

$$\text{déficit}(5\times 8) \;=\; 9 \;>\; 3 \;=\; Q(5,8).$$

Como o déficit é medido exaustivamente sobre os 44.202 tours (valor de
literatura), isto é um **contraexemplo exato à tightness**, não uma
insuficiência de amostragem. O enunciado que ele estabelece:

> **bulk conexo e Q = 3 não implicam tightness.**

A hipótese `min(n,m) ≥ 6` faz trabalho além de tudo o que `Q(n,m)=3` precisa.
Isso muda o que a conjectura de tightness pode assumir: não basta reduzir a
`Q`, é preciso uma hipótese geométrica separada sobre a largura.

Os 6 pontos de déficit extras correspondem exatamente às 6 dimensões de Z_bulk
que não estão em Span(Ham) (rank dos hexágonos realizados = 41, contra
dim Z_bulk = 47).

### Os excepcionais em 6×7 — e a separação das duas órbitas

`rect_unrealized.py` reenumera o milhão de tours e isola os não realizados.
Em 6×7 são **4**, todos falhando no emparelhamento (o mesmo modo de falha do
6×6), e formam **exatamente uma órbita fechada** sob o grupo de Klein — que é
o grupo de simetria de um retângulo não quadrado. Uma por canto.

Isso separa as duas órbitas do 6×6 por posição, e a separação é total:

| órbita (6×6) | tamanho | ancorada em canto | fundo do bulk | sobrevive em 6×7 |
|---|---|---|---|---|
| anel compacto, bbox 3×4 | 4 | 0 | **4** | **não** |
| alongada, bbox 4×5 | 8 | **8** | 0 | **sim** (4, uma por canto) |

Ou seja: **os excepcionais do fundo do bulk desaparecem já no primeiro passo
para cima** (6×6 → 6×7), e as 4 formas que sobram em 6×7 estão todas entre as
10 formas do 6×6, todas encostadas num canto. O que resta é efeito de canto,
não de bulk — e por n=8 as formas já são realizáveis (tabela acima).

**5×6 é o extremo.** Só 8 tours fechados (também valor de literatura), e
**todos com grau de flip 0** — 8 componentes isoladas. Dos 48 pares que passam
no emparelhamento, todos os 48 morrem na conexidade. O movimento de flip é
simplesmente fraco demais em tabuleiros estreitos.

Note que nos casos exaustivos medidos `rank(realizados) = rank{τ ⊕ τ′}` (41 em
6×6 e 41 em 5×8). Isso não é coincidência: **se o grafo de flips é conexo,
todo τ ⊕ τ′ é soma de hexágonos realizados ao longo de um caminho**, logo os
dois spans coincidem. A conectividade é o que converte flips em geradores.

## O que isso estabelece, e o que não

### Estabelecido (verificação exata, sem amostragem na lógica)

Todo flip válido dá `C = τ ⊕ (τ ⊕ C)` com ambos os termos hamiltonianos, logo
`C ∈ Span(Ham)`. Como os hexágonos realizados atingem `rank = dim Z_bulk` e
todos vivem em Z_bulk:

> `Z_bulk ⊆ Span(Ham(n))` para **n ∈ {6, 8, 10, 12, 14}**.

Combinado com `reduction` (Lean, sem `sorry`) e `prop:cota-Q` com Q(n)=3
(provado), isso dá `rank(Ham(n)) = β₁(n) − 3` nesses cinco valores de n.

Para n = 8, 10, 12 isso **reconfirma** o que a Fase 2/3 já havia obtido por
construção; n = 6 e n = 14 são novos, e a rota inteira é ~3 ordens de grandeza
mais barata.

### (U2), separadamente

Os hexágonos do bulk — **todos** eles, não só os realizados — geram Z_bulk com
igualdade exata:

| n | 6 | 8 | 10 | 12 | 14 |
|---|---|---|---|---|---|
| # hexágonos do bulk | 532 | 2.264 | 5.088 | 9.000 | 14.000 |
| rank | 41 | 101 | 185 | 293 | 425 |
| dim Z_bulk | 41 | 101 | 185 | 293 | 425 |

Cinco pontos, todos com igualdade. **(U2) segue sem prova para n geral.**

### Não estabelecido

- **Tightness para todo n.** Cinco valores de n não são uma indução — e o
  5×8 mostra que a hipótese sobre os lados não é decorativa.
- **Conectividade do grafo de flips fora dos casos exaustivos.** Exaustivos:
  6×6, 6×7 e 5×8, todos conexos; 5×6 totalmente desconexo. Para os quadrados
  n ≥ 8 só há bolas, que mostram crescimento exponencial local — compatível
  tanto com conectividade quanto com muitas componentes gigantes.
- **O enunciado forte para n ≥ 10.** O estado é: **falso em n=6** (12 exceções,
  caracterizadas — duas órbitas D₄ fechadas, efeito de tamanho finito),
  **verdadeiro em n=8** (0 exceções, 2.264/2.264 certificados), **aberto de
  n=10 em diante**. Dois pontos não são uma lei.
  *Isto substitui a conclusão anterior de que o enunciado uniforme teria de
  ser o de span; aquela leitura valia só enquanto n=6 era o único dado.*
- **Onde exatamente está o limiar em `min(n,m)`.** Que é em `min` e não em `n`
  está decidido (6×8 tem 4 exceções provadas). Falta separar `min ≥ 7` de
  `min ≥ 8`: o teste é **7×8**, já que 7×7 não tem tours fechados.

## Armadilhas respeitadas

- n=6 usa os 9.862 tours exatos, **não amostragem** (`tours_closed_6x6.npy`).
- O único corte é o teto da bola, e ele está reportado em toda linha, com o
  flag `ball_truncated`. Nenhum corte entra no predicado medido: o critério de
  flip testa **todos** os hexágonos do bulk em **todo** tour visitado.
- Grau 0 foi tratado como resultado, não como bug: `flip_degree_min = 2` em
  n=6 é medida exata sobre o conjunto completo.

---

## Fase 4 — "uma única forma basta como gerador?" (`single_shape_flip.py`)

Pergunta testada: restringindo o conjunto de movimentos a **um único tipo de
hexágono** (uma forma, transladada por todo o bulk), o grafo de flips ainda
conecta os 9.862 tours a partir de uma seed?

**Não — e não chega nem perto.** Os 532 hexágonos do bulk se organizam em
**128 formas** de translação / **23 órbitas D₄**.

| restrição | melhor caso | maior componente | % dos tours |
|---|---|---|---|
| **1 forma** (12 hex, 1.004 flips) | 8.868 componentes | **4** | 0,0 % |
| **1 órbita D₄** (64 hex, 3.200 flips) | 6.802 componentes | **51** | 0,5 % |
| todas as 532 (Fase 1) | **1** componente | 9.862 | 100 % |

**Limiar de geração: 101 das 128 formas** (478 dos 532 hexágonos) são
necessárias para o grafo ficar conexo. A conexão chega tarde e de repente —
comportamento de percolação:

| formas | componentes | maior |
|---|---|---|
| 68–80 | 73 | 9.790 |
| 93 | 11 | 9.847 |
| 99 | 6 | 9.855 |
| 100 | 2 | 9.861 |
| **101** | **1** | **9.862** |

**Leitura.** A conectividade do grafo de flips não é propriedade de um
movimento privilegiado: é propriedade *coletiva* de quase toda a família de
hexágonos. Não existe "hexágono semente" que gere o espaço. O platô 68–80
formas com 73 componentes constantes indica que as últimas formas a entrar
carregam quase toda a conexão — vale investigar quais são.

Dados: `data/single_shape_flip.json`.

---

## Fase 5 — taxonomia das 128 formas (`shape_taxonomy.py`, `shape_taxonomy2.py`)

⚠️ **Correção da Fase 4.** O "101 formas" era **artefato da ordem** (formas
ordenadas por número de hexágonos). Não define limiar, e as 27 restantes não
formam classe alguma. Medidas independentes de ordem:

| categoria | definição | contagem |
|---|---|---|
| **não realizadas** | a forma não habilita nenhum flip válido | **0 / 128** |
| **essenciais** | remover a forma do conjunto completo desconecta | **0 / 128** |
| **redundantes** | remover a forma sozinha mantém conexo | **128 / 128** |

**Nenhuma forma é insubstituível.** A conectividade é sustentada com folga
massiva: qualquer forma pode sair sem consequência.

### Mínimo de formas para conectar

| método | formas |
|---|---|
| ordem por tamanho (Fase 4, enviesado) | 101 |
| **guloso por redução de componentes** | **39** |
| cota inferior por contagem (maior cobertura = 1920 tours) | **≥ 6** |

$6 \le \text{mínimo} \le 39$. O guloso mostra transição abrupta:
1 forma → 8.868 comps; 10 → 2.723; 20 → 429; 30 → 45; 35 → 13; **39 → 1**.
Uma cobertura gulosa de todos os 9.862 tours usa 37 formas — próximo de 39,
sugerindo que o gargalo é **cobrir**, não **costurar**.

### Onde vivem os 12 hexágonos nunca realizados

Espalhados por **10 formas**, sempre como **minoria dos translados** —
nenhuma forma está inteiramente morta:

| bbox | mortos / total | leitura |
|---|---|---|
| 3×4 e 4×3 | 2 / 12 cada | o anel compacto (órbita D₄ de tamanho 4) |
| 4×5 e 5×4 (8 formas) | 1 / 5 cada | a forma alongada (órbita D₄ de tamanho 8) |

Confirma a decomposição 4+8 da Fase 1 e reforça que a obstrução é **posicional**
(qual translado), não da forma em si.

Dados: `data/shape_taxonomy.json`.
