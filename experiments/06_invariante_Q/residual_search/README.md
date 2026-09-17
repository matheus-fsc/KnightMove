# residual_search/: propagação em cascata e busca residual no 6×6

Quantas variáveis de aresta restam realmente livres depois de propagar
**todas** as restrições conhecidas (obrigatórias, pares, triplas, quadras
e XOR duais)? E quanto da contagem dos 9.862 tours fechados sobrevive
a essa propagação?

## Setup

- E = 80 arestas; x_e ∈ {0, 1, FREE}
- 8 obrigatórias: A6-C5, A6-B4, F6-D5, F6-E4, B3-A1, E3-F1, C2-A1, D2-F1
  (todas pares de aresta nos 4 cantos, ambas em deg=2)
- 88 pares excluídos · 1.776 triplas minimais · 17.004 quadras minimais
- 3 cláusulas XOR duais φ₀, φ₁, φ₂ (suportes 22, 2, 20; paridade alvo
  = 0 nos 9.862 tours)

## Regras de propagação (até ponto fixo)

| R  | Disparo                                                                                            | Conclusão                       |
|----|----------------------------------------------------------------------------------------------------|---------------------------------|
| R1 | freq[e] = 1.0 (resp. 0.0)                                                                          | x_e := 1 (resp. 0)              |
| R2 | vértice v com 2 arestas em 1                                                                       | demais incidentes em v := 0     |
| R2 | vértice v com (deg(v) − n_zero) = 2                                                                | as 2 livres restantes := 1      |
| R3 | par excluído (e₁,e₂) e x_{e₁}=1                                                                    | x_{e₂} := 0                     |
| R4 | tripla excluída (e₁,e₂,e₃) com 2 delas = 1                                                         | a terceira := 0                 |
| R5 | quadra excluída (e₁..e₄) com 3 delas = 1                                                           | a quarta := 0                   |
| R6 | cláusula XOR de suporte S e paridade alvo p, com apenas uma livre em S                             | valor determinado pela paridade |

Contradição (e₁ fixada em 0 e 1, ou vértice com 3+ uns, ou par excluído
com ambos 1, etc.) ⇒ levanta erro local (em backtracking: poda).

## Resultado da cascata por nível

| Nível | Regras                         | fixadas = 1 | fixadas = 0 | FREE | iterações |
|------:|--------------------------------|------------:|------------:|-----:|----------:|
|     0 | R1                             |           8 |           0 |   72 |         2 |
|     1 | R1 + R2                        |           8 |           0 |   72 |         2 |
|     2 | + R3 (88 pares)                |           8 |           0 |   72 |         2 |
|     3 | + R4 (1.776 triplas)           |           8 |           0 |   72 |         2 |
|     4 | + R5 (17.004 quadras)          |           8 |           0 |   72 |         2 |
|     5 | + R6 (3 XORs)                  |           8 |           0 |   72 |         2 |

**A propagação em cascata não fixa nenhuma variável adicional** além das 8
obrigatórias dos cantos. A razão é estrutural:

- As 8 obrigatórias estão nas 4 esquinas (cada canto tem deg=2 já
  preenchido). R2 vê esses vértices saturados mas não tem livres
  incidentes para propagar.
- Nenhum dos 88 pares excluídos contém uma obrigatória; idem para
  triplas (precisariam de 2 obrigatórias) e quadras (precisariam de 3).
  Verificado por contagem direta (`grep` por interseção).
- Apenas φ₁ ({F6-D5, B3-A1}) tem suporte contido nas obrigatórias:
  consistente, mas sem livres para fixar.

Visualização: [`data/plots/propagation_cascade.png`](data/plots/propagation_cascade.png)

## Espaço residual

Pós-propagação (idêntico para todos os níveis):

- **n_free = 72** · n_fixadas_1 = 8 · n_fixadas_0 = 0
- Bound ingênuo: **2^72 ≈ 4,72 × 10²¹**
- Ground truth: **9.862 tours** ⇒ densidade ≈ 2,1 × 10⁻¹⁸

Grafo de dependência residual (nó = livre, aresta = constraint que
liga duas livres):

- **1 única componente conexa** com os 72 nós
- O subgrafo do cavalo 6×6 já é conectado; cada vértice acopla suas
  arestas via R2, então o grafo residual herda essa conectividade

Correlações condicionais nos tours reais (entre pares de livres):

| estatística                 | valor  |
|-----------------------------|-------:|
| \|corr\| máx                | 0.771  |
| \|corr\| p99                | 0.477  |
| \|corr\| p90                | 0.177  |
| \|corr\| p50                | 0.051  |
| média                       | 0.081  |
| fração \|corr\| > 0.10      | 23.5%  |
| fração \|corr\| > 0.25      |  6.4%  |

→ existem dependências fortes entre livres (a hipótese de
independência aproximada é falsa: 23,5% dos pares têm |corr| > 0,10).

Visualização: [`data/plots/residual_structure.png`](data/plots/residual_structure.png)

## Enumeração via backtracking com propagação

Como n_free = 72 inviabiliza 2^72 ingênuo, foi rodado backtracking
com R1..R6 a cada nó (vetorizado em NumPy), priorizando arestas com
freq mais decidida (menor min(freq, 1−freq)) e testando primeiro o
valor mais provável.

| métrica                                  | valor       |
|------------------------------------------|------------:|
| nós explorados                           |      29.189 |
| folhas (todas as 72 fixadas)             |      13.422 |
| 2-fatores válidos (= folhas, todas R²)   |      13.422 |
| **tours conexos** (validação final)      |   **9.862** |
| tempo total (NumPy vetorizado)           |     53,4 s  |
| fator de poda vs 2^72                    | ≈ 1,6 × 10¹⁷ |

A contagem coincide exatamente com o ground truth. Diferença
13.422 − 9.862 = **3.560 2-fatores multi-ciclo** que satisfazem todas
as restrições locais R1..R6 mas **não são conexos**, ou seja, são
uniões de ciclos curtos no grafo do cavalo. Nenhuma constraint local
do ideal de grau ≤ 4 (nem as XOR duais) elimina esses pseudo-tours.

## Resposta à pergunta central

> Conhecer todas as restrições locais (obrigatórias + 88 pares +
> 1.776 triplas + 17.004 quadras + 3 XORs) **elimina** a busca?

**Não.** A propagação em cascata, partindo apenas das obrigatórias
dos cantos, deixa 100% do interior (72/72 arestas) livre e numa
única componente acoplada. O ideal local de grau ≤ 4 + as XORs
duais conhecidas **caracterizam exatamente** os 2-fatores do tabuleiro
6×6 que respeitam todas as restrições conhecidas, mas esse conjunto
tem **13.422 elementos**, dos quais 3.560 não são tours porque são
desconexos.

A **conectividade** (= ser um único ciclo de 36 vértices) é o
obstáculo restante. Ela é uma propriedade global do subgrafo de
arestas em 1 e não admite expressão como combinação local fixa de
arestas, explica a diferença `2-fatores_locais − tours = 3.560`.

Para o algoritmo proposto "borda fixada + interior combinatório":

- **A decomposição não é trivialmente realizável**: o grafo de
  dependência residual é uma única componente que mistura interior
  e borda (a única propagação que funcionaria seria a partir de uma
  escolha de borda específica, mas isso já é busca, não propagação).
- **Mas funciona com branching guiado**: o backtracking acima é
  precisamente "decidir uma livre, propagar, repetir", chega a
  9.862 em 53 s sem invocar Z3, com 29.189 nós explorados, todos
  podados por R1..R6 + viabilidade de R2.

## Implicação

As restrições locais conhecidas são **necessárias e quase suficientes**:
elas reduzem 2^72 a 13.422 candidatos (poda de ≈10¹⁷). Para fechar
a conta com 9.862 tours, ou:

1. acoplar a conectividade ao processo de busca (já feito como
   filtro de folha aqui, mas idealmente como constraint global tipo
   "no-subtour" de TSP); ou
2. caracterizar os 3.560 2-fatores extras com novas
   constraints, provavelmente associadas aos **3 ciclos proibidos**
   já identificados em `forbidden_cycles_6x6/` mas que aqui só apareceram
   como φ₀, φ₁, φ₂ duais sobre H¹ (cláusulas XOR de paridade).

A queda mais barata é (1): backtracking + R1..R6 + filtro de
conectividade chega exatamente a 9.862 em ~minuto, sem solver SAT.

## Arquivos

- `propagation_engine.py`: motor R1..R6 (NumPy vetorizado), gera
  `data/propagation_log.json`
- `residual_analysis.py`: estrutura do residual e correlações,
  gera `data/residual_variables.json`
- `enumeration_test.py`: backtracking com propagação, gera
  `data/residual_enumeration.json`
- `plots.py`: gera `data/plots/{propagation_cascade,residual_structure}.png`
