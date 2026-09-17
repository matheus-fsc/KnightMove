# Forbidden cycles 6×6 — caracterização das 3 dimensões de H₁ não atingidas por tours

Pipeline para identificar e explicar geometricamente as **três dimensões
do espaço de ciclos `H₁(G, F₂)` do grafo do cavalo 6×6** que **nenhum
dos 9 862 tours hamiltonianos fechados realiza**. Achado herdado de
`6x6_higher_order/`:

```
rank_{F₂}(T) = 42   <   β₁(G) = E − V + 1 = 80 − 36 + 1 = 45
                deficit = 3
```

## 1. Os 3 ciclos proibidos (representantes via spanning tree BFS-0)

Os representantes específicos dependem da escolha da árvore geradora;
**o subespaço quociente `H₁ / Ham`** é canônico. Com BFS a partir do
vértice 0 (A6) e percorrendo as arestas não-árvore na ordem canônica:

| ciclo | \|S\| | tipo  | descrição | obrigatórias em S |
|---|---|---|---|---|
| f₀ | 6 | ciclo simples | hexágono no quadrante NW: A6 → C5 → A4 → B6 → D5 → B4 → A6 | A6-C5, A6-B4 |
| f₁ | 6 | ciclo simples | hexágono no quadrante NE (reflexo de f₀): A6 → C5 → E4 → F6 → D5 → B4 → A6 | A6-C5, A6-B4, F6-D5, F6-E4 |
| f₂ | 8 | ciclo simples | octógono "diagonal": A6 → C5 → E4 → D2 → F1 → E3 → D5 → B4 → A6 | A6-C5, A6-B4, E3-F1, D2-F1 |

Plots: `data/plots/forbidden_cycle_{0,1,2}.png`. Dados:
`data/results/forbidden_cycles.{npy,json}`.

Propriedades comuns: cada f_i é um **ciclo simples** (componente única,
todos os vértices de grau exatamente 2 no subgrafo) — não uma união de
ciclos disjuntos. Todos passam pelos dois cantos superiores A6 e B6/F6.

## 2. Por que esses ciclos são proibidos?

O caderno do prompt sugeria quatro hipóteses. Os testes em
`obstruction_theorem.py` deram:

| H | descrição | status | nota |
|---|---|---|---|
| H1 | paridade bipartida | **REFUTADA** | cavalo é bipartido pela cor xadrez → toda aresta cruza cores; nenhuma das 5 outras colorações testadas (paridade de linha, de coluna, half-board) distingue f_i de tours |
| H2 | corte impar | **INCONCLUSIVA** | 7 825 cortes testados (`\|S\|≤3` e estruturais); todos os que zeram para tours também zeram para f_i — não exclui H2 mas torna improvável que exista um corte "natural" obstruindo |
| H3 | f_i não cobre 36 vértices | tautológico | os 3 representantes cobrem 6, 6, 8 vértices — mas **qualquer** vetor de H₁ que não seja exatamente um tour deixa vértices descobertos, então este teste não distingue proibidos de outros elementos do quociente |
| H4 | 2-fatores expandem para 45 dim | **REFUTADA** | enumerei **36 236 2-fatores** (uniões disjuntas de ciclos cobrindo V); `rank_{F₂}(2-fatores) = 42`, **idêntico ao rank dos tours**. A obstrução não vem da conectividade — vem da exigência local "grau exatamente 2 em todo vértice". |

### O que realmente explica: as arestas obrigatórias

A análise dual (em `obstruction_theorem.py` ao final) constrói
explicitamente **três funcionais lineares φ₀, φ₁, φ₂ ∈ (F₂^80)\*** tais
que:

- φ_i(t) = 0 para todo tour t (e para todo 2-fator)
- φ_i(f_j) = δ_ij  (i.e., `F · D^T = I₃`)

São os "invariantes lineares" que o Z3 NUNCA inferiria por restrições de
grau-2 + sub-tour elimination, e que portanto poderiam ser injetados como
cláusulas.

O detector **mais simples tem suporte 2**:

```
φ₁ : x_{F6-D5} + x_{B3-A1}  ≡ 0  (mod 2)
```

Tanto F6-D5 quanto B3-A1 são duas das **8 arestas obrigatórias**
(P(x_e = 1) = 1 em todo tour). Para um tour: 1 + 1 = 0 ✓. Para f₁ (que
usa F6-D5 mas não B3-A1): 1 + 0 = 1 ✗.

A leitura combinatória correta: cada um dos 4 cantos do tabuleiro (A1,
A6, F1, F6, grau 2 no grafo do cavalo) força suas 2 arestas a estarem
em todo 2-fator. Isso dá 8 arestas obrigatórias e impõe **7 igualdades
F₂-lineares** `x_{e_i} = x_{e_j}` entre elas. Dessas 7, descontando o
que já é redundante com a estrutura local de grau-2 nos vizinhos dos
cantos, **sobram 3 restrições genuínas** — e essas 3 restrições são
precisamente as três dimensões "proibidas".

Os outros dois detectores têm suporte 22 e 20 arestas; seus padrões
incluem o **anel intermediário do tabuleiro** (linhas 4-5, colunas
B-E) e vários pares de obrigatórias — confirmam que as 3 obstruções
estão ligadas à interação dos cantos com o miolo.

Arquivos: `data/results/dual_detectors.{npy,json}`,
`data/results/obstruction_hypotheses.json`.

### Simetria D₄

O subespaço `Forbidden` (dim 3) é **invariante** sob a ação canônica do
grupo diédrico D₄ no quociente — cada uma das 7 simetrias não-triviais
realiza um automorfismo F₂-linear de `F₂³` (matriz 3×3 com det = 1
sobre F₂ em todos os 7 casos). f₂ é **fixo coset-a-coset por toda D₄**:
para qualquer g ∈ D₄, `g·f₂ ≡ f₂ (mod Ham)`. f₀ e f₁ se trocam por
reflexão vertical (`fy`). Logo o quociente carrega a representação
{trivial} ⊕ {standard 2d}, o que é compatível com a interpretação "uma
obstrução simétrica + um par antissimétrico canto-NW ↔ canto-NE".

## 3. Conjectura para 10×10

Em 10×10 (`board_10x10/`) já está medido:

```
V = 100, E = 288, β₁ = 189
rank(T amostrado, 5k tours) = 186
8 arestas obrigatórias (mesmos 4 cantos × 2 arestas cada)
```

A simetria estrutural sugere: o "deficit" deve ser também governado
pelas 8 obrigatórias dos cantos, e portanto **deficit_{10×10} ≈ 3**
também, a menos de novas obrigatórias internas (que não aparecem no
6×6 pela ausência de vértices interiores de grau baixo).

Predição testável: ao amostrar suficientes tours 10×10,
`rank(T_{10×10})` estabiliza em **186 = 189 − 3**, ou seja, o
deficit é o mesmo. Se aparecer deficit maior, será sinal de **novas
obrigatórias estruturais** específicas do 10×10 (descobríveis pela
análise dual).

## 4. Implicação para busca via Z3

Os 3 funcionais `dual_detectors.npy` podem ser injetados como
cláusulas F₂-lineares:

```python
# para i = 0, 1, 2:
xor_sum(x_e for e in supp(φ_i)) == 0
```

Em particular `x_{F6-D5} ⊕ x_{B3-A1} = 0` (φ₁) é uma cláusula **de
duas variáveis** que nenhuma propagação grau-2 + sub-tour elimination
infere diretamente — o Z3 só descobre essa relação depois de ramificar
profundamente. Adicioná-la a priori pode ser barato e útil para
acelerar a busca em tamanhos maiores.

A leitura mais geral: enumerar todos os pares (e_i, e_j) de
obrigatórias e adicionar `x_{e_i} ⊕ x_{e_j} = 0` para todos eles dá
**C(8,2) = 28 cláusulas binárias**. Conta:

| espaço | dim |
|---|---|
| `Span(28 pares)` em F₂⁸⁰ | 7 (= 8 obrigatórias − 1 redundância) |
| `row-space(∂)` (= anulador de H₁) em F₂⁸⁰ | 35 |
| `Span(28 pares) + row-space(∂)` | 38 = `ker(T)` |
| **`Span(28 pares)` mod `row-space(∂)`** | **3** ← deficit |

Ou seja, ao agir como funcionais em H₁, as 28 cláusulas binárias entre
obrigatórias **geram exatamente o quociente dual** `ker(T)/row(∂)` de
dim 3. **Bastam 3 cláusulas binárias bem escolhidas para fechar o
quociente** — qualquer base do espaço {x_{e_i}⊕x_{e_j} : i,j ∈ mand}
modulo row(∂) serve. Esta é a versão mais econômica das obstruções:
três XORs entre arestas obrigatórias bastam.

## 5. Estrutura de arquivos

```
forbidden_cycles_6x6/
├── README.md
├── find_forbidden.py          ← T0 + T1 (rank/sanity + os 3 vetores)
├── geometry_forbidden.py      ← T2 (plot + topologia + D₄ + relações)
├── obstruction_theorem.py     ← T3 (H1..H4 + funcionais duais)
└── data/
    ├── plots/
    │   ├── forbidden_cycle_0.png
    │   ├── forbidden_cycle_1.png
    │   └── forbidden_cycle_2.png
    └── results/
        ├── forbidden_cycles.npy             # 3 x 80, uint8
        ├── forbidden_cycles.json
        ├── base_h1.npy                      # 45 x 80
        ├── base_ham_rref.npy                # 42 x 80
        ├── incidence_vertex_edge.npy        # 36 x 80
        ├── geometry_summary.json
        ├── obstruction_hypotheses.json
        ├── dual_detectors.npy               # 3 x 80
        └── dual_detectors.json
```

## 6. Como reproduzir

```bash
python3 forbidden_cycles_6x6/find_forbidden.py        # ~1 s
python3 forbidden_cycles_6x6/geometry_forbidden.py    # ~3 s
python3 forbidden_cycles_6x6/obstruction_theorem.py   # ~30 s
                                                      # (enumera 2-fatores)
```

## 7. Achados em uma linha

- **rank(T) = rank(2-fatores) = 42**: a obstrução não é a hamiltonicidade.
- **`Forbidden ⊂ H₁/Ham`** é D₄-invariante; sua decomposição é
  {trivial: f₂} ⊕ {2d: span(f₀, f₁)}.
- **Detector mínimo**: `x_{F6-D5} ⊕ x_{B3-A1} = 0` para todo tour;
  obstrução localizada nos pares de arestas obrigatórias dos cantos.
- **C(8,2) = 28 cláusulas binárias** entre obrigatórias geram, como
  funcionais em H₁, um quociente de dim **exatamente 3** mod o anulador
  de H₁ (verificado numericamente). **3 dessas cláusulas bastam para
  fechar o quociente** — 3 XORs entre obrigatórias substituem qualquer
  outra família de obstruções.
