# Reformulação por *coset* — resultados

Validação de `COSET_PROMPT.md`. Aritmética exata sobre GF(2), enumeração
**exaustiva** de tours em todo tabuleiro onde ela existe. Nenhuma amostragem.

## Comandos

```bash
cd tightness_witnesses

# varredura c(Punc) — é aqui que moram os controles negativos de T2
python coset_tests.py --scan --out data/coset_scan.json

# T1, T2 (rápidos; incluem 3x3, 3x4, 3x5 = bulk desconexo)
python coset_tests.py --test T1,T2 \
  --boards 6x6,8x8,6x7,5x8,5x6,10x10,3x3,3x4,3x5 --out data/coset_T1T2.json

# T3..T7 com enumeração exaustiva
python coset_tests.py --test T3,T4,T5,T6,T7 --boards 6x6,5x8,5x6 --tours \
  --out data/coset_T3T7.json
python coset_tests.py --test T3,T4,T5,T6,T7 --boards 6x7 --tours \
  --out data/coset_6x7.json          # ~20 min (1.067.638 tours)

# Lean
cd ../knight_tour_lean && lake build KnightTour.CornerFiber
```

Código: `coset_tests.py`. Dados: `data/coset_*.json`. Lean:
`knight_tour_lean/KnightTour/CornerFiber.lean`.

---

## Veredicto por teste

| teste | veredicto | controle negativo |
|---|---|---|
| T1 — `y_i` bem-definido | **PASS** (9 tabuleiros) | n/a (predição é incondicional) |
| T2 — posto dos 4 funcionais | **PASS** (9) | **presente**: 3×3, 3×4, 3×5 |
| T3 — cota `rk(Ham) ≤ β₁−3` | **PASS** (4 exaustivos) | 5×6, 5×8 (cota folgada) |
| T4 — tours na fibra `y=1` | **PASS** (4 exaustivos) | n/a (predição é incondicional) |
| T5 — a equivalência | **PASS** (4 exaustivos) | **presente**: 5×8, 5×6 |
| T6 — a classificação (D) | **PASS** — iso em 6×6 **e** 6×7 | **presente**: 5×8, 5×6 |
| T7 — `Q(5,8)` | **PASS** — e é o resultado em destaque | — |

**As duas ressalvas do prompt se confirmam. Nenhuma foi refutada.** Detalhe
abaixo.

---

## O resultado em destaque (T7)

> **Bulk conexo + `Q = 3` NÃO implicam tightness.**

`5×8`: `c(Punc) = 5` (bulk conexo), `Q(5,8) = 3` pela Definição 2.4 do paper,
`β₁ = 51`, **`rk(Ham) = 42`**, portanto **deficit = 9 ≠ 3 = Q**. Sobre a
enumeração **exaustiva** dos 44.202 tours fechados — nada de amostragem.

`5×6` é ainda mais extremo: bulk conexo, `Q = 3`, `β₁ = 33`, só 8 tours,
`rk(Ham) = 6`, **deficit = 27**.

Isso delimita exatamente o que a hipótese `min(n,m) ≥ 6` faz no paper: ela
**não** é necessária para `Q = 3` nem para conexidade do bulk (ambos valem já
em 3×6, 4×4, 5×5 — ver varredura), e sim para a inclusão `Z_bulk ⊆ Span(Ham)`,
que é a metade *inferior* da conta. `Q = 3` é a cota superior; a hipótese de
largura sustenta a inferior.

Os quatro pontos exaustivos separam-se exatamente por `min(n,m)`:

| board | `min(n,m)` | `Q` | bulk | deficit | tight? |
|---|---|---|---|---|---|
| 6×6 | 6 | 3 | conexo | 3 | **sim** |
| 6×7 | 6 | 3 | conexo | 3 | **sim** |
| 5×8 | 5 | 3 | conexo | **9** | não |
| 5×6 | 5 | 3 | conexo | **27** | não |

`Q` e a conexidade do bulk são constantes nas quatro linhas; só `min(n,m)`
muda, e é ele que acompanha o deficit. Quatro pontos não são uma prova — mas
localizam a hipótese na metade certa do argumento.

---

## T1 — `y_i` bem-definido

Todo `z ∈ Z₁` usa as duas arestas de cada canto ou nenhuma. Verificado sobre
uma **base** de `Z₁`, o que basta: "coef(e₁) = coef(e₂)" é condição linear, e
uma condição linear satisfeita numa base vale em todo o espaço.

`6×6, 8×8, 6×7, 5×8, 5×6, 10×10, 3×3, 3×4, 3×5`: 0 violações, cantos de grau
`[2,2,2,2]`, 8 arestas obrigatórias distintas em todos.

## T2 — posto dos 4 funcionais, e a identidade exata

Não só a predição se confirma; ela vira **identidade fechada**. Sendo
`r = rk{y₁..y₄}` em `Z₁`:

```
dim Z₁      = E − V + c(G)
dim ker(y)  = dim Z_bulk = (E − 8) − V + c(Punc)
⇒  r = 8 + c(G) − c(Punc)          [G conexo:  r = 9 − c(Punc)]
```

Como `c(Punc) ≥ 5` sempre (os 4 cantos ficam isolados), segue `r ≤ 4`, com
igualdade **sse** o bulk é conexo. Verificado em todos os 9 tabuleiros.

### Controles negativos (obrigatórios — e disponíveis sem construção artificial)

| board | `c(G)` | `c(Punc)` | `r` previsto | `r` medido | `dim ker(y)` | `β₁−4` | `Q` |
|---|---|---|---|---|---|---|---|
| 3×3 | 2 | 9 | 1 | **1** | 0 | −4 | 0 |
| 3×4 | 1 | 7 | 2 | **2** | 1 | −1 | 1 |
| 3×5 | 1 | 6 | 3 | **3** | 3 | 2 | 2 |
| 3×6 … 12×12 | 1 | 5 | 4 | **4** | β₁−4 | β₁−4 | 3 |

Nos três primeiros o bulk é desconexo, `r < 4` e `dim ker(y) > β₁ − 4` —
exatamente a predição negativa. **T2 discrimina.**

> **Erro meu, corrigido:** a primeira versão do teste marcou 3×3 como FAIL
> usando `r = 9 − c(Punc)`. O grafo do cavalo 3×3 é **desconexo** (o centro é
> isolado, `c(G) = 2`), e a fórmula certa é `r = 8 + c(G) − c(Punc)`. Não era
> falha da reformulação; era hipótese de conexidade não declarada na minha
> fórmula. Corrigido, 3×3 passa.

## T3 — a cota superior

`dim A = dim Z_bulk = β₁ − r`; `A` é afim e não contém `0`, logo
`dim Span(A) ≤ β₁ − r + 1 = β₁ − 3` quando `r = 4`.

| board | `β₁` | `r` | `Q` (def. do paper) | cota `β₁−3` | `rk(Ham)` |
|---|---|---|---|---|---|
| 6×6 | 45 | 4 | 3 | 42 | **42** (atinge) |
| 6×7 | 57 | 4 | 3 | 54 | **54** (atinge) |
| 5×8 | 51 | 4 | 3 | 48 | 42 (folga 6) |
| 5×6 | 33 | 4 | 3 | 30 | 6 (folga 24) |

E `Q = r − 1` em **todos** os tabuleiros da varredura, inclusive os de bulk
desconexo. Combinando com T2:

```
Q(n,m) = 7 + c(G) − c(Punc)      [G conexo:  Q = 8 − c(Punc)]
```

## T4 — tours na fibra `y = 1`

| board | tours (exaustivo) | fora da fibra | `y` mal-definido |
|---|---|---|---|
| 6×6 | 9.862 | 0 | 0 |
| 5×8 | 44.202 | 0 | 0 |
| 5×6 | 8 | 0 | 0 |
| 6×7 | 1.067.638 | 0 | 0 |

As contagens batem com a literatura (9.862 e 1.067.638 são os valores
conhecidos), o que audita o enumerador.

## T5 — a equivalência

| board | `rk(Ham)` | `β₁−3` | `dim span{τ+τ'}` | `dim Z_bulk` | LHS | RHS | equivale? |
|---|---|---|---|---|---|---|---|
| 6×6 | 42 | 42 | 41 | 41 | ✔ | ✔ | **sim** |
| 6×7 | 54 | 54 | 53 | 53 | ✔ | ✔ | **sim** |
| 5×8 | 42 | 48 | 41 | 47 | ✘ | ✘ | **sim** |
| 5×6 | 6 | 30 | 5 | 29 | ✘ | ✘ | **sim** |

Os dois lados foram computados de forma independente. Em 5×8 e 5×6 as
diferenças `τ+τ'` estão **contidas** em `Z_bulk` (como devem) mas não o
esgotam. A equivalência vale nos dois sentidos, inclusive onde falha.

> **Segundo erro meu, corrigido:** eu estava filtrando a *base* de `Z₁` pelo
> critério `y(z) = 0` para obter `Z_bulk`. A fibra é um **subespaço**, não um
> subconjunto da base escolhida — isso dava `dim Z_bulk = 31` em 5×8 (o
> correto é 47) e fazia `Ham_inside_preimage` sair falso. Corrigido com
> eliminação simultânea sobre `F₂⁴` e `F₂^E` (`zbulk_basis`); os valores
> passaram a bater com a fórmula independente `(E−8) − V + c(Punc)`.

## T6 — a classificação (D) — o teste central

### 6×6 — **isomorfismo confirmado**

`dim Z₁ = 45`, `rk(Ham) = 42`, `dim Z₁/Span(Ham) = 3`. O mapa induzido por `y`
é bem definido (`y(Span Ham) = {0000, 1111}`), tem núcleo **0**, e vai sobre
`F₂⁴/{0,1}` de dimensão 3. **É isomorfismo.** A classe de um ciclo depende só
do seu padrão de cantos, módulo complementação.

### 6×7 — **isomorfismo confirmado** (segundo caso positivo, exaustivo)

`dim Z₁ = 57`, `rk(Ham) = 54` sobre **todos** os 1.067.638 tours,
`dim Z₁/Span(Ham) = 3`, `y(Span Ham) = {0000, 1111}`, núcleo do mapa induzido
**0**. Isomorfismo sobre `F₂⁴/{0,1}`. Que a confirmação venha em geometria
não-quadrada é o ponto: (D) não é acidente do 6×6.

Os 3 ciclos proibidos (de `forbidden_cycles_6x6/`, lidos por rótulo de aresta,
conferidos como ciclos de grau par) têm padrões — bits na ordem
`(A6, F6, A1, F1)`:

| ciclo (aresta-fonte) | padrão `y` |
|---|---|
| `B6-D5` | `1000` |
| `F6-D5` | `1100` |
| `E3-F1` | `1001` |

Exatamente os previstos. Independentes módulo `1` (posto 4 junto com `1111`) e
independentes no quociente (`rk(Ham ∪ {3 ciclos}) = 42 + 3 = 45`).

### Controle negativo — 5×8 (deficit 9)

`dim Z₁/Span(Ham) = 9 ≠ 3`, então `y` **não pode** ser isomorfismo — e não é.
O que sobra:

- núcleo do mapa induzido tem dimensão **6** (= 9 − 3), como previsto;
- os 6 geradores extras estão **todos dentro de `Z_bulk`** (verificado);
- padrão de canto de todos eles: **`0000`** — invisíveis a `y`.

Predição do prompt confirmada literalmente. Mesmo padrão em 5×6: 24 dimensões
extras, todas em `Z_bulk`, todas com padrão `0000`.

---

## As duas ressalvas — confirmadas

**1. (C) não é mais barata que a prova existente. CONFIRMADO.**
A independência dos 4 funcionais é `r = 4`, e T2 mostra que
`r = 8 + c(G) − c(Punc)`. Logo `r = 4 ⟺ c(Punc) = 5 ⟺ bulk conexo` — o mesmo
ingrediente geométrico do `Q(n)=3` atual (`bulk_connected_general`). A
reformulação **não elimina** a conexidade do bulk; ela a reexpressa.

O que ela ganha é real, mas é de exposição: a cota superior sai sem tocar no
quociente `F₂^E/row(∂₁)`. Ver a seção Lean.

**2. (D) é consequência da tightness. CONFIRMADO.**
(D) vale em 6×6 (tight) e falha em 5×8 e 5×6 (não-tight), com o modo de falha
exatamente previsto: as dimensões excedentes vivem em `Z_bulk` e têm padrão de
canto trivial. (D) não é derivação independente de nada.

---

## Parte Lean — `KnightTour/CornerFiber.lean`

Compila **sem `sorry`**, integrado à raiz `KnightTour.lean`.

| declaração | axiomas |
|---|---|
| `mandEdge1_ne_mandEdge2` | `propext, Classical.choice, Quot.sound` |
| `mem_corner_iff` | `propext, Classical.choice, Quot.sound` |
| `filter_mem_corner` | `propext, Classical.choice, Quot.sound` |
| `bd_apply` | `propext, Classical.choice, Quot.sound` |
| `bd_at_corner` | `propext, Classical.choice, Quot.sound` |
| `cornerVal` | `propext, Classical.choice, Quot.sound` |
| `cycle_corner_eq` | `propext, Classical.choice, Quot.sound` |
| `mem_ZBulk_iff_cornerVal_zero` | `propext, Classical.choice, Quot.sound` |
| `hitsMandatory_of_fiberOne` | `propext, Classical.choice, Quot.sound` |
| `add_mem_ZBulk_of_fiberOne` | `propext, Classical.choice, Quot.sound` |
| `span_le_sup` | `propext, Classical.choice, Quot.sound` |
| `finrank_le_of_fiberOne` | `propext, Classical.choice, `**`Lean.ofReduceBool`**`, Quot.sound` |

O `Lean.ofReduceBool` **não é introduzido aqui**: entra por `dim_ZBulk`, que
já dependia dele, assim como `reduction` e `Q_abstract_eq_three`. A rota do
coset não piora nem melhora a base axiomática.

### O que foi formalizado

- `cornerVal n hn k : (KEdge n → F) →ₗ[F] F` — o funcional `y_k`.
- `cycle_corner_eq` — para `z ∈ Zcyc n`, `z(e₁(c_k)) = z(e₂(c_k))`. É a
  boa-definição de (A), e sai de `bd_at_corner`, que por sua vez sai de
  `corner_adj_iff` (grau 2) transportado de vértices para arestas
  (`mem_corner_iff`).
- `mem_ZBulk_iff_cornerVal_zero` — a fibra `y = 0` **é** `ZBulk n`. Fecha o
  encaixe com `ZBulk_eq`, como o prompt pedia.
- `add_mem_ZBulk_of_fiberOne` — a fibra `y = 1` é afim.
- `span_le_sup` / **`finrank_le_of_fiberOne`** — a cota
  `rk ≤ β₁ − 3` para qualquer conjunto de ciclos na fibra `y = 1`.

### Valor honesto disto

`finrank_le_of_fiberOne` dá **a mesma cota** que `Q_abstract_eq_three` já dava
sem `sorry`. **Nada fica provado que já não estivesse.** O ganho é que a cota
agora se obtém a partir de `dim_ZBulk` + grau 2 dos cantos, **sem** construir o
quociente `F₂^E/row(∂₁)` — que é o objeto cuja exposição em prosa foi criticada
por confusão com `H₁/row(∂₁)` (mal definido: o espaço de cortes não é subespaço
do de ciclos). A prova nova tem ~200 linhas contra o desenvolvimento de
`QAbstract.lean`.

### O que **não** foi formalizado, e por quê

- `tour_corner_functional_eq_one` (T4 em Lean) — **não é enunciável**: `Ham`
  não existe como definição em Lean neste projeto. O que existe é o predicado
  abstrato `HitsMandatory`, e provei a ponte
  `hitsMandatory_of_fiberOne : z ∈ Zcyc → CornerFiberOne z → HitsMandatory z`.
  Aplicar isso a tours de verdade continua bloqueado pelo mesmo motivo de
  sempre — `reduction` também nunca foi aplicado.
- `rank_Ham_le` com `Ham` literal — mesma razão. `finrank_le_of_fiberOne` é a
  forma abstrata, parametrizada por um conjunto `S` de ciclos na fibra.

Seguindo a convenção do projeto, não há `theorem ... := sorry` no módulo.

---

## Padrão de honestidade

- **Isto validou uma reformulação, não provou um teorema novo.** T1–T6 apenas
  reescrevem, de forma mais legível, conteúdo que o paper já tinha.
- **A única afirmação genuinamente nova é T7**, e ela é *negativa*: bulk conexo
  e `Q = 3` não bastam para tightness.
- **Verificado para** — `T1/T2`: `{6×6, 8×8, 6×7, 5×8, 5×6, 10×10, 3×3, 3×4,
  3×5}`; varredura `c(Punc)`/`Q`: todos `R×C` com `3 ≤ R ≤ C ≤ 12`; `T3–T7`:
  `{6×6, 6×7, 5×8, 5×6}` por enumeração exaustiva. **Não provado para todo
  `n`** — exceto a identidade `r = 8 + c(G) − c(Punc)`, que é derivação
  algébrica e vale em geral (dadas as hipóteses de cantos de grau 2 e 8 arestas
  obrigatórias distintas, ambas verificadas caso a caso acima).
- Dois erros meus de implementação foram encontrados **pelos próprios testes** e
  estão registrados acima (T2/3×3 e T5/base-vs-subespaço), não silenciados.
