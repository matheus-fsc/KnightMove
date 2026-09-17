# Gates 1 e 2 — resultados

**Data:** 2026-08-07 · **Plano:** `PLAN.md` · **Scripts:** `gate1_templates.py`, `gate2_certify_templates.py`
**Dados:** `data/gate1_templates.json`, `data/gate1_greedy_shapes.json`, `data/gate2_certify_templates.json`

---

## Gate 1 — templates de translação · **PASSOU**

### G1.1 Classificação das formas

| m | hexágonos | formas (translação) | formas (translação + D₄) | translações por forma |
|---|---|---|---|---|
| 8 | 2.352 | **136** | **25** | 30 = (m−2)(m−3) |
| 10 | 5.176 | **136** | **25** | 56 = (m−2)(m−3) |
| 12 | 9.088 | **136** | **25** | 90 = (m−2)(m−3) |

O número de formas **estabiliza** — é o mesmo em m ∈ {8,10,12}. Era a condição
necessária para existir template finito, e ela vale.

Além disso, as translações de **uma única forma** são linearmente independentes
entre si: rank = número de translações, em todos os m testados.

### G1.2 Família geradora mínima

**k = 5 formas, e 5 é exatamente mínimo.**

- Guloso: 5 formas geram `Z₁(G_m)` para m = 8, 10, 12 (rank 105/105, 189/189, 297/297).
- Cota inferior por contagem: `dim Z₁ / translações por forma` = 105/30 = 3,5 · 189/56 = 3,375 · 297/90 = 3,3 ⇒ k ≥ 4.
- **k = 4 esgotado exaustivamente em m = 8**: todas as C(136,4) = 12.919.190 combinações testadas, nenhuma gera. Logo k = 5.

As **mesmas 5 formas** servem nos três tabuleiros:

| # | forma (offsets de célula) | caixa | translações (m=12) |
|---|---|---|---|
| T0 | (0,0)(0,1)(1,2)(1,3)(2,0)(2,1) | 3×4 | 90 |
| T1 | (0,0)(0,2)(1,0)(1,2)(2,1)(3,1) | 4×3 | 90 |
| T2 | (0,1)(0,2)(1,0)(1,3)(2,1)(2,2) | 3×4 | 90 |
| T3 | (0,0)(0,2)(1,2)(2,1)(2,3)(3,1) | 4×4 | 81 |
| T4 | (0,2)(0,3)(1,0)(1,1)(2,2)(2,3) | 3×4 | 90 |

Todas cabem em **caixa 4×4**.

> ⚠️ **ALVO ERRADO — ver G1.3.** Este `k = 5` é para `Z₁(G_m)` com **todos** os
> hexágonos, inclusive os que tocam canto. O alvo relevante é `Z_bulk`, onde
> hexágonos de canto não existem — e aí T0…T4 **não bastam**.

---

## Gate 2 (parte 1) — a certificação depende só da forma?

Para cada translação de cada template mediu-se certificação × `inset`
(distância da caixa à borda do tabuleiro).

### n = 10

| template | translações | certificadas | inset 0 | inset 1 | inset 2 | inset 3 | d₀ |
|---|---|---|---|---|---|---|---|
| T0 | 54 | 52 (96,3%) | 22/24 | 18/18 | 10/10 | 2/2 | 1 |
| T1 | 54 | 52 (96,3%) | 22/24 | 18/18 | 10/10 | 2/2 | 1 |
| T2 | 56 | 48 (85,7%) | 22/26 | 14/18 | 10/10 | 2/2 | 2 |
| T3 | 48 | 44 (91,7%) | 19/23 | 16/16 | 8/8 | 1/1 | 1 |
| T4 | 54 | 52 (96,3%) | 22/24 | 18/18 | 10/10 | 2/2 | 1 |

### n = 12

| template | translações | certificadas | inset 0 | 1 | 2 | 3 | 4 | d₀ |
|---|---|---|---|---|---|---|---|---|
| T0 | 88 | 86 (97,7%) | 30/32 | 26/26 | 18/18 | 10/10 | 2/2 | 1 |
| T1 | 88 | 86 (97,7%) | 30/32 | 26/26 | 18/18 | 10/10 | 2/2 | 1 |
| T2 | 90 | 82 (91,1%) | 30/34 | 22/26 | 18/18 | 10/10 | 2/2 | 2 |
| T3 | 80 | 76 (95,0%) | 27/31 | 24/24 | 16/16 | 8/8 | 1/1 | 1 |
| T4 | 88 | 86 (97,7%) | 30/32 | 26/26 | 18/18 | 10/10 | 2/2 | 1 |

**Toda falha ocorre em `inset ≤ 1`.** A partir de `inset ≥ 2` a certificação é
100% nos cinco templates, em n=10 e n=12. Isto é o padrão previsto: certificação
determinada pela forma, com efeito de borda de espessura constante.

### G2.1 — rigidez: é o mesmo **conjunto**, não só a mesma contagem

Contagem igual não é a lei. A lei é: as posições que falham, em coordenadas
relativas ao canto mais próximo, são as mesmas. Assinatura registrada:
`(min(top,bot), lado_vertical, min(left,right), lado_horizontal)` — o lado
importa porque as formas não são D₄-simétricas.

| template | assinaturas de falha | nº | n=10 | n=12 | n=14 |
|---|---|---|---|---|---|
| T0 | (1,0,0,1) (1,1,0,1) | 2 | ✅ | ✅ | ✅ |
| T1 | (0,1,1,0) (0,1,1,1) | 2 | ✅ | ✅ | ✅ |
| T2 | (0,0,0,0) (0,0,0,1) (0,1,0,0) (0,1,0,1) (1,0,1,0) (1,0,1,1) (1,1,1,0) (1,1,1,1) | 8 | ✅ | ✅ | ✅ |
| T3 | (0,0,1,1) (0,1,0,1) (1,0,0,0) (1,1,0,1) | 4 | ✅ | ✅ | ✅ |
| T4 | (1,0,0,0) (1,1,0,0) | 2 | ✅ | ✅ | ✅ |

**Conjuntos idênticos nos três tabuleiros**, com bijeção (cada assinatura ocorre
exatamente uma vez por `n`). Não é coincidência numérica: o efeito de borda é
**rígido**, ancorado nos cantos, e independente de `n`.

> Este é o enunciado de que a indução precisa: a borda deixa de ser "um caso novo
> a cada `n`" e vira um conjunto finito enumerável de 18 posições canto-relativas.

Contraste que valida a métrica: a família de 4 formas de G1.5 **falha** neste
teste (assinatura `(0,0,2,1)` só em n=10). A rigidez é propriedade da família
escolhida, não automática.

### G1.3 — correção: T0…T4 é família para `Z₁`, não para `Z_bulk`

Medindo os ranks contra `dim Z_bulk` (hexágonos que tocam canto **excluídos**):

| conjunto de geradores | n=10 (dim 185) | n=12 (dim 293) | n=14 (dim 425) |
|---|---|---|---|
| todos os hexágonos sem canto | 185 ✅ | 293 ✅ | — |
| só T0…T4, sem canto | 184 | 292 | — |
| **T0…T4 certificados, sem canto** | **183** | **291** | **423** |

O déficit é **exatamente 2 e constante em `n`**, e decompõe-se em duas causas de
tamanho 1 cada: **−1** por restringir às 5 formas, **−1** pelas falhas de
certificação.

### G1.4 — a correção de escopo do interior, medida

Restringindo a `inset ≥ 2` (onde a translação valeria sem ressalva):

| n | rank do interior | `dim Z_bulk` | `dim Z₁(G_{n−4})` | déficit |
|---|---|---|---|---|
| 10 | **45** | 185 | **45** | 140 |
| 12 | **105** | 293 | **105** | 188 |
| 14 | — | 425 | 189 | 236 |

O interior gera exatamente `Z₁` do tabuleiro `(n−4)×(n−4)` — nem uma dimensão a
mais. Em forma fechada:

```
dim Z_bulk(n)   = 3(n−1)(n−3) − 4 = 3n² − 12n + 5
dim Z₁(G_{n−4}) = 3(n−5)(n−7)     = 3n² − 36n + 105
déficit         = 24n − 100
```

**A borda carrega `24n − 100` dimensões, para sempre.** Não é um resto que
some assintoticamente.

### G1.4-bis — a insuficiência era do índice, não da construção

Conclusão anterior — *"o catálogo é necessário e insuficiente por construção"* —
**estava errada**, e o erro foi pressupor índice binário interior/borda. Nesse
índice o catálogo de fato nunca fecha, porque a borda é Θ(n).

Mas as posições de borda são **Θ(n) em número e O(1) em tipo**. O índice certo é
a assinatura de G2.1, `(min(top,bot), lado_v, min(left,right), lado_h)`: capando
as distâncias em ~3 (além disso é interior por translação), o número de classes
de assinatura é **finito**. Cada classe se certifica uma vez no catálogo e vale
para todo `n`.

Indexado por assinatura, o catálogo cobre interior **e** borda, e é finito e
completo. O dimensionamento (passo 2) deve contar **classes de assinatura**, não
configurações perfuradas interiores.

### G1.5 — família de 4 formas para `Z_bulk`: **resultado negativo**

Refazendo o guloso com o alvo certo (`Z_bulk`, sem hexágonos de canto), 4 formas
bastam — as mesmas em n=10 e n=12:

```
S0 (0,1)(0,2)(1,0)(1,3)(2,1)(2,2)     S2 (0,0)(0,1)(1,2)(1,3)(2,0)(2,1)
S1 (0,1)(1,0)(1,2)(2,0)(2,2)(3,1)     S3 (0,0)(1,2)(1,3)(2,0)(2,1)(3,2)
```

Mas **sob certificação essa família colapsa**: rank 173/185 (n=10) e 282/293
(n=12) — déficit 12 e 11, **não constante**. Pior, as assinaturas de falha de S3
**não são rígidas**: `(0,0,2,1)` aparece em n=10 e some em n=12. Família
descartada. Menos formas não é melhor se elas certificam pior.

### G1.6 — família de 6 que funciona · **T5**

Das 136 formas, 24 fecham o déficit 2 de G1.3 (sem filtro de certificação).
Certificando as 10 de menor caixa, várias fecham; a melhor é

```
T5 = (0,1)(0,3)(1,2)(2,0)(2,2)(2,4)      caixa 3×5
```

| n | T5 certificadas | falhas | rank T0…T5 certificados | `dim Z_bulk` |
|---|---|---|---|---|
| 10 | 46/46 | **nenhuma** | 185 | 185 ✅ |
| 12 | 78/78 | **nenhuma** | 293 | 293 ✅ |
| 14 | 118/118 | **nenhuma** | 425 | 425 ✅ |

**T5 certifica em 100% das posições, inclusive encostada na borda** — não tem
efeito de canto algum. E `{T0,…,T5}` restrito às translações **certificadas**
gera `Z_bulk` exatamente, em n = 10, 12 e 14.

> **Alvo corrigido de (U1)+(U2):** 6 formas, todas em caixa ≤ 3×5 ou 4×4, cujas
> translações certificadas geram `Z_bulk`. O efeito de borda vive em 5 delas e é
> um conjunto finito rígido (G2.1); a sexta não tem borda.

### G1.7 — formas livres de borda: **existem, são 42, e geram**

Certificação de **todas** as posições de **todas** as 136 formas em n=10
(≈ 7 mil chamadas, 2617 s):

- **42 de 136 formas são livres de borda** (taxa de certificação 100% em toda
  posição, inclusive encostadas na borda).
- **O span delas gera `Z_bulk` exatamente: rank 185/185.**

Isto **dissolve o conjunto excepcional de 18 posições** de G2.1: existe família
geradora em que toda posição de toda forma certifica. Sem enumeração de borda,
sem classes de assinatura, sem catálogo de borda.

#### O preço: compacidade é o que causa sensibilidade à borda

| lado máx. da caixa | livres de borda |
|---|---|
| ≤ 4 | **0/14 = 0%** |
| 5 | 17/56 = 30% |
| 6 | 24/62 = 39% |
| 7 | 1/4 = 25% |

**Nenhuma** forma de caixa ≤ 4×4 é livre de borda — T0…T4 falham todas. A
hipótese de que a caixa 3×5 do T5 explicaria a propriedade está **refutada**: das
2 formas (3,5), só uma é livre; o padrão é o oposto, formas *espalhadas* é que
são livres.

Explicação consistente com a causa de falha já documentada (canto com < 2
vizinhos restantes): um hexágono compacto perto do canto remove vários vizinhos
do canto de uma vez; um espalhado toca um por vez.

| família livre de borda | formas | rank (dim 185) |
|---|---|---|
| caixa ≤ 4×4 | 0 | 0 |
| caixa ≤ 5×5 | 17 | 176 (falta 9) |
| **caixa ≤ 6×6** | **41** | **185 ✅** |
| caixa ≤ 7×7 | 42 | 185 ✅ |

#### G1.8 — o trade-off dissolve: 18 formas dão as duas coisas

`bbox(W)` medido (ordenação por caixa) sobre as 41 formas livres de caixa ≤ 6×6:
pior caso **8×8**, histograma `{6: 18, 7: 21, 8: 2}`. Isoladamente isso pareceria
exigir bloco D&C 8×8 — bem mais caro que o catálogo 6×6 de 25,6k caminhos que já
existe.

Mas as 18 de `bbox(W) ≤ 6` **bastam sozinhas**:

| família livre de borda | formas | rank (dim 185) |
|---|---|---|
| **`bbox(W) ≤ 6×6`** | **18** | **185 ✅** |
| `bbox(W) ≤ 7×7` | 39 | 185 ✅ |
| `bbox(W) ≤ 8×8` | 41 | 185 ✅ |

> **Melhor cenário disponível, em n=10.** Família de 18 formas com
> **(a)** certificação 100% em toda posição — sem conjunto excepcional de borda —
> e **(b)** gadget dentro de 6×6 — o catálogo D&C existente é do tamanho certo.
> As duas metades de (U1) fecham com o mesmo objeto.

Todas as 18 têm ao menos um lado ≥ 5 (nenhuma compacta), coerente com G1.7.
Lista em `data/borderfree_bboxW6.json`.

#### G1.8-bis — replicação em n=12 · **passou**

Comparação como **conjunto** (não "as 18 sobreviveram?"), que distingue rigidez
de facilidade-crescente:

| | n=10 | n=12 |
|---|---|---|
| formas totais | 136 | 136 (idênticas) |
| livres de borda | **42** | **60** |
| `f₁₀ ∖ f₁₂` | — | **1** |
| `f₁₂ ∖ f₁₀` | — | **19** |
| rank do span das livres | 185/185 ✅ | 293/293 ✅ |

Quase-supraconjunto: ser livre de borda é **mais fácil** em tabuleiro maior, logo
n=10 é o caso mais restritivo e testar o pequeno basta. Não é o caso "cruzado"
que invalidaria a lista.

**As 18 de `bbox(W) ≤ 6`: 18/18 continuam livres em n=12, e geram nos dois:**

| n | rank das 18 | `dim Z_bulk` |
|---|---|---|
| 10 | 185 | 185 ✅ |
| 12 | 293 | 293 ✅ |

**A única forma que regride** — `(0,2)(0,3)(1,0)(1,5)(2,2)(2,3)`, caixa 3×6 — não
está entre as 18, e falha em **1 de 70** posições em n=12. Repetindo essa posição
com orçamento 20× maior (8M nós): continua sem certificado, mas com
`exhausted = False`. Ou seja **não encontrado dentro do orçamento**, *não*
provado inexistente.

#### G1.8-ter — n=14 · a exceção era artefato de orçamento

| n | livres de borda | rank do span | rank das 18 |
|---|---|---|---|
| 10 | 42 | 185/185 ✅ | 185/185 ✅ |
| 12 | 60 | 293/293 ✅ | 293/293 ✅ |
| 14 | **61** | 425/425 ✅ | **425/425 ✅** |

`f₁₂ ⊆ f₁₄` (inclusão estrita, `f₁₄ ∖ f₁₂` = 1) e `f₁₀ ⊆ f₁₄`. E a forma
acrescentada em n=14 é **exatamente** a que "regredia" em n=12 — a única quebra
de monotonia da cadeia era o artefato de orçamento já sinalizado.

Com isso a família é **monótona crescente em `n`**: n=10 é o caso mais
restritivo, e as 18 de `bbox(W) ≤ 6` geram `Z_bulk` em **n = 10, 12 e 14**.

### G1.9 — busca de um critério para "livre de borda"

Uma *lista* de 18 formas tem de ser reverificada a cada `n` e nunca vira teorema;
um *critério* vale para todo `n` por construção. Duas hipóteses testadas em n=10,
sobre as 136 formas:

**(i) Alongamento da caixa — REFUTADO como critério.**

| critério | separa | necessário | suficiente |
|---|---|---|---|
| `max(lado) ≥ 5` | 41% | ✅ | ❌ (42 livres × 80 não-livres) |
| `max(lado) ≥ 6` | 57% | ❌ | ❌ |
| `min(lado) ≥ 5` | 69% | ❌ | ❌ |
| área ≥ 25 | 68% | ❌ | ❌ |

`max(lado) ≥ 5` é **necessário** (nenhuma forma de caixa ≤ 4×4 é livre) mas está
longe de suficiente. Tamanho de caixa não é o mecanismo.

**(ii) Pressão sobre o canto — necessário, não suficiente.**

Pressão = `max` sobre os 4 cantos de `|hexágono ∩ N(canto)|`; cantos têm grau 2.

| pressão máx. | livres | não-livres | % livres |
|---|---|---|---|
| 1 | 42 | 64 | 40% |
| **2** | **0** | **30** | **0%** |

**Pressão 2 ⟹ nunca livre de borda**, sem exceção: se o hexágono pode cobrir os
*dois* vizinhos de um canto, o canto fica com grau 0 e (S3) morre. Necessário,
mas 42/106 entre as de pressão 1 — o discriminante restante está no gadget, não
no hexágono.

**(iii) Causa dominante da rejeição — 90% é o canto.** Sobre todas as posições
que falham em n=10:

| causa | n |
|---|---|
| `canto_grau < 2` | **180** |
| `ham_path` **exausto** (provado inexistente) | 20 |
| cor/paridade | 0 |
| orçamento estourado | **0** |

Zero falhas por orçamento: em n=10, "não é livre de borda" está **provado**, não
limitado por budget.

**(iv) O critério combinatório — necessário, sem busca de caminho.**

Como cantos têm grau exatamente 2, "o canto mantém ≥2 vizinhos" equivale a "o
conjunto removido não toca nenhum vizinho de canto". Isso dá um teste que **não
usa `ham_path` nenhum**:

> **Critério C.** Uma forma passa se, em **toda** posição, existe `W` válido
> (cor/paridade ok) com `V(W) ∖ {polos}` disjunto de `⋃ N(canto)`.

| n=10 | livre | não-livre |
|---|---|---|
| C verdadeiro | **42** | 11 |
| C falso | **0** | 83 |

**Necessário** (zero falsos negativos), separa 91,9%. Os 11 falsos positivos são
exatamente as formas cujas falhas eram `ham_path` exausto (taxas 0,94–0,98).

**(v) O critério é `n`-independente — e a versão anterior deste item estava
contaminada.**

⚠️ **CORREÇÃO.** A primeira medição usava `build_switchers(..., max_w=20)`, que
trunca a enumeração de candidatos `W` em 20 por rotação. Como a ordem de
enumeração depende dos índices de vértice, e esses dependem de `n`, o corte
introduzia dependência em `n` que **não é do critério**. Os números reportados
antes — 53 formas em n=10, 73 a partir de n=12, com "limiar n₀=12" — são
artefatos. Nenhum deles sobrevive.

Sem o corte (`max_w = 10 000`, saturado — `200`, `1 000` e `10 000` dão o mesmo):

| n | formas que passam em C | conjunto |
|---|---|---|
| 10 | **94** | idêntico |
| 12 | **94** | idêntico |
| 14 | **94** | idêntico |
| 16 | **94** | idêntico |

**Não há limiar.** O conjunto é o mesmo desde n=10, o menor tabuleiro par em que
a medição faz sentido. `max_internal = 5` devolve os mesmos 94, logo esse corte
também não morde.

Consequência para o poder discriminante, que cai bastante:

| n=10 | livre de borda | não-livre |
|---|---|---|
| C verdadeiro | 42 | **52** |
| C falso | **0** | 42 |

C continua **necessário** (zero falsos negativos — um predicado mais permissivo
não pode perder nenhum livre-de-borda), mas separa **61,8%**, não os 91,9%
reportados antes.

### G1.11 — Lema de localidade (por que C não depende de `n`)

> **Lema (localidade de C).** Seja `n` par. Para uma forma `S`, a validade de
> `C(S, n)` não depende de `n`.

*Argumento.* As duas cláusulas de C são locais e invariantes por translação:

**(1) A condição de cor.** C exige `cor(v₁) ≠ cor(v₄)` e que o conjunto
permitido seja balanceado em cor. Num tabuleiro `n × n` com `n` par o conjunto
de vértices já é balanceado (`n²/2` de cada cor); remover `R = V(W) ∖ {polos}`
preserva o balanço **se e só se** `R` é balanceado. Ora, `R` é determinado por
`W` a menos de translação, e translação preserva balanço de cor. Igualmente,
`cor(v₁) ≠ cor(v₄)` depende só do deslocamento `v₁ − v₄`. Ambas as cláusulas
são portanto propriedades de `W` como configuração, sem referência a `n`.

**(2) A condição de canto.** `⋃_c N(c)` são 8 vértices em deslocamentos fixos a
partir dos 4 cantos. `V(W)` está contido numa janela de tamanho limitado — os
vértices internos de `P₂`/`P₃` estão a no máximo `max_internal + 1` movimentos
de cavalo de um vértice do hexágono, e o `bbox(W)` medido nunca passa de 8. A
condição `R ∩ ⋃_c N(c) = ∅` depende, portanto, apenas da posição da janela
relativa a cada canto — não de `n`.

De (1) e (2), o conjunto de configurações locais realizáveis é o mesmo para todo
`n` maior que o diâmetro da janela, e o critério estabiliza. ∎

⚠️ **Estatuto.** Isto é um *argumento*, não uma prova formal: falta explicitar a
janela e verificar que toda configuração local realizável num `n` o é em
qualquer outro. O que está **verificado** é a conclusão, para
`n ∈ {10, 12, 14, 16}` com enumeração saturada. A previsão do argumento — que
o limiar seja pequeno — é consistente com não haver limiar nenhum no intervalo
testado.

---

⚠️ **O que a parte 1 não é.** É medida de *existência*, não prova de invariância por
translação. O caminho de (S3) é **global** — cobre o tabuleiro inteiro — logo
não pode ser translation-invariant no sentido ingênuo. O dado diz que a
existência é governada pela forma; a prova ainda precisa de um argumento global
de caminho hamiltoniano. É exatamente o que o Gate 2 propriamente dito (bloco
6×6 perfurado + catálogo D&C) deve fornecer.

---

## Estado dos gates

| Gate | Estado |
|---|---|
| 1 — templates | ✅ **passou**: k=5 mínimo exato, formas estáveis, caixa 4×4 |
| 2 (parte 1) — certificação por forma | ✅ passou: d₀ ≤ 2, idêntico em n=10 e n=12 |
| 2 (parte 2) — premissa do bloco 6×6 (`bbox(W)`) | ✅ **verificada** (G2.2) |
| 2 (parte 2) — catálogo D&C perfurado | ⏳ não executado |
| 3 — (U2) Rota A | ⏳ não iniciado (agora com alvo T0…T4) |
| 0 — Lean da redução | ⏳ não iniciado |

---

## G2.2 — `bbox(W)`: o gadget cabe no bloco 6×6 · **verificado**

Preocupação levantada: `|V(W)| = 8` = 6 vértices do hexágono + 2 internos de
`P₂`/`P₃`; os hexágonos cabem em 4×4, mas os internos entram por movimento de
cavalo e podem sair a distância 2 — `bbox(W)` poderia passar de 6×6. De fato,
com os candidatos ordenados por `|V(W)|` (como em `switcher_search.py`),
aparece `bbox` (6,7) no T2.

Mas há folga: `build_switchers` devolve até 20 candidatos por rotação.
Reordenando os candidatos por **`bbox` em vez de `|V(W)|`** e recertificando
(n=12, `inset ≥ 2`, 145 posições):

| template | cabem em 6×6 | `bbox` observados |
|---|---|---|
| T0 | 30/30 | (4,4) |
| T1 | 30/30 | (4,4) |
| T2 | 30/30 | (4,5), (5,5) |
| T3 | 25/25 | (4,4), (5,5) |
| T4 | 30/30 | (4,4) |
| **total** | **145/145** | máximo **5×5** |

**Todo gadget certificado cabe em 5×5**, folgado dentro de um bloco 6×6. A
premissa do Gate 2 parte 2 está verificada e o catálogo D&C existente é do
tamanho certo. Custo da mudança: trocar a chave de ordenação em
`build_switchers` de `VW` para `(max(bbox), VW)`.


---

## G1.12 — Alcance da contaminação por `max_w`: o que sobrevive

O corte `max_w = 20` estava em **todas** as medições de certificação desta
sessão, não só no Critério C. Ele torna a certificação **mais difícil** (menos
candidatos `W` testados), logo o viés tem direção conhecida — e isso decide,
item a item, o que continua válido.

**Direção do viés.** Certificado com corte ⟹ certificado sem corte. Portanto:

| tipo de afirmação | efeito do corte | estatuto |
|---|---|---|
| "esta forma/posição **certifica**" | não pode ser falso positivo | ✅ **sobrevive** |
| "esta forma/posição **não certifica**" | pode ser falso negativo | ⚠️ **não sobrevive** |

### Sobrevive intacto

- **G1.6, G1.8, G1.8-bis, G1.8-ter — as 18 formas.** São certificados
  positivos. Se certificaram com corte, certificam sem. O rank do span delas
  é o rank de um conjunto de vetores exibidos, e `185/185`, `293/293`,
  `425/425` continuam de pé.
- **G1.3** — os ranks de certificados são cotas inferiores; o déficit 2 é
  cota superior. A conclusão "T0…T4 certificados não geram" **precisaria** ser
  reverificada sem corte, mas a conclusão que se usou dela — que T5 fecha o
  gap — é positiva e sobrevive.
- **Fase 3 (`FASE3_RESULTS.md`)** — os certificados de n ∈ {8,10,12} são pares
  explícitos de tours auditados aresta a aresta. Independentes disto.

### Não sobrevive sem reverificação

- **G1.7 — "42 de 136 são livres de borda".** É **cota inferior**: as 42 são
  genuinamente livres, mas algumas das outras 94 podem ser livres também.
  A tabela de "livre de borda por tamanho de caixa" e a conclusão de que
  *compacidade causa sensibilidade à borda* dependem das não-livres, logo
  ficam **em suspenso**.
- **G2.1 — as 18 posições de falha rígidas.** Conjuntos de falha são
  **supraconjuntos** das falhas verdadeiras. A coincidência exata entre
  n=10, 12 e 14 é sugestiva mas foi medida sob um corte cuja ordem depende de
  `n`; precisa de reverificação sem corte antes de virar enunciado.
- **G1.9 (iv) — "C é necessário, 0 falsos negativos".** Verificado contra as
  42 livres-de-borda com corte. Se o conjunto verdadeiro de livres for maior,
  a necessidade só está testada nas 42.

### Lição de método

O corte não era um parâmetro de desempenho: era um parâmetro que **entrava no
predicado medido**, porque a ordem de enumeração depende de `n`. Qualquer
truncamento cuja ordem dependa da variável do experimento é uma variável
oculta. Antes de reportar dependência (ou independência) em `n`, saturar todo
corte de busca e mostrar a saturação — foi o que expôs isto (`max_w` ∈ {200,
1 000, 10 000} dão o mesmo 94, `20` dá 73).

Este é o mesmo modo de falha da `sec:vies` do paper (viés do amostrador
produzindo rank falso), numa roupa diferente: lá a variável oculta era a ordem
do DFS, aqui é a ordem de enumeração de candidatos.
