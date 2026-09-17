# As 6 direções perdidas do 5×8

Pergunta: por que o `5×8` perde exatamente 6 direções, e a hipótese
"os hexágonos que as gerariam exigem gadgets que não cabem em largura 5".

**Veredicto:**

1. A hipótese na forma **bbox** está **REFUTADA** — hexágonos de altura 3 num
   tabuleiro de 5 linhas já geram as 6 direções (§3).
2. A hipótese de substituição por **contato com a borda** também está
   **REFUTADA como mecanismo** (§6): 2 das 6 direções vêm de hexágonos
   inteiramente interiores, e a linha extra conserta esses tanto quanto os de
   borda. Contato com borda é **correlato**, não causa.
3. **A largura 5 não é a obstrução** (§6b): as 6 formas realizam em `5×10`,
   e duas delas na **posição idêntica** onde eram `no_match` provado em `5×8`.
   A obstrução é **global**, efeito de tamanho finito do `5×8` — não local ao
   hexágono nem devida ao número de linhas.
4. **Fato estabelecido:** acrescentar uma linha (`5×8 → 6×8`) conserta as **6**
   direções, 12/12 colocações certificadas — mas acrescentar **colunas**
   também conserta, então isto é um caso de "aumentar o tabuleiro", não uma
   explicação do limiar 6. **Por que**, segue em aberto.
4. Subproduto independente: **`deficit(6×8) = 3` provado por certificados**
   (§7), sem enumerar os tours de um tabuleiro de 48 casas.

> Registro de superdeclaração: uma versão anterior deste arquivo dizia
> "hipótese revisada CONFIRMADA" e derivava o limiar 6 de "uma linha de folga
> de cada lado". Ambos estavam errados — a aritmética dava 5, e o teste tinha
> sido rodado só nos geradores de borda, que eram justamente onde a explicação
> já funcionava. Corrigido nas §5 e §6.

## Comandos

```bash
cd tightness_witnesses
python missing_directions.py --rows 5 --cols 8     # ~2 min
python missing_directions.py --rows 6 --cols 6     # controle tight
python border_slack.py --rows 5 --cols 8 --exhaustive   # valida o buscador (~17 min)
python border_slack.py --embed 5x8 --into 6x8      # o teste da linha extra
python border_slack.py --span 6x8 --budget 3000000 # prova deficit(6x8)=3 (~4 min)
```

---

## 1. As 6 direções são exatamente a codimensão dos hexágonos realizados

A conta fecha sem folga:

```
rk(Ham) = dim span{τ+τ'} + 1
deficit = β₁ − rk(Ham) = 3 + (dim Z_bulk − dim span{τ+τ'})
```

e no `5×8` o grafo de flips é **conexo** (1 componente, `flip_graph_rect_5x8`),
logo `span{τ+τ'} = span(hexágonos realizados)` — verificado como igualdade de
subespaços, não só de dimensões (`span_diffs_eq_span_realized: true`).

| | 5×8 |
|---|---|
| `β₁` | 51 |
| `dim Z_bulk` (via núcleo dos funcionais de canto) | 47 |
| `rank(todos os 588 hexágonos do bulk)` | 47 = `dim Z_bulk` |
| `rank(hexágonos realizados)` | **41** |
| `dim span{τ+τ'}` | 41 (= igual, como subespaço) |
| `rk(Ham)` | 42 |
| **deficit** | **9** = 3 + **6** |

**As 6 direções = `Z_bulk / span(realizados)`.** Os 588 hexágonos do bulk geram
`Z_bulk` inteiro; só 435 são flip-realizáveis, e esses geram 41.

## 2. O modo de falha é `no_match`, não desconexão

Dos 153 não-realizados:

| modo | n | não-nulos no quociente |
|---|---|---|
| `no_match` — nenhum tour encontra o hexágono num emparelhamento alternado | **151** | **151** |
| `disconnect` — há tour compatível, mas todo XOR quebra em vários ciclos | 2 | **0** |

Partição limpa: os 2 casos de desconexão já estão dentro de
`span(realizados)` e não contribuem nada. **As 6 direções vêm inteiramente de
hexágonos que nenhum tour sequer alcança.** Isso é "o gadget não existe", não
"o gadget existe e desconecta" — o que é favorável ao espírito da hipótese,
ainda que não à sua forma métrica.

## 3. A hipótese bbox — REFUTADA

Filtração: quantas das 6 direções são recuperadas se só se admitem hexágonos
não-realizados de uma dada classe?

| filtro | direções recuperadas (de 6) |
|---|---|
| bbox **altura ≤ 2** | 0 |
| bbox **altura ≤ 3** (de 5 linhas disponíveis) | **6** |
| bbox **largura ≤ 3** | **6** |
| bbox **lado máx ≤ 3** | 0 |
| bbox **lado máx ≤ 4** | **6** |
| todos os não-realizados | 6 |

Hexágonos de **altura 3 num tabuleiro de 5 linhas** já geram as **6** direções.
Cabem com duas linhas de folga. A hipótese "exigem bbox com um lado ≥ 5 ou 6"
é falsa.

Os histogramas confirmam que bbox não separa realizado de não-realizado:

| bbox | realizados | não-realizados |
|---|---|---|
| 4×5 | 88 | 32 |
| 5×5 | 56 | 18 |
| 4×4 | 52 | 20 |
| 3×4 | 27 | 14 |

A forma mais comum (`4×5`) é a mais comum nos dois grupos. E entre os 6
geradores escolhidos gulosamente há um de bbox **3×4**.

## 4. O que discrimina: borda do lado curto

| filtro | direções recuperadas (de 6) |
|---|---|
| só hexágonos que **não tocam linha 0 nem linha 4** (as bordas de 8 casas) | **2** |
| só hexágonos que **não tocam coluna 0 nem coluna 7** (as bordas de 5 casas) | **6** |
| só hexágonos **inteiramente interiores** | 2 |

Evitar as bordas do lado **curto** custa 4 das 6 direções; evitar as do lado
**longo** custa zero. A assimetria é total, e é a única medida testada que
separa.

Os 6 geradores gulosos são todos `no_match`, e todos tocam a linha 0 ou a
linha 4:

| gerador | bbox | bordas tocadas |
|---|---|---|
| 0 | 3×6 | topo |
| 1 | 4×6 | topo |
| 2 | 3×4 | topo, esquerda |
| 3 | 4×4 | topo, esquerda |
| 4 | 5×4 | topo, base, esquerda |
| 5 | 5×6 | topo, base |

## 5. Hipótese revisada — e o resíduo que ela NÃO cobre

> Não é a **forma** que não cabe — é o **switcher** que não tem folga do lado
> de fora quando o hexágono está encostado na borda do lado curto.

⚠️ **Esta hipótese cobre 4 das 6 direções, não as 6.** A própria filtração da
§4 mostra: hexágonos que **não tocam** linha 0 nem linha 4 ainda recuperam
**2** direções, e `interior_only` (não toca borda nenhuma) também dá 2. Existem
hexágonos não-realizados **inteiramente interiores**, com folga acima e abaixo,
e eles respondem por 2 das 6 dimensões perdidas.

Escolhendo os geradores com prioridade para interiores (a escolha por ordem de
enumeração escondia isso — dava 6 geradores todos encostados numa borda):

| gerador | tipo | bbox | linhas ocupadas | bordas |
|---|---|---|---|---|
| 0 | **INTERIOR** | 3×4 | 1–3 | nenhuma |
| 1 | **INTERIOR** | 3×6 | 1–3 | nenhuma |
| 2 | borda | 3×6 | 0–2 | topo |
| 3 | borda | 3×4 | 0–2 | topo |
| 4 | borda | 4×3 | 1–4 | base |
| 5 | borda | 4×5 | 1–4 | base |

Os dois interiores têm altura 3 nas linhas 1–3: **uma linha livre acima e uma
abaixo**. "Sem folga do lado de fora" não se aplica a eles. Por que falham
segue **em aberto**.

## 6. O teste da linha extra — e o que ele realmente mostra

Cada gerador foi embutido em `6×8` — **mesma forma, mesma posição relativa**,
uma linha a mais — nos dois deslocamentos verticais (`dr = 0, 1`), com busca
dirigida por certificado positivo.

**12 de 12 colocações realizam — incluindo as dos dois geradores INTERIORES.**

| gerador | tipo | `dr=0` | nós | `dr=1` | nós |
|---|---|---|---|---|---|
| 0 | **interior** | realizado | 24.854 | realizado | 92.731 |
| 1 | **interior** | realizado | 5.950 | realizado | 112.649 |
| 2 | borda | realizado | 421 | realizado | 31.543 |
| 3 | borda | realizado | 316.757 | realizado | 2.000.586 |
| 4 | borda | realizado | 140.324 | realizado | 232.806 |
| 5 | borda | realizado | 22.666 | realizado | 99.397 |

Todos são `no_match` **provado** (busca saturada) em `5×8`.

> A marca "não saturou" aparece numa linha, mas é irrelevante: ela só limitaria
> conclusões **negativas**. Estes são certificados **positivos** — um tour
> concreto foi exibido e `τ XOR C` verificado hamiltoniano.

### O que isto REFUTA — inclusive a hipótese da §5

Que os dois geradores **interiores** também sejam consertados pela linha extra
mata a explicação por borda como *mecanismo*. Eles já tinham uma linha livre
acima e uma abaixo em `5×8`; ganharam uma terceira e passaram a realizar.
Contato com a borda é **correlato**, não causa.

Também não sobrevive a formulação "6 deixa uma linha de folga de cada lado
para um hexágono de altura 3": essa conta dá `3+1+1 = 5`, e 5 não basta. O
gerador 2 ocupa as linhas 0–2 e tem **zero** linhas acima tanto em `5×8`
quanto em `6×8`; o que mudou foi 2 → 3 linhas abaixo. Não houve simetria de
folga em momento algum.

### O que fica

Estabelecido: **acrescentar uma linha conserta todas as 6 direções**, de borda
e interiores igualmente. Por que — se é o número de linhas, a área total, ou
outra coisa — **não está estabelecido**.

O controle que separa "mais linhas" de "mais espaço" é embutir os mesmos
geradores em `5×10`: mesmas 5 linhas, 20 casas a mais. Resultados na §6b.

## 6b. O controle `5×10` — e por que ele é caro

Um acidente numérico torna `5×10` o controle ideal para `6×8`:

| | 5×10 | 6×8 |
|---|---|---|
| casas | 50 | 48 |
| `E` | 118 | 116 |
| `β₁` | **69** | **69** |
| `dim Z_bulk` | **65** | **65** |
| `Q` | 3 | 3 |

`β₁` e `dim Z_bulk` **idênticos**, e `5×10` é ligeiramente *maior* em casas e
arestas. Se `deficit(6×8) = 3` (§7, provado) e `deficit(5×10) > 3`, a diferença
não pode ser tamanho, área, número de arestas nem dimensão do espaço de
ciclos — só a proporção.

### Resultado: largura 5 NÃO é a obstrução

`border_slack.py --embed 5x8 --into 5x10 --budget 20000000` — os 6 geradores,
em cada deslocamento de coluna `dc ∈ {0,1,2}`, num tabuleiro que **continua com
5 linhas**:

| gerador | tipo | `dc=0` | `dc=1` | `dc=2` |
|---|---|---|---|---|
| 0 | interior | — | **real.** (10.451) | **real.** (6.987) |
| 1 | interior | — | **real.** (196) | — |
| 2 | borda | **real. (3.061)** | **real.** (3.314) | — |
| 3 | borda | — | **real.** (3.500) | **real.** (402) |
| 4 | borda | — | **real.** (92.851) | **real.** (639) |
| 5 | borda | **real. (8.719)** | **real.** (16.598) | — |

(`—` = `not_found_in_budget` a 40 M nós, **não** prova de inexistência.)

**Todos os 6 realizam em pelo menos uma posição, com 5 linhas.** E os
geradores 2 e 5 realizam em `dc=0`, onde as células são **idênticas** às do
`5×8` — mesmas linhas, mesmas colunas, mesma distância às bordas superior e
inferior. Lá são `no_match` **provado** (busca saturada); aqui realizam.

Consequências:

1. **A largura 5 não proíbe nenhuma dessas formas.** A leitura "5 linhas não
   dão folga" está morta em todas as suas versões.
2. **A obstrução não é local.** A mesma configuração local, com a mesma
   vizinhança de linhas, muda de estado só porque o tabuleiro ficou mais
   comprido. É efeito global de tamanho finito do `5×8`.
3. **"Acrescentar uma linha conserta" perde o papel de explicação.** É um caso
   particular de "aumentar o tabuleiro conserta" — acrescentar **colunas**
   também conserta.

O que continua **não** estabelecido: que `deficit(5×10) = 3`. Realizar estas 6
formas não implica que os hexágonos realizados do `5×10` gerem `Z_bulk(5×10)` —
são 916 hexágonos e `dim Z_bulk = 65`. Ver abaixo.

### A pergunta em aberto: estado das medições de `deficit(5×10)`

Três tentativas, nenhuma decisiva:

| via | resultado |
|---|---|
| `deficit_direct.py --rows 5 --cols 10` (enumeração exaustiva) | **inviável** — extrapolando `6×6→6×7` (≈2,2× por casa), `5×10` deve ter ≈10⁸ tours; abortado aos 39 min ainda enumerando |
| `border_slack.py --span 5x10 --budget 2000000` | cortado em 200/916 hexágonos (2227 s): `rank = 44/65`, **48 `not_found_in_budget`** |
| `border_slack.py --embed 5x8 --into 5x10` | em curso |

⚠️ **`not_found_in_budget` não é `no_match`.** O contraste com `6×8` é grande —
lá o mesmo procedimento atingiu `rank = 65/65` no hexágono 427 em 241 s, com
**zero** falhas de orçamento — mas contraste de custo de busca não é prova de
não-existência. Registrado como sugestivo, não como resultado.

### Por que a resposta não é óbvia

Os dados de largura 5 que existem mostram a codimensão **caindo** com colunas:

| board | `dim Z_bulk` | `dim span{τ+τ'}` | codim | deficit |
|---|---|---|---|---|
| 5×6 | 29 | 5 | 24 | 27 |
| 5×8 | 47 | 41 | **6** | 9 |

`24 → 6`. Se a tendência continuar, `5×N` poderia atingir codim 0 para `N`
grande — o que faria o limiar **não** ser em `min(n,m)`. Nada nos dados atuais
exclui isso. (Ressalva: `5×6` tem só 8 tours e grafo de flips totalmente
desconexo, então os dois pontos não estão no mesmo regime.)

Esta é hoje a pergunta em aberto mais direta do assunto:
**`deficit(5×N)` estabiliza em algum valor > 3, ou converge para 3?**

## 7. `deficit(6×8) = 3` — PROVADO sem enumerar tours

O mesmo maquinário dá um resultado, não só uma explicação. O argumento é
unilateral na direção segura:

```
span{τ+τ'} ⊇ span(hexágonos realizados)      [todo flip é diferença de tours]
span{τ+τ'} ⊆ Z_bulk                          [τ, τ' têm y=1 ⟹ τ+τ' tem y=0]
```

Logo, se os hexágonos **certificados** já geram `Z_bulk`, as duas inclusões
fecham em igualdade e o deficit está determinado — sem conhecer o conjunto de
tours.

`python border_slack.py --span 6x8 --budget 3000000` (241 s):

| | 6×8 |
|---|---|
| `V`, `E`, `β₁` | 48, 116, 69 |
| hexágonos do bulk | 1.126 |
| `dim Z_bulk` | 65 |
| **`rank(certificados)`** | **65** = `dim Z_bulk` |
| ⟹ `rk(Ham)` | 66 |
| ⟹ **`deficit`** | **3** = `Q(6,8)` |

Posto pleno atingido no hexágono 427 de 1.126; 65 certificados individuais,
cada um com um tour concreto e `τ XOR C` verificado hamiltoniano. Houve 1
`not_found_in_budget`, que **não entra na conta** — o posto pleno foi atingido
sem ele. (Se entrasse, só poderia aumentar o rank, nunca diminuí-lo.)

`6×8` tem 48 casas: enumeração exaustiva de tours é inviável. Este é o primeiro
ponto de tightness do projeto obtido **por certificados**, não por enumeração.

### O quadro do limiar, agora

| board | `min(n,m)` | `Q` | bulk | deficit | tight? | método |
|---|---|---|---|---|---|---|
| 5×6 | 5 | 3 | conexo | **27** | não | exaustivo |
| 5×8 | 5 | 3 | conexo | **9** | não | exaustivo |
| 6×6 | 6 | 3 | conexo | 3 | **sim** | exaustivo |
| 6×7 | 6 | 3 | conexo | 3 | **sim** | exaustivo |
| 6×8 | 6 | 3 | conexo | 3 | **sim** | **certificados** |

`Q` e conexidade do bulk são constantes nas cinco linhas. Só `min(n,m)`
acompanha o deficit.

---

## Ressalvas metodológicas

- A escolha dos 6 geradores é **gulosa** e portanto arbitrária. Todas as
  conclusões acima vêm da **filtração**, que é invariante à ordem — não da
  bbox dos 6 geradores específicos.
- `border_slack.py` produz **certificados positivos**. Quando não encontra
  dentro do orçamento de nós, reporta `not_found_in_budget`, que **não é**
  prova de inexistência; só o modo `--exhaustive` (busca saturada, saturação
  reportada) distingue `no_match` de `disconnect` como fato provado.
- O buscador **foi validado** contra a classificação exaustiva independente do
  `missing_directions.py` antes de ser usado no `6×8`: DFS restrito (força `m₀`,
  proíbe `m₁`) reproduz **435 / 151 / 2** exatamente, com `n_unsaturated = 0`
  — busca saturada em todos os 588 hexágonos, então a distinção
  `no_match` vs `disconnect` é fato provado, não limite de orçamento.
  São dois métodos independentes: varrer os 44.202 tours vs. busca dirigida.

## Invariante barata que qualquer script desta família deve checar

```
rank(qualquer conjunto de ciclos do bulk)  ≤  dim Z_bulk
```

Todo hexágono do bulk **é** um elemento de `Z_bulk`, logo o rank de um conjunto
deles nunca pode exceder `dim Z_bulk`. O bug do `zbulk_basis` (filtrar a *base*
de `Z₁` por `y=0` em vez de tomar o núcleo — a fibra é subespaço, não
subconjunto da base) produzia `dim Z_bulk = 31` contra `rank(realizados) = 41`:
os hexágonos gerariam **mais** que o espaço que os contém. Absurdo detectável
**sem nenhum valor de referência externo** — e eu não o vi; quem pegou o erro
foi o controle com valor conhecido. Agora está como `assert` em
`missing_directions.py`.
