# Estado do paper

O paper não fica neste repositório (ver `docs/README.md`). Este documento
registra o que ele reivindica, para que o código seja auditável contra as
afirmações que sustenta.

**Versão corrente:** `knight_tour_tightness`, v5, 15/08/2026.
Substitui a linha `knight_tour_complete` (v1–v3), que é anterior à auditoria
de citações e **não deve ser usada como referência**.

## Sistema de selos

O paper marca cada enunciado com um de três selos, e mantém uma regra
explícita entre eles:

| selo | significado |
|---|---|
| **provado** | tem demonstração válida para toda a família indicada — ou, quando é sobre tabuleiros específicos, estabelecida por construção explicitamente auditada |
| **verificado** | verificado por computação nos casos listados; afirma o que foi medido e nada além |
| **aberto** | conjectura |

> **Regra:** nenhum enunciado marcado apenas como *verificado* é usado como
> premissa de um enunciado *provado*.

Essa disciplina é o que torna o documento auditável, e vale a pena preservá-la
em qualquer reescrita.

## O que está provado

| resultado | escopo |
|---|---|
| `Q(n) = 3` | todo `n >= 6`, mais o caso retangular `Q(n,m) = 3` para `min(n,m) >= 6`. Formalizado em Lean 4 sem `sorry` (`Q_abstract_eq_three`) |
| conexidade do bulk | provada e formalizada sem `sorry`; é o único ingrediente geométrico de `Q(n)=3` |
| lema da redução | `Z_bulk(n) ⊆ Ham(n)` implica `dim Ham(n) = beta_1(n) - 3` |
| **tightness por certificados** | quadrados `n ∈ {6, 8, 10, 12, 14}` e retângulos `6×7`, `6×8`. Por construção, via certificados `C = H_A ⊕ H_B` auditados um a um, **sem amostragem** |

### A tabela de certificados

| `n` | hexágonos | certificados | sem cert. | rank | gera `Z_bulk`? |
|---|---|---|---|---|---|
| 8 | 2.264 | 1.966 | 298 (13,2%) | 101/101 | sim |
| 10 | 5.088 | 4.888 | 200 (3,9%) | 185/185 | sim |
| 12 | 9.000 | 8.848 | 152 (1,7%) | 293/293 | sim |

Rank cheio, com folga: os hexágonos sem certificado não são necessários. O
argumento é unilateral na direção segura — certificados a menos só impediriam
de *provar* a cota inferior, jamais produziriam deficit maior.

**Duas rotas independentes.** Para `n ∈ {8,10,12}` o resultado saiu duas vezes:
por certificação dirigida e por flips hexagonais. A segunda rota acrescenta
`n ∈ {6,14}`. As duas erram de formas diferentes, o que torna a concordância
informativa.

**O caso `6×8` é o mais forte metodologicamente:** 48 casas, enumeração
inviável, e a tightness foi estabelecida **sem conhecer o conjunto de tours** —
bastou certificar hexágonos até o rank atingir `dim Z_bulk = 65`.

## O que está aberto

### Conjectura principal (geração por hexágonos)

> Para todo `n >= 8` par, os hexágonos certificados do bulk geram `Z_bulk(n)`.

Pelos lemas de dimensão e redução, isso implica **tightness para todo `n` par**.

O ganho em relação à formulação anterior ("requer um argumento de deformação
local não coberto") é de precisão: agora é uma afirmação sobre um subespaço
**fixo e explicitamente descrito**, de dimensão `beta_1 - 4`, e não sobre o
rank de uma matriz com `~10^22` linhas.

Duas propriedades favorecem indução `n -> n+2`:
o gadget que produz os certificados tem tamanho **constante**
(`|V(W)| = 8`, em caixa `5×5`), e no interior o grafo do cavalo é invariante
por translação. O passo indutivo tem alvo numérico exato:
`dim Z_1(n+2) - dim Z_1(n) = 12(n-1)`.

### O que governa o limiar

Os cantos dão **só a cota superior**. Os tabuleiros `5×6` e `5×8` têm
exatamente a mesma estrutura de cantos — quatro vértices de grau 2, bulk
conexo, `Q = 3` — e deficit **27** e **9**.

O que separa os casos *tight* é a cota **inferior**, `Z_bulk ⊆ Ham`. Qual
parâmetro a governa permanece em aberto. A obstrução do `5×8` **desaparece ao
alongar o tabuleiro**, o que é incompatível com qualquer explicação por gadget
local.

Esta é, na minha leitura, a lacuna mais relevante: o fenômeno ainda não está
entendido, só delimitado.

## Posicionamento

A literatura de grafos Hamilton-generated opera no regime denso ou
pseudoaleatório, onde a resposta típica é deficit nulo. O grafo do cavalo é o
extremo oposto — esparso, grau mínimo 2 — e ainda assim admite resposta exata.

Ver `docs/prior-art.md` para o detalhamento do que é e do que não é
contribuição.

## Notas de leitura do código

Afirmações do paper que dependem de scripts deste repositório:

| afirmação | código |
|---|---|
| `Q(n)=3`, verificação `n=8,12` | `experiments/06_invariante_Q/deficit_theorem/` |
| localidade de `Q` | `experiments/06_invariante_Q/Q_locality_theorem.py` |
| certificados de tightness | `experiments/06_invariante_Q/tightness_witnesses/` |
| família híbrida | `experiments/06_invariante_Q/tightness_intermediate.py` (exige `K >= 250k`) |
| deficit toroidal nulo | `experiments/06_invariante_Q/rank_ham_torus.py` (exige muitas sementes) |
| viés do amostrador | documentado como resultado, não como erro: é a razão de a tightness ter sido construída sem amostrar |
