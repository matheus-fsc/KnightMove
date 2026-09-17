# Fase 3 — o parity-switcher do CNP é construtível em `G_n`?

> Executa `FASE3_PROMPT.md`. Data: 2026-08-06.
> Scripts: `bulk_perp_dims.py`, `switcher_search.py`, `hexagon_span.py`,
> `good_hexagon_span.py`, `hexagon_in_tourspace.py`.
> Dados: `data/bulk_perp_dims.json`, `data/switcher_n{6,8,10,12}.json`,
> `data/hexagon_span.json`, `data/good_hexagon_span.json`,
> `data/hexagon_in_tourspace.json`, `data/full_audit.log`.

## Comandos exatos

```bash
python bulk_perp_dims.py                                          # Tarefa 1, <1 s
python hexagon_span.py --ns 6,8,10,12,14 --lengths 6              # extra,  ~3 min
python switcher_search.py --n 6  --limit 400 --max-internal 4 --max-w 20
python switcher_search.py --n 8  --limit 150 --max-internal 4 --max-w 20   # ~60 s
python switcher_search.py --n 10 --limit 100 --max-internal 4 --max-w 20   # ~62 s
python switcher_search.py --n 12 --limit 100 --max-internal 4 --max-w 20   # ~7 s
python good_hexagon_span.py --ns 8,10,12 --max-hex 4000            # §3.5, ~6 min
python hexagon_in_tourspace.py --ns 8,10,12                        # §3.6, enumeração
                                                                   # COMPLETA + auditoria
```

Tempos do §3.6 (single-thread): n=8 475 s, n=12 2.329 s, n=10 37.803 s (esta
máquina ficou ociosa durante parte da execução; o número não é comparável).

---

## Veredicto

**O esquema do CNP é instanciável em `G_n`.** Exibimos parity-switchers
completos — ciclo `C` + caminhos `P₂,P₃` vértice-disjuntos + caminho
hamiltoniano de (S3) — verificados independentemente, para
`n ∈ {6,8,10,12}`. A taxa de sucesso **cresce** com `n` e atinge 100% em
`n=12`.

Preenchendo a tabela de interpretação do prompt:

| Resultado de 2b | |
|---|---|
| **Existe `W` com (S3) satisfeito** | ✅ **é este o caso** |
| (S3) falha para todo `W`, com causa identificada | não |
| (S3) falha sem causa clara | não |

Isto **reverte** a conclusão da Fase 2 (§2.2/§2.3 do `RESULTS.md`), que estava
errada — ver §5 abaixo.

Mais: a auditoria end-to-end (§3.6) transformou isto num resultado mais forte
que o pedido. Com **enumeração completa** dos hexágonos do bulk e **zero
falhas de auditoria**, obtém-se `Z_bulk ⊆ Span(Ham)` por construção explícita
e, pelo teorema de redução do §3.45,

> **`rank(Ham(n)) = β₁ − 3` (tightness) para `n ∈ {8,10,12}`, provada por
> construção** — cada dimensão atingida por um par explícito de passeios do
> cavalo cuja diferença simétrica é um hexágono, verificados aresta a aresta.

⚠️ Isto **não** prova a tightness para todo `n`. O que falta está no §4.

---

## 1. Tarefa 1 — Lemas A e B: dimensões

`dim C^⊥ = |V|−1`, `dim X = |V|+2`, `dim perp(Z_bulk) = |V|+3` em **todos** os
12 tabuleiros testados. `X ⊆ perp(Z_bulk)` e `perp(Z_bulk) = X ⊕ ⟨1_e⟩`.

| tabuleiro | V | E | comp. bulk | dim C^⊥ | dim X | dim perp | Q | codim | Lema A |
|---|---|---|---|---|---|---|---|---|---|
| 6×6 | 36 | 80 | 1 | 35 | 38 | 39 | 3 | 1 | OK |
| 8×8 | 64 | 168 | 1 | 63 | 66 | 67 | 3 | 1 | OK |
| 10×10 | 100 | 288 | 1 | 99 | 102 | 103 | 3 | 1 | OK |
| 12×12 | 144 | 440 | 1 | 143 | 146 | 147 | 3 | 1 | OK |
| 14×14 | 196 | 624 | 1 | 195 | 198 | 199 | 3 | 1 | OK |
| 16×16 | 256 | 840 | 1 | 255 | 258 | 259 | 3 | 1 | OK |
| 6×8, 6×10, 8×10 | | | 1 | V−1 | V+2 | V+3 | 3 | 1 | OK |
| 6×7, 7×9, 8×12 | | | 1 | V−1 | V+2 | V+3 | 3 | 1 | OK |

Predição falseável `codim = 1`: **confirmada em 12/12**, incluindo
retangulares e lados ímpares.

### 1.1 Isto é demonstração, não só verificação

O Lema A segue dos dois teoremas que você já provou:

- `perp(Z_bulk) = ⟨1_e : e ∈ Mand⟩ ⊕ row(∂₁^{bulk})`, pois arestas de `Mand`
  incidem em canto e ciclos do bulk não as tocam. Logo
  `dim perp(Z_bulk) = 8 + (|V|−4) − c`, onde `c` = nº de componentes do bulk.
- **Conexidade do bulk** (provada, `n ≥ 6`) ⟹ `c = 1` ⟹ `dim = |V|+3`.
- **`Q(n)=3` / `Q(n,m)=3`** (provados, formalizados em Lean) ⟹
  `dim X = (|V|−1) + 3 = |V|+2`.

Portanto `codim = 1` **para todo `n ≥ 6`** — não só para os 12 testados. A
tabela é confirmação numérica de um fato já demonstrável.

---

## 2. Tarefa extra — os hexágonos geram `Z_bulk`

Esta não estava no prompt e fecha uma lacuna que o prompt não menciona.

**O problema.** O Lema B dá "algum ciclo do bulk com interseção ímpar", mas
não diz o **comprimento**. E o comprimento importa: `v₁` e `v_{k+1}` estão a
distância `k` em `C`, então têm cores opostas ⟺ `k` ímpar ⟺ `|C| ≡ 2 (mod 4)`.
Se o único ciclo disponível tiver `|C| ≡ 0 (mod 4)`, (S3) é **impossível por
paridade de cor** e o esquema morre.

**O resultado.** Os hexágonos do bulk sozinhos **geram `Z_bulk`**:

| n | dim `Z_bulk` | nº hexágonos no bulk | rank do span | gera? |
|---|---|---|---|---|
| 6 | 41 | 532 | 41 | ✅ |
| 8 | 101 | 2.264 | 101 | ✅ |
| 10 | 185 | 5.088 | 185 | ✅ |
| 12 | 293 | 9.000 | 293 | ✅ |
| 14 | 425 | 14.000 | 425 | ✅ |

Consequência: se `R ∉ X`, então `R ∉ perp(Z_bulk)` (Lema B), e como os
hexágonos geram `Z_bulk`, existe um **hexágono do bulk** `C` com
`⟨C,R⟩ = 1`. O caso `k=3` não é uma escolha de conveniência — é suficiente.

*Status:* verificado computacionalmente para `n ∈ {6,8,10,12,14}`; **não
provado para todo `n`**. A contagem de hexágonos (532, 2.264, 5.088, 9.000,
14.000) cresce como `Θ(n²)`, com folga enorme sobre `dim Z_bulk` — mas folga
de contagem não é prova de geração, e nenhuma construção local foi tentada.

---

## 3. Tarefa 2 — construção do switcher

Def. 2.2 lida no PDF original. Para `k=3`: `C = (v₁…v₆)`, `P₂` liga `v₂↔v₆`,
`P₃` liga `v₃↔v₅`; `v₁,v₄` são os polos (grau 2 em `W`), os demais têm grau 3.

> **Nota sobre o texto do CNP:** a Def. 2.2 diz caminhos **vértice-disjuntos**,
> mas o passo (S2.b) diz *"edge-disjoint (short) paths"*. É uma inconsistência
> do próprio artigo. Usamos vértice-disjunto (a definição), que é a condição
> mais forte — nossos resultados valem *a fortiori* sob a leitura fraca.

### 3.1 (S2.b) e (S3) — resultados

Busca: `max_internal = 4`, até 20 candidatos `W` por par de polos, 3 pares de
polos por hexágono, orçamento de 400k nós por busca hamiltoniana. Amostra
aleatória de hexágonos do bulk (seed 0).

| n | hexágonos amostrados | admitem `W` (S2.b) | **switcher COMPLETO (S3)** | falhas |
|---|---|---|---|---|
| 6 | 400 → 299 c/ W | 299 | **14 (4,7%)** | 285 |
| 8 | 150 | 150 (100%) | **131 (87,3%)** | 19 |
| 10 | 100 | 100 (100%) | **97 (97,0%)** | 3 |
| 12 | 100 | 100 (100%) | **100 (100%)** | 0 |

O `|V(W)|` mínimo é **8** em todos os `n`: hexágono + 1 vértice interno em
`P₂` + 1 em `P₃`. É o switcher mais econômico possível.

Exemplo verificado (n=12):
```
C  = H7 J6 H5 J4 H3 G5      P2 = [I4]   P3 = [F4]
|V(W)| = 8;  (S3): caminho hamiltoniano cobrindo os 138 vértices restantes
```

### 3.2 Causa das falhas — identificada

O filtro de cor **nunca** rejeitou nada (`cor=0` em todos os `n`). Isso não é
sorte: `v₂,v₆` têm a mesma cor e `v₃,v₅` também, então os internos de `P₂`
desbalanceiam +1 numa cor e os de `P₃` +1 na outra — o buraco é
automaticamente balanceado. **A restrição de cor é satisfeita por construção
quando `k` é ímpar.**

A causa dominante de falha é outra: **um canto ficar com < 2 vizinhos
restantes**. Canto tem grau 2; se `W` engole um dos seus dois vizinhos, o
canto não pode ser interno a um caminho hamiltoniano e não é extremo. Isso
rejeitou 269–2225 candidatos `W` conforme o `n`.

Esta causa é **um efeito de tabuleiro pequeno** e desaparece com `n`: no 6×6
o bulk tem só 32 vértices e quase todo hexágono fica perto de um canto; no
12×12, nenhum hexágono amostrado ficou sem solução. A monotonia
4,7% → 87,3% → 97,0% → 100% é consistente com "obstrução de fronteira que
some no limite", não com obstrução estrutural.

### 3.3 O que é "provado ausente" vs "não encontrado"

Distinguimos os dois. Somando todos os `n`: (S3) **provado ausente** (busca
exaurida) em 3–120 candidatos `W` por tabuleiro; **estourou orçamento** em
2–22. Os "sem orçamento" são não-resultados honestos, não impossibilidades.
Nenhum hexágono foi classificado como falha só por estouro de orçamento nos
casos `n ≥ 10`.

### 3.4 Verificação independente

Todo `W` aceito passa por `verify_switcher()`, que confere do zero: `C` é
6-ciclo de `G_n`; `P₂,P₃` são caminhos legítimos, vértice-disjuntos e
internamente disjuntos de `C`; os **dois** caminhos hamiltonianos de `W`
existem, vão de `v₁` a `v₄` e cobrem `V(W)`; cada aresta de `C` está em
exatamente um deles (é o que garante paridades distintas); e o caminho de
(S3) é hamiltoniano no grafo com buraco, sem invadir o interior de `W`.
**Zero falhas de verificação** em todas as execuções.

---

## 3.45 A redução: **tightness ⟺ `Z_bulk ⊆ Span(Ham)`**

Isto saiu da auditoria do §3.6 e reenquadra a Fase 3 inteira.

**Observação dimensional.** `dim Z_bulk = (E−8) − (V−4) + 1 = β₁ − 4` para todo
`n ≥ 6` (usa conexidade do bulk). E `rank(Ham) = β₁ − 3` é a tightness. Ou
seja, `Z_bulk` tem **exatamente uma dimensão a menos** que o espaço dos tours.

**Teorema (redução).** Para `n` par `≥ 6`:
`rank(Ham(n)) = β₁ − 3` ⟺ `Z_bulk ⊆ Span(Ham(n))`.

*Prova.* (⇐) Todo tour `τ` contém as 8 obrigatórias, enquanto todo vetor de
`Z_bulk` é nulo nelas; logo `τ ∉ Z_bulk`. Se `Z_bulk ⊆ Span(Ham)`, então
`dim Span(Ham) ≥ dim Z_bulk + 1 = β₁ − 3`. Como `dim Span(Ham) ≤ β₁ − 3` já
está provado (dual de `Q(n)=3`), vale a igualdade.
(⇒) Se vale a tightness, `Span(Ham) = perp(X)`. E `Z_bulk ⊆ perp(X)`, pois
ciclos são ortogonais a cortes e os ciclos do bulk não tocam `Mand`. ∎

**Verificado em n=6** com os 9.862 tours exatos:

```
rank(Ham) = 42 = β₁−3        dim Z_bulk = 41 = β₁−4
rank(Ham ∪ Z_bulk) = 42      →  Z_bulk ⊆ Span(Ham)  ✓
rank(Z_bulk ∪ {τ}) = 42      →  Span(Ham) = Z_bulk ⊕ ⟨τ⟩  ✓
```

**Por que isto importa.** A conjectura de tightness — uma afirmação sobre o
rank de uma matriz com bilhões de linhas — vira uma afirmação sobre um
subespaço **fixo e explícito**, `Z_bulk`, de dimensão `β₁−4`. E os
certificados do §3.6 exibem elementos de `Z_bulk ∩ Span(Ham)` **um a um, por
construção**.

---

## 3.5 Os hexágonos **bons** também geram `Z_bulk`

Esta é a computação que fecha o buraco mais concreto (era o item 3 da lista
do §4; ficou resolvido).

**O problema.** Numa prova por contradição, `R` não é escolhido por nós: ele
"acende" algum hexágono, e não controlamos qual. Que 87–100% dos hexágonos
admitam switcher não basta se `R` acender só os ruins.

**O resultado.** Restringindo o span apenas aos hexágonos que **passam em
(S3)**:

| n | dim `Z_bulk` | hexágonos bons | ruins | rank do span dos bons | gera? |
|---|---|---|---|---|---|
| 8 | 101 | 137 | 27 | 101 | ✅ |
| 10 | 185 | 501 | 20 | 185 | ✅ |
| 12 | 293 | 1.141 | 20 | 293 | ✅ |

Em cada caso a varredura **parou assim que o rank saturou** — os hexágonos
bons já geram `Z_bulk` sem precisar de todos.

Combinando com o Lema B: se `R ∈ C_n^⊥ ∖ X`, então `R ∉ perp(Z_bulk)`; como os
hexágonos bons geram `Z_bulk`, existe um hexágono `C` que é **simultaneamente**
`⟨C,R⟩ = 1` e portador de um switcher completo. Os hexágonos ruins deixam de
importar.

*Status:* verificado para `n ∈ {8,10,12}`; não provado para todo `n`.
Numa primeira tentativa com 900 hexágonos, `n=12` ficou em 292/293 — era
limite de amostra, não falha; com 4.000 saturou.

---

## 3.6 Auditoria end-to-end + enumeração COMPLETA

O §3.5 rodou sobre amostra. Como a afirmação "os certificados geram `Z_bulk`"
é **monótona** (rank só cresce ao acrescentar vetores, e cada vetor entra com
certificado positivo), amostrar só pode *sub*-reportar — mas o ponto de falha
real seria um **falso positivo**. Por isso: enumeração completa **e**
verificação independente.

### O certificado

Para cada hexágono, em vez de confiar no `verify_switcher`, reconstruímos do
zero:

```
H'   = caminho hamiltoniano de v1 a v4 em G ∖ (V(W)∖{v1,v4})     [passo (S3)]
H_A  = H' ∘ (1º caminho hamiltoniano de W)
H_B  = H' ∘ (2º caminho hamiltoniano de W)
```

e verificamos que `H_A` e `H_B` são **ciclos hamiltonianos genuínos de `G_n`**
(todo vértice com grau 2, todas as arestas existem em `G_n`, `n²` arestas,
**componente única** — o teste de conexidade é o que exclui união de
sub-ciclos) e que

$$H_A \oplus H_B = C.$$

Isto é estritamente mais forte que "o esquema do CNP é instanciável": **exibe
dois passeios do cavalo cuja diferença simétrica é exatamente o hexágono**,
provando construtivamente `C ∈ Span(Ham(n))`.

### Resultados (enumeração completa de todos os hexágonos do bulk)

| n | hexágonos no bulk | certificados | sem certificado | **falhas de auditoria** | rank do span | gera `Z_bulk`? |
|---|---|---|---|---|---|---|
| 8 | 2.264 (**todos**) | 1.966 | 298 | **0** | 101/101 | ✅ |
| 10 | 5.088 (**todos**) | 4.888 | 200 | **0** | 185/185 | ✅ |
| 12 | 9.000 (**todos**) | 8.848 | 152 | **0** | 293/293 | ✅ |

Enumeração completa nos três: nenhum hexágono do bulk ficou de fora. A fração
sem certificado cai monotonicamente (13,2% → 3,9% → 1,7%), e os certificados
sozinhos já geram `Z_bulk` em todos os casos — os não-certificados são
irrelevantes para a conclusão.

Nenhum certificado falhou a auditoria em nenhuma execução.

### O que isto estabelece

Pelo teorema de redução do §3.45, `Z_bulk ⊆ Span(Ham)` **é** a tightness.
Como os certificados geram `Z_bulk` em `n ∈ {8,10,12}`, obtém-se

> **`rank(Ham(n)) = β₁ − 3` para `n ∈ {8,10,12}`, por construção** — sem
> calcular o rank de uma matriz de tours amostrados.

Isto é metodologicamente melhor que a evidência anterior: o `deficit(8)=3` do
`deficit_theorem/` vinha do rank de 3.000 tours amostrados via Z3 (e a Fase 2
mostrou como amostragem enviesada pode falsear justamente esse número). Aqui
não há amostragem de tours: cada dimensão de `Z_bulk` é atingida por um par
explícito de tours verificados aresta a aresta.

---

## 4. O que falta para virar prova

1. **Uniformizar em `n`.** Tudo aqui é amostragem em `n ∈ {6,8,10,12}`. Uma
   prova precisa de uma construção explícita: dado um hexágono qualquer do
   bulk, exibir `P₂,P₃` e o caminho de (S3) por construção. A conexidade do
   bulk e o D&C por blocos são os candidatos naturais.
2. **Provar que hexágonos geram `Z_bulk`** para todo `n ≥ 6` (§2). Hoje é
   empírico.
3. ~~**Tratar os hexágonos ruins.**~~ ✅ **Resolvido no §3.5:** os hexágonos
   que passam em (S3) já geram `Z_bulk` sozinhos (`n ∈ {8,10,12}`), então todo
   `R ∉ X` acende ao menos um hexágono bom. Falta provar isso para todo `n`.
4. **(S4)/(S5)** não foram tocados — são os passos triviais do recipe, mas não
   verificados aqui.

### 4.1 Estado do esqueleto (S1)–(S5)

| passo | estado |
|---|---|
| (S1) tomar `R` do Lema 2.1 | ✅ adaptado ao caso bipartido (Fase 2, §2.1) |
| (S2.a) ciclo par com `⟨C,R⟩` ímpar | ✅ **de graça** pelos Lemas A+B, e pode-se tomar hexágono (§2, §3.5) |
| (S2.b) caminhos `P₂,P₃` | ✅ existem; `|V(W)|=8` mínimo (§3.1) |
| (S3) caminho hamiltoniano com buraco | ✅ verificado, 100% no 12×12 (§3.1) |
| (S4)/(S5) | não verificados (passos formais do recipe) |
| **uniformidade em `n`** | ❌ **é o que falta** — tudo acima é `n ∈ {6..12}` |

Ou seja: o esqueleto inteiro está disponível ponto a ponto nos tabuleiros
testados. O trabalho restante é trocar amostragem por construção.

---

## 5. Correção da Fase 2

O `RESULTS.md` §2.2 concluía que nenhum parity-switcher pode viver no bulk. O
argumento era: *"toda classe não-trivial tem representante suportado em
`Mand`, logo ciclo que evita cantos tem interseção par com todo witness"*.

**Isso está errado e o erro é circular.** "Toda classe tem representante em
`Mand`" descreve as classes de `X`. Mas o `R` de uma prova por contradição é
exatamente um `R ∈ C_n^⊥ ∖ X` — hipotético, fora de `X`, e sobre cujas classes
aquela afirmação nada diz. Assumir que não existe é assumir a tightness, que é
a conclusão.

O Lema B dá o oposto: `C_n^⊥ ∩ perp(Z_bulk) = X`, logo todo `R ∈ C_n^⊥ ∖ X`
tem interseção **ímpar** com algum ciclo do bulk. (S2.a) sai de graça
justamente para os `R` que importam.

O §2.3 (bloqueio de grau) cai junto, e por um segundo motivo independente: ele
supunha `v₁` isolado em `G ∖ (V(W)∖{v₁,v₄})`. Falso — `v₁` perde só seus dois
vizinhos de `C` e mantém os demais (grau até 8 em `G_n`). O argumento só
valeria se `v₁` fosse canto, o que decorria do §2.2 já refutado.

O que **sobrevive** da Fase 2, intacto: o mapa de witnesses de baixo peso, o
lema da complementação (§2.1), as medidas de (C3), e o registro da armadilha
de amostragem.
