# Fase 2 — mapa computacional dos witnesses duais

> Executa a Fase 2 de `paper/notes/tightness_plan.md`.
> Scripts: `knight_gf2.py` (infra), `witness_map.py` (análise).
> Dados: `data/witness_map_n{6,8}.json`.
> Data: 2026-08-06.

---

## 0. Correção de enquadramento da Fase 2

O plano descreve a Fase 2 como o **teste decisivo da tightness**: *"se aparecer
um R inesperado, a tightness pode ser FALSA para n grande"*. Isso não é
exatamente o que a Fase 2 pode entregar, e vale corrigir antes de ler os
números.

`deficit(n) = dim(C_n^⊥ / C^⊥)` já é conhecido e igual a 3 para
`n ∈ {6,8,10,12}` (`deficit_theorem/`). O teorema `Q(n)=3` já exibe 3 classes
independentes. Logo, **para esses n a contagem de dimensão já força
`C_n^⊥ = X`** — não há espaço para um witness "inesperado". Nenhuma
enumeração pode refutar a tightness em `n` onde o rank já foi calculado; ela
só poderia ser refutada em `n ≥ 14`, onde o rank é desconhecido.

O que a Fase 2 entrega de fato — e que é o que a Fase 3 precisa — é o **mapa
estrutural**: quais suportes realizam cada classe, onde eles vivem no
tabuleiro, e como a condição (C3) do CNP se comporta. Foi isso que rodou.

---

## 1. Resultados n=6 (exaustivo, 9.862 tours)

| Quantidade | Valor |
|---|---|
| V, E, β₁ | 36, 80, 45 |
| rank(Ham) | 42 |
| deficit | 3 |
| dim C^⊥ = row(∂₁) | 35 (= V−1) |
| dim C_n^⊥ | 38 (= E − rank) |
| dim(C_n^⊥/C^⊥) | 3 ✓ |
| vetores de C_n^⊥ com peso ≤ 6 | 814 |
| witnesses fora das 8 classes previstas | **0** |

As 8 classes previstas por `Q(n)=3` (a saber `0`, os 6 pares `r_i+r_j`, e
`r₁+r₂+r₃+r₄`) são distintas mod `C^⊥` e esgotam o quociente.

### 1.1 Pesos mínimos por classe

| Classe | \|R\|min | #mínimos | #que evitam cantos | exemplo |
|---|---|---|---|---|
| `r_i+r_j` (6 classes) | **2** | 4 cada | **0** | `A6-C5 ⊕ F6-D5` |
| `r₁+r₂+r₃+r₄` | **4** | 16 | **0** | `A6-C5 ⊕ F6-D5 ⊕ B3-A1 ⊕ E3-F1` |

Confirma que `φ₁ = x_{F6-D5} ⊕ x_{B3-A1}` é um witness de peso mínimo (classe
`r₂+r₃`), e que **não existe detector mais curto** em nenhuma classe.

### 1.2 Localização — o resultado central

Contagem de witnesses de peso ≤ 6 que **evitam toda aresta incidente a canto**:

| Classe | total |
|---|---|
| `0` (cortes, triviais) | 64 |
| todas as 7 classes não-triviais | **0** |

**Nenhum witness não-trivial de peso ≤ 6 vive no bulk.**

---

## 1-bis. Resultados n=8 (4.000 tours amostrados, peso ≤ 4)

`V=64, E=168, β₁=105`, **rank(Ham) = 102, deficit = 3** — bate com
`deficit_theorem/` (que gastou ~12 min de Z3; aqui são 8 s).

Replica o n=6 ponto a ponto:

| | n=6 | n=8 |
|---|---|---|
| witnesses fora das 8 classes | 0 | 0 |
| \|R\|min classes-par | 2 (×4 cada) | 2 (×4 cada) |
| \|R\|min classe 4-cantos | 4 (×16) | 4 (×16) |
| witnesses não-triviais no bulk | **0** | **0** |
| (C3) `min_v e_R(v)/deg(v)` | 0.500 | 0.500 |
| \|R\|max | 78 = 80−2 | 166 = 168−2 |

As contagens de mínimos (4 e 16) são idênticas nos dois tabuleiros — como
esperado, já que vêm de escolher 1 das 2 arestas obrigatórias por canto
(`2²=4` e `2⁴=16`), independente de `n`.

### ⚠ Armadilha de amostragem (vale registrar)

A primeira rodada de n=8 deu `rank = 96, deficit = 9` e **786 witnesses
"desconhecidos"**, incluindo alguns de peso 1 — o que "refutaria" a
tightness. Era artefato: com Warnsdorff, o DFS fecha um tour na primeira
tentativa e as opções seguintes do nível raso nunca são exploradas, deixando
as duas arestas de saída de `C7` com frequência **zero**. Toda aresta nunca
usada vira um witness espúrio de peso 1.

Correção em `knight_gf2.py`: ordem uniforme nos 4 primeiros níveis + jitter
no Warnsdorff abaixo disso. `rank` sobe 96 → 100 → **102**.

**Consequência metodológica:** qualquer futuro teste de tightness em `n ≥ 14`
por amostragem precisa **primeiro** certificar que nenhuma aresta tem
frequência zero; caso contrário o resultado é um falso contraexemplo. Isto é
o mesmo fenômeno já registrado na família híbrida (`K ≥ 250k` para os
intermediários).

---

## 2. Dois lemas que saíram da análise

### 2.1 Lema da complementação

> `G_n` é bipartido (o cavalo alterna de cor) e `|V| = n²` é par, logo
> `1_E = δ(A)` pertence a `C^⊥`. Portanto `R ↦ E∖R` é uma involução de
> `C_n^⊥` que **preserva a classe** mod `C^⊥`, e em cada coset
> `max|R| = |E| − min|R|`.

Verificado: `|R|max = 78 = 80−2` para as classes-par e `76 = 80−4` para a
classe dos 4 cantos, batendo exatamente com a busca local (que não superou o
limite em nenhuma classe). Ou seja, os pesos máximos são **certificados**, não
heurísticos.

Três consequências:

1. **(C1) do CNP sai de graça e a razão é exatamente a do §1 do plano.**
   `R ≠ G` porque `1_E ∈ C^⊥` e witnesses vivem fora de `C^⊥`. A pendência
   *"confirmar que o argumento de maximalidade não usa paridade"* fica
   resolvida no sentido útil: no caso bipartido não é preciso maximizar coisa
   alguma.
2. **O R maximal do Lema 2.1 é o complemento do detector mínimo.** Não há
   conteúdo de otimização em "tome o elemento de maior cardinalidade" aqui —
   é uma bijeção explícita com um objeto de peso 2.
3. **(C3) é uma afirmação sobre o detector mínimo:**
   `e_R(A,B) ≥ e_G(A,B)/2 ⟺ e_{Rmin}(A,B) ≤ e_G(A,B)/2`.
   Medida em cortes de vértice único: `min_v e_R(v)/deg(v) = 0.500` exato,
   **zero violações** — e a igualdade é atingida precisamente nos cantos
   (onde `deg = 2` e `Rmin` usa 1 das 2 arestas). (C3) vale, e vale
   *justo*, no lugar exato onde mora o deficit.

### 2.2 RETRATADO — "o switcher não pode viver no bulk" estava ERRADO

> ⚠️ **Esta seção continha um erro circular. Foi refutada pela Fase 3
> (`FASE3_RESULTS.md`). O texto original está preservado abaixo, riscado,
> porque o erro é instrutivo.**

O argumento era: *"a paridade `|C ∩ R| mod 2` só depende da classe de `R` mod
`C^⊥`; pelo teorema `Q(n)=3` toda classe não-trivial tem representante
suportado em `Mand`; logo todo ciclo que evita os cantos tem interseção par
com todo `R ∈ C_n^⊥`."*

**Onde quebra.** A primeira parte está certa: a paridade de fato só depende da
classe. A segunda não. `Q(n)=3` descreve as classes de
`X = C^⊥ + Span(XOR_pairs)` — e afirmar que essas são *todas* as classes de
`C_n^⊥` **é exatamente a tightness**, que é o que se quer provar. O `R` de uma
prova por contradição é, por hipótese, um `R ∈ C_n^⊥ ∖ X`; sobre a classe dele
o teorema `Q(n)=3` não diz nada.

**O que vale de verdade** (Lema B, Fase 3):

> `C_n^⊥ ∩ perp(Z_bulk) = X`.
> *Prova:* `perp(Z_bulk) = X ⊕ ⟨1_e⟩` para `e ∈ Mand` (Lema A). Se
> `R = x + 1_e` com `x ∈ X`, então `⟨τ,R⟩ = 0 + 1 = 1` para todo tour `τ`,
> pois toda obrigatória está em todo tour — logo `R ∉ C_n^⊥`. ∎

Ou seja: todo `R ∈ C_n^⊥ ∖ X` tem interseção **ímpar** com algum ciclo do
bulk. O passo (S2.a) do CNP sai **de graça** justamente para os `R` que
importam — o oposto do que esta seção afirmava.

A Fase 3 exibiu switchers completos no bulk, verificados, para
`n ∈ {6,8,10,12}`.

### 2.3 RETRATADO — o "bloqueio de grau" também cai

Dependia do §2.2 (para forçar o canto a estar em `C`) e, independentemente,
supunha que `v₁` ficaria isolado em `G ∖ (V(W)∖{v₁,v₄})`. Lendo a Def. 2.2 no
PDF original: em `W`, `v₁` tem grau 2, mas em `G_n` ele mantém todos os seus
outros vizinhos (grau até 8) — só perde os dois de `C`. O isolamento só
ocorreria se `v₁` fosse canto, o que vinha do §2.2 já refutado.

---

## 3. Impacto no plano

| Fase | Status após esta análise |
|---|---|
| Fase 0 (reformulação dual) | inalterada; ganha 2 lemas novos para incluir |
| Fase 1 (Lema 2.1 bipartido) | **simplificada** — §2.1 mostra que no caso bipartido a maximalidade é uma involução explícita, não um argumento de extremalidade |
| Fase 2 | **concluída** para n=6 (exaustiva até peso 6) e n=8 (peso ≤ 4) |
| Fase 3 (switcher no bulk) | ~~inviável~~ → **executada e bem-sucedida**, ver `FASE3_RESULTS.md` |
| Fase 4 (completar S3) | ~~vazia~~ → (S3) satisfeito em 100% dos hexágonos amostrados no 12×12 |

**A recomendação original desta nota ("não iniciar a Fase 3") estava baseada
no §2.2 errado e foi revertida.** A Fase 3 foi executada e o esquema do CNP
é instanciável em `G_n`.

---

## 4. Reprodução

```bash
python witness_map.py --n 6 --max-weight 6                 # exaustivo, ~8 s
python witness_map.py --n 8 --tours 4000 --max-weight 4    # ~8 s
```

`knight_gf2.py` traz um enumerador próprio (DFS + poda de grau + poda de
conexidade; 9.862 tours no 6×6 em ~3 s) porque a Fase 2 precisa de
**cobertura de arestas**, não de throughput. Para estender a n ≥ 10 o
caminho é plugar o motor BT v2 do paper (§ Algoritmo base: R2 + Union-Find
incremental + pressão de vértice, 214–336 tours/s no 10×10) — mas mantendo a
aleatorização rasa, senão a armadilha do §1-bis volta.

### Relação com o paper v2

`knight_tour_complete_v2.tex` §"Tightness generalizada" diz que a rota de
prova *"requer um argumento de deformação local de tour não coberto"*. Esta
nota estreita esse buraco: a deformação candidata da literatura (o
parity-switcher do CNP) está **excluída** pelo §2.2, e a exclusão decorre do
Teorema `Q(n)=3` que o próprio paper já prova.
