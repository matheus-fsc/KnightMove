# Pesca espectral na T_bulk(6) real — veredito

**Script:** `pesca_espectral_run.py` · **Resultado:** `data/pesca_espectral_n6.json`
**Matriz:** `data/transfer_n6_bulk.npz` (T_bulk coluna-a-coluna, estado
`(deg_c, deg_{c+1}) ∈ {0,1,2}⁶×{0,1,2}⁶`).

## PASSO 0 — reconhecimento (o que é o quê)

A premissa "T_bulk(6) construída pelo `knight_transfer.cpp`" está **incorreta** e
foi a primeira coisa a corrigir:

| objeto | arquivo | estado | saída | tamanho |
|---|---|---|---|---|
| broken-profile DP | `knight_transfer.cpp` | `PState{lo,hi}`, 4 bits/célula (0=vazio, 15=fechado, 1–14=plug-id), canonicalizado | **escalar** N(6)=9862 | "187M estados" em n=8 |
| **column-transfer T_bulk** | `transfer_build.py` → `transfer_n6_bulk.npz` | `(deg_c, deg_{c+1})` perfil de graus | **matriz** 33346×33346 | nnz 1.966.428 |

São DPs **diferentes**. O `knight_transfer.cpp` nunca materializa matriz alguma —
faz um sweep célula-a-célula acumulando contagem. A matriz das âncoras (33.346
estados, λ₁=70,48) é a column-transfer. O export dos `.txt` foi feito a partir
do `.npz` real (o leitor `pesca_espectral` espera `(row,col,val)` 0-indexados;
o formato casa). O check de fidelidade do stub original (`T.diagonal().sum()`)
era placeholder **errado** — a contagem de 2-fatores não é o traço; é uma
contração de fronteira (ver PASSO 2).

## PASSO 2 — fidelidade (gate) — TODAS confirmadas

| âncora | obtido | esperado |
|---|---|---|
| K (estados) | 33.346 | 33.346 |
| nnz | 1.966.428 | 1.966.428 |
| 2-fatores `e_{s0}·Tᵀ⁴·T_nc2ᵀ·e_target` | **36.236** | 36.236 |
| λ₁ | **+70,4766** | +70,4766 |
| λ₂ | **−39,0867** | −39,0867 |
| gap = \|λ₁\|−\|λ₂\| | 31,3899 | ~31,39 |
| ρ = \|λ₂\|/\|λ₁\| | 0,5546 | ~0,5546 |

**N(6)=9862 (tours conexos): NÃO extraível desta matriz.** O estado é só o
perfil de *graus* de coluna — é cego à conectividade e conta todos os 36.236
2-fatores, sem distinguir o sub-tour conexo. (Coerente com a tese do projeto:
conectividade é a obstrução global que nenhuma estrutura local captura.)

## PASSO 3 — a pesca

### Há degenerescência além da paridade conhecida? **NÃO.**
- Top-60 autovalores: **0 degenerescências** (nenhum valor complexo repetido com
  multiplicidade > 1, tol 1e-6).
- 24 pares conjugados `λ, λ̄` (mesmo \|λ\|, valores distintos) — **artefato de
  matriz real**, não simetria escondida. Os "empates de módulo" (ex.
  −6,637±38,32i, ambos \|λ\|=38,889) são apenas conjugação.

### Quais operadores comutam com T? Qual fator de redução?

| operador S | ‖ST−TS‖ | comuta? | fator de redução |
|---|---|---|---|
| **(a)** paridade (r+c) diagonal `diag(±1)` | 8,0 | **NÃO** | — |
| **(b/c)** reflexão vertical `R_v: r→n−1−r` | **0,0** | **SIM** | **1,972** (≈2) |

**(a) Paridade diagonal — não comuta, e é estrutural:** T_bulk é irredutível
(Perron λ₁ real positivo simples ⇒ grafo de transição fortemente conexo). Para
`S=diag(σ)`, `ST=TS` exige `σ(s)=σ(s')` em **toda** aresta de transição ⇒ σ
constante. Logo **nenhum** funcional de paridade diagonal não-trivial comuta. A
"bipartição (r+c)" não é simetria estática: é uma graduação ℤ₂ que depende do
índice de coluna — a própria T a desloca a cada passo, e é por isso que λ₂<0
**sem** o espectro ser simétrico sob λ→−λ (não há `S` anticomutante exato).

**(b/c) Reflexão vertical R_v — comuta exatamente:** involução bem-definida nos
33.346 estados (inverte os deg-tuples), `‖P·T−T·P‖ = 0`. Decomposição:
- estados fixos (deg-tuples palíndromos): **468**
- setor simétrico: (K+fix)/2 = **16.907** ← onde vive o Perron
- setor antissimétrico: (K−fix)/2 = **16.439**
- Perron R_v-invariante: `‖P·v₁−v₁‖ = 1,4e-15` ✓ (como Perron-Frobenius exige)
- **fator de redução = K/16.907 = 1,972**

**D₄ vs sweep de coluna:** só `{id, R_v}` preservam a direção do sweep. A
reflexão horizontal (c→n−1−c) e a rotação 180° **invertem** o sweep ⇒ mapeiam
T→Tᵀ (não comutam); as reflexões diagonais e rotações 90° trocam linha↔coluna ⇒
destroem a estrutura coluna-a-coluna. Sobra o subgrupo **ℤ₂**, fator 2 no máximo.

## PASSO 4 — veredito honesto

**Resultado: NEGATIVO LIMPO.** A única simetria de transferência da T_bulk(6) é
a reflexão vertical R_v (ℤ₂, fator 1,97). Não há invariante escondido: zero
degenerescência espectral genuína, paridade-diagonal provadamente impossível, e
nenhum outro elemento de D₄ sobrevive ao sweep de coluna. **Isto encerra a busca
por invariante de transferência não-trivial** — o que reduz o espaço de estados
é exatamente a simetria geométrica óbvia, nada além.

### Estimativa concreta para n=8 (sem otimismo)

- O fator é **constante = 2** (subgrupo ℤ₂), não muda a classe de complexidade.
- O "187M estados" do brief é do **broken-profile** (`knight_transfer.cpp`), não
  da column-transfer. Para a column-transfer, K cresce 342→7.322→33.346
  (n=4,5,6); K(8) cai na faixa de ~10⁶ (perfis de grau reachable ⊂ 3¹⁶≈43M).
- R_v bloco-diagonaliza ⇒ maior bloco ≈ K/2. Isso ajuda o **eigensolve/memória**
  por ~2×, mas **não** ajuda o gargalo real, que é a **enumeração BFS dos
  estados** (já estoura 8 GB em pure-Python para n=8, conforme
  `transfer_matrix/README`). Você precisa enumerar todos os estados antes de
  poder quocientar por R_v.
- **Conclusão:** R_v sozinho **não destrava o n=8**. Ele dá um 2× de folga no
  solver, mas n=8 continua exigindo C/numba para a enumeração. n=8 exato é
  alcançável **com** essa reescrita (e R_v ajuda na margem); **n=12 permanece
  fora de alcance** — fator 2 não toca a explosão combinatória.

### Invariante de transferência ≠ obstrução Q
- **Invariante de transferência** (este relatório): operador S que comuta com
  T_bulk e bloco-diagonaliza o espaço de estados. Encontrado: só R_v (fator 2).
- **Obstrução Q=3** (`cycle_space_hyperplanes/`): restrição linear estática
  sobre F₂ no espaço de ciclos Z₁, já saturada em 3. É um objeto **diferente** —
  não é alvo desta pesca e não interage com o espectro de T_bulk.
