# Item 4 — o que o "perfil" do D&C realmente codifica

**Data:** 2026-08-07 · **Fonte:** `knight_tours_dnc.py`, `dnc_t1_theorem.py`, `dnc_q3_18x18.py`

Pergunta: *o perfil codifica emparelhamento ou só multi-conjunto?*

**Resposta: nenhum dos dois — é mais fraco.** O catálogo é indexado por um
**único par ordenado `(s, e)`** de vértices de fronteira, e cada bloco 6×6
contribui **exatamente um** caminho hamiltoniano, de `s` a `e`
(`build_catalog`, linha 253). Não há emparelhamento porque não há mais de uma
travessia por bloco a emparelhar.

## O que isso implica — e a boa notícia

O teorema T1 (`dnc_t1_theorem.py`) enuncia:

> dados `k²` blocos 6×6, cada um contribuindo 1 caminho hamiltoniano com
> extremos `(s_i, e_i)`, ligados por cross-edges, o resultado global é **um**
> ciclo hamiltoniano **sse** o multigrafo de blocos induzido pelas cross-edges
> forma um único ciclo.

O docstring já contém o argumento estrutural, e **ele está correto**:

- cada cross-edge liga `e_i` a `s_{σ(i)}`, logo `σ` é uma permutação;
- no meta-grafo cada bloco tem grau 2 (uma entra, uma sai) ⇒ união disjunta de
  ciclos, um por ciclo de `σ`;
- a união cobre todos os vértices porque cada caminho de bloco é hamiltoniano
  no seu bloco;
- global é um ciclo só ⟺ `σ` é um `k²`-ciclo.

**Isto não é empírico.** É álgebra de permutações elementar; os 1639 trials
verificam a *implementação*, não a matemática. A preocupação de que a
suficiência do D&C repousasse sobre evidência não se confirma — ela repousa
sobre um argumento de duas linhas que só não foi escrito como teorema.

Soundness é o que (U1) precisa (queremos **exibir** um caminho, não enumerar
todos), então a incompletude — o catálogo só gera tours da forma "uma passagem
por bloco" — não é problema para esta aplicação.

## Os dois problemas de escopo que ninguém tinha levantado

### (A) `N = 6k` apenas

`N_GLOBAL = 12` está fixo em `knight_tours_dnc.py:39`, com grade 2×2 de blocos
6×6; `dnc_q3_18x18.py` é um módulo separado com grade 3×3. A decomposição exige
`N` múltiplo de 6.

Os tabuleiros onde a tightness está verificada são **n ∈ {8, 10, 12, 14}** — e
só o 12 é múltiplo de 6. A família livre de borda foi verificada em
n = 10, 12, 14; o catálogo, como existe, **não se aplica a n = 8, 10 nem 14**.

Para servir a (U1) em `n` par geral, a decomposição precisa de blocos de
tamanhos mistos (p.ex. uma faixa 6×m mais blocos 6×6), o que muda o catálogo:
deixa de ser um catálogo e vira uma família de catálogos indexada pelos
tamanhos de bloco usados.

### (B) (S3) pede caminho, não ciclo

O teorema T1 produz um **ciclo** hamiltoniano global. (S3) precisa de um
**caminho** hamiltoniano de `v₁` a `v₄` no tabuleiro **menos o patch**, com
extremos prescritos e ambos dentro do bloco perfurado.

Estruturalmente:

- o bloco perfurado contribui **dois** segmentos (`v₁ → saída` e
  `entrada → v₄`), não um;
- portanto ele tem duas "portas" e o meta-grafo deixa de ser 2-regular;
- a condição "meta-grafo é um único ciclo" precisa ser substituída por uma
  condição sobre um meta-grafo com um nó desdobrado.

É uma variante do mesmo argumento, provavelmente igualmente elementar, mas
**não é o teorema que existe**, e é o teorema de que o catálogo perfurado
precisa.

## Veredicto

| preocupação | status |
|---|---|
| suficiência do D&C é empírica | ❌ **infundada** — o argumento existe e é correto |
| perfil codifica emparelhamento | ❌ é par único `(s,e)`, mais fraco que ambas as opções |
| aplicável aos `n` verificados | ⚠️ **só n=12** dos quatro |
| cobre (S3) (caminho com extremos, patch) | ⚠️ **não** — precisa da variante com nó desdobrado |

A rota do catálogo **não precisa de outro mecanismo de suficiência**; precisa de
duas extensões de escopo, ambas aparentemente elementares, mas nenhuma escrita.
Custo revisado: a semana do catálogo continua uma semana, **mais** o trabalho de
generalizar a decomposição para `n` não-múltiplo de 6, que é o item de risco
real e não estava no orçamento.
