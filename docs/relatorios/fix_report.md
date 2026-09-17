# Relatório de Correção — `knight_tour_complete.tex`

**Base:** `agent_prompt_fix_review.md` + `review_report.md` (12 achados).
**Método:** cada achado foi reproduzido contra código/dados do repositório
*antes* de qualquer edição. Provas matemáticas necessárias (Teorema 4.3)
foram escritas como argumento combinatório curto, fiel ao que os scripts já
verificam; não foi necessário formalizar em Lean 4 (o caso é uma simples
contagem de órbitas de permutação, proporcional ao estilo do paper).

**Compilação final:** `pdflatex` 2×, exit 0, **28 páginas**, **0 refs/citações
indefinidas**, **0 overfull/underfull**.

---

## [MAIOR-1] Speedups DnC do abstract — **CONFIRMADO** (corrigido)

**Evidência:** busca exaustiva no repo pelos números 8,8 / 18,9 / 31,3 ligados
a speedup/DnC:
```
grep -rn -e "8\.8" -e "18\.9" -e "31\.3" (json/py/txt) | grep -i speed|dnc  → 0 ocorrências
```
Não existem em nenhum dado. A `tab:scaling` (0,8× / ∞ / 30×) e a `rem:dnc12`
(BT mais rápido em N=12) são internamente consistentes e marcadas `\verified`.
Os números do abstract eram *stale*.

**Diff:** abstract — substituída a frase pelos números reais de `tab:scaling`:
> "backtracking direto competitivo em N=12, timeout do backtracking em N=18 e
> ~30× em N=24 (Tabela 4.5)".

---

## [MAIOR-2] `Q(n)=3` `[Provado]` para todo `n≥6` — **CONFIRMADO** (corrigido)

**Evidência:** `Sorries.lean` confirma:
- `Q_eq_three_nX` via `native_decide` (kernel) só para `n∈{4,5,6,7,8,9,10}`;
- `bulk_connected_X` só para `n∈{6,8}`;
- SORRY #2 = `Q_eq_three_general` (n simbólico) permanece aberto.

`tab:grafo` estende a verificação GF(2) direta até n=14 (não kernel).
A prova estrutural (§2) contradizia o §6 (sorry).

**Diff (Opção A do prompt):**
- Teorema 2.5: reescrito o escopo — verificado pelo kernel para n∈{4..10};
  prova estrutural fecha para n∈{6,8,10}; estende a {12,14} via GF(2); **caso
  n≥6 geral marcado como conjectura**, pendente da conexidade do bulk.
- Abstract e sumário (§1.3) alinhados ao mesmo escopo.
- Marcador `\proved\;\verified` mantido (válido para o conjunto finito; o
  enunciado agora declara explicitamente a parte conjectural). §2 e §6 deixaram
  de se contradizer.

---

## [MAIOR-3] Família híbrida nunca alcança o plano — **CONFIRMADO** (corrigido)

**Evidência (computada, determinística):**
```
classify_edges(6): |E_torus|=144, |E_plane|=80, |W_all|=64, |W_corners|=24 (6/canto)
G(Corners) = T − W_corners:  |E|=120, β₁=85, Q=3   (Def 2.12, k=4 verdadeiro)
Plano G_n  = remover todas wraps: |E|=80, β₁=45, Q=3
família cw acumulado: k=0..4 → |E|=144,138,132,126,120 ; Q=0,0,1,2,3
```
A linha k=4 da `tab:tightness` usava o **plano** (80/45), objeto distinto de
G(Corners) (120/85). Confirmado: a interpolação Def 2.12 termina em G(Corners),
não no plano. Os scripts (`Q_locality_theorem.py`, `tightness_intermediate.py`)
constroem G(H) removendo só wraps de canto.

**rank(Ham(G(Corners))) computado** (não estimado), 400.000 tours únicos,
convergência sólida:
```
rank final = 82  → deficit = 85 − 82 = 3 = Q   (tightness vale)
progressão: …,(396000,82),(398000,82),(400000,82)
```

**Diff:**
- `tab:tightness` linha k=4 → `G(Corners)`: |E|=120, β₁=85, Q=3, rk=82,
  deficit=3, 400000† (toda a tabela agora é uma única família consistente,
  −6 arestas/passo).
- Bullets de interpolação (§2.4): k=4 = G(Corners) (4 cantos planarizados),
  **explicitamente distinto** do plano G_n; nota Q(G(Corners))=Q(G_n)=3.
- Cor. 2.17 renomeado/esclarecido: G(Corners)≠G_n como grafos, mas Q coincide.
- Caption e texto pós-tabela atualizados; plano próprio remetido às
  Tabelas `tab:grafo`/`tab:proibidos`.
- Abstract: "interpolando entre o toro e G(Corners) (cujo Q coincide com o do
  plano)".

---

## [MAIOR-4] Teorema 4.3 `[Provado]` sem prova — **CONFIRMADO** (corrigido)

**Evidência:** `dnc_t1_theorem.py` (docstring) já contém o argumento estrutural
e a hipótese decisiva: "o resultado global é UM ciclo hamiltoniano **sse** o
meta-grafo de blocos forma UM ciclo (Hamilton no meta-grafo)". O paper só
exibia "Verificado 100%: 50/50…" e escondia essa hipótese em "necessariamente
1 ciclo" (l.970).

**Diff (Opção A):**
- Enunciado: hipótese "meta-ciclo hamiltoniano único σ sobre os k² blocos"
  tornada explícita.
- Inserido bloco `\begin{proof}`: união 2-regular ⇒ ciclos disjuntos; nº de
  ciclos globais = nº de órbitas de σ; σ ciclo único ⇒ exatamente 1 tour
  hamiltoniano (36k²=N² vértices). (Caso c≥2 órbitas ⇒ c sub-ciclos.)
- l.970: "necessariamente 1 ciclo" → "dado o meta-ciclo perimetral único, …
  formam 1 ciclo hamiltoniano (Teorema 4.3)".

---

## [MAIOR-5] Klein: 1 vs 2 classes — **CONFIRMADO; verdade-base = 1** (corrigido)

**Evidência (rodado `klein_parity_survey.py`):** para n+m par, **1 única**
classe realizável (100% dos tours numa classe):
```
5×5 (V ímpar): só (1,1)   → 1/4 classes  (15.000 tours únicos, z=+211)
5×7 (V ímpar): só (1,1)   → 1/4
7×7 (V ímpar): só (1,1)   → 1/4
6×6 (V par)  : só (0,0)   → 1/4
```
**Isto refuta a premissa do review** (que assumia 2 corretas e l.1074 errada).
O dado mostra **1** classe. Portanto a linha errada é **l.1107** ("apenas 2"),
não a l.1074 ("apenas 1"). Regra do prompt cumprida: alinhar ao script.

**Diff:** §5.3 ("Estrutura de paridade") reescrito para "apenas **1** classe
((1,1)) realizável no 5×5; as 3 demais — incl. off-diagonais (0,1),(1,0) —
proibidas", com remissão explícita à Obs. `rem:klein-q` (que já dizia 1).
As duas passagens agora contam a mesma história.

---

## [MENOR-6] `tab:knuth` n=8: log 12,9 vs exato 13,12 — **CONFIRMADO** (anotado)

**Evidência:** `data/tour_count_estimates.json` → n8 `ci_boot_log10 =
[12.001, 13.256]`; estimativa pontual `log10 = 12.871`. O exato McKay
(13,12) **está dentro** do IC bootstrap 95%.

**Diff:** adicionada nota de rodapé `$^\ast$` à célula n=8: o exato McKay
(13,12) situa-se no IC bootstrap [12,0; 13,3]; a estimativa pontual subestima
levemente, como esperado de estimador de alta variância. (Não-bloqueador.)

---

## [MENOR-7] Abstract generaliza Teorema U — **CONFIRMADO** (corrigido)

**Evidência:** Teorema 2.14 (`thm:U`) é enunciado só para n=6; abstract dava
Q(H)=max(0,k−1) sem restrição.

**Diff:** abstract → "estabelecemos o Teorema U (verificado para n=6)".

---

## [MENOR-8] "Três categorias" vs 4 marcadores — **CONFIRMADO** (corrigido)

**Evidência:** existem 4 macros de status (`proved/verified/empirical/aberto`)
e resultados `\empirical` (Knuth, ratio, martingale) não cabiam nas 3.

**Diff:** §1.3 → "quatro categorias" + adicionado item (4) "Resultados
empíricos `\empirical`: Knuth, razão tours/2-fatores, martingale".

---

## [MENOR-9] Descrição de `N(8)` — **CONFIRMADO** (corrigido)

**Evidência:** 13.267.364.410.532 é o nº de tours fechados **não-direcionados**
(direcionado = 2×); "sem redução por simetria" era ambíguo.

**Diff:** "(ciclos hamiltonianos não-direcionados, sem quociente por D₄)".

---

## [ESTILO-10] `\ref` hard-coded "Lema 2.7" — **CONFIRMADO** (corrigido)

**Evidência:** `lem:bordas` = 2.7 no `.aux` (confirmado pós-compilação).

**Diff:** duas ocorrências "Lema 2.7" → `\ref{lem:bordas}` (§6.3, lista e
tabela de inventário). Renderiza "2.7"; robusto a renumeração.

---

## [ESTILO-11] Inventário "26 teoremas" não soma — **CONFIRMADO** (corrigido)

**Evidência (grep no projeto Lean, excl. .lake):** 26 declarações
theorem/lemma, das quais **4 contêm `sorry`** (`corner_degree_eq_two_general`,
`Q_eq_three_general`, `theorem_U_general`, `Q_rect_eq_three`). Logo
**22 sem sorry**. Quebra real:
```
grau dos cantos = 2 : 11   (n4,n5 + 4×n6 + 4×n8 + n10)
Q(n)=3 (n=4..10)    :  7
conectividade bulk  :  2   (n=6,8)
propriedades n=6    :  2   (lemma_2_7_corner00_n6, num_xor_vectors_n6)
TOTAL               : 22
```
O "Teorema U (k=0..4)" e "Q retangular" do inventário antigo são **`#eval`**
(não teoremas de kernel) — não contam.

**Diff:**
- Box `lake build`: "26 teoremas sem sorry" → "**22**".
- Tabela de inventário: rows refeitos (11+7+2+2=**22**), com linha Total;
  removidas as linhas U(~5) e retangular(3); "~" eliminado.
- Bullets §6 de U e retangular anotados "(via `#eval`)"; nota pós-tabela
  esclarece que não entram na contagem de 22.

---

## [DÚVIDA-12] Reivindicação de primazia — **PENDENTE** (não editado)

Conforme instrução, **não editado**. A frase já contém "ao nosso
conhecimento" (l.236-241). **Decisão do autor**: manter, atenuar mais, ou
buscar trabalhos de knight-tour em Lean/Coq para confirmar. Sem ação.

---

## Verificação de critérios de PRONTO

- [x] `fix_report.md` cobre os 12 itens.
- [x] `.tex` recompila limpo (28 p., 0 ref/cit indefinida, 0 overfull).
- [x] §2 (Teorema 2.5) e §6 (sorry) não se contradizem mais.
- [x] Abstract não cita speedups inexistentes nem Teorema U geral.
- [x] Nenhuma linha de tabela mistura híbrido vs plano sem rótulo
      (k=4 = G(Corners), rotulado; plano tratado à parte).
- [x] Todo número alterado tem origem rastreável: rank=82/deficit=3 de
      `tightness_intermediate.py`+`Q_locality_theorem.py` (400k tours);
      Klein 1-classe de `klein_parity_survey.py`; 22 teoremas de grep Lean;
      IC Knuth de `data/tour_count_estimates.json`. **Nenhum número inventado.**
