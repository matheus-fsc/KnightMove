# Agent Prompt — Verificação e Correção dos Achados da Revisão (`knight_tour_complete.tex`)

## CONTEXTO

Projeto `knight_tour/`: estudo de invariantes topológicos do passeio do
cavalo via espaço de ciclos sobre GF(2). O artigo principal é
`knight_tour_complete.tex` (1704 linhas, pt-BR, compila com
`pdflatex knight_tour_complete.tex`, 2 passadas). Autor: Matheus Coelho.

Uma revisão cética anterior (`review_report.md`) levantou **12 achados**.
Sua missão **não é re-revisar do zero** — é **verificar empiricamente cada
achado contra o código/dados do repositório** e, onde confirmado, **aplicar a
correção mínima no `.tex`**. Onde a verdade-de-base for ambígua ou exigir
decisão do autor, **NÃO edite**: registre como pendência.

Marcadores de status no paper (macros): `\proved` (prova completa),
`\verified` (evidência computacional), `\empirical` (medição ruidosa),
`\openp` (aberto).

---

## REGRAS

1. **Verifique antes de editar.** Cada item abaixo traz uma *fonte-de-verdade*
   (script ou JSON). Rode-a / leia-a e compare com o que o `.tex` afirma.
   Nunca "corrija" um número sem antes reproduzi-lo ou achar sua origem.
2. **Edição mínima e cirúrgica.** Mude só o token/linha necessário. Não
   reescreva parágrafos nem mexa em formatação alheia ao achado.
3. **Recompile ao final** (`pdflatex` 2×) e confirme: 0 refs/citações
   indefinidas, 0 overfull novo, 28 páginas (±1).
4. **Não invente dados.** Se um benchmark não existe no repo, a correção é
   *remover/atenuar a alegação*, não fabricar números.
5. Entregue `fix_report.md`: por item → `[CONFIRMADO/REFUTADO/PENDENTE]`,
   evidência (comando + saída resumida), e o diff aplicado (ou por que não).

---

## ITENS A VERIFICAR E CORRIGIR

Severidade entre colchetes. Itens 1–5 são os que mais importam.

### [MAIOR-1] Speedups do DnC no abstract contradizem o corpo
- **Afirmação suspeita (abstract, l.88–90):** speedup DnC **8,8× (N=12),
  18,9× (N=18), 31,3× (N=24)**.
- **Conflito:** Tabela `tab:scaling` (l.985–987) dá **0,8× / ∞ / 30×**, e a
  Observação `rem:dnc12` (l.992–1000) diz que **BT é mais rápido** no N=12.
- **Fonte-de-verdade:** `dnc_t3_scaling.py`, `dnc_q3_18x18.py`,
  `dnc_t2_24x24.py`; procurar JSON de saída em `data/` e raiz
  (`grep -rn "speedup" --include=*.json .`). Rodar `dnc_t3_scaling.py` se
  barato; senão localizar o JSON que gerou `tab:scaling`.
- **Verificação:** os números 8,8/18,9/31,3 existem em algum dado? (busca
  inicial NÃO os encontrou em `data/*.json` — provavelmente *stale*.)
- **Regra de decisão:** se 8,8/18,9/31,3 não têm fonte reproduzível →
  **substituir no abstract pelos números de `tab:scaling` (0,8×/∞/30×)** ou
  reescrever a frase para "speedup que cresce com N (BT competitivo em N=12,
  timeout em N=18, ~30× em N=24)". Manter consistência com `rem:dnc12`.

### [MAIOR-2] `Q(n)=3` marcado `[Provado]` para todo `n≥6`, mas a prova só fecha até n=10
- **Afirmação suspeita:** Teorema 2.5 (l.320–322) `\proved\;\verified`
  "para todo `n≥6`". Sumário l.167 repete "(provado)".
- **Conflito:** a prova do Lema 2.10(a) (l.380–409) depende da conexidade do
  bulk, **explicitamente deixada para trabalho futuro para `n≥10`**
  (l.392–397). O Lean admite isso: `sorry` #2 = "Teorema Q(n)=3 geral"
  (l.1479–1483). Direta auto-contradição §2 vs §6.
- **Fonte-de-verdade:** `Q_locality_theorem.py`, `deficit_theorem/` (verificar
  até que `n` o `Q=3` é de fato *computado* — `tab:grafo` afirma n=6..14 via
  GF(2)); `knight_tour_lean/.../Sorries.lean` (confirmar que sorry #2 é o
  geral).
- **Verificação:** (a) confirmar que `Q=3` é computado direto (não só inferido)
  para n∈{6,8,10,12,14}; (b) confirmar que a conexidade do bulk só é provada
  p/ n∈{6,8} (Lean) + verificada n=10.
- **Regra de decisão (escolher 1, sem inventar prova):**
  - **Opção A (recomendada, mínima):** trocar o escopo do `\proved`. Enunciar
    "`Q(n)=3` provado para `n∈{6,8,10}` (Lema 2.10 + conexidade do bulk
    verificada); verificado por GF(2) direto até `n=14` (Tab. `tab:grafo`); o
    caso `n` geral permanece **conjectura** pendente da conexidade do bulk."
    Ajustar marcador do Teorema 2.5 e a linha do sumário (l.167)
    correspondentemente.
  - **Opção B:** manter `\proved` mas adicionar hipótese explícita ao
    enunciado ("supondo `G_n∖Mand` conexo, fato verificado p/ `n≤10`").
  - Qualquer opção deve **deixar §2 e §6 coerentes** (parar de afirmar
    "provado p/ todo n" num lugar e "sorry" no outro).

### [MAIOR-3] A família híbrida nunca alcança o plano (contagem de arestas)
- **Afirmação suspeita:** Def 2.12 (l.480–486) — `G(H)` remove só
  *wrap-edges incidentes a cantos*. Tabela `tab:tightness` (l.640–644) lista
  `k=4` com **|E|=80, β₁=45** (= o plano real). Cor 2.17 (l.567) afirma
  `Q(G_n)=Q(G(Corners))`.
- **Conflito aritmético:** padrão −6 arestas/canto dá k=0..4 =
  144,138,132,126,**120** — não 80. O plano remove *todas* as wraps de borda
  (inclusive de vértices não-canto). `G(Corners)` (|E|=120) ≠ `G_n` (|E|=80).
- **Fonte-de-verdade:** `Q_locality_theorem.py` e/ou `tightness_intermediate.py`
  (constroem o grafo híbrido). **Rodar e imprimir `|E|` para H = cada canto
  acumulado** e para H = todos os 4 cantos.
- **Verificação:** confirmar `|E|(G(Corners))`. Se = 120 → a Def gera "toro
  com 4 cantos planarizados", não o plano. Se o script na verdade remove mais
  arestas (chega a 80) → a Def 2.12 no texto está incompleta (descreve mal o
  que o código faz).
- **Regra de decisão:**
  - Se `|E|(G(Corners))=120`: **corrigir a linha k=4 de `tab:tightness`** para
    o grafo híbrido verdadeiro (|E|=120, β₁=85, e recomputar Q/rk/deficit com
    o script) **OU** marcar a linha k=4 explicitamente como "plano (objeto
    distinto, fora da interpolação)" e ajustar Cor 2.17 para não igualar
    `G(Corners)=G_n`. Ajustar também a redação do abstract ("interpola entre
    toro e plano") se a interpolação não chega ao plano.
  - Se o script chega a 80: **corrigir a Def 2.12** para descrever a remoção
    real (todas as wraps, não só as dos cantos).
  - **Não chutar** β₁/Q/deficit do ponto k=4 verdadeiro: computar com o script.

### [MAIOR-4] Teorema 4.3 (Suficiência Local) marcado `[Provado]` sem prova escrita
- **Afirmação suspeita:** Teorema 4.3 (l.925–931) `\proved\;\verified`, mas
  só há "Verificado 100%: 50/50…" (l.933). Linha 970 ("4 caminhos + 4
  cross-edges formam **necessariamente 1** ciclo hamiltoniano") esconde a
  hipótese de meta-ciclo único.
- **Fonte-de-verdade:** `dnc_t1_theorem.py` (nome sugere o teorema),
  `dnc_t4_verify.py`, `dnc_q2_success.py`.
- **Verificação:** existe prova/argumento formal nos scripts ou só verificação
  empírica? A hipótese "meta-grafo dos blocos é um único ciclo" está explícita?
- **Regra de decisão:**
  - **Opção A:** inserir prova curta (3–5 linhas): meta-ciclo hamiltoniano nos
    `k²` blocos + caminho hamiltoniano compatível por bloco ⇒ concatenação é 1
    ciclo hamiltoniano (cada cross-edge liga fim de um a início do próximo na
    ordem do meta-ciclo). Tornar a hipótese "meta = ciclo único" explícita no
    enunciado.
  - **Opção B:** rebaixar para `\verified` e remover "necessariamente 1 ciclo"
    (l.970–971), substituindo por "1 ciclo, *dado* um meta-ciclo hamiltoniano".

### [MAIOR-5] Klein: contradição interna 1 vs 2 classes realizáveis
- **Afirmação suspeita:** l.1074 "Com `n+m` par, **apenas 1** das 4 classes é
  realizável"; l.1107–1108 "**apenas 2** das 4 classes são realizáveis; (0,1)
  e (1,0) sempre proibidas".
- **Fonte-de-verdade:** `klein_parity_survey.py`, `knight_tours_klein.py`,
  `klein_full_search.py`. Memória do projeto diz: off-diag (0,1)(1,0) sempre
  proibido ⇒ **2 classes** realizáveis; `n+m` par ⇒ Q_top=1.
- **Verificação:** rodar/ler o survey de paridade: quantas das 4 classes
  `(w_x mod 2, w_y mod 2)` aparecem para 5×5 (e p/ `n+m` par em geral)?
- **Regra de decisão:** alinhar as duas linhas ao resultado do script. Se 2
  realizáveis (esperado): **corrigir l.1074** ("1" → "2", ou reformular para
  "1 *classe proibida a mais* além das de borda", esclarecendo a relação com
  Q_top=1). Garantir que `rem:klein-q` e §5.3 contem a mesma história.

### [MENOR-6] `tab:knuth` n=8: log 12,9 vs exato 13,12
- **l.1132–1133:** estimativa `n=8 ~10¹³, log 12,9`; valor McKay exato
  `1,327×10¹³`, log10 = **13,12**.
- **Fonte-de-verdade:** `data/tour_count_estimates.json`,
  `tour_count_estimator.py`.
- **Regra:** se o IC bootstrap contém 13,12, manter mas anotar; se a célula é
  só ponto-estimado, atualizar para o IC ou nota de rodapé. Não é bloqueador.

### [MENOR-7] Abstract generaliza o Teorema U além de n=6
- **l.81–83:** abstract enuncia `Q(H)=max(0,k−1)` sem restrição; Teorema 2.14
  (l.494–496) é só `n=6`.
- **Regra:** qualificar o abstract: "(Teorema U, verificado para `n=6`)".
  Edição de texto pura, sem dados.

### [MENOR-8] "Três categorias" vs 4 marcadores
- **l.165:** sumário diz "três categorias" mas há 4 status
  (`proved/verified/empirical/aberto`) e resultados `\empirical`
  (martingale/Knuth/ratio) não cabem.
- **Regra:** mencionar a 4ª categoria (empírica) ou reclassificar. Texto puro.

### [MENOR-9] Descrição imprecisa de `N(8)`
- **l.201–204:** `N(8)=13.267.364.410.532` dito "sem redução por simetria".
  Esse é o número **não-direcionado**; "sem simetria" sugere sem quociente por
  `D₄`.
- **Regra:** precisar a redação ("tours fechados não-direcionados, sem
  quociente por `D₄`"). Confirmar a semântica com a literatura McKay antes.

### [ESTILO-10] `\ref` hard-coded "Lema 2.7"
- **l.1409, 1435:** "Lema 2.7" escrito à mão (hoje casa: `lem:bordas`=2.7).
- **Regra:** trocar por `\ref{lem:bordas}`. Recompilar e confirmar que rende
  "2.7".

### [ESTILO-11] Inventário "26 teoremas" não soma
- **l.1432–1436:** 11+7+2+3+~5+3 = ~31, não 26; e o "~5" é frouxo.
- **Fonte-de-verdade:** contar teoremas reais em `knight_tour_lean/` (grep
  `theorem`/`lemma` sem `sorry`).
- **Regra:** acertar os números para somar ao total real; remover o "~".

### [DÚVIDA-12] Reivindicação de primazia (não confirmar)
- **l.236–241:** "primeira formalização em Lean 4 do espaço de ciclos GF(2) do
  cavalo".
- **Regra:** **NÃO editar para confirmar.** No máximo sugerir atenuação ("ao
  nosso conhecimento" já presente). Registrar como decisão do autor.

---

## ORDEM DE EXECUÇÃO

1. Verificar 1–5 (MAIORES) com os scripts/JSON antes de tocar no `.tex`.
2. Aplicar edições mínimas dos itens **CONFIRMADOS** (1–11). Itens com
   "Opção A/B" → escolher a opção *mínima* que torne o texto **honesto e
   internamente coerente**; se a correção exigir recomputar números, **rodar o
   script e usar o valor real** (jamais estimar).
3. Item 12: só registrar.
4. `pdflatex` 2×; conferir refs/citações/overfull/páginas.
5. Escrever `fix_report.md` com, por item: status, evidência (comando+saída
   curta), diff aplicado (ou motivo de PENDENTE).

## CRITÉRIO DE PRONTO

- `fix_report.md` cobre os 12 itens.
- `.tex` recompila limpo; nenhuma nova ref/citação indefinida.
- §2 (Teorema 2.5) e §6 (Lean/sorry) **não se contradizem mais**.
- Abstract não cita números (speedup, Teorema U geral) que o corpo não
  sustenta.
- Nenhuma linha de tabela mistura objetos distintos (híbrido vs plano) sem
  rótulo.
- Nenhum número novo foi inventado: todo valor alterado tem origem rastreável
  num script/JSON citado no `fix_report.md`.
