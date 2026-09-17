# Relatório de Revisão Crítica — `knight_tour_complete.tex`

**Revisor:** auditoria cética (Claude)
**Alvo:** `knight_tour_complete.tex` (1704 linhas) / `.pdf` (28 páginas, compilação limpa)
**Data:** 2026-05-30

---

## 1. Veredito em uma linha

**Núcleo matemático sólido, mas com 1 furo de prova no resultado central e
exageros de calibração: ~4 achados MAIORES (abstract contradiz o corpo no
DnC; `Q(n)=3` marcado `[Provado]` para todo `n≥6` quando a prova só fecha
para `n∈{6,8,10}`; família híbrida nunca alcança o plano; Teorema de
Suficiência `[Provado]` sem prova escrita) + 1 contradição interna na Klein +
menores.** Sem bloqueador formal de compilação; o "bloqueador" é de
honestidade, não de LaTeX.

---

## 2. Tabela de achados

| sev | linha(s) .tex | seção | descrição | sugestão |
|-----|---------------|-------|-----------|----------|
| **[MAIOR]** | 88–90 vs 985–999 | abstract vs §4.3 | Abstract: speedup DnC **8,8× (N=12), 18,9× (N=18), 31,3× (N=24)**. Tabela `tab:scaling`: **0,8× / ∞ / 30×**. Rem. `dnc12` (l.994) diz explicitamente que **BT é mais rápido** que DnC no N=12. Os números do abstract não existem no corpo. | Reconciliar: ou corrigir o abstract para os números reais (0,8/∞/30), ou explicar de onde vêm 8,8/18,9/31,3 e citar a fonte. |
| **[MAIOR]** | 320–322, 392–397 | Thm 2.5 + Lema 2.10(a) | `Q(n)=3` marcado `\proved` "para todo `n≥6`". A prova depende da conexidade do bulk, **provada só para `n≤10` (computacional)**; para `n≥10` geral o próprio texto diz "*deixada para trabalho futuro*". O Lean confirma: sorry #2 = "Teorema Q(n)=3 geral" (l.1479–1483). | Downgrade do enunciado: `\proved` para `n∈{6,8,10}`, `\verified` (GF(2) direto) até `n=14`, **conjectura** para `n` geral. Ou destacar a hipótese de conexidade do bulk como premissa explícita. |
| **[MAIOR]** | 480–492, 567, 640–644 | Def 2.12, Cor 2.17, `tab:tightness` | `G(H)` remove só wrap-edges **incidentes a cantos** (−6 arestas/canto: 144→138→132→126). O ponto `k=4` deveria ter **|E|=120**, mas a tabela usa **|E|=80, β₁=45 = o plano real**. O plano remove *todas* as wraps de borda (inclusive de vértices não-canto), não só as dos cantos. A interpolação **nunca chega ao plano**; `Q(G_n)=Q(G(Corners))` conflaciona dois grafos distintos. | Ou (a) redefinir `G(H)` para de fato interpolar até o plano, ou (b) admitir que `k=4` no esquema híbrido é "toro com 4 cantos planarizados" (|E|=120) ≠ plano, e calcular tightness nesse grafo. A linha `k=4` atual mistura objetos. |
| **[MAIOR]** | 925–934, 966–971 | Thm 4.3 | "Suficiência Local" marcado `\proved\;\verified` mas **não há bloco de prova** — só "Verificado 100%: 50/50…". Além disso "4 caminhos + 4 cross-edges formam **necessariamente 1** ciclo hamiltoniano" (l.970) só vale se o meta-grafo de blocos for um único ciclo; em geral poderiam formar sub-ciclos (justamente a parte difícil). | Escrever a prova (curta: meta-ciclo hamiltoniano nos blocos + path compatível por bloco ⇒ 1 ciclo) explicitando a hipótese "meta = ciclo único", ou downgrade para `\verified`. |
| **[MAIOR]** | 1074 vs 1107–1108 | Rem 5.2 vs §5.3 | Contradição interna: l.1074 "Com `n+m` par, **apenas 1** das 4 classes é realizável"; l.1107 "**apenas 2** das 4 classes são realizáveis; (0,1) e (1,0) sempre proibidas". 1 ≠ 2. | Reconciliar. Se (0,1),(1,0) são proibidas, sobram 2 realizáveis — corrigir l.1074. |
| **[MENOR]** | 1131–1134 vs 202 | `tab:knuth` vs intro | Estimador Knuth lista `n=8 ~10¹³, log 12,9`. Valor exato de McKay é `1,327×10¹³`, **log 13,12**. A estimativa subestima o exato conhecido. | Notar que é estimativa com IC; ou ajustar para o IC conter 13,12. |
| **[MENOR]** | 81–83, 494–496 | abstract vs Thm 2.14 | Abstract apresenta o Teorema U `Q(H)=max(0,k−1)` **sem restrição de `n`**; o teorema é enunciado/provado **só para `n=6`**. | Qualificar o abstract: "Teorema U (verificado para `n=6`)". |
| **[MENOR]** | 165–176 | §1.3 sumário | Diz "três categorias", mas há **4 marcadores** de status (`proved/verified/empirical/aberto`); resultados `\empirical` (martingale, Knuth, ratio) não cabem nas 3. | Mencionar a 4ª categoria (empírica) ou reclassificar. |
| **[MENOR]** | 201–204 | §1.4 contagem | `N(8)=13.267.364.410.532` descrito como "sem redução por simetria". Esse é o número **não-direcionado** (direcionado = 2×); "sem simetria" costuma significar sem redução por `D₄`. Frase ambígua. | Especificar: "tours fechados não-direcionados (sem quociente por `D₄`)". |
| **[DÚVIDA]** | 236–241 | §1.4 formalização | Reivindicação de primazia: "primeira formalização em Lean 4 do espaço de ciclos GF(2) do cavalo". | Autor deve verificar/atenuar ("ao nosso conhecimento" já ajuda; considerar buscar trabalhos de knight-tour em Lean/Coq). |
| **[ESTILO]** | 1409, 1435 | §6.3 Lean | "Lema 2.7" e "Lema 2.7" escritos **hard-coded** em vez de `\ref{lem:bordas}`. Hoje casa (lem:bordas = 2.7), mas quebra se a numeração mudar. | Trocar por `\ref{lem:bordas}`. |
| **[ESTILO]** | 1436 | tab "26 teoremas" | "Teorema U … `$\sim$5`" — o "~5" num inventário que soma para "26" é frouxo (11+7+2+3+~5+3 = ~31, não 26). | Dar o número exato; a soma não fecha em 26. |

---

## 3. Auditoria prova-a-prova (resultados `\proved`)

**Lema 2.6 — Grau dos cantos = 2 (l.330).** ✔ Fecha. WLOG `(0,0)`; dos 8
offsets de cavalo só `(1,2),(2,1)` caem em `[0,n)²`. Correto para `n≥4`.

**Lema 2.7 — Relação aresta–canto (l.342).** ✔ Fecha. `[∂(δ_c)]_e = [c∈e]`;
pelo Lema 2.6 só `e₁(c),e₂(c)` contêm `c`. Limpo.

**Corolário 2.8 (l.354).** ✔ Imediato do Lema 2.7.

**Lema 2.10 — Independência dos 4 representantes (l.366).** ⚠ **A lógica
fecha, mas há um furo de cobertura.** A parte (a): assume `Σ r_c = 0`, deriva
`∂α = Σ e₁(c)`, usa **conexidade do bulk** ⇒ `α` constante `= a` no bulk; daí
`[∂α]_{e₁(c)} = α[c]+a = 1` força `α[c]=1+a`, e então
`[∂α]_{e₂(c)} = (1+a)+a = 1`, contradizendo `e₂(c) ∉ supp(∂α)` (vizinhanças
disjuntas para `n≥6`). **Esse encadeamento está correto** — verifiquei que vale
para `a∈{0,1}` e que `e₁(c),e₂(c)` vão ambos para vértices do bulk
(`(1,2),(2,1)` são não-canto para `n≥4`). **O furo é a premissa**: "bulk
conexo após remover `Mand`" é provado só para `n≤10` (BFS no Lean); l.394–397
admite que `n≥10` geral "*requer enumeração dos tipos de vértice e é deixada
para trabalho futuro*". Logo o lema (a) **não está provado para `n≥12`**.

**Lema 2.11 — Imagem dos pares XOR (l.411).** ✔ Fecha. Casos (a) mesmo canto
⇒ `v_{ij}∈R(n)`; (b) cantos distintos ⇒ `r_{c_a}+r_{c_b}`. Correto.

**Teorema 2.5 — `Q(n)=3` (l.320).** ⚠ **Fecha condicionalmente.** Dado o Lema
2.10, o funcional `σ(Σ a_k r_{c_k}) = Σ a_k` tem `ker` de codimensão 1 em
`W≅F₂⁴`; os pares estão em `ker(σ)` e `{r₁+r₂, r₁+r₃, r₁+r₄}` são LI ⇒
`Q=3`. **O passo `≥3` e `≤3` fecha sem buraco** *desde que* os 4
representantes sejam independentes — ou seja, herda integralmente o furo do
Lema 2.10. **Conclusão: provado para `n∈{6,8,10}`; verificado por GF(2) direto
até `n=14`; em aberto para `n` geral.** O marcador `\proved` "para todo `n≥6`"
**sobrevende**. (A Observação 2.9 sobre a origem funcional do "−1" está
correta e é uma boa intuição.)

**Teorema 2.14 — Teorema U (l.494).** ⚠ Enunciado **só para `n=6`**, e nesse
escopo a prova por casos fecha (k=0 trivial; k=1 via Lema 2.7; k≥2 herda a
independência, verificada nas 16 configs). **Mas** a Corolário 2.17 ("plano
como extremal", l.567) usa `G(Corners)=G_n`, que é **falso** sob a Def 2.12
(ver achado MAIOR de `|E|`). O teorema em si (n=6) está ok; o corolário
extrapola para um grafo que a definição não produz.

**Proposição 2.18 — `Q(n,m)=3` retangular (l.578).** ✔ Corretamente marcada
`\verified` (não `\proved`); o "esboço" é honesto ao dizer que `min(n,m)∈{4,5}`
depende de verificação computacional. Sem sobrevenda aqui. β₁ de todas as 4
linhas confere (21,33,93,177).

**Teorema 4.3 — Suficiência Local (l.925).** ✘ **Marcado `\proved` sem
prova.** Só há a verificação 100% (50/50, 20/20, 10/10). O enunciado é
plausivelmente verdadeiro, mas "prova" não está escrita; e o complemento em
l.970 ("necessariamente 1 ciclo") esconde a hipótese de meta-ciclo único.

**Lean (`Q_eq_three_n4/n5`, l.1360–1361).** ✔ **Não dependem de `sorry`.** São
teoremas `native_decide` fechados pelo kernel; os 4 `sorry` (l.1475–1488) são
todos os enunciados **gerais** (`n` simbólico), não os casos concretos. A
alegação "verificado via Lean para `n∈{4,5}`" (l.167) está **correta e não
contaminada por sorry**. Sem bloqueador aqui — bom.

---

## 4. Tabela de consistência numérica

| quantidade | onde aparece | confere? |
|-----------|--------------|----------|
| `β₁ = |E|−|V|+1` plano (n=6..14: 45,105,189,297,429) | `tab:grafo` l.273–277 | ✔ recalculado, todos batem |
| `β₁` toro (n=4..10: 17,109,193,301) | `tab:toro`, B.2 | ✔ (`|E|=4n²` p/ `n≥6`; `n=4` especial deg=4) |
| `β₁` retangular (21,33,93,177) | `tab:retangulares` | ✔ todos batem |
| `β₁` híbrido k=0..3 (109,103,97,91) | `tab:tightness` | ✔ seguem `144−6k` |
| **`|E|` híbrido k=4 = 80** | `tab:tightness` l.644 | ✘ **deveria ser 120** (4 cantos × −6); 80 é o plano, objeto diferente |
| `N(8)=13.267.364.410.532` | l.202, `tab:knuth` | ✔ valor McKay; log10=13,12 (mas tab:knuth diz 12,9 — ver MENOR) |
| `N_DnC=3,690×10²¹` | abstract l.92, l.206, `tab:ndnc` | ✔ idêntico nas 3 ocorrências; número cheio = 3.690.310.991.081.969.141.272 |
| `N_DnC/8 = 461.288.873.885.246.142.659` | l.1037 | ✔ verificado (resto 0) |
| soma topologias DnC = 3,690×10²¹ | `tab:ndnc` | ✔ (2,414+0,638+0,638) |
| soma configs = 79.968 | `tab:ndnc` | ✔ (76.832+1.568+1.568) |
| `N(12)≈1,3×10³³` | rem. l.1030, `tab:knuth` | ✔ consistente (log 33,1) |
| `r(n)` e `log₂(1/r)` (6..12) | `tab:ratio` | ✔ todos batem; r(6)=9862/36236=0,272 casa com `tab:transfer` |
| espectro `T_bulk(6)`: λ₁=70,48, λ₂=−39,09, gap=31,39, ρ=0,555 | `tab:transfer` | ✔ aritmética interna confere |
| deficit=3 / rank=102 (n=8) | l.664, l.1292, l.1599 | ✔ consistente (β₁=105−3=102) |
| **speedup DnC (8,8/18,9/31,3)** | abstract l.88–90 | ✘ **contradiz `tab:scaling` (0,8/∞/30) e rem. dnc12** |

---

## 5. Lista de ação priorizada (o que consertar primeiro)

1. **Abstract vs corpo no DnC (MAIOR, mais visível).** Os números 8,8/18,9/31,3
   do abstract contradizem frontalmente a Tabela 4.5 e a Observação 4.6. É o
   primeiro lugar que um leitor/revisor checa. Corrigir já.
2. **Calibrar `Q(n)=3` (MAIOR, é a tese central).** Tornar explícito que o
   `\proved` vale para `n∈{6,8,10}` (e GF(2) até 14), e que `n` geral é
   conjectura pendente da conexidade do bulk. O próprio Lean já admite isso
   (sorry #2) — alinhar §2 com §6 para o paper parar de se contradizer.
3. **Resolver a interpolação híbrida → plano (MAIOR, conceitual).** Decidir se
   `G(H)` realmente interpola até o plano; se não, parar de chamar `k=4` de
   plano e recomputar a linha `k=4` no grafo correto (|E|=120) ou marcá-la
   como o plano *separado*, fora da família.
4. **Escrever (ou rebaixar) o Teorema 4.3.** Uma prova de meia página resolve;
   senão, `\verified`.
5. **Corrigir a contradição 1-vs-2 classes na Klein** (l.1074 vs 1107).
6. Menores: log do Knuth `n=8`; "Teorema U para `n=6`" no abstract; "três
   categorias"; descrição de `N(8)`; `\ref{lem:bordas}` no lugar de "2.7";
   soma "26 teoremas".

---

## Nota de compilação

`pdflatex` (2 passadas) → **28 páginas, 0 referências/citações indefinidas, 0
caixas overfull/underfull**, todas as figuras presentes. Tecnicamente limpo.
Os problemas são de **conteúdo e honestidade**, não de LaTeX. O hard-code
"Lema 2.7" funciona hoje (lem:bordas = 2.7) mas é frágil.
