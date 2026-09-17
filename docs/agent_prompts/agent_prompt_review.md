# Agent Prompt — Revisão Crítica do Artigo (knight_tour_complete.tex / .pdf)

## CONTEXTO

Projeto: `knight_tour/` — estudo de invariantes topológicos do passeio do
cavalo via espaço de ciclos sobre $\GF(2)$.

**Artigo a revisar:** `knight_tour_complete.tex` (1704 linhas, pt-BR,
`pdflatex knight_tour_complete.tex`). O PDF compilado é
`knight_tour_complete.pdf`. Autor: Matheus Coelho.

O artigo usa três marcadores de status, definidos como macros:
- `\proved` → **[Provado]** (prova matemática completa)
- `\verified` → **[Verificado]** (evidência computacional, não prova)
- `\openp` → **[Aberto]** (conjectura / trabalho futuro)

**Tese central do paper:** o invariante algébrico $Q(G)$ (baseado em
vértices de grau 2 do grafo do cavalo) satisfaz $Q(n)=3$ para todo
$n\geq 6$; há um "Teorema U" para a família híbrida toro→plano; e o
método Divide-and-Conquer por blocos $6\times6$ dá hamiltonicidade global
a partir de compatibilidade local.

---

## REGRA ABSOLUTA

**Você é um revisor cético, não um co-autor.** Seu trabalho é encontrar
erros, lacunas e exageros — não elogiar. Cada afirmação marcada `\proved`
deve ser tratada como uma alegação de prova que você precisa auditar linha
a linha. Não confie no marcador de status: verifique se o que está escrito
realmente corresponde ao nível de rigor reivindicado.

**Não invente correções.** Se você suspeita de um erro mas não tem certeza,
classifique como "dúvida a investigar", não como "erro confirmado".

---

## TAREFA

Produza um relatório de revisão em `review_report.md` cobrindo os eixos
abaixo. Para cada achado, cite **linha(s) do .tex** e classifique a
severidade: `[BLOQUEADOR]` / `[MAIOR]` / `[MENOR]` / `[DÚVIDA]` / `[ESTILO]`.

### 1. Correção matemática das provas `\proved`
Audite cada `theorem`/`lemma`/`corollary`/`proposition` marcado `\proved`:
- Teorema $Q(n)=3$ (l.320) e seus lemas de apoio (graus dos cantos l.330,
  relação aresta–canto l.342, independência dos 4 representantes l.366,
  imagem dos pares XOR l.411). A prova de $Q(n)=3$ depende de "4
  representantes independentes + kernel da soma" — verifique se a
  independência é realmente estabelecida e se o passo "$\geq 3$ e $\leq 3$"
  fecha sem buracos.
- Teorema U (l.494) e seus corolários (monotonia, invariância geométrica,
  plano/toro como extremos).
- Teorema de Suficiência Local do DnC (l.925).
- Para cada prova: o enunciado casa com o que é provado? Há quantificadores
  implícitos ("para todo canto", "existe ciclo") que ficam não-justificados?
  Casos de base e indução cobrem o que dizem cobrir?

### 2. Conjectura vs. teorema (calibração de status)
- Algo marcado `\proved` que na verdade é só `\verified` (apenas evidência
  computacional)? Ex.: $Q(n,m)=3$ retangular está como `\verified` (l.578) —
  confira se o texto não escorrega para linguagem de "provado".
- A *tightness* (l.623) e o deficit toroidal (l.667) estão como conjectura/
  evidência — confira se não há frase no corpo que os trate como teorema.
- Lean 4 (Seção l.1330): o paper admite 4 `sorry`'s (l.1470). Verifique se
  as alegações de "verificado via Lean" para $n\in\{4,5\}$ (l.167) NÃO
  dependem de um `sorry`. Se dependerem, é `[BLOQUEADOR]`.

### 3. Consistência numérica entre seções
O paper é denso em números; eles precisam bater em todo lugar:
- Tabela de $\beta_1$ (l.273–276): $n=6\to45$, $n=8\to105$, $n=10\to189$,
  $n=12\to297$. Recalcule $\beta_1=|E|-|V|+1$ a partir dos $|V|,|E|$ dados
  e confirme.
- $N(6)=9862$, $N(8)=13\,267\,364\,410\,532$ (McKay, l.219),
  $N_{\mathrm{DnC}}(12\times12)=3{,}69\times10^{21}$ (l.224 e Seção l.1002):
  o mesmo número aparece igual nas duas ocorrências?
- Estimativas de Knuth (Seção l.1116) vs. números citados na introdução.
- Qualquer número que apareça duas vezes (deficit=3, rank, "4–6 nós/tour",
  contagens) deve ser idêntico nas duas. Liste discrepâncias.

### 4. Definições e notação
- $Q(G)$ tem definição única e usada de forma consistente? O paper alerta
  (l.176+) que existe um "$Q$ topológico" distinto para Klein/superfícies —
  verifique se toda ocorrência de $Q$ deixa claro qual é, e se a Observação
  `rem:klein-q` (Seção l.1043) é coerente com a nota da introdução.
- Símbolos usados antes de definir? Macros de notação (l.~40+) que somem?
- `\ref`/`\label` quebrados (saída "??" no PDF), `\cite` sem entrada na
  bibliografia, figuras/tabelas referenciadas mas ausentes.

### 5. Trabalhos relacionados e citações
- As atribuições estão corretas? (Schwenk 1991 — caracterização de
  existência; McKay 1997 — contagem $N(8)$; Warnsdorff 1823; Diestel 2017
  cap. 1.9; Knuth — estimador.) Confira se o que é atribuído a cada um
  corresponde ao resultado real.
- O paper reivindica originalidade ("primeira formalização em Lean 4 do
  espaço de ciclos $\GF(2)$ do cavalo", l.~237). Sinalize reivindicações de
  primazia como `[DÚVIDA]` para o autor checar — não as confirme.

### 6. Algoritmos (Seção l.730)
- O pseudocódigo do backtracking R2+UF (l.733) e do DnC (l.919) está
  correto e corresponde ao que o texto descreve? Invariantes de loop
  declarados realmente se mantêm?
- Alegações de complexidade ("$\sim$4–6 nós/tour", *scaling* l.973) são
  empíricas ou provadas? Estão rotuladas como tal?

### 7. LaTeX / compilação
- Compile: `pdflatex -interaction=nonstopmode knight_tour_complete.tex`
  (rode 2×). Liste **warnings** relevantes: referências indefinidas,
  citações indefinidas, overfull \hbox graves, fontes faltando.
- NÃO conserte o .tex; apenas reporte. (Se o usuário pedir depois, aí sim.)

### 8. Coerência global e honestidade
- O Sumário de Resultados (l.163) promete exatamente o que o corpo entrega?
  Algo prometido como "provado" no sumário mas só "verificado" no corpo?
- A introdução vende algo que a conclusão não sustenta?
- Tom: o paper se mantém honesto sobre provado-vs-conjecturado? Aponte
  qualquer lugar onde a linguagem exagera.

---

## ENTREGÁVEL

Arquivo `review_report.md` com:

1. **Veredito em uma linha** (ex.: "Sólido, com 2 maiores e 6 menores" /
   "Bloqueador em §X").
2. **Tabela de achados**: `severidade | linha | seção | descrição | sugestão`.
3. **Auditoria prova-a-prova**: um parágrafo por resultado `\proved`,
   dizendo se a prova fecha, e onde aperta.
4. **Tabela de consistência numérica**: cada número-chave, onde aparece, se
   bate.
5. **Lista de ação priorizada** para o autor (o que consertar primeiro).

Comece lendo o `.tex` inteiro **antes** de escrever qualquer achado.
Não pule a auditoria das provas — é a parte mais importante.
