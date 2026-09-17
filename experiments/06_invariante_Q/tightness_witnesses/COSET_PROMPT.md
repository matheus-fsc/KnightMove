# Prompt para o agente — validação da reformulação por *coset*

## Contexto

Repositório `/home/math/Dev/knight_tour`. Material relevante:

- `tightness_witnesses/` — `knight_gf2.py` (infra GF(2), enumerador sem viés),
  `bulk_perp_dims.py`, `flip_graph_rect.py` (enumeração exaustiva de tours em
  retângulos), `FLIP_GRAPH_RESULTS.md`.
- `forbidden_cycles_6x6/data/results/forbidden_cycles.json` — os 3 ciclos
  proibidos com rótulos de aresta.
- `paper/knight_tour_tightness.tex` — Paper A, §2 (Q(n)=3) e §3 (tightness).
- `knight_tour_lean/KnightTour/` — `QAbstract.lean`, `Reduction.lean`,
  `PuncGraph.lean`, `CycleSpaceDim.lean` (todos sem `sorry`).

## A reformulação a validar

Notação: $G$ grafo do cavalo $n\times m$; $Z_1 = \ker \partial_1$ (espaço de
ciclos); $\beta_1 = |E| - |V| + 1$; cantos $c_1..c_4$, todos de grau 2, com
arestas obrigatórias $e_1(c_i), e_2(c_i)$.

**(A) Os funcionais de canto.** Para $z \in Z_1$, o grau em $c_i$ é par e
$\deg_G(c_i)=2$, logo $z$ usa **as duas** arestas de $c_i$ ou **nenhuma**.
Define-se $y_i : Z_1 \to \mathbb{F}_2$, $y_i(z) = [z$ usa $c_i]$.

**(B) Os tours vivem numa fibra.** Todo tour satisfaz $y = (1,1,1,1)$.
Logo $\Ham \subseteq A := \{z \in Z_1 : y(z) = \mathbf{1}\}$, e
$Z_{\mathrm{bulk}} = \{z \in Z_1 : y(z) = \mathbf{0}\}$.

**(C) Cota superior.** Se os 4 funcionais são independentes em $Z_1$, então
$\dim A = \dim Z_{\mathrm{bulk}} = \beta_1 - 4$, e como $A$ é afim sem conter
$0$, $\dim \Span(A) \le \beta_1 - 3$. Daí $\rk(\Ham) \le \beta_1 - 3$.

**(D) O quociente.** *Sob tightness*, $\Span(\Ham) = \langle\tau\rangle \oplus
Z_{\mathrm{bulk}} = \{z : y(z) \in \{\mathbf{0},\mathbf{1}\}\}$, e portanto
$Z_1/\Span(\Ham) \cong \mathbb{F}_2^4/\{\mathbf{0},\mathbf{1}\}$ — a classe de
um ciclo depende **só** do seu padrão de cantos, módulo complementação.

## Duas ressalvas que quem escreveu isto já identificou — confirme ou refute

1. **(C) NÃO é mais barata que a prova existente.** A independência dos 4
   funcionais equivale a $c(\mathrm{Punc}) = 5$, isto é, **conexidade do
   bulk** — o mesmo ingrediente do `Q(n)=3` atual. A reformulação ganha em
   exposição, não em força. *Verifique que a equivalência é mesmo essa.*
2. **(D) é consequência da tightness, não derivação independente.** Só vale
   onde $\span\{\tau - \tau'\} = Z_{\mathrm{bulk}}$. *Verifique que (D) falha
   exatamente onde a tightness falha.*

---

## Testes

Use aritmética exata sobre $\mathbb{F}_2$. Não amostre onde houver enumeração
exaustiva disponível.

### T1 — $y_i$ bem-definido ⟨rápido⟩
Sobre uma base de $Z_1$, confirme que todo elemento tem 0 ou 2 arestas em cada
canto. Tabuleiros: $6\times6$, $8\times8$, $6\times7$, $5\times8$, $5\times6$.
Predição: sempre. Se falhar, (A) está errado.

### T2 — posto dos 4 funcionais ⟨rápido⟩
Calcule $r = \rk\{y_1,\dots,y_4\}$ como funcionais em $Z_1$, e
$\dim\{z : y(z)=\mathbf 0\}$. Predição: $r = 4$ e
$\dim = \beta_1 - 4$, **sse** o bulk é conexo.

⚠️ **Inclua um tabuleiro de bulk desconexo** se conseguir construir um (tente
$m\times n$ pequenos, ou remova arestas artificialmente). Predição: $r < 4$ e
$\dim > \beta_1 - 4$. Sem esse controle negativo, T2 não distingue "vale
sempre" de "vale porque testei só casos bons".

### T3 — cota superior ⟨rápido⟩
Verifique $\rk(\Ham) \le \beta_1 - 3$ e compare com $Q$ calculado pela
definição atual do paper. Predição: coincidem em todos os casos com $r=4$.

### T4 — tours na fibra ⟨exaustivo em 6×6, 5×8, 6×7⟩
Confirme $y(\tau) = \mathbf{1}$ para **todo** tour. Use os conjuntos
exaustivos (`flip_graph_rect.py` já os produz). Predição: sempre.

### T5 — a equivalência ⟨exaustivo⟩
Verifique que $\rk(\Ham) = \beta_1 - 3 \iff \span\{\tau - \tau'\} =
Z_{\mathrm{bulk}}$, computando os dois lados independentemente em $6\times6$
(vale) e $5\times8$ (**falha**: deficit 9).

### T6 — a classificação (D) ⟨o teste central⟩
Em $6\times6$: compute $Z_1/\Span(\Ham)$ e o mapa induzido por $y$. Verifique
que é isomorfismo sobre $\mathbb{F}_2^4/\{\mathbf 0,\mathbf 1\}$, e que os 3
ciclos proibidos têm padrões $(1000), (1100), (1001)$ — independentes módulo
$\mathbf 1$.

⚠️ **Controle negativo obrigatório:** repita em $5\times8$. Lá o deficit é 9,
então $\dim Z_1/\Span(\Ham) = 9 \neq 3$ e o mapa $y$ **não pode** ser
isomorfismo. Identifique o que sobra: quais são as 6 dimensões extras, e elas
têm padrão de canto trivial? (Predição: sim — vivem dentro de
$Z_{\mathrm{bulk}}$, invisíveis a $y$.)

### T7 — $Q(5,8)$ diretamente ⟨decisivo⟩
Compute $Q(5,8)$ pela definição do paper e confirme conexidade do bulk de
$5\times8$. Predição: bulk conexo e $Q = 3$.

Se ambos, então $5\times8$ é um contraexemplo exato mostrando

> **bulk conexo + $Q=3$ NÃO implicam tightness**,

o que é resultado publicável por si e delimita o que a hipótese
$\min(n,m)\ge 6$ está realmente fazendo.

---

## Interpretação — fixe antes de rodar

| resultado | leitura |
|---|---|
| T1–T5 confirmam, T6 confirma com o controle negativo | reformulação válida; vale reescrever §2 do paper e considerar Lean |
| T7 dá bulk conexo e $Q(5,8)=3$ | contraexemplo forte; **reportar em destaque** |
| algum teste falha | a reformulação está errada — **reporte e pare**, não contorne |
| T2 sem controle negativo disponível | reporte como não-testado, não como confirmado |

---

## Armadilhas deste projeto — respeite

1. **⚠️ Enumeração exaustiva onde existir.** $6\times6$ (9.862), $5\times8$
   (44.202), $6\times7$ (1.067.638) são exatos e batem com a literatura.
   Ranks amostrais já produziram um falso deficit aqui (rank 96 em $n=8$ por
   viés do Warnsdorff).
2. **⚠️ Sature todo corte de busca** e exiba a saturação. Um corte no
   predicado medido já invalidou uma sessão inteira (`max_w=20`).
3. **⚠️ Distinga "não encontrado no orçamento" de "provado inexistente".**
4. **⚠️ Controles negativos não são opcionais.** Cada predição positiva deste
   prompt tem um caso onde deve falhar; sem ele o teste não discrimina.

---

## Parte Lean — **condicional**

Só inicie se T1–T6 confirmarem.

**O que formalizar:**
- `cornerFunctional (i) : Z₁ →ₗ[F2] F2`, bem-definido via grau 2 do canto
- `tour_corner_functional_eq_one` — todo tour está na fibra $\mathbf 1$
- `finrank_span_affine` — span linear de conjunto afim sem 0 tem dim = dim + 1
- `rank_Ham_le` — a cota $\rk(\Ham) \le \beta_1 - 3$ por essa rota

**Seja honesto sobre o valor:** o projeto **já tem** `Q_abstract_eq_three` sem
`sorry`. Esta rota não acrescenta conteúdo matemático; o ganho é (i) prova mais
curta e (ii) evitar o quociente $\mathbb{F}_2^E/\mathrm{row}(\partial_1)$, cuja
exposição em prosa foi criticada por confusão com $H_1/\mathrm{row}(\partial_1)$
(que é mal definido, pois o espaço de cortes não é subespaço do de ciclos).

**Não** escreva `theorem ... := sorry`. Se um passo não fechar, documente como
comentário com a assinatura pretendida — é a convenção do projeto e existe
porque `sorry` no namespace contamina dependências em silêncio.

Reaproveite `PuncGraph.lean` e `CycleSpaceDim.lean`; a fibra $y=\mathbf 0$ deve
coincidir com `ZBulk` já definido (`ZBulk_eq`).

---

## Entregáveis

1. `tightness_witnesses/COSET_RESULTS.md` — tabela de interpretação preenchida,
   veredicto explícito por teste, e os controles negativos em destaque.
2. Código reproduzível, comandos no topo do RESULTS.
3. Dados em `tightness_witnesses/data/`.
4. Se a parte Lean rodar: módulo novo, com `#print axioms` de cada declaração.
5. Se algo contradisser as ressalvas 1–2 da seção inicial: **reporte, não
   acomode** — elas foram escritas como suspeitas, não como fatos.

## Padrão de honestidade

- Reporte negativo com o mesmo destaque do positivo.
- Distinga sempre "verificado para $n \in \{...\}$" de "provado para todo $n$".
- Esta validação testa uma **reformulação**, não um resultado novo. Se ela se
  confirmar, nada fica provado que já não estivesse — apenas fica mais legível.
  Não anuncie como avanço matemático.
