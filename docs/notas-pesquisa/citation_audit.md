# Auditoria de citações — `knight_tour_complete.tex`

**Data:** 2026-08-06
**Escopo:** bibliografia atual (10 `\bibitem`), seção *Trabalhos relacionados* (linhas 235–311), e todas as reivindicações de novidade espalhadas pelo texto.
**Fontes da auditoria:** OpenAlex API, busca web, OEIS A001230.
**Status:** apenas diagnóstico — nenhuma edição feita no `.tex`.

---

## TL;DR — o veredicto sobre "redescobri a roda"

Separando em três camadas, porque o resultado é **bem melhor do que você temia**:

| Camada | Existe literatura? | Impacto |
|---|---|---|
| **GF(2) / espaço de ciclos aplicado ao grafo do cavalo** (o núcleo: `Q(n)=3`, deficit, ciclos proibidos) | **Não encontrei nada.** Buscas dirigidas em OpenAlex e web não retornam nenhum trabalho que aplique espaço de ciclos $\mathbb{F}_2$ ao grafo do cavalo. | ✅ **A novidade específica sobrevive.** |
| **O *framework* geral "quando ciclos hamiltonianos geram o espaço de ciclos"** | **SIM, e é uma área ativa com nome próprio.** Chama-se ***Hamilton space* / grafos *Hamilton-generated***. Heinig, Hefetz–Krivelevich, Hou–Yin, Christoph–Nenadov–Petrova, 2024–2026. | 🔴 **Crítico.** Seu `Ham(n) ⊆ Z₁` e seu `deficit` são *exatamente* $\mathcal{C}_n(G) \subseteq \mathcal{C}(G)$ e a codimensão dele. Não citar isso é o único item que pode custar o paper numa revisão. |
| **Topologia de tours do cavalo em superfícies (toro, cilindro, Klein, Möbius)** | **SIM, publicado e recente.** Watkins (livro, 2004), Forrest–Teehan (2015), Forrest–Lague (2024). Classificam tours por classe de homotopia em $\pi_1$. | 🟠 **Aqui você redescobriu a roda de fato.** Sua Seção 6 (hierarquia de $Q$ por superfície, `Q_toro=0`, `Q_Klein`) refaz — com outra linguagem — o que esses papers fazem. |

**Leitura prática:** você *não* perdeu tempo no resultado principal. Perdeu na seção de superfícies, e está exposto por não conhecer a literatura de *Hamilton space* — que, ironicamente, é a que **valida** seu enquadramento e te dá a ferramenta que falta para a *tightness* (ver §1.3 abaixo).

**Duas ressalvas descobertas na 2ª rodada** (leitura dos PDFs originais), ambas com consequência concreta:

- 🔴 **`deficit > 0` é prior art.** Heinig (2013) já observou que $\delta(G) \geq 3$ é *necessário* para Hamilton-generation. O cavalo tem cantos de grau 2. Você precisa dizer isso e reposicionar a contribuição no **valor exato 3** e na **constância em $n$** — ver §1.2 reescrito.
- 🔴 **A lei $N(n) \sim 1{,}82^{n^2}$ não se sustenta nos seus próprios dados.** É artefato de OLS com resíduos de fator ~7 e constantes locais decrescentes — ver §2.4. Não é erro de citação, é erro estatístico, mas é o item mais fácil de um parecerista atacar.

---

## 1. Lacunas CRÍTICAS (bloqueiam publicação)

### 1.1 A literatura de *Hamilton space* — seu enquadramento inteiro já tem nome

Você define (linhas 166–171 e Def. de `deficit`):

$$\Ham(n) \subseteq \text{2-fatores}(n) \subseteq \ker(\partial_1) = Z_1(G_n;\mathbb{F}_2), \qquad \mathrm{deficit}(n) := \beta_1(n) - \rank_{\mathbb{F}_2}(\Ham(n))$$

Na literatura isto é: $\mathcal{C}_n(G)$ = **Hamilton space** (subespaço de $\mathcal{C}(G)$ gerado pelos vetores de incidência dos ciclos hamiltonianos); um grafo é **Hamilton-generated** quando $\mathcal{C}_n(G) = \mathcal{C}(G)$; e seu `deficit` é a **codimensão de $\mathcal{C}_n(G)$ em $\mathcal{C}(G)$**.

Referências a adicionar (ordem de importância):

| Ref | Trabalho | Por quê |
|---|---|---|
| `Heinig2014` | S. Heinig, resultados iniciais: $\delta(G) \ge n/2 + 37$ implica $\mathcal{C}_n(G)=\mathcal{C}(G)$ | Origem da pergunta. |
| `HouYin2025` | Hou & Yin, *Dirac-type condition for Hamilton-generated graphs*, arXiv:2503.15950 | Resolve completamente o caso Dirac ($\delta > n/2$). Introduz o termo **Hamilton-generated**. |
| `HefetzKrivelevich2025` | Hefetz & Krivelevich, *The Hamilton cycle space of random graphs*, arXiv:2506.19731; e *...random regular and randomly perturbed graphs*, arXiv:2507.04488 | Estado da arte 2025. |
| `ChristophNenadovPetrova2024` | M. Christoph, R. Nenadov, K. Petrova, *The Hamilton space of pseudorandom graphs*, arXiv:2402.01447 | Introduz **parity-switchers** — ver §1.3. |
| `HamGenerated2026` | *On graphs whose cycle space is spanned by their Hamilton cycles*, arXiv:2606.05835 | Survey de fato da área na introdução; fortalece Chvátal–Erdős e McDiarmid–Yolov. |

**Onde inserir:** novo parágrafo em *Trabalhos relacionados* (§1.5), imediatamente após o parágrafo `\paragraph{$\GF(2)$ e homologia.}` (linha 264). Sugestão de conteúdo:

> A pergunta "quando os ciclos hamiltonianos geram todo o espaço de ciclos?" tem literatura própria sob o nome *Hamilton space* / grafos *Hamilton-generated* [Heinig, Hou–Yin, Hefetz–Krivelevich, Christoph–Nenadov–Petrova]. Esses trabalhos operam no regime **denso/pseudoaleatório** (condições tipo Dirac, $\delta(G) \gtrsim n/2$), onde a resposta típica é $\mathcal{C}_n(G) = \mathcal{C}(G)$ — deficit zero. O grafo do cavalo está no extremo oposto: é **esparso** ($|E_n| \approx 4|V_n|$), tem grau máximo 8 e grau mínimo 2. Nossa contribuição é exibir uma família esparsa natural com **deficit positivo, exato e constante** ($=3$), determinado por geometria local (os quatro cantos), em contraste com o regime denso onde o deficit desaparece.

**Isto é bom para você, não ruim.** Reposiciona `Q(n)=3` de "invariante ad hoc" para "primeiro exemplo estruturado de deficit constante numa família esparsa" — que é uma contribuição *mais* forte e claramente complementar à literatura existente.

### 1.2 A obstrução de grau 2 É prior art documentada — ✅ RESOLVIDO

**Atualização (verificado no texto completo de arXiv:2606.05835).** Minha suposição inicial de que isso seria apenas "folclore não escrito" estava errada. Está publicado, explicitamente:

> "if $p$ is large enough to ensure that $\delta(G) \geq 3$ holds a.a.s., then it also ensures that a.a.s. $\mathcal{C}_n(G) = \mathcal{C}(G)$ **(the necessity of this condition was observed by Heinig in [22])**."
> — Hefetz & Krivelevich, arXiv:2606.05835, §1

Ou seja: $\delta(G) \geq 3$ é **necessário** para Hamilton-generation, e a necessidade é atribuída a **Heinig, *When Hamilton circuits generate the cycle space of a random graph*, arXiv:1303.0026, 2013.**

O argumento é o elementar: se $\deg(v)=2$ com arestas $e_1,e_2$, todo ciclo hamiltoniano contém ambas, logo todo $\mathbf{1}_\tau$ satisfaz $x_{e_1}=x_{e_2}$ — um hiperplano próprio contendo $\mathcal{C}_n(G)$.

**Consequência direta para o seu paper:** o grafo do cavalo tem $\delta(G_n) = 2$ (os quatro cantos). Logo, **pela observação de Heinig, $G_n$ não é Hamilton-generated — isto é, $\mathrm{deficit}(n) > 0$ é prior art.**

🔴 **Isto precisa ser dito explicitamente no paper.** Hoje o texto apresenta a existência da obstrução como parte da descoberta. Não é.

**O que continua sendo seu (e é o conteúdo real):**

1. O **valor exato** $\mathrm{deficit} = 3$, não apenas $> 0$. O colapso $4 \to 3$ (quatro cantos, três dimensões, via o núcleo do funcional soma $\sigma$) é onde está a matemática.
2. A **constância em $n$** — o deficit não cresce com o tabuleiro.
3. A **conexidade do bulk**, que é o ingrediente geométrico não-trivial e o que você formalizou em Lean.
4. A **tightness** ($\rank(\Ham) = \beta_1 - 3$ exatamente), ainda conjectural, que diz que *não há outras obstruções*.

**Ação (reescrita sugerida da observação no §2.2):**

> Observação. Que $G_n$ não seja *Hamilton-generated* segue de um fato conhecido: $\delta(G) \geq 3$ é necessário para $\mathcal{C}_n(G) = \mathcal{C}(G)$ [Heinig 2013], e os quatro cantos de $G_n$ têm grau 2. O conteúdo do Teorema~\ref{thm:Qn3} não é, portanto, a *existência* da obstrução, mas (i) o **valor exato** da dimensão do espaço gerado pelas quatro obstruções simultâneas — que colapsa de 4 para 3 pelo núcleo do funcional soma — e (ii) sua **constância em $n$**. A tightness (Conjectura~\ref{conj:deficit}) afirma adicionalmente que estas são as *únicas* obstruções.

Isso te protege de um parecerista dizendo "a obstrução é trivial e conhecida" — porque ele estaria certo sobre a existência, e você precisa ter dito isso primeiro.

**Referência adicional relacionada, também da bibliografia de [2606.05835]:**
> I. B.-A. Hartman, *Long cycles generate the cycle space of a graph*, European J. Combin. **4** (1983), 237–246. — Ciclos de comprimento $\geq d+1$ geram $\mathcal{C}(G)$ para $G$ 2-conexo com grau mínimo $d$. Ancestral direto de toda a linha.

### 1.3 *Parity-switchers* — a ferramenta que falta para a `tightness`

Sua §12.1 (*Tightness generalizada*, linha 1695) descreve como estratégia plausível um "lema de deformação algébrica via *twin-tour expansion*". A literatura de Hamilton space já construiu exatamente esse objeto: Christoph–Nenadov–Petrova (arXiv:2402.01447) usam **parity-switchers**, descritos como "análogos de *absorbers* no contexto de espaços de ciclos" — gadgets locais que permitem trocar a paridade de uma coordenada mantendo hamiltonicidade.

**Ação:** citar em §12.1 e reescrever a conjectura de *tightness* em termos de "existência de parity-switchers suficientes no grafo do cavalo". Isso transforma uma conjectura vaga numa pergunta técnica bem-posta e conectada — muito mais forte num paper.

### 1.4 Superfícies: você foi escooped (parcialmente)

Sua §6 (*Superfícies alternativas*, hierarquia de $Q$: plano 3, cilindro 0, toro 0; garrafa de Klein) não cita nada. Existe:

| Ref | Trabalho |
|---|---|
| `Watkins2004` | J. J. Watkins, *Across the Board: The Mathematics of Chessboard Problems*, Princeton UP, 2004. Capítulos dedicados a "The Torus and the Cylinder" e "The Klein Bottle and Other Variations". Watkins & Hoenigman caracterizaram os toros que admitem tour; Watkins fez o análogo para cilindros, Möbius e Klein. |
| `ForrestTeehan2015` | B. Forrest, K. Teehan, *The Topology of Knight's Tours on Surfaces*, arXiv:1507.02917. Caracterizam dimensões de tabuleiros cilíndricos e toroidais que admitem tours fechados realizando **a identidade de $\pi_1$** versus **um gerador de $\pi_1$**. |
| `ForrestLague2024` | B. Forrest, Lague, *Nullhomotopic and Generating Knight's Tours on Non-Orientable Surfaces*, arXiv:2406.05226. Mesmo programa para Möbius e garrafa de Klein. |

**Sobreposição concreta.** Sua memória de projeto registra:
- `Q_toro = 0` "4 classes (par_y, par_x) ~25% cada" → isso é a afirmação de que tours toroidais realizam **todas** as classes de $H_1(T^2;\mathbb{F}_2)$, i.e. existem tours nullhomotópicos *e* geradores. Forrest–Teehan.
- Klein: "off-diag (0,1),(1,0) SEMPRE proibido (H₁=ℤ não ℤ²)" → é a estrutura de $\pi_1$(Klein) sendo detectada empiricamente. Forrest–Lague.

**Atualização — li o texto completo de Forrest–Teehan (arXiv:1507.02917, 24 pp.).** A sobreposição é real mas **mais estreita do que eu temia**, e há uma diferenciação legítima e defensável.

O que eles provam exatamente (Teoremas 4.1, 5.1, 4.2, 6.1): para cilindros $C_{m,n}$ e toros $T_{m,n}$, caracterizam **para quais $(m,n)$ existe** um tour realizando (a) a identidade de $\pi_1$, e (b) um gerador de $\pi_1$ (no toro: a classe de homotopia de uma longitude). Exemplos: Teorema 4.2 — "quando $m+n>3$, $T_{m,n}$ suporta um tour nullhomotópico"; Teorema 6.1 — "quando $m,n$ não são ambos 1, $T_{m,n}$ suporta um tour que [realiza a longitude]".

**A diferenciação que te salva:** eles respondem uma pergunta de **existência** ("existe tour na classe $X$?"). Você responde uma pergunta de **distribuição** ("que fração dos tours cai em cada classe?"). Sua memória registra `Q_toro=0` como "4 classes $(\text{par}_y,\text{par}_x)$ com $\approx 25\%$ cada em 20k tours $6\times6$" e `Q_cyl=0` como "paridade do winding 50/50 em 18k tours únicos". Isso é **equidistribuição**, que é estritamente mais forte que existência — e não está nos papers deles.

**Ação recomendada (revisada):** manter a seção, encolher para ~1 página, e enquadrar assim:

> Forrest e Teehan [FT2015] caracterizaram, para cilindros e toros, quais dimensões admitem tours realizando a identidade ou um gerador de $\pi_1$; Forrest e Lague [FL2024] estenderam para Möbius e Klein. Nossa contribuição aqui é **complementar e de natureza distinta**: onde eles estabelecem *existência* por classe de homotopia, nós medimos a *distribuição* dos tours entre as classes. O achado é de **equidistribuição** — no toro $6\times6$, as quatro classes de $H_1(T^2;\mathbb{F}_2)$ recebem $\approx 25\%$ dos tours cada; no cilindro, as duas classes recebem $\approx 50\%$ cada. Em termos do nosso invariante algébrico, isso é exatamente $Q = 0$: nenhuma restrição linear separa as classes. O contraste com $Q = 3$ no plano isola a obstrução dos cantos como fenômeno **de fronteira**, ausente quando a fronteira desaparece.

Esse enquadramento é honesto, cita o prior art, e ainda assim deixa a seção com conteúdo próprio. Melhor que cortar.

**Nota adicional:** sua *"Nota sobre notação"* (linhas 224–233) já sinaliza que existe um "$Q$ topológico" distinto para a garrafa de Klein baseado em winding number. **É exatamente o invariante desses papers.** Cite ali — é o lugar natural.

**Referências extras colhidas da bibliografia do Forrest–Teehan** (úteis se você mantiver a seção):
- **J. J. Watkins, R. Hoenigman**, *Knight's Tour on a Torus*, Mathematics Magazine **70**(3) (1997), 175–184. — A origem do estudo de tours toroidais.
- **J. J. Watkins**, *Across the Board*, Princeton UP, 2004, **pp. 65–77** (páginas exatas do capítulo sobre cilindros/toros). Teorema 2.2 do FT é citado como "Watkins [12]; pg 71": $C_{m,n}$ admite ciclo hamiltoniano exceto se $m=1, n>1$, ou $m \in \{2,4\}$ e $n$ par.
- **Cannon, Dolan**, *The Knight's Tour*, The Mathematical Gazette **70** (1986), 91–100.
- **Ralston** — boards com $m,n$ ambos ímpares são "odd-tourable".
- **Miller, Farnsworth** — tours fechados em cilindros e toros com um quadrado removido.

---

## 2. Lacunas IMPORTANTES (credibilidade)

### 2.1 Contagem: falta a metade da história

`\bibitem{McKay1997}` está sozinho. A história padrão exige o par:

- **Löbbing & Wegener (1996)**, *The Number of Knight's Tours Equals 33,439,123,484,294 — Counting with Binary Decision Diagrams*, Electron. J. Combin. **3**(1), R5. DOI `10.37236/1229`. **O número no título está errado** (erro de fator/off-by-one no BDD).
- **McKay (1997)** corrigiu para $N(8) = 13\,267\,364\,410\,532$ não-direcionados ($= 26\,534\,728\,821\,064$ direcionados).
- **Wegener (2000)**, *Branching Programs and Binary Decision Diagrams*, SIAM, p. 369 — onde o valor corrigido é registrado.
- **OEIS A001230** — "Number of undirected closed knight's tours on a 2n × 2n chessboard": `0, 0, 9862, 13267364410532`. **Cite isto para $N(6)=9862$**, que hoje aparece no seu paper sem fonte (§5.1 e Tabela de validação).

**Ação:** reescrever `\paragraph{Contagem.}` (linha 246) narrando Löbbing–Wegener → McKay, e citar OEIS A001230 na primeira aparição de 9862 (§3, validação do DP C++, e §5.1).

### 2.2 Estimador de Knuth: referência errada

Você cita `\bibitem{Knuth2022}` = TAOCP Vol. 4B (2022). A referência primária é:

> D. E. Knuth, **"Estimating the efficiency of backtrack programs"**, *Mathematics of Computation* **29**(129), 121–136, 1975. DOI `10.1090/S0025-5718-1975-0373371-6`. (255 citações)

**Ação:** citar Knuth 1975 como primária; manter TAOCP 4B como secundária se quiser.

Complementos relevantes ao seu §5.1:
- **P. C. Chen (1992)**, *Heuristic Sampling: A Method for Predicting the Performance of Tree Searching Programs*, SIAM J. Comput. **21**(2). DOI `10.1137/0221022` — estratificação que reduz a variância; **diretamente aplicável ao seu CV≈150 em $n=10$**, que hoje você apenas reporta como limitação.

### 2.3 Você não é o primeiro a estimar $N(n)$ do cavalo por amostragem

🟠 **Prior art direto na sua §5.**

- **H. Cancela, E. Mordecki**, *Counting Knight's Tours through the Randomized Warnsdorff Rule*, arXiv:math/0609009, 2006. Usam **importance sampling** com Warnsdorff randomizado num esquema de backtracking. Reportam $\approx 1{,}22 \times 10^{15}$ tours abertos geometricamente distintos.
- **H. Cancela, E. Mordecki**, *On the number of open knight's tours*, arXiv:1507.03642, 2015.
- **P. Hingston, G. Kendall**, *Enumerating Knight's Tours using an Ant Colony Algorithm*, CEC 2005. DOI `10.1109/CEC.2005.1554800`.

**Ação:** citar em §5.1 e diferenciar explicitamente: o estimador deles é *importance sampling calibrado por Warnsdorff*; o seu é *Knuth puro sobre a árvore R2+UF*, sem calibração prévia — e cobre tours **fechados**, $n \in \{6,8,10,12\}$. A diferenciação é fácil de fazer e te protege.

### 2.4 Bounds provados para $N(n)$ — ✅ RESOLVIDO, e a lei $1{,}82^{n^2}$ não se sustenta

Extraí o texto completo de Kyek–Parberry–Wegener 1997. Os expoentes são em $\boldsymbol{c^{n^2}}$ (confirmado), e os resultados são:

| Resultado | Enunciado |
|---|---|
| Teorema 4 | $\mathcal{Y}_n = \Omega(1{,}3535^{n^2})$ — **melhor lower bound provado** para tours no $n \times n$ |
| Teorema 3 | tours *estruturados*: $\Omega(1{,}2862^{n^2})$ |
| Teorema 2 | para todo $n \geq 12$: $\mathcal{Y}_n > 1{,}1646^{n^2}$ (versão construtiva/efetiva) |
| Observação 6 | upper bound geral: $\mathcal{Y}_n \leq 4^{n^2}$ (o grafo tem "um pouco menos que $4n^2$" arestas) |
| Teorema 7 | $n=8$: $\mathcal{Y}_8 \leq 3{,}019 \times 10^{22}$ |
| — | $n=8$, melhor lower bound anterior: $122\,802\,512$, devido a **Kraitchik** |

Portanto o intervalo provado para a constante de crescimento é $\boldsymbol{[1{,}3535,\; 4]}$.

**Boa notícia:** $1{,}82 \in [1{,}3535,\ 4]$. Sua estimativa é consistente com o que está provado, e é **muito mais precisa** — você estreita um intervalo de largura 2,6 para um ponto. Isso é vendável e você deve dizer.

**Má notícia — e aqui está o problema real.** Refiz o ajuste com seus quatro pontos:

| $n$ | $n^2$ | $\log_{10} N$ | $c_n = N^{1/n^2}$ |
|---|---|---|---|
| 6 | 36 | 3,994 | 1,2911 |
| 8 | 64 | 13,123 | 1,6034 |
| 10 | 100 | 22,380 | 1,6742 |
| 12 | 144 | 33,111 | 1,6980 |

Constantes implicadas por pares consecutivos (i.e. $10^{\Delta \log N / \Delta n^2}$):

| intervalo | constante local |
|---|---|
| $6 \to 8$ | **2,119** |
| $8 \to 10$ | **1,808** |
| $10 \to 12$ | **1,753** |

⚠️ **As constantes locais decrescem monotonicamente.** O ajuste OLS de $\log_{10} N = \alpha n^2 + \beta$ nos 4 pontos dá $c = 10^\alpha = 1{,}848$, mas com **resíduos de $-0{,}82,\ +0{,}84,\ +0{,}49,\ -0{,}51$ em $\log_{10}$** — ou seja, erros de até um **fator 7** na contagem, com sinais alternados. Resíduo alternante e de grande magnitude é a assinatura clássica de **modelo mal especificado**: falta um termo de ordem inferior (efeito de borda, $\log N = \alpha n^2 + \beta n + \gamma$).

Testei o modelo de 3 parâmetros: é instável com 4 pontos (ajustando em $n=8,10,12$ dá $c = 1{,}528$ e erra $n=6$ por $+1{,}34$ em $\log_{10}$; ajustando em $n=6,8,10$ dá $c = 1{,}038$). **Os dados não determinam a constante assintótica.**

**Ação — reescrever a afirmação.** O enunciado atual "$N(n) \sim 1{,}82^{n^2}$" (linhas 151 e 666) é um artefato do ajuste e deve sair. Substituir por algo como:

> A sequência $c_n = N(n)^{1/n^2}$ é crescente em nossos dados ($1{,}29,\ 1{,}60,\ 1{,}67,\ 1{,}70$ para $n=6,8,10,12$), com constantes locais decrescentes ($2{,}12,\ 1{,}81,\ 1{,}75$), sugerindo um limite na faixa $1{,}7$–$1{,}85$ — consistente com, e substancialmente mais estreito que, o intervalo provado $[1{,}3535,\ 4]$ de Kyek–Parberry–Wegener. Ressalvamos que quatro pontos (dois deles estimativas com CV $\approx 150$) não permitem separar $\alpha n^2$ de um termo de borda $\beta n$; não afirmamos uma constante assintótica.

Isso é mais honesto **e mais forte**: você passa de uma extrapolação frágil para um estreitamento defensável de um intervalo provado.

**Referências adicionais colhidas da bibliografia do KPW** (todas relevantes e ausentes do seu paper):
- **A. Conrad, T. Hindrichs, H. Morsy, I. Wegener**, *Solution of the knight's Hamiltonian path problem on a chessboard*, Discrete Appl. Math. **50** (1994), 125–134. — Caracteriza existência de caminhos hamiltonianos $s \to t$; **relevante à sua §8** (enumeração de caminhos) e à sua interface `knight_path`.
- **M. Kraitchik**, *Mathematical Recreations*, Cap. 11, Dover, 1953. — Lower bound histórico para $n=8$.
- **P. Cull, J. De Curtins**, *Knight's Tour Revisited*, Fibonacci Quarterly **16** (1978), 276–285. — Prova que todo $m \times n$ com $m,n \geq 5$ admite tour aberto e caracteriza os fechados; anterior a Schwenk.

Relacionado: **I. Parberry (1997)**, *An efficient algorithm for the Knight's tour problem*, Discrete Appl. Math. **73**(3), 251–260. DOI `10.1016/S0166-218X(96)00010-8`. É um algoritmo **divide-and-conquer** para tours do cavalo — **sua §3.2 (DnC via blocos $6\times6$) precisa citar e diferenciar.** Parberry constrói *um* tour em tempo linear por D&C; você **enumera/conta** a subclasse decomponível. Diferença real, mas precisa ser dita.

### 2.5 Matriz de transferência / broken-profile DP

Sua §9 diz "tema clássico em combinatória computacional" sem citar ninguém. Referência central:

> J. L. Jacobsen, **"Exact enumeration of Hamiltonian circuits, walks and chains in two and three dimensions"**, *J. Phys. A: Math. Theor.* **40**(49), 14667, 2007. DOI `10.1088/1751-8113/40/49/003`.

Também:
- **V. H. Pettersson**, *Enumerating Hamiltonian Cycles*, Electron. J. Combin. **21**(4), 2014. DOI `10.37236/4510`.
- **S. Minato**, *Zero-suppressed BDDs and their applications*, STTT **3**(2), 2001 — a técnica ZDD que Löbbing–Wegener usaram; contextualiza por que seu C++ estourou em $n=8$ com 187M estados.

### 2.6 Conectividade como obstrução: *Cut & Count* é o análogo exato

Este é conceitualmente o mais interessante que encontrei fora do seu radar. Todo o seu paper argumenta que **conectividade é a obstrução global que GF(2) não captura** (Lema `xor-overlapzero`, razão tours/2-fatores decaindo, `residual_search`). Existe uma técnica canônica que ataca exatamente isso:

> M. Cygan, J. Nederlof, M. Pilipczuk, M. Pilipczuk, J. M. M. van Rooij, J. O. Wojtaszczyk, **"Solving Connectivity Problems Parameterized by Treewidth in Single Exponential Time"**, FOCS 2011, DOI `10.1109/FOCS.2011.23`; versão de journal ACM TALG **18**(2), 2022, DOI `10.1145/3506707`.
> H. L. Bodlaender, M. Cygan, S. Kratsch, J. Nederlof, *Deterministic single exponential time algorithms for connectivity problems parameterized by treewidth*, Inf. Comput. **243**, 2015. DOI `10.1016/j.ic.2014.12.008`.

**Cut & Count** conta soluções **módulo 2** justamente para que as soluções *desconexas* se cancelem aos pares (cada uma é contada um número par de vezes pelos cortes consistentes), deixando só as conexas. É literalmente "usar $\mathbb{F}_2$ para capturar conectividade" — a coisa que você conclui ser impossível na sua formulação.

**Ação:** parágrafo novo em *Trabalhos relacionados* + menção nas conclusões da §8 (XOR). Não invalida nada seu — o Cut&Count opera sobre uma decomposição em árvore com o *Isolation Lemma*, não sobre o espaço de ciclos do grafo. Mas mostra que a barreira que você identificou tem uma rota conhecida de contorno, e **não conhecer isso é o tipo de coisa que um parecerista de TCS nota imediatamente**.

### 2.7 Complexidade: falta a referência de NP-completude

> A. Itai, C. H. Papadimitriou, J. L. Szwarcfiter, **"Hamilton Paths in Grid Graphs"**, *SIAM J. Comput.* **11**(4), 676–686, 1982. DOI `10.1137/0211056`. (515 citações)

Necessária para justificar por que backtracking/heurísticas são a abordagem certa. Citar na §1.1 ou §3.1.

### 2.8 Heurísticas: Warnsdorff sozinho é insuficiente

Você cita `Warnsdorff1823` mas não a análise moderna:
- **I. Pohl**, *A method for finding Hamilton paths and Knight's tours*, CACM **10**(7), 446–449, 1967. DOI `10.1145/363427.363463` — o desempate de Pohl para a regra de Warnsdorff.
- **S. Marateck**, *How good is the Warnsdorff's knight's tour heuristic?*, arXiv:0803.4321, 2008.
- **Squirrel & Cull (2000)**, análise de quando Warnsdorff falha.

Relevante porque sua heurística $f_\infty(L)$ é uma **regra de ordenação estática por camada**, e Warnsdorff é uma **regra dinâmica por grau residual**. Você deveria comparar as duas explicitamente na §3.1 — e, dado que sua memória registra que o ganho real vem de *pressão de vértice* (≈ Warnsdorff) e **não** de $f_\infty$ (benchmark refutado, speedup flat ~1,0), essa comparação é honesta e informativa.

### 2.9 Retangulares: sua conjectura $Q(n,m)=3$ tem vizinhança

> G. L. Chia, S.-H. Ong, *Generalized knight's tours on rectangular chessboards*, Discrete Appl. Math. **150**(1–3), 80–98, 2005. DOI `10.1016/j.dam.2004.11.008`.

E, para as extensões que você menciona em trabalho futuro (§12.3, "grafos cavalo-tipo"):
- J. DeMaio, *Which Chessboards have a Closed Knight's Tour within the Cube?*, Electron. J. Combin. **14**, 2007. DOI `10.37236/950`.
- J. Erde, B. Golénia, S. Golénia, *The Closed Knight Tour Problem in Higher Dimensions*, Electron. J. Combin. **19**, 2012. DOI `10.37236/2272`.

### 2.10 Cultural / survey

> N. D. Elkies, R. P. Stanley, *The mathematical knight*, *Math. Intelligencer* **25**(1), 22–34, 2003. DOI `10.1007/BF02985635`.

Survey de alto nível. Barato de citar, sinaliza domínio da área.

---

## 3. ERRO FACTUAL a corrigir

### 3.1 Enunciado de Schwenk está incorreto (linhas 239–244)

Seu texto:

> "todos exceto aqueles com (a) $m$ e $n$ ambos ímpares, (b) $m \in \{1,2\}$, (c) $m = 3$ e $n \in \{4,6,8\}$, ou (d) $m = n = 4$"

Enunciado correto de Schwenk (1991), para $m \le n$ — não existe tour fechado sse:

> (a) $m$ e $n$ ambos ímpares; (b) $m \in \{1, 2, \mathbf{4}\}$; (c) $m = 3$ e $n \in \{4,6,8\}$.

**O erro:** você restringiu a exclusão do caso $m=4$ a $m=n=4$. Está errado — **nenhum** tabuleiro $4 \times n$ admite tour fechado, para qualquer $n$ (é o argumento clássico de coloração de Pósa em 4 cores). Sua cláusula (d) deve ser eliminada e (b) deve virar $m \in \{1,2,4\}$.

Impacto baixo no resto do paper (você só trabalha com $n \ge 6$), mas é o tipo de erro num enunciado clássico que destrói confiança imediatamente.

---

## 4. Reivindicações de novidade a moderar

| Local | Texto atual | Problema | Sugestão |
|---|---|---|---|
| L. 267–274 | *"Nossa contribuição é aplicar sistematicamente este framework ao passeio do cavalo"* | Ok, mas incompleto: omite que o framework tem nome (*Hamilton space*) e literatura ativa. | Manter, mas emendar com o parágrafo de §1.1 acima e posicionar o cavalo como o caso **esparso** que falta na literatura. |
| L. 287–292 | *"é, ao nosso conhecimento, a primeira formalização de resultados sobre o espaço de ciclos GF(2) do passeio do cavalo"* | Cuidadosamente hedged e **provavelmente verdadeira** — não achei nada contrário. | ✅ Manter como está. É a formulação certa. |
| §6 (superfícies) | Hierarquia $Q$: plano 3, cilindro 0, toro 0; Klein | Prior art não citado (Watkins, Forrest–Teehan, Forrest–Lague). | Ver §1.4 — reposicionar como validação cruzada. |
| §5 (estimador) | Estimador de Knuth para $N(n)$ | Prior art não citado (Cancela–Mordecki). | Ver §2.3 — diferenciar por método e por tours fechados. |
| §3.2 (DnC) | DnC por blocos $6\times6$ | Parberry 1997 é DnC para o mesmo problema. | Ver §2.4 — diferenciar construção-de-um vs. contagem-de-classe. |
| L. 151, 666 | $N(n) \sim 1{,}82^{n^2}$ | Ajuste a 4 pontos, dois deles estimados com CV≈150, apresentado sem barra de erro na constante e sem confronto com bounds provados. | Reportar IC para a constante; confrontar com Kyek–Parberry–Wegener (verificar constantes — §2.4). |
| §8 (XOR paths) | *"a desconexão é a obstrução universal"* | Verdadeiro no seu setup, mas Cut&Count mostra rota conhecida de contorno. | Citar Cygan et al. e qualificar: "no espaço de ciclos; cf. Cut&Count para uma rota alternativa via decomposição em árvore". |

---

## 5. Bibliografia proposta — deltas

**Atual:** 10 itens. **Proposta:** ~28. Para um paper de 36 páginas com reivindicações de novidade em 4 frentes, 10 é baixo demais e por si só sinaliza revisão bibliográfica insuficiente.

### Adicionar (ordenado por prioridade)

```
P0 — bloqueiam publicação                                    [✅ = referência verificada no texto-fonte]
✅[Heinig2013]     P. Heinig. When Hamilton circuits generate the cycle space of a random graph.
                     arXiv:1303.0026, 2013.   <-- OBSERVA A NECESSIDADE DE delta(G) >= 3. Cite no Teorema Q(n)=3.
✅[Heinig2014]     P. Heinig. On prisms, Möbius ladders and the cycle space of dense graphs.
                     European Journal of Combinatorics 36 (2014), 503-530.
✅[Hartman1983]    I. B.-A. Hartman. Long cycles generate the cycle space of a graph.
                     European Journal of Combinatorics 4 (1983), 237-246.
✅[HouYin2025]     X. Hou, Z. Yin. Dirac-type condition for Hamilton-generated graphs.
                     arXiv:2503.15950v1, 2025.
✅[CNP2026]        M. Christoph, R. Nenadov, K. Petrova. The Hamilton space of pseudorandom graphs.
                     J. Combin. Theory Ser. B 176 (2026), 254-267.   <-- PUBLICADO; nao cite o arXiv
✅[HK2025a]        D. Hefetz, M. Krivelevich. The Hamilton cycle space of random graphs.
                     arXiv:2506.19731, 2025.
✅[HK2025b]        D. Hefetz, M. Krivelevich. The Hamilton cycle space of random regular graphs and
                     randomly perturbed graphs. arXiv:2507.04488, 2025.
✅[HamGen2026]     D. Hefetz, M. Krivelevich. On graphs whose cycle space is spanned by their
                     Hamilton cycles. arXiv:2606.05835, 2026.
✅[Watkins2004]    J. J. Watkins. Across the Board: The Mathematics of Chessboard Problems.
                     Princeton UP, 2004, pp. 65-77.
✅[WatkinsHoenigman1997] J. J. Watkins, R. Hoenigman. Knight's Tour on a Torus.
                     Mathematics Magazine 70(3) (1997), 175-184.
✅[ForrestTeehan2015] B. Forrest, K. Teehan. The Topology of Knight's Tours on Surfaces.
                     arXiv:1507.02917, 2015.
 [ForrestLague2024]  B. Forrest, Lague. Nullhomotopic and Generating Knight's Tours on
                     Non-Orientable Surfaces. arXiv:2406.05226, 2024.

P1 — credibilidade
[LobbingWegener1996] Löbbing, Wegener. Electron. J. Combin. 3(1), R5, 1996. doi:10.37236/1229. [contagem incorreta; ver McKay1997]
[Wegener2000]     Wegener. Branching Programs and Binary Decision Diagrams. SIAM, 2000, p. 369.
[OEIS_A001230]    OEIS Foundation. Sequence A001230. https://oeis.org/A001230
[Knuth1975]       Knuth. Estimating the efficiency of backtrack programs. Math. Comp. 29(129):121-136, 1975.
[Chen1992]        Chen. Heuristic Sampling. SIAM J. Comput. 21(2), 1992. doi:10.1137/0221022
[CancelaMordecki2006] Cancela, Mordecki. Counting Knight's Tours through the Randomized Warnsdorff Rule. arXiv:math/0609009, 2006.
✅[KyekParberryWegener1997] Kyek, Parberry, Wegener. Bounds on the number of knight's tours.
                     DAM 74(2):171-181, 1997.  [Thm 4: Omega(1.3535^{n^2}); Rmk 6: <= 4^{n^2};
                     Thm 7: N(8) <= 3.019e22]
✅[Conrad1994]     A. Conrad, T. Hindrichs, H. Morsy, I. Wegener. Solution of the knight's
                     Hamiltonian path problem on a chessboard. DAM 50 (1994), 125-134.
✅[Kraitchik1953]  M. Kraitchik. Mathematical Recreations, Ch. 11. Dover, 1953.
✅[CullDeCurtins1978] P. Cull, J. De Curtins. Knight's Tour Revisited.
                     Fibonacci Quarterly 16 (1978), 276-285.
[Parberry1997]    Parberry. An efficient algorithm for the Knight's tour problem. DAM 73(3):251-260, 1997.
[Jacobsen2007]    Jacobsen. Exact enumeration of Hamiltonian circuits... J. Phys. A 40(49):14667, 2007.
[Cygan2011]       Cygan et al. Solving Connectivity Problems Parameterized by Treewidth... FOCS 2011.
[Itai1982]        Itai, Papadimitriou, Szwarcfiter. Hamilton Paths in Grid Graphs. SIAM J. Comput. 11(4), 1982.
[Pohl1967]        Pohl. A method for finding Hamilton paths and Knight's tours. CACM 10(7), 1967.

P2 — completude
[Pettersson2014]  Pettersson. Enumerating Hamiltonian Cycles. Electron. J. Combin. 21(4), 2014.
[Minato2001]      Minato. Zero-suppressed BDDs and their applications. STTT 3(2), 2001.
[ChiaOng2005]     Chia, Ong. Generalized knight's tours on rectangular chessboards. DAM 150, 2005.
[ElkiesStanley2003] Elkies, Stanley. The mathematical knight. Math. Intelligencer 25(1), 2003.
[Marateck2008]    Marateck. How good is the Warnsdorff's knight's tour heuristic? arXiv:0803.4321, 2008.
[Erde2012]        Erde, Golénia, Golénia. The Closed Knight Tour Problem in Higher Dimensions. E-JC 19, 2012.
[Bodlaender2015]  Bodlaender, Cygan, Kratsch, Nederlof. Inf. Comput. 243, 2015.
```

### Corrigir os existentes

- `Schwenk1991` — enunciado errado no corpo do texto (§3.1 acima).
- `Knuth2022` — despromover a secundária; primária é Knuth 1975.
- `Reconfig2022` — Ito et al., *Rerouting Planar Curves and Disjoint Paths*, arXiv:2210.11778. ⚠️ **Verificar se saiu versão publicada** (provavelmente ICALP 2023); se sim, citar a publicada.
- `McKay1997` — adicionar link estável: `http://hdl.handle.net/1885/40759`, e mencionar que corrige Löbbing–Wegener.
- `Grinberg1968` — verificar a referência; a grafia usual é *Latvian Math. Yearbook* **4**, 51–58, 1968, em russo. Confirmar que você realmente usa a fórmula (hoje aparece só como "resultado relacionado" na linha 310 e 1677) — se for só menção de contexto, tudo bem; se sustenta algum argumento, precisa ser explicitada.

---

## 6. Plano de ação sugerido

**Ordem de execução** (do maior retorno por esforço):

0. **Reescrever a lei $1{,}82^{n^2}$** (30 min) — §2.4. Não é citação, é estatística; é o alvo mais fácil de um parecerista e a correção só te fortalece (estreitar $[1{,}3535,\ 4]$ é um resultado melhor que uma extrapolação frágil).
1. **Corrigir Schwenk** (5 min, erro factual puro).
2. **Escrever o parágrafo de *Hamilton space*** em Trabalhos Relacionados + reposicionar o abstract em torno de "deficit constante em família esparsa" (2h). Maior ganho isolado do paper.
3. **Reposicionar §6 (superfícies)** como validação cruzada, citando Watkins/Forrest (1h). Maior redução de risco.
4. **Corrigir a atribuição do estimador**: Knuth 1975 + Cancela–Mordecki com diferenciação explícita (30 min).
5. **Adicionar Löbbing–Wegener/McKay/OEIS** na narrativa de contagem (30 min).
6. **Citar Parberry 1997** na §3.2 e diferenciar do seu DnC (30 min).
7. **Adicionar parágrafo Cut&Count** (45 min). Alto valor intelectual — e potencialmente uma direção de pesquisa nova para a `tightness`.
8. **Verificar as constantes de Kyek–Parberry–Wegener** e confrontar com $1{,}82^{n^2}$ (pendente — ver §2.4).
9. Preencher P2 (1h).

---

## 7. Pendências

### ✅ Resolvidas (2ª rodada — via `pdftotext` sobre os PDFs originais)

- [x] **Constantes de Kyek–Parberry–Wegener** → §2.4. Expoente é $c^{n^2}$. Intervalo provado $[1{,}3535,\ 4]$. **Consequência inesperada: a lei $1{,}82^{n^2}$ do paper é artefato de ajuste e precisa ser reescrita.**
- [x] **Referência exata de Heinig** → duas entradas ([21] EJC 2014, [22] arXiv:1303.0026 2013). **Consequência inesperada: Heinig [22] já observa que $\delta \geq 3$ é necessário — a obstrução de grau 2 é prior art documentada.** Ver §1.2 reescrito.
- [x] **Forrest–Teehan texto completo** (24 pp.) → §1.4. Sobreposição é sobre *existência* por classe de homotopia; sua contribuição de *equidistribuição* sobrevive.

### ⬜ Ainda abertas

- [x] **Versão publicada de Ito et al.** → ✅ resolvido via DBLP. Existem **três** versões; cite a de journal:
  > T. Ito, Y. Iwamasa, N. Kakimura, Y. Kobayashi, S. Maezawa, Y. Nozaki, Y. Okamoto, K. Ozeki, *Rerouting Planar Curves and Disjoint Paths*, **ACM Trans. Algorithms 21 (2025)**, DOI `10.1145/3715694`.
  > (versão de conferência: ICALP 2023, DOI `10.4230/LIPIcs.ICALP.2023.81`; preprint: arXiv:2210.11778)

  Seu `\bibitem{Reconfig2022}` cita o preprint de 2022 — **desatualizado em duas gerações**. Renomeie para `Reconfig2025`.
- [ ] **Forrest–Lague 2024** (arXiv:2406.05226, Klein/Möbius) — só li o abstract. Baixar e comparar com seus achados de Klein (`n+m` par ⇒ $Q=1$; off-diagonal $(0,1),(1,0)$ sempre proibido). É o único ponto onde ainda pode haver sobreposição não mapeada.
- [ ] Confirmar se algum trabalho da linha *Hamilton space* já trata explicitamente grafos **esparsos com deficit positivo exato**. Busquei e não achei — a literatura toda vive no regime denso/pseudoaleatório onde o deficit é 0. Mas a área se move rápido (JCTB 2026, arXivs de 2025–2026); vale um alerta salvo no arXiv/Google Scholar.
- [ ] **Parberry, *Algorithms for touring knights*** (Tech. Report, Univ. North Texas, 1994) — é a ref [6] do KPW e a base do D&C. Verificar relação com o `Parberry1997` publicado antes de citar os dois.

---

## 8. Acesso: o que é grátis e o que precisa da UNIFEI

### 8.1 Sobre a chave que você passou

A chave `vFE9FBWmkXNS9pObK8DTyV` **não foi usada** — o **OpenAlex não exige API key**. Ele usa "polite pool" por e-mail: basta `?mailto=seu@email`, que dá rate limit melhor sem autenticação. Toda a auditoria rodou assim. (A chave tem cara de CORE ou Semantic Scholar; guarde para esses.) Como você mesmo disse que vai rotacionar, sem problema — mas registro que era desnecessária.

### 8.2 As APIs que resolveram tudo — todas sem key ou com key trivial

| Serviço | Key? | Para quê |
|---|---|---|
| **OpenAlex** `api.openalex.org` | ❌ só `mailto=` | Descoberta ampla, contagem de citações, metadados. Foi o cavalo de batalha aqui. |
| **DBLP** `dblp.org/search/publ/api` | ❌ | **Achar a versão publicada de um preprint.** Resolveu o Ito et al. (arXiv → ICALP'23 → TALG'25). Use sempre antes de citar arXiv. |
| **Unpaywall** `api.unpaywall.org/v2/{doi}?email=` | ❌ só e-mail | **Cópia OA legal de artigo pago.** Achou o KPW e o Jacobsen em PDF aberto — nenhum dos dois precisou de CAPES. |
| **Crossref** `api.crossref.org` | ❌ | Metadados canônicos p/ montar `.bib` correto. |
| **arXiv API** `export.arxiv.org/api/query` | ❌ | Preprints. |
| **OEIS** `oeis.org/search?fmt=json` | ❌ | A001230 (contagens de tours). |
| **zbMATH Open** `api.zbmath.org` | ❌ (aberto desde 2021) | **Específico de matemática** — reviews e classificação MSC. Ideal p/ este projeto. |
| **Semantic Scholar Graph** | 🔸 opcional, grátis | Rate limit maior; grafo de citações. Talvez sua chave sirva aqui. |
| **CORE** `core.ac.uk/services/api` | 🔸 grátis c/ cadastro | Agregador de repositórios; fallback do Unpaywall. |

**Fluxo recomendado:** OpenAlex (descobre) → DBLP (versão publicada) → Unpaywall (PDF legal) → só então CAPES.

Ferramenta local que fez a diferença: **`pdftotext -layout`** (poppler, já instalado). As três pendências caíram porque o `WebFetch` falha em converter PDFs acadêmicos, mas `curl | pdftotext` extrai perfeitamente. Todas as constantes do KPW vieram daí.

### 8.3 O que realmente precisa do acesso UNIFEI

Ponto importante: **nenhuma das três pendências precisou de paywall.** Todas estavam em arXiv ou na página pessoal do autor. O que sobra que de fato exige assinatura:

| Item | Onde | Precisa de acesso? |
|---|---|---|
| Kyek–Parberry–Wegener 1997 | DAM/Elsevier | ❌ OA via Unpaywall + cópia em `ianparberry.com/pubs/boundsknight.pdf` |
| Jacobsen 2007 | J. Phys. A/IOP | ❌ OA via Unpaywall |
| Forrest–Teehan, Forrest–Lague, todos os *Hamilton space* | arXiv | ❌ grátis |
| **Conrad et al. 1994** | DAM/Elsevier | ✅ **sim** — Unpaywall diz `is_oa=False` |
| **Schwenk 1991** | Mathematics Magazine/JSTOR | ✅ sim (você já cita; vale conferir o enunciado — §3.1) |
| **Watkins–Hoenigman 1997** | Mathematics Magazine/JSTOR | ✅ sim |
| **Watkins 2004** (livro) | Princeton UP | ✅ biblioteca / e-book institucional |
| **Chia–Ong 2005**, **Hartman 1983**, **Heinig 2014** | Elsevier | ✅ sim |
| Cull–De Curtins 1978 | Fibonacci Quarterly | ❌ o FQ é aberto: `fq.math.ca` |
| Kraitchik 1953 | Dover (livro) | 🔸 tentar `archive.org` |

**Rota de acesso no Brasil:** **Portal de Periódicos CAPES** com login federado **CAFe** usando credenciais UNIFEI — funciona fora do campus, é o caminho certo (eduroam só te dá acesso por IP dentro do campus). Cobre Elsevier, JSTOR e Taylor & Francis, que é exatamente o conjunto acima.

Vale checar também se a UNIFEI assina **MathSciNet** (AMS): para este projeto seria útil, porque os *reviews* do MathSciNet frequentemente dizem explicitamente "este resultado generaliza X" — que é justamente o tipo de conexão que uma auditoria de prior art precisa e que buscas por palavra-chave não pegam.

### 8.4 Sugestão de higiene

Vale montar um `paper/refs.bib` e migrar de `thebibliography` manual para BibTeX — com ~28 referências, manter à mão vira fonte de erro (o enunciado errado do Schwenk e o Ito desatualizado são sintomas disso). Crossref entrega BibTeX pronto:

```bash
curl -sLH "Accept: application/x-bibtex" "https://doi.org/10.1145/3715694"
```
