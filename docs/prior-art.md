# Prior art

Registro honesto do que já existia e do que sobrou. Este documento existe
porque a descoberta de que boa parte do trabalho já estava publicada foi o
evento que mudou o rumo do projeto, e escondê-la tornaria o repositório
enganoso.

A auditoria completa, com citações e trechos dos originais, está em
`docs/notas-pesquisa/citation_audit.md` (agosto de 2026) e nos dois índices
de bibliografia em `docs/notas-pesquisa/bibliografia-{teoria,algoritmos}.md`.

---

## 1. O algoritmo: totalmente prior art

O solver de `core/knight_tours.py` (backtracking em variáveis-aresta +
propagação de grau 2 + union-find incremental de sub-ciclos) é o
**multi-path method**, com correspondência termo a termo:

| neste repositório | literatura |
|---|---|
| `FREE` / `ACTIVE` / `INACTIVE` | *undecided* / *required* / *deleted*, rotulação de **Rubin (1974)** |
| `avail == 2` força aresta | *forced vertex* → *forced edge* |
| `_process_queue` até esvaziar | a *chain reaction* até consistência |
| union-find de segmentos | **Kocay (1992)** |

Referências: Rubin (1974), Christofides (1975), Kocay (1992), Chalaturnyk (2008).

## 2. O XOR de ciclos: prior art, e o resultado negativo também

A ideia de combinar ciclos fundamentais por XOR é o **circuit vector space
method**, de Welch (1966) e Mateti–Deo (1976).

Mais duro: **o resultado negativo já era teorema em 1976.** Mateti e Deo
provaram que a razão tours / 2-fatores tende a zero e que enumerar por uniões
de circuitos é desperdício assintótico. A medição independente feita aqui
(`r ≈ 0.84 · e^(-0.17n)`) confirma um teorema de 49 anos atrás.

E ambos os métodos **já foram aplicados ao grafo do cavalo especificamente**.
O 6×6 com 9862 tours é linha de tabela publicada.

Outros itens da mesma categoria:

- flips por face sobre tours = *Z-transformation*, **Sheffield (2000)**;
- cláusulas XOR em solvers = **DPLL(XOR)**, técnica padrão — o speedup medido
  de 1,00x era previsível.

## 3. O enquadramento tem nome próprio: Hamilton space

Este foi o achado mais importante da auditoria, e não estava no radar do
projeto.

O que aqui se chamou `Ham(n) ⊆ Z_1` e `deficit = beta_1 - rank(Ham)` é, na
literatura, o **Hamilton space** `C_n(G)` e sua **codimensão** em `C(G)`. Um
grafo com `C_n(G) = C(G)` é dito **Hamilton-generated**. É uma área ativa:

| trabalho | relevância |
|---|---|
| Heinig (2013), *When Hamilton circuits generate the cycle space of a random graph* | origem da pergunta |
| Hou & Yin (2025), *Dirac-type condition for Hamilton-generated graphs* | resolve o caso Dirac; cunha o termo |
| Hefetz & Krivelevich (2025), duas notas sobre grafos aleatórios e regulares | estado da arte |
| Christoph, Nenadov & Petrova (2024), *The Hamilton space of pseudorandom graphs* | introduz **parity-switchers** |
| *On graphs whose cycle space is spanned by their Hamilton cycles* (2026) | survey da área |

### 3.1 `deficit > 0` também é prior art

Heinig observou que `delta(G) >= 3` é **necessário** para Hamilton-generation.
O argumento é elementar: se `deg(v) = 2` com arestas `e_1, e_2`, todo ciclo
hamiltoniano contém ambas, logo todo tour satisfaz `x_{e1} = x_{e2}` — um
hiperplano próprio contendo `C_n(G)`.

O grafo do cavalo tem `delta = 2` nos quatro cantos. Portanto **a existência
da obstrução é conhecida**, não é descoberta deste projeto.

## 4. Superfícies: parcialmente escooped

A hierarquia de `Q` por superfície (plano 3, cilindro 0, toro 0, Klein) não
foi construída no vácuo, embora tenha sido derivada sem conhecer a literatura:

| trabalho | o que cobre |
|---|---|
| Watkins (2004), *Across the Board* | capítulos sobre toro, cilindro, Möbius e Klein |
| Forrest & Teehan (2015), *The Topology of Knight's Tours on Surfaces* | caracterizam quais dimensões de cilindro e toro admitem tours realizando a identidade de `pi_1` versus um gerador |
| Forrest & Lague (2024), *Nullhomotopic and Generating Knight's Tours on Non-Orientable Surfaces* | mesmo programa para Möbius e Klein |

**A diferença que sobrevive:** esses trabalhos respondem uma pergunta de
**existência** ("existe tour na classe X?"). As medições deste projeto
respondem uma pergunta de **distribuição** ("que fração dos tours cai em cada
classe?"): 4 classes com ~25% cada no toro em 20k tours, winding 50/50 no
cilindro em 18k tours. **Equidistribuição é estritamente mais forte que
existência**, e não aparece nesses papers.

---

## 5. O que sobra como contribuição

Com as ressalvas acima já descontadas:

1. **O valor exato `deficit = 3`**, não apenas `> 0`. O colapso `4 -> 3`
   (quatro cantos, três dimensões, via o núcleo do funcional soma) é onde
   está a matemática.
2. **A constância em `n`** — o deficit não cresce com o tabuleiro.
3. **A localidade**: `Q(G_T \ S) = 3 <=> W_corners ⊆ S`, e
   `Q = max(0, k_deg2 - 1)`.
4. **A conexidade do bulk** (`n >= 6`, indução `n -> n+2`), ingrediente
   geométrico não-trivial.
5. **A formalização em Lean 4**, sem `sorry`.
6. **Tightness** (`rank(Ham) = beta_1 - 3` exatamente), provada por construção
   em `n ∈ {8,10,12}`.
7. **Equidistribuição** por classe de homologia nas superfícies (§4).
8. **Os resultados negativos documentados** — em especial o diagnóstico de que
   conectividade não é uma condição linear sobre `F_2^E`, testado e confirmado
   em cinco domínios.

### O reposicionamento que a auditoria sugere

A literatura de Hamilton space opera no regime **denso/pseudoaleatório**
(condições tipo Dirac, `delta(G) >~ n/2`), onde o deficit tipicamente é zero.
O grafo do cavalo está no extremo oposto: **esparso** (`|E| ≈ 4|V|`), grau
máximo 8, grau mínimo 2.

Lido assim, o resultado deixa de ser "um invariante ad hoc" e passa a ser
**o primeiro exemplo estruturado de deficit positivo, exato e constante numa
família esparsa natural, determinado por geometria local** — complementar à
literatura existente, não concorrente com ela.

## 6. Estado das pendências da auditoria

A auditoria é de 06/08/2026. O paper corrente (`knight_tour_tightness`, v5,
15/08/2026) já a incorporou. Registro do que foi feito:

| pendência | estado |
|---|---|
| citar a literatura de Hamilton space | **resolvido** — Heinig 2013/2014, CNP2026, Hou–Yin, Hefetz–Krivelevich a/b, HamGen2026 na bibliografia (31 itens) |
| dizer que `deficit > 0` é prior art | **resolvido** — os cantos passam a dar explicitamente só a *cota superior* |
| adotar parity-switchers | **resolvido** — e com um ganho: provada a **dicotomia canto/bulk** (um parity-switcher ancorado em canto é impossível em `G_n`, porque o canto teria grau 2 e seus dois vizinhos são removidos no passo S3) |
| ler Hartman (1983) | **resolvido** — citado |
| superfícies: Watkins, Forrest–Teehan | **resolvido**, com distinção explícita: o `Q` deste trabalho é o **invariante algébrico** de vértices de grau 2; o `Q` topológico da literatura é por classes de winding. Coincidem no plano e divergem alhures |
| lei `N(n) ~ 1,82^(n²)` | **retirada** do texto, com registro no changelog |

### Um item que a auditoria não previu

A auditoria não pedia, mas o paper encontrou: o deficit nulo no toro **não é
acidente da superfície**. O grafo do cavalo toroidal é um grafo de Cayley sobre
`Z_n × Z_n`, e deficit zero é corolário de um teorema de **Alspach, Locke e
Witte (1990)**. O que elimina o deficit é a **transitividade por vértices** —
toda aresta equivalente a toda outra — e não a topologia. No plano, a fronteira
quebra a transitividade e cria exatamente os vértices de grau 2 que geram
`Q = 3`.

Plano versus toro é, em essência, **com fronteira versus sem fronteira**.

### E um erro da minha leitura anterior deste repositório

Uma versão anterior deste documento afirmava que a hierarquia de `Q` por
superfície era contribuição própria. Não é, na forma em que estava escrita —
e o paper corrente já a reposiciona corretamente. O que sobrevive ali é a
distinção algébrico/topológico e a medição de **distribuição** (equidistribuição
por classe) onde a literatura estabelece **existência**.

## Referências

- Rubin, F. (1974). *A search procedure for Hamilton paths and circuits*.
- Christofides, N. (1975). *Graph Theory: An Algorithmic Approach*.
- Welch, J. T. (1966). Circuit vector space method.
- Mateti, P., Deo, N. (1976). *On algorithms for enumerating all circuits of a graph*.
- Hartman, I. B.-A. (1983). *Long cycles generate the cycle space of a graph*. European J. Combin. 4, 237–246.
- Kocay, W. (1992). *An algorithm for finding Hamiltonian cycles*.
- Sheffield, S. (2000). Z-transformation.
- Watkins, J. J. (2004). *Across the Board: The Mathematics of Chessboard Problems*. Princeton UP.
- Chalaturnyk, A. (2008). *A fast algorithm for finding Hamilton cycles* (tese).
- Heinig, S. (2013). *When Hamilton circuits generate the cycle space of a random graph*. arXiv:1303.0026.
- Forrest, B., Teehan, K. (2015). *The Topology of Knight's Tours on Surfaces*. arXiv:1507.02917.
- Christoph, M., Nenadov, R., Petrova, K. (2024). *The Hamilton space of pseudorandom graphs*. arXiv:2402.01447.
- Forrest, B., Lague (2024). *Nullhomotopic and Generating Knight's Tours on Non-Orientable Surfaces*. arXiv:2406.05226.
- Hou, Yin (2025). *Dirac-type condition for Hamilton-generated graphs*. arXiv:2503.15950.
- Hefetz, D., Krivelevich, M. (2025). *The Hamilton cycle space of random graphs*. arXiv:2506.19731; *...random regular and randomly perturbed graphs*. arXiv:2507.04488.
