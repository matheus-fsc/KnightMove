# Prior art **algorítmico** — `knight_tours.py` e as frentes XOR

Biblioteca separada de `../referencia/` (que guarda a linha **teórica**:
espaço de Hamilton, superfícies, contagem).

Montada em 2026-09-01 para responder duas perguntas:

1. O método de otimização de `knight_tours.py` (backtracking em variáveis-aresta
   + propagação de grau 2 + Union-Find incremental de sub-ciclos) é inovador?
2. E o **XOR de loops** (espaço de ciclos GF(2))?

> **Veredito: não, em nenhuma das duas.** O algoritmo é o *multi-path method*
> (Rubin 1974 / Christofides 1975 / Kocay 1992); o XOR de ciclos é o *circuit
> vector space method* (Welch 1966 → Mateti–Deo 1976), cujo **resultado
> negativo já era teorema em 1976**. E ambos **já foram aplicados ao grafo do
> cavalo especificamente** — o 6×6 com 9 862 tours é linha de tabela publicada.
>
> Sobra como contribuição: $Q(n)=3$, tightness, formalização Lean. Ver `../referencia/` §1.1.

---

## 1. O algoritmo de `knight_tours.py` = *multi-path method*

### 1.1 Correspondência 1:1

| `knight_tours.py` | Multi-path (Rubin 74 / Christofides 75 / Kocay 92) |
|---|---|
| `FREE` / `ACTIVE` / `INACTIVE` | *undecided* / *required* / *deleted* (rotulação de Rubin) |
| `avail == 2 and dw < 2` → forçar `ACTIVE` | **forced vertex** → **forced edge** (Chalaturnyk §3.4) |
| `_process_queue` iterando até esvaziar a fila | a "**chain reaction**" até `M` ficar *consistent* (§3.4) |
| `CONTRADICTION` (`avail < 2`) | **Stop Condition 1A / 1B** (§3.4) |
| `uf_find(u) == uf_find(v)` → `SUBTOUR` / `COMPLETE_TOUR` | **Stop Condition 2** (§3.4) |
| `uf_parent` / `uf_size`, union by size | "*a **merge-find** like structure is maintained in the array*" — §5.1.3, atribuído a **[Koc92]**, detecção em **O(1)** |
| `_choose_next_edge` (prioriza grau 1 = ponta de segmento) | `ChooseVertex` / `ExtendSegments` (§6.2) |

### 1.2 Arquivos

| Arquivo | Trabalho | Papel |
|---|---|---|
| `Chalaturnyk2008_fast_algorithm_hamilton_cycles_thesis.pdf` | A. Chalaturnyk, *A Fast Algorithm for Finding Hamilton Cycles*, MSc, U. Manitoba, 2008 | ⭐⭐ **A referência central.** §3.3–3.4 = o algoritmo; §5.1.3 = o merge-find de Kocay; **Cap. 11 = testes em grafos do cavalo** (ver §3 abaixo). |
| `VandegriendCulberson1998_gnm_phase_transition_hc.pdf` | JAIR **9**, 219–245 (= arXiv:1105.5443) | §3: poda por vizinhos de grau 2 + caminhos forçados iterada até ponto fixo, + conectividade e cut-points. **§6 testa cavalo generalizado e o acha DIFÍCIL** (13 % estouraram 30 min). |
| `Vandegriend1998_finding_hamiltonian_cycles_thesis_SCAN.pdf` | B. Vandegriend, *Finding Hamiltonian Cycles: Algorithms, Graphs and Performance*, MSc, U. Alberta, 1998, 150 p. | O "[Van98]". Capítulo dedicado a cavalo generalizado $(A,B)_{n\times m}$. ⚠️ **Digitalização sem camada de texto** — `pdftotext` devolve vazio, precisa de `ocrmypdf`. |
| `CaseauLaburthe1997_solving_small_tsps_constraints.pdf` | Caseau & Laburthe, Bouygues | §3.1.3 `nocycle`: mesmo propagador de sub-ciclo, versão CP. Hoje é o global constraint `circuit`. |
| `SleegersVandenBerg2021_backtracking_algorithms_hcp.pdf` | Int. J. Adv. Intelligent Systems **14** (= arXiv:2107.00314) | Benchmark moderno: Vandegriend–Culberson segue sendo o baseline exato de referência. |

### 1.3 Não obtidos (paywall) — rota CAPES/CAFe

| Trabalho | Venue | Por que importa |
|---|---|---|
| **Rubin 1974**, *A search procedure for Hamilton paths and circuits* | JACM **21**(4) 576–580 | Origem da rotulação required/deleted/undecided. |
| **Kocay 1992**, *An extension of the multi-path algorithm for finding Hamilton cycles* | Discrete Math. **101**, 171–188 | O merge-find de segmentos + separadores e bipartições. |
| **Christofides 1975**, *Graph Theory: An Algorithmic Approach* | Academic Press | Exposição original do multi-path. |

### 1.4 Onde o nosso código está **atrás** de 1992/2008
- `State.snapshot()` copia 6 arrays por nó, O(V+E); Kocay/Chalaturnyk usam *in-place
  reduction* com trilha de desfazer (Chalaturnyk: estado total O(n+e)).
- Sem checagem de conectividade / pontos de articulação (VC98 tem).
- Sem detecção de conjuntos separadores nem de bipartições (Kocay 92 tem) —
  justamente a poda que ataca a barreira de conectividade que medimos.

---

## 2. XOR de loops

### 2.1 Como enumerador — o resultado negativo é de 1976

`MatetiDeo1975_enumerating_all_circuits_TR585.txt`
— P. Mateti & N. Deo, *On Algorithms for Enumerating All Circuits of a Graph*,
**SIAM J. Comput. 5(1) 90–99, 1976** (TR nº 585, OA em archive.org).

§2.1 é o método do `path_decomposition_experiment/`: base de ciclos
fundamentais → gerar os $2^\mu-1$ *circs* por **ring-sum** (XOR) → testar quais
são circuitos. E eles **já provaram** o que medimos:

> §4: "*the main source of **inefficiency** in the algorithms using the circuit
> vector space approach is in the **wasted computation of edge-disjoint unions
> of circuits**.*"

= nosso `bad_count` 73–84 % e os 2-fatores multi-ciclo do `residual_search`.

> §2.1 / Fig. 5: "*The **ratio of circuits to all vectors** in the vector space
> for numerous classes of graphs tends to **zero**.*"

= nosso $r(n) \approx 0{,}84\,e^{-0{,}17n}$ do `ratio_analysis`.

Mais: existem **apenas 4 grafos** cujo espaço de ciclos consiste só de circuitos.
Ancestrais comparados no mesmo survey: Welch (1966), Gibbs (1969), Rao–Murti.

### 2.2 Flips por face = **Z-transformation**

| Arquivo | Trabalho | Correspondência |
|---|---|---|
| `Sheffield2000_sampling_hamiltonian_cycles_Ztransform.pdf` | S. Sheffield, *Computing and Sampling Restricted Vertex Degree Subgraphs and Hamiltonian Cycles*, arXiv:math/0008231 | **Z-transformation** = "*symmetric difference consists of the boundary edges of a **single face***" = nosso flip hexagonal. Estuda $S_\phi$ com $\phi\equiv 2$ (2-fatores); **Teorema 3 caracteriza a conexidade do grafo de flips**; função-altura; hamiltonicidade de poliominós em $O(\lvert V\rvert^2)$; **cadeias de Markov para amostrar e contar tours**. |
| `TratnikYe2017_resonance_graphs_surfaces.pdf` | Tratnik & Ye, arXiv:1710.00761 | Survey da linha **resonance graph / Z-transformation graph** (Zhang–Guo–Chen, Discrete Math. **72**, 1988), incl. superfícies. Conexidade do grafo de flips é o tema central há ~40 anos. |

⚠️ **A brecha honesta:** Sheffield e a linha *resonance graph* exigem grafo
**bipartido mergulhado no plano ou no toro** (toda a máquina é a função-altura).
O grafo do cavalo **não é planar** → os teoremas não se aplicam diretamente.
O *método* não é novo; a conexidade do flip graph do cavalo não é corolário de ninguém.

### 2.3 XOR dentro do SAT

| Arquivo | Trabalho | O que já dizia |
|---|---|---|
| `LaitinenJunttilaNiemela2012_complete_parity_reasoning.pdf` | Laitinen, Junttila & Niemelä, arXiv:1207.0988 | Framework DPLL(XOR): propagação unitária → equivalência → **Gauss–Jordan incremental completo**. CryptoMiniSat faz Gauss durante a busca desde ~2009. Explica o *speedup* ≈ 1,00 do `xor_clauses_benchmark`: sem paridade densa, ganho nulo é o esperado. |

---

## 3. "Aplicaram ao cavalo, ou só a grafos gerais?" — ao cavalo

Chalaturnyk 2008 §11.2.1:

> "*The **first, and most extensive part**, demonstrates the performance of the
> new algorithm versus the old using a set of graphs known as **knight's
> graphs** … As mentioned in [Van98], graphs of this type are suitable for
> comparative analysis of Hamilton cycle algorithms because they tend to be
> **difficult to fully explore**.*"

Tabela 11.1 + 11.2 — **busca exaustiva**, contagem completa:

| Grafo | n | e | ciclos hamiltonianos | $t_{desc}$ (s) | $t_{old}$ = Kocay 92 (s) |
|---|---|---|---|---|---|
| kn4x8 | 32 | 64 | 0 | 0,000077 | 0,00015 |
| kn5x6 | 30 | 62 | 8 | 0,000026 | 0,00010 |
| **kn6x6** | **36** | **80** | **9 862** | **0,0077** | 0,022 |
| kn5x8 | 40 | 90 | 44 202 | 0,026 | 0,087 |
| kn6x7 | 42 | 98 | 1 067 638 | 0,74 | 2,42 |
| kn7x4 | 24 | 86 | 207 360 000 | 61,18 | 943,01 |

O **9 862** do nosso 6×6 já era benchmark em 2008, pelo mesmo algoritmo, em
**7,7 ms** — contra 53 s do `residual_search` em Python. O `kn5x8` = 44 202 é o
mesmo 5×8 do `missing_directions_5x8`.

**Uso construtivo:** esta tabela é *ground truth* independente e publicado para
4×8, 5×6, 5×8, 6×6, 6×7 e 7×4 — vale um teste de regressão do nosso enumerador
contra números de terceiros.

---

## 4. Relacionados que ficaram em `../referencia/`

Não movidos para não quebrar `../notes/citation_audit.md` e os `.tex`:

- `Parberry1997_efficient_algorithm_knights_tour.pdf` — D&C construtivo.
- `Marateck2008_warnsdorff_how_good.pdf`, `CancelaMordecki2006/2015` — Warnsdorff.
- `Pettersson2014_enumerating_hamiltonian_cycles.pdf`, `Jacobsen2007`, `Minato2001`,
  `LobbingWegener1996`, `McKay1997` — contagem por DP / matriz de transferência / BDD.
- `Knuth1975_estimating_backtrack.pdf` — o estimador.
- `Cygan2011_cut_and_count.pdf`, `Bodlaender2015`, `Ito2025` — conectividade parametrizada
  e reconfiguração.
- `Heinig2013`, `CNP2026`, `HamGen2026`, `HouYin2025`, `HefetzKrivelevich2025a/b` —
  **a linha teórica onde $Q(n)=3$ é contribuição.**

---

## 5. BibTeX

```bibtex
@article{Rubin1974,
  author  = {Frank Rubin},
  title   = {A Search Procedure for {H}amilton Paths and Circuits},
  journal = {Journal of the ACM}, volume = {21}, number = {4},
  pages   = {576--580}, year = {1974}, doi = {10.1145/321850.321854}}

@article{Kocay1992,
  author  = {William Kocay},
  title   = {An extension of the multi-path algorithm for finding {H}amilton cycles},
  journal = {Discrete Mathematics}, volume = {101}, number = {1-3},
  pages   = {171--188}, year = {1992}, doi = {10.1016/0012-365X(92)90600-K}}

@mastersthesis{Chalaturnyk2008,
  author = {Andrew Chalaturnyk},
  title  = {A Fast Algorithm for Finding {H}amilton Cycles},
  school = {University of Manitoba}, year = {2008}}

@mastersthesis{Vandegriend1998,
  author = {Basil Vandegriend},
  title  = {Finding {H}amiltonian Cycles: Algorithms, Graphs and Performance},
  school = {University of Alberta}, year = {1998}}

@article{VandegriendCulberson1998,
  author  = {Basil Vandegriend and Joseph Culberson},
  title   = {The $G_{n,m}$ Phase Transition is Not Hard for the {H}amiltonian Cycle Problem},
  journal = {Journal of Artificial Intelligence Research}, volume = {9},
  pages   = {219--245}, year = {1998}, doi = {10.1613/jair.512}}

@article{MatetiDeo1976,
  author  = {Prabhaker Mateti and Narsingh Deo},
  title   = {On Algorithms for Enumerating All Circuits of a Graph},
  journal = {SIAM Journal on Computing}, volume = {5}, number = {1},
  pages   = {90--99}, year = {1976}, doi = {10.1137/0205007}}

@misc{Sheffield2000,
  author = {Scott Sheffield},
  title  = {Computing and Sampling Restricted Vertex Degree Subgraphs and {H}amiltonian Cycles},
  year   = {2000}, eprint = {math/0008231}, archivePrefix = {arXiv}}

@inproceedings{CaseauLaburthe1997,
  author = {Yves Caseau and Fran\c{c}ois Laburthe},
  title  = {Solving Small {TSP}s with Constraints},
  booktitle = {ICLP}, year = {1997}}

@article{ZhangGuoChen1988,
  author  = {Fuji Zhang and Xiaofeng Guo and Rongsi Chen},
  title   = {{Z}-transformation graphs of perfect matchings of hexagonal systems},
  journal = {Discrete Mathematics}, volume = {72}, pages = {405--415}, year = {1988}}
```
