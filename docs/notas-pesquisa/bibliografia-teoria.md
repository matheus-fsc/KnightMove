# Biblioteca local de referências

Cópias locais dos trabalhos relacionados a `knight_tour_complete_v2.tex`.
Montado em 2026-08-06 via arXiv, Unpaywall e páginas dos autores.

**Somente material de acesso aberto.** Nada aqui contorna paywall — os itens
pagos estão listados na §2 com o caminho de acesso institucional.

- Auditoria de citações: `../notes/citation_audit.md`
- Plano de *tightness*: `../notes/tightness_plan.md`

---

## 1. Disponível localmente (23 PDFs, 12,6 MB)

### 1.1 Espaço de Hamilton — o núcleo teórico ⭐
A linha que enquadra o seu `deficit`. Leia nesta ordem.

| Arquivo | Trabalho | Por que importa |
|---|---|---|
| `Heinig2013_hamilton_circuits_cycle_space.pdf` | Heinig, arXiv:1303.0026 | **Observa que $\delta \geq 3$ é necessário** — é a fonte do fato de que `deficit > 0` é prior art. Cite no Teorema Q(n)=3. |
| `CNP2026_hamilton_space_pseudorandom.pdf` | Christoph–Nenadov–Petrova, JCTB 176 (2026) | ⭐ **O paper central do plano.** Lema 2.1 (dualidade), Def. 2.2 (parity-switcher), Lema 2.3 (Cycle Lemma). |
| `HamGen2026_cycle_space_spanned.pdf` | Hefetz–Krivelevich, arXiv:2606.05835 | Melhor survey da área na introdução; bibliografia completa. |
| `HouYin2025_dirac_hamilton_generated.pdf` | Hou–Yin, arXiv:2503.15950 | Resolve o caso Dirac; cunha o termo *Hamilton-generated*. |
| `HefetzKrivelevich2025a_random_graphs.pdf` | arXiv:2506.19731 | Estado da arte, grafos aleatórios. |
| `HefetzKrivelevich2025b_random_regular.pdf` | arXiv:2507.04488 | Regulares e perturbados. |

### 1.2 Superfícies — prior art da sua §6
| Arquivo | Trabalho |
|---|---|
| `ForrestTeehan2015_topology_knights_surfaces.pdf` | Forrest–Teehan, arXiv:1507.02917. Cilindro e toro; Teoremas 4.1, 5.1, 4.2, 6.1 (identidade vs. gerador de $\pi_1$). |
| `ForrestLague2024_nullhomotopic_nonorientable.pdf` | Forrest–Lague, arXiv:2406.05226. Möbius e Klein. **Ainda não comparado com os seus achados de Klein — pendência aberta.** |

### 1.3 Contagem e enumeração
| Arquivo | Trabalho |
|---|---|
| `KyekParberryWegener1997_bounds_knights_tours.pdf` | Bounds provados: Thm 4 $\Omega(1{,}3535^{n^2})$, Rmk 6 $\leq 4^{n^2}$, Thm 7 $N(8) \leq 3{,}019\times10^{22}$. Bibliografia dá Conrad, Kraitchik, Cull. |
| `LobbingWegener1996_BDD_count.pdf` | A contagem BDD **incorreta** do 8×8. |
| `McKay1997_knights_tours_8x8.pdf` | A correção: $N(8) = 13\,267\,364\,410\,532$. |
| `Knuth1975_estimating_backtrack.pdf` | Referência **primária** do estimador (não o TAOCP 4B). |
| `CancelaMordecki2006_counting_warnsdorff.pdf` | Prior art da sua §5: importance sampling via Warnsdorff randomizado. |
| `CancelaMordecki2015_open_knights_tours.pdf` | Continuação, tours abertos. |
| `Jacobsen2007_exact_enumeration_hamiltonian.pdf` | Matriz de transferência 2D/3D — referência da sua §9. |
| `Pettersson2014_enumerating_hamiltonian_cycles.pdf` | Enumeração de ciclos hamiltonianos. |
| `Minato2001_zero_suppressed_bdds.pdf` | ZDDs; contextualiza o estouro do seu C++ em $n=8$. |

### 1.4 Algoritmos e heurísticas
| Arquivo | Trabalho |
|---|---|
| `Parberry1997_efficient_algorithm_knights_tour.pdf` | **D&C para tours** — prior art da sua §3.2. Ele *constrói* um tour; você *conta* uma subclasse. |
| `Parberry_knights_tour_3d.pdf` | Bônus: versão 3D. |
| `Marateck2008_warnsdorff_how_good.pdf` | Análise da taxa de sucesso de Warnsdorff. |

### 1.5 Conectividade parametrizada
| Arquivo | Trabalho |
|---|---|
| `Cygan2011_cut_and_count.pdf` | ⭐ **Cut & Count** — conta mod 2 para que desconexas se cancelem. A rota conhecida de contorno da barreira que você identificou. |
| `Bodlaender2015_deterministic_connectivity.pdf` | Versão determinística. |
| `Ito2025_rerouting_planar_curves.pdf` | Reconfiguração de caminhos (ACM TALG 21, 2025). |

---

## 2. Não disponível — precisa de acesso institucional

Rota: **Portal de Periódicos CAPES** com login federado **CAFe** usando
credenciais UNIFEI (funciona fora do campus; eduroam só dá acesso por IP
dentro do campus).

| Trabalho | Venue | Prioridade |
|---|---|---|
| **Conrad et al. 1994**, Solution of the knight's Hamiltonian path problem | Discrete Appl. Math. 50, 125–134 | 🔴 **Alta** — resolve o passo (S3) do plano de tightness. DOI `10.1016/0166-218X(92)00170-Q` |
| **Schwenk 1991**, Which Rectangular Chessboards Have a Knight's Tour? | Math. Magazine 64(5) / JSTOR | 🔴 Alta — você cita e o enunciado estava errado na v1 |
| **Hartman 1983**, Long cycles generate the cycle space | European J. Combin. 4, 237–246 | 🟠 Média — ancestral da linha Hamilton space |
| **Heinig 2014**, On prisms, Möbius ladders… | European J. Combin. 36, 503–530 | 🟠 Média |
| **Watkins–Hoenigman 1997**, Knight's Tour on a Torus | Math. Magazine 70(3) / JSTOR | 🟠 Média — prior art de superfícies |
| **Watkins 2004**, *Across the Board*, pp. 65–77 | Princeton UP (livro) | 🟠 Média — biblioteca |
| **Itai–Papadimitriou–Szwarcfiter 1982**, Hamilton Paths in Grid Graphs | SIAM J. Comput. 11(4) | 🟡 Baixa — só NP-completude |
| **Chia–Ong 2005**, Generalized knight's tours | Discrete Appl. Math. 150 | 🟡 Baixa |
| **Pohl 1967** | CACM 10(7) | 🟡 Baixa |
| **Chen 1992**, Heuristic Sampling | SIAM J. Comput. 21(2) | 🟡 Baixa — redução de variância |
| **Yen 1971**, **Eppstein 1998** | Manag. Sci. / SIAM J. Comput. | 🟡 Baixa — já citados |
| **Cull–De Curtins 1978** | Fibonacci Quarterly 16 | 🟢 **Aberto** em `fq.math.ca` — só não automatizei |
| **MacLane 1937** | Fundamenta Math. 28 | 🟢 Aberto em EuDML/matwbn — o link direto falhou (403) |
| **Kraitchik 1953**, *Mathematical Recreations* | Dover (livro) | 🟢 Tentar `archive.org` |
| **Diestel 2017**, *Graph Theory* | Springer | 🟢 **Grátis** em `diestel-graph-theory.com` |
| **Alspach–Locke–Witte** | ref [2] do CNP2026 | ❓ **Referência completa ainda não obtida** — extrair da bibliografia do CNP |

---

## 3. Como reproduzir / estender

Fluxo que funcionou: **OpenAlex** (descobre) → **DBLP** (acha versão
publicada) → **Unpaywall** (acha PDF legal) → CAPES só no que sobrar.
Nenhum precisa de API key; Unpaywall e OpenAlex só pedem `email=`/`mailto=`.

```bash
# versão publicada de um preprint
curl -s "https://dblp.org/search/publ/api?q=TITULO&format=json" | python3 -m json.tool

# todas as cópias OA de um DOI (inclui repositórios, não só o publisher)
curl -s "https://api.unpaywall.org/v2/DOI?email=SEU@EMAIL" | python3 -m json.tool

# BibTeX canônico
curl -sLH "Accept: application/x-bibtex" "https://doi.org/DOI"
```

⚠️ Duas lições da montagem desta pasta:
1. **Não confie no `best_oa_location`** — para Jacobsen e Bodlaender ele
   apontava para o publisher (que devolve HTML ou 403). A lista completa
   `oa_locations` tinha os repositórios (arXiv, TU/e) que funcionaram.
2. **Verifique o que baixou.** Dois arXiv IDs que chutei por busca vieram
   com artigos completamente diferentes (um de física de férmions). Sempre
   confira: `pdftotext -f 1 -l 1 arquivo.pdf - | head -3`.

---

## 4. Prior art **algorítmico** → pasta separada

O levantamento sobre `knight_tours.py` (multi-path / Rubin 1974 / Kocay 1992) e
sobre o **XOR de loops** (Mateti–Deo 1976, Z-transformation de Sheffield 2000,
DPLL(XOR)) fica em **[`../referencia_algoritmo/`](../referencia_algoritmo/README.md)**
— 9 itens, com tabela de correspondência linha-a-linha e BibTeX.

**TL;DR:** nem o algoritmo nem o XOR são novos, e ambos já foram aplicados ao
grafo do cavalo (Chalaturnyk 2008 enumera o 6×6 = 9 862 tours em 7,7 ms).
A contribuição defensável é a linha teórica da §1.1 desta pasta.
