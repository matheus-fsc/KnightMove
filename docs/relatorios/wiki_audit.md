# Auditoria do Projeto Knight Tour

> Saída da **Tarefa 0** de `agent_prompt_wiki.md`.
> Gerado em 2026-05-20.
> Nenhum arquivo de wiki foi criado — aguarda aprovação.

---

## 1. Site existente

### Tecnologia
- **HTML puro + ES Modules + Canvas 2D.** Sem build, sem framework.
- Estilos: um único `css/styles.css` (7.5 KB).
- Lógica: módulos ES (`type="module"`) em `js/` + `js/modes/`.
  Pacote total ≈ 1.660 linhas de JS, sem dependências externas.
- Sem `package.json`, `mkdocs.yml`, `_config.yml`, `next.config.*`.

### Arquivos do site
```
/
├── index.html              ← landing page (links para os 2 visualizadores)
├── cavalo_viz.html         ← visualizador principal (8×8, dados do JSON)
├── unicornios_viz.html     ← visualizador secundário (6×6, "unicórnios topológicos")
├── cavalo_data.json        ← dados do 8×8 consumidos por cavalo_viz
├── assinaturas_6x6.json    ← dados do 6×6 consumidos por unicornios_viz
├── css/
│   └── styles.css          ← único stylesheet
└── js/
    ├── app.js              ← orquestração: carrega DATA, monta board + tree, gerencia modos
    ├── core.js             ← loadData, edgeKey, getLoopEdges
    ├── board.js            ← canvas do tabuleiro
    ├── loop-tree.js        ← canvas da árvore de loops
    └── modes/              ← 7 modos do visualizador
        ├── grafo.js        (47 linhas) — exibe o grafo do cavalo
        ├── loops.js        (143)      — ciclos fundamentais
        ├── solucao.js      (141)      — tour hamiltoniano heurístico
        ├── arvore.js       (127)      — árvore geradora BFS
        ├── hierarquia.js   (145)      — hierarquia de loops
        ├── simetria.js     (133)      — órbitas D₄ + arestas obrigatórias
        └── xor.js          (286)      — restrições XOR/cláusulas
```

### O que o visualizador faz
1. `cavalo_viz.html` é o painel principal. Layout em 3 colunas:
   - **Coluna 1:** canvas 560×560 com o tabuleiro do cavalo + controles
     de velocidade.
   - **Coluna 2:** painel de estatísticas globais (V, E, H¹, obrigatórias,
     impossíveis, loops determinados, loops livres) + painel do modo ativo.
   - **Coluna 3:** canvas 420×460 com a árvore de loops + contador.
2. Modos selecionáveis: grafo, loops, solução, árvore BFS, hierarquia,
   XOR, simetria. Cada modo é um módulo ES com hooks `enter()`/`exit()`/
   `draw(board)`.
3. `unicornios_viz.html` é um one-pager menor: visualização das caudas
   da distribuição de assinaturas booleanas no 6×6.
4. Tudo é estático — `cavalo_data.json` (138 KB) e `assinaturas_6x6.json`
   (375 KB) são carregados via `fetch`, processados client-side.

### Como o site é servido
- **Local.** Comando recomendado pelo próprio `app.js`:
  `python3 -m http.server` (precisa ser HTTP, não `file://`).
- Sem deploy configurado (sem `gh-pages`, sem CI). É um site estático
  que funciona em qualquer host HTTP.

---

## 2. READMEs encontrados (18 arquivos `.md`)

| # | Caminho | Módulo | Status | Resultado principal |
|---|---------|--------|--------|---------------------|
| 1 | `README.md` | raiz | desatualizado | descrição básica + pipeline 8×8 antigo |
| 2 | `deficit_theorem/README.md` | deficit | **completo** | `Q(n)=3` provado p/ todo `n≥4`; deficit=3 verificado p/ `n∈{6,8,10,12}` |
| 3 | `deficit_theorem/conjecture_proof.md` | deficit | **completo** | prova formal Lemas 1–4 + Teorema |
| 4 | `forbidden_cycles_6x6/README.md` | forbidden | **completo** | 3 ciclos proibidos identificados; `φ₁: x_{F6-D5}⊕x_{B3-A1}=0` (detector mínimo) |
| 5 | `6x6_higher_order/README.md` | higher-order 6×6 | **completo** | 88 pares + 1.776 triplas + 17.004 quadras minimais; `rank(T)=42<β₁=45` |
| 6 | `board_10x10/README.md` | board 10×10 | **completo** | `V=100, E=288, β₁=189`; rank=186; 8 obrigatórias; bifurcação top-1 `B10-D9 ↔ B10-C8` |
| 7 | `board_10x10/higher_order_10x10/README.md` | h.o. 10×10 | **completo (negativo)** | 0 triplas genuínas; `NOT(A∧B∧C)` gera slowdown 0.71× — transição qualitativa 6×6→10×10 |
| 8 | `xor_clauses_benchmark/README.md` | XOR bench | **completo** | speedup ≈ 1.00×; deficit 10×10 = 3 confirmado |
| 9 | `residual_search/README.md` | residual 6×6 | **completo** | R1..R6 fixa 0 livres; BT+R1..R6 → 9.862 tours em 53s; 3.560 2-fatores extras |
| 10 | `residual_search_10x10/README.md` | residual 10×10 | **completo** | BT+pressão acha K=50 em 1.63s (2.04× vs Z3); só 5.7% dos 2-fatores são conexos |
| 11 | `incremental_subtour/README.md` | UF detector | **completo** | BT v2: 25× vs Z3, 5.5× vs v1 no 10×10; razão 2-fat/tour 18.1×→1.00× |
| 12 | `local_phase_heuristic/README.md` | f∞(L) | **completo** | tabela teórica de 6 valores; ≤1.14× v2 amostrado; 623× mais rápido em t_total_justo |
| 13 | `martingale_analysis/README.md` | martingale | **completo** | H_MART forte refutada; versão fraca confirmada (`μ_min ≥ 0.173`) |
| 14 | `transfer_matrix/README.md` | transfer matrix | parcial | n=6 OK (count=36.236, λ₁=70.48); n=8 bloqueado por RAM |
| 15 | `complex_orbit/README.md` | complex orbit | **completo** | DFT/winding/CTQW; rank(T)=42 sobre GF(2); 18 invariantes geométricos |
| 16 | `patch_lut_results.md` | patch LUT | **completo (negativo)** | LUT correta mas taxa de resgate = 0% (R2+UF funde componentes cedo) |
| 17 | `agent_prompt_8x8.md` | prompt arquivado | histórico | brief antigo — ignorar para wiki |
| 18 | `agent_prompt_benchmark.md` | prompt arquivado | histórico | brief antigo — ignorar para wiki |
| 19 | `cavalo_8x8/agent_prompt_fase2.md` | prompt arquivado | histórico | brief antigo — ignorar para wiki |
| 20 | `agent_prompt_wiki.md` | prompt atual | — | este brief |

### Dependências entre módulos
```
deficit_theorem ────────► usa rank/β₁ medidos em forbidden_cycles_6x6 e board_10x10
forbidden_cycles_6x6 ───► usa enumeração exaustiva 9.862 tours (complex_orbit/data/tours_closed_6x6.npy)
6x6_higher_order ───────► usa mesma incidência aresta×tour
board_10x10 ────────────► amostragem SAT (Z3); base p/ higher_order_10x10, xor_clauses_benchmark, residual_search_10x10
incremental_subtour ────► motor v2; base p/ local_phase_heuristic, martingale_analysis, knight_tours.py
local_phase_heuristic ──► motor "theory" autocontido; codificado em knight_tours.py (F_INF)
knight_tours.py ────────► destino final do pipeline: standalone, sem deps, sem amostragem
tour_count_estimator.py ► importa knight_tours.py
patch_lut_results ──────► investigação fechada (negativa)
transfer_matrix ────────► explora autovetor de T; conexão com f∞ confirmada na marginal exata
complex_orbit ──────────► linha lateral (ℂ, DFT, CTQW); confirma rank(T)=42
```

---

## 3. Plots disponíveis (todos `.png`)

Total: **38 plots** em 14 diretórios.

| Arquivo | Localização | Assunto |
|---------|-------------|---------|
| `ratio_analysis.png` | `data/plots/` | r(n) = tours/2-fatores, decaimento exponencial |
| `tour_count_estimates.png` | `data/plots/` | N(n) por Knuth, IC95% |
| `forbidden_cycle_0.png` | `forbidden_cycles_6x6/data/plots/` | hexágono NW (f₀) |
| `forbidden_cycle_1.png` | `forbidden_cycles_6x6/data/plots/` | hexágono NE (f₁) |
| `forbidden_cycle_2.png` | `forbidden_cycles_6x6/data/plots/` | octógono diagonal (f₂) |
| `ideal_structure.png` | `6x6_higher_order/data/plots/` | histogramas H₁ por ordem |
| `vertex_coverage.png` | `6x6_higher_order/data/plots/` | heatmap 6×6 de cobertura |
| `propagation_cascade.png` | `residual_search/data/plots/` | cascata R1..R6 no 6×6 |
| `residual_structure.png` | `residual_search/data/plots/` | grafo residual e correlações |
| `propagation_cascade_10x10.png` | `residual_search_10x10/data/plots/` | cascata 10×10 |
| `residual_structure_10x10.png` | `residual_search_10x10/data/plots/` | grafo residual 10×10 |
| `benchmark_10x10.png` | `residual_search_10x10/data/plots/` | BT vs Z3 (n=10) |
| `speedup_v2.png` | `incremental_subtour/data/plots/` | barras de speedup BT v2 |
| `poda_breakdown.png` | `incremental_subtour/data/plots/` | eventos v1 vs v2 |
| `subtour_size_hist.png` | `incremental_subtour/data/plots/` | histograma de tamanhos |
| `subtour_depth_hist.png` | `incremental_subtour/data/plots/` | histograma de profundidades |
| `scaling_full.png` | `incremental_subtour/data/plots/` | escala n∈{6..14} |
| `edge_freq_convergence.png` | `incremental_subtour/data/plots/` | convergência de f∞ por anel |
| `nodes_per_tour_comparison.png` | `local_phase_heuristic/data/plots/` | A/B/C/D barras |
| `scaling_theory.png` | `local_phase_heuristic/data/plots/` | theory vs v2 em n∈{6..14} |
| `sensitivity_heatmap.png` | `local_phase_heuristic/data/plots/` | sensibilidade por nível |
| `p_distribution_by_n.png` | `martingale_analysis/data/plots/` | p(b) × n com IC Wilson |
| `mu_decomposition.png` | `martingale_analysis/data/plots/` | μ_total / interior / borda |
| `chi2_heatmap.png` | `martingale_analysis/data/plots/` | -log₁₀(p) por bin |
| `convergence_evidence.png` | `martingale_analysis/data/plots/` | 4 testes formais |
| `martingale_process.png` | `martingale_analysis/data/plots/` | X_t, M_t, Var, QQ |
| `spectral_gap.png` | `transfer_matrix/data/plots/` | gap λ₁ − λ₂ |
| `eigenvector_vs_finf.png` | `transfer_matrix/data/plots/` | autovetor vs f∞ |
| `tour_count_scaling.png` | `transfer_matrix/data/plots/` | 3 modelos de calibração |
| `speedup_comparison.png` | `xor_clauses_benchmark/data/plots/` | speedups XOR 6×6/10×10 |
| `benchmark_comparison.png` | `board_10x10/higher_order_10x10/data/plots/` | benchmark NOT(A∧B∧C) |
| `fourier_mean_spectrum.png` | `complex_orbit/data/plots/` | espectro DFT médio dos 9.862 tours |
| `invariant_classification.png` | `complex_orbit/data/plots/` | heatmaps 6×6 de invariantes |
| `quantum_walk_snapshots.png` | `complex_orbit/data/plots/` | CTQW em t=0,10,36,100 |
| `quantum_walk_uniform_distance.png` | `complex_orbit/data/plots/` | TV vs t |
| `winding_linearity_basis_free.png` | `complex_orbit/data/plots/` | linearidade GF(2) winding |
| `winding_mean_heatmap.png` | `complex_orbit/data/plots/` | mean \|w(p)\| 6×6 |

Plots **não disponíveis** mencionados no brief:
- `deficit_theorem/data/plots/` está vazio.
- `board_10x10/data/plots/` está vazio.

### JSONs de resultado relevantes (não dados brutos)
- `data/ratio_main_results.json`, `data/tour_count_estimates.json`, `data/benchmark_8x8.json`
- `forbidden_cycles_6x6/data/results/forbidden_cycles.{npy,json}`, `dual_detectors.{npy,json}`
- `6x6_higher_order/data/exclusions_{pairs,triples,quads}.json`, `ideal_summary.json`
- `deficit_theorem/data/results/verify_n{6,8,10,12}.json`, `structural_n{4..12}.json`
- `board_10x10/data/invariants/{summary,candidate_clauses,mandatory_edges}.json`
- `incremental_subtour/data/{benchmark_results,subtour_analysis,scaling_minimal_v2}.json`
- `local_phase_heuristic/data/{benchmark_theory,calibration_results,benchmark_comparison}.json`
- `martingale_analysis/data/{p_distribution,convergence_tests}.json`
- `transfer_matrix/data/{eigenvalues_n6,exact_2factor_marginals_n6,state_count}.json`
- `xor_clauses_benchmark/data/{benchmark_6x6_xor,benchmark_10x10_xor,xor_10x10_clauses}.json`
- `complex_orbit/data/{tours_closed_6x6.npy,fourier_coefficients.npy,winding_numbers.npy,invariant_classification.json}`
- `analysis_6x6.json`, `summary.json`

---

## 4. Entry points de código

### Standalone (sem amostragem, sem solver, sem dados externos)
- **`knight_tours.py`** (raiz, 472 linhas, dep: numpy). Funções públicas:
  - `knight_tours(n, K, seed=None) -> list[np.ndarray]`
  - `verify_tour(tour, n) -> bool`
  - Internals reutilizáveis: `build_graph`, `State`, `_propagate_initial`,
    `_choose_next_edge`, `fix_and_propagate`, `_is_complete_tour`.
  - Constantes: `F_INF = {0:0.528, 1:0.193, 2:0.198, 3:0.294, 4:0.247, 5:0.261}`,
    `F_INF_DEFAULT = 0.25`.

- **`tour_count_estimator.py`** (raiz, ~10 KB, dep: numpy + matplotlib opcional).
  Importa de `knight_tours`. Funções públicas:
  - `random_path_weight(n, ctx, rng) -> float`
  - `estimate_tour_count(n, M, seed, ...)`
  - `validate_estimator()` (valida vs `n=6` exato e `n=8` literatura)
  - `run_all_estimates(M6, M8, M10, M12, seed)`

### Pacote `knight_tours_optimized/`
- `__init__.py` exporta:
  - de `tours.py`: `knight_tours`, `verify_tour` (mesma API)
  - de `paths.py`: `knight_path`, `verify_path` (caminhos abertos)
  - de `symmetry.py`: `d4_symmetries`, `canonical_form`, `is_canonical`,
    `expand_d4`, `verify_d4_decomposition`, `knight_tours_canonical`
  - de `parallel.py`: `knight_tours_parallel`, `knight_path_parallel`
- Opcionais (require numba):
  - `core_numba.py`: `knight_tours_numba`, `NUMBA_AVAILABLE`
  - `parallel_prefix.py`: `generate_prefixes`, `knight_tours_prefix_parallel`,
    `full_benchmark_optimized`
- Outros: `core.py`, `benchmark.py`, `cli.py`, `divide_conquer.py`, `tests.py`.

### Como um usuário novo rodaria o código do zero
```bash
git clone <repo>
cd knight_tour
pip install numpy      # única dependência obrigatória
# (opcional) pip install matplotlib z3-solver numba
python3 -c "from knight_tours import knight_tours, verify_tour; \
            tours = knight_tours(n=10, K=10, seed=42); \
            print(len(tours), verify_tour(tours[0], 10))"
```

### Scripts auxiliares na raiz (não fazem parte da API pública)
- `cavalo_engine.py`, `cavalo_engine_z3`, `cavalo_engine_8x8_path*.py`:
  geradores legados que produzem `cavalo_data.json` para o visualizador.
- `analise_6x6.py`, `b6_grouping_6x6.py`, `b6_interior_correlations_6x6.py`:
  análises do 6×6 que alimentam `analysis_6x6.json`, `b6_*.json`.
- `cavalo_loop_destruicao_6x6.py`: backtracking exaustivo + destruição H₁.
- `knight_tours_patch.py`, `knight_tours_deformation_gf2.py`: variantes
  experimentais (patch LUT — encerrada negativamente).
- `ratio_analysis.py`, `run_ratio_experiment.py`, `ratio_analyze_results.py`:
  experimento r(n) = tours/2-fatores.
- `benchmark/` (sub-dir): benchmarks históricos (`baseline_bt_6x6.py`,
  `benchmark_z3.py`, `gt_paths_6x6.py`, `kl_divergence.py`, `report.py`).

### Relatório acadêmico já existente
- `relatorio_professor.tex` (37 KB) + PDF (1.1 MB) compilam um documento
  formal com Teoremas e Conjecturas, cobrindo o mesmo material que a
  wiki vai apresentar — fonte cruzada útil para validar números.

---

## 5. Estrutura de wiki proposta

Mantenho o esqueleto sugerido pelo brief (11 páginas). Justificativa:
cada página tem material concreto (números + plots) e não há sobreposição
significativa entre temas.

```
wiki/
├── index.html              ← landing: mapa de resultados + 3 destaques + timeline
├── style.css               ← classes com prefixo .w- (não conflita com css/styles.css)
├── nav.js                  ← sidebar dinâmica
├── 01-introducao.html      ← problema, contexto histórico, escopo
├── 02-teoria-gf2.html      ← GF(2), Z₁, β₁, deficit
├── 03-teorema-q3.html      ← cantos, 4 representantes, σ, 3 ciclos proibidos
├── 04-algoritmo.html       ← BT + R2 + UF + f∞; API knight_tours()
├── 05-fase-local.html      ← f∞(L), tabela universal, sensibilidade
├── 06-benchmark.html       ← Z3 vs BT v1 vs BT v2 vs theory; 48× no 10×10
├── 07-estimativas.html     ← Knuth, validação, N(10)≈2.4×10²², N(12)≈1.3×10³³
├── 08-topologia.html       ← r(n) decai exp; transição de fase n=10→12
├── 09-extensoes.html       ← D₄ simetria; caminhos abertos; divide-and-conquer
├── 10-resultados.html      ← tabela completa proven/empirical/open/negative
└── 11-codigo.html          ← instalação, API, estrutura do projeto
```

Mapeamento README → página wiki:
| Página | Fontes READMEs |
|--------|----------------|
| 02 | `deficit_theorem/`, `forbidden_cycles_6x6/`, `6x6_higher_order/` |
| 03 | `deficit_theorem/conjecture_proof.md`, `forbidden_cycles_6x6/` |
| 04 | `incremental_subtour/`, `residual_search/`, `local_phase_heuristic/`, `knight_tours.py` |
| 05 | `local_phase_heuristic/`, `incremental_subtour/` |
| 06 | `xor_clauses_benchmark/`, `incremental_subtour/`, `residual_search_10x10/`, `board_10x10/` |
| 07 | `tour_count_estimator.py`, `data/tour_count_estimates.json`, `transfer_matrix/` |
| 08 | `ratio_analysis.py`, `data/ratio_main_results.json`, `martingale_analysis/` |
| 09 | `knight_tours_optimized/symmetry.py`, `.../paths.py`, `.../divide_conquer.py`, `complex_orbit/` |
| 10 | consolidação de todas as fontes |
| 11 | `knight_tours.py`, `knight_tours_optimized/__init__.py`, `tour_count_estimator.py` |

---

## 6. Pontos de atenção

### 6.1 Integração com o site existente
- **Não há conflito de CSS** desde que o brief seja seguido (prefixo `.w-`).
  O `css/styles.css` atual usa seletores como `body`, `#board-col`,
  `#info-col`, `.gs-item`, `#tree-col` — nada com prefixo `w-`.
- **Não há conflito de JS** desde que a wiki seja servida em arquivos
  separados (`wiki/nav.js` não importa nada de `js/app.js`).
- **MathJax e highlight.js via CDN** não interferem nos visualizadores
  (eles não usam essas libs).

### 6.2 Como integrar (Tarefa 8)
- O `index.html` atual da raiz é uma landing com 2 links. Adicionar
  um terceiro link `<a href="./wiki/index.html">📖 Wiki</a>` é trivial
  e mantém o estilo dark já existente.
- **Não modificar `cavalo_viz.html`, `unicornios_viz.html`, `js/*` ou
  `css/styles.css`.**

### 6.3 Dados disponíveis vs faltantes para a wiki
- ✅ `forbidden_cycle_{0,1,2}.png` existem — podem ser embebidos diretos.
- ✅ Plots de benchmark, scaling, martingale, transfer matrix existem.
- ⚠️ `deficit_theorem/data/plots/` está vazio. A página 03 pode usar
  diagrama ASCII ou referenciar `forbidden_cycle_*.png` em vez disso.
- ⚠️ `board_10x10/data/plots/` está vazio. A página 06 usa
  `xor_clauses_benchmark/.../speedup_comparison.png` e
  `residual_search_10x10/.../benchmark_10x10.png` no lugar.
- ⚠️ Não há plot pronto para a "transição de fase" do `ratio_analysis`
  (`ratio_analysis.png` em `data/plots/` cobre o decaimento exponencial,
  mas não o cruzamento com Q(n)=3). Página 08 mostra a tabela e o plot
  existente.

### 6.4 Conflitos numéricos entre fontes
- `forbidden_cycles_6x6/README.md` lista `rank(2-fatores)=42` (mesma do
  rank dos tours). O brief de wiki menciona "razão tours/2-fatores =
  N_t/N_{2f}". Distinguir: **rank de incidência GF(2)** vs **contagem
  combinatória**. Os dois números coexistem.
- `tour_count_estimates.json` reporta `N(10)=2.4×10²²` com IC95% bem
  largo (`[5.6×10²⁰, 6.8×10²²]`). A página 07 deve mostrar o IC
  honestamente, não só o ponto-estimado.
- `local_phase_heuristic` reporta tempo "623× mais rápido em t_total_justo"
  que inclui amostragem de 417 s; o `incremental_subtour` reporta 25.2×
  vs Z3 no benchmark direto. Ambos são corretos para definições
  diferentes — a página 06 deve diferenciá-los.

### 6.5 Material já consolidado (cross-check)
- `relatorio_professor.tex` (já compilado) é uma versão acadêmica do
  mesmo conteúdo. Usar como referência cruzada para enunciados formais.
  Não embutir o PDF na wiki por enquanto (1.1 MB).

### 6.6 Item do brief que muda diante da auditoria
- O brief original lista o caminho `knight_tour/site/` — **não existe**.
  O site existente fica na raiz do repositório. A integração (Tarefa 8)
  deve usar `index.html` da raiz.
- O brief presume `knight_tour/README.md` informativo — o atual é
  desatualizado (descreve o pipeline 8×8 antigo). Não é fonte primária
  para a wiki, mas a página 01 deve substituir esse texto.

---

## 7. Próximo passo

→ **PARAR e aguardar aprovação** antes de criar qualquer arquivo dentro
de `wiki/`. Após o ok, executar Tarefa 1: criar esqueleto de `wiki/` com
arquivos `.html` vazios + `style.css` + `nav.js`.
