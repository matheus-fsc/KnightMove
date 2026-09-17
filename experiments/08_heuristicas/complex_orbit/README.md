# complex_orbit — Knight's Tour as a Complex Curve

Research log for a new branch of the Knight's Tour project. The aim is to
represent the knight movement on a 6×6 board using complex numbers and to
investigate what the real/imaginary decomposition reveals about the global
topology of closed tours.

This work plugs into the main project's existing findings on the 6×6 board:

- **9,862** closed (Hamiltonian) tours, computed exhaustively
- **1,232** D₄-canonical equivalence classes
- **8** mandatory edges
- **7** topological invariants under D₄ symmetry (GF(2) / H₁ homology basis)
- Conservation law: total loop destruction = H₁ = 45 across all directed solutions

(See the parent project at `../analysis_6x6.json`, `../cavalo_loop_destruicao_6x6.py`,
and `../b6_grouping_6x6.json`.)

---

## 1. Mathematical framing

### 1.1 Complex position representation

A board cell `(row, col)` is encoded as

    z = col + i·row   ∈   ℤ[i]   ⊂   ℂ

The eight knight moves are complex displacements

    Δz ∈ { ±1 ± 2i ,  ±2 ± i }

They factor as the orbit of `δ₀ = 1 + 2i` under the cyclic group
`⟨i⟩ ≅ ℤ/4`, together with the orbit of its complex conjugate `δ̄₀ = 1 − 2i`:

    { iᵏ · δ₀  : k = 0,1,2,3 }  ∪  { iᵏ · δ̄₀  : k = 0,1,2,3 }

Equivalently the knight move set is closed under the action of the
dihedral group **D₄**, acting on ℂ by `z ↦ iᵏ z` (rotations) and
`z ↦ z̄` (reflection). Every move satisfies

    |Δz|² = 5

so the move locus is the discrete circle of radius √5 around the current
position.

The minimal monic polynomial vanishing on all eight moves is

    P(z) = (z² − 1)² · (z² − 4)² − …   = z⁸ + 14 z⁴ + 625  (mod sign)

(Concretely: the eight points are roots of `(z⁴+a)(z⁴+b)` for a careful
choice of `a, b`; expanded → `z⁸ + 14 z⁴ + 625`.) Task A verifies this
numerically.

### 1.2 Tour as a closed discrete complex curve

A closed Hamiltonian tour `T` of `N = n² = 36` steps gives a sequence

    z(0), z(1), …, z(N−1), z(N) = z(0)

This is a closed polygon in ℂ. The piecewise-linear interpolation never
crosses any cell-center `(c + ½) + i (r + ½)` because every knight segment
goes between *integer* points by displacements `(±1, ±2)` or `(±2, ±1)`;
a short algebraic check shows that no such segment can hit a half-integer
point. This lets us define winding numbers around all 36 cell centers
without singularities.

### 1.3 DFT decomposition

For each closed tour we compute the Discrete Fourier Transform of the
position sequence:

    c_k  =  (1/N) · Σ_{t=0..N−1}   z(t) · exp( −2π i · k · t / N ),
    k = 0, …, N − 1

Interpretation of the spectral lines:

| k          | meaning                                              |
|------------|------------------------------------------------------|
| `k = 0`    | centroid of the tour                                 |
| `k = 1`    | fundamental "rotation" — winding around the centroid |
| `k` small  | global structure / large-scale geometry              |
| `k` large  | local zig-zag detail                                 |

Under the D₄ action on the board, a rotation by `i` sends `z ↦ i·z`,
which multiplies every `c_k` by `i` (a global phase, identical for every
`k`). A reflection `z ↦ z̄` sends `c_k ↦ c_{−k}` (i.e. `c_{N−k}`) and
conjugates. So magnitudes `|c_k|` are invariant under the full D₄ action
on the board *up to* the involution `k ↔ N − k`. Task B reports this.

### 1.4 Winding number

For each closed tour `T` and each cell-center `p = (c+½) + i(r+½)` of the
board, the discrete winding number is

    w(p)  =  (1 / 2π) · Σ_{t=0..N−1}   Im( log( (z(t+1) − p) / (z(t) − p) ) )

where the principal branch of `log` is used and the imaginary parts are
each in `(−π, π]`. This is the standard "signed lap counter" — the path
is closed so `w(p) ∈ ℤ`.

### 1.5 Dependency tree (t-indexed state tree)

At time `t` the state is `S(t) = (z(t), V_t)` with `V_t ⊆ board` the
set of visited squares. The branching factor at depth `t` is

    |{ Δz : z(t) + Δz ∈ board,  z(t) + Δz ∉ V_t }|

The tree is doubly non-Markovian: the valid moves depend on the *full
history* `V_t`, not just on `z(t)`. This is why the problem is hard.

---

## 2. Modules

| File                    | Task | Purpose                                                                     |
|-------------------------|------|-----------------------------------------------------------------------------|
| `complex_moves.py`      | A    | Core complex representation, 8 moves, factorization check, P(z) check.      |
| `fourier_tours.py`      | B    | DFT of all 9,862 closed tours, mean spectrum plot.                          |
| `winding_number.py`     | C    | Winding around all 36 cell centers, correlation with H₁ invariant bits.     |
| `quantum_walk.py`       | D    | Continuous-time quantum walk on the 6×6 knight graph.                       |
| `h1_invariants.py`      | §5   | Reconstrução da base GF(2) de H₁ dos tours; classificação geom × algébrico. |
| `_tours.py`             | —    | Helper: backtracking exaustivo + edge orbits + tour ↔ complex conversion.   |

Outputs go to `data/` (`.npy`, `.json`) and `data/plots/` (`.png`).

---

## 3. Findings

### Task A — complex moves (`complex_moves.py`)

All four formal checks pass:

| Check                                                       | Result |
|-------------------------------------------------------------|--------|
| `|Δz|² = 5` for all 8 moves                                 | ✓      |
| Move set = `⟨i⟩·δ₀  ∪  ⟨i⟩·δ̄₀`, with `δ₀ = 1 + 2i`         | ✓      |
| `P(z) = z⁸ + 14 z⁴ + 625` vanishes on every Δ               | ✓ (machine-zero) |
| `valid_moves` respects board bounds and visited set         | ✓      |

The factorisation makes the D₄ symmetry of the move set explicit:
`i^k · δ₀` are 4 of the 8 moves; the other 4 are `i^k · δ̄₀`.

### Task B — Fourier decomposition (`fourier_tours.py`)

- DFT computed for all 9,862 canonical closed tours, shape `(9862, 36, 2)` saved to `data/fourier_coefficients.npy`.
- `|c_0 − (2.5 + 2.5i)|_max = 0` exactly: every closed tour has centroid equal to the geometric centre of the 6×6 board, as it must (Hamiltonian → visits every cell once).
- Mean magnitudes across all tours order as: `k = 0 ≫ k = 3, 4 ≫ k = 5, 6 ≫ …`. By the symmetry `|c_{N−k}| = |c_k|` after reflection of the board, the modes `k = 32, 33` mirror `k = 3, 4`.
- The dominant non-trivial mode is `k = 3`, **not** `k = 1`. The "fundamental rotation" `|c_1|` averages to `0.4681` — *smaller* than the third harmonic. Interpretation: on average a closed knight tour does not look like a single loop around the centroid; it looks like a 3-lobed / 4-lobed pattern. This is consistent with the geometric fact that knight tours are visually "weavy".
- Spectral floor (smallest mean |c_k| for k > 0) occurs at `k = 18 = N/2` — the highest Nyquist-like frequency.
- Plot: `data/plots/fourier_mean_spectrum.png` (mean ± 1σ across tours).

### Task C — winding numbers vs. H₁ invariants (`winding_number.py`)

- All 9,862 winding maps computed and saved to `data/winding_numbers.npy` of shape `(9862, 6, 6)`.
- Per-cell winding range observed: `[−4, 5]`.
- **Geometric observation.** Choosing test points `p = (c+½, r+½)` for `r, c ∈ {0,…,5}` gives 36 candidates, but only the 25 cells with `r, c ∈ {0,…,4}` lie in the convex hull of the integer board. The 11 boundary points (last row, last column of the grid) have `w(p) = 0` for *every* tour — the path never winds around something outside its own convex hull.
- Mean `|w(p)|` is highest on the central cells (≈ 1.28 at the four off-centre cells around the middle) and drops smoothly towards 0.9 at the inner-corner cells, then sharply to 0 at the boundary.
- **D₄ self-consistency check** (random sample of 500 tours × 7 non-identity D₄ symmetries = 3,500 mapped pairs): the multiset of `|w(p)|` is preserved in 100% of pairs — confirming the winding maps respect D₄ up to grid orientation / direction reversal.
- **Correlation with the 7 H₁ invariants — caveat.** The exact 7-bit basis of the parent project is not stored per-tour in any of the existing JSONs. The natural D₄-quotient invariant per tour that *is* readily computable is the GF(2)-parity of the tour's edge-incidence inside each of the 10 D₄ edge orbits. Empirically these 10 parity bits have GF(2)-rank **4** across the 9,862 tours, so they capture *less* than the project's 7-bit basis — there is genuine homology information not captured by orbit-parities. **Pearson correlation between winding numbers and these 4 effective bits is weak:** max |r| = 0.145 (cell `(3,4)` vs. bit 8). A more interesting cross-correlation will likely require recovering the actual 7-bit cycle basis from `../cavalo_loop_destruicao_6x6.py` directly — see Open Question 2 below.
- Plot: `data/plots/winding_mean_heatmap.png` (mean `|w(p)|` over all tours).

### Task D — quantum walk (`quantum_walk.py`)

- Continuous-time quantum walk on the adjacency `A` of the 6×6 knight graph, starting from `|A6⟩`. Step `τ = 1.0`. 100 steps. Unitarity holds to machine precision (max norm deviation `1.1e-15`).
- Spectral radius `‖A‖₂ ≈ 4.987`.
- The walk does **not** mix to the uniform distribution. TV-distance to uniform stays > 0.36 throughout 100 steps; minimum 0.36 at `t = 96`, far from zero.
- At the "tour length" `t = N = 36` the probability concentrates **back on the starting cell A6** (max prob = 0.32), with TV-distance to uniform 0.568. This is a clean signature of **quantum recurrence** on a structured graph — the knight graph is small enough and regular enough that interference brings amplitude back to the start.
- Plots: `data/plots/quantum_walk_snapshots.png` (heatmaps at t=0, 10, 36, 100), `data/plots/quantum_walk_uniform_distance.png` (TV vs. t).

---

## 5. Classificação dos invariantes H₁ (`h1_invariants.py`, 2026-05-17)

Esta seção fecha o gap deixado pelo `winding_number.py`, que usou uma proxy
(paridade-por-órbita D₄ de arestas) cobrindo só rank-4 do subespaço de
invariantes. Aqui o cálculo é feito diretamente sobre a matriz de
incidência aresta–tour `T ∈ GF(2)^{9862 × |E|}`.

### 5.1 — Inventário linear

| Quantidade                                                                       | Valor | Comentário                                                                                                |
|----------------------------------------------------------------------------------|------:|-----------------------------------------------------------------------------------------------------------|
| `V`, `E`, β₁ do grafo do cavalo 6×6                                              | 36, 80, 45 | β₁ = E − V + 1 (grafo conexo).                                                                       |
| `rank(T)` sobre GF(2) (arestas **não-direcionadas**)                              | **42** | Os 9862 tours geram um subespaço próprio de H₁ (3 ciclos de H₁ não são expressáveis como ⊕ de tours).  |
| `rank(T_dir)` sobre GF(2) (arestas **direcionadas**, 2·E = 160 colunas)           | **78** | Tours direcionados não preenchem GF(2)^{160} mesmo aproximadamente.                                  |
| `n_distinct_signatures` em `T @ B^T mod 2` com base pivô de T                     | 9862  | Cada tour tem signature única — esperado dado que rank=42 ≫ log₂(9862) ≈ 13.3.                       |
| `rank({integer winding W[:, p] | p ∈ 36 cells})` sobre ℝ                          | **21** | Só 21 das 36 séries de winding inteiro são ℝ-independentes.                                          |
| `rank({W[:, p] mod 2 | p ∈ 36 cells})` sobre GF(2)                                | **18** | Só 18 invariantes geométricos GF(2)-LI extraíveis de winding.                                        |
| `# D₄-orbits` no quadrante interior 5×5 das células                               | 6     | Burnside no 5×5 com D₄ ⇒ 6 órbitas (4 + 8 + 4 + 4 + 4 + 1).                                          |
| `rank({Σ_{p∈orb} W[:, p] mod 2 | orb})` sobre GF(2)                               | **4** | XOR sobre órbita D₄ dá só 4 funcionais LI (uma órbita é sempre 0, outra é dependente).               |

**Surpresa central** (Task 2b, basis-independent):

> ∀ célula `p`, a parida winding(tour, p) mod 2 é uma função **GF(2)-linear das
> arestas não-direcionadas** do tour. As 36 funcionais (uma por célula) vivem
> num subespaço 18-dimensional de GF(2)^9862.

Isto é um teorema topológico: o "linking number mod 2" entre uma curva
fechada e um ponto é a paridade do número de cruzamentos com qualquer raio
emergindo do ponto, e cada cruzamento depende apenas da aresta atravessada
(não da direção). Portanto, mod 2 é insensível à direção do passeio — daí a
linearidade em `T` (não em `T_dir`). Sobre ℤ, isso falha: 24/36 células
exigem `T_dir` (R² < 1 em `T` solo).

### 5.2 — Por que a tabela literal de "42 invariantes" do `invariant_classification.json` é toda algébrica

A pergunta literal do enunciado ("invariantes da base pivô × winding por
célula, threshold de Pearson 0.3") foi computada e respondida: **todas as
42 invariantes ficaram classificadas como `ALGÉBRICAS`** (max |r| em GF(2)
foi 0.116, no invariante #39).

O motivo é estrutural, não vazio: bits de uma base arbitrária (linhas-pivô
de `T`) misturam ~36 arestas cada um. Cada bit é uma "soma de produtos"
sem semântica geométrica local. Como a winding-por-célula só "vê" arestas
locais a `p`, a sobreposição entre dois objetos não-alinhados é diluída.
Pearson de 0.1 com 9862 amostras é tudo, menos ruído (`p`-valor ≪ 10⁻³⁰),
mas o sinal por bit é genuinamente fraco.

O teste **basis-independent** (Task 2b) responde à pergunta correta: é
*toda* informação winding-mod-2 capturável por GF(2)-funcionais lineares
sobre `T`. **Resposta: sim**, com 18 dimensões de informação independente.

### 5.3 — Sobre o "7" do projeto principal

O texto do enunciado supunha `rank(T) = 7` e 7 invariantes H₁. **Empiricamente
isso não se confirmou em nenhum dos cálculos acima** — os números relevantes
são 42 (rank de T), 18 (rank de winding mod 2), e 4 (rank D₄-orbit-quotient).

Hipótese mais provável para o "7" original:

- O conjunto de **arestas obrigatórias** identificadas em `../analysis_6x6.json`
  está em uma única órbita D₄ de tamanho 8; combinado com a relação "todo
  tour tem 36 arestas", isso pode reduzir um "rank 8" a "rank 7" de algum
  subconjunto invariante.
- Os "7 invariantes" do projeto principal podem ser **não-lineares**: por
  exemplo, paridades de destruição de loops em **passos** específicos do
  caminho, que dependem do estado dinâmico, não da incidência. Esse tipo
  de invariante NÃO é captável por nenhuma análise linear sobre `T` ou `T_dir`,
  por construção.

→ Para investigar: importar o destruction tracking de
`../cavalo_loop_destruicao_6x6.py` e calcular, por tour, o vetor binário
"destruiu loop O_k no passo t_k?" para k = 0..6.

### 5.4 — Tabela final dos invariantes (literal, basis-de-pivô de T)

42 bits, todos `ALGEBRAIC` pelo threshold |r| > 0.3. JSON completo em
`data/invariant_classification.json`. Top três por max |r| (mod 2):

| bit | n₁/N (%) | max |r| (W) | max |r| (W mod 2) | célula testemunha mod 2 |
|----:|---------:|------------:|------------------:|--------------------------|
| 39  |   49.7%  | 0.018       | **0.116**         | (provável (1,3) ou simétrica) |
| 12  |   51.3%  | 0.043       | 0.052             | —                        |
| 25  |   51.0%  | 0.019       | 0.050             | —                        |

Plot por bit: `data/plots/invariant_classification.png` (heatmaps 6×6).
Plot basis-independent: `data/plots/winding_linearity_basis_free.png`.

### 5.5 — Implicação para o salto 10×10

No 10×10 cresce `V = 100`, `E = ?` (≈ 384 = 8 × 100 / 2 corrigido nas
bordas), β₁ ≈ 285. O subespaço gerado pelos tours em GF(2)^{384} terá
rank previsivelmente próximo de 200–280 (extrapolando a relação
rank(T)/β₁ ≈ 42/45 ≈ 0.93 do 6×6).

Resultado prático: tentar enumerar **invariantes GF(2)-lineares** sobre
arestas (estilo 5.1 acima) gera ordens de magnitude mais quantidades do
que os "9 órbitas / médias 18.6 invariantes" reportados na [[project-fase-c-done]].
**Conclusão estratégica**: a Fase C/D do projeto principal deve continuar
usando os invariantes **não-lineares** (correlações de destruição,
[[project-fase-b-final]]) — eles capturam estrutura genuinamente
dinâmica que a análise GF(2)-linear não captura. A linearização sobre
GF(2) é útil principalmente como **upper-bound** ("se um invariante é
GF(2)-linear em arestas, ele cabe num subespaço de dim ≤ 42 no 6×6 e
≤ ~280 no 10×10"), não como gerador novo de invariantes.

---

## 6. Open questions / next steps

1. **DFT vs. canonical classes.** Do the magnitudes `{|c_k|}` separate the
   1,232 D₄-canonical classes? If two D₄-inequivalent tours have identical
   `(|c_0|, …, |c_{N-1}|)` they would be a *spectral collision* — easy
   experiment using the saved npy.
2. **Recover the true 7-bit H₁ basis.** The 10 orbit-parity bits used in
   Task C only span a rank-4 subspace of the H₁ invariant space, so the
   weak winding↔bit correlation (max 0.145) is on a degraded proxy.
   Construct the genuine 7-element GF(2) basis from `../cavalo_loop_destruicao_6x6.py`'s
   per-tour destruction signatures and rerun the correlation matrix.
3. **Winding ↔ loop destruction.** The conservation law `Σ destroyed = 45`
   holds across all tours and equals `H₁ = 45`. The mean `|w(p)|` over
   interior cells (≈ `Σ_p |w(p)| ≈ 27` per tour, from `data/winding_correlation.json`)
   should encode part of the same homological information. Is there a
   linear functional `f(w) ≈ #(destroyed loops at step t)`?
4. **CTQW recurrence at `t = N`.** The walk concentrates back on the
   starting vertex at `t = 36`. The 36 step-times of a closed tour and
   the eigenvalue spacings of the knight adjacency may be in resonance.
   Verify: compute the IPR `(Σ_v |ψ_v|⁴)⁻¹` at every `t` and look for
   minima at `t ∈ {N, 2N, 3N, …}`.
5. **Discrete vs. continuous quantum walks.** Replace CTQW with a
   discrete-time coined walk (Grover coin on the 8-regular knight graph
   restricted to interior cells, with reflective boundary). Compare
   mixing behaviour — the discrete walk is closer in spirit to the
   classical knight tour itself.

---

## 7. Links back to the main project

- `../analysis_6x6.json` — canonical counts and edge orbits
- `../b6_grouping_6x6.json` — orbit groupings
- `../b6_interior_correlations_6x6.json` — interior loop correlations
- `../cavalo_loop_destruicao_6x6.py` — exhaustive backtracking & H₁ tracking
- `../destruction_map.json` — per-square destruction map

---

## Run log

### 2026-05-17 — segundo turno (h1_invariants)

- `h1_invariants.py`: criou ∂₁ (36×80) e T (9862×80) em GF(2); confirmou
  `T · ∂₁ᵀ ≡ 0` (todo tour é ciclo); **rank(T) = 42** sobre GF(2).
- Base pivô e coordenadas salvas em `data/h1_basis_gf2.npy` (42×80) e
  `data/h1_coordinates_gf2.npy` (9862×42). 9862 signatures distintas.
- D₄-invariância das coordenadas pivô: 0.14% (98/69034) — esperado, base
  não é D₄-invariante por construção.
- Pearson per-bit × winding: **todas as 42 bits classificadas como
  ALGÉBRICAS** com threshold 0.3 (max |r| = 0.116).
- Análise basis-independent: **mod-2 winding é GF(2)-linear em T para todas
  as 36 células** (resultado topológico); 18 dimensões geométricas LI.
- D₄-orbit quotient (5×5 interior): 6 órbitas, 4 funcionais LI sob XOR.
- Conclusão registrada no README §5: o "7 H₁ invariants under D₄" da
  memória provavelmente referencia invariantes **não-lineares** (parities
  de destruição), não a estrutura linear que `h1_invariants.py` mede.

### 2026-05-17 — first end-to-end run

- `_tours.py` enumeration: 19,724 directed Hamiltonian cycles anchored at
  v=0 → 9,862 canonical (undirected) tours, exhaustive backtracking with
  Warnsdorff ordering, ~70 minutes wall time. Cached to
  `data/tours_closed_6x6.npy` (355 KB, int8).
- `complex_moves.py`: all Task-A self-tests pass.
- `fourier_tours.py`: DFT saved (5.4 MB); centroid sanity check exact;
  mean spectrum dominated by `k = 3, 4` (and mirror `k = 32, 33`), not by
  the fundamental `k = 1`.
- `winding_number.py`: winding maps saved (348 KB); D₄ multiset
  invariance holds for 3,500 / 3,500 sampled pairs (rate 1.0); rank-4
  orbit-parity proxy has weak correlation (`max |r| = 0.145`) with
  windings. 11 outer cell-centres consistently have `w = 0`.
- `quantum_walk.py`: 100 CTQW steps, unitarity to 1e-15; at `t = N = 36`
  the walk recurs to the starting vertex (max prob 0.32 at A6, TV to
  uniform 0.568).
