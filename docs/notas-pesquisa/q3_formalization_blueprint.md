# Q(n)=3 General — Formalization Blueprint

> Target: a separate Lean agent should be able to close `Q_eq_three_general`
> (`KnightTour/TheoremQ3.lean:60`) for all `n ≥ 6` by mechanically following
> this document, **without re-deriving any mathematics**. The maths is settled
> (paper Thm 2.5, Lemmas 2.6/2.7/2.10/2.11, Remark 2.12). The only ingredient
> that was geometric and open — bulk connectivity — is already Lean-proved
> sorry-free as `bulk_connected_general` (`KnightTour/BulkConnectivity.lean:235`).
> What remains is **pure linear algebra over F₂ = ZMod 2**.

All Mathlib names below were checked against the project's pinned source tree
`knight_tour_lean/.lake/packages/mathlib` at `inputRev = v4.16.0`
(`lean-toolchain = leanprover/lean4:v4.16.0`). Line references are to that tree.

---

## Exact paper statements (transcribed)

Notation from the paper (Def 2.2, Def 2.4, Lemmas 2.6/2.7/2.10/2.11, Thm 2.5,
Remark 2.12). The Lean source labels differ slightly (it calls the independence
lemma 2.9); this blueprint uses the **paper** numbering and cross-references the
`.tex` label keys.

- **Graph** (`def:cantos`, `Def Grafo do cavalo`): `G_n = (V_n, E_n)`,
  `V_n = {(i,j) : 0≤i,j<n}`, `E_n = { {u,v} : |u₁−v₁|·|u₂−v₂| = 2 }`.
- **Boundary / cycle space** (`Def 2.2`): `∂₁ ∈ F₂^{|V|×|E|}` the
  vertex–edge incidence matrix; `Z₁ = ker(∂₁) ⊆ F₂^{|E|}`;
  `β₁ = |E| − |V| + 1`.
- **Corners & mandatory edges** (`def:cantos`): `Corners(n) =
  {(0,0),(0,n−1),(n−1,0),(n−1,n−1)}`. Each corner `c` has degree 2
  (Lemma 2.6 = `lem:grau`), its two edges `mand(c)={e₁(c),e₂(c)}`;
  `Mand(n)=⋃_c mand(c)`, `|Mand|=8`.
- **Def 2.4 (`def:Q`)** — for `i,j ∈ Mand(n)`, `i<j`, let `v_ij ∈ F₂^{E_n}` be
  the indicator of `{i,j}`. Let `R(n)=row(∂₁)`, `π: F₂^{E_n} → F₂^{E_n}/R(n)`.
  Then
  `Q(n) := dim_{F₂}( π( Span{ v_ij : i,j∈Mand, i<j } ) )`.
- **Lemma 2.7 (`lem:bordas`)** — for every corner `c`,
  `e₁(c) + e₂(c) = ∂(δ_c) ∈ R(n)`, where `δ_c ∈ F₂^{V_n}` is the indicator of
  `c`. ⇒ **Cor 2.8 (`cor:rc`)**: `[e₁(c)] = [e₂(c)] =: r_c` in `F₂^{E_n}/R(n)`.
- **Lemma 2.10 (independence, `lem:indep`)** — `r_{c₁},…,r_{c₄}` are linearly
  independent in `F₂^{E_n}/R(n)`. Proof (a), `n≥6`: if `Σ_{c∈S} r_c = 0` for
  nonempty `S⊆Corners`, then `Σ_{c∈S} e₁(c) ∈ R(n)`, so `∃ α∈F₂^{V_n}` with
  `∂(α) = Σ_{c∈S} e₁(c)`; by **bulk connectivity** `α` is constant `=a` off the
  corners; `[∂α]_{e₁(c)} = α[c]+a = 1 ⇒ α[c]=1+a ⇒ [∂α]_{e₂(c)} = (1+a)+a = 1`,
  contradicting `e₂(c) ∉ supp(∂α) = {e₁(c'):c'∈S}` (disjoint corner
  neighbourhoods, `n≥6`).
- **Lemma 2.11 (XOR image, `lem:xor`)** —
  `π(Span{v_ij}) = Span{ r_{c_i}+r_{c_j} : 1≤i<j≤4 }`.
  Proof: same-corner pair ⇒ `v_ij∈R(n)` ⇒ `π=0`; cross-corner pair
  `i∈mand(c_a), j∈mand(c_b), a≠b` ⇒ `π(v_ij)=r_{c_a}+r_{c_b}`.
- **Thm 2.5 (`thm:Qn3`)** — `W := Span(r_{c₁},…,r_{c₄}) ≅ F₂⁴`. With
  `σ: W→F₂`, `σ(Σ a_k r_{c_k})=Σ a_k`, every `r_{c_i}+r_{c_j} ∈ ker σ`,
  `dim ker σ = 3`, and `{r_{c₁}+r_{c₂}, r_{c₁}+r_{c₃}, r_{c₁}+r_{c₄}}` are LI in
  `ker σ`, so `Q(n)=3`.
- **Remark 2.12 (`rem:sigma`)** — the “−1” is functional: `Σ_c r_c ≠ 0` in the
  quotient (`ker σ` has codim 1 in `W`), not geometric.

---

## Route decision (Part A)

### The two routes

**ROUTE 1 — Abstract-only.** Define everything as Mathlib objects:
`δ := coboundary LinearMap (Fin n×Fin n → F₂) → (E_n → F₂)`;
`R(n) := LinearMap.range δ`; quotient `Qt := (F₂^E ⧸ R(n))` via
`Submodule.Quotient`; `Q_abstract(n) := finrank F₂ ((Submodule.span F₂ {v_ij}).map R(n).mkQ)`.
Prove `Q_abstract(n) = 3` directly. `computeQ` is untouched by the general proof
(it survives only for concrete-`n` `native_decide`).

**ROUTE 2 — Bridge `computeQ`.** Prove `computeQ n hn = Q_abstract n` for all
`n` (i.e. verify the imperative `gf2Rank` Gaussian elimination computes the
quotient dimension), then prove `Q_abstract n = 3`.

### Cost of ROUTE 2 (rejected)

`gf2Rank` (`GF2Space.lean:61`) is an `Id.run do` loop with `Array.set!`,
`break`, nested `for ... in [rank:m]`, and in-place row reduction on
`Array (Array GF2)`. To prove `computeQ n = Q_abstract n` generally one must:

- give a loop invariant tying the mutable `mat`/`rank` to the row space and the
  pivot structure after column `col`;
- discharge every `set!`/`!`-index obligation (`getElem!` on `Array`, default
  values) — these are partial and reason about bounds the loop maintains;
- prove the final `rank` equals `finrank` of the span of the input rows
  (a correctness theorem for GF(2) Gaussian elimination), then *separately* do
  all of Route 1 anyway to get the value `3`.

This is a self-contained verified-imperative-linear-algebra project on top of
the entire Route-1 work. **It strictly dominates Route 1 in cost and buys
nothing** (`computeQ` already gives `native_decide` theorems for `n∈{4..10}`,
which is all we need it for). Reject.

### Cost of ROUTE 1 (chosen)

New objects to build (all standard Mathlib constructions): one `LinearMap`
(`δ`), one `Submodule` (`R(n) = range δ`), the quotient (already a `Submodule`
mechanism), one finite spanning set (`{v_ij}`), four explicit quotient vectors
`r_c`. Every finrank step is a named Mathlib lemma (table below). The hardest
single sub-step is **not** linear algebra — it is the new “locally-constant on a
connected graph ⇒ globally constant” lemma feeding Lemma 2.10 (Part C, Hard
Core #1), which is a clean `SimpleGraph.Walk` induction on top of the already
proved `bulk_connected_general`.

### A genuinely better variant of Route 1 (adopted): drop the σ functional

Remark 2.12’s `σ` is the *conceptual* reason for the “−1”, but **formalizing
`σ` is unnecessary and harder than the alternative.** Constructing `σ` as a
well-defined linear functional on `W` requires that `{r_{c_k}}` be a *basis* of
`W` (to read off coordinates `a_k`) — i.e. it presupposes Lemma 2.10 anyway —
and then needs `LinearMap.finrank_range_add_finrank_ker` to get `dim ker σ = 3`.

Instead, use the elementary identity over F₂:

> `Span{ r_{c_i}+r_{c_j} : i<j } = Span{ w₂, w₃, w₄ }` where `w_k := r_{c₁}+r_{c_k}`,

because `r_{c_i}+r_{c_j} = w_i + w_j` (with `w₁ = 0`). And `{w₂,w₃,w₄}` is
linearly independent **directly from** Lemma 2.10 (any F₂-dependence among the
`w_k` rearranges to a nonempty-subset dependence among the `r_{c_k}`, excluded
by Lemma 2.10). Then `finrank = 3` by `finrank_span_eq_card`. **No quotient
functional, no rank–nullity.** This is the adopted endgame (Steps 7–9 below).

**ROUTE CHOSEN: 1 (abstract-only), σ-free endgame.**

---

## Step decomposition (Part B)

Types: `V := Fin n × Fin n` (= `Square n`, `Basic.lean:18`); `F := ZMod 2`
(= `GF2`, `Basic.lean:15`). For the edge index type use a `Fintype` of
undirected edges — see Step 0 for the representation choice, which is the one
real design decision feeding all later steps.

### Step 0 — Edge index type `Edge n` (foundational design choice)

- **Statement.** Fix a `Fintype`/`DecidableEq` type `Edge n` indexing the
  undirected knight edges, with an incidence predicate
  `inc : V → Edge n → Prop` (decidable) s.t. each edge has exactly two
  incident vertices. Define the F₂ vector spaces `CV := V → F` and
  `CE := Edge n → F` (both `= Fintype.card`-dimensional via
  `Module.finrank_fintype_fun_eq_card`).
- **Paper source.** Def 2.2 (`∂₁ ∈ F₂^{|V|×|E|}`).
- **Mathlib tools.** `Module.finrank_fintype_fun_eq_card`
  (`Dimension/Constructions.lean:305`); `Pi.module`; `ZMod.instField`
  (`Data/ZMod/Basic.lean:1062`) with `Nat.fact_prime_two`
  (`Data/Nat/Prime/Defs.lean:439`) to get `Field F`.
- **Recommended representation.** `Edge n := { e : Sym2 V // KAdj' e }` *or* a
  `Finset (V×V)` of canonical edges. **Recommendation: `Sym2 V`-based**, because
  `Sym2` gives symmetric incidence (`Sym2.mem`) for free and matches Mathlib’s
  `SimpleGraph.edgeSet`. The existing `KAdj` (`BulkConnectivity.lean:37`) already
  provides the adjacency predicate; reuse it to define the knight `SimpleGraph`
  on all of `V` (the bulk graph is its induced subgraph). Avoid the
  `List`-based `undirectedEdges` (`GF2Space.lean:19`) for the abstract proof —
  it is for `computeQ`/`native_decide` only.
- **Dependencies.** none.
- **Difficulty.** ROUTINE (boilerplate Fintype/incidence; the only care is
  proving “exactly two incident vertices” = restatement of Lemma 2.6 for the two
  corner edges and not needed in general for interior edges).

### Step 1 — Coboundary `δ : CV →ₗ[F] CE` and `R(n) := range δ`

- **Statement.** `δ : (V → F) →ₗ[F] (Edge n → F)`,
  `δ α e = ∑_{v inc e} α v` (over F₂: `α u + α w` for `e={u,w}`). Define
  `R : Submodule F CE := LinearMap.range δ`.
- **Paper source.** Def 2.2 / `lem:bordas` (`∂(δ_v)` = incidence row of `v`;
  here `δ` is `∂₁ᵀ`, matching `boundaryVector`, `GF2Space.lean:28`). Note the
  orientation: the paper’s `row(∂₁) ⊆ F₂^E` is the **range of `∂₁ᵀ = δ`**, i.e.
  the column space of the transpose = row space of `∂₁`. Match this.
- **Mathlib tools.** `LinearMap.mk` / build via `Finset.sum`; `LinearMap.range`
  (`Submodule`); `LinearMap.range_eq_map`. Linearity of `δ` is `map_add`/
  `map_smul` from `Finset.sum` linearity.
- **Dependencies.** Step 0.
- **Difficulty.** ROUTINE. (`δ` is a finite signed-incidence sum; building it
  as a `LinearMap` is standard. Key sanity lemma to prove here:
  `δ α e = α u + α w` for the two endpoints — a `simp` over `Sym2`.)

### Step 2 — Quotient `Qt := CE ⧸ R(n)` and projection `π := R.mkQ`

- **Statement.** `Qt := CE ⧸ R`; `π := (R).mkQ : CE →ₗ[F] Qt`; it is surjective.
- **Paper source.** Def 2.4 (`π : F₂^{E_n} → F₂^{E_n}/R(n)`).
- **Mathlib tools.** `Submodule.Quotient` (instance `AddCommGroup`/`Module`);
  `Submodule.mkQ` (`Quotient/Basic.lean`); `Submodule.range_mkQ`
  (`Quotient/Basic.lean:157`, `= ⊤`, gives surjectivity);
  `Submodule.Quotient.mk_eq_zero` / `Submodule.mkQ_apply`;
  membership `x - y ∈ R ↔ π x = π y` via `Submodule.Quotient.eq`.
- **Dependencies.** Step 1.
- **Difficulty.** TRIVIAL (entirely library).

### Step 3 — Mandatory edges, `v_ij`, and `Span{v_ij}` as a Submodule

- **Statement.** Define `mandEdges n : Finset (Edge n)` = the 8 corner-incident
  edges (two per corner, via Lemma 2.6). For `i,j ∈ mandEdges`, `v i j : CE` :=
  `Pi.single i 1 + Pi.single j 1`. Define
  `S := Submodule.span F (↑{ v i j | i,j ∈ mandEdges, i<j } : Set CE)`.
- **Paper source.** Def 2.4.
- **Mathlib tools.** `Pi.single`; `Submodule.span`; `Finset.image`/`Finset.offDiag`
  to enumerate `i<j`. To connect to `computeQ`’s `xorVector`/`allXorVectors`
  *only at concrete n* (not needed generally).
- **Dependencies.** Step 0, plus Lemma 2.6 (`knightDegree n c = 2`, already in
  `Lemma1.lean`) to know each corner contributes exactly two edges.
- **Difficulty.** ROUTINE.

### Step 4 — Corner representatives `r_c` and Cor 2.8

- **Statement.** For each corner `c` with edges `e₁(c),e₂(c)`:
  (i) `Pi.single (e₁ c) 1 + Pi.single (e₂ c) 1 = δ (Pi.single c 1)`  (Lemma 2.7);
  hence (ii) `π (Pi.single (e₁ c) 1) = π (Pi.single (e₂ c) 1) =: r c` in `Qt`.
  Define `r : Corners → Qt`, `r c := π (Pi.single (e₁ c) 1)`.
- **Paper source.** Lemma 2.7 (`lem:bordas`) + Cor 2.8 (`cor:rc`).
- **Mathlib tools.** (i) is a finite computation: `δ (Pi.single c 1) e =
  [c inc e]`, and by Lemma 2.6 only `e₁ c, e₂ c` are incident to `c` — a
  `Finset.sum`/`Sym2` `decide`-free unfolding (or `Finset.sum_eq` with the
  two-element incident set). (ii): `Submodule.Quotient.eq` +
  `LinearMap.mem_range_self`.
- **Dependencies.** Steps 1, 2; Lemma 2.6 (`Lemma1.lean`).
- **Difficulty.** ROUTINE. The content is “the incidence row of a corner is
  supported exactly on its two edges”, i.e. a restatement of Lemma 2.6 at the
  vector level.

### Step 5 — Lemma 2.10: `r` is linearly independent in `Qt`  **(HARD — see Part C #1)**

- **Statement.** `LinearIndependent F (r : Corners → Qt)`. Equivalent F₂
  formulation actually proved: `∀ S : Finset Corners, S.Nonempty →
  ∑_{c∈S} r c ≠ 0`.
- **Paper source.** Lemma 2.10(a) (`lem:indep`), the bulk-connectivity argument.
- **Mathlib tools.** `Fintype.linearIndependent_iff` (reduces F₂ LI to
  “no nonempty subset sums to 0”); the new lemma `bulk_const` (Part C #1) built
  on `bulk_connected_general` (`BulkConnectivity.lean:235`) via
  `SimpleGraph.Reachable.elim` + `Walk` induction; `Submodule.Quotient.eq` /
  `LinearMap.mem_range` to turn `∑ r c = 0` into `∃ α, δ α = ∑ e₁(c)`.
- **Dependencies.** Steps 1, 2, 4; `bulk_connected_general`; disjointness facts
  (Step 5b).
- **Difficulty.** HARD. The linear algebra is trivial; the work is the
  graph-constancy lemma and the F₂ support bookkeeping. Full detail in Part C.

### Step 5b — Disjointness / distinctness of the 8 mandatory edges (support facts)

- **Statement.** For `n ≥ 6`: (D1) the two knight-neighbours of each corner are
  **non-corner** (bulk) vertices; (D2) the 8 vectors `e₁(c),e₂(c)` over the four
  corners are pairwise distinct edges; equivalently corner neighbourhoods are
  disjoint, so `e₂(c) ∉ {e₁(c') : c'∈Corners}` and each `e₁(c)` appears with
  multiplicity exactly 1 in `∑_{c∈S} e₁(c)`.
- **Paper source.** Lemma 2.10(a), the “separação ≥ n−3 ≥ 3 para n≥6” clause.
- **Mathlib tools.** explicit coordinates `(1,2),(2,1)` near `(0,0)`,
  `(1,n−3),(2,n−2)` near `(0,n−1)`, etc.; `omega` after unfolding `KAdj`/corner
  coordinates. No deep lemma — purely arithmetic over `Fin n` with `n≥6`.
- **Dependencies.** Step 0.
- **Difficulty.** ROUTINE but **tedious** (4 corners × 2 edges, pairwise: ~28
  `omega` goals). Genuinely general; not n=6-specific. Budget real time here.

### Step 6 — Lemma 2.11: `S.map π = Span{ r_{c_i} + r_{c_j} : i<j }`

- **Statement.** `(S).map π = Submodule.span F (↑{ r a + r b | a b : Corners,
  a≠b } : Set Qt)` (over F₂ the cross-corner pairs give exactly these; same-corner
  pairs map to 0).
- **Paper source.** Lemma 2.11 (`lem:xor`).
- **Mathlib tools.** `Submodule.map_span`
  (`Algebra/Lie/IdealOperations.lean:204` usage; lemma `Submodule.map_span` in
  `Mathlib/LinearAlgebra/Span.lean`) to push `π` through the span:
  `(span F T).map π = span F (π '' T)`. Then a set-image computation:
  `π (v i j) = 0` if `i,j` same corner (Lemma 2.7 ⇒ `v_ij ∈ R`), and
  `π (v i j) = r(c_a)+r(c_b)` if cross-corner (additivity of `π` + Step 4).
- **Dependencies.** Steps 2, 3, 4; Lemma 2.7 (Step 4(i)).
- **Difficulty.** ROUTINE. (Image of a span = span of image; the only content is
  the per-pair case split, which is finite over the 8 mandatory edges.)

### Step 7 — Reduce pairwise-sum span to `{w₂,w₃,w₄}` (σ-free)

- **Statement.** With `w k := r c₁ + r c_k`, prove
  `Span{ r a + r b : a≠b } = Span{ w₂, w₃, w₄ }`, using `r a + r b = w a + w b`
  and `w₁ = 0`.
- **Paper source.** Endgame of Thm 2.5 (replaces the `σ`/`ker σ` step of
  Remark 2.12 with the equivalent explicit basis).
- **Mathlib tools.** `Submodule.span_le` / `Submodule.span_mono` both directions,
  or `le_antisymm` with `Submodule.mem_span_pair` / `Submodule.subset_span`;
  arithmetic `r a + r b = (r c₁ + r a) + (r c₁ + r b)` over `AddCommGroup`
  (char 2: `x+x=0`).
- **Dependencies.** Step 4 (the `r c`), Step 6.
- **Difficulty.** ROUTINE.

### Step 8 — `{w₂,w₃,w₄}` is linearly independent  ⇒ finrank = 3

- **Statement.** `LinearIndependent F ![w₂,w₃,w₄]`, hence
  `finrank F (Span{w₂,w₃,w₄}) = 3`.
- **Paper source.** Thm 2.5 (“`{r_{c₁}+r_{c₂},…}` LI in `ker σ`, dim 3”);
  here phrased without `σ`.
- **Mathlib tools.** `finrank_span_eq_card`
  (`Dimension/Constructions.lean:447`, hyp: `[Nontrivial F]` (✓ for `ZMod 2`),
  `Fintype ι`, `LinearIndependent F b` ⇒ `finrank (span (range b)) = card ι`).
  The LI of `{w₂,w₃,w₄}` follows from Step 5 (LI of `r`): an F₂-dependence
  `∑_{k∈T} w_k = 0` (`T ⊆ {2,3,4}`) expands to
  `(|T| mod 2)·r c₁ + ∑_{k∈T} r c_k = 0`, a nonempty-subset dependence among
  `{r c₁,…,r c₄}`, excluded by Step 5. Encode via
  `Fintype.linearIndependent_iff` on both ends, or `LinearIndependent.map'`
  pushing the standard basis of F₂³ through the injective map `e_k ↦ w_k`.
- **Dependencies.** Steps 5, 7.
- **Difficulty.** ROUTINE (the algebra is the “rearrange to subset of `r`”
  bookkeeping; mechanically a `Fin 3`/`Finset` case analysis or a clean
  `LinearIndependent.map'`).

### Step 9 — Assemble `Q_abstract n = 3`

- **Statement.** `finrank F ((S).map π) = 3`; define `Q_abstract n` as this and
  conclude.
- **Paper source.** Thm 2.5 conclusion.
- **Mathlib tools.** chain: Step 6 (`(S).map π = Span{r a+r b}`) ▸ Step 7
  (`= Span{w₂,w₃,w₄}`) ▸ Step 8 (`finrank = 3`). Pure `rw`/`calc`.
- **Dependencies.** Steps 6, 7, 8.
- **Difficulty.** TRIVIAL.

### Step 10 (optional, for parity with `computeQ`) — concrete `n` bridge

- **Statement.** `Q_eq_three_general` as currently typed proves
  `computeQ n hn = 3`. With Route 1 we instead prove `Q_abstract n = 3`. To keep
  the *exact existing signature* (`computeQ n hn = 3`), one still needs
  `computeQ n hn = Q_abstract n` — the rejected Route-2 bridge.
- **Recommendation.** **Restate the target.** Change the public theorem to
  `theorem Q_eq_three_general (n) (hn : 6 ≤ n) : Q_abstract n = 3`, and keep
  `computeQ`-based `native_decide` theorems for `n ∈ {4..10}` as the concrete
  cross-check (they already exist for `n=6`, `TheoremQ3.lean:23`). Document that
  `Q_abstract` is the paper’s `Q(n)` verbatim (Def 2.4) and `computeQ` is a
  decidable shadow agreeing on all tested `n`. This is honest and avoids
  verifying imperative Gaussian elimination. *(Note `n≥6`: for `n∈{4,5}` keep the
  `native_decide` route; the structural proof needs `n≥6` for bulk connectivity.)*
- **Difficulty.** TRIVIAL (it is a definitional/relabelling decision), but it is
  a **decision the Lean agent must be told to make** — flagged in Risk register.

---

## Hard core (Part C)

### Hard Core #1 — Step 5: independence via bulk connectivity

**What makes it hard.** Two sub-parts: (i) a *new* Mathlib-absent lemma turning
`(Bulk n).Connected` into “`α` constant off the corners”; (ii) F₂ support
bookkeeping connecting `δ α = ∑ e₁(c)` to a contradiction. Neither is deep, but
(i) requires a `Walk` induction and (ii) requires the Step 5b distinctness facts
to be in hand.

**Mathematics is settled** — paper Lemma 2.10(a), and the geometric core
(`bulk_connected_general`) is already Lean-proved sorry-free. No gap (see the
explicit n-generality trace below).

**Proof sketch for a Lean agent.**

*New lemma `bulk_const`* (the load-bearing absent lemma):
```
lemma edgewise_const_of_connected
    {W : Type*} (H : SimpleGraph W) (hH : H.Connected)
    (α : W → F) (hedge : ∀ a b, H.Adj a b → α a = α b) :
    ∀ a b, α a = α b
```
Proof: `hH.preconnected a b` gives `Reachable`; `Reachable.elim` extracts a
`Walk a b`; induct on the walk (`Walk.rec` / `induction p`): `nil` ⇒ `rfl`;
`cons (h : Adj a x) (q : Walk x b)` ⇒ `hedge a x h ▸ (ih)`. (Mathlib has
`Reachable.elim` `Path.lean:745`, `Walk` and its recursor in `Path.lean`.) This
is ~10 lines, the only genuinely new lemma in the whole blueprint.

*Apply it to the bulk.* `Bulk n` = `SimpleGraph` on `{v // ¬ isCorner n v}`
(`BulkConnectivity.lean:43`), already proved `Connected` for `n≥6`. Given
`α : V → F` from `δ α = ∑_{c∈S} e₁(c)`:
- For every **bulk edge** `e = {u,w}` (both `u,w` non-corner): `e ≠ e₁(c)` for
  all `c` (an `e₁(c)` is incident to a corner; `u,w` are not), so
  `(∑_{c∈S} e₁(c)) e = 0`, i.e. `δ α e = α u + α w = 0`, i.e. `α u = α w`
  (char 2). This is exactly `hedge` for `Bulk n` restricted to `α ∘ (↑)`.
- `edgewise_const_of_connected (Bulk n) (bulk_connected_general n hn) (α∘↑) …`
  ⇒ `α` is constant `= a` on all non-corner vertices.

*Derive the contradiction.* Fix `c ∈ S` (nonempty). Its neighbours `p,q`
(`= e₁(c) = {c,p}`, `e₂(c) = {c,q}`) are bulk (Step 5b D1), so `α p = α q = a`.
- `δ α (e₁ c) = α c + α p = α c + a`. RHS `= (∑_{c'∈S} e₁(c')) (e₁ c) = 1`
  (the term `c'=c` contributes 1; all others 0 by Step 5b D2 distinctness).
  ⇒ `α c = 1 + a`.
- `δ α (e₂ c) = α c + α q = (1+a) + a = 1`. But
  `(∑_{c'∈S} e₁(c')) (e₂ c) = 0` because `e₂ c ∉ {e₁(c') : c'∈S}` (Step 5b D2).
  ⇒ `1 = 0` in F₂. Contradiction.

This closes `∑_{c∈S} r c = 0 → False`, hence `LinearIndependent` via
`Fintype.linearIndependent_iff`.

**n-generality trace (gap check).** Every quantity above is a uniform function
of `n`: `bulk_connected_general` is `∀ n≥6`; D1/D2 are `omega`-provable from the
explicit corner-neighbour coordinates for all `n≥6`; the support computation
`(∑ e₁(c)) e ∈ {0,1}` uses only distinctness, not the value of `n`. The argument
is **NOT** n=6-specific. The only place `n=6` enters elsewhere is the
`native_decide` base cases `bulk_conn_6/7` *inside* `bulk_connected_general` —
already discharged. **No hidden assumption.**

### Hard Core #2 — Step 5b distinctness, and the “exactly two incident vertices”

**What makes it hard.** Not conceptually hard — it is the *volume* of explicit
`Fin n` coordinate arithmetic. Eight mandatory edges, their endpoints, pairwise
distinctness, and “neighbours are non-corner” all for symbolic `n≥6`. Mathlib
gives no shortcut; it is `omega` after unfolding, repeated ~28×. Risk is not
soundness but *boilerplate fatigue* / `Fin` casting friction (`Fin.val`,
`Nat.sub`, `n-1`, `n-3` underflow for small `n` — but guarded by `n≥6`).

**Settled?** Yes — finite explicit geometry, the same kind already handled in
`BulkConnectivity.lean`’s `witness_of_target` (which does exactly this style of
`omega`-on-`KAdj` reasoning). Reuse that file’s idioms (`nmul21/nmul12`,
`isCorner` unfolding) verbatim.

**Mitigation.** Define the corner-neighbour pairs as an explicit function
`cornerEdges : Corners → Edge n × Edge n` and prove the 8-distinctness as a
single `Finset.card = 8` (`decide` won’t work symbolically; use
`Finset.card_eq_of_injOn`-style with `omega` per pair). Budget this as its own
sub-file, comparable in size to one layer of `BulkConnectivity.lean`.

---

## Cascade (Part D)

### Theorem U (hybrid torus→plane graph, `TheoremU.lean`)

The paper’s Theorem U is the meta-cycle / single-orbit decomposition argument
(`thm` around `.tex:1046`), conceptually downstream of `Q(n)=3` but **not a
re-run of the same quotient machinery** — it is about counting orbits of a
permutation `σ` over `k²` blocks (`.tex:1063`). It does *not* reuse Steps 1–9.
**Verdict: needs-new-work** (different combinatorial object); Q(n)=3 is a *named
input* to it, not a template. Out of scope for this blueprint; scope separately.

### Q(n,m)=3 rectangular general (`RectangularQ3.lean`)

This **is** mostly the same blueprint. Deltas:
- Vertex set `Fin n × Fin m`, corners still the 4 extremal cells, Lemma 2.6 still
  gives degree 2 for `n,m ≥ 4` (same knight-offset argument).
- Steps 1–4, 6–9 are **verbatim reusable** (replace `n` by the pair `(n,m)`;
  `δ`, `R`, quotient, `r_c`, σ-free endgame are dimension-agnostic).
- Step 5 needs a **rectangular bulk connectivity** theorem
  `Bulk(n,m).Connected`. The current `bulk_connected_general` is square-only.
  This is the *one* genuinely new ingredient: extend the `n→n+2` induction to a
  2-parameter induction (`(n,m)→(n+2,m)` and `(n,m)→(n,m+2)`) with rectangular
  base cases. The translation-embedding core (`emb_adj`) is already
  parameter-symmetric and reuses directly; only the witness tables and base
  cases grow.
- Step 5b distinctness re-proved with `(n,m)` corners (same `omega` style).

**Verdict: mostly-reuse.** New work = rectangular bulk connectivity
(~one `BulkConnectivity.lean`-sized effort) + re-parameterising Step 5b. The
abstract Steps 1–4,6–9 carry over with a `(n,m)` substitution.

---

## Mathlib dependency table

| # | Lemma / def (exact Mathlib name) | Where used | Source (v4.16.0) | Exists |
|---|----------------------------------|-----------|------------------|:--:|
| 1 | `Module.finrank_fintype_fun_eq_card` | Step 0,8 | `LinearAlgebra/Dimension/Constructions.lean:305` | Y |
| 2 | `ZMod.instField` | Step 0 (Field F₂) | `Data/ZMod/Basic.lean:1062` | Y |
| 3 | `Nat.fact_prime_two` (`Fact (Prime 2)`) | Step 0 | `Data/Nat/Prime/Defs.lean:439` | Y |
| 4 | `LinearMap.range` (Submodule) | Step 1 | `LinearAlgebra/Basic` | Y |
| 5 | `Submodule.Quotient` (Module instance) | Step 2 | `LinearAlgebra/Quotient/Basic.lean` | Y |
| 6 | `Submodule.mkQ`, `Submodule.range_mkQ` | Step 2 | `Quotient/Basic.lean:157` | Y |
| 7 | `Submodule.Quotient.eq` / `.mk_eq_zero` | Steps 2,4,5 | `Quotient/Basic.lean` | Y |
| 8 | `Submodule.span`, `Submodule.subset_span`, `span_le` | Steps 3,7 | `LinearAlgebra/Span.lean` | Y |
| 9 | `Pi.single` | Steps 3,4 | `Data/Pi` | Y |
| 10 | `LinearMap.mem_range_self` | Step 4 | `LinearAlgebra/Basic` | Y |
| 11 | `Fintype.linearIndependent_iff` | Steps 5,8 | `LinearAlgebra/LinearIndependent.lean` | Y |
| 12 | `LinearIndependent.map'` (ker = ⊥) | Step 8 | `LinearAlgebra/LinearIndependent.lean:920` | Y |
| 13 | `Submodule.map_span` | Step 6 | `LinearAlgebra/Span.lean` (used `Lie/IdealOperations.lean:204`) | Y |
| 14 | `finrank_span_eq_card` | Step 8 | `LinearAlgebra/Dimension/Constructions.lean:447` | Y |
| 15 | `SimpleGraph.Connected.preconnected` | Step 5 | `Connectivity/Subgraph.lean:58` | Y |
| 16 | `SimpleGraph.Reachable.elim`, `Walk` recursor | Step 5 (new lemma) | `SimpleGraph/Path.lean:745` | Y |
| 17 | `SimpleGraph.Reachable.map` | (already used in bulk proof) | `SimpleGraph/Path.lean:787` | Y |
| 18 | `bulk_connected_general` (project) | Step 5 | `KnightTour/BulkConnectivity.lean:235` | Y |
| 19 | `knightDegree n c = 2` (project, Lemma 2.6) | Steps 3,4 | `KnightTour/Lemma1.lean` | Y |
| — | `Submodule.finrank_quotient_add_finrank` | (σ-route only; **not needed** in chosen σ-free endgame) | `LinearAlgebra/FiniteDimensional.lean` | Y |
| — | `LinearMap.finrank_range_add_finrank_ker` | (σ-route only; **not needed**) | `LinearAlgebra/FiniteDimensional.lean:130` | Y |
| **NEW** | `edgewise_const_of_connected` | Step 5 core | **must be proved** (Walk induction, ~10 lines) | **N** |

All load-bearing Mathlib lemmas for the chosen route **exist**. Exactly **one**
new lemma must be proved from scratch (`edgewise_const_of_connected`), and it is
short and routine.

---

## Risk register

1. **Edge representation choice (Step 0)** ripples through every later step.
   `Sym2 V` is recommended; if the agent picks the `List`-based `undirectedEdges`
   to reuse `computeQ`, it inherits `List.eraseDups`/index reasoning that fights
   the abstract proof. *Mitigation:* commit to `Sym2`/`Finset` incidence up front;
   keep `computeQ` strictly for `native_decide`.
2. **Target restatement (Step 10).** The chosen route proves `Q_abstract n = 3`,
   **not** `computeQ n hn = 3`. Closing the *existing* `sorry` verbatim would
   require the rejected Route-2 bridge. The Lean agent must be authorised to
   restate the public theorem in terms of `Q_abstract` (the paper’s actual Def
   2.4) and demote `computeQ` to a concrete-`n` cross-check. *If forbidden from
   restating, the project balloons by a full verified-Gaussian-elimination
   proof.* **Decision required before Lean work starts.**
3. **Step 5b volume.** ~28 pairwise `omega` distinctness goals over symbolic
   `Fin n`, `n≥6`; `Fin.val`/`Nat.sub` underflow friction. Soundness is not at
   risk; effort is underestimated if treated as “trivial”. Reuse
   `BulkConnectivity.lean` idioms.
4. **`δ α e = α u + α w` packaging.** Getting the coboundary to compute cleanly on
   a `Sym2` edge (summing over the two endpoints without double counting / with
   the right `Finset`) is the fiddliest *routine* step; a bad definition makes
   Steps 4 and 5 painful. *Mitigation:* prove a single clean
   `δ_apply_edge : δ α e = α e.out.1 + α e.out.2` (or via `Sym2.lift`) early and
   use it everywhere.
5. **LI bookkeeping in Step 8** (`{w_k}` LI from `{r_c}` LI). Over F₂ this is a
   subset-parity rearrangement; doing it by hand is error-prone. *Mitigation:*
   prefer the `LinearIndependent.map'` route (inject the F₂³ standard basis via
   `e_k ↦ r c₁ + r c_{k}` with explicit `ker = ⊥`) over a raw
   `Fintype.linearIndependent_iff` Finset sum.
6. **No mathematical gap found.** The σ argument is genuinely general; the
   adopted σ-free endgame is strictly more elementary and removes the only place
   rank–nullity would have been needed.

---

=== Q3 FORMALIZATION BLUEPRINT ===
Route chosen: 1 (abstract-only) — define Q via Submodule/quotient/finrank;
  never verify the imperative `gf2Rank`. Adopt a σ-free endgame (explicit
  3-vector basis `{r_c1+r_ck}`) that removes the only rank–nullity dependence.
Total steps: 10 (Steps 0–9; Step 10 = a target-restatement decision).
Hard steps: Step 5 (independence via bulk connectivity) and its support
  Step 5b (8-mandatory-edge distinctness).
Mathlib lemmas relied on: 19 verified to exist in v4.16.0 (+2 σ-route lemmas
  available but not needed) + 1 new lemma to prove (`edgewise_const_of_connected`,
  ~10-line Walk induction). → 19 verified, 1 new.
Highest risk: the Step 10 target-restatement decision — proving the *existing*
  `computeQ n hn = 3` signature verbatim would force a verified-Gaussian-
  elimination project; the route instead proves `Q_abstract n = 3` (paper Def 2.4)
  and keeps `computeQ` for concrete-n `native_decide`. Must be authorised.
Hidden-assumption flags: none. Step 5 traced to be uniformly general in n
  (bulk_connected_general is ∀n≥6; distinctness is omega for all n≥6); the only
  n=6 input is the already-discharged base case inside bulk connectivity.
Cascade (Thm U, Q(n,m)): Q(n,m)=3 = mostly-reuse (Steps 1–4,6–9 verbatim with
  (n,m); new work = rectangular bulk connectivity + Step 5b reparam). Theorem U =
  needs-new-work (distinct meta-cycle/orbit argument; uses Q(n)=3 as input only).
Estimated Lean effort: comparable to bulk connectivity overall — roughly
  (a) one linear-algebra layer (Steps 1–4,6–9, mostly named lemmas: light),
  (b) one Step-5b distinctness layer (~BulkConnectivity-sized: medium),
  (c) the `edgewise_const_of_connected` lemma + Step 5 assembly (small but the
  conceptual crux). Net: similar magnitude to BulkConnectivity.lean, front-loaded
  on Step 5b boilerplate rather than on any deep proof.

Honest go/no-go: **GO.** The abstract proof is genuinely general and
formalizable on stock Mathlib v4.16.0 plus one short new graph lemma. No
mathematical gap surfaced; the σ functional is general and is in fact replaceable
by a more elementary explicit-basis endgame, de-risking the hardest-sounding
step. The only non-mathematical decision needed is whether the Lean agent may
restate the public target as `Q_abstract n = 3` (recommended) rather than
re-verifying the imperative `computeQ`.
