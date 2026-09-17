/-
  Gate 0 — the REDUCTION theorem.

      rank(Ham n) = β₁(n) − 3   ⟺   Z_bulk(n) ⊆ Span(Ham n)

  This is the spine of the tightness argument: it replaces a statement about
  the span of *all* Hamiltonian cycles (a set nobody can enumerate) by a
  statement about a single fixed, explicit subspace.

  Architecture mirrors QAbstract.lean deliberately: the theorem is proved for
  an ABSTRACT submodule `H` satisfying four hypotheses, and the bridge to the
  combinatorial `Ham` (spanned by actual tour vectors) is a separate layer —
  the same split that `Q_abstract` vs `computeQ` already uses in this project.
  Doing it the other way round is what made the computeQ bridge painful.

  All dimension statements use `beta1_abstract` (counts via `Fintype.card`),
  NOT the computational `Basic.beta1` (counts via `edgeListOriented`); see L0'
  for why, and for the intended bridge between them.

  Layers:
    L0' beta1_abstract — the primary β₁ for this chain
    L0  boundary ∂ : F₂^E → F₂^V, cycle space Z = ker ∂
    L1  Z_bulk = cycles vanishing on the 8 corner-incident (mandatory) edges
    L2  dim Z_bulk = β₁ − 4        (uses bulk_connected_general)
    L3  abstract reduction theorem (SORRY-FREE given L2)
    L4  bridge: the combinatorial Ham satisfies the hypotheses   [open]

  STATUS: L0, L1, L3 complete. L2 and L4 carry documented `sorry`s — see the
  block comment above each. No claim of tightness is made here; this file
  proves only the equivalence, which is what makes tightness attackable.
-/
import KnightTour.QAbstract
import KnightTour.CycleSpaceDim
import KnightTour.PuncGraph
import KnightTour.PuncCount
import Mathlib.LinearAlgebra.Dimension.Finrank
import Mathlib.LinearAlgebra.Finsupp.LinearCombination

namespace KnightTour

open scoped BigOperators

variable {n : ℕ}

/-! ## L0': the first Betti number, abstractly

  `Basic.beta1` computes `β₁` from `edgeListOriented`, a `List`-based
  definition. Everything in this file is stated instead against
  `beta1_abstract`, which reads the counts off `Fintype.card`.

  **`beta1_abstract` is the PRIMARY definition for the reduction chain.**
  `Basic.beta1` remains the computational one, used by the `native_decide`
  layers. Naming follows the `Q_abstract` / `computeQ` convention already
  established in this project: the `_abstract` form carries the proofs, the
  computational form carries the evaluation, and they are bridged only where
  a bridge is actually needed.

  The bridge is deliberately NOT stated as a `sorry`d theorem — an unproved
  claim sitting in the namespace contaminates silently, the same reason
  `tightness` is a comment in L4 rather than a declaration. Intended signature,
  for whoever needs it:

      theorem beta1_abstract_eq_beta1 (n : ℕ) (hn : 4 ≤ n) :
          beta1_abstract n = beta1 n

  Proving it means relating `(edgeListOriented n).length / 2` to
  `Fintype.card (KEdge n)` for symbolic `n` — the same shape of bridge that
  made `computeQ` painful. Nothing downstream of `reduction` needs it. -/

/-- `β₁ = |E| − |V| + 1`, read off `Fintype.card`. Primary for this file. -/
noncomputable def beta1_abstract (n : ℕ) : Int :=
  (Fintype.card (KEdge n) : Int) - (Fintype.card (Square n) : Int) + 1

/-! ## L0: boundary map and cycle space -/

/-- The boundary `∂ : F₂^E → F₂^V`, `∂(x)(v) = Σ_{e ∋ v} x e`.
    Dual to `coboundary`; the cycle space is its kernel. -/
noncomputable def boundary (n : ℕ) : (KEdge n → F) →ₗ[F] (Square n → F) where
  toFun x v := ∑ e : KEdge n, if v ∈ e.val then x e else 0
  map_add' x y := by
    funext v
    simp only [Pi.add_apply]
    rw [← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl fun e _ => ?_
    by_cases h : v ∈ e.val <;> simp [h]
  map_smul' c x := by
    funext v
    simp only [Pi.smul_apply, RingHom.id_apply, smul_eq_mul, Finset.mul_sum]
    refine Finset.sum_congr rfl fun e _ => ?_
    by_cases h : v ∈ e.val <;> simp [h]

/-- The cycle space `Z₁ = ker ∂`.

    Defined as `CycleSpaceDim.cycleSpace` of the knight graph rather than via
    the local `boundary` above, so that the general dimension formula and the
    extension-by-zero transfer apply verbatim. `boundary` is kept only as the
    readable spelling of what `bd` computes. -/
noncomputable def Zcyc (n : ℕ) : Submodule F (KEdge n → F) :=
  cycleSpace (KnightGraph n)

/-! ## L1: the bulk cycle space -/

/-- Vectors vanishing on a prescribed set of coordinates. -/
def vanishOn (n : ℕ) (S : Set (KEdge n)) : Submodule F (KEdge n → F) where
  carrier := {x | ∀ e ∈ S, x e = 0}
  zero_mem' := by intro e _; rfl
  add_mem' := by
    intro a b ha hb e he
    simp [ha e he, hb e he]
  smul_mem' := by
    intro c a ha e he
    simp [ha e he]

/-- `Z_bulk` — the cycle space of `G_n` minus the four corners, realised inside
    `F₂^E` as the cycles supported away from the corner-incident edges. -/
noncomputable def ZBulk (n : ℕ) : Submodule F (KEdge n → F) :=
  Zcyc n ⊓ vanishOn n {e | IsCornerEdge n e}

theorem zbulk_le_zcyc (n : ℕ) : ZBulk n ≤ Zcyc n := inf_le_left

theorem zbulk_vanishes {n : ℕ} {x : KEdge n → F} (hx : x ∈ ZBulk n)
    {e : KEdge n} (he : IsCornerEdge n e) : x e = 0 :=
  hx.2 e he

/-! ## L2: `dim Z_bulk = β₁ − 4`

  `Z_bulk` is the cycle space of the graph `Bulk n`, which has `n² − 4`
  vertices and `|E| − 8` edges (each corner has degree exactly 2, and the four
  corners are pairwise non-adjacent, so no edge is double-counted).
  `bulk_connected_general` gives connectivity for `n ≥ 6`, hence

      dim Z_bulk = (|E| − 8) − (n² − 4) + 1 = β₁ − 4.

  STATUS of the three ingredients:

  (a) `dim_{𝔽₂}(cycle space) = |E| − |V| + #components` for ANY finite simple
      graph — **DONE, sorry-free**: `KnightTour.finrank_cycleSpace_general` in
      `CycleSpaceDim.lean`. This was the only piece missing from Mathlib; it is
      stated for an arbitrary `SimpleGraph`, independent of knights, and the
      general (#components) form is what `Z_bulk` actually needs, since the
      punctured graph has FIVE components — the bulk plus four isolated
      corners.

  (c) the identification of `ZBulk n` with a cycle space living in `𝔽₂^{E(G)}`
      — **DONE, sorry-free**: `cycleSpaceOn` / `finrank_cycleSpaceOn_eq` in
      `CycleSpaceDim.lean`, via extension by zero along `G' ≤ G`.

  (b) the knight-specific counting — **THE ONLY THING LEFT**:
        · `ZBulk n = cycleSpaceOn (KnightGraph n) (Punc n)`, where `Punc n` is
          the spanning subgraph with corner-incident edges deleted; this is the
          statement that "vanishes on corner edges" and "is in the image of
          extension by zero" agree;
        · `|E(Punc n)| = |E(G_n)| − 8`, from `corner_adj_iff` (QAbstract, gives
          corner degree exactly 2) plus the four corners being pairwise
          non-adjacent, so the 8 edges are distinct;
        · `#components(Punc n) = 5`, from `bulk_connected_general` plus the
          four corners being isolated.

      Then `dim Z_bulk = (|E| − 8) − n² + 5 = β₁ − 4` by arithmetic.

      (b.3) IS DONE, sorry-free: `card_connectedComponent_punc` in
      `PuncGraph.lean` proves `#components(Punc n) = 5` for `n ≥ 6`, via an
      explicit `Equiv` with `Fin 5` rather than any quotient cardinality
      computation. What is left of (b) is (b.1) the identification and (b.2)
      the edge count.

  So the mathematics is done. The part that was NOT mere bookkeeping — the
  component count over a quotient type — is now also done; what remains, (b.1)
  and (b.2), genuinely is bookkeeping. Note for the record: calling all of (b)
  "mechanical" in an earlier pass was wrong, and (b.3) was the reason. -/
/-- **(b.1)** `Z_bulk` is exactly the cycle space of the punctured graph, seen
    inside `𝔽₂^{E(G_n)}`: "is a cycle and vanishes on corner edges" and "is an
    extension by zero of a cycle of `Punc n`" are the same condition. -/
theorem ZBulk_eq (n : ℕ) :
    ZBulk n = cycleSpaceOn (KnightGraph n) (Punc n) (punc_le n) := by
  apply le_antisymm
  · rintro x ⟨hxZ, hxV⟩
    have hz : bd (KnightGraph n) x = 0 := hxZ
    set f : (Punc n).edgeSet → F :=
      fun e => x ⟨e.val, SimpleGraph.edgeSet_mono (punc_le n) e.property⟩ with hfdef
    have hext : extEdge (KnightGraph n) (Punc n) (punc_le n) f = x := by
      funext e
      show (if h' : e.val ∈ (Punc n).edgeSet then f ⟨e.val, h'⟩ else 0) = x e
      by_cases h : e.val ∈ (Punc n).edgeSet
      · rw [dif_pos h]
      · rw [dif_neg h]
        symm
        apply hxV
        by_contra hc
        exact h (mem_punc_edgeSet_iff.mpr ⟨e.property, hc⟩)
    refine ⟨f, ?_, hext⟩
    show bd (Punc n) f = 0
    funext v
    rw [← bd_extEdge (KnightGraph n) (Punc n) (punc_le n) f v, hext]
    exact congrFun hz v
  · rintro x ⟨f, hf, rfl⟩
    refine ⟨cycleSpaceOn_le_cycleSpace (KnightGraph n) (Punc n) (punc_le n)
              ⟨f, hf, rfl⟩, ?_⟩
    intro e he
    have hnot : e.val ∉ (Punc n).edgeSet := by
      intro hmem
      exact (mem_punc_edgeSet_iff.mp hmem).choose_spec he
    show (if h' : e.val ∈ (Punc n).edgeSet then f ⟨e.val, h'⟩ else 0) = 0
    rw [dif_neg hnot]

theorem dim_ZBulk (n : ℕ) (hn : 6 ≤ n) :
    (Module.finrank F (ZBulk n) : Int) = beta1_abstract n - 4 := by
  classical
  letI : Fintype (Punc n).ConnectedComponent := puncCompFintype hn
  have hE := card_punc_edgeSet (n := n) hn
  have hC : Fintype.card (Punc n).ConnectedComponent = 5 :=
    card_connectedComponent_punc hn
  -- transport the dimension along (b.1), keeping ONE syntactic form of the
  -- `cycleSpaceOn` term so that `omega` sees a single atom
  have hEq : Module.finrank F (ZBulk n)
      = Module.finrank F2 (cycleSpaceOn (KnightGraph n) (Punc n) (punc_le n)) := by
    rw [ZBulk_eq n]
  have hdim := finrank_cycleSpaceOn_eq (KnightGraph n) (Punc n) (punc_le n)
  rw [hC] at hdim
  rw [hEq]
  unfold beta1_abstract
  omega

/-! ## L3: the abstract reduction theorem

  Stated for an abstract submodule `H` so that the combinatorial content
  (what a tour *is*) stays in L4. The four hypotheses are exactly what the
  set of Hamiltonian cycles supplies. -/

section Abstract

variable (H : Submodule F (KEdge n → F))

/-- A tour vector is 1 on every corner-incident edge: corners have degree 2,
    so both their edges are forced. This is what separates `H` from `Z_bulk`. -/
def HitsMandatory (τ : KEdge n → F) : Prop :=
  ∀ e : KEdge n, IsCornerEdge n e → τ e = 1

/-- **Key separation.** A vector hitting the mandatory edges is not in `Z_bulk`,
    because `Z_bulk` vanishes there. Needs at least one corner edge to exist. -/
theorem not_mem_ZBulk_of_hitsMandatory {τ : KEdge n → F}
    (hτ : HitsMandatory τ) (e₀ : KEdge n) (he₀ : IsCornerEdge n e₀) :
    τ ∉ ZBulk n := by
  intro hmem
  have h1 : τ e₀ = 1 := hτ e₀ he₀
  have h0 : τ e₀ = 0 := zbulk_vanishes hmem he₀
  rw [h1] at h0
  exact one_ne_zero h0

/-- `Z_bulk ⊓ span{τ} = ⊥` when `τ` hits the mandatory edges.
    Over `F₂` the span is `{0, τ}`, so this is immediate from the separation. -/
theorem inf_span_eq_bot {τ : KEdge n → F}
    (hτ : HitsMandatory τ) (e₀ : KEdge n) (he₀ : IsCornerEdge n e₀) :
    ZBulk n ⊓ Submodule.span F {τ} = ⊥ := by
  rw [Submodule.eq_bot_iff]
  rintro x ⟨hxZ, hxS⟩
  rw [SetLike.mem_coe, Submodule.mem_span_singleton] at hxS
  rw [SetLike.mem_coe] at hxZ
  obtain ⟨c, rfl⟩ := hxS
  -- over F₂ either c = 0 (done) or c = 1 (contradicts the separation)
  have hF2 : ∀ d : F, d = 0 ∨ d = 1 := by decide
  rcases hF2 c with rfl | rfl
  · simp
  · exfalso
    rw [one_smul] at hxZ
    exact not_mem_ZBulk_of_hitsMandatory hτ e₀ he₀ hxZ

/-- **THE REDUCTION THEOREM.**  If a submodule `H` of the edge space
    (i) contains `Z_bulk`, (ii) contains some vector hitting every mandatory
    edge, and (iii) has rank at most `β₁ − 3`, then its rank is exactly
    `β₁ − 3`.

    Hypothesis (iii) is `Q(n) = 3`, already proved sorry-free in this project
    (`Q_abstract_eq_three`). Hypothesis (ii) holds for any single tour.
    Hypothesis (i) is the whole content of the tightness problem. -/
theorem reduction (hn : 6 ≤ n)
    (hZ : ZBulk n ≤ H)
    {τ : KEdge n → F} (hτH : τ ∈ H) (hτ : HitsMandatory τ)
    (e₀ : KEdge n) (he₀ : IsCornerEdge n e₀)
    (hub : (Module.finrank F H : Int) ≤ beta1_abstract n - 3) :
    (Module.finrank F H : Int) = beta1_abstract n - 3 := by
  refine le_antisymm hub ?_
  -- lower bound: Z_bulk ⊔ span{τ} ≤ H, and that sup has rank (β₁ − 4) + 1
  have hsub : ZBulk n ⊔ Submodule.span F {τ} ≤ H :=
    sup_le hZ (Submodule.span_le.mpr (Set.singleton_subset_iff.mpr hτH))
  have hinf : ZBulk n ⊓ Submodule.span F {τ} = ⊥ := inf_span_eq_bot hτ e₀ he₀
  have hspan : Module.finrank F (Submodule.span F ({τ} : Set (KEdge n → F))) = 1 := by
    apply finrank_span_singleton
    intro h
    subst h
    have : (0 : KEdge n → F) e₀ = 1 := hτ e₀ he₀
    simp at this
  have hsup : Module.finrank F (ZBulk n ⊔ Submodule.span F {τ} : Submodule F _)
      = Module.finrank F (ZBulk n) + 1 := by
    have := Submodule.finrank_sup_add_finrank_inf_eq
      (ZBulk n) (Submodule.span F ({τ} : Set (KEdge n → F)))
    rw [hinf] at this
    simp [hspan, finrank_bot] at this
    omega
  have hmono := Submodule.finrank_mono hsub
  have hdim := dim_ZBulk n hn
  omega

end Abstract

/-! ## L4: bridge to the combinatorial `Ham`   [OPEN]

  What remains is to define `Ham n` as the span of the indicator vectors of
  actual Hamiltonian cycles of `KnightGraph n`, and to discharge the three
  hypotheses of `reduction`:

    (ii) `HitsMandatory τ` for a tour τ — from `corner` having degree exactly 2
         in `KnightGraph n`, so both incident edges are forced. This is the
         easy one and is essentially `corner_adj_iff` from QAbstract.
    (iii) `finrank Ham ≤ β₁ − 3` — this is `Q_abstract_eq_three`, but stated
         against `Q_abstract`; the bridge is the same shape as the existing
         `computeQ` shadow lemmas.
    (i)  `ZBulk n ≤ Ham n` — OPEN. This is (U1)+(U2) and is exactly what the
         `tightness_witnesses/` experiments are attacking. Verified by explicit
         construction for n ∈ {8,10,12}; not proved for general n.

  Deliberately NOT stated as a `sorry`d theorem here: writing
  `theorem tightness : ... := sorry` would put an unproved headline claim in
  the namespace, and this project has been bitten by exactly that before. -/

end KnightTour
