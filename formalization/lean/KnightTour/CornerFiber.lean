/-
  # The corner functionals and the fibre `y = 1`

  A second route to the upper bound `rk(Ham) ≤ β₁ − 3`, phrased through four
  linear functionals on the cycle space instead of through the quotient
  `𝔽₂^E / row(∂₁)`.

  Setup. Each corner `c_k` has degree exactly 2 (`corner_adj_iff`), so a cycle
  `z ∈ Z₁` uses BOTH of its mandatory edges or NEITHER. Hence

      cornerVal k (z) := z (mandEdge1 k)

  already determines `z` on both edges of `c_k` (`cycle_corner_eq`), and the
  joint fibre `y = 0` is exactly `Z_bulk` (`mem_ZBulk_iff_cornerVal_zero`).
  Every tour lies in the fibre `y = 1`, which is an AFFINE subspace not
  containing `0`; so its linear span exceeds `Z_bulk` by at most one dimension:

      span(fibre 1)  ≤  Z_bulk ⊔ span{τ}     (`span_le_sup`)

  and `dim (Z_bulk ⊔ span{τ}) = (β₁ − 4) + 1 = β₁ − 3` (`finrank_le_of_fiberOne`).

  ## What this does and does not add

  It does NOT add mathematical content. `Q_abstract_eq_three` (QAbstract.lean)
  already gives the same bound, sorry-free. What is gained is that the bound
  here is reached WITHOUT the quotient `𝔽₂^E / row(∂₁)` — only the degree-2
  fact and `dim_ZBulk` — which is a shorter and more readable route.

  It also does not close the tightness problem. The matching LOWER bound still
  requires `Z_bulk ≤ span(Ham)`, exactly as in `Reduction.reduction`; nothing
  here supplies it. Boards such as `5 × 8` have connected bulk and `Q = 3` yet
  deficit 9 — the fibre picture is equally valid there and equally silent.

  Following project convention, no `theorem ... := sorry` appears here.
-/
import KnightTour.Reduction

namespace KnightTour

open Finset

variable {n : ℕ}

/-! ## The two edges at a corner -/

/-- The two mandatory edges of one corner are distinct. -/
theorem mandEdge1_ne_mandEdge2 (hn : 6 ≤ n) (k : Fin 4) :
    mandEdge1 n hn k ≠ mandEdge2 n hn k := by
  intro h
  have : ((k, true) : Fin 4 × Bool) = (k, false) :=
    mandE_injective hn (by simpa [mandE] using h)
  simp at this

/-- **Degree 2, edgewise.** An edge contains corner `k` iff it is one of that
    corner's two mandatory edges. This is `corner_adj_iff` transported from
    vertices to edges. -/
theorem mem_corner_iff (hn : 6 ≤ n) (k : Fin 4) (e : KEdge n) :
    corner n (by omega) k ∈ e.val ↔
      e = mandEdge1 n hn k ∨ e = mandEdge2 n hn k := by
  constructor
  · intro hmem
    obtain ⟨w, hw⟩ := Sym2.mem_iff_exists.mp hmem
    have hadj : (KnightGraph n).Adj (corner n (by omega) k) w := by
      have := e.property
      rw [hw] at this
      exact (KnightGraph n).mem_edgeSet.mp this
    rcases (corner_adj_iff n hn k w).mp hadj with rfl | rfl
    · exact Or.inl (Subtype.ext (by simp [mandEdge1, mkEdge, hw]))
    · exact Or.inr (Subtype.ext (by simp [mandEdge2, mkEdge, hw]))
  · rintro (rfl | rfl) <;> simp [mandEdge1, mandEdge2, mkEdge]

/-- The corner-incident edges of corner `k`, as a `Finset`. -/
theorem filter_mem_corner (hn : 6 ≤ n) (k : Fin 4) :
    Finset.univ.filter
        (fun e : KEdge n => corner n (by omega) k ∈ e.val) =
      ({mandEdge1 n hn k, mandEdge2 n hn k} : Finset (KEdge n)) := by
  classical
  ext e
  simp [mem_corner_iff hn k e]

/-! ## The boundary at a corner -/

theorem bd_apply (x : KEdge n → F) (v : Square n) :
    bd (KnightGraph n) x v = ∑ e : KEdge n, if v ∈ e.val then x e else 0 := by
  simp only [bd, Matrix.mulVecLin_apply, Matrix.mulVec, dotProduct, inc]
  refine Finset.sum_congr rfl fun e _ => ?_
  by_cases h : v ∈ e.val <;> simp [h]

/-- The boundary at a corner reads off exactly the two mandatory coordinates. -/
theorem bd_at_corner (hn : 6 ≤ n) (k : Fin 4) (x : KEdge n → F) :
    bd (KnightGraph n) x (corner n (by omega) k)
      = x (mandEdge1 n hn k) + x (mandEdge2 n hn k) := by
  classical
  rw [bd_apply, ← Finset.sum_filter, filter_mem_corner hn k]
  exact Finset.sum_pair (mandEdge1_ne_mandEdge2 hn k)

/-! ## The corner functionals -/

/-- `cornerVal k : 𝔽₂^E →ₗ 𝔽₂`, evaluation at the corner's first mandatory
    edge. On `Z₁` it does not matter which of the two is chosen
    (`cycle_corner_eq`). -/
noncomputable def cornerVal (n : ℕ) (hn : 6 ≤ n) (k : Fin 4) :
    (KEdge n → F) →ₗ[F] F :=
  LinearMap.proj (mandEdge1 n hn k)

@[simp] theorem cornerVal_apply (hn : 6 ≤ n) (k : Fin 4) (x : KEdge n → F) :
    cornerVal n hn k x = x (mandEdge1 n hn k) := rfl

/-- **`y_k` is well defined.** A cycle takes the same value on both edges of a
    corner: either it uses both or it uses neither. -/
theorem cycle_corner_eq (hn : 6 ≤ n) (k : Fin 4) {x : KEdge n → F}
    (hx : x ∈ Zcyc n) :
    x (mandEdge1 n hn k) = x (mandEdge2 n hn k) := by
  have h0 : bd (KnightGraph n) x = 0 := hx
  have := bd_at_corner hn k x
  rw [h0] at this
  have hsum : x (mandEdge1 n hn k) + x (mandEdge2 n hn k) = 0 := this.symm
  -- over 𝔽₂, `a + b = 0` is `a = b`
  have := congrArg (· + x (mandEdge2 n hn k)) hsum
  simpa [add_assoc, CharTwo.add_self_eq_zero] using this

/-! ## The fibres -/

/-- The joint fibre `y = 0` inside `Z₁` is exactly `Z_bulk`. -/
theorem mem_ZBulk_iff_cornerVal_zero (hn : 6 ≤ n) {x : KEdge n → F}
    (hx : x ∈ Zcyc n) :
    x ∈ ZBulk n ↔ ∀ k : Fin 4, cornerVal n hn k x = 0 := by
  constructor
  · intro hmem k
    exact zbulk_vanishes hmem (isCornerEdge_mandE hn (k, true))
  · intro hzero
    refine ⟨hx, ?_⟩
    intro e he
    obtain ⟨p, rfl⟩ := exists_mandE_of_isCornerEdge hn he
    cases hb : p.2
    · -- second mandatory edge: equal to the first by `cycle_corner_eq`
      have : mandE n hn p = mandEdge2 n hn p.1 := by simp [mandE, hb]
      rw [this, ← cycle_corner_eq hn p.1 hx]
      exact hzero p.1
    · have : mandE n hn p = mandEdge1 n hn p.1 := by simp [mandE, hb]
      rw [this]
      exact hzero p.1

/-- The fibre `y = 1`: every corner functional evaluates to 1. For cycles this
    is equivalent to `HitsMandatory` (`hitsMandatory_of_fiberOne`). -/
def CornerFiberOne (n : ℕ) (hn : 6 ≤ n) (x : KEdge n → F) : Prop :=
  ∀ k : Fin 4, cornerVal n hn k x = 1

/-- On cycles, the fibre `y = 1` is `HitsMandatory`: a cycle through all four
    corners uses all eight mandatory edges. -/
theorem hitsMandatory_of_fiberOne (hn : 6 ≤ n) {x : KEdge n → F}
    (hx : x ∈ Zcyc n) (h1 : CornerFiberOne n hn x) : HitsMandatory x := by
  intro e he
  obtain ⟨p, rfl⟩ := exists_mandE_of_isCornerEdge hn he
  cases hb : p.2
  · have : mandE n hn p = mandEdge2 n hn p.1 := by simp [mandE, hb]
    rw [this, ← cycle_corner_eq hn p.1 hx]
    exact h1 p.1
  · have : mandE n hn p = mandEdge1 n hn p.1 := by simp [mandE, hb]
    rw [this]
    exact h1 p.1

/-- **The fibre is affine.** The difference of two elements of the fibre `y = 1`
    lies in `Z_bulk`. (Over `𝔽₂`, difference = sum.) -/
theorem add_mem_ZBulk_of_fiberOne (hn : 6 ≤ n) {x y : KEdge n → F}
    (hx : x ∈ Zcyc n) (hy : y ∈ Zcyc n)
    (hx1 : CornerFiberOne n hn x) (hy1 : CornerFiberOne n hn y) :
    x + y ∈ ZBulk n := by
  have hxy : x + y ∈ Zcyc n := Submodule.add_mem _ hx hy
  refine (mem_ZBulk_iff_cornerVal_zero hn hxy).mpr ?_
  intro k
  rw [map_add, hx1 k, hy1 k]
  decide

/-! ## The upper bound -/

/-- **The affine-span step.** Any set inside the fibre `y = 1` spans no more
    than `Z_bulk` plus one dimension, witnessed by any single member `τ`. -/
theorem span_le_sup (hn : 6 ≤ n) {S : Set (KEdge n → F)}
    (hSZ : ∀ x ∈ S, x ∈ Zcyc n) (hS1 : ∀ x ∈ S, CornerFiberOne n hn x)
    {τ : KEdge n → F} (hτS : τ ∈ S) :
    Submodule.span F S ≤ ZBulk n ⊔ Submodule.span F {τ} := by
  rw [Submodule.span_le]
  intro x hx
  -- x = (x + τ) + τ, with x + τ ∈ Z_bulk
  have hmem : x + τ ∈ ZBulk n :=
    add_mem_ZBulk_of_fiberOne hn (hSZ x hx) (hSZ τ hτS) (hS1 x hx) (hS1 τ hτS)
  have hτmem : τ ∈ Submodule.span F ({τ} : Set (KEdge n → F)) :=
    Submodule.mem_span_singleton_self τ
  have hback : (x + τ) + τ = x := by
    ext e
    simp only [Pi.add_apply]
    have : τ e + τ e = 0 := CharTwo.add_self_eq_zero _
    rw [add_assoc, this, add_zero]
  have hmem2 : (x + τ) + τ ∈ ZBulk n ⊔ Submodule.span F ({τ} : Set (KEdge n → F)) :=
    Submodule.add_mem _ (Submodule.mem_sup_left hmem)
      (Submodule.mem_sup_right hτmem)
  rwa [hback] at hmem2

/-- **THE UPPER BOUND, by the fibre route.** A set of cycles all lying in the
    corner fibre `y = 1` — in particular any set of Hamiltonian tours — spans a
    subspace of dimension at most `β₁ − 3`.

    This is the same bound as `Q_abstract_eq_three` supplies, obtained here
    without the quotient `𝔽₂^E / row(∂₁)`: only corner degree 2 and
    `dim_ZBulk` are used. -/
theorem finrank_le_of_fiberOne (hn : 6 ≤ n) {S : Set (KEdge n → F)}
    (hSZ : ∀ x ∈ S, x ∈ Zcyc n) (hS1 : ∀ x ∈ S, CornerFiberOne n hn x)
    {τ : KEdge n → F} (hτS : τ ∈ S) :
    (Module.finrank F (Submodule.span F S) : Int) ≤ beta1_abstract n - 3 := by
  classical
  have hτZ : τ ∈ Zcyc n := hSZ τ hτS
  have hτ1 : HitsMandatory τ := hitsMandatory_of_fiberOne hn hτZ (hS1 τ hτS)
  have he₀ : IsCornerEdge n (mandE n hn (0, true)) := isCornerEdge_mandE hn _
  -- `Z_bulk ⊓ span{τ} = ⊥`, so the sup has dimension (β₁ − 4) + 1
  have hinf : ZBulk n ⊓ Submodule.span F {τ} = ⊥ :=
    inf_span_eq_bot hτ1 _ he₀
  have hτ0 : τ ≠ 0 := by
    intro h
    have h1 : τ (mandEdge1 n hn 0) = 1 := hS1 τ hτS 0
    rw [h] at h1
    simp at h1
  have hspan : Module.finrank F (Submodule.span F ({τ} : Set (KEdge n → F))) = 1 :=
    finrank_span_singleton hτ0
  have hsup : Module.finrank F (ZBulk n ⊔ Submodule.span F {τ} : Submodule F _)
      = Module.finrank F (ZBulk n) + 1 := by
    have h := Submodule.finrank_sup_add_finrank_inf_eq
      (ZBulk n) (Submodule.span F ({τ} : Set (KEdge n → F)))
    rw [hinf] at h
    simp [hspan, finrank_bot] at h
    omega
  have hle := Submodule.finrank_mono (span_le_sup hn hSZ hS1 hτS)
  have hZB := dim_ZBulk n hn
  omega

end KnightTour
