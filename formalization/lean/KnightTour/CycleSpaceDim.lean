/-
  The dimension of the cycle space of a finite connected graph over 𝔽₂.

      dim_{𝔽₂} ker ∂  =  |E| − |V| + 1        (G connected)

  This is the one genuinely missing ingredient for the reduction theorem
  (`KnightTour.Reduction.dim_ZBulk`), and it is not about knights at all — it
  is a general statement about finite simple graphs, stated here for an
  arbitrary `SimpleGraph V` with `V` finite. Upstreamable as-is.

  Proof strategy (the trick is to route everything through
  `Matrix.rank_transpose` instead of proving `range ∂ = {f | Σ f = 0}` by hand):

    1. incidence matrix `inc : Matrix V E 𝔽₂`, `inc v e = [v ∈ e]`
    2. `∂ = inc.mulVecLin`, `δ = incᵀ.mulVecLin`; note `δ α e = α u + α v`
    3. `ker δ = constants`, hence `finrank (ker δ) = 1`  (connectivity)
       — this is where `edgewise_const_of_connected` enters
    4. rank-nullity on `δ` gives `rank δ = |V| − 1`; `Matrix.rank_transpose`
       gives `rank ∂ = rank δ`; rank-nullity on `∂` finishes.
-/
import KnightTour.EdgewiseConst
import Mathlib.Data.Matrix.Rank
import Mathlib.LinearAlgebra.FiniteDimensional
import Mathlib.Combinatorics.SimpleGraph.Finite
import Mathlib.Data.ZMod.Basic

namespace KnightTour

open Matrix Module

/-- The two-element field. -/
abbrev F2 := ZMod 2

variable {V : Type*} [Fintype V] [DecidableEq V]
variable (G : SimpleGraph V) [DecidableRel G.Adj] [Fintype G.edgeSet]

/-- The incidence matrix of `G` over `𝔽₂`: rows indexed by vertices, columns by
    edges, entry 1 iff the vertex lies on the edge. -/
noncomputable def inc : Matrix V G.edgeSet F2 :=
  fun v e => if v ∈ e.val then 1 else 0

/-- The boundary `∂ : 𝔽₂^E → 𝔽₂^V`. -/
noncomputable def bd : (G.edgeSet → F2) →ₗ[F2] (V → F2) := (inc G).mulVecLin

/-- The coboundary `δ : 𝔽₂^V → 𝔽₂^E`. -/
noncomputable def cbd : (V → F2) →ₗ[F2] (G.edgeSet → F2) := (inc G)ᵀ.mulVecLin

/-- The cycle space `Z₁ = ker ∂`. -/
noncomputable def cycleSpace : Submodule F2 (G.edgeSet → F2) :=
  LinearMap.ker (bd G)

/-! ### Step 3: `ker δ` is exactly the constants -/

theorem cbd_apply (α : V → F2) (e : G.edgeSet) :
    cbd G α e = ∑ v : V, if v ∈ e.val then α v else 0 := by
  simp only [cbd, Matrix.mulVecLin_apply, Matrix.mulVec, dotProduct, inc,
             Matrix.transpose_apply]
  refine Finset.sum_congr rfl fun v _ => ?_
  by_cases h : v ∈ e.val <;> simp [h]

/-- For an edge `s(u,w)` the sum collapses to `α u + α w`. -/
theorem cbd_apply_pair {u w : V} (huw : G.Adj u w) (α : V → F2) :
    cbd G α ⟨s(u, w), G.mem_edgeSet.mpr huw⟩ = α u + α w := by
  have hne : u ≠ w := G.ne_of_adj huw
  rw [cbd_apply]
  rw [Finset.sum_eq_add_of_mem u w (Finset.mem_univ u) (Finset.mem_univ w) hne]
  · simp
  · intro x _ hx
    obtain ⟨h1, h2⟩ := hx
    have hx' : x ∉ (s(u, w) : Sym2 V) := by
      simp only [Sym2.mem_iff]
      push_neg
      exact ⟨h1, h2⟩
    simp [hx']

private theorem f2_add_self (x : F2) : x + x = 0 := by revert x; decide

private theorem f2_of_add_eq_zero {x y : F2} (h : x + y = 0) : x = y := by
  revert h; revert x y; decide

/-- **Kernel of the coboundary = constant functions**, for `G` connected. -/
theorem mem_ker_cbd_iff (hG : G.Connected) (α : V → F2) :
    α ∈ LinearMap.ker (cbd G) ↔ ∀ a b, α a = α b := by
  constructor
  · intro hα
    have hz : cbd G α = 0 := LinearMap.mem_ker.mp hα
    apply edgewise_const_of_connected hG
    intro u w hadj
    have h0 : cbd G α ⟨s(u, w), G.mem_edgeSet.mpr hadj⟩ = 0 := by
      rw [hz]; rfl
    rw [cbd_apply_pair G hadj] at h0
    exact f2_of_add_eq_zero h0
  · intro hconst
    rw [LinearMap.mem_ker]
    funext e
    obtain ⟨e, he⟩ := e
    induction e using Sym2.ind with
    | _ u w =>
      have hadj : G.Adj u w := G.mem_edgeSet.mp he
      show cbd G α ⟨s(u, w), _⟩ = (0 : G.edgeSet → F2) _
      rw [cbd_apply_pair G hadj, hconst u w]
      exact f2_add_self _

/-! ### `ker δ` is one-dimensional -/

/-- Evaluation at a base point, as a linear map out of `ker δ`. -/
noncomputable def evalAt (v₀ : V) : LinearMap.ker (cbd G) →ₗ[F2] F2 where
  toFun α := α.val v₀
  map_add' _ _ := rfl
  map_smul' _ _ := rfl

theorem finrank_ker_cbd (hG : G.Connected) (v₀ : V) :
    finrank F2 (LinearMap.ker (cbd G)) = 1 := by
  classical
  have hinj : Function.Injective (evalAt G v₀) := by
    rw [injective_iff_map_eq_zero]
    intro α hα
    ext v
    have hconst := (mem_ker_cbd_iff G hG α.val).mp α.property
    show α.val v = 0
    rw [hconst v v₀]
    exact hα
  have hsurj : Function.Surjective (evalAt G v₀) := by
    intro c
    refine ⟨⟨fun _ => c, ?_⟩, rfl⟩
    exact (mem_ker_cbd_iff G hG _).mpr fun _ _ => rfl
  have hiso : LinearMap.ker (cbd G) ≃ₗ[F2] F2 :=
    LinearEquiv.ofBijective (evalAt G v₀) ⟨hinj, hsurj⟩
  rw [hiso.finrank_eq]
  simp

/-! ### General version: `#components` instead of connectivity

  The connected case is what the knight application needs for `G_n` itself, but
  `Z_bulk` is the cycle space of a graph with FIVE components (the bulk plus the
  four isolated corners), so the general count is required downstream. -/

/-- Pull a function on components back to a function on vertices. -/
noncomputable def liftComp :
    (G.ConnectedComponent → F2) →ₗ[F2] (V → F2) where
  toFun f v := f (G.connectedComponentMk v)
  map_add' _ _ := rfl
  map_smul' _ _ := rfl

theorem liftComp_mem_ker (f : G.ConnectedComponent → F2) :
    liftComp G f ∈ LinearMap.ker (cbd G) := by
  rw [LinearMap.mem_ker]
  funext e
  obtain ⟨e, he⟩ := e
  induction e using Sym2.ind with
  | _ u w =>
    have hadj : G.Adj u w := G.mem_edgeSet.mp he
    show cbd G (liftComp G f) ⟨s(u, w), _⟩ = (0 : G.edgeSet → F2) _
    rw [cbd_apply_pair G hadj]
    have : G.connectedComponentMk u = G.connectedComponentMk w :=
      SimpleGraph.ConnectedComponent.sound hadj.reachable
    show f (G.connectedComponentMk u) + f (G.connectedComponentMk w) = 0
    rw [this]
    exact f2_add_self _

/-- `ker δ ≃ 𝔽₂^{components}`: an edgewise-constant function is exactly a
    function on the set of connected components. -/
noncomputable def kerCbdEquiv :
    (G.ConnectedComponent → F2) ≃ₗ[F2] LinearMap.ker (cbd G) := by
  refine LinearEquiv.ofBijective
    ((liftComp G).codRestrict (LinearMap.ker (cbd G)) (liftComp_mem_ker G)) ⟨?_, ?_⟩
  · rw [injective_iff_map_eq_zero]
    intro f hf
    funext c
    obtain ⟨v, rfl⟩ := c.exists_rep
    have : liftComp G f v = 0 := congrFun (congrArg Subtype.val hf) v
    exact this
  · rintro ⟨α, hα⟩
    have hz : cbd G α = 0 := LinearMap.mem_ker.mp hα
    have hedge : ∀ u w, G.Adj u w → α u = α w := by
      intro u w hadj
      have h0 : cbd G α ⟨s(u, w), G.mem_edgeSet.mpr hadj⟩ = 0 := by rw [hz]; rfl
      rw [cbd_apply_pair G hadj] at h0
      exact f2_of_add_eq_zero h0
    refine ⟨fun c => Quot.liftOn c α (fun a b h => edgewise_const_of_reachable α hedge h), ?_⟩
    apply Subtype.ext
    funext v
    rfl

theorem finrank_ker_cbd_general [Fintype G.ConnectedComponent] :
    finrank F2 (LinearMap.ker (cbd G)) = Fintype.card G.ConnectedComponent := by
  rw [← (kerCbdEquiv G).finrank_eq, finrank_pi]

/-- **Cycle space dimension, general form.** For any finite simple graph,
    `dim_{𝔽₂} ker ∂ = |E| − |V| + #components`. -/
theorem finrank_cycleSpace_general [Fintype G.ConnectedComponent] :
    (finrank F2 (cycleSpace G) : Int)
      = (Fintype.card G.edgeSet : Int) - (Fintype.card V : Int)
        + (Fintype.card G.ConnectedComponent : Int) := by
  classical
  have hd : finrank F2 (LinearMap.range (cbd G))
      + finrank F2 (LinearMap.ker (cbd G)) = Fintype.card V := by
    have h := (cbd G).finrank_range_add_finrank_ker
    simpa [finrank_pi] using h
  have hbdy : finrank F2 (LinearMap.range (bd G))
      + finrank F2 (LinearMap.ker (bd G)) = Fintype.card G.edgeSet := by
    have h := (bd G).finrank_range_add_finrank_ker
    simpa [finrank_pi] using h
  have hrk : finrank F2 (LinearMap.range (bd G))
      = finrank F2 (LinearMap.range (cbd G)) := by
    have h : (inc G)ᵀ.rank = (inc G).rank := Matrix.rank_transpose _
    simpa [Matrix.rank, bd, cbd] using h.symm
  have hone := finrank_ker_cbd_general G
  have : finrank F2 (cycleSpace G) = finrank F2 (LinearMap.ker (bd G)) := rfl
  omega

/-! ### The main theorem -/

/-- **Cycle space dimension.** For a finite connected simple graph,
    `dim_{𝔽₂} ker ∂ = |E| − |V| + 1`. -/
theorem finrank_cycleSpace (hG : G.Connected) (v₀ : V) :
    (finrank F2 (cycleSpace G) : Int)
      = (Fintype.card G.edgeSet : Int) - (Fintype.card V : Int) + 1 := by
  classical
  have hδ : finrank F2 (LinearMap.range (cbd G))
      + finrank F2 (LinearMap.ker (cbd G)) = Fintype.card V := by
    have h := (cbd G).finrank_range_add_finrank_ker
    simpa [finrank_pi] using h
  have hbdy : finrank F2 (LinearMap.range (bd G))
      + finrank F2 (LinearMap.ker (bd G)) = Fintype.card G.edgeSet := by
    have h := (bd G).finrank_range_add_finrank_ker
    simpa [finrank_pi] using h
  have hrk : finrank F2 (LinearMap.range (bd G))
      = finrank F2 (LinearMap.range (cbd G)) := by
    have h : (inc G)ᵀ.rank = (inc G).rank := Matrix.rank_transpose _
    simpa [Matrix.rank, bd, cbd] using h.symm
  have hone := finrank_ker_cbd G hG v₀
  have : finrank F2 (cycleSpace G) = finrank F2 (LinearMap.ker (bd G)) := rfl
  omega

/-! ### Transfer: a spanning subgraph's cycle space, seen inside `𝔽₂^{E(G)}`

  The downstream need (`Reduction.ZBulk`) is a subspace of `𝔽₂^{E(G)}`, not of
  `𝔽₂^{E(G')}`. This section builds the extension-by-zero isomorphism so that
  the dimension formula can be applied to `G'` and read off in `G`'s edge
  space. -/

section Transfer

variable {V : Type*} [Fintype V] [DecidableEq V]
variable (G G' : SimpleGraph V) [DecidableRel G.Adj] [DecidableRel G'.Adj]
variable [Fintype G.edgeSet] [Fintype G'.edgeSet]

/-- The inclusion of edge sets induced by `G' ≤ G`. -/
def edgeIncl (hsub : G' ≤ G) : G'.edgeSet → G.edgeSet :=
  fun e => ⟨e.val, SimpleGraph.edgeSet_mono hsub e.property⟩

theorem edgeIncl_injective (hsub : G' ≤ G) :
    Function.Injective (edgeIncl G G' hsub) := by
  intro a b h
  simpa [edgeIncl, Subtype.ext_iff] using h

/-- Extension by zero, `𝔽₂^{E(G')} → 𝔽₂^{E(G)}`. -/
noncomputable def extEdge (hsub : G' ≤ G) :
    (G'.edgeSet → F2) →ₗ[F2] (G.edgeSet → F2) where
  toFun f e := if h : e.val ∈ G'.edgeSet then f ⟨e.val, h⟩ else 0
  map_add' f g := by
    funext e; by_cases h : e.val ∈ G'.edgeSet <;> simp [h]
  map_smul' c f := by
    funext e; by_cases h : e.val ∈ G'.edgeSet <;> simp [h]

theorem extEdge_incl (hsub : G' ≤ G) (f : G'.edgeSet → F2) (e : G'.edgeSet) :
    extEdge G G' hsub f (edgeIncl G G' hsub e) = f e := by
  simp only [extEdge, edgeIncl, LinearMap.coe_mk, AddHom.coe_mk]
  rw [dif_pos e.property]

theorem extEdge_injective (hsub : G' ≤ G) :
    Function.Injective (extEdge G G' hsub) := by
  intro f g h
  funext e
  have := congrFun h (edgeIncl G G' hsub e)
  rwa [extEdge_incl, extEdge_incl] at this

/-- **Extension by zero commutes with the boundary.** Both sides count, for a
    fixed vertex `v`, the `G'`-edges at `v` carrying a 1. -/
theorem bd_extEdge (hsub : G' ≤ G) (f : G'.edgeSet → F2) (v : V) :
    bd G (extEdge G G' hsub f) v = bd G' f v := by
  classical
  simp only [bd, Matrix.mulVecLin_apply, Matrix.mulVec, dotProduct, inc]
  set emb : G'.edgeSet ↪ G.edgeSet :=
    ⟨edgeIncl G G' hsub, edgeIncl_injective G G' hsub⟩ with hemb
  set F : G.edgeSet → F2 :=
    fun e => (if v ∈ e.val then (1 : F2) else 0) * extEdge G G' hsub f e with hF
  -- terms of `G` outside the image of `G'` vanish, so the sum collapses
  have h1 : ∑ e : G.edgeSet, F e = ∑ e ∈ Finset.univ.map emb, F e := by
    refine (Finset.sum_subset (Finset.subset_univ _) ?_).symm
    intro e _ he
    have hnot : e.val ∉ G'.edgeSet := by
      intro hmem
      exact he (Finset.mem_map.mpr ⟨⟨e.val, hmem⟩, Finset.mem_univ _, rfl⟩)
    simp [hF, extEdge, hnot]
  rw [h1, Finset.sum_map]
  refine Finset.sum_congr rfl fun e _ => ?_
  simp only [hF, hemb, Function.Embedding.coeFn_mk,
             extEdge_incl G G' hsub f e]
  rfl

/-- The image of the cycle space of `G'` under extension by zero: cycles of `G`
    supported on `E(G')`. This is exactly the shape of `ZBulk`. -/
noncomputable def cycleSpaceOn (hsub : G' ≤ G) : Submodule F2 (G.edgeSet → F2) :=
  (cycleSpace G').map (extEdge G G' hsub)

theorem cycleSpaceOn_le_cycleSpace (hsub : G' ≤ G) :
    cycleSpaceOn G G' hsub ≤ cycleSpace G := by
  rintro x ⟨f, hf, rfl⟩
  have hf' : bd G' f = 0 := hf
  show bd G ((extEdge G G' hsub) f) = 0
  funext v
  rw [bd_extEdge G G' hsub f v]
  exact congrFun hf' v

/-- **Dimension transfers.** -/
theorem finrank_cycleSpaceOn (hsub : G' ≤ G) :
    finrank F2 (cycleSpaceOn G G' hsub) = finrank F2 (cycleSpace G') :=
  (Submodule.equivMapOfInjective _ (extEdge_injective G G' hsub) _).symm.finrank_eq

/-- **The formula, in `G`'s edge space.** -/
theorem finrank_cycleSpaceOn_eq (hsub : G' ≤ G) [Fintype G'.ConnectedComponent] :
    (finrank F2 (cycleSpaceOn G G' hsub) : Int)
      = (Fintype.card G'.edgeSet : Int) - (Fintype.card V : Int)
        + (Fintype.card G'.ConnectedComponent : Int) := by
  rw [finrank_cycleSpaceOn]
  exact finrank_cycleSpace_general G'

end Transfer

end KnightTour
