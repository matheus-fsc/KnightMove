/-
  Q_abstract(n) = 3 for all n ≥ 6: the general theorem.

  This file defines Q via Submodule/finrank (the paper's Def 2.4) and proves
  Q_abstract n = 3 using the abstract linear algebra argument from Thm 2.5,
  WITHOUT bridging to the imperative computeQ / gf2Rank.

  Structure (following the validated blueprint):
    Layer 1: Knight SimpleGraph, coboundary LinearMap δ, R = range δ, quotient
    Layer 2: Mandatory edges, v_ij vectors, Q_abstract definition
    Layer 3: Corner representatives r_c, Corollary 2.8
    Layer 4: Independence of r_c₁..r_c₄ via bulk connectivity (Lemma 2.10)
    Layer 5: XOR image, σ-free endgame, Q_abstract = 3 (Thm 2.5)
-/
import KnightTour.Basic
import KnightTour.BulkConnectivity
import KnightTour.EdgewiseConst
import Mathlib.LinearAlgebra.Quotient.Basic
import Mathlib.LinearAlgebra.Dimension.Constructions
import Mathlib.LinearAlgebra.Span.Basic
import Mathlib.LinearAlgebra.LinearIndependent
import Mathlib.LinearAlgebra.FreeModule.Finite.Basic
import Mathlib.Combinatorics.SimpleGraph.Finite
import Mathlib.Data.ZMod.Basic

namespace KnightTour

/-! ## Layer 1: Knight SimpleGraph, coboundary LinearMap, quotient -/

/-- The field F₂ = ZMod 2. -/
abbrev F := ZMod 2

/-- The knight graph as a SimpleGraph on Square n. -/
def KnightGraph (n : ℕ) : SimpleGraph (Square n) where
  Adj u v := KAdj n u v ∧ u ≠ v
  symm := by
    intro u v ⟨h, hne⟩
    refine ⟨?_, hne.symm⟩
    unfold KAdj at h ⊢
    have h1 : ((v.1.val : Int) - u.1.val).natAbs
            = ((u.1.val : Int) - v.1.val).natAbs := by omega
    have h2 : ((v.2.val : Int) - u.2.val).natAbs
            = ((u.2.val : Int) - v.2.val).natAbs := by omega
    rw [h1, h2]; exact h
  loopless := by intro u ⟨_, hne⟩; exact hne rfl

instance knightGraphDecRel (n : ℕ) : DecidableRel (KnightGraph n).Adj := by
  intro u v; unfold KnightGraph; simp only; infer_instance

/-- The edge type: elements of the edge set of the knight graph. -/
abbrev KEdge (n : ℕ) := (KnightGraph n).edgeSet

instance kedgeFintype (n : ℕ) : Fintype (KEdge n) := SimpleGraph.fintypeEdgeSet _
instance kedgeDecEq (n : ℕ) : DecidableEq (KEdge n) := Subtype.instDecidableEq

/-- Build an edge from an adjacency proof. -/
def mkEdge (n : ℕ) {u v : Square n} (h : (KnightGraph n).Adj u v) : KEdge n :=
  ⟨s(u, v), (KnightGraph n).mem_edgeSet.mpr h⟩

/-- The coboundary map δ : F₂^V → F₂^E, defined by δ(α)(e) = α(u) + α(v)
    for e = {u,v}. Its range is the row space R(n) = row(∂₁). -/
noncomputable def coboundary (n : ℕ) : (Square n → F) →ₗ[F] (KEdge n → F) where
  toFun α e := Sym2.lift ⟨fun u v => α u + α v, fun _ _ => by ring⟩ e.val
  map_add' α β := by
    ext ⟨e, he⟩
    induction e using Sym2.ind with
    | _ u v => simp [Pi.add_apply]; ring
  map_smul' c α := by
    ext ⟨e, he⟩
    induction e using Sym2.ind with
    | _ u v => simp [Pi.smul_apply, RingHom.id_apply]; ring

/-- R(n) = range(δ) = row space of ∂₁. -/
noncomputable def R (n : ℕ) : Submodule F (KEdge n → F) :=
  LinearMap.range (coboundary n)

/-- The quotient module F₂^E / R(n). -/
noncomputable abbrev Qt (n : ℕ) := (KEdge n → F) ⧸ R n

/-- The quotient projection π : F₂^E → Qt. -/
noncomputable def π (n : ℕ) : (KEdge n → F) →ₗ[F] Qt n :=
  (R n).mkQ

/-! ## Layer 2: Corners, mandatory edges, Q_abstract -/

/-- The 4 corners of the n×n board, indexed by Fin 4. -/
noncomputable def corner (n : ℕ) (hn : 4 ≤ n) : Fin 4 → Square n
  | ⟨0, _⟩ => (⟨0, by omega⟩, ⟨0, by omega⟩)
  | ⟨1, _⟩ => (⟨0, by omega⟩, ⟨n-1, by omega⟩)
  | ⟨2, _⟩ => (⟨n-1, by omega⟩, ⟨0, by omega⟩)
  | ⟨3, _⟩ => (⟨n-1, by omega⟩, ⟨n-1, by omega⟩)

/-- The two knight-neighbours of each corner. For n ≥ 6, each corner has
    exactly degree 2 in the knight graph with these specific neighbours. -/
noncomputable def cornerNbr (n : ℕ) (hn : 6 ≤ n) : Fin 4 → Square n × Square n
  | ⟨0, _⟩ => ((⟨1, by omega⟩, ⟨2, by omega⟩), (⟨2, by omega⟩, ⟨1, by omega⟩))
  | ⟨1, _⟩ => ((⟨1, by omega⟩, ⟨n-3, by omega⟩), (⟨2, by omega⟩, ⟨n-2, by omega⟩))
  | ⟨2, _⟩ => ((⟨n-3, by omega⟩, ⟨1, by omega⟩), (⟨n-2, by omega⟩, ⟨2, by omega⟩))
  | ⟨3, _⟩ => ((⟨n-3, by omega⟩, ⟨n-2, by omega⟩), (⟨n-2, by omega⟩, ⟨n-3, by omega⟩))

-- Reuse the helper from BulkConnectivity for natAbs multiplication
private theorem nmul21' (a b : Int) (ha : a.natAbs = 2) (hb : b.natAbs = 1) :
    a.natAbs * b.natAbs = 2 := by rw [ha, hb]
private theorem nmul12' (a b : Int) (ha : a.natAbs = 1) (hb : b.natAbs = 2) :
    a.natAbs * b.natAbs = 2 := by rw [ha, hb]

/-- The only factorizations of 2 over ℕ. -/
private theorem nat_mul_eq_two {x y : ℕ} (h : x * y = 2) :
    (x = 1 ∧ y = 2) ∨ (x = 2 ∧ y = 1) := by
  have hx : x ≤ 2 := Nat.le_of_dvd (by norm_num) ⟨y, h.symm⟩
  interval_cases x <;> omega

/-- Corner k is KAdj to its first neighbour. -/
theorem corner_kadj_fst (n : ℕ) (hn : 6 ≤ n) (k : Fin 4) :
    KAdj n (corner n (by omega) k) (cornerNbr n hn k).1 := by
  fin_cases k <;> simp only [corner, cornerNbr, KAdj] <;>
    (first | exact nmul12' _ _ (by omega) (by omega)
           | exact nmul21' _ _ (by omega) (by omega))

/-- Corner k is not equal to its first neighbour. -/
theorem corner_ne_fst (n : ℕ) (hn : 6 ≤ n) (k : Fin 4) :
    corner n (by omega) k ≠ (cornerNbr n hn k).1 := by
  fin_cases k <;> simp [corner, cornerNbr, Prod.ext_iff, Fin.ext_iff] <;> omega

/-- Corner k is adjacent to its first neighbour in KnightGraph. -/
theorem corner_adj_fst (n : ℕ) (hn : 6 ≤ n) (k : Fin 4) :
    (KnightGraph n).Adj (corner n (by omega) k) (cornerNbr n hn k).1 :=
  ⟨corner_kadj_fst n hn k, corner_ne_fst n hn k⟩

/-- Corner k is KAdj to its second neighbour. -/
theorem corner_kadj_snd (n : ℕ) (hn : 6 ≤ n) (k : Fin 4) :
    KAdj n (corner n (by omega) k) (cornerNbr n hn k).2 := by
  fin_cases k <;> simp only [corner, cornerNbr, KAdj] <;>
    (first | exact nmul12' _ _ (by omega) (by omega)
           | exact nmul21' _ _ (by omega) (by omega))

/-- Corner k is not equal to its second neighbour. -/
theorem corner_ne_snd (n : ℕ) (hn : 6 ≤ n) (k : Fin 4) :
    corner n (by omega) k ≠ (cornerNbr n hn k).2 := by
  fin_cases k <;> simp [corner, cornerNbr, Prod.ext_iff, Fin.ext_iff] <;> omega

/-- Corner k is adjacent to its second neighbour in KnightGraph. -/
theorem corner_adj_snd (n : ℕ) (hn : 6 ≤ n) (k : Fin 4) :
    (KnightGraph n).Adj (corner n (by omega) k) (cornerNbr n hn k).2 :=
  ⟨corner_kadj_snd n hn k, corner_ne_snd n hn k⟩

/-- **Lemma 2.6 (abstract).** For n ≥ 6, a corner's only knight-neighbours are
    its two designated ones. This is the "degree exactly 2" characterisation. -/
theorem corner_adj_iff (n : ℕ) (hn : 6 ≤ n) (k : Fin 4) (w : Square n) :
    (KnightGraph n).Adj (corner n (by omega) k) w ↔
      w = (cornerNbr n hn k).1 ∨ w = (cornerNbr n hn k).2 := by
  constructor
  · rintro ⟨hadj, -⟩
    obtain ⟨⟨a, ha⟩, ⟨b, hb⟩⟩ := w
    fin_cases k <;>
      simp only [corner, cornerNbr, KAdj] at hadj ⊢ <;>
      rcases nat_mul_eq_two hadj with ⟨h1, h2⟩ | ⟨h1, h2⟩ <;>
      (first
        | (left; simp only [Prod.mk.injEq, Fin.mk.injEq]; omega)
        | (right; simp only [Prod.mk.injEq, Fin.mk.injEq]; omega))
  · rintro (rfl | rfl)
    · exact corner_adj_fst n hn k
    · exact corner_adj_snd n hn k

/-- The two designated neighbours of a corner are distinct. -/
theorem cornerNbr_distinct (n : ℕ) (hn : 6 ≤ n) (k : Fin 4) :
    (cornerNbr n hn k).1 ≠ (cornerNbr n hn k).2 := by
  fin_cases k <;> simp only [cornerNbr, ne_eq, Prod.mk.injEq, Fin.mk.injEq, not_and] <;> omega

/-- The first mandatory edge of corner k. -/
noncomputable def mandEdge1 (n : ℕ) (hn : 6 ≤ n) (k : Fin 4) : KEdge n :=
  mkEdge n (corner_adj_fst n hn k)

/-- The second mandatory edge of corner k. -/
noncomputable def mandEdge2 (n : ℕ) (hn : 6 ≤ n) (k : Fin 4) : KEdge n :=
  mkEdge n (corner_adj_snd n hn k)

/-- Indicator function: 1 at edge e, 0 elsewhere (as F₂-valued function). -/
noncomputable def edgeInd (n : ℕ) (e : KEdge n) : KEdge n → F :=
  Pi.single e 1

/-- XOR vector v_{e₁,e₂} = indicator(e₁) + indicator(e₂). -/
noncomputable def vij (n : ℕ) (e₁ e₂ : KEdge n) : KEdge n → F :=
  edgeInd n e₁ + edgeInd n e₂

/-- The set of all XOR vectors for mandatory edge pairs. -/
noncomputable def mandXorSet (n : ℕ) (hn : 6 ≤ n) : Set (KEdge n → F) :=
  { v | ∃ (k₁ k₂ : Fin 4) (b₁ b₂ : Bool),
    v = vij n (if b₁ then mandEdge1 n hn k₁ else mandEdge2 n hn k₁)
              (if b₂ then mandEdge1 n hn k₂ else mandEdge2 n hn k₂) }

/-- The span of XOR vectors as a Submodule. -/
noncomputable def S (n : ℕ) (hn : 6 ≤ n) : Submodule F (KEdge n → F) :=
  Submodule.span F (mandXorSet n hn)

/-- Q_abstract(n) := finrank F₂ (π(Span{v_ij})) — the paper's Definition 2.4. -/
noncomputable def Q_abstract (n : ℕ) (hn : 6 ≤ n) : ℕ :=
  Module.finrank F ((S n hn).map (π n))

/-! ## Layer 3: Corner representatives and Corollary 2.8 -/

/-- Lemma 2.7: e₁(c) + e₂(c) = δ(δ_c) ∈ R(n). -/
theorem mandEdge_xor_in_R (n : ℕ) (hn : 6 ≤ n) (k : Fin 4) :
    edgeInd n (mandEdge1 n hn k) + edgeInd n (mandEdge2 n hn k) ∈ R n := by
  rw [R, LinearMap.mem_range]
  refine ⟨Pi.single (corner n (by omega) k) 1, ?_⟩
  funext ε
  obtain ⟨e, he⟩ := ε
  revert he
  induction e using Sym2.ind with
  | _ u v =>
    intro he
    have hadj : (KnightGraph n).Adj u v := (KnightGraph n).mem_edgeSet.mp he
    have hc1 : corner n (by omega) k ≠ (cornerNbr n hn k).1 := corner_ne_fst n hn k
    have hc2 : corner n (by omega) k ≠ (cornerNbr n hn k).2 := corner_ne_snd n hn k
    have h12 : (cornerNbr n hn k).1 ≠ (cornerNbr n hn k).2 := cornerNbr_distinct n hn k
    simp only [coboundary, edgeInd, mandEdge1, mandEdge2, mkEdge, LinearMap.coe_mk,
      AddHom.coe_mk, Sym2.lift_mk, Pi.add_apply, Pi.single_apply, Subtype.mk.injEq,
      Sym2.eq_iff]
    by_cases hu : u = corner n (by omega) k <;> by_cases hv : v = corner n (by omega) k
    · exact absurd (hu.trans hv.symm) hadj.ne
    · subst hu
      rcases (corner_adj_iff n hn k v).mp hadj with hv1 | hv2
      · subst hv1
        simp [hc1, hc2, h12, Ne.symm hc1, Ne.symm hc2, Ne.symm h12, hv]
      · subst hv2
        simp [hc1, hc2, h12, Ne.symm hc1, Ne.symm hc2, Ne.symm h12, hv]
    · subst hv
      rcases (corner_adj_iff n hn k u).mp hadj.symm with hu1 | hu2
      · subst hu1
        simp [hc1, hc2, h12, Ne.symm hc1, Ne.symm hc2, Ne.symm h12, hu]
      · subst hu2
        simp [hc1, hc2, h12, Ne.symm hc1, Ne.symm hc2, Ne.symm h12, hu]
    · simp [hu, hv]

/-- Corollary 2.8: π(indicator(e₁(c))) = π(indicator(e₂(c))) in Qt. -/
theorem mandEdge_eq_in_quotient (n : ℕ) (hn : 6 ≤ n) (k : Fin 4) :
    π n (edgeInd n (mandEdge1 n hn k)) =
    π n (edgeInd n (mandEdge2 n hn k)) := by
  have h := mandEdge_xor_in_R n hn k
  have h0 : π n (edgeInd n (mandEdge1 n hn k) + edgeInd n (mandEdge2 n hn k)) = 0 := by
    rw [π]; exact (Submodule.Quotient.mk_eq_zero _).mpr h
  rw [map_add] at h0
  -- In a ZMod 2-module: a + b = 0 ↔ a = b (since -x = x)
  have neg_self : -(π n (edgeInd n (mandEdge2 n hn k))) = π n (edgeInd n (mandEdge2 n hn k)) :=
    ZModModule.neg_eq_self _
  rw [← neg_self]
  exact eq_neg_of_add_eq_zero_left h0

/-- The corner representative r(k) in the quotient. -/
noncomputable def r (n : ℕ) (hn : 6 ≤ n) (k : Fin 4) : Qt n :=
  π n (edgeInd n (mandEdge1 n hn k))

/-! ## Layer 4: Independence (Lemma 2.10) -/

/-- Corner neighbours are non-corner. -/
theorem cornerNbr_not_corner (n : ℕ) (hn : 6 ≤ n) (k : Fin 4) :
    ¬ isCorner n (cornerNbr n hn k).1 ∧ ¬ isCorner n (cornerNbr n hn k).2 := by
  fin_cases k <;> simp [cornerNbr, isCorner] <;> omega

/-- Each corner is in fact a corner. -/
theorem corner_isCorner (n : ℕ) (hn : 4 ≤ n) (k : Fin 4) :
    isCorner n (corner n hn k) := by
  fin_cases k
  · exact ⟨Or.inl rfl, Or.inl rfl⟩
  · exact ⟨Or.inl rfl, Or.inr rfl⟩
  · exact ⟨Or.inr rfl, Or.inl rfl⟩
  · exact ⟨Or.inr rfl, Or.inr rfl⟩

/-- The 8 mandatory edges are pairwise distinct (Step 5b). -/
theorem mandEdges_distinct (n : ℕ) (hn : 6 ≤ n) :
    Function.Injective (fun p : Fin 4 × Bool =>
      if p.2 then mandEdge1 n hn p.1 else mandEdge2 n hn p.1) := by
  rintro ⟨k1, b1⟩ ⟨k2, b2⟩ heq
  have hval := congrArg Subtype.val heq
  simp only [mandEdge1, mandEdge2, mkEdge, apply_ite Subtype.val] at hval
  fin_cases k1 <;> fin_cases k2 <;> cases b1 <;> cases b2 <;>
    simp_all only [corner, cornerNbr, Sym2.eq_iff, Prod.mk.injEq, Fin.mk.injEq,
      Bool.false_eq_true, reduceIte, ite_true, ite_false, Prod.ext_iff] <;>
    omega

/-- In F₂, `x + y = 0` forces `x = y`. -/
private theorem f2_add_eq_zero {x y : F} (h : x + y = 0) : x = y := by
  have hyy : y + y = 0 := by
    have h2 : (2 : F) * y = 0 := by rw [show (2 : F) = 0 from by decide, zero_mul]
    rwa [two_mul] at h2
  calc x = x + (y + y) := by rw [hyy, add_zero]
    _ = (x + y) + y := by rw [add_assoc]
    _ = y := by rw [h, zero_add]

/-- Lemma 2.10: r(0)..r(3) are linearly independent in Qt.
    Proof: a vanishing F₂-combination ∑ gₖ rₖ = 0 lifts to ∑ gₖ e₁(cₖ) = δ(α);
    α is constant (= a) off the corners by bulk connectivity + edgewise_const;
    then gⱼ = δα(e₁(cⱼ)) = α(cⱼ)+a = δα(e₂(cⱼ)) = 0 since e₂(cⱼ) is not a mandatory
    first-edge. -/
theorem corner_reps_linearIndependent (n : ℕ) (hn : 6 ≤ n) :
    LinearIndependent F (r n hn) := by
  rw [Fintype.linearIndependent_iff]
  intro g hg
  -- mandEdge1 is injective in k, and mandEdge1 k ≠ mandEdge2 j for all j,k.
  have mE1_inj : ∀ {a b : Fin 4}, mandEdge1 n hn a = mandEdge1 n hn b → a = b := by
    intro a b h
    have h2 : ((a, true) : Fin 4 × Bool) = (b, true) := mandEdges_distinct n hn (by simpa using h)
    exact (Prod.ext_iff.mp h2).1
  have mE12_ne : ∀ (a b : Fin 4), mandEdge1 n hn b ≠ mandEdge2 n hn a := by
    intro a b h
    have h2 : ((b, true) : Fin 4 × Bool) = (a, false) := mandEdges_distinct n hn (by simpa using h)
    simp at h2
  -- v* := ∑ gₖ e₁(cₖ); its value is gⱼ at e₁(cⱼ) and 0 at any non-first-edge.
  have hvstar_one : ∀ j, (∑ k, g k • edgeInd n (mandEdge1 n hn k)) (mandEdge1 n hn j) = g j := by
    intro j
    rw [Finset.sum_apply, Finset.sum_eq_single j]
    · rw [Pi.smul_apply, edgeInd, Pi.single_apply, if_pos rfl, smul_eq_mul, mul_one]
    · intro k _ hkj
      rw [Pi.smul_apply, edgeInd, Pi.single_apply, if_neg (fun hc => hkj (mE1_inj hc).symm),
        smul_zero]
    · intro hju; exact absurd (Finset.mem_univ j) hju
  have hvstar_zero : ∀ (e : KEdge n), (∀ k, e ≠ mandEdge1 n hn k) →
      (∑ k, g k • edgeInd n (mandEdge1 n hn k)) e = 0 := by
    intro e hne
    rw [Finset.sum_apply]
    apply Finset.sum_eq_zero
    intro k _
    rw [Pi.smul_apply, edgeInd, Pi.single_apply, if_neg (hne k), smul_zero]
  -- ∑ gₖ rₖ = 0 ⇒ v* ∈ R n ⇒ ∃ α, δ α = v*.
  have hv : (π n) (∑ k, g k • edgeInd n (mandEdge1 n hn k)) = 0 := by
    rw [map_sum]; simp only [map_smul]; exact hg
  rw [π, Submodule.mkQ_apply, Submodule.Quotient.mk_eq_zero, R, LinearMap.mem_range] at hv
  obtain ⟨α, hα⟩ := hv
  -- α is constant off the corners (bulk connectivity).
  have hbulkconst : ∀ (a b : {v : Square n // ¬ isCorner n v}), α a.val = α b.val := by
    refine edgewise_const_of_connected (bulk_connected_general n hn) (fun w => α w.val) ?_
    intro a b hab
    have hne : a.val ≠ b.val := fun h => hab.ne (Subtype.ext h)
    have hknight : (KnightGraph n).Adj a.val b.val := ⟨hab, hne⟩
    have hedge_ne : ∀ k, mkEdge n hknight ≠ mandEdge1 n hn k := by
      intro k hcontra
      rw [mkEdge, mandEdge1, mkEdge, Subtype.mk.injEq, Sym2.eq_iff] at hcontra
      rcases hcontra with ⟨h1, _⟩ | ⟨_, h2⟩
      · exact a.property (by rw [h1]; exact corner_isCorner n (by omega) k)
      · exact b.property (by rw [h2]; exact corner_isCorner n (by omega) k)
    have hcob : coboundary n α (mkEdge n hknight) = α a.val + α b.val := by
      simp [coboundary, mkEdge, Sym2.lift_mk]
    have hz : α a.val + α b.val = 0 := by
      rw [← hcob, hα]; exact hvstar_zero (mkEdge n hknight) hedge_ne
    exact f2_add_eq_zero hz
  -- Per corner j: gⱼ = δα(e₁) = α(cⱼ)+a = δα(e₂) = 0.
  intro j
  set cj := corner n (by omega) j with hcjdef
  have hnn : α (cornerNbr n hn j).1 = α (cornerNbr n hn j).2 :=
    hbulkconst ⟨_, (cornerNbr_not_corner n hn j).1⟩ ⟨_, (cornerNbr_not_corner n hn j).2⟩
  have cval1 : coboundary n α (mandEdge1 n hn j) = α cj + α (cornerNbr n hn j).1 := by
    simp [coboundary, mandEdge1, mkEdge, Sym2.lift_mk, hcjdef]
  have cval2 : coboundary n α (mandEdge2 n hn j) = α cj + α (cornerNbr n hn j).2 := by
    simp [coboundary, mandEdge2, mkEdge, Sym2.lift_mk, hcjdef]
  have e1 : coboundary n α (mandEdge1 n hn j) = g j := by rw [hα]; exact hvstar_one j
  have e2 : coboundary n α (mandEdge2 n hn j) = 0 := by
    rw [hα]; exact hvstar_zero _ (fun k => (mE12_ne j k).symm)
  rw [← e1, cval1, hnn, ← cval2, e2]

/-! ## Layer 5: XOR image + σ-free endgame -/

/-- Lemma 2.11: π(S) = Span{r(i) + r(j) : i < j}. -/
theorem xor_image_span (n : ℕ) (hn : 6 ≤ n) :
    (S n hn).map (π n) = Submodule.span F
      (Set.range (fun p : {p : Fin 4 × Fin 4 // p.1 < p.2} =>
        r n hn p.val.1 + r n hn p.val.2)) := by
  -- π sends the indicator of either mandatory edge of corner k to r(k).
  have hpi : ∀ (k : Fin 4) (b : Bool),
      π n (edgeInd n (if b then mandEdge1 n hn k else mandEdge2 n hn k)) = r n hn k := by
    intro k b; cases b
    · exact (mandEdge_eq_in_quotient n hn k).symm
    · rfl
  -- r(a)+r(b) lies in the span of the i<j pairs (handle a=b, a<b, a>b).
  have hrr : ∀ a b : Fin 4, r n hn a + r n hn b ∈
      Submodule.span F (Set.range (fun p : {p : Fin 4 × Fin 4 // p.1 < p.2} =>
        r n hn p.val.1 + r n hn p.val.2)) := by
    intro a b
    rcases lt_trichotomy a b with h | h | h
    · exact Submodule.subset_span ⟨⟨(a, b), h⟩, rfl⟩
    · subst h; rw [ZModModule.add_self]; exact Submodule.zero_mem _
    · rw [add_comm]; exact Submodule.subset_span ⟨⟨(b, a), h⟩, rfl⟩
  simp only [S, Submodule.map_span]
  apply le_antisymm
  · rw [Submodule.span_le]
    rintro y hy
    simp only [Set.mem_image, mandXorSet, Set.mem_setOf_eq] at hy
    obtain ⟨v, ⟨k1, k2, b1, b2, rfl⟩, rfl⟩ := hy
    simp only [vij, map_add, hpi]
    exact hrr k1 k2
  · rw [Submodule.span_le]
    rintro x hx
    simp only [Set.mem_range] at hx
    obtain ⟨⟨⟨i, j⟩, hij⟩, rfl⟩ := hx
    apply Submodule.subset_span
    rw [Set.mem_image]
    refine ⟨vij n (mandEdge1 n hn i) (mandEdge1 n hn j), ?_, ?_⟩
    · simp only [mandXorSet, Set.mem_setOf_eq]
      exact ⟨i, j, true, true, rfl⟩
    · simp only [vij, map_add]; rfl

/-- σ-free identity: Span{r(i)+r(j) : i<j} = Span{w₂,w₃,w₄} where w_k = r(0)+r(k). -/
theorem pairwise_span_eq_three (n : ℕ) (hn : 6 ≤ n) :
    Submodule.span F (Set.range (fun p : {p : Fin 4 × Fin 4 // p.1 < p.2} =>
        r n hn p.val.1 + r n hn p.val.2)) =
    Submodule.span F (Set.range (fun k : Fin 3 => r n hn 0 + r n hn (Fin.succ k))) := by
  -- r(0)+r(a) lies in the target span for every a (a=0 gives 0).
  have hw : ∀ a : Fin 4, r n hn 0 + r n hn a ∈
      Submodule.span F (Set.range (fun k : Fin 3 => r n hn 0 + r n hn (Fin.succ k))) := by
    intro a
    by_cases ha : a = 0
    · subst ha; rw [ZModModule.add_self]; exact Submodule.zero_mem _
    · have heq : r n hn 0 + r n hn a
          = (fun k : Fin 3 => r n hn 0 + r n hn (Fin.succ k)) (a.pred ha) := by
        simp only [Fin.succ_pred]
      rw [heq]; exact Submodule.subset_span ⟨a.pred ha, rfl⟩
  apply le_antisymm
  · rw [Submodule.span_le]
    rintro x ⟨⟨⟨i, j⟩, hij⟩, rfl⟩
    show r n hn i + r n hn j ∈ _
    rw [(ZModModule.add_add_add_cancel (r n hn i) (r n hn 0) (r n hn j)).symm]
    refine Submodule.add_mem _ ?_ (hw j)
    rw [add_comm]; exact hw i
  · rw [Submodule.span_le]
    rintro x ⟨k, rfl⟩
    exact Submodule.subset_span ⟨⟨(0, Fin.succ k), Fin.succ_pos k⟩, rfl⟩

/-- The three w-vectors are linearly independent (from independence of the 4 reps). -/
theorem three_vecs_linearIndependent (n : ℕ) (hn : 6 ≤ n) :
    LinearIndependent F (fun k : Fin 3 => r n hn 0 + r n hn (Fin.succ k)) := by
  rw [Fintype.linearIndependent_iff]
  intro c hc
  have hLI := corner_reps_linearIndependent n hn
  rw [Fintype.linearIndependent_iff] at hLI
  have key : (∑ k, c k) • r n hn 0 + ∑ i, c i • r n hn (Fin.succ i)
      = ∑ k, c k • (r n hn 0 + r n hn (Fin.succ k)) := by
    rw [Finset.sum_smul, ← Finset.sum_add_distrib]
    exact Finset.sum_congr rfl
      (fun k _ => (smul_add (c k) (r n hn 0) (r n hn (Fin.succ k))).symm)
  have hsum : ∑ j, (Fin.cons (∑ k, c k) c : Fin 4 → F) j • r n hn j = 0 := by
    rw [Fin.sum_univ_succ]
    simp only [Fin.cons_zero, Fin.cons_succ]
    rw [key]; exact hc
  have hz := hLI _ hsum
  intro i
  have hi := hz (Fin.succ i)
  rwa [Fin.cons_succ] at hi

/-- **Theorem 2.5: Q_abstract(n) = 3 for all n ≥ 6.** -/
theorem Q_abstract_eq_three (n : ℕ) (hn : 6 ≤ n) : Q_abstract n hn = 3 := by
  unfold Q_abstract
  rw [xor_image_span, pairwise_span_eq_three]
  have hli := three_vecs_linearIndependent n hn
  rw [finrank_span_eq_card hli, Fintype.card_fin]

end KnightTour
