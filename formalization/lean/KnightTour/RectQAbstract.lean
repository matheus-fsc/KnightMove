/-
  Q_abstract_rect(n,m) = 3 for all n,m ≥ 6: the rectangular general theorem
  (paper Prop 2.18 / 2.19).

  This is the rectangular twin of QAbstract.lean. The quotient machinery is
  square-specific there (not generic over the board), so per the blueprint's
  Cascade section this file DUPLICATES the abstract argument with rectangular
  coordinates. The ONE genuinely-new ingredient — rectangular bulk
  connectivity — is supplied by `bulkR_connected` (RectBulkConnectivity.lean).
  Everything else (Steps 1–4, 6–9) is a verbatim adaptation with the column
  size `n` replaced by `m`.
-/
import KnightTour.RectBulkConnectivity
import KnightTour.EdgewiseConst
import KnightTour.QAbstract
import Mathlib.LinearAlgebra.Quotient.Basic
import Mathlib.LinearAlgebra.Dimension.Constructions
import Mathlib.LinearAlgebra.Span.Basic
import Mathlib.LinearAlgebra.LinearIndependent
import Mathlib.Combinatorics.SimpleGraph.Finite
import Mathlib.Data.ZMod.Basic

namespace KnightTour

/-! ## Layer 1: rectangular knight SimpleGraph, coboundary, quotient -/

/-- The rectangular knight graph as a SimpleGraph on RSquare n m. -/
def KnightGraphR (n m : ℕ) : SimpleGraph (RSquare n m) where
  Adj u v := KAdjR n m u v ∧ u ≠ v
  symm := by
    intro u v ⟨h, hne⟩
    refine ⟨?_, hne.symm⟩
    unfold KAdjR at h ⊢
    have h1 : ((v.1.val : Int) - u.1.val).natAbs = ((u.1.val : Int) - v.1.val).natAbs := by omega
    have h2 : ((v.2.val : Int) - u.2.val).natAbs = ((u.2.val : Int) - v.2.val).natAbs := by omega
    rw [h1, h2]; exact h
  loopless := by intro u ⟨_, hne⟩; exact hne rfl

instance knightGraphRDecRel (n m : ℕ) : DecidableRel (KnightGraphR n m).Adj := by
  intro u v; unfold KnightGraphR; simp only; infer_instance

abbrev KEdgeR (n m : ℕ) := (KnightGraphR n m).edgeSet

instance kedgeRFintype (n m : ℕ) : Fintype (KEdgeR n m) := SimpleGraph.fintypeEdgeSet _
instance kedgeRDecEq (n m : ℕ) : DecidableEq (KEdgeR n m) := Subtype.instDecidableEq

def mkEdgeR (n m : ℕ) {u v : RSquare n m} (h : (KnightGraphR n m).Adj u v) : KEdgeR n m :=
  ⟨s(u, v), (KnightGraphR n m).mem_edgeSet.mpr h⟩

/-- The coboundary map δ : F₂^V → F₂^E for the rectangle. -/
noncomputable def coboundaryR (n m : ℕ) : (RSquare n m → F) →ₗ[F] (KEdgeR n m → F) where
  toFun α e := Sym2.lift ⟨fun u v => α u + α v, fun _ _ => by ring⟩ e.val
  map_add' α β := by
    ext ⟨e, he⟩
    induction e using Sym2.ind with
    | _ u v => simp [Pi.add_apply]; ring
  map_smul' c α := by
    ext ⟨e, he⟩
    induction e using Sym2.ind with
    | _ u v => simp [Pi.smul_apply, RingHom.id_apply]; ring

noncomputable def RR (n m : ℕ) : Submodule F (KEdgeR n m → F) :=
  LinearMap.range (coboundaryR n m)

noncomputable abbrev QtR (n m : ℕ) := (KEdgeR n m → F) ⧸ RR n m

noncomputable def πR (n m : ℕ) : (KEdgeR n m → F) →ₗ[F] QtR n m := (RR n m).mkQ

/-! ## Layer 2: corners, mandatory edges, Q_abstract_rect -/

noncomputable def cornerR (n m : ℕ) (hn : 4 ≤ n) (hm : 4 ≤ m) : Fin 4 → RSquare n m
  | ⟨0, _⟩ => (⟨0, by omega⟩, ⟨0, by omega⟩)
  | ⟨1, _⟩ => (⟨0, by omega⟩, ⟨m-1, by omega⟩)
  | ⟨2, _⟩ => (⟨n-1, by omega⟩, ⟨0, by omega⟩)
  | ⟨3, _⟩ => (⟨n-1, by omega⟩, ⟨m-1, by omega⟩)

noncomputable def cornerNbrR (n m : ℕ) (hn : 6 ≤ n) (hm : 6 ≤ m) :
    Fin 4 → RSquare n m × RSquare n m
  | ⟨0, _⟩ => ((⟨1, by omega⟩, ⟨2, by omega⟩), (⟨2, by omega⟩, ⟨1, by omega⟩))
  | ⟨1, _⟩ => ((⟨1, by omega⟩, ⟨m-3, by omega⟩), (⟨2, by omega⟩, ⟨m-2, by omega⟩))
  | ⟨2, _⟩ => ((⟨n-3, by omega⟩, ⟨1, by omega⟩), (⟨n-2, by omega⟩, ⟨2, by omega⟩))
  | ⟨3, _⟩ => ((⟨n-3, by omega⟩, ⟨m-2, by omega⟩), (⟨n-2, by omega⟩, ⟨m-3, by omega⟩))

private theorem rnmul21 (a b : Int) (ha : a.natAbs = 2) (hb : b.natAbs = 1) :
    a.natAbs * b.natAbs = 2 := by rw [ha, hb]
private theorem rnmul12 (a b : Int) (ha : a.natAbs = 1) (hb : b.natAbs = 2) :
    a.natAbs * b.natAbs = 2 := by rw [ha, hb]

private theorem rnat_mul_eq_two {x y : ℕ} (h : x * y = 2) :
    (x = 1 ∧ y = 2) ∨ (x = 2 ∧ y = 1) := by
  have hx : x ≤ 2 := Nat.le_of_dvd (by norm_num) ⟨y, h.symm⟩
  interval_cases x <;> omega

theorem corner_kadj_fstR (n m : ℕ) (hn : 6 ≤ n) (hm : 6 ≤ m) (k : Fin 4) :
    KAdjR n m (cornerR n m (by omega) (by omega) k) (cornerNbrR n m hn hm k).1 := by
  fin_cases k <;> simp only [cornerR, cornerNbrR, KAdjR] <;>
    (first | exact rnmul12 _ _ (by omega) (by omega)
           | exact rnmul21 _ _ (by omega) (by omega))

theorem corner_ne_fstR (n m : ℕ) (hn : 6 ≤ n) (hm : 6 ≤ m) (k : Fin 4) :
    cornerR n m (by omega) (by omega) k ≠ (cornerNbrR n m hn hm k).1 := by
  fin_cases k <;> simp [cornerR, cornerNbrR, Prod.ext_iff, Fin.ext_iff] <;> omega

theorem corner_adj_fstR (n m : ℕ) (hn : 6 ≤ n) (hm : 6 ≤ m) (k : Fin 4) :
    (KnightGraphR n m).Adj (cornerR n m (by omega) (by omega) k) (cornerNbrR n m hn hm k).1 :=
  ⟨corner_kadj_fstR n m hn hm k, corner_ne_fstR n m hn hm k⟩

theorem corner_kadj_sndR (n m : ℕ) (hn : 6 ≤ n) (hm : 6 ≤ m) (k : Fin 4) :
    KAdjR n m (cornerR n m (by omega) (by omega) k) (cornerNbrR n m hn hm k).2 := by
  fin_cases k <;> simp only [cornerR, cornerNbrR, KAdjR] <;>
    (first | exact rnmul12 _ _ (by omega) (by omega)
           | exact rnmul21 _ _ (by omega) (by omega))

theorem corner_ne_sndR (n m : ℕ) (hn : 6 ≤ n) (hm : 6 ≤ m) (k : Fin 4) :
    cornerR n m (by omega) (by omega) k ≠ (cornerNbrR n m hn hm k).2 := by
  fin_cases k <;> simp [cornerR, cornerNbrR, Prod.ext_iff, Fin.ext_iff] <;> omega

theorem corner_adj_sndR (n m : ℕ) (hn : 6 ≤ n) (hm : 6 ≤ m) (k : Fin 4) :
    (KnightGraphR n m).Adj (cornerR n m (by omega) (by omega) k) (cornerNbrR n m hn hm k).2 :=
  ⟨corner_kadj_sndR n m hn hm k, corner_ne_sndR n m hn hm k⟩

theorem corner_adj_iffR (n m : ℕ) (hn : 6 ≤ n) (hm : 6 ≤ m) (k : Fin 4) (w : RSquare n m) :
    (KnightGraphR n m).Adj (cornerR n m (by omega) (by omega) k) w ↔
      w = (cornerNbrR n m hn hm k).1 ∨ w = (cornerNbrR n m hn hm k).2 := by
  constructor
  · rintro ⟨hadj, -⟩
    obtain ⟨⟨a, ha⟩, ⟨b, hb⟩⟩ := w
    fin_cases k <;>
      simp only [cornerR, cornerNbrR, KAdjR] at hadj ⊢ <;>
      rcases rnat_mul_eq_two hadj with ⟨h1, h2⟩ | ⟨h1, h2⟩ <;>
      (first
        | (left; simp only [Prod.mk.injEq, Fin.mk.injEq]; omega)
        | (right; simp only [Prod.mk.injEq, Fin.mk.injEq]; omega))
  · rintro (rfl | rfl)
    · exact corner_adj_fstR n m hn hm k
    · exact corner_adj_sndR n m hn hm k

theorem cornerNbr_distinctR (n m : ℕ) (hn : 6 ≤ n) (hm : 6 ≤ m) (k : Fin 4) :
    (cornerNbrR n m hn hm k).1 ≠ (cornerNbrR n m hn hm k).2 := by
  fin_cases k <;> simp only [cornerNbrR, ne_eq, Prod.mk.injEq, Fin.mk.injEq, not_and] <;> omega

noncomputable def mandEdge1R (n m : ℕ) (hn : 6 ≤ n) (hm : 6 ≤ m) (k : Fin 4) : KEdgeR n m :=
  mkEdgeR n m (corner_adj_fstR n m hn hm k)

noncomputable def mandEdge2R (n m : ℕ) (hn : 6 ≤ n) (hm : 6 ≤ m) (k : Fin 4) : KEdgeR n m :=
  mkEdgeR n m (corner_adj_sndR n m hn hm k)

noncomputable def edgeIndR (n m : ℕ) (e : KEdgeR n m) : KEdgeR n m → F := Pi.single e 1

noncomputable def vijR (n m : ℕ) (e₁ e₂ : KEdgeR n m) : KEdgeR n m → F :=
  edgeIndR n m e₁ + edgeIndR n m e₂

noncomputable def mandXorSetR (n m : ℕ) (hn : 6 ≤ n) (hm : 6 ≤ m) : Set (KEdgeR n m → F) :=
  { v | ∃ (k₁ k₂ : Fin 4) (b₁ b₂ : Bool),
    v = vijR n m (if b₁ then mandEdge1R n m hn hm k₁ else mandEdge2R n m hn hm k₁)
                 (if b₂ then mandEdge1R n m hn hm k₂ else mandEdge2R n m hn hm k₂) }

noncomputable def SR (n m : ℕ) (hn : 6 ≤ n) (hm : 6 ≤ m) : Submodule F (KEdgeR n m → F) :=
  Submodule.span F (mandXorSetR n m hn hm)

/-- Q_abstract_rect(n,m) := finrank F₂ (π(Span{v_ij})), the paper's Def 2.4. -/
noncomputable def Q_abstract_rect (n m : ℕ) (hn : 6 ≤ n) (hm : 6 ≤ m) : ℕ :=
  Module.finrank F ((SR n m hn hm).map (πR n m))

/-! ## Layer 3: corner representatives and Cor 2.8 -/

theorem mandEdge_xor_in_RR (n m : ℕ) (hn : 6 ≤ n) (hm : 6 ≤ m) (k : Fin 4) :
    edgeIndR n m (mandEdge1R n m hn hm k) + edgeIndR n m (mandEdge2R n m hn hm k) ∈ RR n m := by
  rw [RR, LinearMap.mem_range]
  refine ⟨Pi.single (cornerR n m (by omega) (by omega) k) 1, ?_⟩
  funext ε
  obtain ⟨e, he⟩ := ε
  revert he
  induction e using Sym2.ind with
  | _ u v =>
    intro he
    have hadj : (KnightGraphR n m).Adj u v := (KnightGraphR n m).mem_edgeSet.mp he
    have hc1 : cornerR n m (by omega) (by omega) k ≠ (cornerNbrR n m hn hm k).1 :=
      corner_ne_fstR n m hn hm k
    have hc2 : cornerR n m (by omega) (by omega) k ≠ (cornerNbrR n m hn hm k).2 :=
      corner_ne_sndR n m hn hm k
    have h12 : (cornerNbrR n m hn hm k).1 ≠ (cornerNbrR n m hn hm k).2 :=
      cornerNbr_distinctR n m hn hm k
    simp only [coboundaryR, edgeIndR, mandEdge1R, mandEdge2R, mkEdgeR, LinearMap.coe_mk,
      AddHom.coe_mk, Sym2.lift_mk, Pi.add_apply, Pi.single_apply, Subtype.mk.injEq, Sym2.eq_iff]
    by_cases hu : u = cornerR n m (by omega) (by omega) k <;>
      by_cases hv : v = cornerR n m (by omega) (by omega) k
    · exact absurd (hu.trans hv.symm) hadj.ne
    · subst hu
      rcases (corner_adj_iffR n m hn hm k v).mp hadj with hv1 | hv2
      · subst hv1; simp [hc1, hc2, h12, Ne.symm hc1, Ne.symm hc2, Ne.symm h12, hv]
      · subst hv2; simp [hc1, hc2, h12, Ne.symm hc1, Ne.symm hc2, Ne.symm h12, hv]
    · subst hv
      rcases (corner_adj_iffR n m hn hm k u).mp hadj.symm with hu1 | hu2
      · subst hu1; simp [hc1, hc2, h12, Ne.symm hc1, Ne.symm hc2, Ne.symm h12, hu]
      · subst hu2; simp [hc1, hc2, h12, Ne.symm hc1, Ne.symm hc2, Ne.symm h12, hu]
    · simp [hu, hv]

theorem mandEdge_eq_in_quotientR (n m : ℕ) (hn : 6 ≤ n) (hm : 6 ≤ m) (k : Fin 4) :
    πR n m (edgeIndR n m (mandEdge1R n m hn hm k)) =
    πR n m (edgeIndR n m (mandEdge2R n m hn hm k)) := by
  have h := mandEdge_xor_in_RR n m hn hm k
  have h0 : πR n m (edgeIndR n m (mandEdge1R n m hn hm k)
      + edgeIndR n m (mandEdge2R n m hn hm k)) = 0 := by
    rw [πR]; exact (Submodule.Quotient.mk_eq_zero _).mpr h
  rw [map_add] at h0
  have neg_self : -(πR n m (edgeIndR n m (mandEdge2R n m hn hm k)))
      = πR n m (edgeIndR n m (mandEdge2R n m hn hm k)) := ZModModule.neg_eq_self _
  rw [← neg_self]
  exact eq_neg_of_add_eq_zero_left h0

noncomputable def rR (n m : ℕ) (hn : 6 ≤ n) (hm : 6 ≤ m) (k : Fin 4) : QtR n m :=
  πR n m (edgeIndR n m (mandEdge1R n m hn hm k))

/-! ## Layer 4: independence (Lemma 2.10) via rectangular bulk connectivity -/

theorem cornerNbr_not_cornerR (n m : ℕ) (hn : 6 ≤ n) (hm : 6 ≤ m) (k : Fin 4) :
    ¬ isCornerR n m (cornerNbrR n m hn hm k).1 ∧ ¬ isCornerR n m (cornerNbrR n m hn hm k).2 := by
  fin_cases k <;> simp [cornerNbrR, isCornerR] <;> omega

theorem corner_isCornerR (n m : ℕ) (hn : 4 ≤ n) (hm : 4 ≤ m) (k : Fin 4) :
    isCornerR n m (cornerR n m hn hm k) := by
  fin_cases k
  · exact ⟨Or.inl rfl, Or.inl rfl⟩
  · exact ⟨Or.inl rfl, Or.inr rfl⟩
  · exact ⟨Or.inr rfl, Or.inl rfl⟩
  · exact ⟨Or.inr rfl, Or.inr rfl⟩

theorem mandEdges_distinctR (n m : ℕ) (hn : 6 ≤ n) (hm : 6 ≤ m) :
    Function.Injective (fun p : Fin 4 × Bool =>
      if p.2 then mandEdge1R n m hn hm p.1 else mandEdge2R n m hn hm p.1) := by
  rintro ⟨k1, b1⟩ ⟨k2, b2⟩ heq
  have hval := congrArg Subtype.val heq
  simp only [mandEdge1R, mandEdge2R, mkEdgeR, apply_ite Subtype.val] at hval
  fin_cases k1 <;> fin_cases k2 <;> cases b1 <;> cases b2 <;>
    simp_all only [cornerR, cornerNbrR, Sym2.eq_iff, Prod.mk.injEq, Fin.mk.injEq,
      Bool.false_eq_true, reduceIte, ite_true, ite_false, Prod.ext_iff] <;>
    omega

private theorem rf2_add_eq_zero {x y : F} (h : x + y = 0) : x = y := by
  have hyy : y + y = 0 := by
    have h2 : (2 : F) * y = 0 := by rw [show (2 : F) = 0 from by decide, zero_mul]
    rwa [two_mul] at h2
  calc x = x + (y + y) := by rw [hyy, add_zero]
    _ = (x + y) + y := by rw [add_assoc]
    _ = y := by rw [h, zero_add]

theorem corner_reps_linearIndependentR (n m : ℕ) (hn : 6 ≤ n) (hm : 6 ≤ m) :
    LinearIndependent F (rR n m hn hm) := by
  rw [Fintype.linearIndependent_iff]
  intro g hg
  have mE1_inj : ∀ {a b : Fin 4}, mandEdge1R n m hn hm a = mandEdge1R n m hn hm b → a = b := by
    intro a b h
    have h2 : ((a, true) : Fin 4 × Bool) = (b, true) :=
      mandEdges_distinctR n m hn hm (by simpa using h)
    exact (Prod.ext_iff.mp h2).1
  have mE12_ne : ∀ (a b : Fin 4), mandEdge1R n m hn hm b ≠ mandEdge2R n m hn hm a := by
    intro a b h
    have h2 : ((b, true) : Fin 4 × Bool) = (a, false) :=
      mandEdges_distinctR n m hn hm (by simpa using h)
    simp at h2
  have hvstar_one : ∀ j, (∑ k, g k • edgeIndR n m (mandEdge1R n m hn hm k))
      (mandEdge1R n m hn hm j) = g j := by
    intro j
    rw [Finset.sum_apply, Finset.sum_eq_single j]
    · rw [Pi.smul_apply, edgeIndR, Pi.single_apply, if_pos rfl, smul_eq_mul, mul_one]
    · intro k _ hkj
      rw [Pi.smul_apply, edgeIndR, Pi.single_apply, if_neg (fun hc => hkj (mE1_inj hc).symm),
        smul_zero]
    · intro hju; exact absurd (Finset.mem_univ j) hju
  have hvstar_zero : ∀ (e : KEdgeR n m), (∀ k, e ≠ mandEdge1R n m hn hm k) →
      (∑ k, g k • edgeIndR n m (mandEdge1R n m hn hm k)) e = 0 := by
    intro e hne
    rw [Finset.sum_apply]
    apply Finset.sum_eq_zero
    intro k _
    rw [Pi.smul_apply, edgeIndR, Pi.single_apply, if_neg (hne k), smul_zero]
  have hv : (πR n m) (∑ k, g k • edgeIndR n m (mandEdge1R n m hn hm k)) = 0 := by
    rw [map_sum]; simp only [map_smul]; exact hg
  rw [πR, Submodule.mkQ_apply, Submodule.Quotient.mk_eq_zero, RR, LinearMap.mem_range] at hv
  obtain ⟨α, hα⟩ := hv
  have hbulkconst : ∀ (a b : {v : RSquare n m // ¬ isCornerR n m v}), α a.val = α b.val := by
    refine edgewise_const_of_connected (bulkR_connected n m hn hm) (fun w => α w.val) ?_
    intro a b hab
    have hne : a.val ≠ b.val := fun h => hab.ne (Subtype.ext h)
    have hknight : (KnightGraphR n m).Adj a.val b.val := ⟨hab, hne⟩
    have hedge_ne : ∀ k, mkEdgeR n m hknight ≠ mandEdge1R n m hn hm k := by
      intro k hcontra
      rw [mkEdgeR, mandEdge1R, mkEdgeR, Subtype.mk.injEq, Sym2.eq_iff] at hcontra
      rcases hcontra with ⟨h1, _⟩ | ⟨_, h2⟩
      · exact a.property (by rw [h1]; exact corner_isCornerR n m (by omega) (by omega) k)
      · exact b.property (by rw [h2]; exact corner_isCornerR n m (by omega) (by omega) k)
    have hcob : coboundaryR n m α (mkEdgeR n m hknight) = α a.val + α b.val := by
      simp [coboundaryR, mkEdgeR, Sym2.lift_mk]
    have hz : α a.val + α b.val = 0 := by
      rw [← hcob, hα]; exact hvstar_zero (mkEdgeR n m hknight) hedge_ne
    exact rf2_add_eq_zero hz
  intro j
  set cj := cornerR n m (by omega) (by omega) j with hcjdef
  have hnn : α (cornerNbrR n m hn hm j).1 = α (cornerNbrR n m hn hm j).2 :=
    hbulkconst ⟨_, (cornerNbr_not_cornerR n m hn hm j).1⟩ ⟨_, (cornerNbr_not_cornerR n m hn hm j).2⟩
  have cval1 : coboundaryR n m α (mandEdge1R n m hn hm j) = α cj + α (cornerNbrR n m hn hm j).1 := by
    simp [coboundaryR, mandEdge1R, mkEdgeR, Sym2.lift_mk, hcjdef]
  have cval2 : coboundaryR n m α (mandEdge2R n m hn hm j) = α cj + α (cornerNbrR n m hn hm j).2 := by
    simp [coboundaryR, mandEdge2R, mkEdgeR, Sym2.lift_mk, hcjdef]
  have e1 : coboundaryR n m α (mandEdge1R n m hn hm j) = g j := by rw [hα]; exact hvstar_one j
  have e2 : coboundaryR n m α (mandEdge2R n m hn hm j) = 0 := by
    rw [hα]; exact hvstar_zero _ (fun k => (mE12_ne j k).symm)
  rw [← e1, cval1, hnn, ← cval2, e2]

/-! ## Layer 5: XOR image + σ-free endgame -/

theorem xor_image_spanR (n m : ℕ) (hn : 6 ≤ n) (hm : 6 ≤ m) :
    (SR n m hn hm).map (πR n m) = Submodule.span F
      (Set.range (fun p : {p : Fin 4 × Fin 4 // p.1 < p.2} =>
        rR n m hn hm p.val.1 + rR n m hn hm p.val.2)) := by
  have hpi : ∀ (k : Fin 4) (b : Bool),
      πR n m (edgeIndR n m (if b then mandEdge1R n m hn hm k else mandEdge2R n m hn hm k))
        = rR n m hn hm k := by
    intro k b; cases b
    · exact (mandEdge_eq_in_quotientR n m hn hm k).symm
    · rfl
  have hrr : ∀ a b : Fin 4, rR n m hn hm a + rR n m hn hm b ∈
      Submodule.span F (Set.range (fun p : {p : Fin 4 × Fin 4 // p.1 < p.2} =>
        rR n m hn hm p.val.1 + rR n m hn hm p.val.2)) := by
    intro a b
    rcases lt_trichotomy a b with h | h | h
    · exact Submodule.subset_span ⟨⟨(a, b), h⟩, rfl⟩
    · subst h; rw [ZModModule.add_self]; exact Submodule.zero_mem _
    · rw [add_comm]; exact Submodule.subset_span ⟨⟨(b, a), h⟩, rfl⟩
  simp only [SR, Submodule.map_span]
  apply le_antisymm
  · rw [Submodule.span_le]
    rintro y hy
    simp only [Set.mem_image, mandXorSetR, Set.mem_setOf_eq] at hy
    obtain ⟨v, ⟨k1, k2, b1, b2, rfl⟩, rfl⟩ := hy
    simp only [vijR, map_add, hpi]
    exact hrr k1 k2
  · rw [Submodule.span_le]
    rintro x hx
    simp only [Set.mem_range] at hx
    obtain ⟨⟨⟨i, j⟩, hij⟩, rfl⟩ := hx
    apply Submodule.subset_span
    rw [Set.mem_image]
    refine ⟨vijR n m (mandEdge1R n m hn hm i) (mandEdge1R n m hn hm j), ?_, ?_⟩
    · simp only [mandXorSetR, Set.mem_setOf_eq]
      exact ⟨i, j, true, true, rfl⟩
    · simp only [vijR, map_add]; rfl

theorem pairwise_span_eq_threeR (n m : ℕ) (hn : 6 ≤ n) (hm : 6 ≤ m) :
    Submodule.span F (Set.range (fun p : {p : Fin 4 × Fin 4 // p.1 < p.2} =>
        rR n m hn hm p.val.1 + rR n m hn hm p.val.2)) =
    Submodule.span F (Set.range (fun k : Fin 3 => rR n m hn hm 0 + rR n m hn hm (Fin.succ k))) := by
  have hw : ∀ a : Fin 4, rR n m hn hm 0 + rR n m hn hm a ∈
      Submodule.span F (Set.range (fun k : Fin 3 =>
        rR n m hn hm 0 + rR n m hn hm (Fin.succ k))) := by
    intro a
    by_cases ha : a = 0
    · subst ha; rw [ZModModule.add_self]; exact Submodule.zero_mem _
    · have heq : rR n m hn hm 0 + rR n m hn hm a
          = (fun k : Fin 3 => rR n m hn hm 0 + rR n m hn hm (Fin.succ k)) (a.pred ha) := by
        simp only [Fin.succ_pred]
      rw [heq]; exact Submodule.subset_span ⟨a.pred ha, rfl⟩
  apply le_antisymm
  · rw [Submodule.span_le]
    rintro x ⟨⟨⟨i, j⟩, hij⟩, rfl⟩
    show rR n m hn hm i + rR n m hn hm j ∈ _
    rw [(ZModModule.add_add_add_cancel (rR n m hn hm i) (rR n m hn hm 0) (rR n m hn hm j)).symm]
    refine Submodule.add_mem _ ?_ (hw j)
    rw [add_comm]; exact hw i
  · rw [Submodule.span_le]
    rintro x ⟨k, rfl⟩
    exact Submodule.subset_span ⟨⟨(0, Fin.succ k), Fin.succ_pos k⟩, rfl⟩

theorem three_vecs_linearIndependentR (n m : ℕ) (hn : 6 ≤ n) (hm : 6 ≤ m) :
    LinearIndependent F (fun k : Fin 3 => rR n m hn hm 0 + rR n m hn hm (Fin.succ k)) := by
  rw [Fintype.linearIndependent_iff]
  intro c hc
  have hLI := corner_reps_linearIndependentR n m hn hm
  rw [Fintype.linearIndependent_iff] at hLI
  have key : (∑ k, c k) • rR n m hn hm 0 + ∑ i, c i • rR n m hn hm (Fin.succ i)
      = ∑ k, c k • (rR n m hn hm 0 + rR n m hn hm (Fin.succ k)) := by
    rw [Finset.sum_smul, ← Finset.sum_add_distrib]
    exact Finset.sum_congr rfl
      (fun k _ => (smul_add (c k) (rR n m hn hm 0) (rR n m hn hm (Fin.succ k))).symm)
  have hsum : ∑ j, (Fin.cons (∑ k, c k) c : Fin 4 → F) j • rR n m hn hm j = 0 := by
    rw [Fin.sum_univ_succ]
    simp only [Fin.cons_zero, Fin.cons_succ]
    rw [key]; exact hc
  have hz := hLI _ hsum
  intro i
  have hi := hz (Fin.succ i)
  rwa [Fin.cons_succ] at hi

/-- **Proposition 2.18 / 2.19: Q_abstract_rect(n,m) = 3 for all n,m ≥ 6.** -/
theorem Q_abstract_rect_eq_three (n m : ℕ) (hn : 6 ≤ n) (hm : 6 ≤ m) :
    Q_abstract_rect n m hn hm = 3 := by
  unfold Q_abstract_rect
  rw [xor_image_spanR, pairwise_span_eq_threeR]
  have hli := three_vecs_linearIndependentR n m hn hm
  rw [finrank_span_eq_card hli, Fintype.card_fin]

end KnightTour
