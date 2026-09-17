/-
  Bulk connectivity of the rectangular knight graph — GENERAL min(n,m) ≥ 6.

  Target (paper Prop 2.19, the open ingredient of Q(n,m)=3 general):
    For all n,m ≥ 6, BulkR(n,m) := G_{n×m}[V \ Corners] is connected.

  Proof: a 2-parameter induction. The vertical embedding τ(i,j)=(i+1,j) maps the
  n×m board into the interior rows {1..n} of the (n+2)×m board; knight adjacency
  is row-translation invariant, so τ is a graph hom and the embedded block is
  connected by IH. Every frame vertex (rows 0, n+1) and embedded corner has an
  explicit knight neighbour in the block. The horizontal direction is handled by
  the coordinate-swap isomorphism BulkR(n,m) ≅ BulkR(m,n). Base cases {6,7}².
-/
import KnightTour.Basic
import Mathlib.Combinatorics.SimpleGraph.Path
import Mathlib.Combinatorics.SimpleGraph.Connectivity.WalkCounting
import Mathlib.Combinatorics.SimpleGraph.Maps
import Mathlib.Tactic

namespace KnightTour

/-- A cell on the rectangular n×m board. -/
abbrev RSquare (n m : ℕ) := Fin n × Fin m

/-- Corner of the n×m board: both coordinates extremal. -/
def isCornerR (n m : ℕ) (v : RSquare n m) : Prop :=
  (v.1.val = 0 ∨ v.1.val = n - 1) ∧ (v.2.val = 0 ∨ v.2.val = m - 1)

instance (n m) : DecidablePred (isCornerR n m) := by
  intro v; unfold isCornerR; infer_instance

/-- Knight adjacency on the rectangle. -/
def KAdjR (n m : ℕ) (u v : RSquare n m) : Prop :=
  (((u.1.val : Int) - v.1.val).natAbs) * (((u.2.val : Int) - v.2.val).natAbs) = 2

instance (n m u v) : Decidable (KAdjR n m u v) := by unfold KAdjR; infer_instance

/-- The rectangular bulk graph: knight graph induced on the non-corner cells. -/
def BulkR (n m : ℕ) : SimpleGraph {v : RSquare n m // ¬ isCornerR n m v} where
  Adj a b := KAdjR n m a.val b.val
  symm := by
    intro a b h
    unfold KAdjR at h ⊢
    have h1 : (((b.val).1.val : Int) - (a.val).1.val).natAbs
            = (((a.val).1.val : Int) - (b.val).1.val).natAbs := by omega
    have h2 : (((b.val).2.val : Int) - (a.val).2.val).natAbs
            = (((a.val).2.val : Int) - (b.val).2.val).natAbs := by omega
    rw [h1, h2]; exact h
  loopless := by intro a h; unfold KAdjR at h; simp at h

instance (n m) : DecidableRel (BulkR n m).Adj := by
  intro a b; unfold BulkR; simp only; infer_instance

/-! ## Base cases (both parities per dimension), kernel-verified. -/

theorem bulkR_conn_66 : (BulkR 6 6).Connected := by native_decide
theorem bulkR_conn_67 : (BulkR 6 7).Connected := by native_decide
theorem bulkR_conn_76 : (BulkR 7 6).Connected := by native_decide
theorem bulkR_conn_77 : (BulkR 7 7).Connected := by native_decide

/-! ## Coordinate-swap isomorphism BulkR(n,m) ≅ BulkR(m,n). -/

/-- The swap map on bulk vertices. -/
def swapR (n m : ℕ) (a : {v : RSquare n m // ¬ isCornerR n m v}) :
    {w : RSquare m n // ¬ isCornerR m n w} :=
  ⟨(a.val.2, a.val.1), by
    intro hc; exact a.property (by simp only [isCornerR] at hc ⊢; tauto)⟩

theorem swapR_adj (n m : ℕ) (a b : {v : RSquare n m // ¬ isCornerR n m v})
    (h : (BulkR n m).Adj a b) : (BulkR m n).Adj (swapR n m a) (swapR n m b) := by
  have h' : KAdjR n m a.val b.val := h
  show KAdjR m n (swapR n m a).val (swapR n m b).val
  unfold KAdjR at h' ⊢
  simp only [swapR]
  rw [mul_comm]; exact h'

def swapHomR (n m : ℕ) : BulkR n m →g BulkR m n where
  toFun := swapR n m
  map_rel' := by intro a b h; exact swapR_adj n m a b h

theorem swapR_swapR (n m : ℕ) (a : {v : RSquare n m // ¬ isCornerR n m v}) :
    swapR m n (swapR n m a) = a := by
  apply Subtype.ext; simp only [swapR]

/-- Connectedness is swap-invariant. -/
theorem bulkR_swap_connected (n m : ℕ) (hC : (BulkR m n).Connected) :
    (BulkR n m).Connected := by
  rw [SimpleGraph.connected_iff_exists_forall_reachable]
  obtain ⟨c, hc⟩ := (SimpleGraph.connected_iff_exists_forall_reachable _).mp hC
  refine ⟨swapR m n c, ?_⟩
  intro w
  have h : (BulkR n m).Reachable (swapR m n c) (swapR m n (swapR n m w)) :=
    (hc (swapR n m w)).map (swapHomR m n)
  rwa [swapR_swapR] at h

/-! ## Vertical embedding τ(i,j) = (i+1,j) : (n×m) ↪ ((n+2)×m). -/

def embR (n m : ℕ) (a : {v : RSquare n m // ¬ isCornerR n m v}) :
    {w : RSquare (n+2) m // ¬ isCornerR (n+2) m w} :=
  ⟨(⟨a.val.1.val + 1, by have := a.val.1.isLt; omega⟩, a.val.2),
   by
     have h1 := a.val.1.isLt
     intro hc
     simp only [isCornerR] at hc
     omega⟩

@[simp] theorem embR_fst (n m a) : (embR n m a).val.1.val = a.val.1.val + 1 := rfl
@[simp] theorem embR_snd (n m a) : (embR n m a).val.2.val = a.val.2.val := rfl

theorem embR_adj (n m : ℕ) (a b : {v : RSquare n m // ¬ isCornerR n m v})
    (h : (BulkR n m).Adj a b) : (BulkR (n+2) m).Adj (embR n m a) (embR n m b) := by
  have h' : KAdjR n m a.val b.val := h
  show KAdjR (n+2) m (embR n m a).val (embR n m b).val
  unfold KAdjR at h' ⊢
  have r1 : (((embR n m a).val.1.val : Int) - (embR n m b).val.1.val).natAbs
          = ((a.val.1.val : Int) - b.val.1.val).natAbs := by
    rw [embR_fst, embR_fst]; omega
  have r2 : (((embR n m a).val.2.val : Int) - (embR n m b).val.2.val).natAbs
          = ((a.val.2.val : Int) - b.val.2.val).natAbs := by
    rw [embR_snd, embR_snd]
  rw [r1, r2]; exact h'

def embHomR (n m : ℕ) : BulkR n m →g BulkR (n+2) m where
  toFun := embR n m
  map_rel' := by intro a b h; exact embR_adj n m a b h

theorem reach_embR (n m : ℕ) (hC : (BulkR n m).Connected) (a b) :
    (BulkR (n+2) m).Reachable (embR n m a) (embR n m b) :=
  (hC.preconnected a b).map (embHomR n m)

/-! ## Witnesses for the vertical step. -/

private theorem rnmul21 (a b : Int) (ha : a.natAbs = 2) (hb : b.natAbs = 1) :
    a.natAbs * b.natAbs = 2 := by rw [ha, hb]
private theorem rnmul12 (a b : Int) (ha : a.natAbs = 1) (hb : b.natAbs = 2) :
    a.natAbs * b.natAbs = 2 := by rw [ha, hb]

/-- Given an interior target `(ti,tj)` of the embedded block (row `1≤ti≤n`,
    column `tj<m`, not an embedded corner) knight-adjacent to `w`, exhibit a
    neighbour of `w` in the image of τ. -/
theorem witness_of_targetR (n m : ℕ)
    (w : {v : RSquare (n+2) m // ¬ isCornerR (n+2) m v}) (wi0 wj0 : ℕ)
    (hwi : w.val.1.val = wi0) (hwj : w.val.2.val = wj0) (ti tj : ℕ)
    (hti1 : 1 ≤ ti) (htin : ti ≤ n) (htjm : tj < m)
    (hnc : ¬ ((ti = 1 ∨ ti = n) ∧ (tj = 0 ∨ tj = m - 1)))
    (hk : ((wi0 : Int) - ti).natAbs * ((wj0 : Int) - tj).natAbs = 2) :
    ∃ u, (BulkR (n+2) m).Adj w (embR n m u) := by
  have hu : ¬ isCornerR n m ((⟨ti - 1, by omega⟩ : Fin n), (⟨tj, by omega⟩ : Fin m)) := by
    simp only [isCornerR]; omega
  refine ⟨⟨(⟨ti - 1, by omega⟩, ⟨tj, by omega⟩), hu⟩, ?_⟩
  show KAdjR (n+2) m w.val (embR n m ⟨(⟨ti - 1, by omega⟩, ⟨tj, by omega⟩), hu⟩).val
  unfold KAdjR
  rw [show (embR n m ⟨(⟨ti - 1, by omega⟩, ⟨tj, by omega⟩), hu⟩).val.1.val = ti from by
        simp only [embR_fst]; omega,
      show (embR n m ⟨(⟨ti - 1, by omega⟩, ⟨tj, by omega⟩), hu⟩).val.2.val = tj from by
        simp only [embR_snd],
      hwi, hwj]
  exact hk

/-- Every vertex of `BulkR (n+2) m` is in the image of τ or has a neighbour there. -/
theorem mem_image_or_witnessR (n m : ℕ) (hn : 6 ≤ n) (hm : 6 ≤ m)
    (w : {v : RSquare (n+2) m // ¬ isCornerR (n+2) m v}) :
    (∃ u, embR n m u = w) ∨ (∃ u, (BulkR (n+2) m).Adj w (embR n m u)) := by
  obtain ⟨⟨wi, wj⟩, hw⟩ := w
  have hbi : wi.val < n + 2 := wi.isLt
  have hbj : wj.val < m := wj.isLt
  have hncw : ¬ ((wi.val = 0 ∨ wi.val = n + 1) ∧ (wj.val = 0 ∨ wj.val = m - 1)) := by
    intro hcc; exact hw (by simp only [isCornerR]; omega)
  by_cases hImg : (1 ≤ wi.val ∧ wi.val ≤ n ∧
                   ¬ ((wi.val = 1 ∨ wi.val = n) ∧ (wj.val = 0 ∨ wj.val = m - 1)))
  · -- in image of τ
    left
    obtain ⟨hi1, hin, hnc⟩ := hImg
    have hu : ¬ isCornerR n m ((⟨wi.val - 1, by omega⟩ : Fin n), (⟨wj.val, by omega⟩ : Fin m)) := by
      simp only [isCornerR]; omega
    refine ⟨⟨(⟨wi.val - 1, by omega⟩, ⟨wj.val, by omega⟩), hu⟩, ?_⟩
    apply Subtype.ext
    rw [Prod.ext_iff]
    refine ⟨Fin.ext ?_, Fin.ext ?_⟩
    · rw [embR_fst]; change wi.val - 1 + 1 = wi.val; omega
    · rw [embR_snd]
  · right
    by_cases hi0 : wi.val = 0
    · -- top row: wj ∈ [1, m-2]
      exact witness_of_targetR n m ⟨(wi, wj), hw⟩ wi.val wj.val rfl rfl 2 (wj.val + 1)
        (by omega) (by omega) (by omega) (by omega)
        (rnmul21 _ _ (by omega) (by omega))
    · by_cases hin1 : wi.val = n + 1
      · -- bottom row
        exact witness_of_targetR n m ⟨(wi, wj), hw⟩ wi.val wj.val rfl rfl (n - 1) (wj.val + 1)
          (by omega) (by omega) (by omega) (by omega)
          (rnmul21 _ _ (by omega) (by omega))
      · -- embedded corner: 1 ≤ wi ≤ n, (wi=1∨wi=n) ∧ (wj=0∨wj=m-1)
        have hEC : (wi.val = 1 ∨ wi.val = n) ∧ (wj.val = 0 ∨ wj.val = m - 1) := by
          by_contra hcon
          exact hImg ⟨by omega, by omega, hcon⟩
        obtain ⟨hi', hj'⟩ := hEC
        rcases hi' with h1 | h1 <;> rcases hj' with h2 | h2
        · -- (1,0) → (3,1)
          exact witness_of_targetR n m ⟨(wi, wj), hw⟩ wi.val wj.val rfl rfl 3 1
            (by omega) (by omega) (by omega) (by omega)
            (rnmul21 _ _ (by omega) (by omega))
        · -- (1,m-1) → (3,m-2)
          exact witness_of_targetR n m ⟨(wi, wj), hw⟩ wi.val wj.val rfl rfl 3 (m - 2)
            (by omega) (by omega) (by omega) (by omega)
            (rnmul21 _ _ (by omega) (by omega))
        · -- (n,0) → (n-2,1)
          exact witness_of_targetR n m ⟨(wi, wj), hw⟩ wi.val wj.val rfl rfl (n - 2) 1
            (by omega) (by omega) (by omega) (by omega)
            (rnmul21 _ _ (by omega) (by omega))
        · -- (n,m-1) → (n-2,m-2)
          exact witness_of_targetR n m ⟨(wi, wj), hw⟩ wi.val wj.val rfl rfl (n - 2) (m - 2)
            (by omega) (by omega) (by omega) (by omega)
            (rnmul21 _ _ (by omega) (by omega))

/-- Vertical inductive step: `BulkR n m` connected ⇒ `BulkR (n+2) m` connected. -/
theorem bulkR_step (n m : ℕ) (hn : 6 ≤ n) (hm : 6 ≤ m) (hC : (BulkR n m).Connected) :
    (BulkR (n+2) m).Connected := by
  rw [SimpleGraph.connected_iff_exists_forall_reachable]
  refine ⟨embR n m ⟨(⟨1, by omega⟩, ⟨1, by omega⟩), by simp only [isCornerR]; omega⟩, ?_⟩
  intro w
  rcases mem_image_or_witnessR n m hn hm w with ⟨u, hu⟩ | ⟨u, hadj⟩
  · rw [← hu]; exact reach_embR n m hC _ u
  · exact (reach_embR n m hC _ u).trans hadj.symm.reachable

/-- **Main theorem.** The rectangular bulk is connected for all n, m ≥ 6.
    Closes the open ingredient of Prop 2.19 / Q(n,m)=3 general.
    Double induction on `n+m`: reduce the larger side by 2 (vertical step,
    via the swap iso when it is the second coordinate); base cases {6,7}². -/
theorem bulkR_connected_general :
    ∀ N n m, n + m = N → 6 ≤ n → 6 ≤ m → (BulkR n m).Connected := by
  intro N
  induction N using Nat.strong_induction_on with
  | _ N ih =>
    intro n m hN hn hm
    rcases Nat.lt_or_ge n 8 with hn8 | hn8
    · rcases Nat.lt_or_ge m 8 with hm8 | hm8
      · -- both sides in {6,7}: base cases
        interval_cases n <;> interval_cases m
        · exact bulkR_conn_66
        · exact bulkR_conn_67
        · exact bulkR_conn_76
        · exact bulkR_conn_77
      · -- m ≥ 8: reduce the second coordinate via the swap iso
        apply bulkR_swap_connected
        have hprev : (BulkR (m - 2) n).Connected :=
          ih ((m - 2) + n) (by omega) (m - 2) n rfl (by omega) (by omega)
        have hstep := bulkR_step (m - 2) n (by omega) (by omega) hprev
        rwa [Nat.sub_add_cancel (by omega)] at hstep
    · -- n ≥ 8: reduce the first coordinate (vertical step)
      have hprev : (BulkR (n - 2) m).Connected :=
        ih ((n - 2) + m) (by omega) (n - 2) m rfl (by omega) (by omega)
      have hstep := bulkR_step (n - 2) m (by omega) hm hprev
      rwa [Nat.sub_add_cancel (by omega)] at hstep

/-- Convenience form. -/
theorem bulkR_connected (n m : ℕ) (hn : 6 ≤ n) (hm : 6 ≤ m) :
    (BulkR n m).Connected :=
  bulkR_connected_general (n + m) n m rfl hn hm

end KnightTour
