/-
  Bulk connectivity of the knight graph — GENERAL n ≥ 6.

  Target (paper Lemma 2.10(a), the open ingredient of Q(n)=3 general):
    For all n ≥ 6, Bulk(n) := G_n[V_n \ Corners(n)] is connected.

  Proof (validated on paper + adversarial review, see
  paper/notes/bulk_connectivity_proof_draft.tex):
    Induction n → n+2, base cases {6,7}.
    Step: embed the n-board as the interior block B = {1..n}² of the
    (n+2)-board via τ(i,j)=(i+1,j+1). Knight adjacency is translation
    invariant, so τ is a graph hom and B\EC (= image of τ) is connected
    by IH. Every remaining vertex (the 4 embedded corners EC and the 4n
    frame cells) has an explicit knight neighbour in B\EC, so the whole
    (n+2)-bulk is connected.

  This file builds the clean SimpleGraph formulation (independent of the
  BFS `bulkReachableCount` in Connectivity.lean) and proves the general
  theorem with no `sorry`.
-/
import KnightTour.Basic
import Mathlib.Combinatorics.SimpleGraph.Path
import Mathlib.Combinatorics.SimpleGraph.Connectivity.WalkCounting
import Mathlib.Combinatorics.SimpleGraph.Maps
import Mathlib.Tactic

namespace KnightTour

/-- A square is a corner of the n×n board iff both coordinates are extremal. -/
def isCorner (n : ℕ) (v : Square n) : Prop :=
  (v.1.val = 0 ∨ v.1.val = n - 1) ∧ (v.2.val = 0 ∨ v.2.val = n - 1)

instance (n) : DecidablePred (isCorner n) := by
  intro v; unfold isCorner; infer_instance

/-- Knight adjacency via the integer move signature |Δi|·|Δj| = 2. -/
def KAdj (n : ℕ) (u v : Square n) : Prop :=
  (((u.1.val : Int) - v.1.val).natAbs) * (((u.2.val : Int) - v.2.val).natAbs) = 2

instance (n u v) : Decidable (KAdj n u v) := by unfold KAdj; infer_instance

/-- The bulk graph: knight graph induced on the non-corner vertices. -/
def Bulk (n : ℕ) : SimpleGraph {v : Square n // ¬ isCorner n v} where
  Adj a b := KAdj n a.val b.val
  symm := by
    intro a b h
    unfold KAdj at h ⊢
    have h1 : (((b.val).1.val : Int) - (a.val).1.val).natAbs
            = (((a.val).1.val : Int) - (b.val).1.val).natAbs := by omega
    have h2 : (((b.val).2.val : Int) - (a.val).2.val).natAbs
            = (((a.val).2.val : Int) - (b.val).2.val).natAbs := by omega
    rw [h1, h2]; exact h
  loopless := by
    intro a h; unfold KAdj at h; simp at h

instance (n) : DecidableRel (Bulk n).Adj := by
  intro a b; unfold Bulk; simp only; infer_instance

/-! ## Base cases (n ∈ {6,7}, one per parity), kernel-verified. -/

theorem bulk_conn_6 : (Bulk 6).Connected := by native_decide
theorem bulk_conn_7 : (Bulk 7).Connected := by native_decide
-- extra safety margin
theorem bulk_conn_8 : (Bulk 8).Connected := by native_decide
theorem bulk_conn_10 : (Bulk 10).Connected := by native_decide

/-! ## Layer 1+2: translation embedding τ(i,j)=(i+1,j+1) and reachability
    transport (Step 1 of the induction). -/

/-- The translation embedding of the n-board into the interior block
    `{1..n}²` of the (n+2)-board, as a map of bulk vertex types. -/
def emb (n : ℕ) (a : {v : Square n // ¬ isCorner n v}) :
    {w : Square (n+2) // ¬ isCorner (n+2) w} :=
  ⟨(⟨a.val.1.val + 1, by have := a.val.1.isLt; omega⟩,
    ⟨a.val.2.val + 1, by have := a.val.2.isLt; omega⟩),
   by
     have h1 := a.val.1.isLt
     have h2 := a.val.2.isLt
     intro hc
     simp only [isCorner] at hc
     omega⟩

@[simp] theorem emb_fst (n a) : (emb n a).val.1.val = a.val.1.val + 1 := rfl
@[simp] theorem emb_snd (n a) : (emb n a).val.2.val = a.val.2.val + 1 := rfl

/-- Knight adjacency is translation-invariant: τ preserves edges (Step 1 core). -/
theorem emb_adj (n : ℕ) (a b : {v : Square n // ¬ isCorner n v})
    (h : (Bulk n).Adj a b) : (Bulk (n+2)).Adj (emb n a) (emb n b) := by
  have h' : KAdj n a.val b.val := h
  show KAdj (n+2) (emb n a).val (emb n b).val
  unfold KAdj at h' ⊢
  have r1 : (((emb n a).val.1.val : Int) - (emb n b).val.1.val).natAbs
          = ((a.val.1.val : Int) - b.val.1.val).natAbs := by
    rw [emb_fst, emb_fst]; omega
  have r2 : (((emb n a).val.2.val : Int) - (emb n b).val.2.val).natAbs
          = ((a.val.2.val : Int) - b.val.2.val).natAbs := by
    rw [emb_snd, emb_snd]; omega
  rw [r1, r2]; exact h'

/-- The translation embedding as a graph homomorphism Bulk n →g Bulk (n+2). -/
def embHom (n : ℕ) : Bulk n →g Bulk (n+2) where
  toFun := emb n
  map_rel' := by intro a b h; exact emb_adj n a b h

/-- Step 1: reachability transports along τ. If Bulk n is connected, all
    images of τ (the block B\EC) are mutually reachable in Bulk (n+2). -/
theorem reach_emb (n : ℕ) (hC : (Bulk n).Connected) (a b) :
    (Bulk (n+2)).Reachable (emb n a) (emb n b) :=
  (hC.preconnected a b).map (embHom n)

/-! ## Layer 3: explicit witnesses (Steps 2,3) and inductive assembly. -/

-- small product helpers (omega cannot multiply two variable natAbs)
private theorem nmul21 (a b : Int) (ha : a.natAbs = 2) (hb : b.natAbs = 1) :
    a.natAbs * b.natAbs = 2 := by rw [ha, hb]
private theorem nmul12 (a b : Int) (ha : a.natAbs = 1) (hb : b.natAbs = 2) :
    a.natAbs * b.natAbs = 2 := by rw [ha, hb]

/-- Given interior target coordinates `(ti,tj)` (a non-embedded-corner cell of
    the block, `1≤ti,tj≤n`) knight-adjacent to `w`, exhibit a neighbour of `w`
    in the image of τ (i.e. in B\EC). Steps 2 and 3 both reduce to this. -/
theorem witness_of_target (n : ℕ)
    (w : {v : Square (n+2) // ¬ isCorner (n+2) v}) (wi0 wj0 : ℕ)
    (hwi : w.val.1.val = wi0) (hwj : w.val.2.val = wj0) (ti tj : ℕ)
    (hti1 : 1 ≤ ti) (htin : ti ≤ n) (htj1 : 1 ≤ tj) (htjn : tj ≤ n)
    (hnc : ¬ ((ti = 1 ∨ ti = n) ∧ (tj = 1 ∨ tj = n)))
    (hk : ((wi0 : Int) - ti).natAbs * ((wj0 : Int) - tj).natAbs = 2) :
    ∃ u, (Bulk (n+2)).Adj w (emb n u) := by
  have hu : ¬ isCorner n ((⟨ti - 1, by omega⟩ : Fin n), (⟨tj - 1, by omega⟩ : Fin n)) := by
    simp only [isCorner]; omega
  refine ⟨⟨(⟨ti - 1, by omega⟩, ⟨tj - 1, by omega⟩), hu⟩, ?_⟩
  show KAdj (n+2) w.val (emb n ⟨(⟨ti - 1, by omega⟩, ⟨tj - 1, by omega⟩), hu⟩).val
  unfold KAdj
  rw [show (emb n ⟨(⟨ti - 1, by omega⟩, ⟨tj - 1, by omega⟩), hu⟩).val.1.val = ti from by
        simp only [emb_fst]; omega,
      show (emb n ⟨(⟨ti - 1, by omega⟩, ⟨tj - 1, by omega⟩), hu⟩).val.2.val = tj from by
        simp only [emb_snd]; omega,
      hwi, hwj]
  exact hk

/-- Witness lemma (Steps 2+3 unified): every vertex of `Bulk (n+2)` is either in
    the image of τ (the block B\EC) or has a knight neighbour there. -/
theorem mem_image_or_witness (n : ℕ) (hn : 6 ≤ n)
    (w : {v : Square (n+2) // ¬ isCorner (n+2) v}) :
    (∃ u, emb n u = w) ∨ (∃ u, (Bulk (n+2)).Adj w (emb n u)) := by
  obtain ⟨⟨wi, wj⟩, hw⟩ := w
  have hbi : wi.val < n + 2 := wi.isLt
  have hbj : wj.val < n + 2 := wj.isLt
  have hncw : ¬ ((wi.val = 0 ∨ wi.val = n + 1) ∧ (wj.val = 0 ∨ wj.val = n + 1)) := by
    intro hcc; exact hw (by simp only [isCorner]; omega)
  by_cases hImg : (1 ≤ wi.val ∧ wi.val ≤ n ∧ 1 ≤ wj.val ∧ wj.val ≤ n ∧
                   ¬ ((wi.val = 1 ∨ wi.val = n) ∧ (wj.val = 1 ∨ wj.val = n)))
  · -- in image of τ
    left
    obtain ⟨hi1, hin, hj1, hjn, hnc⟩ := hImg
    have hu : ¬ isCorner n ((⟨wi.val - 1, by omega⟩ : Fin n), (⟨wj.val - 1, by omega⟩ : Fin n)) := by
      simp only [isCorner]; omega
    refine ⟨⟨(⟨wi.val - 1, by omega⟩, ⟨wj.val - 1, by omega⟩), hu⟩, ?_⟩
    apply Subtype.ext
    simp only [emb, Prod.mk.injEq, Fin.ext_iff]
    omega
  · -- not in image: build an explicit neighbour in B\EC
    right
    by_cases hi0 : wi.val = 0
    · -- top edge: 1 ≤ wj.val ≤ n
      by_cases hjt : wj.val ≤ n - 1
      · exact witness_of_target n ⟨(wi, wj), hw⟩ wi.val wj.val rfl rfl 2 (wj.val + 1)
          (by omega) (by omega) (by omega) (by omega) (by omega)
          (nmul21 _ _ (by omega) (by omega))
      · exact witness_of_target n ⟨(wi, wj), hw⟩ wi.val wj.val rfl rfl 2 (n - 1)
          (by omega) (by omega) (by omega) (by omega) (by omega)
          (nmul21 _ _ (by omega) (by omega))
    · by_cases hin1 : wi.val = n + 1
      · -- bottom edge
        by_cases hjt : wj.val ≤ n - 1
        · exact witness_of_target n ⟨(wi, wj), hw⟩ wi.val wj.val rfl rfl (n - 1) (wj.val + 1)
            (by omega) (by omega) (by omega) (by omega) (by omega)
            (nmul21 _ _ (by omega) (by omega))
        · exact witness_of_target n ⟨(wi, wj), hw⟩ wi.val wj.val rfl rfl (n - 1) (n - 1)
            (by omega) (by omega) (by omega) (by omega) (by omega)
            (nmul21 _ _ (by omega) (by omega))
      · by_cases hj0 : wj.val = 0
        · -- left edge
          by_cases hit : wi.val ≤ n - 1
          · exact witness_of_target n ⟨(wi, wj), hw⟩ wi.val wj.val rfl rfl (wi.val + 1) 2
              (by omega) (by omega) (by omega) (by omega) (by omega)
              (nmul12 _ _ (by omega) (by omega))
          · exact witness_of_target n ⟨(wi, wj), hw⟩ wi.val wj.val rfl rfl (n - 1) 2
              (by omega) (by omega) (by omega) (by omega) (by omega)
              (nmul12 _ _ (by omega) (by omega))
        · by_cases hjn1 : wj.val = n + 1
          · -- right edge
            by_cases hit : wi.val ≤ n - 1
            · exact witness_of_target n ⟨(wi, wj), hw⟩ wi.val wj.val rfl rfl (wi.val + 1) (n - 1)
                (by omega) (by omega) (by omega) (by omega) (by omega)
                (nmul12 _ _ (by omega) (by omega))
            · exact witness_of_target n ⟨(wi, wj), hw⟩ wi.val wj.val rfl rfl (n - 1) (n - 1)
                (by omega) (by omega) (by omega) (by omega) (by omega)
                (nmul12 _ _ (by omega) (by omega))
          · -- interior block (1≤wi,wj≤n) but not in image ⇒ embedded corner
            have hEC : (wi.val = 1 ∨ wi.val = n) ∧ (wj.val = 1 ∨ wj.val = n) := by
              by_contra hcon
              exact hImg ⟨by omega, by omega, by omega, by omega, hcon⟩
            obtain ⟨hi', hj'⟩ := hEC
            rcases hi' with h1 | h1 <;> rcases hj' with h2 | h2
            · -- (1,1) → (2,3)
              exact witness_of_target n ⟨(wi, wj), hw⟩ wi.val wj.val rfl rfl 2 3
                (by omega) (by omega) (by omega) (by omega) (by omega)
                (nmul12 _ _ (by omega) (by omega))
            · -- (1,n) → (3,n-1)
              exact witness_of_target n ⟨(wi, wj), hw⟩ wi.val wj.val rfl rfl 3 (n - 1)
                (by omega) (by omega) (by omega) (by omega) (by omega)
                (nmul21 _ _ (by omega) (by omega))
            · -- (n,1) → (n-2,2)
              exact witness_of_target n ⟨(wi, wj), hw⟩ wi.val wj.val rfl rfl (n - 2) 2
                (by omega) (by omega) (by omega) (by omega) (by omega)
                (nmul21 _ _ (by omega) (by omega))
            · -- (n,n) → (n-1,n-2)
              exact witness_of_target n ⟨(wi, wj), hw⟩ wi.val wj.val rfl rfl (n - 1) (n - 2)
                (by omega) (by omega) (by omega) (by omega) (by omega)
                (nmul12 _ _ (by omega) (by omega))

/-- Inductive step (Steps 1–3): `Bulk n` connected ⇒ `Bulk (n+2)` connected. -/
theorem bulk_step (n : ℕ) (hn : 6 ≤ n) (hC : (Bulk n).Connected) :
    (Bulk (n+2)).Connected := by
  rw [SimpleGraph.connected_iff_exists_forall_reachable]
  refine ⟨emb n ⟨(⟨1, by omega⟩, ⟨2, by omega⟩), by simp only [isCorner]; omega⟩, ?_⟩
  intro w
  rcases mem_image_or_witness n hn w with ⟨u, hu⟩ | ⟨u, hadj⟩
  · rw [← hu]; exact reach_emb n hC _ u
  · exact (reach_emb n hC _ u).trans hadj.symm.reachable

/-- **Main theorem.** The bulk of the knight graph is connected for all n ≥ 6.
    Closes the open ingredient of Lemma 2.10(a) / Q(n)=3 general. -/
theorem bulk_connected_general : ∀ n, 6 ≤ n → (Bulk n).Connected := by
  intro n
  induction n using Nat.strong_induction_on with
  | _ n ih =>
    intro hn
    rcases Nat.lt_or_ge n 8 with h8 | h8
    · interval_cases n
      · exact bulk_conn_6
      · exact bulk_conn_7
    · have hprev : (Bulk (n - 2)).Connected := ih (n - 2) (by omega) (by omega)
      have hstep := bulk_step (n - 2) (by omega) hprev
      rwa [Nat.sub_add_cancel (by omega)] at hstep

end KnightTour
