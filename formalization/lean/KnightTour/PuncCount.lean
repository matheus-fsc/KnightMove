/-
  (b.2) — the edge count `|E(Punc n)| = |E(G_n)| − 8`, and (b.1) — the
  identification of the corner-incident edges with the 8 mandatory edges.

  The count goes through `mandEdges_distinct` (QAbstract), which already
  supplies the injection `Fin 4 × Bool ↪ KEdge n`. What is added here is
  SURJECTIVITY onto the corner-incident edges: every edge touching a corner is
  one of the eight. That direction is where `corner_adj_iff` (corner degree
  exactly 2) does the work.

  Note on why the four corners being pairwise non-adjacent is not needed as a
  separate step: distinctness of the eight edges is exactly `mandEdges_distinct`,
  which compares the edges directly rather than reasoning about neighbourhoods.
  (Reasoning via neighbourhoods would in fact be wrong — at n = 5 the corner
  neighbourhoods intersect, at (1,2), and the edges are still distinct.)
-/
import KnightTour.PuncGraph

namespace KnightTour

open Finset

variable {n : ℕ}

/-- An edge incident to a corner: the *mandatory* edges, since corners have
    degree exactly 2. -/
def IsCornerEdge (n : ℕ) (e : KEdge n) : Prop :=
  ∃ v, v ∈ e.val ∧ isCorner n v

noncomputable instance : DecidablePred (IsCornerEdge n) := Classical.decPred _

/-- The eight mandatory edges, indexed by corner and side. -/
noncomputable def mandE (n : ℕ) (hn : 6 ≤ n) (p : Fin 4 × Bool) : KEdge n :=
  if p.2 then mandEdge1 n hn p.1 else mandEdge2 n hn p.1

theorem mandE_injective (hn : 6 ≤ n) : Function.Injective (mandE n hn) :=
  mandEdges_distinct n hn

/-- Every corner is `corner n hn k` for some `k`. -/
theorem exists_corner_index (hn : 4 ≤ n) {v : Square n} (hv : isCorner n v) :
    ∃ k : Fin 4, v = corner n hn k := by
  obtain ⟨h1, h2⟩ := hv
  rcases h1 with a | a <;> rcases h2 with b | b
  · exact ⟨⟨0, by omega⟩, Prod.ext (Fin.ext a) (Fin.ext b)⟩
  · exact ⟨⟨1, by omega⟩, Prod.ext (Fin.ext a) (Fin.ext b)⟩
  · exact ⟨⟨2, by omega⟩, Prod.ext (Fin.ext a) (Fin.ext b)⟩
  · exact ⟨⟨3, by omega⟩, Prod.ext (Fin.ext a) (Fin.ext b)⟩

/-- The mandatory edges are corner-incident. -/
theorem isCornerEdge_mandE (hn : 6 ≤ n) (p : Fin 4 × Bool) :
    IsCornerEdge n (mandE n hn p) := by
  refine ⟨corner n (by omega) p.1, ?_, corner_isCorner n (by omega) p.1⟩
  unfold mandE mandEdge1 mandEdge2 mkEdge
  cases p.2 <;> simp

/-- **Surjectivity.** Every corner-incident edge is one of the eight. -/
theorem exists_mandE_of_isCornerEdge (hn : 6 ≤ n) {e : KEdge n}
    (he : IsCornerEdge n e) : ∃ p : Fin 4 × Bool, mandE n hn p = e := by
  obtain ⟨v, hv, hcv⟩ := he
  obtain ⟨k, rfl⟩ := exists_corner_index (by omega) hcv
  obtain ⟨w, hw⟩ := Sym2.mem_iff_exists.mp hv
  have hadj : (KnightGraph n).Adj (corner n (by omega) k) w := by
    have := e.property
    rw [hw] at this
    exact (KnightGraph n).mem_edgeSet.mp this
  rcases (corner_adj_iff n hn k w).mp hadj with rfl | rfl
  · exact ⟨(k, true), Subtype.ext (by simp [mandE, mandEdge1, mkEdge, hw])⟩
  · exact ⟨(k, false), Subtype.ext (by simp [mandE, mandEdge2, mkEdge, hw])⟩

/-- There are exactly 8 corner-incident edges. -/
theorem card_cornerEdges (hn : 6 ≤ n) :
    Fintype.card {e : KEdge n // IsCornerEdge n e} = 8 := by
  classical
  have hbij : Function.Bijective
      (fun p : Fin 4 × Bool => (⟨mandE n hn p, isCornerEdge_mandE hn p⟩ :
        {e : KEdge n // IsCornerEdge n e})) := by
    constructor
    · intro a b h
      exact mandE_injective hn (congrArg Subtype.val h)
    · rintro ⟨e, he⟩
      obtain ⟨p, hp⟩ := exists_mandE_of_isCornerEdge hn he
      exact ⟨p, Subtype.ext hp⟩
  have := Fintype.card_of_bijective hbij
  simpa using this.symm

/-! ### From corner edges to `Punc` edges -/

theorem mem_punc_edgeSet_iff {e : Sym2 (Square n)} :
    e ∈ (Punc n).edgeSet ↔ ∃ he : e ∈ (KnightGraph n).edgeSet,
      ¬ IsCornerEdge n ⟨e, he⟩ := by
  induction e using Sym2.ind with
  | _ u w =>
    constructor
    · intro h
      obtain ⟨hadj, hu, hw⟩ := (Punc n).mem_edgeSet.mp h
      refine ⟨(KnightGraph n).mem_edgeSet.mpr hadj, ?_⟩
      rintro ⟨v, hv, hcv⟩
      rcases Sym2.mem_iff.mp hv with rfl | rfl
      · exact hu hcv
      · exact hw hcv
    · rintro ⟨he, hnc⟩
      have hadj : (KnightGraph n).Adj u w := (KnightGraph n).mem_edgeSet.mp he
      refine (Punc n).mem_edgeSet.mpr ⟨hadj, ?_, ?_⟩
      · intro hu; exact hnc ⟨u, by simp, hu⟩
      · intro hw; exact hnc ⟨w, by simp, hw⟩

/-- `E(Punc n)` is in bijection with the non-corner edges of `G_n`. -/
noncomputable def puncEdgeEquiv :
    (Punc n).edgeSet ≃ {e : KEdge n // ¬ IsCornerEdge n e} where
  toFun e := ⟨⟨e.val, (mem_punc_edgeSet_iff.mp e.property).choose⟩,
              (mem_punc_edgeSet_iff.mp e.property).choose_spec⟩
  invFun e := ⟨e.val.val, mem_punc_edgeSet_iff.mpr ⟨e.val.property, e.property⟩⟩
  left_inv e := by ext; rfl
  right_inv e := by ext; rfl

/-- **(b.2)** `|E(Punc n)| + 8 = |E(G_n)|`. -/
theorem card_punc_edgeSet (hn : 6 ≤ n) :
    Fintype.card (Punc n).edgeSet + 8 = Fintype.card (KEdge n) := by
  classical
  have h1 : Fintype.card (Punc n).edgeSet
      = Fintype.card {e : KEdge n // ¬ IsCornerEdge n e} :=
    Fintype.card_congr puncEdgeEquiv
  have h2 : Fintype.card {e : KEdge n // ¬ IsCornerEdge n e}
      = Fintype.card (KEdge n) - Fintype.card {e : KEdge n // IsCornerEdge n e} :=
    Fintype.card_subtype_compl _
  have h3 := card_cornerEdges (n := n) hn
  have h4 : Fintype.card {e : KEdge n // IsCornerEdge n e}
      ≤ Fintype.card (KEdge n) := Fintype.card_subtype_le _
  omega

end KnightTour
