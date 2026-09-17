/-
  The punctured knight graph `Punc n`: the knight graph with every
  corner-incident edge deleted, so the four corners become isolated vertices.

  This is the graph whose cycle space is `Z_bulk`, realised as a SPANNING
  subgraph of `KnightGraph n` (same vertex type) rather than as a graph on a
  subtype. That choice is deliberate: it lets `CycleSpaceDim.cycleSpaceOn`
  apply directly, with no transport across vertex types.

  Main result here: `card_connectedComponent_punc`, that `Punc n` has exactly
  **5** connected components for `n ≥ 6` — the bulk plus the four corners.

  The counting deliberately AVOIDS computing a cardinality of a quotient type.
  Instead it exhibits an explicit `Equiv` with `Fin 5` built from a colouring
  `compIdx : Square n → Fin 5` that is constant on components. Both halves of
  the bijectivity then live in `Prop` over concrete vertices.

  ## A tactic note worth generalising

  `cidx` encodes the corner index as `2 * (row bit) + (column bit)` rather than
  as nested `if`s. That is not cosmetic. With the arithmetic encoding plus the
  bounds `bit ≤ 1`, `omega` recovers `bit = bit'` from `2a + b = 2c + d` by
  itself; with nested `if`s it must branch, and the corner-separation argument
  turns into sixteen hand-written cases.

  **General rule: when a finite index has to be separated later, encode it
  arithmetically with bounds, not by nested conditionals.** This applies to any
  finite colouring or tag introduced for a subsequent case analysis.
-/
import KnightTour.BulkConnectivity
import KnightTour.QAbstract

namespace KnightTour

/-- The knight graph with all corner-incident edges removed. -/
def Punc (n : ℕ) : SimpleGraph (Square n) where
  Adj u v := (KnightGraph n).Adj u v ∧ ¬ isCorner n u ∧ ¬ isCorner n v
  symm := by
    rintro u v ⟨h, hu, hv⟩
    exact ⟨h.symm, hv, hu⟩
  loopless := by
    rintro u ⟨h, _, _⟩
    exact (KnightGraph n).loopless u h

instance puncDecRel (n : ℕ) : DecidableRel (Punc n).Adj := by
  intro u v; unfold Punc; simp only; infer_instance

theorem punc_le (n : ℕ) : Punc n ≤ KnightGraph n := fun _ _ h => h.1

/-! ### Corners are isolated -/

theorem punc_corner_isolated {n : ℕ} {u v : Square n} (hu : isCorner n u) :
    ¬ (Punc n).Adj u v := fun h => h.2.1 hu

/-- From a corner, nothing is reachable but the corner itself. -/
theorem punc_reachable_corner {n : ℕ} {u v : Square n} (hu : isCorner n u)
    (h : (Punc n).Reachable u v) : u = v := by
  obtain ⟨w⟩ := h
  cases w with
  | nil => rfl
  | cons hadj _ => exact absurd hadj (punc_corner_isolated hu)

/-! ### The bulk is one component -/

/-- The inclusion `Bulk n →g Punc n`. -/
def puncHom (n : ℕ) : Bulk n →g Punc n where
  toFun := Subtype.val
  map_rel' := by
    intro a b h
    have hk : KAdj n (a : Square n) (b : Square n) := h
    have hne : (a : Square n) ≠ (b : Square n) := by
      intro heq
      rw [heq] at hk
      unfold KAdj at hk
      simp at hk
    exact ⟨⟨hk, hne⟩, a.property, b.property⟩

theorem punc_reachable_of_noncorner {n : ℕ} (hn : 6 ≤ n) {u v : Square n}
    (hu : ¬ isCorner n u) (hv : ¬ isCorner n v) : (Punc n).Reachable u v := by
  have hC := bulk_connected_general n hn
  have hr : (Bulk n).Reachable ⟨u, hu⟩ ⟨v, hv⟩ := hC.preconnected _ _
  exact hr.map (puncHom n)

/-! ### The colouring -/

/-- Component index as a natural number: corners get `0,1,2,3` according to
    which side they sit on, everything else gets `4`. The additive encoding
    (rather than nested `if`s) is what makes the corner-separation argument
    fall to `omega`. -/
noncomputable def cidx (n : ℕ) (v : Square n) : ℕ :=
  if isCorner n v then
    2 * (if v.1.val = 0 then 0 else 1) + (if v.2.val = 0 then 0 else 1)
  else 4

theorem cidx_lt (n : ℕ) (v : Square n) : cidx n v < 5 := by
  unfold cidx
  split
  · split <;> split <;> omega
  · omega

/-- Component index, as an element of `Fin 5`. -/
noncomputable def compIdx (n : ℕ) (v : Square n) : Fin 5 :=
  ⟨cidx n v, cidx_lt n v⟩

theorem compIdx_noncorner {n : ℕ} {v : Square n} (hv : ¬ isCorner n v) :
    compIdx n v = 4 := by
  apply Fin.ext
  simp [compIdx, cidx, hv]

theorem compIdx_corner_ne_four {n : ℕ} {v : Square n} (hv : isCorner n v) :
    compIdx n v ≠ 4 := by
  intro h
  have : cidx n v = 4 := congrArg Fin.val h
  unfold cidx at this
  rw [if_pos hv] at this
  split at this <;> split at this <;> omega

/-- Corners with the same index are equal (needs `n ≥ 2`, so `0 ≠ n−1`). -/
theorem corner_eq_of_compIdx {n : ℕ} (hn : 2 ≤ n) {u v : Square n}
    (hu : isCorner n u) (hv : isCorner n v) (h : compIdx n u = compIdx n v) :
    u = v := by
  have hn1 : n - 1 ≠ 0 := by omega
  have h' : cidx n u = cidx n v := congrArg Fin.val h
  unfold cidx at h'
  rw [if_pos hu, if_pos hv] at h'
  obtain ⟨hu1, hu2⟩ := hu
  obtain ⟨hv1, hv2⟩ := hv
  -- the two "bits" are recovered from the value by linear arithmetic
  have b1 : (if u.1.val = 0 then 0 else 1 : ℕ) ≤ 1 := by split <;> omega
  have b2 : (if u.2.val = 0 then 0 else 1 : ℕ) ≤ 1 := by split <;> omega
  have b3 : (if v.1.val = 0 then 0 else 1 : ℕ) ≤ 1 := by split <;> omega
  have b4 : (if v.2.val = 0 then 0 else 1 : ℕ) ≤ 1 := by split <;> omega
  have e1 : (if u.1.val = 0 then 0 else 1 : ℕ)
          = (if v.1.val = 0 then 0 else 1 : ℕ) := by omega
  have e2 : (if u.2.val = 0 then 0 else 1 : ℕ)
          = (if v.2.val = 0 then 0 else 1 : ℕ) := by omega
  refine Prod.ext (Fin.ext ?_) (Fin.ext ?_)
  · rcases hu1 with a | a <;> rcases hv1 with c | c <;>
      simp only [a, c] <;> simp [a, c, hn1] at e1 ⊢
  · rcases hu2 with a | a <;> rcases hv2 with c | c <;>
      simp only [a, c] <;> simp [a, c, hn1] at e2 ⊢

/-- `compIdx` is constant on reachability classes. -/
theorem compIdx_const {n : ℕ} {u v : Square n} (h : (Punc n).Reachable u v) :
    compIdx n u = compIdx n v := by
  by_cases hu : isCorner n u
  · rw [punc_reachable_corner hu h]
  · have hv : ¬ isCorner n v := by
      intro hv
      exact hu ((punc_reachable_corner hv h.symm) ▸ hv)
    rw [compIdx_noncorner hu, compIdx_noncorner hv]

/-- The induced map on components. -/
noncomputable def compIdxLift (n : ℕ) : (Punc n).ConnectedComponent → Fin 5 :=
  Quot.lift (compIdx n) (fun _ _ h => compIdx_const h)

/-! ### The bijection with `Fin 5` -/

theorem compIdxLift_surjective {n : ℕ} (hn : 6 ≤ n) :
    Function.Surjective (compIdxLift n) := by
  intro i
  have h4 : (4 : ℕ) ≤ n := by omega
  -- the four corners and one interior vertex
  have hcorner : ∀ k : Fin 4, isCorner n (corner n h4 k) := fun k =>
    corner_isCorner n h4 k
  have hz : (0 : ℕ) < n := by omega
  have hm : n - 1 < n := by omega
  have hn1 : n - 1 ≠ 0 := by omega
  have hcor : ∀ (a b : Fin n), (a.val = 0 ∨ a.val = n - 1) →
      (b.val = 0 ∨ b.val = n - 1) → isCorner n (a, b) := by
    intro a b ha hb; exact ⟨ha, hb⟩
  fin_cases i
  · refine ⟨(Punc n).connectedComponentMk (⟨0, hz⟩, ⟨0, hz⟩), ?_⟩
    show compIdx n _ = _
    apply Fin.ext
    simp [compIdx, cidx, hcor ⟨0, hz⟩ ⟨0, hz⟩ (Or.inl rfl) (Or.inl rfl)]
  · refine ⟨(Punc n).connectedComponentMk (⟨0, hz⟩, ⟨n - 1, hm⟩), ?_⟩
    show compIdx n _ = _
    apply Fin.ext
    simp [compIdx, cidx, hcor ⟨0, hz⟩ ⟨n - 1, hm⟩ (Or.inl rfl) (Or.inr rfl), hn1]
  · refine ⟨(Punc n).connectedComponentMk (⟨n - 1, hm⟩, ⟨0, hz⟩), ?_⟩
    show compIdx n _ = _
    apply Fin.ext
    simp [compIdx, cidx, hcor ⟨n - 1, hm⟩ ⟨0, hz⟩ (Or.inr rfl) (Or.inl rfl), hn1]
  · refine ⟨(Punc n).connectedComponentMk (⟨n - 1, hm⟩, ⟨n - 1, hm⟩), ?_⟩
    show compIdx n _ = _
    apply Fin.ext
    simp [compIdx, cidx, hcor ⟨n - 1, hm⟩ ⟨n - 1, hm⟩ (Or.inr rfl) (Or.inr rfl), hn1]
    rfl
  · have hone : (1 : ℕ) < n := by omega
    refine ⟨(Punc n).connectedComponentMk (⟨1, hone⟩, ⟨1, hone⟩), ?_⟩
    have hnc : ¬ isCorner n ((⟨1, hone⟩ : Fin n), (⟨1, hone⟩ : Fin n)) := by
      unfold isCorner; simp only; omega
    show compIdx n _ = _
    rw [compIdx_noncorner hnc]
    rfl

theorem compIdxLift_injective {n : ℕ} (hn : 6 ≤ n) :
    Function.Injective (compIdxLift n) := by
  refine SimpleGraph.ConnectedComponent.ind₂ ?_
  intro u v h
  have h' : compIdx n u = compIdx n v := h
  apply SimpleGraph.ConnectedComponent.sound
  by_cases hu : isCorner n u
  · by_cases hv : isCorner n v
    · rw [corner_eq_of_compIdx (by omega) hu hv h']
    · exact absurd (h' ▸ compIdx_corner_ne_four hu) (by rw [compIdx_noncorner hv]; simp)
  · by_cases hv : isCorner n v
    · exact absurd ((compIdx_noncorner hu).symm.trans h' ▸ compIdx_corner_ne_four hv)
        (by simp)
    · exact punc_reachable_of_noncorner hn hu hv

/-- **`Punc n` has exactly five components** for `n ≥ 6`: the bulk plus the
    four isolated corners. -/
noncomputable def puncCompEquiv {n : ℕ} (hn : 6 ≤ n) :
    (Punc n).ConnectedComponent ≃ Fin 5 :=
  Equiv.ofBijective (compIdxLift n) ⟨compIdxLift_injective hn, compIdxLift_surjective hn⟩

noncomputable instance puncCompFintype {n : ℕ} (hn : 6 ≤ n) :
    Fintype (Punc n).ConnectedComponent :=
  Fintype.ofEquiv (Fin 5) (puncCompEquiv hn).symm

theorem card_connectedComponent_punc {n : ℕ} (hn : 6 ≤ n) :
    @Fintype.card (Punc n).ConnectedComponent (puncCompFintype hn) = 5 := by
  rw [Fintype.ofEquiv_card]
  simp

end KnightTour
