/-
  ConjectureXOR.lean — the XOR path-validity conjecture.

  Status:
    • LEMMA 3 (overlap_zero_implies_invalid) — PROVED (no sorry).
    • weak_xor_monotone — REFUTED computationally; stated as the corrected
      threshold fact (overlap = 0 ⇒ never valid), not as monotonicity.
    • strong_xor_conjecture — OPEN; stated with `sorry`, documented with the
      experimental evidence.

  Experimental support (pathfinding_xor_experiment/results/):
    - P(valid | C∩P contiguous)      = 81.8%
    - P(valid | C∩P not contiguous)  =  2.6%   (ratio ≈ 32×)
    - overlap = 0  ⇒  0.0% valid     (n = 14325 candidates) ← LEMMA 3
    - AUC(overlap+contiguous) = 0.969 vs AUC(overlap+length) = 0.862
    - 'monotone increasing in overlap' is FALSE: 0%→80%→72%→62%→27%.
-/
import KnightTour.GF2Path
import KnightTour.FaceCycles

namespace KnightTour

open SimpleGraph Finset

variable {V : Type*} [DecidableEq V]

/-! ### LEMMA 3 — disjoint overlap ⇒ invalid (the deterministic lower bound) -/

/-- If the cycle `C` is edge-disjoint from the path `P` (`overlap = 0`), then
    `C ⊆ edgeXor P C`. -/
theorem subset_edgeXor_of_disjoint
    {P C : Finset (Sym2 V)} (hdis : C ∩ P = ∅) :
    C ⊆ edgeXor P C := by
  intro e he
  have heP : e ∉ P := by
    intro hp
    have : e ∈ C ∩ P := Finset.mem_inter.mpr ⟨he, hp⟩
    rw [hdis] at this
    exact absurd this (Finset.not_mem_empty e)
  simp only [edgeXor, Finset.mem_union, Finset.mem_sdiff]
  exact Or.inr ⟨he, heP⟩

/-- **Lemma 3.**  If `P` is a valid path and `C` is a cycle edge-disjoint from
    `P` (overlap = 0), then `edgeXor P C` is NOT a valid path.

    Proof: `C ⊆ edgeXor P C`, so the graph cycle carried by `C` transfers
    (via `IsCycle.mapLe` along `fromEdgeSet_mono`) to a cycle in
    `fromEdgeSet (edgeXor P C)`, contradicting its acyclicity.  This is the
    deterministic core matching the experimental `overlap = 0 ⇒ 0% valid`. -/
theorem overlap_zero_implies_invalid
    {P C : Finset (Sym2 V)} {s t : V}
    (_hP : isValidPath P s t)
    (hC : isCycle C)
    (hdis : C ∩ P = ∅) :
    ¬ isValidPath (edgeXor P C) s t := by
  intro hV
  obtain ⟨_, u, w, hw⟩ := hC
  have hsub : C ⊆ edgeXor P C := subset_edgeXor_of_disjoint hdis
  have hcoe : (C : Set (Sym2 V)) ⊆ (edgeXor P C : Set (Sym2 V)) :=
    Finset.coe_subset.mpr hsub
  have hle : fromEdgeSet (C : Set (Sym2 V))
      ≤ fromEdgeSet (edgeXor P C : Set (Sym2 V)) := fromEdgeSet_mono hcoe
  have hcyc : (w.mapLe hle).IsCycle := hw.mapLe hle
  exact hV.2.2.2 (w.mapLe hle) hcyc

/-! ### Conjecture statements (open; documented with experimental evidence) -/

/-- WEAK form — **corrected**.  The literal "P(valid) monotone increasing in
    overlap" is REFUTED (0%→80%→72%→62%→27%).  What is actually true and
    deterministic is the overlap-0 lower bound, which is exactly `LEMMA 3`
    above.  We therefore record the corrected weak fact as a theorem rather
    than a conjecture: overlap = 0 forbids validity. -/
theorem weak_xor_threshold
    {P C : Finset (Sym2 V)} {s t : V}
    (hP : isValidPath P s t) (hC : isCycle C)
    (hdis : C ∩ P = ∅) :
    ¬ isValidPath (edgeXor P C) s t :=
  overlap_zero_implies_invalid hP hC hdis

/-- STRONG conjecture (OPEN).  Conditioned on `overlap ≥ 1`, validity of
    `edgeXor P C` is equivalent to `C ∩ P` forming a single contiguous segment
    of `P`.  Experimental support: P(valid|contiguous)=81.8% vs 2.6%, AUC 0.969.
    Left as `sorry` deliberately — this is open. -/
theorem strong_xor_conjecture
    {P C : Finset (Sym2 V)} {s t : V}
    (hP : isValidPath P s t) (hC : isCycle C)
    (path_order : List V)
    (h_overlap : 1 ≤ overlap C P) :
    isValidPath (edgeXor P C) s t ↔ contiguousOverlap C P path_order := by
  sorry

/-! ### Face-cycle specialization (Mac Lane basis) — PART 2 & PART 3 -/

/-- **PART 3 (proved, no sorry).**  Mac Lane local-perturbation bound: every
    edge lies in ≤ 2 face boundaries, so a face flip touches a bounded
    neighborhood.  Immediate from the `PlaneGraph` sparse-basis axiom. -/
theorem face_cycle_local_perturbation
    {V : Type*} [DecidableEq V] [Fintype V]
    (pg : PlaneGraph V) (_C : Finset (Sym2 V))
    (_hC : isFaceCycle pg _C) (e : Sym2 V) :
    (pg.faces.filter (fun f => e ∈ f)).card ≤ 2 :=
  pg.sparse_basis e

/-- **PART 2 corollary (proved, no sorry).**  A face XOR with overlap 0 is
    never a valid path — direct from LEMMA 3. -/
theorem face_xor_invalid_iff_disjoint
    {V : Type*} [DecidableEq V] [Fintype V]
    (pg : PlaneGraph V) (P C : Finset (Sym2 V)) (s t : V)
    (hP : isValidPath P s t) (hC : isFaceCycle pg C)
    (h_disjoint : C ∩ P = ∅) :
    ¬ isValidPath (edgeXor P C) s t :=
  overlap_zero_implies_invalid hP (isFaceCycle_isCycle hC) h_disjoint

/-- **PART 2 main theorem.**  For a FACE cycle with overlap ≥ 1, `edgeXor P C`
    is a valid path.  The `face_contiguous_crossing` axiom supplies contiguity
    and `xor_path_cycle_degrees` (LEMMA 2b) supplies exact interior degrees;
    the remaining step — deriving global connectivity / acyclicity from
    contiguous overlap — is the topological core and is left as `sorry`. -/
theorem face_xor_valid_iff_overlap_pos
    {V : Type*} [DecidableEq V] [Fintype V]
    (pg : PlaneGraph V) (P C : Finset (Sym2 V)) (s t : V)
    (hP : isValidPath P s t) (hC : isFaceCycle pg C)
    (hAlt : alternating C P) (h_overlap : 1 ≤ overlap C P) :
    isValidPath (edgeXor P C) s t := by
  -- Step 1: face cycles cross paths contiguously (combinatorial Jordan curve)
  have _hCont := face_contiguous_crossing pg C P s t hC hP
  -- Step 2: alternating + LEMMA 2b ⇒ exact degree 2 at interior overlap vertices
  have _hCyc : isCycle C := isFaceCycle_isCycle hC
  have _hDeg : ∀ v, vertexDegree P v = 2 → vertexDegree C v = 2 →
      0 < vertexDegree (C ∩ P) v → vertexDegree (edgeXor P C) v = 2 :=
    fun v hPv hCv hov => xor_path_cycle_degrees hP _hCyc hAlt v hPv hCv hov
  -- Step 3: contiguous overlap + correct degrees ⇒ connected & acyclic.
  sorry  -- connectivity from contiguous overlap — topological core, requires
         -- Jordan curve theorem; open in this formalization.

end KnightTour
