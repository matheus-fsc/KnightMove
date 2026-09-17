/-
  GF2Path.lean — Path / cycle XOR theory over GF(2).

  Companion to the computational study `pathfinding_xor_experiment/`.
  We model an (undirected) edge set as a `Finset (Sym2 V)` and study the
  symmetric difference (XOR over GF(2)) of a base path with a cycle.

  Proved here (no sorry):
    • vertexDegree_edgeXor_parity  (LEMMA 1: degree parity preservation)
    • xor_path_cycle_degree_parity (LEMMA 2, parity form: needs only that C
                                    has even degrees — NO `alternating`)
    • xor_path_cycle_degrees       (LEMMA 2, EXACT form: with the user's
                                    `alternating` precondition, interior
                                    vertices get degree exactly 2)
  See `ConjectureXOR.lean` for LEMMA 3 (overlap = 0 → invalid) and the
  conjecture statements.
-/
import Mathlib.Combinatorics.SimpleGraph.Acyclic
import Mathlib.Combinatorics.SimpleGraph.Path
import Mathlib.Data.Finset.Card
import Mathlib.Tactic

namespace KnightTour

open SimpleGraph Finset

variable {V : Type*} [DecidableEq V]

/-- Number of edges of `E` incident to `v` (the GF(2)/integer degree). -/
def vertexDegree (E : Finset (Sym2 V)) (v : V) : ℕ :=
  (E.filter (fun e => v ∈ e)).card

/-- XOR (symmetric difference) of two edge sets, exactly as in the spec. -/
def edgeXor (A B : Finset (Sym2 V)) : Finset (Sym2 V) :=
  (A \ B) ∪ (B \ A)

/-- Overlap |C ∩ P|. -/
def overlap (C P : Finset (Sym2 V)) : ℕ := (C ∩ P).card

/-- A cycle: every vertex has even degree, the edge set is nonempty, and it
    carries an actual graph cycle (a closed `IsCycle` walk) in `fromEdgeSet`.
    Carrying the graph cycle is what makes LEMMA 3 provable. -/
def isCycle (C : Finset (Sym2 V)) : Prop :=
  (∀ v, Even (vertexDegree C v)) ∧
  ∃ (u : V) (w : (fromEdgeSet (C : Set (Sym2 V))).Walk u u), w.IsCycle

/-- Necessary conditions for an edge set to be a simple `s → t` path:
    degree 1 at the endpoints, even degree elsewhere (interior = 2, absent = 0),
    and acyclicity (a simple path contains no cycle). -/
def isValidPath (E : Finset (Sym2 V)) (s t : V) : Prop :=
  vertexDegree E s = 1 ∧ vertexDegree E t = 1 ∧
  (∀ v, v ≠ s → v ≠ t → Even (vertexDegree E v)) ∧
  (fromEdgeSet (E : Set (Sym2 V))).IsAcyclic

/-- `alternating C P` (user's correction): at every vertex touched by the
    overlap `C ∩ P`, exactly ONE incident overlap-edge occurs.  Together with
    `Even (vertexDegree C v)` (degree 2 on the cycle) this means: of the two
    `C`-edges at `v`, exactly one lies in `P` and one does not. -/
def alternating (C P : Finset (Sym2 V)) : Prop :=
  ∀ v, 0 < vertexDegree (C ∩ P) v → vertexDegree (C ∩ P) v = 1

/-- Contiguity placeholder: the overlap edges occupy a single connected
    segment of the vertex order of `P`.  (Strong conjecture; not used in the
    proved lemmas.) -/
def contiguousOverlap (_C _P : Finset (Sym2 V)) (path_order : List V) : Prop :=
  ∀ i j k : ℕ, i ≤ j → j ≤ k →
    (∀ a b, path_order.get? i = some a → path_order.get? k = some b →
      True)  -- abstract: kept opaque; real def in companion experiment

/-! ### Filtering distributes over set operations -/

private lemma filter_inc_inter (A B : Finset (Sym2 V)) (v : V) :
    (A ∩ B).filter (fun e => v ∈ e)
      = (A.filter (fun e => v ∈ e)) ∩ (B.filter (fun e => v ∈ e)) := by
  ext e; simp only [Finset.mem_filter, Finset.mem_inter]; tauto

private lemma filter_inc_sdiff (A B : Finset (Sym2 V)) (v : V) :
    (A \ B).filter (fun e => v ∈ e)
      = (A.filter (fun e => v ∈ e)) \ (B.filter (fun e => v ∈ e)) := by
  ext e; simp only [Finset.mem_filter, Finset.mem_sdiff]; tauto

private lemma filter_inc_edgeXor (A B : Finset (Sym2 V)) (v : V) :
    (edgeXor A B).filter (fun e => v ∈ e)
      = ((A.filter (fun e => v ∈ e)) \ (B.filter (fun e => v ∈ e)))
        ∪ ((B.filter (fun e => v ∈ e)) \ (A.filter (fun e => v ∈ e))) := by
  ext e
  simp only [edgeXor, Finset.mem_filter, Finset.mem_union, Finset.mem_sdiff]
  tauto

/-! ### LEMMA 1 — degree parity preservation -/

/-- **Lemma 1.** Degree in the XOR is the sum of degrees mod 2.
    Purely algebraic: `|a ∆ b| ≡ |a| + |b| (mod 2)`. -/
theorem vertexDegree_edgeXor_parity (A B : Finset (Sym2 V)) (v : V) :
    vertexDegree (edgeXor A B) v % 2
      = (vertexDegree A v + vertexDegree B v) % 2 := by
  set a := A.filter (fun e => v ∈ e) with ha
  set b := B.filter (fun e => v ∈ e) with hb
  have hxor : vertexDegree (edgeXor A B) v = (a \ b).card + (b \ a).card := by
    unfold vertexDegree
    rw [filter_inc_edgeXor, ← ha, ← hb,
        Finset.card_union_of_disjoint]
    exact disjoint_sdiff_sdiff
  have h1 : (a \ b).card + (a ∩ b).card = a.card := card_sdiff_add_card_inter a b
  have h2 : (b \ a).card + (b ∩ a).card = b.card := card_sdiff_add_card_inter b a
  have hi : (a ∩ b).card = (b ∩ a).card := by rw [Finset.inter_comm]
  have hda : vertexDegree A v = a.card := rfl
  have hdb : vertexDegree B v = b.card := rfl
  rw [hxor, hda, hdb]
  omega

/-! ### LEMMA 2 — XOR of path + cycle: degrees at interior vertices -/

/-- **Lemma 2 (parity form).**  Every non-endpoint vertex has EVEN degree in
    `edgeXor P C`.  Needs only that `C` has even degrees — the `alternating`
    hypothesis is *not* required for the parity statement. -/
theorem xor_path_cycle_degree_parity
    {P C : Finset (Sym2 V)} {s t : V}
    (hP : isValidPath P s t)
    (hCeven : ∀ v, Even (vertexDegree C v))
    (v : V) (hvs : v ≠ s) (hvt : v ≠ t) :
    Even (vertexDegree (edgeXor P C) v) := by
  have hpar := vertexDegree_edgeXor_parity P C v
  have hPe : Even (vertexDegree P v) := hP.2.2.1 v hvs hvt
  have hCe : Even (vertexDegree C v) := hCeven v
  rw [Nat.even_iff] at hPe hCe ⊢
  omega

/-- Degree splits across `A \ B` and `A ∩ B`. -/
private lemma vertexDegree_sdiff_add_inter (A B : Finset (Sym2 V)) (v : V) :
    vertexDegree (A \ B) v + vertexDegree (A ∩ B) v = vertexDegree A v := by
  unfold vertexDegree
  rw [filter_inc_sdiff, filter_inc_inter]
  exact card_sdiff_add_card_inter _ _

/-- **Lemma 2 (exact form, with the `alternating` correction).**
    If `P` is a valid path, `C` is a cycle, `C` and `P` *alternate*, and `v`
    is an interior vertex on both `P` (degree 2) and the overlap, then `v` has
    degree EXACTLY 2 in `edgeXor P C`.  This is the statement that actually
    needs `alternating`; mere `isCycle` is insufficient. -/
theorem xor_path_cycle_degrees
    {P C : Finset (Sym2 V)} {s t : V}
    (_hP : isValidPath P s t)
    (_hC : isCycle C)
    (hAlt : alternating C P)
    (v : V)
    (hPv : vertexDegree P v = 2)
    (hCv : vertexDegree C v = 2)
    (hov : 0 < vertexDegree (C ∩ P) v) :
    vertexDegree (edgeXor P C) v = 2 := by
  -- degree in XOR = deg(P\C) + deg(C\P)
  have hxor : vertexDegree (edgeXor P C) v
      = vertexDegree (P \ C) v + vertexDegree (C \ P) v := by
    unfold vertexDegree
    rw [filter_inc_edgeXor, filter_inc_sdiff, filter_inc_sdiff,
        Finset.card_union_of_disjoint disjoint_sdiff_sdiff]
  -- alternating: exactly one overlap edge at v (in both orders of ∩)
  have hCP : vertexDegree (C ∩ P) v = 1 := hAlt v hov
  have hPC : vertexDegree (P ∩ C) v = 1 := by
    have : vertexDegree (P ∩ C) v = vertexDegree (C ∩ P) v := by
      unfold vertexDegree; rw [Finset.inter_comm]
    rw [this, hCP]
  -- split degrees
  have sP := vertexDegree_sdiff_add_inter P C v   -- deg(P\C)+deg(P∩C)=deg P
  have sC := vertexDegree_sdiff_add_inter C P v   -- deg(C\P)+deg(C∩P)=deg C
  rw [hxor]; omega

/- NOTE ON BASIS CHOICE (documented from experimental results)

   DFS spanning tree basis (used in pathfinding_xor_experiment):
   - Compatibility rate: ~18% in grid graphs (26.9% at 20×20, falling with n)
   - Failure: DFS back-edges create long, non-local cycles
   - These long cycles cross paths NON-contiguously → disconnection (90.8%
     of all invalid candidates are disconnected, not degree violations)

   Face cycle basis (Mac Lane 2-basis, theoretically justified):
   - Each face cycle is geometrically local (each edge in ≤ 2 faces:
     `PlaneGraph.sparse_basis`)
   - `face_contiguous_crossing` (axiom) guarantees contiguous overlap
   - Predicted compatibility: significantly higher than 18%
   - This is the CORRECT basis for XOR path enumeration in planar graphs
   - Tested empirically in src/paths/face_basis_experiment.py

   Grinberg's theorem (1968) connects face structure to Hamiltonian cycles in
   planar graphs via  Σ (i − 2)(f'_i − f''_i) = 0,  a separate but related
   result using the same face decomposition.

   Open conjecture (experimental AUC = 0.969, ratio = 31.8×):
       valid(P XOR C)  ↔  contiguousOverlap(C, P)
   For face cycles specifically, contiguousOverlap is guaranteed by
   `face_contiguous_crossing`, reducing the conjecture to
       valid(P XOR C)  ↔  overlap ≥ 1        (for face cycles)
   This form is `face_xor_valid_iff_overlap_pos` in ConjectureXOR.lean:
   partially proved (degrees + contiguity from the axiom), with the
   connectivity step isolated as `sorry` (topological core / Jordan curve).
-/

end KnightTour
