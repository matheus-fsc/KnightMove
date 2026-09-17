/-
  FaceCycles.lean — Planar face-cycle basis (Mac Lane 2-basis).

  Mathlib v4.16.0 has NO planar-graph / face / Jordan-curve theory
  (only an incidental mention in Coloring.lean).  We therefore AXIOMATIZE a
  plane-graph embedding via its face boundaries, following Mac Lane's
  planarity criterion (1937): a graph is planar iff its cycle space has a
  2-basis in which each edge lies in ≤ 2 basis vectors.  The bounded faces of
  a planar embedding are exactly such a basis.

  The combinatorial Jordan-curve fact (a face cycle meets any path in a single
  contiguous segment) is stated as `face_contiguous_crossing` — an axiom, not
  proved here.  Proving it would require a full planarity library.
-/
import Mathlib.Combinatorics.SimpleGraph.Basic
import Mathlib.Combinatorics.SimpleGraph.Path
import KnightTour.GF2Path

namespace KnightTour

open SimpleGraph Finset

variable {V : Type*} [DecidableEq V] [Fintype V]

/-- A plane-graph embedding given by its face boundaries.  Axiomatized:
    we assume the embedding axioms as fields rather than constructing them. -/
structure PlaneGraph (V : Type*) [DecidableEq V] [Fintype V] where
  /-- the underlying simple graph -/
  graph : SimpleGraph V
  /-- each face = its set of boundary edges -/
  faces : Finset (Finset (Sym2 V))
  /-- Axiom 1: every face boundary is a cycle -/
  face_cycles : ∀ f ∈ faces, isCycle f
  /-- Axiom 2 (Mac Lane 2-basis): each edge lies in ≤ 2 face boundaries -/
  sparse_basis : ∀ e : Sym2 V, (faces.filter (fun f => e ∈ f)).card ≤ 2
  /-- Axiom 3: faces generate the cycle space.  "C = XOR of a face subset S"
      is expressed in GF(2) as: e ∈ C ↔ an odd number of chosen faces contain e
      (avoids `Finset.fold`, which would need comm/assoc instances). -/
  face_generating : ∀ C : Finset (Sym2 V), isCycle C →
    ∃ S : Finset (Finset (Sym2 V)), S ⊆ faces ∧
      ∀ e : Sym2 V, e ∈ C ↔ Odd ((S.filter (fun f => e ∈ f)).card)

/-- A face cycle is one of the embedding's face boundaries. -/
def isFaceCycle (pg : PlaneGraph V) (C : Finset (Sym2 V)) : Prop :=
  C ∈ pg.faces

/-- A face cycle is a cycle (Axiom 1). -/
theorem isFaceCycle_isCycle {pg : PlaneGraph V} {C : Finset (Sym2 V)}
    (hC : isFaceCycle pg C) : isCycle C :=
  pg.face_cycles C hC

/-- Placeholder vertex order of a path.  `contiguousOverlap` (GF2Path.lean) is
    currently abstract, so the concrete order is not needed; kept for the
    axiom's signature. -/
def pathVertexOrder (_P : Finset (Sym2 V)) (_s _t : V) : List V := []

/-- **AXIOM (combinatorial Jordan curve).**  A bounded face cycle meets any
    simple path in a single contiguous segment.  This is the topological core
    of the strong conjecture; out of scope to prove without a planarity
    library, so we axiomatize it (per task spec). -/
axiom face_contiguous_crossing
    {V : Type*} [DecidableEq V] [Fintype V]
    (pg : PlaneGraph V)
    (C P : Finset (Sym2 V)) (s t : V)
    (hC : isFaceCycle pg C)
    (hP : isValidPath P s t) :
    contiguousOverlap C P (pathVertexOrder P s t)

end KnightTour
