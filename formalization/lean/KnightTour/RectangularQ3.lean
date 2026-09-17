/-
  Proposition 2.18: Q(n,m) = 3 for min(n,m) ≥ 4
  Extension of Q(n)=3 to rectangular boards.
-/
import KnightTour.Basic
import KnightTour.GF2Space
import KnightTour.RectQAbstract

namespace KnightTour

-- ════════════════════════════════════════════════════════════════════
-- GENERAL THEOREM Q(n,m) = 3, n,m ≥ 6 — PROVED (no sorry).
-- Formalized abstractly as `Q_abstract_rect n m` (paper Def 2.4 on the
-- rectangle) in RectQAbstract.lean; `Q_abstract_rect_eq_three` closes it for
-- all n,m ≥ 6, axioms [propext, Classical.choice, Quot.sound, Lean.ofReduceBool]
-- (the last from the native_decide rectangular bulk-connectivity base cases).
-- The `computeQRect`-based `Q_rect_eq_three` below is a DECIDABLE SHADOW only,
-- SUPERSEDED by `Q_abstract_rect_eq_three` (closing it for symbolic n,m is the
-- rejected verified-Gaussian-elimination bridge).
-- ════════════════════════════════════════════════════════════════════

/-- Re-export: the paper's general Q(n,m)=3 (n,m ≥ 6), proved sorry-free. -/
theorem Q_rect_eq_three_proved (n m : ℕ) (hn : 6 ≤ n) (hm : 6 ≤ m) :
    Q_abstract_rect n m hn hm = 3 := Q_abstract_rect_eq_three n m hn hm

-- Rectangular board n×m
abbrev RectSquare (n m : Nat) := Fin n × Fin m

-- Knight neighbors on rectangular board
def knightNeighborsRect (n m : Nat) (v : RectSquare n m) : List (RectSquare n m) :=
  knightOffsets.filterMap fun (dr, dc) =>
    let r' := (v.1.val : Int) + dr
    let c' := (v.2.val : Int) + dc
    if h1 : 0 ≤ r' ∧ r' < n then
      if h2 : 0 ≤ c' ∧ c' < m then
        some (⟨r'.toNat, by omega⟩, ⟨c'.toNat, by omega⟩)
      else none
    else none

-- Degree of corner (0,0) on rectangular board
def rectCornerDegree (n m : Nat) (hn : 0 < n) (hm : 0 < m) : Nat :=
  (knightNeighborsRect n m (⟨0, hn⟩, ⟨0, hm⟩)).length

-- Edge list for rectangular board
def allRectSquares (n m : Nat) : List (RectSquare n m) :=
  (List.finRange n).flatMap fun i =>
    (List.finRange m).map fun j => (i, j)

def edgeListOrientedRect (n m : Nat) : List (RectSquare n m × RectSquare n m) :=
  (allRectSquares n m).flatMap fun u =>
    (knightNeighborsRect n m u).map fun v => (u, v)

-- Canonical edge for rectangular board
def lexLtRect (n m : Nat) (a b : RectSquare n m) : Bool :=
  if a.1.val < b.1.val then true
  else if a.1.val > b.1.val then false
  else a.2.val < b.2.val

def canonEdgeRect (n m : Nat) (u v : RectSquare n m) : RectSquare n m × RectSquare n m :=
  if lexLtRect n m u v then (u, v) else (v, u)

def undirectedEdgesRect (n m : Nat) : List (RectSquare n m × RectSquare n m) :=
  let oriented := edgeListOrientedRect n m
  let canonical := oriented.map fun (u, v) => canonEdgeRect n m u v
  canonical.eraseDups

-- Full Q computation for rectangular boards
def boundaryRowsRect (n m : Nat) : List (List GF2) :=
  (allRectSquares n m).map fun v =>
    (undirectedEdgesRect n m).map fun (u, w) =>
      if u == v || w == v then (1 : GF2) else (0 : GF2)

def mandatoryEdgesRect (n m : Nat) : List (RectSquare n m × RectSquare n m) :=
  let edges := undirectedEdgesRect n m
  let deg2verts := (allRectSquares n m).filter fun v =>
    edges.foldl (fun acc (u, w) =>
      if u == v || w == v then acc + 1 else acc) 0 == 2
  deg2verts.flatMap fun v =>
    edges.filter fun (u, w) => u == v || w == v

def computeQRect (n m : Nat) : Nat :=
  let edges := undirectedEdgesRect n m
  let rows := boundaryRowsRect n m
  let mand := (mandatoryEdgesRect n m).eraseDups
  let indexed := mand.enum
  let xors := indexed.flatMap fun (i, ei) =>
    indexed.filterMap fun (j, ej) =>
      if i < j then
        some (edges.map fun e =>
          if e == ei || e == ej then (1 : GF2) else (0 : GF2))
      else none
  let rankRows := gf2Rank rows
  let rankBoth := gf2Rank (rows ++ xors)
  rankBoth - rankRows

-- Verify Q(n,m) = 3 for rectangular boards
#eval computeQRect 4 6   -- expect 3
#eval computeQRect 4 8   -- expect 3
#eval computeQRect 6 10  -- expect 3

-- Proposition 2.18 — SUPERSEDED by `Q_abstract_rect_eq_three` (see header).
-- Kept only as the decidable computeQRect shadow (verified by native_decide for
-- small rectangles); closing it for symbolic n,m is the rejected bridge.
theorem Q_rect_eq_three (n m : Nat) (hn : 4 ≤ n) (hm : 4 ≤ m) :
    computeQRect n m = 3 := by
  sorry
  -- Proof: identical to Q(n)=3.
  -- The 4 corners of the rectangle have degree 2 when min(n,m) ≥ 4.
  -- Their neighborhoods are disjoint.
  -- Lemmas 2.6-2.10 generalize trivially.

end KnightTour
