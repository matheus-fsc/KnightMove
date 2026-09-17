/-
  Cycle space over GF(2).
  Incidence matrix, boundary operator, and cycle space Z₁ = ker(∂₁).
-/
import KnightTour.Basic

namespace KnightTour

-- Canonical edge ordering (lexicographic)
def lexLt (n : Nat) (a b : Square n) : Bool :=
  if a.1.val < b.1.val then true
  else if a.1.val > b.1.val then false
  else a.2.val < b.2.val

def canonEdge (n : Nat) (u v : Square n) : Square n × Square n :=
  if lexLt n u v then (u, v) else (v, u)

-- All undirected edges (canonical form, deduped)
def undirectedEdges (n : Nat) : List (Square n × Square n) :=
  let oriented := edgeListOriented n
  let canonical := oriented.map fun (u, v) => canonEdge n u v
  canonical.eraseDups

#eval (undirectedEdges 6).length   -- expect 80
#eval (undirectedEdges 8).length   -- expect 168

-- Boundary vector ∂₁(δ_v): 1 at all edges incident to v
def boundaryVector (n : Nat) (v : Square n) : List GF2 :=
  let edges := undirectedEdges n
  edges.map fun (u, w) =>
    if u == v || w == v then (1 : GF2) else (0 : GF2)

-- GF(2) vector ops
def gf2VecAdd (a b : List GF2) : List GF2 :=
  List.zipWith (· + ·) a b

def gf2Dot (a b : List GF2) : GF2 :=
  (List.zipWith (· * ·) a b).foldl (· + ·) (0 : GF2)

-- Check membership in cycle space: every vertex has even degree
def inCycleSpace (n : Nat) (vec : List GF2) : Bool :=
  (allSquares n).all fun v =>
    gf2Dot vec (boundaryVector n v) == (0 : GF2)

-- Vertex degree in an edge subset
def vertexDegreeInSubset (n : Nat) (vec : List GF2) (v : Square n) : Nat :=
  let edges := undirectedEdges n
  (List.zip edges vec).foldl (fun acc ((u, w), b) =>
    if b == (1 : GF2) && (u == v || w == v) then acc + 1 else acc) 0

-- 2-factor check
def is2Factor (n : Nat) (vec : List GF2) : Bool :=
  (allSquares n).all fun v => vertexDegreeInSubset n vec v == 2

-- Row space generators: {∂₁(δ_v) : v ∈ V}
def boundaryRows (n : Nat) : List (List GF2) :=
  (allSquares n).map fun v => boundaryVector n v

-- GF(2) Gaussian elimination to compute rank
-- Uses Id monad for mutable state
def gf2Rank (vecs : List (List GF2)) : Nat := Id.run do
  let arr := vecs.toArray.map (·.toArray)
  let m := arr.size
  if m == 0 then return 0
  let ncols := if h : 0 < m then (arr[0]).size else 0
  let mut mat := arr
  let mut rank := 0
  for col in [:ncols] do
    let mut pivotRow : Option Nat := none
    for row in [rank:m] do
      if (mat[row]!)[col]! == (1 : GF2) then
        pivotRow := some row
        break
    match pivotRow with
    | none => pure ()
    | some pr =>
      let tmp := mat[rank]!
      mat := mat.set! rank (mat[pr]!)
      mat := mat.set! pr tmp
      for row in [:m] do
        if row != rank && (mat[row]!)[col]! == (1 : GF2) then
          let newRow := Array.zipWith mat[row]! mat[rank]! (· + ·)
          mat := mat.set! row newRow
      rank := rank + 1
  return rank

-- rank(∂₁) = |V| - 1 for a connected graph
#eval gf2Rank (boundaryRows 6)  -- expect 35

-- XOR vector: 1 at positions of edges i and j
def xorVector (n : Nat) (i j : Square n × Square n) : List GF2 :=
  let edges := undirectedEdges n
  let ci := canonEdge n i.1 i.2
  let cj := canonEdge n j.1 j.2
  edges.map fun e =>
    if e == ci || e == cj then (1 : GF2) else (0 : GF2)

-- Mandatory edges of a corner
def mandatoryEdgesOf (n : Nat) (c : Square n) : List (Square n × Square n) :=
  (knightNeighborsList n c).map fun nb => canonEdge n c nb

-- All mandatory edges
def allMandatoryEdges (n : Nat) (hn : 4 ≤ n) : List (Square n × Square n) :=
  (cornersList n hn).flatMap fun c => mandatoryEdgesOf n c

-- Generate all distinct XOR pairs {v_ij} for i < j in Mand
def allXorVectors (n : Nat) (hn : 4 ≤ n) : List (List GF2) :=
  let mand := allMandatoryEdges n hn
  let indexed := mand.enum
  let pairs := indexed.flatMap fun (i, ei) =>
    (indexed.filterMap fun (j, ej) =>
      if i < j then some (xorVector n ei ej)
      else none)
  pairs.eraseDups

-- Q(n) = rank(rows ∪ xors) - rank(rows)
def computeQ (n : Nat) (hn : 4 ≤ n) : Nat :=
  let rows := boundaryRows n
  let xors := allXorVectors n hn
  let rankRows := gf2Rank rows
  let rankBoth := gf2Rank (rows ++ xors)
  rankBoth - rankRows

-- KEY VERIFICATION: Q(n) = 3
#eval computeQ 6 (by omega)   -- expect 3

-- rank(Ham) computation requires actual tours, done in Verification.lean

end KnightTour
