/-
  Knight Tour: Basic Definitions
  Formalization of "Topologia Algébrica e Algoritmos para o Passeio do Cavalo"
  With Mathlib imports.
-/
import Mathlib.Data.Fin.Basic
import Mathlib.Data.Finset.Basic
import Mathlib.Data.ZMod.Basic
import Mathlib.Data.List.Basic
import Mathlib.Tactic

namespace KnightTour

-- GF(2) as ZMod 2 from Mathlib (has CommRing, DecidableEq, etc.)
abbrev GF2 := ZMod 2

-- A square on the n×n board
abbrev Square (n : ℕ) := Fin n × Fin n

-- The 8 knight move offsets as (dr, dc)
def knightOffsets : List (Int × Int) :=
  [(1,2),(1,-2),(-1,2),(-1,-2),(2,1),(2,-1),(-2,1),(-2,-1)]

-- Compute knight neighbors of a square on the n×n board
def knightNeighborsList (n : ℕ) (v : Square n) : List (Square n) :=
  knightOffsets.filterMap fun (dr, dc) =>
    let r' := (v.1.val : Int) + dr
    let c' := (v.2.val : Int) + dc
    if h1 : 0 ≤ r' ∧ r' < n then
      if h2 : 0 ≤ c' ∧ c' < n then
        some (⟨r'.toNat, by omega⟩, ⟨c'.toNat, by omega⟩)
      else none
    else none

-- Knight neighbors as a Finset
def knightNeighbors (n : ℕ) (v : Square n) : Finset (Square n) :=
  (knightNeighborsList n v).toFinset

-- Degree of a vertex in the knight graph
def knightDegree (n : ℕ) (v : Square n) : ℕ :=
  (knightNeighbors n v).card

-- The 4 corners of the n×n board (n ≥ 4)
def corner_00 (n : ℕ) (h : 4 ≤ n) : Square n := (⟨0, by omega⟩, ⟨0, by omega⟩)
def corner_0n (n : ℕ) (h : 4 ≤ n) : Square n := (⟨0, by omega⟩, ⟨n-1, by omega⟩)
def corner_n0 (n : ℕ) (h : 4 ≤ n) : Square n := (⟨n-1, by omega⟩, ⟨0, by omega⟩)
def corner_nn (n : ℕ) (h : 4 ≤ n) : Square n := (⟨n-1, by omega⟩, ⟨n-1, by omega⟩)

def cornersList (n : ℕ) (h : 4 ≤ n) : List (Square n) :=
  [corner_00 n h, corner_0n n h, corner_n0 n h, corner_nn n h]

-- List all squares on the n×n board
def allSquares (n : ℕ) : List (Square n) :=
  (List.finRange n).flatMap fun i =>
    (List.finRange n).map fun j => (i, j)

-- List all oriented edges (u, v) where v is a knight neighbor of u
def edgeListOriented (n : ℕ) : List (Square n × Square n) :=
  (allSquares n).flatMap fun u =>
    (knightNeighborsList n u).map fun v => (u, v)

-- Count edges (oriented count / 2 = undirected count)
def numEdges (n : ℕ) : ℕ := (edgeListOriented n).length / 2
def numVertices (n : ℕ) : ℕ := n * n

-- First Betti number: β₁ = |E| - |V| + 1
def beta1 (n : ℕ) : Int := (numEdges n : Int) - (numVertices n : Int) + 1

-- Verify basic graph parameters
#eval numVertices 6    -- expect 36
#eval numEdges 6       -- expect 80
#eval beta1 6          -- expect 45
#eval numVertices 8    -- expect 64
#eval numEdges 8       -- expect 168
#eval beta1 8          -- expect 105
#eval numVertices 10   -- expect 100
#eval numEdges 10      -- expect 288
#eval beta1 10         -- expect 189

-- Verify corner degrees
#eval knightDegree 6 (⟨0, by omega⟩, ⟨0, by omega⟩)  -- expect 2
#eval knightDegree 8 (⟨0, by omega⟩, ⟨0, by omega⟩)  -- expect 2

end KnightTour
