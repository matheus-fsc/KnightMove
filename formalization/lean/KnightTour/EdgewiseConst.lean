/-
  edgewise_const_of_connected: if α : V → F₂ is constant across every edge
  of a connected graph, then α is globally constant.

  This is the ONE new lemma that Mathlib is missing for the Q(n)=3 proof.
  It feeds Lemma 2.10 (independence of corner representatives) via
  bulk_connected_general.
-/
import Mathlib.Combinatorics.SimpleGraph.Path

namespace KnightTour

/-- If `f` agrees across every edge of a connected graph, then `f` is constant.
    Proof: from `Connected.preconnected` we get `Reachable a b` for any `a b`;
    `reachable_iff_reflTransGen` recasts this as `ReflTransGen G.Adj a b`;
    induction on `ReflTransGen` propagates the edge equality. -/
theorem edgewise_const_of_connected
    {V : Type*} {G : SimpleGraph V} (hG : G.Connected)
    {α : Type*} (f : V → α) (hedge : ∀ u v, G.Adj u v → f u = f v) :
    ∀ a b, f a = f b := by
  intro a b
  have hr := (G.reachable_iff_reflTransGen a b).mp (hG.preconnected a b)
  induction hr with
  | refl => rfl
  | tail _ hadj ih => exact ih.trans (hedge _ _ hadj)

/-- The component-local version: an edgewise-constant function is constant on
    every *reachability class*, with no connectivity hypothesis on `G`.
    `edgewise_const_of_connected` is the special case where there is one class;
    this version is what the cycle-space dimension formula needs in order to
    count `#components` instead of assuming 1. -/
theorem edgewise_const_of_reachable
    {V : Type*} {G : SimpleGraph V}
    {α : Type*} (f : V → α) (hedge : ∀ u v, G.Adj u v → f u = f v)
    {a b : V} (h : G.Reachable a b) : f a = f b := by
  have hr := (G.reachable_iff_reflTransGen a b).mp h
  clear h
  induction hr with
  | refl => rfl
  | tail _ hadj ih => exact ih.trans (hedge _ _ hadj)

end KnightTour
