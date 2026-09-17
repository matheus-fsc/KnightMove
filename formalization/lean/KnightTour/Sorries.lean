/-
  Inventory of all sorry's with justification.
  Documents what is proved, what is sketched, and what remains.
-/

namespace KnightTour

/-
═══════════════════════════════════════════════════════
STATUS OF THE LEAN 4 FORMALIZATION (with Mathlib v4.16.0)
═══════════════════════════════════════════════════════

PROVED (no sorry, verified by type-checker):
  ✓ Lemma 2.6: corner degree = 2 for n ∈ {4,5,6,8,10}
    (11 theorems via native_decide, all 4 corners for n=6,8)
  ✓ Q(n) = 3 for n ∈ {4,5,6,7,8,9,10}
    (7 theorems via native_decide on computeQ — full GF(2)
     Gaussian elimination checked by Lean's kernel)
  ✓ Lemma 2.7 for n=6: boundary of corner = XOR of mandatory edges
  ✓ |XOR vectors| = 28 for n=6
  ✓ rank(∂₁) = 35 for n=6
  ✓ Theorem U for n=6: all 5 values of k verified computationally
    k=0→Q=0, k=1→Q=0, k=2→Q=1, k=3→Q=2, k=4→Q=3
  ✓ Q(n,m) = 3 for (4,6), (4,8), (6,10) rectangles
  ✓ Torus degree = 8 for n ∈ {6,8,10}
  ✓ |Mand| = 8 for n ∈ {6,8,10}
  ✓ Bulk connectivity (G_n \ Mand) for n ∈ {6,8}:
    BFS reaches all 32/60 non-corner vertices respectively
  ✓ Graph parameters match paper: V, E, β₁ for n ∈ {6,8,10}
  ✓ Q(n) = 3 GENERAL: ∀ n ≥ 6, Q_abstract n = 3
    (QAbstract.lean: Q_abstract_eq_three — NO sorry; paper Thm 2.5 / Def 2.4).
    Abstract ROUTE 1: coboundary LinearMap + cycle-space quotient, corner
    representatives, independence of the 4 reps via edgewise_const_of_connected
    + bulk_connected_general, XOR image, σ-free 3-vector endgame.
    Axioms: propext, Classical.choice, Quot.sound, Lean.ofReduceBool.
    The ONE Mathlib-absent lemma proved: edgewise_const_of_connected
    (EdgewiseConst.lean — locally-constant on a connected graph ⇒ constant).
  ✓ Bulk connectivity GENERAL: ∀ n ≥ 6, (Bulk n).Connected
    (BulkConnectivity.lean: bulk_connected_general — NO sorry).
    SimpleGraph formulation; induction n→n+2, bases {6,7} (native_decide),
    step via translation-invariance hom (emb_adj, pure kernel) + explicit
    knight-move witnesses (mem_image_or_witness). Axioms: propext,
    Classical.choice, Quot.sound, Lean.ofReduceBool (last from native_decide
    base cases). This is the paper's open Lemma 2.10(a) ingredient, now closed.

SORRY 1: corner_degree_eq_two_general (Lemma1.lean:62)
  WHAT: knightDegree n (0,0) = 2 for arbitrary n ≥ 4
  WHY: filterMap with dite over 8 offsets is opaque for symbolic n
  MITIGATION: proved for n ∈ {4,5,6,8,10} via native_decide

SORRY 2: Q_eq_three_general (TheoremQ3.lean) — SUPERSEDED (no longer the
    general theorem; kept only as the computeQ decidable shadow).
  WHAT (this sorry): computeQ n = 3 for symbolic n (imperative GF(2) elim).
  REAL GENERAL THEOREM (PROVED, no sorry): Q_abstract_eq_three (QAbstract.lean),
    ∀ n ≥ 6, Q_abstract n = 3, where Q_abstract n is the paper's Def 2.4
    (= dim_F₂ of π(Span{v_ij}) in the cycle-space quotient F₂^E / row(∂₁)).
    Re-exported as Q_eq_three_general_proved (TheoremQ3.lean).
    #print axioms: [propext, Classical.choice, Quot.sound, Lean.ofReduceBool]
    (Lean.ofReduceBool from the native_decide bulk-connectivity base cases).
    The full abstract pipeline (ROUTE 1 of the blueprint) is now built:
      (a) coboundary δ : F₂^V →ₗ F₂^E as a LinearMap; R = range δ; quotient. ✓
      (b) Lemma 2.6 abstract (corner_adj_iff), Lemma 2.7 (mandEdge_xor_in_R),
          Cor 2.8, Lemma 2.10 independence (corner_reps_linearIndependent) via
          edgewise_const_of_connected + bulk_connected_general, Lemma 2.11
          (xor_image_span), σ-free endgame, finrank = 3. ✓
  WHY this sorry remains: (c) computeQ n = Q_abstract n (concrete Gaussian elim
    equals quotient dim) is the rejected Route-2 bridge — a verified-imperative-
    linear-algebra project that buys nothing. computeQ stays as a decidable
    cross-check, proved for n ∈ {4,5,6,7,8,9,10} via native_decide.

SORRY 3: theorem_U_general (TheoremU.lean:119)
  WHAT: Q(G(H)) = max(0, k(H)-1) for all n ≥ 6
  WHY: depends on Sorry 2 + rowspace monotonicity
  MITIGATION: all 5 values of k verified for n=6

SORRY 4: Q_rect_eq_three (RectangularQ3.lean) — SUPERSEDED.
  WHAT (this sorry): computeQRect n m = 3 for symbolic n,m (decidable shadow).
  REAL GENERAL THEOREM (PROVED, no sorry): Q_abstract_rect_eq_three
    (RectQAbstract.lean), ∀ n,m ≥ 6, Q_abstract_rect n m = 3 (paper Prop 2.18).
    Full duplication of the square abstract pipeline with rectangular
    coordinates, fed by bulkR_connected (RectBulkConnectivity.lean — the one
    genuinely-new ingredient: rectangular bulk connectivity, n,m ≥ 6).
    Re-exported as Q_rect_eq_three_proved (RectangularQ3.lean).
    #print axioms: [propext, Classical.choice, Quot.sound, Lean.ofReduceBool].
  WHY this sorry remains: same computeQRect = quotient-dim bridge as Sorry 2,
    deliberately not built; computeQRect stays as native_decide cross-check.

GF2 PATH XOR THEORY (GF2Path.lean, ConjectureXOR.lean):
  ✓ LEMMA 1  vertexDegree_edgeXor_parity — deg(A⊕B) ≡ deg A + deg B (mod 2).
             Fully proved (|a∆b| ≡ |a|+|b| via card_sdiff_add_card_inter).
  ✓ LEMMA 2a xor_path_cycle_degree_parity — interior vertices of P⊕C have
             EVEN degree. Needs only Even-degree of C (NO `alternating`).
  ✓ LEMMA 2b xor_path_cycle_degrees — EXACT degree 2 at interior overlap
             vertices. Requires the `alternating C P` precondition (the
             user's correction: `isCycle` alone is insufficient). Proved.
  ✓ LEMMA 3  overlap_zero_implies_invalid — C ∩ P = ∅ ⇒ ¬ isValidPath(P⊕C).
             Proved via: C ⊆ P⊕C ⇒ the graph cycle of C transfers
             (IsCycle.mapLe ∘ fromEdgeSet_mono) into P⊕C, contradicting the
             acyclicity in isValidPath. Matches experiment: overlap=0 ⇒ 0%.
  (all four depend only on propext/Classical.choice/Quot.sound — no sorryAx)

  SORRY (intentional, open): strong_xor_conjecture (ConjectureXOR.lean)
    valid(P⊕C) ↔ C∩P contiguous segment of P, for overlap ≥ 1.
    Experimental support: P(valid|contig)=81.8% vs 2.6%, AUC 0.969.
    NOTE: the literal WEAK 'monotone increasing in overlap' is REFUTED
    computationally (0%→80%→72%→62%→27%); only the overlap=0 lower bound
    (= Lemma 3) is a theorem.

OUTSIDE LEAN SCOPE:
  □ Tightness: rank(Ham(G(H))) = β₁ − Q(H)
  □ deficit_torus = 0
  □ N_DnC(12×12)
  □ Klein bottle tours
  □ Knuth estimator
-/

end KnightTour
