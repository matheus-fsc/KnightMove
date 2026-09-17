# XOR pathfinding section — writing / integration notes

## Deviations from the task brief (inputs that did not exist)

The brief assumed a `paper/draft/` layout with a pre-written
`section_xor_pathfinding.tex` and a `.bib` file. **None existed.** Actual state:

- Main paper: `paper/knight_tour_complete.tex` (monolithic, pt-BR).
- No `section_xor_pathfinding.tex` → the section was **authored here** from the
  established experimental + Lean results of the XOR-pathfinding series
  (`pathfinding_xor_experiment/`, `knight_tour_lean/KnightTour/GF2Path.lean`).
- No `.bib` → new references added as `\bibitem` in the existing
  `thebibliography`.
- Paper is **Portuguese** → the section is written in Portuguese to match.
- A Lean section already exists (`sec:lean`) → updated, not created.
- Paper compiles **from repo root** (figures use root-relative paths), so the
  section was inserted **directly** into the `.tex` (no `\input`).

## Cross-reference resolution (Phase 5)

| brief placeholder | concept | resolved to | status |
|---|---|---|---|
| `thm:local_sufficiency` | Local Sufficiency Theorem | `thm:suficiencia` | EXISTS (L950) |
| `thm:theorem_u` / `thm:tightness_deficit` | Theorem U | `thm:U` | EXISTS (L512) |
| `thm:h1_conservation` | H₁(G;F₂)=45 cycle space | **`sec:grafo` + `tab:grafo`** | no dedicated theorem; β₁=45 is established in the cycle-space subsection and Table `tab:grafo`. Referenced those instead of a (non-existent) "conservation theorem". No `\todo` needed. |
| `sec:dc_algorithm` | Divide-and-Conquer | `sec:dnc` | EXISTS (L944) |

All four concepts exist; **no `\todo{}` markers were required**.

## Honesty constraints honored in the prose

- η is presented as a **qualitative** driver only (corridor-likeness); the
  cross-validated R² of every simple formula is < 0.5, and the `frac_degree_2`
  hypothesis is reported as **falsified** (negative CV R²). No formula is
  claimed.
- XOR is **not** claimed competitive with Yen for diverse enumeration; Yen wins
  diversity at all k ≥ 5. XOR's niche (high-throughput local variations,
  corridor-like graphs) is stated as such.
- Lemma X.3 (overlap=0 ⇒ invalid) is the one fully-general, Lean-checked claim
  and is presented as the section's anchor result.

## Lean (Phase 7)

Added subsection `sec:lean-xor` to `sec:lean` with a table of the three proved
lemmas (`vertexDegree_edgeXor_parity`, `xor_path_cycle_degrees`,
`overlap_zero_implies_invalid`) + the two open `sorry`s (`strong_xor_conjecture`,
`face_xor_valid_iff_overlap_pos`), file `GF2Path.lean` / `ConjectureXOR.lean`.
The pre-existing build-state box (`22 teoremas, 4 sorry`) was **not edited**
(hard constraint: add only); the new subsection states the additional counts
locally.
