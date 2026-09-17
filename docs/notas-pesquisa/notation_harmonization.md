# Notation harmonization — XOR pathfinding section

The existing paper (`paper/knight_tour_complete.tex`) is in **Brazilian
Portuguese** (`babel brazilian`), uses `\bibitem`/`thebibliography` (no `.bib`),
and defines its own macros. The new section adapts to the paper (paper wins).

## Established paper notation (from preamble + body)

| concept | paper symbol | source |
|---|---|---|
| GF(2) | `\GF` → `\mathrm{GF}`, and `\Fdois` → `$\mathbb{F}_2$` | preamble L41/L47 |
| XOR / symmetric difference | `\oplus` | Fig. `fig:xor-gf2`, L160, L479 |
| cycle space | `Z_1 = \ker(\partial_1)`, `\beta_1` | `def` in `sec:grafo` |
| boundary matrix | `\partial_1` | `sec:grafo` |
| rank | `\rk` | preamble L49 |
| generic cycle (specific) | lowercase `c_1, c_2` | L160 |
| Hamiltonian path endpoints | `(s_i, e_i)` | `thm:suficiencia` |
| status tags | `\proved \verified \empirical \openp` | preamble L52–55 |

## Substitutions applied to the new section

| new-section draft used | conflict? | resolution |
|---|---|---|
| `\oplus` (XOR of edge sets) | none — paper uses `\oplus` for the same GF(2) sum | **keep `\oplus`** |
| `\mathbb{F}_2` | paper macro `\Fdois` is the house style | **`\mathbb{F}_2` → `\Fdois`** |
| `\eta` (crossing efficiency) | `\eta` is **not used** anywhere in the paper | **adopt `\eta`**, defined on first use |
| `P` (the base path) | paper does **not** use `P` for probability | **keep `P` for the path** |
| probability `P(valid \mid \dots)` | would clash with `P`=path | **use `\Pr[\cdot]`** instead of `P(\cdot)` |
| `C` (generic cycle) | paper uses lowercase `c_i` for specific cycles only | **keep `C`** for a generic cycle, no clash |
| `overlap(C,P)` | new term | `\mathrm{overlap}(C,P) := |C \cap P|`, defined on first use |
| `compat` | new term | `\mathrm{compat}`, defined on first use |
| `coverage` | new term | `\mathrm{cov}` (abbrev.), defined on first use |
| edge set of a path | — | treated as a subset of `E_n` (consistent with `Z_1 ⊆ \Fdois^{|E|}`) |

No `\renewcommand` scoping was needed: the only true clash (P-as-path vs
P-as-probability) is resolved by writing probabilities as `\Pr[\cdot]`.

## Theorem numbering

The paper's `theorem/lemma/corollary/definition/conjecture` share one counter,
numbered `[section]`. New results therefore auto-number within the new section
(Lemma `\thesection`.1, .2, …). No existing number is touched; all references
use `\label`/`\ref`.
