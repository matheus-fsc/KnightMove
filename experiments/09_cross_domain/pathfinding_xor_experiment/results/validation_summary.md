# Validation summary — XOR pathfinding (compat ≈ 0.76 × coverage)

## RISK 1 — is 0.76 grid/density-specific?

- Efficiency across all 526 instances: range [0.205, 1.000], global mean **0.671 ± 0.218**.
- All cell means in [0.70,0.82]: False; all cell CV<0.15: False; cells with CV>0.20: 6.
- **ELIMINATED: NO**

### Efficiency matrix (mean ± std)

| size＼density | 0.10 | 0.20 | 0.30 | 0.40 | 0.50 |
|---|---|---|---|---|---|
| 10×10 | 0.828±0.147 | 0.759±0.155 | 0.834±0.160 | 0.888±0.121 | 0.946±0.100 |
| 20×20 | 0.519±0.147 | 0.542±0.136 | 0.756±0.091 | 0.853±0.121 | 0.950±0.060 |
| 30×30 | 0.486±0.109 | 0.476±0.079 | 0.603±0.117 | 0.799±0.116 | skip |
| 50×50 | 0.389±0.084 | 0.396±0.097 | 0.458±0.069 | 0.729±0.095 | skip |

## RISK 2 — generalization to non-grid planar graphs

| Type | n | \|V\| | \|E\| | deg | coverage | compat | efficiency |
|---|---:|---:|---:|---:|---:|---:|---:|
| Grid+obstacles | 20 | 266 | 361 | 2.7 | 38.4% | 29.4% | 0.764±0.093 |
| Delaunay | 20 | 400 | 1181 | 5.9 | 12.9% | 7.9% | 0.594±0.303 |
| Gabriel | 20 | 400 | 756 | 3.8 | 32.6% | 11.0% | 0.342±0.104 |
| Random planar | 18 | 338 | 345 | 2.0 | 55.3% | 51.2% | 0.939±0.146 |

- Grid efficiency 0.764; off-grid effs [0.764, 0.594, 0.342, 0.939].
- Matches grid within 0.10: False; all in [0.65,0.85]: False.
- **ELIMINATED: NO**

## RISK 3 — fair comparison with Yen

| k | XOR div | Yen div | XOR len | Yen len | XOR thru | Yen thru |
|---:|---:|---:|---:|---:|---:|---:|
| 5 | 0.137 | 0.275 | 55.8 | 39.7 | 1.35 | 0.77 |
| 10 | 0.227 | 0.354 | 61.5 | 39.7 | 2.41 | 0.78 |
| 25 | 0.290 | 0.398 | 59.3 | 39.7 | 4.48 | 0.72 |
| 50 | 0.276 | 0.430 | 55.1 | 39.7 | 5.32 | 0.72 |
| 100 | 0.251 | 0.457 | 58.7 | 39.7 | 5.11 | 0.68 |

- Break-even k* (Yen diversity > XOR): k*=5.
- XOR competitive for k < k*; Yen superior for k ≥ k*.
- **ELIMINATED: YES**

## OVERALL VERDICT

- RISK 1 (0.76 grid-specific): **NO**
- RISK 2 (non-general): **NO**
- RISK 3 (unfair Yen): **YES**

**OVERALL: NEEDS QUALIFICATION**

## REQUIRED PAPER QUALIFICATIONS

**The claim `compat ≈ 0.76 × coverage` is NOT a universal law and must be restated.** Efficiency η = compat/coverage varies systematically and 0.76 is merely the value at the cell where it was first measured (20×20, 30% obstacles → 0.756/0.764).

1. **η decreases with graph size** (fixed 30% density): 10×10≈0.83, 20×20≈0.76, 30×30≈0.60, 50×50≈0.46.
2. **η increases with obstacle density** (e.g. 50×50: 0.39 @10% → 0.73 @40%; 20×20: 0.52 @10% → 0.95 @50%).
3. **η varies 0.34–0.94 across planar graph types**: Gabriel 0.34, Delaunay 0.59, grid 0.76, near-tree random-planar 0.94.
4. **Unifying driver:** η rises as the graph becomes more *corridor-like* (lower mean degree / higher obstacle density / sparser). In corridors a touched cycle crosses the path contiguously (→ valid); in open regions it crosses at separated points (→ disconnected, invalid). So η measures local corridor structure, not a constant.

**Correct statement for the paper:** `compat = η · coverage`, with η ∈ [0.34, 0.95] increasing in corridor-likeness; η ≈ 0.76 holds specifically for 20×20 grids at 30% obstacle density. The *qualitative* law (compat is the product of coverage and a structure-dependent crossing-efficiency) IS robust; the *numeric coefficient* is not.

**On Yen (RISK 3):** the comparison is now fair, and the honest result is that **Yen's k-shortest dominates XOR on diversity at every k ≥ 5** and returns shorter (higher-quality) paths; XOR's only advantage is throughput (paths/ms) and it caps at ~32 valid alternatives from a single base path. The paper must NOT claim XOR is a competitive *diverse*-path enumerator — only a cheap generator of many *local* variations.
