# kinight_tours_shared: sheared-torus utilities

> Note: directory name keeps the original (mis)spelling "kinight" to avoid
> breaking imports/paths that reference it.

## Purpose
Helper scripts for the **sheared-torus** topology line of work (knight tours on
a sheared/deformed torus and their winding invariants).

## Contents
| file | role |
|------|------|
| `knight_tours_sheared.py` | knight-tour generation on the sheared torus |
| `sheared_winding.py` | winding-number / GF(2) winding computation for sheared tours |
| `bias_check_6x6.py` | sampling-bias check on the 6×6 sheared case |

## Related
- Root-level sheared scripts: `sheared_full_search.py`,
  `sheared_heuristic_search.py`, `sheared_symmetry_analysis.py`,
  `sheared_loop_extractor.py`, `sheared_benchmark_g_infinity.py`.
- Topology siblings at root: `knight_tours_torus.py`, `knight_tours_cylinder.py`,
  `knight_tours_klein.py`.
- Docs: `../SHEARED_TORUS_README.md`.
