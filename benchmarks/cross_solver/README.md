# benchmark: cross-solver comparison

## Purpose
Cross-solver comparison harness: runs the baseline backtracking and Z3 solvers
on small boards, compares tour distributions (KL divergence) against
ground-truth, and produces reports. This is the place for **comparative**
benchmarks (one solver vs another), as opposed to a single experiment's own
results.

## Contents
| file | role |
|------|------|
| `baseline_bt_6x6.py` | naive/baseline backtracking solver on 6×6 |
| `benchmark_z3.py` | Z3-based solver benchmark |
| `gt_paths_6x6.py` | ground-truth path generation for 6×6 |
| `kl_divergence.py` | KL divergence between sampled and ground-truth distributions |
| `report.py` | aggregates the above into a report |
| `results/` | generated JSON + logs (see `results/`) |

## Usage
```bash
python benchmark/baseline_bt_6x6.py
python benchmark/benchmark_z3.py
python benchmark/report.py
```

## Results
Written to `results/` (e.g. `baseline_bt_6x6.json`,
`baseline_bt_6x6_solution_times.json`, run logs).

## Related
- Flagship solver: `../knight_tours_optimized/`
- Aggregated 8×8 parallel benchmark (flagship vs minimal_v2 vs naive):
  `../data/parallel_8x8/`
