# results — cross-solver 8×8 benchmark

## Contents
One JSON per solver strategy, all run on the same 8×8 setup by the root-level
`benchmark_8x8_parallel.py`:

| file | strategy | description |
|------|----------|-------------|
| `flagship.json` | flagship | cycle-space + union-find solver (`../../knight_tours_optimized/`) |
| `minimal_v2.json` | minimal v2 | scaling-minimal backtracking variant |
| `naive.json` | naive | plain backtracking baseline |
| `warnsdorff.json` | Warnsdorff | Warnsdorff degree heuristic |
| `uniforme.json` | uniform | uniform-random move ordering |
| `theory.json` | theory | f∞ local-phase theoretical ordering |
| `z3_puro.json` | Z3 pure | Z3 with no extra clauses |
| `z3_mandatory.json` | Z3 + mandatory | Z3 with mandatory-edge clauses |
| `summary.json` | — | aggregated comparison across all strategies |
| `run.log` | — | run log |

## How to regenerate
```bash
python ../../benchmark_8x8_parallel.py
```
