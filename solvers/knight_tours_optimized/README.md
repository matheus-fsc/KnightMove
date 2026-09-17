# knight_tours_optimized — flagship solver

## Purpose
The flagship Knight's-Tour enumeration/counting engine: cycle-space + union-find
incremental sub-tour detection, with a Numba-accelerated core and parallel
prefix decomposition. This is the fastest solver in the repo (the
`incremental_subtour` work, ~25× over Z3 on 10×10, is built on this package).

## Layout (Python package)
| file | role |
|------|------|
| `__init__.py` | package exports |
| `core.py` | pure-Python solver core (graph, backtracking, sub-tour detection) |
| `core_numba.py` | Numba-JIT accelerated hot loop |
| `divide_conquer.py` | block / divide-and-conquer decomposition |
| `parallel.py` | multiprocessing driver |
| `parallel_prefix.py` | parallel prefix-enumeration strategy |
| `cli.py` | command-line entry point |
| `benchmark.py` | self-benchmark harness |

## Usage
```bash
python -m knight_tours_optimized.cli --help
python -m knight_tours_optimized.benchmark        # run the internal benchmark
```

## Results
This package writes aggregated runs to the shared cross-solver location
`../data/parallel_8x8/flagship.json` (compared there against `minimal_v2.json`
and `naive.json`). See `../data/README.md` and `../benchmark/README.md`.

## Key findings
- Incremental union-find sub-tour detection collapses the 2-factor/tour ratio
  from ~18× to ~1× and gives ~25× speedup over Z3 on 10×10.
- Scales with near-constant nodes/tour across n∈{6,8,10,12,14}.

## Dependencies
- `numba` (for `core_numba.py`; pure-Python `core.py` is the fallback).
- Imported by several top-level benchmark scripts.
