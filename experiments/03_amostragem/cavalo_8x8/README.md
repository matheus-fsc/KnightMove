# cavalo_8x8: 8×8 sampling engine + datasets

## Purpose
The 8×8 (and 6×6 control) **tour-sampling and correlation** pipeline. This is a
self-contained experiment unit: the engine scripts live here alongside their
large generated datasets (~7,400 JSON files). It produced the Phase A/B data
behind the orbit-preservation results.

## Scripts (top of dir)
The `.py` files here run sampling, D₄ expansion, correlation, and destruction
analysis. Run them from this directory so their relative `data_*/` output paths
resolve.

## Data layout (`data_*/` subtrees)
Each dataset directory follows the same internal structure:

| subdir | contents |
|--------|----------|
| `samples/` | sampled tours (JSON) |
| `correlations/` | edge/orbit correlation outputs |
| `logs/` | run logs |
| `results/` | aggregated results (where present) |
| `destruction/` | destruction-catalogue outputs (parallel runs) |

Dataset variants present:
- `data/`, `data_parallel/`, main 8×8 runs
- `data_6x6_sampled/`, `data_6x6_sym/`, `data_6x6_smoke/`, 6×6 controls (incl. symmetry & smoke-test)
- `data_8x8_phase_b/`, `data_8x8_phase_b_smoke/`, Phase B 8×8 re-collection
- `configs/`: run configurations

## Key findings (see project memory)
- Phase A: Z3 sym=False VALID (ρ=0.96, 7/7 orbits); sym=True PARTIAL (ρ=0.89).
- Phase B: 100% of top-13 orbits preserved; B8-D7↔B8-C6 universal top-1.

## Note
These ~7.4k JSON files are kept in place (read by sibling scripts via relative
paths). Per-dataset leaf dirs (`samples/`, `logs/`, …) are documented here in
the parent rather than with individual READMEs, to avoid noise.
