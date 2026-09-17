# cavalo_10x10 — 10×10 sampling engine + datasets

## Purpose
The 10×10 counterpart of `../cavalo_8x8/`: a self-contained Phase-C sampling /
correlation / destruction pipeline with its large generated dataset (~5,000
JSON files).

## Scripts (top of dir)
The `.py` files run the 10×10 parallel sampling and analysis. Run from this
directory so the relative `data_parallel/` paths resolve.

## Data layout (`data_parallel/`)
| subdir | contents |
|--------|----------|
| `samples/` | sampled 10×10 tours (JSON) |
| `correlations/` | edge/orbit correlations |
| `destruction/` | destruction-catalogue outputs |
| `results/` | aggregated results |
| `logs/` | run logs |

## Key findings (see project memory)
- Phase C complete: 500k samples × 8 D₄ = 4M; 9 pools, mean 18.6 orbits,
  7.94× compaction; universal top-1 B10-D9 ↔ B10-C8 — hypothesis confirmed.

## Note
The ~5k JSON files are kept in place (read by sibling scripts via relative
paths). Leaf data dirs are documented here rather than individually.
