# η characterization — crossing efficiency vs graph structure

Dataset: 584 training + 80 held-out (graphs regenerated deterministically; η = valid_XOR / cycles_touching_path).

## Top predictors

| rank | feature | Pearson r | Spearman | p |
|---:|---|---:|---:|---:|
| 1 | treewidth_width | -0.663 | -0.786 | 3.9e-75 |
| 2 | mean_cycle_length | -0.630 | -0.702 | 5.3e-66 |
| 3 | frac_degree_2 | +0.624 | +0.657 | 3.0e-64 |
| 4 | cycle_space_dim | -0.621 | -0.768 | 1.2e-63 |
| 5 | frac_degree_4plus | -0.570 | -0.675 | 1.6e-51 |
| 6 | avg_vertex_cut | -0.545 | -0.626 | 1.5e-46 |
| 7 | frac_local_cycles | +0.541 | +0.571 | 1.0e-45 |
| 8 | mean_path_coverage | +0.514 | +0.520 | 9.7e-41 |
| 9 | mean_degree | -0.466 | -0.677 | 9.3e-33 |
| 10 | cycle_density | -0.457 | -0.670 | 1.6e-31 |

- |r|>0.5: ['treewidth_width', 'mean_cycle_length', 'frac_degree_2', 'cycle_space_dim', 'frac_degree_4plus', 'avg_vertex_cut', 'frac_local_cycles', 'mean_path_coverage']
- |r|<0.1 (useless): ['frac_degree_3', 'std_degree']

## Model comparison (5-fold CV R², 20% test RMSE)

| model | features | CV R² | test RMSE |
|---|---|---:|---:|
| frac_degree_2 (hypothesis) | frac_degree_2 | -0.067 | 0.188 |
| best single | mean_cycle_length | 0.224 | 0.192 |
| best-2 | mean_cycle_length,frac_degree_2 | 0.336 | 0.166 |
| 5-feature | top5 | -0.157 | 0.158 |
| poly deg-2 | top3 | -0.213 | 0.138 |
| random forest | all | 0.423 | 0.140 |

## Best interpretable formula

η(G) ≈ **0.8805 + -0.0097*mean_cycle_length**  [single-feature CV R²=0.224, RMSE=0.192]

Held-out validation (80 grids): MAE=0.104, RMSE=0.126, max_error=0.329.

## RF feature importances (top 6)

- treewidth_width: 0.316
- cycle_space_dim: 0.307
- mean_cycle_length: 0.084
- frac_local_cycles: 0.035
- algebraic_connectivity: 0.032
- mean_degree: 0.029

## Claim strength: **WEAK** (CV R²=0.336)

### PAPER CLAIM (exact text)

> η is **not** well predicted by simple structural features (best single-feature CV R²=0.22). The qualitative driver is real — the strongest correlate is `treewidth_width` (r=-0.66) — but no simple formula captures η reliably; we report the qualitative trend only.
