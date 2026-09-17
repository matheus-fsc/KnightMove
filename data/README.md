# data

| caminho | conteúdo | versionado |
|---|---|---|
| `archives/` | amostras brutas das Fases A/B/C, `tar.gz` | sim, via **Git LFS** |
| `results/` | resultados agregados pequenos (JSON) | sim |
| `parallel_8x8/` | benchmark cruzado 8×8 por solver | sim |
| `plots/` | gráficos gerados | sim |
| `raw/` | saídas regeneráveis e logs de execução | **não** |

## archives/

2,4 GB de JSON comprimidos para 66 MB (razão ~35x).

| arquivo | conteúdo | regenerado por |
|---|---|---|
| `fase_a_6x6_sampled.tar.gz` | Fase A, 6×6, amostragem sem simetria | `experiments/03_amostragem/cavalo_8x8/` |
| `fase_a_6x6_sym.tar.gz` | Fase A, 6×6, com quebra de simetria | idem |
| `fase_b_8x8.tar.gz` | Fase B, 8×8, re-coleta sem viés | idem |
| `fase_b_8x8_parallel.tar.gz` | Fase B, 8×8, execução paralela (4010 arquivos) | idem |
| `fase_c_10x10_parallel.tar.gz` | Fase C, 10×10: 500k amostras × 8 simetrias D4 = 4M | `experiments/03_amostragem/cavalo_10x10/` |

```bash
tar xzf data/archives/fase_c_10x10_parallel.tar.gz
```

## Não versionado

`data/raw/destruction_catalogue.json` (352 MB) é regenerável por
`experiments/01_backtracking/cavalo_loop_destruicao_6x6.py`. O resumo dele,
`destruction_map.json` (277 KB), está versionado no mesmo diretório do script.
