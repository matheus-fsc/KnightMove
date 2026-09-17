# tests — regressão

```bash
for t in tests/*.py; do ./venv/bin/python "$t"; done
```

| teste | verifica |
|---|---|
| `test_gf2.py`, `test_gf2_classic.py` | álgebra do espaço de ciclos sobre GF(2) |
| `test_open_path.py` | caminhos abertos com início/fim fixos, incluindo o caso impossível de mesma cor |
| `test_6x6_z3.py` | modelo Z3 do 6×6 |
| `test_z3_dof.py` | graus de liberdade do modelo Z3 |

O teste de contagem canônico (9862 tours no 6×6) está em
`benchmarks/benchmark_exhaustive_6x6.py`.
