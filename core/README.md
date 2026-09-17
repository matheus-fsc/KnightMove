# core: engines de enumeração

Módulos importáveis por todo o resto do repositório. Cada um é autocontido:
constrói o grafo do cavalo sobre uma superfície, extrai o espaço de ciclos e
enumera tours fechados.

| módulo | superfície / técnica |
|---|---|
| `knight_tours.py` | plano `n×n`. Engine de referência: backtracking + pressão de vértice + union-find. É o que os benchmarks chamam de *baseline honesto*. |
| `knight_tours_torus.py` | toro (ambas as direções modulares). Grau uniforme 8, `Q = 0`. |
| `knight_tours_cylinder.py` | cilindro (X modular, Y rígido). `Q = 0`. |
| `knight_tours_klein.py` | garrafa de Klein. `H1 = Z`, não `Z²`: a paridade off-diagonal é sempre proibida. |
| `knight_tours_sheared.py` | toro cisalhado, parâmetro `s`. |
| `knight_tours_dnc.py` | divide-and-conquer com blocos 6×6. |
| `knight_tours_hybrid.py` | variante híbrida (plano com bordas parcialmente identificadas). |
| `knight_tours_patch.py` | `knight_tours.py` + LUT de deformação GF(2). Resultado negativo: taxa de resgate 0%. |
| `knight_tours_deformation_gf2.py` | deformação de tours por XOR de ciclos de face. |

Scripts fora de `core/` importam estes módulos pelo nome curto (`import
knight_tours`) graças ao bootstrap inserido no topo de cada um.
