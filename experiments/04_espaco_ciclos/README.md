# 04. Espaço de ciclos sobre GF(2)

A virada conceitual. Contando os "loops" da árvore de recursão, o número bateu
com a dimensão do espaço de ciclos de um grafo:

```
beta_1 = |E| - |V| + 1
```

Tours passam a ser vetores em `F_2^E`, e o conjunto de tours vive num coset do
espaço de ciclos. Cada ciclo fundamental é um grau de liberdade.

| diretório / arquivo | conteúdo |
|---|---|
| `board_10x10/` | pipeline GF(2) completo no 10×10: `V=100, E=288, beta_1=189`, rank 186 de 189, 8 arestas obrigatórias |
| `6x6_higher_order/` | exclusões de ordem 3 e superior: 88 pares, 1776 triplas, 17004 quadras minimais; `rank(T) = 42 < beta_1 = 45` |
| `forbidden_cycles_6x6/` | os 3 ciclos proibidos, e o detector dual mínimo `x_{F6-D5} XOR x_{B3-A1}` |
| `cycle_census_6x6/` | censo de ciclos do 6×6 |
| `cycle_space_hyperplanes/` | as restrições GF(2) lidas como hiperplanos em `Z_1` |
| `face_basis/` | base de ciclos de face |
| `b6_*`, `edge_compatibility_*` | correlações entre arestas e compatibilidade local |

**Descoberta central:** o rank do espaço gerado pelos tours fica sistematicamente
3 abaixo de `beta_1`. Esse deficit é o assunto de `06_invariante_Q/`.

**Limite descoberto:** a obstrução `rank < beta_1` **não** é hamiltonicidade: o rank dos 2-fatores é o mesmo 42 no 6×6. O que falta é conectividade.
