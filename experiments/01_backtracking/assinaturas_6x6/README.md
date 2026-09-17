# assinaturas_6x6: linha exploratória encerrada

Tentativa inicial de caracterizar tours do 6x6 por uma **assinatura booleana**
(um vetor de predicados por solução) e procurar estrutura na distribuição
desses vetores. Os casos raros nas caudas da distribuição foram apelidados de
"unicórnios".

**Veredito: abandonada.** A distribuição de pesos das assinaturas não levou a
nenhuma restrição ou invariante que sobrevivesse. O que a substituiu foi o
espaço de ciclos sobre GF(2) (`experiments/04_espaco_ciclos/`), onde a mesma
pergunta ("que estrutura os tours compartilham?") tem resposta algébrica exata:
o invariante `Q(n) = 3`.

Arquivado aqui porque faz parte da história do projeto, não porque seja útil.

| arquivo | o que é |
|---|---|
| `assinaturas_6x6.json` | as assinaturas dos 9862 tours do 6x6 |
| `visualizar_.py` | gera `assinatura_grafico.png` |
| `plot_curva.py` | gera `curva_pesos_assinaturas.png`, distribuição por peso |
| `plot_heatmap.py` | gera `heatmap_correlacao.png`, correlação entre predicados |
| `unicornios_viz.html` | página que plotava a curva e destacava as caudas raras |

O JSON é regenerado por `tests/test_6x6_z3.py`, que continua servindo como
teste do modelo Z3 do 6x6.

A página `unicornios_viz.html` foi removida do
[visualizador](https://github.com/matheus-fsc/knight-tour-visualizer) em
setembro de 2026, por não refletir mais o entendimento atual do problema.
