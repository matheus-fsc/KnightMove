# 01 — Backtracking e primeiras engines

O ponto de partida: backtracking com podas, poda inversa e Warnsdorff, mais a
primeira observação estrutural útil (paridade das cores do tabuleiro, que
descarta de imediato caminhos abertos entre casas de mesma cor).

| arquivo | o que faz |
|---|---|
| `cavalo_engine_8x8_path.py`, `cavalo_engine_8x8_path_v3.py` | engines 8×8 para caminhos com início/fim fixos |
| `cavalo_loop_destruicao_6x6.py` | catálogo de destruição: quais arestas, ao serem removidas, matam quais loops. Gera `data/raw/destruction_catalogue.json` (352 MB, não versionado); o resumo é `destruction_map.json` |
| `analise_6x6.py` | análise das assinaturas e correlações do 6×6 |
| `bishop_queen_6x6.py` | controle: a mesma análise para bispo e dama, para separar o que é do cavalo do que é do tabuleiro |
| `resultados_8x8*.json`, `summary*` | saídas das execuções |

**Veredito:** otimizações locais de backtracking não mudam a classe de
complexidade. O 8×8 continuou fora de alcance por enumeração exaustiva. O que
sobrou de valioso foi o hábito de desenhar a árvore de recursão — que levou
diretamente à descoberta dos loops (etapa 5).
