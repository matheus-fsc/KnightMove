# Path selection for XOR enumeration — final experiment

Grid 20×20, 20 instances (mesmas do série, via used_seed). XOR avaliado sobre a base DFS.

## Resultados por estratégia

| Estratégia | length | cobertura% | compat% | valid_paths | diversidade | t(ms) |
|---|---:|---:|---:|---:|---:|---:|
| A* shortest | 40±2 | 37.0±9% | 26.9±6% | 25±6 | 0.344 | 3 |
| Random DFS | 54±8 | 44.0±10% | 30.7±3% | 29±4 | 0.255 | 4 |
| Coverage-greedy | 54±13 | 44.6±12% | 33.7±8% | 32±12 | 0.293 | 9 |
| Centrality A* | 40±2 | 38.0±9% | 27.0±7% | 25±6 | 0.343 | 4 |

- Correlação cobertura→compat%: **+0.694**
- Correlação cobertura→valid_paths: **+0.719**
- Correlação length→cobertura: +0.515
- Hipótese (cobertura-máxima >> A* em compat%): **REFUTADA** (melhor=Coverage-greedy 33.7% vs A* 26.9%)

## STEP 4 — limite teórico e eficiência

| Estratégia | upper_bound (cobertura) | eficiência (compat/cobertura) |
|---|---:|---:|
| A* shortest | 37.0% | 73.6% |
| Random DFS | 44.0% | 72.9% |
| Coverage-greedy | 44.6% | 76.5% |
| Centrality A* | 38.0% | 72.1% |

## STEP 5 — caracterização geométrica

| Estratégia | mean betweenness no caminho | length |
|---|---:|---:|
| A* shortest | 0.1349 | 40 |
| Random DFS | 0.1078 | 54 |
| Coverage-greedy | 0.1290 | 54 |
| Centrality A* | 0.1562 | 40 |

## STEP 6 — VEREDITO FINAL

1. **Gargalo overlap=0 resolvível via seleção de caminho?** PARCIALMENTE — a cobertura correlaciona com compat (+0.694), mas o ganho da melhor estratégia sobre o A* é modesto (<1.5×).
2. **Compat máxima alcançável (20×20, 30% obstáculos):** upper_bound médio (cobertura da melhor estratégia) = 44.6%, eficiência 76.5% → compat real da melhor = 33.7%. Mesmo com seleção ótima de caminho, a compat é limitada pela cobertura (nem todo ciclo cabe num único caminho) E pela fração de cruzamentos contíguos entre os tocados.
3. **XOR competitivo com A* repetido para enumeração diversa?** Ver diversidade: XOR continua produzindo variações locais; A* penalizado (série anterior ~0.48) ainda é mais diverso. XOR ganha em VAZÃO (dezenas de alternativas de um DFS único), não em diversidade.
4. **Algoritmo recomendado:** para muitas alternativas locais baratas, use um caminho-base de alta cobertura (greedy ou A* central) e faça XOR com a base de ciclos; para poucas rotas globalmente distintas, use A* repetido / Yen k-shortest.

## Conexão com a teoria

- **"o gargalo é o par (caminho,grafo), não a base":** SUPORTADO — trocar a estratégia de CAMINHO move a cobertura e a compat (corr cobertura→compat +0.694), enquanto trocar a BASE (DFS↔face) não movia (26.9%↔23.0%). A alavanca é o caminho.
- **Maximizar cobertura muda valid_paths?** corr cobertura→valid_paths = +0.719; comparar valid_paths das estratégias na tabela.
- **Novo algoritmo prático:** "A* com peso de centralidade + XOR com base de faces" é uma alternativa razoável ao Yen quando se quer muitas alternativas locais; mas NÃO supera o Yen em diversidade estrutural.
