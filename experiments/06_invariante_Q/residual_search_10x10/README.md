# residual_search_10x10/ — propagação em cascata e busca residual no 10×10

Replica em 10×10 o experimento de [`residual_search/`](../residual_search/)
do 6×6. A pergunta a responder:

1. O **ideal local** (mandatory + pares + triplas + quadras + XORs) captura
   completamente os 2-fatores? (no 6×6 capturou: 13.422 2-fatores
   reproduzidos exatamente.)
2. A **conectividade** é o único obstáculo residual? Em qual razão?
3. **Backtracking + R1..R6** é competitivo com Z3 puro / Z3 + mandatory?

Hipótese inicial (à luz de [`project_10x10_higher_order_done`]): no 10×10
o gap entre 2-fatores e tours deve ser **muito maior** (transição
qualitativa já observada em ordens superiores).

## Setup

- E = 288 arestas no grafo do cavalo (vs 80 no 6×6); |V| = 100
- 8 obrigatórias (mesmas dos 4 cantos, agora rotuladas A10/J10/A1/J1):
  A10-C9, A10-B8, J10-H9, J10-I8, B3-A1, I3-J1, C2-A1, H2-J1
- 84 pares estritos (P_amostra=0 com ambas freq>0; modo `todos`)
- 3 cláusulas XOR duais (de `xor_10x10_clauses.json`)
- 5.000 tours amostrados para frequências e correlações
- R4 omitido: equivalente a R2 em 10×10 (cf. [`project_10x10_higher_order_done`])
- Amostragem ao invés de ground truth: **não existe enumeração exaustiva
  10×10**; toda verificação é por amostra ou Z3.

## Regras (idênticas ao 6×6)

| R  | Disparo                                                       | Conclusão                       |
|----|---------------------------------------------------------------|---------------------------------|
| R1 | freq[e] = 1.0 (resp. 0.0)                                     | x_e := 1 (resp. 0)              |
| R2 | vértice v com 2 arestas em 1                                  | demais incidentes em v := 0     |
| R2 | vértice v com (deg(v) − n_zero) = 2                           | as 2 livres restantes := 1      |
| R3 | par excluído (e₁,e₂) e x_{e₁}=1                               | x_{e₂} := 0                     |
| R6 | XOR de suporte S com apenas uma livre em S                    | valor determinado pela paridade |

(R4/R5 omitidas — vide acima.)

## Resultado da cascata por nível

| Nível | Regras                       | fixadas = 1 | fixadas = 0 | FREE | iterações |
|------:|------------------------------|------------:|------------:|-----:|----------:|
|     0 | R1                           |           8 |           0 |  280 |         2 |
|     1 | R1 + R2                      |           8 |           0 |  280 |         2 |
|     2 | + R3 (84 pares)              |           8 |           0 |  280 |         2 |
|     3 | + R6 (3 XORs)                |           8 |           0 |  280 |         2 |

**Idêntico ao 6×6**: a propagação não fixa nenhuma variável além das
8 obrigatórias dos cantos. A razão estrutural é a mesma: cada canto
tem deg=2 já saturado e não há livres incidentes; nenhum par/tripla
contém duas obrigatórias suficientemente próximas para disparar; as
XORs duais têm suporte parcialmente sobre obrigatórias mas restam
livres demais para fechar.

Comparação 6×6 vs 10×10 (após toda a cascata):

| tabuleiro | E | mandatory | FREE | fração FREE |
|-----------|---|-----------|------|-------------|
| 6×6       |  80 |  8 |  72 | 90,0 % |
| 10×10     | 288 |  8 | 280 | 97,2 % |

→ a propagação local é **estritamente menos informativa** em 10×10
em termos relativos (e em absoluto, dado que `n_free` quadruplicou).

Visualização: [`data/plots/propagation_cascade_10x10.png`](data/plots/propagation_cascade_10x10.png)

## Espaço residual

Pós-propagação (idêntico em todos os níveis):

- **n_free = 280** · n_fixadas_1 = 8 · n_fixadas_0 = 0
- Bound ingênuo: **2^280 ≈ 1,94 × 10⁸⁴**
- Tours estimados (literatura, ordem de magnitude): ~10¹⁵
- razão tours/bound ≈ 5 × 10⁻⁷⁰ → o ideal local é praticamente
  sem efeito sobre a contagem global

Estrutura espacial das 280 livres:

| classe                 | nº arestas |
|------------------------|-----------:|
| interior-interior      |        168 |
| borda-interior + reverso |       104 |
| borda-borda            |          8 |

Grafo de dependência residual:

- **1 única componente** com os 280 nós (borda e interior
  acoplados via R2 nos vértices internos)

Distribuição das frequências nas livres (n=280):

| stat | valor |
|------|------:|
| min  | 0.079 |
| max  | 0.902 |
| mean | 0.329 |

Correlações condicionais nos 5.000 tours amostrados (entre pares de livres):

| estatística                 | valor |
|-----------------------------|------:|
| \|ρ\| máx                   | 0.804 |
| \|ρ\| p99                   | 0.216 |
| \|ρ\| p90                   | 0.083 |
| \|ρ\| p50                   | 0.031 |
| média                       | 0.041 |
| fração \|ρ\| > 0.10         |  6.2 % |
| fração \|ρ\| > 0.25         |  0.7 % |

Top pares (todos com correlação fortemente negativa, da ordem
de −0.80):
`34↔35` (B10-..), `266↔283`, `2↔4`, `207↔243`, `66↔67`.
Todos refletem o trade-off local entre as duas únicas saídas de um
vértice de grau-3 ou grau-4: se a aresta a é tomada, a vizinha b
costuma não ser (e vice-versa).

→ Comparado ao 6×6: o residual é **mais esparso em correlação**
(6,2% > 0,10 vs 23,5% no 6×6) — coerente com o n_free 4× maior.
As correlações fortes seguem locais (pares de arestas no mesmo
vértice de baixo grau).

Visualização: [`data/plots/residual_structure_10x10.png`](data/plots/residual_structure_10x10.png)

## Backtracking com propagação (K=50 tours)

Estratégia (chave para fazer funcionar em 10×10):

- **Ordenação dinâmica por pressão de vértice**: a cada nó escolhe-se
  a aresta livre incidente ao vértice de maior `n_um*100 + (deg-n_zero-2)`
  — i.e., vértice mais perto de saturar grau-2.
- Entre as livres desse vértice, escolhe a com freq mais distante de 0,5
  e testa primeiro o valor mais provável.
- Esse heurístico **concentra** o branching onde a R2 cascateia,
  evitando a difusão sem propagação observada na ordenação global.

| métrica                        | valor |
|--------------------------------|------:|
| tours encontrados (alvo=50)    | **50** |
| t_total                        | **1,63 s** |
| t_first (1º tour)              | 0,11 s |
| nós explorados                 | 1.908 |
| 2-fatores válidos              | 872 |
| **fração 2-fatores conexos**   | **5,7 %** |
| nós por tour                   | 38,2 |
| nós por 2-fator                |  2,2 |
| podas R2 (viabilidade)         | 5 |
| podas R3, R6                   | 0 |

**Observação crítica**: 5,7% dos 2-fatores ach­ados são conexos
(vs 73,5% no 6×6 = 9.862 / 13.422). Em 10×10 a hipótese
"2-fatores que respeitam constraints locais ≈ tours" **falha**
drasticamente: a conectividade é o obstáculo dominante.

## Benchmark vs Z3 (K=50 tours fechados)

Todos sem hint de tour inicial; sub-tour eliminado por bloqueio
de componente (Z3) ou por filtro de folha (backtracking).

| método                           | t_first | t_total | tentativas | nós explorados |
|----------------------------------|--------:|--------:|-----------:|---------------:|
| A) Z3 puro (grau-2 + no-subtour) |  0,48 s |  3,33 s |       165 |  — |
| B) Z3 + R1 (8 mandatory)         |  0,23 s |  3,36 s |       158 |  — |
| **C) Backtracking + R1..R6**     | **0,11 s** | **1,63 s** | — | 1.908 |

Speedups vs Z3 puro:

|                   | t_first | t_total |
|-------------------|--------:|--------:|
| Z3 + mandatory    | 2,11×   | 0,99×   |
| Backtracking      | 4,27×   | **2,04×** |

- **Z3 + mandatory** acelera o 1º tour (≈2×) mas o tempo total é
  idêntico ao Z3 puro: o gargalo é a eliminação dos sub-tours,
  não a busca da 1ª solução.
- **Backtracking + R1..R6** ganha consistentemente (≈4× no 1º
  tour, ≈2× no total) porque (i) propagação é mais barata que
  PbEq + BFS; (ii) ordenação por pressão prune cedo as ramificações
  ruins; (iii) não há overhead de bloqueio incremental de soluções.

Visualização: [`data/plots/benchmark_10x10.png`](data/plots/benchmark_10x10.png)

## Resposta às hipóteses

1. **O ideal local captura os 2-fatores no 10×10?**
   Captura-os como conjunto-superset, igual ao 6×6 — o motor encontra
   872 2-fatores em 1908 nós com 5 podas. Mas a captura **não é
   restritiva**: as constraints locais conhecidas (mandatory + 84
   pares + 3 XORs) não eliminam praticamente nenhuma das ramificações,
   o que faz o backtracking degenerar para enumeração explícita do
   tipo "branch+check" no interior.

2. **Conectividade é o obstáculo residual? Em qual razão?**
   Sim. Razão 2-fatores/tours observada em 50 amostras = **17,4×**
   (vs 1,36× no 6×6). A propagação local não consegue ver
   sub-tours de comprimento ≪ 100 — coerente com a impossibilidade
   teórica de expressar conexão como combinação local fixa de
   arestas.

3. **Backtracking + R1..R6 é competitivo com Z3?**
   Sim — vence em ~2× no total e ~4× no 1º tour. A vantagem vem
   inteiramente da heurística de pressão de vértice + custo baixo
   da propagação NumPy-vetorizada, não do poder discriminativo
   das constraints (que é fraco no 10×10).

## Escalonamento projetado para 12×12

Extrapolação grosseira (apenas para orientar trabalho futuro):

- E(12×12) = 384, mandatory = 8 (cantos invariantes), n_free
  esperado ≈ 376 (≈ 97,9 % livre)
- razão 2-fatores conexos/2-fatores deve cair mais ainda
  (estimativa < 2 %), tornando o filtro de conectividade ainda
  mais caro proporcionalmente
- backtracking + pressão deve continuar funcionando, mas
  speedup vs Z3 pode encolher por aumento de gargalo de
  conectividade — caminho mais promissor é integrar
  no-subtour incremental (estilo MTZ ou eliminação de ciclos
  por DFS dentro do solver)

Conclusão para o programa de pesquisa: a hipótese
"borda fixada + interior combinatório" **falha em 10×10**
no sentido que o interior é uma única componente residual
dominada por conectividade global. Investigação subsequente
deve focar em (i) caracterização local de **sub-ciclos curtos**
proibidos (análogo às 3 órbitas proibidas do 6×6, mas escaladas)
ou (ii) uma constraint de fluxo/no-subtour incorporada à propagação.

## Arquivos

- `propagation_engine_10x10.py` — motor R1..R6 vetorizado, carrega
  amostras de `board_10x10/data/samples/`, gera
  `data/propagation_log.json`
- `residual_analysis_10x10.py` — estrutura espacial, componentes,
  correlações; gera `data/residual_variables.json`
- `backtracking_10x10.py` — DPLL com ordenação por pressão; gera
  `data/backtracking_results.json`
- `benchmark_vs_z3.py` — comparação Z3 puro / Z3+R1 / backtracking;
  gera `data/benchmark_comparison.json`
- `plots_10x10.py` — gera `data/plots/{propagation_cascade,residual_structure,benchmark}_10x10.png`
