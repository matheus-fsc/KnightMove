# 10×10: Higher-order exclusões e benchmark NOT(A∧B∧C)

Teste da hipótese motivada pelo achado do 6×6 (96% dos geradores até
ordem 4 são triplas/quadras): se `NOT(A∧B)` não acelera o 10×10,
talvez `NOT(A∧B∧C)` recupere o speedup.

**Resultado:** hipótese **refutada**. Nenhuma tripla genuína foi
encontrada; todas as 118 verificadas são estruturais (já implicadas
pelo grau-2). Speedup com `NOT(A∧B∧C)` é negativo (0.70×).

## 1. Pipeline

### Tarefa 0: Reconstrução de T
- 50 batches × 100 = 5.000 amostras carregadas de `board_10x10/data/samples/`
- T ∈ {0,1}^{5000×288}, dtype uint8
- Cada linha soma 100 arestas (Hamiltoniano fechado) ✓
- `rank(T) GF(2) = 186` (confirmando achado anterior, ≤ β₁ = 189)
- freq média 0.347 = 100/288 (o brief mencionava 0.694 mas é 2×100/288 contando duplas)

### Tarefa 1: Busca dirigida de triplas

Estratégia híbrida (não força bruta dos 3,9M):

| Método | Candidatas geradas |
|---|---:|
| 1.1 pair-extension (a partir dos 17 pares r<-0.5) | 2.132 |
| 1.2 low-degree vertex (grau ≤ 4) | 95 |
| **Total deduplicado** | **2.227** |

Distribuição de P_triple nas candidatas (com filtro `freq_c > 0.05`,
sub-pares vivos `> 0.02`):

| P_triple | contagem |
|---|---:|
| = 0 (estrito) | 118 |
| (0, 0.001) | 5 |
| [0.001, 0.005) | 222 |
| [0.005, 0.01) | 656 |
| [0.01, 0.02) | 1.226 |

### Tarefa 2: Verificação Z3

Verificou todas as 118 com P=0 + amostra de 50 com P>0:

| Veredicto | Total |
|---|---:|
| **proven** (UNSAT confirmado) | **118** |
| false_pos (SAT, tour existe) | 46 |
| unverified (timeout/max-iters) | 4 |

Tempo total: 24s (0.14 s/cand). Todas as 118 com P=0 nas amostras se
confirmaram como exclusões formais.

### Tarefa 3: Benchmark (K=200)

| Config | t_total | t_first | sol/s | speedup |
|---|---:|---:|---:|---:|
| **A** Z3 puro | 15.07 s | 0.853 s | 13.28 | 1.00× |
| **B** + 8 mandatory | 13.60 s | 0.049 s | 14.71 | **1.11×** |
| **C** + mandatory + 8 NOT(A∧B) | 15.24 s | 0.023 s | 13.13 | 0.99× |
| **D** + mandatory + 118 NOT(A∧B∧C) | **21.36 s** | 0.381 s | 9.36 | **0.71×** |
| **E** combined (B + C + D) | 17.71 s | 0.025 s | 11.29 | 0.85× |

Plot: `data/plots/benchmark_comparison.png`.

## 2. Diagnóstico do resultado negativo

Categorização das 118 triplas proven (`analysis.py`):

| Categoria | Definição | Total |
|---|---|---:|
| **estrutural** | 3 arestas compartilham um vértice | **118 (100%)** |
| genuína | sem vértice comum (exclusão topológica real) | **0** |

Distribuição estrutural por grau do vértice pivô:

| Grau do pivô | Vértices no grafo | C(grau,3) | Triplas |
|---:|---:|---:|---:|
| 3 | 8 | 1 | 8 |
| 4 | 28 | 4 | 110 |
| 5+ | 64 | 10-56 | 0 (filtro freq>0.05) |

Esperado teórico de estruturais (somando todos os vértices grau≥3 com
`C(grau, 3)`): **2.616**. Encontradas **118**. As ~2.500 perdidas são
candidatas onde alguma das 3 arestas tem `freq < 0.05` (raras) e
foram cortadas pelo filtro.

**Conclusão essencial:**

*Todas as triplas com P=0 no 10×10 derivam diretamente da restrição
de grau-2.* Para qualquer vértice de grau ≥ 3, escolher 3 arestas
incidentes força esse vértice a ter ≥ 3 arestas no tour, violando o
grau-2 que exige exatamente 2. Z3 já bloqueia isso implicitamente via
PbEq. **Adicionar `NOT(A∧B∧C)` explicitamente apenas sobrecarrega o
solver com cláusulas redundantes** sem ganho funcional.

## 3. Por que o 6×6 funciona e o 10×10 não

| | 6×6 | 10×10 |
|---|---:|---:|
| Tours hamiltonianos fechados | 9.862 | ~10¹⁵ (estimado) |
| Triplas minimais | 1.776 | ?? (≥0 genuínas, ≥118 estruturais) |
| Triplas / C(E,3) | 2.16% | 0.003% (só estruturais) |
| Speedup NOT(A∧B) | 3.8× | ~1.0× |
| Speedup NOT(A∧B∧C) | n/a (não testado) | 0.71× (slowdown) |

No 6×6 a estrutura combinatória pequena cria correlações fortes, as 8
bifurcações `r ≈ -0.77` excluem ~30% dos tours, e *cláusulas explícitas
ajudam o solver a podar agressivamente*. No 10×10 o espaço de busca é
~10¹¹× maior; cláusulas extras adicionam custo de propagação por SAT
check sem aproximadamente reduzir o tamanho do espaço de busca.

**Os geradores genuínos (não estruturais) do ideal de infactibilidade
do 10×10 (se existem) não são detectáveis com 5.000 amostras
heurísticas, e mesmo se fossem detectáveis, podem não estar
codificáveis como cláusulas locais NOT(A∧…∧X).**

## 4. Implicações

1. **Mandatory edges** continuam sendo o único ganho confirmado
   no 10×10 (~10% speedup), pois eliminam variáveis (não constraints).
2. **NOT(A∧B) e NOT(A∧B∧C)** não escalam para 10×10. A complexidade do
   solver cresce mais rápido que o ganho de poda.
3. Caminhos alternativos:
   - **Restrições simbólicas estruturais** em vez de pontuais
     (e.g., perfis de degree-3 inteiros)
   - **Aprender warm-starts** de tours observados (mais semelhante a
     LNS / large-neighborhood search)
   - **Decomposição do tabuleiro** (e.g., 10×10 = 6×6 + bordas
     L-shaped, herdar ideal do 6×6)
   - Aceitar o slowdown e melhorar a *qualidade* das amostras (e.g.,
     reduzindo viés via MCMC) em vez de a velocidade.

## 5. Resultado negativo é resultado

Este experimento confirma uma transição qualitativa entre 6×6 e 10×10:
geradores de ordem baixa que dominam o ideal no 6×6 **simplesmente não
existem** (ou são apenas estruturais) no 10×10. A hipótese motivada
pelo 6×6 está refutada.

## 6. Arquivos

```
higher_order_10x10/
├── README.md
├── triples_search.py
├── z3_verify_triples.py
├── benchmark_triples.py
├── analysis.py
└── data/
    ├── T_10x10.npy                (5000 × 288, uint8)
    ├── freq_singles.npy           (288, float32)
    ├── freq_pairs.npy             (288 × 288, float32)
    ├── triples_candidates.json    (2.227 candidatas)
    ├── triples_verified.json      (118 proven, 46 false_pos, 4 unverified)
    ├── benchmark_triples.json     (5 configs)
    ├── analysis.json              (categorização e escalonamento)
    └── plots/
        └── benchmark_comparison.png
```

## 7. Reproduzir

```bash
# T0 + T1 (rebuilt T + busca dirigida; ~30 s)
python board_10x10/higher_order_10x10/triples_search.py

# T2 (Z3 verify; ~30 s)
python board_10x10/higher_order_10x10/z3_verify_triples.py

# T3 (benchmark 5 configs; ~1.5 min)
python board_10x10/higher_order_10x10/benchmark_triples.py

# T4 (análise + categorização)
python board_10x10/higher_order_10x10/analysis.py
```
