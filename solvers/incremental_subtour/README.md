# incremental_subtour/

Detecção incremental de sub-ciclos via Union-Find durante o backtracking
do passeio fechado do cavalo.

> Status: **T0–T6 concluídos.** No 10×10 com K=200:
> v2 alcança **5.47× speedup sobre v1** e **25.2× sobre Z3 puro**;
> razão 2-fatores/tours cai de **18.1× → 1.00×** (zero 2-fatores
> desconexos chegam à folha). Escala mantida em n ∈ {6, 8, 10, 12}.

---

## 1. Problema

O `backtracking_10x10` (v1) usa propagação R1–R6 (grau-2, pares
estritos, cláusulas XOR) para podar a árvore de busca, mas **só
detecta sub-ciclos na folha** via uma verificação BFS final.

Diagnóstico observado em v1 @ K=200, 10×10:

| métrica           | v1 (com R1..R6)     |
|-------------------|---------------------|
| 200 tours         | 7.18s, 7441 nós     |
| 2-fatores na folha| 3623                |
| razão 2-fat/tour  | **18.1×**           |
| podas R2 internas | 23                  |

Ou seja: para cada tour válido, o backtracking explorava
≈18 ramos completos até a folha que terminavam em 2-fatores
desconexos (subgrafos onde todo vértice tem grau 2 mas o grafo
parte em múltiplos ciclos). Esse era o gargalo dominante,
absorvendo ~94% do trabalho.

## 2. Solução: detector incremental por Union-Find

A cada aresta `e=(u,v)` que vira `1` durante a busca,
mantemos um Union-Find com:

- `degree[v]`  — grau atual de `v` no subgrafo de arestas `=1`
- `size[root]` — tamanho da componente conexa
- `n_deg2_of_root[root]` — quantos vértices da componente já têm grau 2

A regra de detecção é exata:

> Uma componente com `size = k < V` em que **todos os `k` vértices
> têm grau 2** é um ciclo de comprimento `k < V` — um **sub-ciclo
> definitivo**: nenhuma extensão pode evoluí-lo para um tour
> hamiltoniano.

O detector retorna um de quatro status em cada `fix(e)`:

| status         | significado                                |
|----------------|---------------------------------------------|
| `OK`           | aresta absorvida, sem ciclo fechado         |
| `SUBTOUR`      | ciclo de tamanho `<V` → podar imediatamente |
| `COMPLETE_TOUR`| ciclo de tamanho `V` → candidato a tour     |
| `CONTRADICTION`| algum vértice atingiria grau >2             |

### 2.1 Variantes de undo

Backtracking requer reverter o detector ao desfazer um branch.
Implementamos duas variantes (`subtour_detector.py`):

- **`SubtourDetectorRollback`** — pilha de operações; `rollback(cp)`
  desfaz em ordem reversa. Não usa path compression.
- **`SubtourDetectorCopy`** — `checkpoint()` retorna tupla de cópias
  dos arrays internos; `rollback(cp)` restaura.

Micro-bench (2000 nós × 60 fixes cada):

| variante  | ops/ms  |
|-----------|---------|
| Rollback  | 240     |
| Copy      | **528** |

No backtracking real `K=200` os tempos finais ficam empatados
(Copy: 1.31s, Rollback: 1.28s); ambos visitam exatamente os
mesmos 1385 nós e disparam os mesmos 413 sub-ciclos —
prova cruzada de correção.

## 3. Testes unitários (`test_detector.py`)

| teste | descrição                                              |
|-------|--------------------------------------------------------|
| T1    | 4-ciclo fechado em grafo com vértices extras → SUBTOUR |
| T2    | caminho linear nunca dispara SUBTOUR                   |
| T3    | tour completo → COMPLETE_TOUR na última aresta         |
| T4    | rollback restaura estado idêntico ao inicial           |
| T5    | 100 tours aleatórios em ciclos `C_V` (6..30 vértices)  |
| T5b   | replay de 20 tours reais 10×10 em ordem aleatória      |
| T6    | grau 3 → CONTRADICTION                                 |

Todas as 7 verificações passam em ambas as variantes.

## 4. Backtracking v2 (`backtracking_v2.py`)

Mesma estrutura do `backtracking_10x10` (propagação engine,
ordenação por pressão, heurística freq-mais-extrema), com:

- Snapshot adicional: `last_seen[E]` + `detector.checkpoint()`
- Após cada `propagar(R2,R3,R6)`, `sync_detector()` empurra ao
  detector toda aresta nova `=1` desde o último snapshot. Se
  qualquer `fix()` retorna `SUBTOUR` ou `CONTRADICTION`, podar.
- Logging por categoria: `R2`, `R3`, `R6`, `SUBTOUR_EARLY`,
  `CONTRADICTION_DET`.

Cada `SUBTOUR_EARLY` é registrado em `subtour_log.json` com
`{nó, profundidade, edge_idx, tamanho do ciclo}` para análise.

## 5. Resultados @ 10×10, K=200

### 5.1 Comparativo v1 vs v2

| métrica                  | v1     | v2 Copy | v2 Rollback |
|--------------------------|--------|---------|-------------|
| `t_total`                | 7.18s  | 1.31s   | 1.28s       |
| `t_first`                | ≈0.06s | 0.053s  | 0.051s      |
| nós explorados           | 7441   | 1385    | 1385        |
| nós/tour                 | 37.2   | 6.9     | 6.9         |
| 2-fatores na folha       | 3623   | 200     | 200         |
| **razão 2-fat/tour**     | **18.1×** | **1.00×** | **1.00×** |
| `SUBTOUR_EARLY`          | —      | 413     | 413         |
| podas R2                 | 23     | 5       | 5           |

### 5.2 Speedups vs Z3 puro

| método              | t (s)  | speedup |
|---------------------|--------|---------|
| Z3 puro             | 33.11  | 1.00×   |
| Z3 + mandatory      | 54.42  | 0.61×   |
| BT v1               | 7.18   | 4.61×   |
| BT v2 (Copy)        | 1.31   | **25.24×** |
| BT v2 (Rollback)    | 1.28   | **25.89×** |

(Z3 + mandatory PIORA por adicionar restrições redundantes — observação
consistente com `xor_clauses_benchmark`.)

### 5.3 Análise dos sub-ciclos detectados (`subtour_analysis.json`)

- **Cobertura**: 100% (zero 2-fatores desconexos na folha)
- **Tamanho dos sub-ciclos**: mediana 18 vértices, média 32, distribuição
  bimodal (cluster pequeno 5-20 + cluster grande 60-100)
- **Profundidade do disparo**: mediana `n_fixed=288/288` (i.e. quase
  todas as detecções ocorrem nos últimos passos da árvore — sub-ciclos
  fecham TARDE na maioria dos ramos)

Implicação: o ganho de v2 **não vem de cortar profundidade**, vem de
**eliminar verificações BFS na folha** (3423 a menos). Em termos
estruturais: a propagação R2 cascateia rapidamente, fixando quase
todas as 288 arestas; o detector "captura" o momento exato em que
um sub-ciclo se forma e abandona a propagação restante.

### 5.4 Plots

- `data/plots/speedup_v2.png` — barras de `t_total` e `t_first`
  para os 5 métodos (escala log)
- `data/plots/poda_breakdown.png` — categorias de evento em v1 vs v2
- `data/plots/subtour_size_hist.png` — histograma de tamanhos
- `data/plots/subtour_depth_hist.png` — histograma de profundidades

## 6. Escalonamento — n ∈ {6, 8, 10, 12} (`scaling_minimal_v2.py`)

Variante mínima: grafo do cavalo + 8 mandatórias dos cantos + R2 +
detector incremental, **sem R3/R6**. Constrói tudo do zero — não
depende de amostras pré-computadas.

K=500 tours alvo, timeout 300s:

| n  | V   | E   | β₁   | n_free | nós/tour | razão | t (s) |
|----|-----|-----|------|--------|----------|-------|-------|
| 6  | 36  | 80  | 45   | 72     | 5.3      | 1.00× | 0.76  |
| 8  | 64  | 168 | 105  | 160    | 4.4      | 1.00× | 1.20  |
| 10 | 100 | 288 | 189  | 280    | 3.9      | 1.00× | 1.65  |
| 12 | 144 | 440 | 297  | 432    | 4.2      | 1.00× | 2.79  |

Observações:

1. **`nós/tour` é praticamente constante em `n`** (≈ 4–5).
   Comparar com v1 @ 10×10 (37.2) — ganho de **~8–9×** em estrutura
   de busca apenas pela troca R3/R6 → detector.
2. **`razão = 1.00×` em todos os tamanhos** — a hipótese inicial
   ("razão cairá de 17.4× para próximo de 1×") foi conservadora.
3. Tempo cresce moderadamente (~3.7× indo de 6 → 12) — dominado por
   custo de propagação R2 (∝ V).
4. **Sem R3/R6, o detector sozinho já supera o pipeline v1 completo
   no 10×10**: 1.65s @ K=500 minimal vs 7.18s @ K=200 v1.

## 7. Cobertura do detector — por que não há resíduo

O detector dispara `SUBTOUR` no momento exato em que uma aresta fecha
um ciclo dentro de uma componente onde todos os vértices têm grau 2.

Caso teoricamente possível em que escaparia: ciclo "inevitável" mas
ainda não fechado — ex. cadeia A-B-C-D com A-B, B-C, C-D fixadas,
sendo D-A a única livre incidente a D. R2 fixa D-A em 1 → próximo
sync detector fecha o ciclo → `SUBTOUR`. R2 já força esse fechamento
na propagação seguinte.

Por isso `n_subtour_leaf = 0` consistentemente. Não é necessário
implementar look-ahead.

## 8. Lacunas honestas

- **Heurística de variável**: ainda é `argmax(n_um*100 + (grau-n_zero-2))`
  herdada da v1. Pode haver ganho marginal com heurísticas que prefiram
  arestas adjacentes a componentes pequenas (potencialmente fechando
  sub-ciclos cedo).
- **Detector é específico para grau-2 + conectividade**. Não generaliza
  para problemas com restrições mais complexas (ex. caminhos hamiltonianos
  abertos com restrições adicionais).
- **Pipeline R3/R6 ficou redundante**: as 84 cláusulas R3 e 3 XORs R6
  perdem valor quando o detector está ativo. v2 (full) preserva por
  conservadorismo; v2 minimal mostra que podem ser removidos sem perda.
- **Não testamos n > 12**. Custos da propagação R2 crescem com V e E;
  para n=14 (V=196, E=624) o tempo extrapolado seria ~4–6s para K=500.

## 9. Arquivos

```
incremental_subtour/
├── README.md                       ← este arquivo
├── subtour_detector.py             ← Union-Find Rollback e Copy
├── test_detector.py                ← 7 testes unitários (todos passam)
├── bench_detector_overhead.py      ← micro-bench Rollback vs Copy
├── backtracking_v2.py              ← v2 full (R1..R6 + detector)
├── analyse_subtours.py             ← T3: distribuição dos sub-ciclos
├── benchmark_v2.py                 ← T4: comparativo + plots
├── scaling_minimal_v2.py           ← T5: n ∈ {6,8,10,12} sem R3/R6
└── data/
    ├── v2_K50.json                 ← v2 full @ K=50
    ├── v2_K200.json                ← v2 full @ K=200
    ├── v2_K200_rollback.json       ← v2 rollback @ K=200 (sanity)
    ├── subtour_log.json            ← 413 sub-ciclos detectados
    ├── subtour_analysis.json       ← histogramas e estatísticas
    ├── benchmark_results.json      ← comparativo 5 métodos
    ├── scaling_minimal_v2.json     ← K=500 nos 4 tamanhos
    └── plots/
        ├── speedup_v2.png
        ├── poda_breakdown.png
        ├── subtour_size_hist.png
        └── subtour_depth_hist.png
```
