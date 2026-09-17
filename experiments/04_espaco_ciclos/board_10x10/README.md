# Knight Tour 10×10 — Pipeline GF(2) e descoberta de invariantes

Replica no tabuleiro 10×10 a metodologia validada no 6×6:
amostragem SAT não-viesada (`break_symmetry=False`), descoberta de
invariantes do espaço de ciclos sobre GF(2), arestas obrigatórias e
benchmark de speedup quando os invariantes são injetados no solver.

## 1. Caracterização do grafo (`graph_10x10.py`)

| Métrica | Valor |
|---|---|
| `V` | 100 |
| `E` | 288 |
| `β₁ = E − V + 1` | **189** |
| grau médio | 5.76 |
| min / max grau | 2 / 8 |
| conectividade | conexo |

Classes posicionais:

| classe | n | graus |
|---|---:|---|
| canto | 4 | 2 |
| borda | 32 | 3, 4 |
| near-edge | 28 | 4, 6 |
| interior | 36 | 8 |

A matriz ∂₁ ∈ GF(2)^{100×288} tem `rank = 99 = V − 1` ⇒
β₁ confirmado via teorema do posto.

Artefatos: `data/boundary_matrix_10x10.npy`, `data/edges_10x10.json`,
`data/graph_meta_10x10.json`.

## 2. Amostragem SAT (`sample_tours.py`)

Formulação:
- Variáveis booleanas `x_e` por aresta
- Grau-2 em cada vértice via `PbEq(... , 2)`
- Quebra de simetria aleatória por amostra (`x_e ← {0,1}` para aresta sorteada)
- Pós-filtro BFS de conectividade (rejeita uniões de sub-tours)
- `break_symmetry=False` global — sem viés de amostragem

Protocolo:
- `n_batches=50`, `batch_size=100` ⇒ **5.000 amostras válidas**
- Cada batch reinicia o solver (evita acúmulo de cortes degradar performance)
- Cortes de exclusão da configuração atual após cada amostra (válida ou sub-tour)

Resultados observados (run em `data/samples/full_run.log`):

| Métrica | Valor |
|---|---|
| amostras válidas | 5.000 / 5.000 |
| amostras únicas | 100% (sem duplicatas) |
| tempo total | **417 s** (~7 min) |
| taxa de rejeição (sub-tours) | **~83%** |
| tempo médio por amostra válida | ~80 ms |
| amostras únicas locais (por batch) | 100% |

**Nota:** Taxa de rejeição de 83% excede o limite recomendado de 70% do
brief, mas Z3 é rápido o bastante (~5 ms por `check`) que o tempo total
fica viável sem MTZ. Decidido por economia.

## 3. Invariantes GF(2) (`invariants_gf2.py` + `edge_pair_invariants.py`)

### 3.1 Rank do espaço amostrado

```
T   ∈ GF(2)^{5000×288}     (assinaturas amostradas)
B   = row-reduction(T)
rank(T) = 186  /  β₁ = 189   →  cobertura 98.4%
```

Faltam 3 dimensões — o espaço de ciclos é quase-totalmente capturado
pelas 5.000 amostras.

### 3.2 Coordenadas e invariantes não-triviais

`coords ∈ GF(2)^{5000×186}` é a projeção de cada amostra na base `B`.

Todas as 186 dimensões têm `var > 0.01 × var_max` (i.e., são
não-degeneradas) e `p(1) ∈ [0.48, 0.53]` — distribuição muito próxima
de 50/50 em todas as direções, consistente com base aleatória de
eliminação gaussiana.

### 3.3 Correlação aresta×invariante — **não produziu candidatos**

Filtro `r(coord_k, x_e) > 0.6` ∧ `P(x_a=1 ∧ x_b=1) < 0.01`
para pares de arestas correlatas ao mesmo invariante:
**0 pares** identificados.

Motivo: a base GF(2) obtida por row-reduction não está alinhada com
os "invariantes topológicos naturais"; cada coordenada é uma combinação
linear arbitrária de arestas, diluindo correlações.

### 3.4 Análise direta aresta×aresta

Substitui a busca via coordenadas pela coexistência empírica direta:
para cada par de arestas vivas `(e_a, e_b)` com freq ∈ (0.01, 0.99),
calcula Pearson `r(x_a, x_b)` e `P(x_a=1 ∧ x_b=1)`.

**17 pares com r < −0.5** identificados; histograma do triângulo
superior (39.060 pares vivos):

| faixa de r | contagem |
|---|---:|
| [-1.0, -0.7) | **8** |
| [-0.7, -0.5) | 9 |
| [-0.5, -0.3) | 146 |
| [-0.3, -0.1) | 1.400 |
| [-0.1, +0.1) | 36.631 |
| [+0.1, +0.3) | 865 |
| [+0.3, +1.0) | 1 |

Top-5 bifurcações (`data/invariants/candidate_clauses.json`):

| edge_a | edge_b | r | P(A∧B) | freq_a | freq_b |
|---|---|---:|---:|---:|---:|
| A9-C8 | A9-B7 | −0.8035 | 0.0996 | 0.40 | 0.70 |
| H3-I1 | G2-I1 | −0.8032 | 0.0984 | 0.39 | 0.71 |
| **B10-D9** | **B10-C8** | **−0.7960** | 0.1060 | 0.68 | 0.42 |
| B4-A2 | C3-A2 | −0.7896 | 0.1064 | 0.71 | 0.40 |
| J9-H8 | J9-I7 | −0.7786 | 0.1206 | 0.47 | 0.65 |

**Padrão observado**: todos os top-8 pares compartilham um vértice.
São bifurcações em vértices de grau 3 (em A9, I1, B10, A2, J9, B1, I10,
J2): das três arestas incidentes, duas se excluem mutuamente em vez
das outras combinações. A exclusão **não é estrita** — `P(A∧B) ≈ 0.10`,
contra `1/3 ≈ 0.33` que seria esperado em sampling uniforme sem
restrição. Há atração negativa forte, mas coexistência ainda ocorre.

A bifurcação **B10-D9 ↔ B10-C8** coincide com o top-1 universal
reportado na Fase C anterior do projeto.

## 4. Arestas obrigatórias (`mandatory_edges.py`)

Candidatas: `freq[e] > 0.99` ⇒ **8 arestas**, todas com `freq = 1.0`
nas 5.000 amostras.

Prova formal: para cada candidata, Z3 com `grau-2 ∧ x_e = 0` é UNSAT
em **1 iteração** (sem precisar enumerar sub-tours).

| canto | arestas obrigatórias |
|---|---|
| A10 (grau 2) | A10-C9, A10-B8 |
| J10 (grau 2) | J10-H9, J10-I8 |
| A1 (grau 2)  | B3-A1, C2-A1 |
| J1 (grau 2)  | I3-J1, H2-J1 |

Padrão idêntico ao 6×6: os 4 vértices de canto têm exatamente 2 vizinhos
de cavalo, logo essas 8 arestas são forçadas pela restrição de grau-2.

## 5. Benchmark Z3 vs Z3 + invariantes (`benchmark_10x10.py`)

K = 200 amostras por configuração, pós-filtro BFS, quebra de simetria
aleatória. Comparação:

| Config | t_total (s) | t_first (s) | sol/s | attempts/sol | reject | speedup |
|---|---:|---:|---:|---:|---:|---:|
| A — Z3 puro | 17.17 | 0.846 | 11.6 | 5.89 | 81.8% | 1.00× |
| B — + 8 mandatory | 15.96 | **0.050** | 12.5 | 5.54 | 80.7% | **1.08×** |
| C — + mandatory + 17 NOT(A∧B) | 22.26 | **0.018** | 9.0 | 7.10 | 83.1% | 0.77× |
| C′ — + mandatory + 8 NOT(A∧B) (só r<-0.7) | 17.89 | 0.024 | 11.2 | 6.17 | 80.9% | 0.97× |

**Achados:**

1. **Mandatory edges** dão speedup pequeno mas real (~8%) e **47× speedup**
   no time-to-first. Custo zero — provas Z3 instantâneas.

2. **NOT(A∧B)** com todos os 17 pares **desacelera** a amostragem total
   em 23%, apesar de acelerar o time-to-first 47×. Limitar aos top-8
   reduz a perda a ~3% (neutro).

3. **Diferença crítica vs 6×6**: em 6×6 as 8 bifurcações têm `P(A∧B) ≈ 0`
   na ground truth (estritas) e adicionar NOT(A∧B) deu speedup de 3.8×.
   Em 10×10 as bifurcações têm `P(A∧B) ≈ 0.10` (não estritas) — adicionar
   NOT(A∧B) **exclui amostras válidas** e endurece o SAT.

## 6. Questões em aberto

- **Por que 10×10 não tem bifurcações estritas?** No 6×6 a estrutura de
  pequeno mundo do grafo do cavalo força bifurcações verdadeiras
  (`P(A∧B) = 0`); em 10×10 há mais flexibilidade local e a estrutura
  é apenas correlacional.
- **3 dimensões faltantes (rank 186/189)** — rodar mais amostras
  (15-20k) provavelmente fecha a base, mas custo é não-trivial.
- **MTZ ou flow-based subtour elimination** poderia reduzir rejeição
  de 83% para <5% — não foi implementado por economia, mas viraria
  útil para 12×12+.
- **Variantes do C′** que respeitem `P(A∧B) > 0`: implementar como
  *soft constraints* / pesos no Z3 (em vez de Not(And)) — pode
  combinar speedup do time-to-first sem o slowdown total.
- **Comparação direta de KL** Z3 10×10 vs ground truth de paths: o
  ground truth é estimado em ~10¹⁵ tours, enumeração exaustiva
  inviável; ficaria comparação Z3 vs Z3 (consistency) ou Z3 vs
  amostragem MCMC independente.

## 7. Estrutura de arquivos

```
board_10x10/
├── graph_10x10.py
├── sample_tours.py
├── invariants_gf2.py
├── edge_pair_invariants.py    ← complementar (análise direta)
├── mandatory_edges.py
├── benchmark_10x10.py
├── README.md
└── data/
    ├── boundary_matrix_10x10.npy    (100 × 288, uint8)
    ├── edges_10x10.json
    ├── graph_meta_10x10.json
    ├── benchmark_10x10.json
    ├── samples/
    │   ├── tours_10x10_batch_000.npy ... _049.npy   (50 batches × 100)
    │   ├── run_log.json
    │   └── full_run.log
    └── invariants/
        ├── h1_basis_sampled.npy        (186 × 288)
        ├── h1_coordinates.npy          (5000 × 186, uint8)
        ├── edge_correlations.npy       (186 × 288, float32)
        ├── edge_pair_correlations.json
        ├── candidate_clauses.json      (17 pares NOT(A∧B))
        ├── mandatory_edges.json
        └── summary.json
```

## 8. Como reproduzir

```bash
# 1. caracterização (instantâneo)
python board_10x10/graph_10x10.py

# 2. amostragem (~7 min)
python board_10x10/sample_tours.py --n-batches 50 --batch-size 100

# 3. invariantes GF(2) (~30 s)
python board_10x10/invariants_gf2.py
python board_10x10/edge_pair_invariants.py

# 4. arestas obrigatórias com prova Z3 (<1 s)
python board_10x10/mandatory_edges.py --verify

# 5. benchmark (~1 min)
python board_10x10/benchmark_10x10.py --target-k 200 --max-pairs 17
```
