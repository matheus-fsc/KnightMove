# martingale_analysis/

Hipótese de martingale para a busca v2 do passeio do cavalo: a taxa de
descoberta de tours é constante em `n` porque o processo de busca é,
empiricamente, uma sequência de Bernoulli i.i.d. com parâmetro `μ`
independente do tamanho do tabuleiro.

> Origem: `incremental_subtour/scaling_minimal_v2.py` mostrou que
> `nós/tour ≈ 4–5` é **constante** em n ∈ {6,8,10,12}. Se a probabilidade
> de o próximo nó da busca terminar em TOUR é `≈ μ` independentemente do
> estado e de `n`, então o processo é um martingale com incremento médio
> `μ`, e `E[nós/tour] = 1/μ = O(1)` segue trivialmente.
>
> Este diretório verifica a hipótese empiricamente (T1–T4) e produz a
> base para uma prova formal via análise de G∞ (T5).

---

## 0. Definição formal do processo (T0)

### 0.1 Sequência de eventos

A busca v2 (`incremental_subtour/scaling_minimal_v2.py`) explora a
árvore de backtracking em DFS. Cada chamada de `branch()` corresponde
a UM **nó visitado** e gera UM **evento**.

Seja a sequência:

    e_1, e_2, e_3, ...

onde cada `e_t` é o desfecho do nó `t` na ordem DFS:

| outcome             | significado                                                |
|---------------------|------------------------------------------------------------|
| `TOUR`              | nó-folha em que (estado == 1).sum() == V e o grafo é conexo |
| `SUBTOUR_EARLY`     | detector incremental fechou um ciclo de tamanho < V        |
| `CONTRADICTION_DET` | tentativa de fixar uma aresta cujo vértice já tem grau 2   |
| `R2`                | propagação R2 detectou `n_um > 2` ou grau efetivo < 2      |
| `RECURSE_INTERNAL`  | nó não terminal — propagação OK, há aresta livre, recursão |

### 0.2 Variável de Bernoulli e processo de soma

Definimos a indicadora de "sucesso" (descoberta de tour) no passo `t`:

    X_t = 1   se e_t = TOUR
    X_t = 0   caso contrário

A taxa global de sucesso é:

    μ(n) = E[X_t] = (nº de TOUR) / (nº total de eventos)

Para a base `scaling_minimal_v2.json` (K=500 tours, n=10): 500 tours em
1927 nós ⇒ `μ(10) ≈ 0.26`, e `1/μ ≈ 3.85 ≈ nós/tour` — consistência
trivial pela definição.

### 0.3 Estado do processo

Em cada evento registramos o estado:

    s_t = (depth, n_components, max_degree_free, L_current, score_pressure)

| variável           | definição                                                        |
|--------------------|------------------------------------------------------------------|
| `depth`            | `n_fixadas / E` (proporção de arestas com estado != FREE)        |
| `n_components`     | nº de componentes conexas no Union-Find sobre arestas == 1       |
| `max_degree_free`  | `max_v (grau_v - n_arestas_zero_v)` para vértices ainda livres   |
| `L_current`        | nível (anel) do vértice escolhido = `min(r, c, n-1-r, n-1-c)`    |
| `score_pressure`   | `n_um(v) * 100 + (grau_v - n_zero(v) - 2)` para o vértice alvo   |

Para eventos podados ANTES de `escolher_var()` (R2, SUBTOUR/CONTRA na
sincronização do detector), não há vértice escolhido:
`L_current = -1` e `score_pressure = NaN`.

### 0.4 Hipótese de martingale

**H_MART**: condicionado em qualquer realização do estado `s_t`,

    E[X_{t+1} | s_t] ≈ μ

ou seja, o estado `s_t` **não carrega informação preditiva** sobre o
próximo sucesso além do que `μ` já captura.

Definindo o processo acumulado:

    M_t = Σ_{k=1}^{t} (X_k - μ)

Sob H_MART:

1. `E[M_t] = 0` para todo `t`
2. `Var[M_t] = t · μ (1−μ)` (linear em `t`)
3. `M_t / √t  →  N(0, μ(1−μ))` (TCL)

### 0.5 Por que H_MART implica `nós/tour = O(1)`

Se `μ(n) ≥ μ_min > 0` independente de `n`, e o processo é i.i.d.
Bernoulli(μ), então o número esperado de eventos até o próximo `TOUR`
é geométrico com média `1/μ ≤ 1/μ_min = O(1)`, e portanto

    E[nós/tour] = 1/μ(n) ≤ 1/μ_min = O(1)

O conteúdo empírico que este diretório precisa estabelecer:

(a) μ(n) ≈ constante em n ∈ {6,8,10,12,14}      (T2.1 + T3.1)
(b) X_t aproximadamente i.i.d. (sem correlação)  (T3.2)
(c) `p(b, n) = P(TOUR | s ∈ b)` ≈ `μ` para todos os bins de estado e
    todos os `n` (T2.2 + T2.4)
(d) `M_t` sem drift e variância linear (T3.3 + T3.4)

Se TODOS passarem, a hipótese está sustentada empiricamente e o
trabalho formal restante é provar `μ(n) ≥ μ_min` analiticamente
(via análise de G∞ — frequência estacionária por anel).

Se algum FALHAR, a variável de estado que quebra a hipótese aponta
onde o processo carrega memória — informação direta sobre o que a
prova formal precisaria modelar.

---

## 1. Estrutura

```
martingale_analysis/
├── README.md                ← este arquivo (T0 + resultados T2..T5)
├── logged_backtracking.py   ← T1: v2 mínimo com event_log por nó
├── martingale_analysis.py   ← T2: μ(n), p(b,n), distribuição condicional
├── convergence_test.py      ← T3: 4 testes formais (homogeneidade, i.i.d., drift, TCL)
└── data/
    ├── event_log_n{6,8,10,12,14}.json
    ├── p_distribution.json
    ├── convergence_tests.json
    └── plots/
        ├── p_distribution_by_n.png
        ├── martingale_process.png
        └── convergence_evidence.png
```

## 2. Resultados

### 2.1 T1 — taxa μ(n) por tabuleiro

K=500 tours, timeout 600s, todos os logs em `data/event_log_n{n}.json`:

| n  | eventos | tours | **μ(n)** | nós/tour | t (s) |
|----|---------|-------|----------|----------|-------|
| 6  | 2673    | 500   | **0.1871** | 5.35   | 3.51  |
| 8  | 2214    | 500   | **0.2258** | 4.43   | 4.07  |
| 10 | 1927    | 500   | **0.2595** | 3.85   | 5.97  |
| 12 | 2124    | 500   | **0.2354** | 4.25   | 9.40  |
| 14 | 2036    | 500   | **0.2456** | 4.07   | 12.23 |

Spread global = 27.9%; para n ≥ 8, spread = 13.0%.
Forma forte de H_MART (μ constante) **REFUTADA**.

### 2.2 T2 — análise condicional dupla

**(A) p(b, n) por variável de estado** (target = X_{t+1} = é o próximo evento TOUR?):
- depth: bins iniciais com p≈0 (tour só aparece em depth>0.95)
- L_parent: variação substancial entre n; L=2 varia de 0.20 (n=6) a 0.29 (n=10)
- n_components: nc=1 sempre dispara TOUR; nc>1 → p≈0

**(B) μ por região do PAI** (L_parent ≥ 2 vs < 2):
```
n=6: μ_int=0.158  μ_bor=0.229   (n=6 INVERTE — borda > interior!)
n=8: μ_int=0.253  μ_bor=0.154
n=10: μ_int=0.266  μ_bor=0.178
n=12: μ_int=0.252  μ_bor=0.147
n=14: μ_int=0.246  μ_bor=0.246
```
- spread μ_interior = 40.5% (PIOR que μ_total)
- n=6 mostra **inversão qualitativa** — efeito de borda não explica divergência

**(C) μ padronizado** (composição n=10 como referência):
| padronização | spread |
|--------------|--------|
| depth        | 26.5%  |
| L_parent     | 52.5%  |
| n_components | 26.1%  |
| (μ_total)    | 27.9%  |

Nenhuma padronização reduz spread ⇒ **variação é genuinamente intrínseca**, não composição.

**(D) χ² de homogeneidade**: rejeita H₀ em múltiplos bins mesmo restringindo a n ≥ 8 (depth alto, L_parent=2/3, nc=2/4-5/6-10, pressure Q3).

### 2.3 T3 — testes formais de convergência

| critério | resultado | detalhe |
|----------|-----------|---------|
| **μ(n) não decai** | ✓ SIM | slope OLS = +0.0063, p=0.16 (não rejeita H₀:slope=0); R²=0.53 |
| **Var[M_t] linear em t (R²>0.9)** | ✓ SIM | R²_lin médio = 0.958 |
| **α em Var ∝ t^α ≈ 1** | ✗ NÃO | α médio = 1.38 (super-linear; α=1.5 em n=6 e n=14) |
| **μ_min ≥ 0.10 (95% conf.)** | ✓ SIM | min LB Wilson = 0.173 ⇒ E[nós/tour] ≤ 5.79 |
| **TCL martingale** | ✓ SIM | Shapiro p médio = 0.28, KS p médio = 0.57 |

`σ²_estimated / σ²_iid` por n: 7.74, 2.08, 1.18, 2.49, 4.32 — variância empírica
SUPERA o predito i.i.d. ⇒ incrementos têm **dependência positiva fraca**
(eventos próximos na sequência DFS estão na mesma subárvore e tendem a ter
outcomes correlacionados).

## 3. Conclusão

### 3.1 Veredito sobre H_MART

**H_MART forte (X_t i.i.d., μ constante em n) — REFUTADA**:
- spread μ(n) = 28% globalmente
- χ² rejeita homogeneidade em múltiplos bins
- α=1.38 indica correlação positiva nos incrementos
- μ_interior(L_p≥2) tem spread ainda maior (40%)

**Versão fraca SUSTENTADA** (suficiente para O(1)):
- μ(n) não decai com n (slope p=0.16, ligeiramente positivo)
- μ_min ≥ 0.173 com 95% confiança nos n testados
- E[nós/tour] ≤ 1/μ_min ≤ 5.79 (limite empírico)
- Var[M_t] permanece linear (R²=0.96), com TCL aproximado

### 3.2 Teorema-candidato e lacuna formal

Empiricamente sustentável:

> **Conjectura (O(1) por p_min):** Para o backtracking v2 mínimo no
> grafo do cavalo n×n com R2 + detector incremental, existe `μ_min > 0`
> independente de `n` tal que `μ(n) ≥ μ_min` para todo `n ≥ 6`.
> Como consequência, `E[nós/tour] ≤ 1/μ_min = O(1)`.

Lacuna que falta provar analiticamente:

> **Lema (a provar):** `lim inf_{n→∞} μ(n) > 0`.

Os dados sugerem `μ(n) → ` algo entre 0.23 e 0.27 quando n cresce
(slope ligeiramente positivo, não decay). Não há evidência de
`μ(n) → 0`.

### 3.3 Onde a estrutura está concentrada

O χ² indicou que `L_parent` (nível do vértice escolhido pelo pai)
carrega a maior variação preditiva entre n. Em particular:
- **n=6 inverte** o padrão interior > borda (borda mais frutífera)
- **n≥8** tem padrão consistente: interior gera mais tours por evento

Interpretação: em n=6 (V=36), o grafo é pequeno demais para o anel
interno ser "produtivo"; em n≥8 a propagação R2 a partir do interior
encaminha rapidamente para um tour. Para a prova formal: o limite
inferior `μ_min` deve ser derivado a partir do regime interior (que
domina assintoticamente: `n_ev_int / n_total → 1` quando n cresce).

### 3.4 Conexão com a heurística local f∞

A análise de `local_phase_heuristic/` mostrou que a frequência
estacionária por anel `f∞(L)` é fixa e independente de n. Se o
autovetor dominante da matriz de transferência da busca coincide
com `f∞(L)`, então o limite assintótico de `μ(n)` é uma média
ponderada por anel:

    μ_∞ = Σ_L f∞(L) · P(TOUR | L_parent = L, n=∞)

Os números desta análise dão estimativas para P(TOUR | L_parent = L):
- L=0: ≈ 0.10–0.18
- L=1: ≈ 0.04–0.30 (alta variância em n=6)
- L=2: ≈ 0.20–0.29
- L=3: ≈ 0.36 (n=8), 0.36 (n=10), 0.40 (n=12), 0.40 (n=14)

A ESTABILIDADE em L grande (3, 4) é a base candidata para a prova
formal: se `P(TOUR | L)` é constante para L ≥ 2 quando n → ∞, então
`μ_min` é determinado pelos pesos `f∞(L)`.

## 4. Plots gerados

- `data/plots/p_distribution_by_n.png` — grid p(b) × n com IC Wilson
  e marcadores `*` onde χ² rejeita homogeneidade.
- `data/plots/mu_decomposition.png` — μ_total, μ_interior/μ_borda,
  μ_padronizado lado a lado.
- `data/plots/chi2_heatmap.png` — -log₁₀(p) por bin para todos n
  e só n≥8.
- `data/plots/convergence_evidence.png` — resumo dos 4 testes.
- `data/plots/martingale_process.png` — X_t, M_t, Var[M_t]×t e
  QQ-plot para n=10.

## 5. Como reproduzir

```bash
cd martingale_analysis/
python3 logged_backtracking.py --ns 6,8,10,12,14 --alvo 500
python3 martingale_analysis.py
python3 convergence_test.py
```

Tempo total: ~40s em CPU single-thread (Python 3.14, numpy 2.x).
