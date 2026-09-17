# xor_clauses_benchmark — cláusulas XOR de paridade derivadas de H₁ no Z3

Testa empiricamente se as obstruções F₂-lineares descobertas em
`forbidden_cycles_6x6/` (e replicadas estruturalmente no 10×10)
aceleram o solver Z3 quando injetadas como cláusulas de paridade.

## Resumo de uma linha

**Speedup empírico ≈ 1.00× em ambos os tamanhos. A teoria é correta e
útil para caracterização, mas as cláusulas XOR não aceleram este Z3 no
nosso pipeline de amostragem.** Por que: o solver já propaga grau-2 +
sub-tour elimination rápido o bastante, e o ganho marginal das XOR é
absorvido pelo overhead de avaliar essas cláusulas extras.

## 1. Os 3 detectores duais do 6×6 (T0)

Carregados de `forbidden_cycles_6x6/data/results/dual_detectors.json`.
Para cada φᵢ ∈ (F₂⁸⁰)\* válido em todo tour (verificado: zera em
**9 862/9 862** tours conhecidos):

| φᵢ | \|supp\| | exemplo de suporte | interpretação |
|---|---|---|---|
| φ₀ | 22 | A6-C5, B6-D5, B6-A4, … | mistura cantos + anel central |
| φ₁ | **2** | F6-D5, B3-A1 | **detector mínimo**, dois pares de obrigatórias |
| φ₂ | 20 | B6-D5, B6-A4, D6-C4, … | rotação/reflexão de φ₀ |

Como cláusula Z3 (Σ x_e ≡ 0 mod 2):

```python
# suporte 2 — equivalência direta (mais rápida)
solver.add(x_a == x_b)

# suporte > 2 — soma inteira % 2
solver.add(Sum([If(x_e, 1, 0) for e in supp]) % 2 == 0)
```

## 2. Benchmark 6×6 (T1, K=500)

5 configurações comparadas via amostragem de 500 tours hamiltonianos
sem viés:

| cfg | restrições | t_total | speedup vs A | sol/s | attempts/sol |
|---|---|---|---|---|---|
| **A** | grau-2 + sub-tour elim | 6.50 s | 1.00× | 76.9 | 3.58 |
| **B** | A + 8 mandatórias | 6.63 s | 0.98× | 75.5 | 3.58 |
| **C** | B + 3 XOR (φ₀,φ₁,φ₂) | 6.91 s | 0.94× | 72.3 | 3.62 |
| **D** | B + 28 XOR (C(8,2)) | 6.73 s | 0.97× | 74.3 | 3.58 |
| **E** | B + 3 XOR + NOT(A∧B) top-8 | 7.64 s | 0.85× | 65.5 | 3.87 |

Observações:
- **A ≈ B** — adicionar as 8 obrigatórias como fatos não muda nada
  perceptível. O Z3 já infere `x_e=1` em 1 passo de propagação grau-2
  (cada canto tem 2 vizinhos, soma 2 ⇒ ambas saturadas).
- **C levemente pior que B** — as 3 cláusulas XOR (uma com soma de 22
  variáveis, outra com 20) custam mais para avaliar do que poupam.
- **D ≈ B** — usar as 28 cláusulas XOR (todas pares de obrigatórias)
  é redundante mas barato porque cada uma é `x_a == x_b` (suporte 2).
- **E pior** — adicionar `NOT(A∧B)` para top-8 pares minimal-excluded
  duplica restrição que o Z3 já vai descobrir via sub-tour cuts.

KL Bernoulli vs ground-truth: 0.45-0.60 em todas as configs — nenhuma
introduz viés mensurável.

Dados: `data/benchmark_6x6_xor.json`.

## 3. Identificação das XOR no 10×10 (T2)

Mesma análise que `forbidden_cycles_6x6/obstruction_theorem.py` faz
no 6×6, agora sobre as 5 000 amostras de `board_10x10/`:

```
V=100  E=288  β₁=189
rank(T amostrado, 5k tours) = 186   ⇒  deficit empírico = 3
rank(∂)                     = 99
8 obrigatórias (dos 4 cantos):
  A10-C9, A10-B8, J10-H9, J10-I8,
  B3-A1, I3-J1, C2-A1, H2-J1
```

C(8,2) = **28 vetores XOR** entre pares de obrigatórias.
**Todos os 28 são empiricamente válidos** nas 5 000 amostras (zera
para todo tour). Dimensões em F₂²⁸⁸:

| espaço | dim |
|---|---|
| Span(28 XORs)                    | 7  (= 8 mand − 1 redundância) |
| row_space(∂)                     | 99 |
| Span(28 XORs) + row_space(∂)     | 102 |
| **quociente (deficit detectado)**| **3** |

**Confirmação direta da conjectura formulada em
`forbidden_cycles_6x6/README.md`**: o deficit do 10×10 é **exatamente
3**, mesma dimensão que o 6×6, sem novas obstruções estruturais.

### Base mínima escolhida

3 pares linearmente independentes mod row_space(∂):

```
XOR_1 : x_{A10-C9}  ⊕  x_{J10-H9}  ≡ 0   (NW canto-superior ↔ NE)
XOR_2 : x_{A10-C9}  ⊕  x_{B3-A1}   ≡ 0   (NW ↔ SW)
XOR_3 : x_{A10-C9}  ⊕  x_{I3-J1}   ≡ 0   (NW ↔ SE)
```

Interpretação geométrica: a obstrução é a “correlação canto-a-canto”.
A10-C9 é uma das 2 obrigatórias do canto NW; cada XOR liga-a a uma
obrigatória de outro canto. Os 4 cantos estão acoplados em modo
{trivial ⊕ duplo} pelo grupo D₄ que age sobre o tabuleiro.

Dados: `data/xor_10x10_clauses.json`.

## 4. Benchmark 10×10 (T3, K=200)

| cfg | restrições | t_total | speedup vs A | sol/s | attempts/sol |
|---|---|---|---|---|---|
| **A** | grau-2 + sub-tour elim | 16.21 s | 1.00× | 12.34 | 5.82 |
| **B** | A + 8 mandatórias | 16.40 s | 0.99× | 12.19 | 5.82 |
| **C** | **B + 3 XOR mínimas** | **15.91 s** | **1.02×** | 12.57 | 5.82 |
| **D** | B + 28 XOR pares | 16.16 s | 1.00× | 12.38 | 5.82 |
| **E** | B + 3 XOR + 1 NOT(A∧B) | 22.66 s | 0.71× | 8.82 | 7.46 |

- **A ≈ B ≈ C ≈ D** — ruído dentro de 2 %. O Z3 não tira proveito
  prático das obstruções F₂-lineares no nosso pipeline. As 3 XOR
  mínimas dão uma leve melhora (1.02×) mas dentro da margem.
- **E pior** — só **1 par** com `p_coexist = 0` existe em
  `edge_pair_correlations.json` do 10×10 (o resto tem `p_coexist > 0`,
  i.e., não são exclusões estritas). Esse 1 par força UNSAT em alguns
  ramos cedo, gerando mais retries (7.46 vs 5.82 attempts/sol).

Dados: `data/benchmark_10x10_xor.json`.

## 5. Por que o speedup esperado não apareceu?

**A teoria está correta** (verificada três vezes: rank GF(2), validação
empírica, e Z3 confirmou UNSAT ao violar os detectores no 6×6 em
trabalho prévio). Mas três fatores estruturais do nosso pipeline
matam o ganho prático:

1. **Já há corte global por amostra**. O loop adiciona
   `Or([x_i != sig[i] for all i])` após cada tour encontrado para
   evitar duplicatas. Esse corte é forte e domina o tempo total
   bem mais que as 3 XOR.
2. **Sub-tour elimination iterativo é caro** mas converge rápido —
   ele já implicitamente força paridade local (cada vértice tem grau
   par no subgrafo). Nada nas XOR é mais profundo do que isso para
   o Z3 inferir via propagação.
3. **As obrigatórias dos cantos** são triviais para o Z3:
   `PbEq(2,2)` em vértices de grau 2 implica saturação em 1 passo de
   simplificação. Logo a "barreira" que a teoria descreve não existe
   para o solver — ela só é visível em formulações puramente
   simbólicas/algébricas.

Onde **as XOR teriam mais chance de ajudar**:

- problemas de **decisão** (existe um tour com aresta forçada x_e=0?)
  sem o overhead de geração K-completa
- formulações **sem sub-tour elimination iterativa** (e.g., só grau-2)
- variantes **constrangidas**, e.g., "tour com k arestas vetadas e
  3 mandatórias removidas" — aí o quociente pode realmente importar

## 6. Conclusão prática

- ✅ O resultado **teórico** vale: `deficit_{10×10} = deficit_{6×6} = 3`,
  estrutura totalmente governada pelas 8 obrigatórias dos cantos.
- ✅ As 3 XOR mínimas do 10×10 são as mais econômicas:
  cada uma é uma cláusula booleana de 2 variáveis (`x_a == x_b`).
- ❌ O **speedup prático** em amostragem Z3 com sub-tour elimination
  é ≤ 2 %. Não compensa engenharia.
- 🟡 Vale a pena re-testar em formulações alternativas (CSP puro,
  ASP, MaxSAT) onde a estrutura algébrica é menos diluída pelos
  cortes de simetria do solver.

## 7. Comparação 6×6 vs 10×10

| dimensão | 6×6 | 10×10 |
|---|---|---|
| V, E, β₁ | 36, 80, 45 | 100, 288, 189 |
| nº de obrigatórias | 8 | 8 |
| pares válidos / 28 | 28/28 | 28/28 |
| **deficit (dim quociente)** | **3** | **3** ← idêntico |
| t_total Config A | 6.5 s | 16.2 s |
| Config C speedup | 0.94× | 1.02× |
| Config E speedup | 0.85× | 0.71× |

O padrão de speedup é **qualitativamente igual** em ambos os tamanhos:
configs com XOR ficam dentro do ruído (~5 %) de A; configs com NOT(A∧B)
pioram. **A escala não muda a história.**

## 8. Estrutura de arquivos

```
xor_clauses_benchmark/
├── README.md
├── xor_6x6.py             ← T0 + T1
├── find_xor_10x10.py      ← T2 (identifica as 3 XOR via álgebra GF(2))
├── xor_10x10.py           ← T3
├── plot_speedup.py        ← T4
└── data/
    ├── benchmark_6x6_xor.json
    ├── benchmark_10x10_xor.json
    ├── xor_10x10_clauses.json
    ├── bench_6x6.log
    ├── bench_10x10.log
    └── plots/
        └── speedup_comparison.png
```

## 9. Como reproduzir

```bash
# T0 + T1 — benchmark 6×6 (~30 s)
python3 xor_clauses_benchmark/xor_6x6.py --K 500

# T2 — identifica as 3 XOR mínimas do 10×10 (~10 s)
python3 xor_clauses_benchmark/find_xor_10x10.py

# T3 — benchmark 10×10 (~2 min)
python3 xor_clauses_benchmark/xor_10x10.py --K 200

# T4 — gerar plot comparativo
python3 xor_clauses_benchmark/plot_speedup.py
```

Z3 precisa estar disponível (este projeto usa o venv em `.venv/`).
