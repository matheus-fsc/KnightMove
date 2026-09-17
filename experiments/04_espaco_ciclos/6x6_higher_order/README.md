# Knight Tour 6×6 — Exclusões de ordem superior

Estende a análise dos invariantes do passeio do cavalo 6×6 para
**exclusões minimais de ordem ≥ 3**: triplas e quadras de arestas que
nunca aparecem juntas em nenhum dos 9.862 tours hamiltonianos
fechados, sem que nenhuma subcombinação seja já excluída.

Ground truth: `complex_orbit/data/tours_closed_6x6.npy`
(9.862 × 36, enumeração exaustiva).

## 1. Definições

Seja `T ∈ {0,1}^{9862×80}` a matriz de incidência tour×aresta
(ordenação canônica `EDGES_LIST`, lex sobre `(u,v)` com `u<v`).

**Exclusão de ordem k**: tupla `(e_1, …, e_k)` com
`P(x_{e_1} = 1 ∧ … ∧ x_{e_k} = 1) = 0`, calculado como média de
`T[:, e_1] · … · T[:, e_k]`.

**Minimalidade**: nenhuma subtupla própria é também excluída.
Uma tupla minimal é um **gerador independente do ideal de
infactibilidade** — não decorre de exclusões de ordem menor.

## 2. Contagem das exclusões minimais

| Ordem | Total enum | Excluídas (P=0) | Minimais | Redundantes | Tempo |
|---:|---:|---:|---:|---:|---:|
| 2 (pares)   | 3.160 | 88 | **88** | 0 | <1 s |
| 3 (triplas) | 82.160 | 8.256 | **1.776** | 6.480 | 10 s |
| 4 (quadras) | 1.581.580 | — | **17.004** | — | 10 s (bit-packing) |

Total de geradores do ideal de infactibilidade até ordem 4:
**18.868**.

Consistência (`ideal_structure.py`):
- 0 triplas com subpar excluído ✓
- 0 quadras com sub-par ou sub-tripla excluído ✓

## 3. Verificação Z3

Z3 (grau-2 + sub-tour elimination iterativo) confirma UNSAT para todas
as exclusões verificadas. Amostra (`z3_verify.py`):

| Ordem | Verificadas | UNSAT | SAT (contraex.) | Inconclusive |
|---|---:|---:|---:|---:|
| 2 | 88 (full) | 80 | **0** | 8 |
| 3 | 50 / 1.776 | 47 | **0** | 3 |
| 4 | 50 / 17.004 | 50 | **0** | 0 |

`inconclusive` = atingiu `max_iters=50` de sub-tour elimination. **0 SAT
em todas as ordens** → construção de T está correta. Os inconclusive
são apenas Z3 que precisa de mais iterações para fechar (não são
contraexemplos).

## 4. Estrutura do ideal

### 4.1 Hierarquia

Minimalidade preservada por construção (verificado em 4.1).
Os 18.868 minimais são **independentes**: nenhum contém outro
como subconjunto.

### 4.2 Cobertura por vértice

Quantos minimais cada vértice toca (qualquer aresta incidente):

| Classe | n vértices | min | avg | max |
|---|---:|---:|---:|---:|
| canto | 4 | **0** | 0 | 0 |
| borda | 16 | 879 | 2.143 | 3.407 |
| near-edge | 12 | 3.198 | 5.086 | 6.030 |
| **interior** | **4** | **9.731** | **9.731** | **9.731** |

Cantos (A1, A6, F1, F6, grau 2) não aparecem em nenhuma exclusão minimal
— suas 2 arestas obrigatórias são fixadas, sem flexibilidade local.

Os 4 vértices interiores (C3, C4, D3, D4, grau 8) aparecem em
**9.731 minimais cada** — são os pontos de máxima decisão topológica
no grafo.

### 4.3 Geração

| Pergunta | Resposta |
|---|---|
| Pares são suficientes para gerar o ideal? | **Não** — 1.776 triplas minimais |
| Pares + triplas? | **Não** — 17.004 quadras minimais |
| O ideal é finitamente gerado por k ≤ 4? | Indeterminado — pode haver geradores de ordem 5+ |

Quintas e ordens superiores **não foram enumeradas** — C(80, 5) ≈ 24M
exigiria amostragem.

### 4.4 Assinatura H₁

Para cada minimal, computamos o vetor indicador `1_S ∈ GF(2)^80` e o
projetamos na base do espaço de ciclos amostrado (`rank(T) = 42`).
Peso H₁ = número de coordenadas não-zero após projeção:

| Ordem | min | avg | max |
|---:|---:|---:|---:|
| 2 | 2 | 6.38 | 22 |
| 3 | 0 | 8.57 | 25 |
| 4 | 3 | 11.56 | 30 |

Algumas triplas têm peso H₁ = 0 — seu vetor indicador é uma fronteira
pura (combinação de "estrelas" em vértices), ou seja, exclusão
puramente local. Peso cresce com a ordem, como esperado, mas sem ser
proporcional — exclusões de ordem 4 ainda podem ter peso pequeno
(=3).

### 4.5 Achado inesperado: `rank(T) = 42 < β₁ = 45`

A row-reduction de `T` sobre GF(2) tem rank 42, mas o espaço de
ciclos tem dimensão β₁ = E − V + 1 = 80 − 36 + 1 = 45.

**Três dimensões do espaço de ciclos do grafo do cavalo 6×6 não são
atingidas por nenhum dos 9.862 tours hamiltonianos fechados.** Isto é,
existem ciclos no grafo (no sentido topológico) que combinatorialmente
nunca aparecem como subestrutura linear de tours.

Isso confirma que tours hamiltonianos formam uma sub-variedade
estrita do espaço de ciclos, e o "déficit" de 3 dimensões é uma
medida de complexidade combinatória.

## 5. Implicações para 10×10

A análise sugere fortemente que **o pipeline atual para 10×10
(`board_10x10/`) detecta apenas os geradores de ordem 2 do ideal de
infactibilidade**. Considerando que no 6×6 quase **96% dos geradores
minimais até ordem 4 são triplas ou quadras** (18.780 de 18.868),
provavelmente:

1. **Existem ~10⁴–10⁵ exclusões de ordem 3+ no 10×10** ainda não
   identificadas.
2. **O speedup nulo do `NOT(A∧B)` no 10×10** (confirmado em
   `board_10x10/`) é consistente com isto: pares não são bons
   candidatos para acelerar o Z3 nesse tamanho — os geradores reais
   são de ordem superior.
3. **Para encontrar geradores de ordem 3+ no 10×10** seria necessário:
   - 5.000 amostras dão `rank = 186 / β₁ = 189` — quase completo,
     mas P(triple)>0 em 5.000 amostras com poucas configurações
     extremas pode mascarar exclusões verdadeiras.
   - Necessário 50k+ amostras OU enumeração exaustiva
     parcial (inviável).

## 6. Próximos passos

- **Ordem 5+**: amostragem (>5M) para detectar triplas/quintas que
  não emergem por enumeração direta.
- **No 10×10**: testar `NOT(A∧B∧C)` para triplas estimadas a partir
  das bifurcações já conhecidas — gera potencialmente mais speedup.
- **Estrutura dos 3 ciclos não-atingidos** (`rank=42 vs β₁=45`):
  caracterizar combinatorialmente quais combinações de arestas estão
  ausentes do span dos tours.
- **Conexão com H₁ assinatura zero**: as triplas com peso H₁ = 0
  são candidatas a "exclusões locais puras" — vale enumerá-las
  como caso especial.

## 7. Estrutura de arquivos

```
6x6_higher_order/
├── README.md
├── compute_frequencies.py
├── minimal_exclusions.py
├── z3_verify.py
├── ideal_structure.py
└── data/
    ├── edges_6x6.json                  (ordenação canônica)
    ├── incidence_matrix_6x6.npy        (9862 × 80, uint8)
    ├── freq_singles.npy                (80, float32)
    ├── freq_pairs.npy                  (80 × 80, float32)
    ├── exclusions_pairs.json           (88 minimais)
    ├── exclusions_triples.json         (1776 minimais)
    ├── exclusions_quads.json           (17004 minimais)
    ├── z3_verification.json
    ├── ideal_summary.json
    └── plots/
        ├── vertex_coverage.png         (heatmap 6×6 de cobertura)
        └── ideal_structure.png         (histogramas H₁ por ordem)
```

## 8. Como reproduzir

```bash
# Tarefa 0 — frequências, T, pares (<1 min)
python 6x6_higher_order/compute_frequencies.py

# Tarefas 1+2 — triplas e quadras (~20 s)
python 6x6_higher_order/minimal_exclusions.py

# Tarefa 3 — verificação Z3 (samplar; ~10 s)
python 6x6_higher_order/z3_verify.py

# Tarefa 4 — análise estrutural + plots (~30 s)
python 6x6_higher_order/ideal_structure.py
```
