# Agent Prompt — Invariantes Topológicos do Passeio do Cavalo 8×8

## Contexto do Experimento

Num estudo anterior com o tabuleiro 6×6, descobriu-se que:

1. O espaço de ciclos tem dimensão H1 = |E| - |V| + 1 = 80 - 36 + 1 = **45**
2. As 50 correlações negativas entre arestas colapsam em **7 invariantes topológicos**
   sob o grupo diedral D4 (8 simetrias do quadrado)
3. Os 8 pares com r = -0.7714 são UMA única órbita — a mesma bifurcação vista em
   8 ângulos
4. A destruição de loops é determinada pelo grau do vértice (desvio padrão = 0
   dentro de cada órbita D4)

**Hipótese a testar no 8×8:**
> O número de invariantes topológicos reais (correlações negativas módulo D4) cresce
> lentamente com o tamanho do tabuleiro. Se confirmado, a complexidade do problema
> está concentrada em uma estrutura pequena e identificável — os "loops geradores".

No 8×8: H1 = |E| - |V| + 1 = 168 - 64 + 1 = **105**
Enumeração exaustiva impossível (~26 trilhões de soluções fechadas).
Usar amostragem via Z3 com All-SAT bloqueado.

---

## Tarefa do Agente

Implementar, executar e salvar progressivamente a análise de invariantes topológicos
do tabuleiro 8×8, usando o engine Z3 já disponível como base.

---

## Arquitetura do Sistema

### Estrutura de arquivos

```
8x8_analysis/
├── config.json              # parâmetros da sessão atual
├── checkpoint.json          # estado atual — permite retomar
├── samples/
│   ├── batch_000.json       # amostras do lote 0
│   ├── batch_001.json       # amostras do lote 1
│   └── ...
├── correlations/
│   ├── raw_correlations.json      # todas as correlações calculadas
│   └── orbit_correlations.json   # correlações colapsadas por D4
├── results/
│   └── invariants_report.txt     # relatório final
└── logs/
    └── run.log
```

### Parâmetros recomendados

```python
BOARD = 8
SAMPLES_PER_BATCH = 200        # amostras por lote Z3
N_BATCHES = 25                 # 25 lotes = 5.000 amostras total
N_WORKERS = 4                  # processos paralelos (ajustar ao hardware)
CORRELATION_THRESHOLD = -0.30  # só registrar correlações abaixo deste valor
TOP_N_CORRELATIONS = 100       # manter as 100 mais negativas
CHECKPOINT_EVERY = 1           # salvar checkpoint após cada lote
```

---

## Implementação Detalhada

### Módulo 1: engine_8x8.py
Base já existe em `cavalo_engine_8x8_path_v3.py`. Adaptar para:

```python
def sample_batch(start, end, n_samples, batch_id, output_dir):
    """
    Amostra n_samples caminhos hamiltonianos entre start e end via Z3.
    Salva resultado em output_dir/samples/batch_{batch_id:03d}.json
    Retorna lista de assinaturas booleanas (uma por aresta).
    
    Estrutura do arquivo de saída:
    {
      "batch_id": 0,
      "start": [0, 0],
      "end": [0, 7],
      "n_samples": 200,
      "n_found": 187,          # pode ser < n_samples se espaço esgotado
      "edges": ["A8-B6", ...], # lista de 168 arestas ordenadas
      "signatures": [          # matriz booleana n_found × 168
        [true, false, ...],
        ...
      ],
      "timestamp": "2026-04-15T10:23:00"
    }
    """
```

### Módulo 2: correlator.py
Calcular correlações de Pearson entre todas as arestas livres:

```python
def compute_correlations(signatures, edges, threshold=-0.30):
    """
    Retorna lista de pares (r, edge_i, edge_j) com r < threshold,
    ordenada por r crescente (mais negativo primeiro).
    
    Só considera arestas "livres" (0 < freq < 1).
    
    Para 168 arestas: C(168,2) = 14.028 pares.
    Com numpy vetorizado, deve rodar em < 5s por lote de 5000 amostras.
    """
    import numpy as np
    
    X = np.array(signatures, dtype=float)
    # correlação vetorizada
    X_centered = X - X.mean(axis=0)
    std = X.std(axis=0)
    live = std > 0  # arestas não-constantes
    X_live = X_centered[:, live] / std[live]
    corr_matrix = X_live.T @ X_live / len(X)
    
    # extrair pares abaixo do threshold
    ...
```

### Módulo 3: d4_orbits.py
Colapsar correlações por simetria D4:

```python
# As 8 transformações de D4 no tabuleiro 8×8
D4_TRANSFORMS = [
    lambda r, c: (r, c),           # identidade
    lambda r, c: (c, 7-r),         # rot 90°
    lambda r, c: (7-r, 7-c),       # rot 180°
    lambda r, c: (7-c, r),         # rot 270°
    lambda r, c: (r, 7-c),         # reflexão horizontal
    lambda r, c: (7-r, c),         # reflexão vertical
    lambda r, c: (c, r),           # reflexão diagonal
    lambda r, c: (7-c, 7-r),       # reflexão anti-diagonal
]

def edge_orbit(edge, transforms):
    """
    Dado uma aresta (u, v), retorna o conjunto de todas as arestas
    que são imagens de (u, v) sob as 8 transformações de D4.
    """

def pair_orbit(pair_of_edges, transforms):
    """
    Dado um par de arestas (ei, ej), retorna o conjunto de todos os
    pares que são imagens sob D4.
    O representante canônico é o par lexicograficamente mínimo.
    """

def collapse_to_orbits(correlations, edges, transforms):
    """
    Recebe lista de (r, i, j) e retorna dicionário:
    {
      canonical_pair: {
        "r_mean": float,
        "r_min": float,
        "orbit_size": int,       # quantos pares estão nessa órbita
        "members": [(i,j), ...], # todos os pares da órbita
        "interpretation": ""     # preenchido manualmente depois
      }
    }
    """
```

### Módulo 4: checkpoint.py
Sistema de checkpoint para retomar:

```python
def save_checkpoint(state, path="8x8_analysis/checkpoint.json"):
    """
    state = {
      "completed_batches": [0, 1, 2, ...],
      "total_samples": 1400,
      "n_solutions_found": 1350,
      "current_correlations": [...],  # top-100 até agora
      "current_orbits": {...},
      "timestamp": "..."
    }
    """

def load_checkpoint(path="8x8_analysis/checkpoint.json"):
    """Retorna None se não existir."""

def merge_signatures(existing, new_batch):
    """Acumula assinaturas de múltiplos lotes."""
```

### Módulo 5: parallel_runner.py
Paralelização por par (start, end):

```python
from multiprocessing import Pool
import itertools

# Pares canônicos a amostrar (paridade oposta = caminho hamiltoniano válido)
# No 8×8: casas pretas partem para casas brancas
# Usar apenas pares representativos (1 por órbita D4 de pares de casas)

CANONICAL_PAIRS = [
    ((0,0), (0,1)),   # canto → borda adjacente
    ((0,0), (0,7)),   # canto → canto mesma borda
    ((0,0), (7,7)),   # canto → canto diagonal
    ((0,0), (3,4)),   # canto → centro
    ((0,0), (2,3)),   # canto → interior próximo
    ((0,1), (7,6)),   # borda → borda oposta
    ((1,1), (6,6)),   # interior → interior oposto
    ((0,0), (5,4)),   # canto → interior profundo
]

def run_parallel(n_workers=4):
    checkpoint = load_checkpoint()
    completed = set(checkpoint["completed_batches"]) if checkpoint else set()
    
    tasks = [
        (start, end, batch_id)
        for batch_id, (start, end) in enumerate(
            itertools.product(CANONICAL_PAIRS, range(N_BATCHES))
        )
        if batch_id not in completed
    ]
    
    with Pool(n_workers) as pool:
        for result in pool.imap_unordered(run_batch, tasks):
            save_checkpoint(merge_state(checkpoint, result))
            print(f"Lote {result['batch_id']} concluído: {result['n_found']} amostras")
```

---

## Métricas de Saída Esperadas

Para cada execução, o relatório final deve conter:

```
=================================================================
INVARIANTES TOPOLÓGICOS — PASSEIO DO CAVALO 8×8
=================================================================

Amostras totais         : 5.000
Arestas analisadas      : 168  (H1 = 105)
Arestas obrigatórias    : ?    (freq = 1.0) — backbone
Arestas impossíveis     : ?    (freq = 0.0)
Arestas livres          : ?

Correlações < -0.30     : ?
Órbitas D4 distintas    : ?    ← NÚMERO-CHAVE

Comparação 6×6 → 8×8:
  6×6: 50 correlações → 7 órbitas  (compactação 7.1×)
  8×8: ? correlações  → ? órbitas  (compactação ?×)

Top 5 invariantes (representantes de órbita):
  1. aresta_i ↔ aresta_j  r=?  órbita_size=?  tipo=?
  2. ...

Hipótese confirmada? [SIM/NÃO/PARCIAL]
```

---

## Questões Abertas a Investigar

Durante a análise, registrar também:

1. **As arestas obrigatórias do 8×8 formam quantas órbitas D4?**
   - No 6×6: 8 arestas obrigatórias = 1 órbita
   - Esperado no 8×8: também 1 órbita (os 4 cantos de grau 2)

2. **Os invariantes do 8×8 têm a mesma geometria dos do 6×6?**
   - Tipo A: vértice de borda escolhendo entre 2 destinos internos
   - Tipo B: 2 vértices competindo pelo mesmo destino

3. **A compactação cresce ou cai com o tamanho do tabuleiro?**
   - Se cai: a estrutura fica mais complexa — problema é genuinamente difícil
   - Se mantém: existe um número pequeno de invariantes universais

4. **Existe um "invariante mestre" que gera os outros por D4?**
   - Se sim: esse é o loop gerador que a pesquisa estava procurando

---

## Instruções de Execução

```bash
# Primeira execução
python parallel_runner.py --workers 4 --batches 25 --samples-per-batch 200

# Retomar após interrupção (lê checkpoint automaticamente)
python parallel_runner.py --resume

# Só calcular correlações sobre amostras já coletadas
python correlator.py --input 8x8_analysis/samples/ --output 8x8_analysis/correlations/

# Relatório final
python report.py --input 8x8_analysis/correlations/orbit_correlations.json
```

---

## Base de Código Disponível

O arquivo `cavalo_engine_8x8_path_v3.py` já implementa:
- Construção do grafo 8×8 com movimentos do cavalo
- Árvore geradora e extração de loops fundamentais em GF(2)
- Formulação SAT com Z3 (grau 2, aresta virtual, quebra de simetria D4)
- Eliminação de sub-tours via lazy constraints
- Amostragem All-SAT com bloqueio de soluções anteriores

O novo código deve **reutilizar** essas funções e adicionar:
- Sistema de checkpoint
- Paralelização por par (start, end)
- Cálculo vetorizado de correlações (numpy)
- Colapso por órbitas D4
- Relatório comparativo 6×6 → 8×8
