# Agent Prompt — Comparação de Eficiência: Z3 vs Backtracking

## Contexto

Esta pesquisa usa dois métodos para encontrar ciclos hamiltonianos
no grafo do cavalo:

1. **Backtracking** com ordenação Warnsdorff (profundidade)
2. **Z3 SAT solver** via espaço de ciclos GF(2)

Já temos dados parciais de ambos, mas não uma comparação direta
e controlada de eficiência. Este experimento preenche essa lacuna.

---

## Dados já disponíveis (não repetir)

### Backtracking 6×6 (já computado)
- 710.064 soluções dirigidas encontradas
- ~38.5 milhões de dead-ends registrados
- Tempo total: 1.360 segundos (inclui catalogação de destruição)
- Arquivo: `destruction_catalogue.json`

### Z3 6×6 (já computado)
- 80k amostras para validação (sym=False)
- Localização: `data_6x6_sampled/`
- Tempo por batch: registrado em logs mas não sumarizado

**O que falta:** comparação direta de tempo/recurso para
encontrar K soluções com cada método, em condições controladas.

---

## Experimento a implementar

### Objetivo
Comparar Z3 e backtracking em termos de:
1. Tempo para encontrar N soluções distintas
2. Tentativas desperdiçadas (dead-ends no BT, attempts no Z3)
3. Qualidade das amostras (uniformidade)
4. Escalonamento com o tamanho do tabuleiro

### Condições controladas
- Mesmo hardware, mesmo processo, sem paralelismo
- Mesmo alvo: K = {10, 100, 1.000, 10.000} soluções
- Mesmos tabuleiros: 6×6 e 8×8
- Backtracking: com e sem ordenação Warnsdorff
- Z3: com break_symmetry=False (método validado)

---

## Implementação

### Módulo 1: benchmark_backtracking.py

```python
"""
Backtracking com medição de eficiência.
Registra: tempo, dead-ends, soluções encontradas, memória pico.
"""
import time, tracemalloc

def backtrack_timed(board_size, target_k, use_warnsdorff=True):
    """
    Roda backtracking até encontrar target_k soluções.
    
    Retorna:
    {
      "method": "backtracking_warnsdorff" ou "backtracking_puro",
      "board": board_size,
      "target_k": target_k,
      "found_k": int,
      "time_s": float,
      "dead_ends": int,
      "dead_ends_per_solution": float,
      "nodes_visited": int,           # total de vértices expandidos
      "peak_memory_mb": float,
      "solutions_per_second": float,
      "time_to_first_s": float,       # tempo para a 1ª solução
      "time_to_k_s": float,           # tempo para a K-ésima solução
    }
    """
    solutions = []
    dead_ends = [0]
    nodes_visited = [0]
    t0 = time.perf_counter()
    time_to_first = [None]
    
    # ... implementação do backtracking ...
    # Usar o código de cavalo_loop_destruicao_6x6.py como base
    # Adicionar contadores de dead_ends e nodes_visited
    
    return results
```

### Módulo 2: benchmark_z3.py

```python
"""
Z3 com medição de eficiência.
Registra: tempo, attempts, soluções encontradas, memória pico.
"""

def z3_timed(board_size, start, end, target_k):
    """
    Roda Z3 até encontrar target_k soluções para o par (start,end).
    
    Retorna:
    {
      "method": "z3_gf2",
      "board": board_size,
      "pair": (start, end),
      "target_k": target_k,
      "found_k": int,
      "time_s": float,
      "n_attempts": int,             # chamadas a solver.check()
      "failed_attempts": int,        # calls que não geraram solução válida
      "failed_per_solution": float,
      "peak_memory_mb": float,
      "solutions_per_second": float,
      "time_to_first_s": float,
      "time_to_k_s": float,
      "exhausted": bool,
      "timed_out": bool,
    }
    """
    # Usar engine_6x6_sampler.py / engine_8x8.py como base
    # Adicionar medição de tempo e memória
    
    return results
```

### Módulo 3: benchmark_runner.py

```python
"""
Roda todos os benchmarks e gera tabela comparativa.
"""

CONFIGS = [
    # (board_size, target_k)
    (6, 10),
    (6, 100),
    (6, 1_000),
    (6, 9_862),    # todas as soluções do 6×6
    (8, 10),
    (8, 100),
    (8, 1_000),
]

# Para cada config, rodar:
# 1. backtracking puro (sem Warnsdorff)
# 2. backtracking + Warnsdorff
# 3. Z3 (3 pares canônicos diferentes, média)
# 4. Z3 + constraints de borda (usando invariantes descobertos)

# Repetir 3 vezes cada experimento, reportar média ± desvio
```

---

## Métricas a reportar

### Tabela principal

```
┌─────────────┬──────┬────────┬─────────────┬─────────────────┬──────────────┐
│   Método    │  n   │  K sol │  Tempo (s)  │ Desperdício/sol │  Mem (MB)    │
├─────────────┼──────┼────────┼─────────────┼─────────────────┼──────────────┤
│ BT puro     │ 6×6  │  100   │    ?        │ dead-ends/sol   │    ?         │
│ BT Warnsd.  │ 6×6  │  100   │    ?        │ dead-ends/sol   │    ?         │
│ Z3 GF(2)    │ 6×6  │  100   │    ?        │ attempts/sol    │    ?         │
│ Z3+borda    │ 6×6  │  100   │    ?        │ attempts/sol    │    ?         │
│ BT Warnsd.  │ 8×8  │  100   │    ?        │ dead-ends/sol   │    ?         │
│ Z3 GF(2)    │ 8×8  │  100   │    ?        │ attempts/sol    │    ?         │
└─────────────┴──────┴────────┴─────────────┴─────────────────┴──────────────┘
```

### Curva de escalonamento

Para cada método, plotar:
- Tempo em função de K (quantas soluções pedidas)
- Tempo em função de n (tamanho do tabuleiro)

### Análise de uniformidade das amostras

Para o 6×6 onde temos ground truth:
```python
def uniformidade(amostras, ground_truth_9862):
    """
    Mede quão uniformemente o método amostra o espaço.
    
    Métrica: KL-divergência entre distribuição amostrada
    e distribuição uniforme sobre as 9.862 soluções.
    
    KL = 0: amostragem perfeitamente uniforme
    KL > 0: viés de amostragem
    """
```

Isso responde se Z3 e backtracking exploram o mesmo
espaço ou fatias diferentes.

---

## Hipóteses a testar

**H1: Z3 tem menos desperdício por solução que backtracking**
- BT explora dead-ends explicitamente
- Z3 descarta internamente via SAT — sem dead-end visível
- Esperado: Z3 tem menos "trabalho inútil" por solução

**H2: Backtracking é mais rápido para encontrar a PRIMEIRA solução**
- BT com Warnsdorff encontra a primeira solução em ~O(n²)
- Z3 precisa construir o modelo SAT antes de começar
- Esperado: BT mais rápido para k=1, Z3 mais eficiente para k grande

**H3: Z3 tem viés de amostragem, BT tem viés diferente**
- Z3 explora o espaço em ordem determinada pelo solver
- BT com Warnsdorff favorece soluções "bem comportadas"
- Esperado: ambos têm viés, mas em direções diferentes
- Verificável via KL-divergência contra ground truth 6×6

**H4: Z3 + constraints de borda reduz attempts por solução**
- Injetar as bifurcações descobertas como constraints iniciais
- Esperado: menos attempts desperdiçados com constraints de borda

---

## Experimento adicional: Z3 com invariantes como constraints

Testar se injetar os invariantes descobertos acelera o Z3:

```python
def z3_com_invariantes(board_size, start, end, target_k,
                       invariants):
    """
    Adiciona os invariantes como soft constraints antes de amostrar.
    
    Para cada bifurcação (A, B) com r < -0.5:
      - Adicionar: NOT(A AND B)
        (as duas arestas não podem coexistir)
    
    Isso não elimina soluções válidas (correlação ≠ impossibilidade)
    mas guia o solver para regiões mais densas do espaço.
    """
```

Se H4 for confirmada, você tem uma contribuição prática:
os invariantes topológicos descobertos podem ser usados como
guia de busca em solvers SAT para problemas hamiltonianos.

---

## Saídas esperadas

```
benchmark_results.json     # dados brutos de todos os experimentos
benchmark_report.txt       # tabela comparativa legível
uniformidade_6x6.json      # KL-divergência por método
scaling_curves.json        # tempo × K e tempo × n por método
```

---

## Instruções de execução

```bash
# Rodar todos os benchmarks (estimativa: 30-60 min)
python benchmark_runner.py --repeats 3

# Só 6×6 (rápido, ~5 min)
python benchmark_runner.py --board 6 --repeats 3

# Só comparação de uniformidade
python benchmark_runner.py --uniformidade-only

# Gerar relatório dos dados já existentes
python benchmark_runner.py --report-only
```

---

## Nota importante

Se os dados do backtracking 6×6 já existem em
`destruction_catalogue.json`, extrair dali:
- Tempo total já registrado: 1.360s para 710.064 soluções
- Dead-ends: ~38.5 milhões
- Dead-ends por solução: 38.500.000 / 710.064 ≈ 54

Isso já é um baseline para comparação com Z3.
Só precisar medir o Z3 no mesmo tabuleiro com
o mesmo alvo para ter a comparação direta.
