# Agent Prompt — Fase 2: Correção de Viés e Escalonamento

## Contexto da Pesquisa

Esta pesquisa investiga a estrutura topológica interna do espaço de soluções
do Passeio do Cavalo, usando o espaço de ciclos GF(2) como linguagem formal.

### Resultados já obtidos

**6×6 (ground truth exaustivo):**
- H₁ = 45 loops fundamentais
- 710.064 soluções dirigidas = 9.862 ciclos não-direcionados = 1.232 formas D₄-canônicas
- Lei de conservação: destruição total = H₁ = 45 em TODA solução (exato, não aproximado)
- 50 correlações negativas → 7 órbitas D₄ (compactação 7.14×)
- 2 tipos geométricos de bifurcação identificados:
  - Tipo A: vértice de borda escolhendo entre 2 destinos internos (vértice compartilhado)
  - Tipo B: 2 vértices competindo pelo mesmo destino (vértice compartilhado no outro extremo)
- 8 arestas obrigatórias (backbone) = 1 órbita D₄ = 4 cantos grau 2 × 2 arestas cada

**8×8 (amostragem Z3, ~800k amostras por par, 8 pares canônicos):**
- H₁ = 105 loops fundamentais
- Compactação D₄ média: 7.88× (similar ao 6×6 → hipótese suportada)
- ~13 órbitas distintas por par (vs 7 no 6×6, crescimento modesto)
- Top invariante em TODOS os 8 pares: B8-D7 ↔ B8-C6 (análogo direto ao B6-D5 ↔ B6-C4 do 6×6)
- Mesma geometria de bifurcação Tipo A preservada

### Hipótese central
> O número de invariantes topológicos reais (correlações negativas módulo D₄)
> cresce sub-linearmente com o tamanho do tabuleiro. Os tipos geométricos de
> bifurcação são universais — aparecem em qualquer tabuleiro n×n.

### Problemas identificados que precisam ser resolvidos
1. **Viés de amostragem Z3**: Z3 não amostra uniformemente; prefere regiões
   "fáceis" do espaço SAT. Os resultados do 8×8 podem ser artefatos.
2. **Quebra de simetria no engine**: `break_symmetry=True` força a aresta
   A8-C7 ativa para pares com start=(0,0), enviesando 4 dos 8 pares canônicos.
3. **Apenas 2 pontos de dados**: 6×6 e 8×8 não são suficientes para
   estabelecer a forma da curva de crescimento.

### Base de código disponível
- `engine_8x8.py`: engine Z3 para caminhos hamiltonianos abertos
- `correlator.py`: correlações de Pearson vetorizadas (numpy) + colapso D₄
- `expand_d4.py`: expansão D₄ de assinaturas sem re-amostrar
- `d4_orbits.py`: grupo D₄ sobre vértices, arestas e pares de arestas
- `parallel_runner.py`: multiprocessing.Pool com checkpoint retomável
- `checkpoint.py`: estado persistente (atomic write)
- `destruction_8x8.py`: análise de destruição de loops por casa
- `cavalo_loop_destruicao_6x6.py`: enumeração exaustiva 6×6
- `analysis_6x6.json`: baseline completo do 6×6

---

## Tarefa do Agente — 3 Fases Sequenciais

---

## FASE A: Validação do Método (Prioridade Máxima)

**Objetivo**: confirmar que a amostragem Z3 reproduz os invariantes corretos,
comparando com o ground truth exaustivo do 6×6.

**Por quê fazer primeiro**: se o método Z3 for enviesado no 6×6 (onde temos
a verdade), não faz sentido confiar nos resultados do 8×8.

### A1 — Engine Z3 para o 6×6

Adaptar `engine_8x8.py` para o tabuleiro 6×6:

```python
# engine_6x6_sampler.py
# Idêntico ao engine_8x8.py mas com BOARD=6
# CRÍTICO: break_symmetry=False por padrão nesta versão
# Motivo: queremos comparar com ground truth sem viés adicional

def sample_signatures_6x6(start, end, n_target,
                           break_symmetry=False,  # False por padrão
                           ...):
    ...
```

### A2 — Coleta de amostras Z3 no 6×6

Rodar o parallel_runner adaptado para 6×6 com os mesmos pares canônicos
usados na análise exaustiva. Alvo: 10.000 amostras por par (rápido no 6×6).

```bash
python parallel_runner_6x6.py \
  --workers 4 \
  --samples-per-batch 200 \
  --n-batches-per-pair 50 \
  --break-symmetry false \
  --out-dir data_6x6_sampled
```

### A3 — Comparação ground truth vs Z3

Implementar `validate_sampler.py` que:

```python
def compare_invariants(ground_truth_json, sampled_json):
    """
    Compara os invariantes do ground truth (analysis_6x6.json)
    com os invariantes obtidos por amostragem Z3.

    Métricas de comparação:
    1. Órbitas encontradas: são as mesmas 7?
    2. Ranking por |r|: a ordem dos invariantes bate?
    3. Valores de r: diferença absoluta média entre GT e Z3
    4. Arestas obrigatórias: as mesmas 8?
    5. Distribuição de frequências de arestas: KL-divergence GT vs Z3

    Saída:
    {
      "orbits_match": bool,          # mesmas 7 órbitas?
      "ranking_correlation": float,  # Spearman rank corr dos r-values
      "r_mae": float,                # mean absolute error dos r-values
      "mandatory_match": bool,       # mesmas 8 arestas obrigatórias?
      "freq_kl_divergence": float,   # KL-div das frequências de arestas
      "verdict": "VALID" | "BIASED" | "PARTIAL"
    }
    """
```

**Critério de aprovação**: `ranking_correlation > 0.90` e `orbits_match = True`.
Se não passar: identificar a fonte do viés e corrigir antes de prosseguir.

### A4 — Testar com e sem break_symmetry

Rodar com `break_symmetry=True` e `break_symmetry=False` e comparar.
Quantificar o viés introduzido pela quebra de simetria.

**Resultado esperado**: com `break_symmetry=False`, os resultados devem
ser mais próximos do ground truth. Isso confirmará que a quebra de simetria
é a principal fonte de viés nos resultados do 8×8.

---

## FASE B: 8×8 Sem Viés (após aprovação da Fase A)

**Objetivo**: re-coletar dados do 8×8 com o método validado.

### B1 — Engine corrigido

Modificar `engine_8x8.py`:

```python
# Remover break_symmetry hardcoded para start=(0,0)
# Substituir por amostragem sem viés direcional

# OPÇÃO 1 (simples): simplesmente desligar break_symmetry
sample_signatures(..., break_symmetry=False, ...)

# OPÇÃO 2 (melhor): randomizar a simetria quebrada
# Em vez de sempre forçar A8-C7, sortear aleatoriamente qual das
# arestas do start será forçada. Isso distribui o viés uniformemente.
if break_symmetry:
    edges_from_start = [e for e in edges_aug if start in e and e != virtual]
    forced_edge = random.choice(edges_from_start)
    solver.add(ea(forced_edge) == True)
```

A Opção 2 mantém a eficiência computacional da quebra de simetria
(reduz o espaço de busca do Z3 à metade) sem introduzir viés direcional.

### B2 — Pares canônicos mais abrangentes

Ampliar para cobrir mais órbitas D₄ de pares de casas:

```python
CANONICAL_PAIRS_EXTENDED = [
    # Originais (manter para comparação)
    [[0,0],[0,7]],  [[0,0],[7,6]],  [[0,0],[4,3]],  [[0,0],[2,5]],
    [[0,1],[6,4]],  [[0,1],[7,5]],  [[0,1],[4,4]],  [[1,2],[7,5]],

    # Novos: start no centro (testa se invariantes são do grafo, não do par)
    [[3,3],[0,1]],  # centro → borda
    [[3,3],[0,0]],  # centro → canto (paridade?)
    [[2,3],[5,4]],  # interior → interior oposto
    [[3,4],[4,3]],  # vizinhos do centro
]
```

**Pergunta-chave**: os mesmos ~13 invariantes aparecem quando o start
está no centro? Se sim, são propriedades do grafo. Se não, são do par.

### B3 — Volume de amostras

Por par canônico, alvo de 200.000 amostras (antes eram 800.000 com
expansão D₄; aqui 200.000 originais × 8 transformações D₄ = 1.6M).

```python
CONFIG_B3 = {
    "samples_per_batch": 200,
    "n_batches_per_pair": 1000,  # 200k amostras por par
    "break_symmetry": False,      # método validado
}
```

### B4 — Análise de convergência

Verificar se os invariantes estabilizam com menos amostras:

```python
def convergence_analysis(samples_dir, checkpoints=[1000, 5000, 20000,
                                                    50000, 200000]):
    """
    Para cada tamanho de amostra, calcula os invariantes e verifica
    se as órbitas D₄ e rankings estabilizaram.

    Plota: N_amostras × {n_orbitas, r_top1, r_top2, r_top3}
    Determina N_mínimo para estabilização (critério: variação < 5%
    entre dois checkpoints consecutivos).
    """
```

Isso informa o volume mínimo necessário para o 10×10.

---

## FASE C: Escalonamento para 10×10

**Objetivo**: terceiro ponto de dados para testar a forma da curva de crescimento.

**Pré-requisito**: Fase A aprovada + Fase B concluída.

### C1 — Engine 10×10

```python
# engine_10x10.py
BOARD = 10
# H1 = |E| - |V| + 1
# Para 10×10: |V|=100, |E|=?
# Arestas do cavalo: fórmula analítica
# E = 4(n-1)(n-2) = 4×9×8 = 288  → H1 = 288 - 100 + 1 = 189

# Atenção: o grafo do cavalo 10×10 tem loops menores que 8×8
# Verificar se o spanning tree do engine escala bem
```

### C2 — Pares canônicos 10×10

```python
CANONICAL_PAIRS_10x10 = [
    [[0,0],[0,9]],  # canto → canto mesma borda
    [[0,0],[9,8]],  # canto → canto diagonal
    [[0,0],[5,4]],  # canto → centro
    [[0,1],[9,6]],  # borda → borda oposta
    [[0,1],[5,5]],  # borda → centro
    [[4,4],[5,5]],  # centro → centro adjacente
]
```

### C3 — Volume e infraestrutura

Com H₁ = 189 e grafo maior, o Z3 fica mais lento por amostra.
Estimativa: ~5-10× mais lento que 8×8 por amostra.

Para VPS múltiplos, configurar distribuição de pares por servidor:

```bash
# VPS 1: pares 0-1 do 10×10
python parallel_runner.py --pairs-subset 0,1 --workers 8

# VPS 2: pares 2-3 do 10×10
python parallel_runner.py --pairs-subset 2,3 --workers 8

# VPS 3: mais amostras do 8×8 (Fase B)
python parallel_runner.py --config config_8x8_phase_b.json --workers 8
```

Implementar `--pairs-subset` no parallel_runner:

```python
ap.add_argument("--pairs-subset", default=None,
    help="índices dos pares a processar, ex: '0,1,3' (default: todos)")
```

### C4 — Análise comparativa final

Com 3 pontos de dados (6×6, 8×8, 10×10), rodar regressão sobre
o crescimento das órbitas:

```python
def fit_growth_curve(data_points):
    """
    data_points = [(n, n_orbits), (6,7), (8,13), (10,?)]

    Ajustar modelos:
    - Linear:       f(n) = a*n + b
    - Quadrático:   f(n) = a*n² + b*n + c
    - n²/8:         f(n) = n²/8 + c    ← hipótese do paper
    - n*log(n):     f(n) = a*n*log(n) + b

    Reportar R², AIC, BIC para cada modelo.
    Não afirmar qual é o "certo" com 3 pontos — apenas qual ajusta melhor.
    """
```

---

## Métricas de Sucesso por Fase

### Fase A — Validação
- [ ] `ranking_correlation > 0.90` entre GT e Z3 no 6×6
- [ ] Mesmas 7 órbitas encontradas
- [ ] Viés de `break_symmetry` quantificado
- [ ] Método aprovado ou corrigido

### Fase B — 8×8 Sem Viés
- [ ] Resultados com `break_symmetry=False` reproduzem padrão anterior
- [ ] Invariantes com start no centro identificados
- [ ] N_mínimo de convergência determinado
- [ ] Confirmação (ou refutação) dos ~13 invariantes

### Fase C — 10×10
- [ ] H₁ = 189 confirmado
- [ ] N_orbitas medido (ponto 3 da curva)
- [ ] Forma do crescimento ajustada com 3 modelos
- [ ] Tipo geométrico das bifurcações verificado (mesmo padrão?)

---

## Infraestrutura VPS

### Configuração recomendada para múltiplos servidores

```
VPS Principal (maior RAM):
  - Roda correlator.py (consome ~2-4GB RAM para 1.6M amostras)
  - Consolida resultados de todos os VPS
  - Roda report.py e validate_sampler.py

VPS Workers (2-4 servidores):
  - Rodam parallel_runner.py com --pairs-subset
  - Escrevem batches em diretório local
  - Rsync periódico para VPS Principal:
    rsync -avz data_parallel/ vps_main:/data/knight_tour/8x8/

Script de sincronização (rodar a cada 30min via cron):
  #!/bin/bash
  rsync -avz --progress \
    data_parallel/samples/ \
    user@vps_main:/data/knight_tour/samples/
  ssh user@vps_main "cd /data/knight_tour && python correlator.py --expand-d4"
```

### Estimativas de tempo por configuração

| Tabuleiro | Amostras/par | Workers | Tempo estimado |
|-----------|-------------|---------|----------------|
| 6×6 validação | 10k | 4 | ~30 min |
| 8×8 Fase B | 200k/par × 12 pares | 8 × 2 VPS | ~48h |
| 10×10 Fase C | 50k/par × 6 pares | 8 × 2 VPS | ~72-120h |

---

## Perguntas de Pesquisa a Responder

Ao final das 3 fases, o relatório deve responder:

1. **O método Z3 é confiável?**
   Resposta quantificada pela comparação com GT do 6×6.

2. **Os ~13 invariantes do 8×8 são reais ou artefatos?**
   Resposta: reproduced com método corrigido + start no centro.

3. **Os tipos geométricos de bifurcação são universais?**
   Resposta: verificado em 3 tamanhos (6×6, 8×8, 10×10).

4. **Qual é a forma da curva de crescimento?**
   Resposta: regressão com 3 pontos, modelos comparados por AIC/BIC.

5. **A complexidade topológica é O(n²/8)?**
   Resposta: confirmada, refutada, ou inconclusiva com os dados disponíveis.

---

## Notas para o Agente

- **Não pular a Fase A**. É a mais importante. Sem validação do método,
  os resultados das fases seguintes não têm base.

- **Manter compatibilidade de formato** com os JSONs existentes
  (analysis_6x6.json, raw_correlations.json, orbit_correlations.json).
  O report.py deve funcionar com os novos dados sem modificação.

- **Checkpoint sempre**. Cada fase pode ser interrompida e retomada.
  O parallel_runner.py já implementa isso — usar consistentemente.

- **Não tirar conclusões com 2 pontos**. Toda afirmação sobre curvas
  de crescimento deve esperar o 10×10.

- **Documentar o viés quantitativamente**. Se break_symmetry introduz
  viés, medir o tamanho do viés em termos de KL-divergência ou
  diferença nos r-values. Isso vai para o paper como metodologia.
