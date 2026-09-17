# Viabilidade do n=8 — tamanho do espaço de estados da T_bulk de graus

**Script:** `measure_states_n8.py` · **Dados:** `data/measure_states.json`
**Medição pura de estados** (sem matriz, sem espectro, sem contagem). Estado =
`(deg_c, deg_{c+1})`, deg ∈ {0,1,2}ⁿ, bitpacked em int (2 bits/célula).

> Objeto medido: T_bulk **de graus** (cega à conectividade). NÃO é o
> broken-profile de `knight_transfer.cpp` (que rastreia conexão e estoura
> >187M em n=8). São objetos distintos.

## Gate de fidelidade (PASSO 1) ✓

n=6 reproduziu **33.346 estados** (esperado 33.346) em 10,7 s pela mesma BFS de
`transfer_build.transitions`. Enumeração confiável.

## Curva de crescimento (PASSO 2)

| n | estados | status | ratio K(n)/K(n−1) | ratio mesma paridade K(n)/K(n−2) | tempo |
|---|---|---|---|---|---|
| 4 | 342 | completo | — | — | — |
| 5 | 7.322 | completo | ×21,4 | — | — |
| 6 | 33.346 | completo | ×4,55 | ×97,5 (n6/n4) | 10,7 s |
| 7 | 516.417 | completo | ×15,5 | ×70,5 (n7/n5) | 479 s |
| 8 | **≥ 1.631.530** | **time_cap** | ×≥3,2 (parcial) | — | 793 s (cortado) |

**Alternância par/ímpar nítida:** saltos para n ímpar são grandes (×21, ×15);
saltos para n par são menores (×4,55 em 5→6). O passo 7→8 é par ⇒ esperado
menor que os saltos ímpares.

### Curva interna do n=8 (cortado em 793 s / cap 14 GB nunca atingido)

| t (s) | estados descobertos | expandidos | fila | RSS |
|---|---|---|---|---|
| 19 | 514.066 | 10.883 | 503.183 | 0,10 GB |
| 266 | 1.014.076 | 80.577 | 933.499 | 0,15 GB |
| 742 | 1.514.336 | 316.796 | 1.197.540 | 0,22 GB |
| **793 (corte)** | **1.631.530** | 320.000 | 1.311.530 | **0,22 GB** |

No corte só **20%** dos estados descobertos foram expandidos (320 k de 1,63 M;
fila de 1,31 M ainda não processada). O número real é bem maior que 1,63 M.

## Número do n=8: estimativa

- **Cota inferior dura (medida):** **1.631.530 estados.**
- **Extrapolação:** os ratios de mesma paridade decaem (97,5 → 70,5); projetando
  K(8)/K(6) ≈ 55–65 daria ~1,8–2,2 M. Mas como já há 1,63 M descobertos com 80%
  por expandir, e pela saturação tardia da BFS observada no n=7 (a 30% expandido
  havia 77% dos estados finais), o total realista é **~2,5–3,5 milhões de
  estados** (central ~2,5 M).

## Memória estimada (n=8 a ~2,5 M estados)

| representação | por estado | total | nota |
|---|---|---|---|
| **tuplas Python** (dict `{state:i}`) | ~350 B | **~0,9 GB** | + fila/visited na BFS |
| **bitpacked em RAM** (int 32 bits) | 4 B cru | **~10 MB cru** / ~200 MB (set Python) / ~30–50 MB (hashset C) | n=8 cabe em 32 bits |
| **matriz esparsa** (nnz ≈ estados × out-deg) | — | **~4–13 GB** (central ~7 GB) | **o verdadeiro gargalo de memória** |

Out-degree médio: n=6 = 1.966.428/33.346 = **59**; cresce ~×2 por linha extra,
projetando ~150–250 para n=8. Logo nnz(8) ≈ 2,5 M × ~230 ≈ **~575 M triplas** →
~7 GB a 12 B/tripla (a construção em Python usaria transitoriamente bem mais).

## Ambiente e gargalo

- **Python puro.** O `transitions` recursivo (combinações por linha) custa
  ~1 ms/estado em n=7 e ~3,3 ms/estado em n=8.
- **O gargalo é VELOCIDADE, não memória.** RSS ficou em **0,22 GB** para 1,63 M
  estados — o cap de 14 GB nunca chegou perto. A 14 GB caberiam ~100 M estados.
  O que mata é o tempo: a ~300 estados/s, expandir ~2,5 M estados levaria
  **~2,3 horas** só de enumeração em Python.

## Veredito de caixa (PASSO 3)

**A pergunta de hardware está mal-posta: o n=8 de 2-fatores é limitado por
SOFTWARE (Python lento), não por hardware.**

1. **Contagem de 2-fatores (o alvo):** NÃO precisa materializar a matriz. O DP
   coluna-a-coluna (`count_2factors_dp`) faz streaming das transições e guarda só
   o vetor de estados (~2,5 M doubles ≈ **20 MB**). Memória trivial. O único
   bloqueio é a velocidade do `transitions` em Python.
   → **Viável numa única máquina (até o ROG sozinho, ou um laptop), SEM GPU e SEM
   cluster, após reescrever `transitions` em C/numba.** A tabela de estados
   (~10–200 MB bitpacked) cabe em <8 GB com folga.

2. **Matriz explícita (só se quiser ESPECTRO):** ~7 GB (faixa 4–13 GB).
   - **<8 GB:** não — a matriz sozinha passa de 8 GB ⇒ **não cabe na VRAM da
     4060** (8 GB), nem com folga para o workspace do eigensolver.
   - **8–16 GB:** sim — **cabe na RAM do ROG** com bitpacking + C/numba.
   - **>16 GB:** só no pior caso (≥4 M estados e out-deg ~300 ⇒ nnz ~1,5 B ⇒
     ~18 GB); não é o cenário central.

3. **Cluster heterogêneo (32 GB):** desnecessário. Fator de paralelismo não ajuda
   um BFS sequencial intrinsecamente serial; e a memória nunca foi o problema.

### Resumo de uma linha

O espaço de estados da T_bulk de graus em n=8 é **~2,5 M (cota inferior medida
1,63 M)** — pequeno em memória (tabela bitpacked ~10–200 MB). O n=8 de 2-fatores
é alcançável **numa única máquina**, mas o gargalo é a **velocidade da
enumeração em Python puro**: a solução é reescrever `transitions` em C/numba, não
comprar hardware. GPU/cluster só seriam relevantes para a matriz explícita do
espectro (~7 GB, cabe no ROG, não na 4060) — e mesmo aí não para a contagem.
