# Resultados — Path-finding via Loop-XOR vs. heurísticas clássicas

Gerado por `path_decomposition.py` (dados brutos em
`results/path_decomposition_experiment.json`). Reportado honestamente,
**favorável ou não** à hipótese.

## Setup efetivo

| | 6×6 | 8×8 |
|---|---|---|
| V, E | 36, 80 | 64, 168 |
| Ciclos fundamentais (β₁ = E−V+1) | 45 | 105 |
| s → t **usado** | (0,0)→(0,1) | (0,0)→(0,1) |
| Caminho-base achado | sim | sim |

**Paridade:** `(0,0)→(n−1,n−1)` é **impossível** em ambos: têm a mesma cor,
mas um caminho hamiltoniano tem `n²−1` arestas (ímpar), exigindo extremos de
cores opostas. O script detecta e cai para `(0,0)→(0,1)`, como previsto.

## Números principais

| Métrica | 6×6 | 8×8 |
|---|---|---|
| **Ground truth** (caminhos s→t, exaustivo) | **46 666** (13.3 s) | — (8×8 inviável) |
| XOR — caminhos válidos distintos | 4 | 10 |
| XOR — cobertura | **0,0086 %** | — |
| XOR — tempo | 51 ms | 63 ms |
| Ciclos individualmente compatíveis | **0 / 45 (0 %)** | **4 / 105 (3,8 %)** |
| Backtracking (60 s cap no 8×8) | 46 666 exaustivo | **0** completados |
| Warnsdorff — taxa de sucesso (t fixo) | 2,35 % (47/2000) | 0,55 % (11/2000) |
| Warnsdorff — caminhos distintos | 3 | 11 |
| **Diversidade** XOR | 9,33 | 9,07 |
| Diversidade Warnsdorff | 10,67 | 43,38 |
| Diversidade população real (amostra bt) | **19,43** | — |

### Modo de falha do XOR (de ~10 000 subsets testados)

| | 6×6 | 8×8 |
|---|---|---|
| `valid` | 0,1 % | 0,7 % |
| `bad_count` (|arestas| ≠ V−1) | **73,2 %** | **84,0 %** |
| `bad_degree` (algum vértice grau ≠ esperado) | 24,8 % | 14,5 % |
| `disconnected` (graus ok, mas caminho + ciclos) | 1,9 % | 0,8 % |

## Respostas às perguntas-chave

**1. O XOR encontra caminhos válidos, ou a maioria das combinações produz
subgrafos inválidos/desconexos?**
Esmagadoramente **inválidos**. Apenas **0,1 % (6×6) / 0,7 % (8×8)** dos subsets
de ciclos produzem um caminho hamiltoniano válido. O modo de falha dominante
**não é desconexão**, e sim **`bad_count`** (73–84 %): o XOR do caminho-base
com ciclos fundamentais arbitrários **muda o número de arestas**. Razão: para
que `P ⊕ C` continue tendo `V−1` arestas é preciso que `C` *alterne* com `P`
(metade das arestas de `C` em `P`); os ciclos da base de `cycle_basis` quase
nunca alternam, então o resultado tem `≈ (V−1)+|C|` arestas. Em seguida vem
`bad_degree` (vértices de grau 4); desconexão pura é rara (< 2 %).

**2. A diversidade dos caminhos XOR é maior que a do Warnsdorff repetido?**
**Não.** A diversidade XOR (≈ 9 arestas de diferença média) é a **mais baixa**:
abaixo do Warnsdorff (10,7 no 6×6; 43,4 no 8×8) e bem abaixo da **população
verdadeira** (19,4 no 6×6). Os caminhos XOR ficam **agrupados em torno do
caminho-base** — diferem dele apenas por um ciclo alternante curto. Ressalva:
amostras pequenas (4–11 caminhos), mas o sinal é consistente.

**3. Que fração dos ciclos é "compatível" com o caminho-base?**
Pouquíssima: **0/45 (6×6)** e **4/105 (3,8 %, 8×8)** ciclos isolados geram
caminho válido ao serem XOR-ados. Confirma que a base de ciclos fundamentais
está **desalinhada** com o caminho — só ciclos *alternantes* funcionam, e a
base do DFS quase não contém nenhum.

## Conclusão (honesta)

A decomposição em ciclos fundamentais + XOR **falha como enumerador** de
caminhos hamiltonianos `s→t` no grafo do cavalo:

- cobre **0,0086 %** do conjunto-verdade no 6×6 e satura em ~10 caminhos no
  8×8 (limitada pelos poucos ciclos compatíveis — não melhora com mais tempo);
- gera caminhos **menos diversos** que tanto Warnsdorff quanto a população real;
- a operação XOR só preserva a estrutura de caminho para **ciclos alternantes**,
  ausentes na base `cycle_basis` — então 99,3–99,9 % das combinações são lixo.

O baseline **backtracking com poda de grau + conectividade** é o vencedor claro
no 6×6 (enumera os 46 666 caminhos em 13 s, exato). No 8×8 nenhum método
enumera bem com extremo fixo: o DFS exaustivo não fecha nenhum caminho em 60 s
e o Warnsdorff tem ~0,5 % de sucesso — o XOR produz uns poucos caminhos
baratos, mas é um teto, não uma vantagem.

Alinha com o tema do projeto: **conectividade/estrutura global é a obstrução**;
operações algébricas locais (XOR de ciclos) não capturam a restrição de "um
único caminho gerador". Para o XOR funcionar seria preciso restringir a
**ciclos alternantes** em relação a `P` (ciclos aumentantes), não à base de DFS.
