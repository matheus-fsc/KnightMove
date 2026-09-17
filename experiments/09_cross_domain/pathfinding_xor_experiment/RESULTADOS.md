# Resultados — GF(2) XOR vs A* repetido em mapas-grade

60 instâncias (3 tamanhos × 20 seeds), 30% de obstáculos.

## Tabela agregada

| Grade | Compat. single-XOR | Caminhos XOR (méd) | Caminhos A* (méd) | t_XOR | t_A* | Diversidade XOR | Diversidade A* |
|------:|:------------------:|:------------------:|:-----------------:|------:|-----:|:---------------:|:--------------:|
| 20×20 | **26.9% ± 6.0%**   | 66.8 (20–112)      | 5.6               | 23.5ms | 2.0ms | 0.298 | 0.453 |
| 30×30 | **15.1% ± 4.1%**   | 53.5 (32–121)      | 9.3               | 49.3ms | 6.3ms | 0.262 | 0.521 |
| 50×50 | **11.7% ± 2.0%**   | 84.8 (57–108)      | 10.0              | 115.7ms | 15.5ms | 0.170 | 0.456 |

Compatibilidade global média: **17.9%**. Todas as 60 instâncias produziram
≥10 caminhos válidos via XOR. O A* repetido achou os 10 em apenas 30/60
(falha quando a componente é pequena e as arestas se esgotam).

## Respostas às perguntas-chave (honestas)

**1. Taxa de compatibilidade single-XOR em grades?**
**11.7%–26.9%** (média 17.9%) — claramente **>> 3%** do TSP. Hipótese de que
grafos esparsos são mais "amigáveis" ao XOR: **CONFIRMADA**. Mas note: a taxa
**cai** quando a grade cresce (26.9% → 11.7%), porque os ciclos da árvore DFS
ficam mais longos e menos locais.

**2. XOR é mais rápido que A* repetido para 10 caminhos?**
**NÃO**, do jeito medido (t_XOR cobre enumerar *todos* os ~50–100 caminhos +
500 pares; t_A* cobre só 10). Per-caminho a leitura muda: o A* é ~8× mais
rápido por caminho, **mas** só consegue 10 caminhos em metade dos casos,
enquanto o XOR entrega 50–100 alternativas distintas no mesmo orçamento de
DFS único. O XOR ganha em **vazão de alternativas**, não em latência.

**3. XOR produz caminhos mais diversos que a penalização A*?**
**NÃO.** A* penalizado é consistentemente mais diverso (0.45–0.52 vs
0.17–0.30). Faz sentido: o XOR com ciclo fundamental faz edições **locais**
no caminho-base (troca um corredor), mantendo a maior parte das arestas; a
penalização A* força desvios **globais**. Os caminhos XOR são variações
locais; os caminhos A* são rotas estruturalmente distintas.

**4. Compatibilidade correlaciona com densidade ou comprimento de ciclo?**
- densidade efetiva de obstáculos: **+0.343** (mais buracos → ciclos menores → mais compatível)
- comprimento médio do ciclo: **−0.665** (ciclos mais longos → menos compatível)
- fração de ciclos locais (≤10 arestas): **+0.499**

Confirma a intuição: **ciclos curtos/locais → maior compatibilidade**. O
comprimento do ciclo é o preditor mais forte.

**5. Modo de falha dominante?**
**DESCONEXÃO**, esmagadoramente: de 16.847 candidatos inválidos →
**90.8% desconexos**, 5.5% violação de grau, 3.7% ambos. Quando o ciclo DFS
não compartilha um trecho contíguo com o caminho-base, o XOR gera o caminho +
um laço solto separado → componente desconexa.

## Conexão com o trabalho existente

**vs. resultado TSP (~3%):** a compatibilidade é ~4–9× maior em grades
esparsas. A esparsidade planar realmente favorece o XOR — confirma que o
fracasso no TSP era específico de grafos densos, não uma limitação universal
do GF(2)/XOR.

**vs. achado do Knight's Tour (conexão = obstrução global):** o modo de falha
casa **exatamente**. No Knight's Tour, a conectividade era a obstrução global
restante (poucos 2-fatores são conexos: 73.5% no 6×6 caindo para 5.7% no
10×10 — ver `[[project_residual_search_10x10_done]]`). Aqui, **90.8%** das
falhas do XOR são desconexão, não grau. A restrição de **grau** (local) é
quase sempre satisfeita de graça (a paridade de grau é preservada pelo XOR
com um ciclo); o que mata o candidato é a **conectividade** (global). Mesmo
padrão em dois domínios independentes.

**Conclusão — domínio-específico, não universalmente limitado:**
o XOR de ciclos GF(2) é **viável em grafos planares esparsos** (compat ~12–27%)
e **inviável em grafos densos** (TSP ~3%). Em ambos os domínios a obstrução
que sobra é a **conectividade global** — a estrutura de grau local nunca é o
gargalo. Isso unifica os dois experimentos: GF(2)/XOR resolve trivialmente as
restrições locais; o desafio que persiste, em qualquer densidade, é manter o
grafo conexo.
