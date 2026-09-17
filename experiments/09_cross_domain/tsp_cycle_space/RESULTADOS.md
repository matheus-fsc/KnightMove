# Resultados — TSP Cycle-Space vs 2-opt

Benchmark: `product([15,20,25], range(10))` = 30 instâncias, seed mestre 42.

## Tabela agregada

```
=== TSP CYCLE-SPACE vs 2-OPT ===
n=15: cycle_avg= 370.57  2opt_avg= 333.81  gap=+11.44%  time_ratio=1.41x
n=20: cycle_avg= 473.86  2opt_avg= 401.43  gap=+17.89%  time_ratio=1.36x
n=25: cycle_avg= 516.70  2opt_avg= 442.68  gap=+16.72%  time_ratio=1.39x

Best ordering strategy: (empate — todas idênticas)
Cases where cycle-method beats 2-opt: 0 / 30
Gap médio (melhor ciclo vs 2-opt): +15.35%
```

## Resposta à pergunta-chave

> Ordenar os ciclos fundamentais por delta estimado (estratégia *d*) torna a
> busca guiada por ciclos competitiva com o 2-opt?

**Não.** As quatro ordenações (random, length, min_edge, delta) produziram
**custo final idêntico em todas as 30 instâncias**. A ordenação é irrelevante
porque quase não há movimentos válidos para ordenar.

**2-opt vence por ~15,4% em média e é mais rápido** (cycle-method gasta
1,4× o tempo só para varrer ciclos que não pode aplicar).

## Por que o método de ciclos falha (achado central, mecânica verificada)

Apenas **3/30** instâncias tiveram **um único** movimento aplicado; as demais
ficaram presas no tour do vizinho-mais-próximo (`final == init`).

Diagnóstico (n=15, sanity test): dos **91 ciclos fundamentais** da MST, apenas
**3** preservam a Hamiltonicidade quando somados (XOR) ao tour. A mecânica está
correta — confirmamos que um movimento 2-opt, escrito como XOR de um 4-ciclo,
**é** aceito pelo validador. O problema é estrutural:

- Um movimento 2-opt corresponde ao XOR com um **4-ciclo específico**
  `{(a,b),(c,d),(a,c),(b,d)}`.
- Esse 4-ciclo **quase nunca** é um ciclo fundamental da MST.
- Os ciclos fundamentais da MST somados ao tour produzem, na esmagadora
  maioria, grafos com vértices de grau 4 e/ou múltiplas componentes —
  inválidos como tour Hamiltoniano.

Ou seja: **a base do espaço de ciclos derivada da MST não está alinhada com a
vizinhança de movimentos que melhoram o TSP.** A estrutura `H₁(G; F₂)` é rica,
mas o subconjunto de elementos do espaço de ciclos que mapeiam tour→tour é
minúsculo, e a base MST praticamente não o intersecta.

## Conexão com o projeto do passeio do cavalo

Isto ecoa o achado de [[project_ratio_analysis]]: a razão entre tours válidos e
2-fatores (elementos de grau-par do espaço de ciclos) decai exponencialmente.
Aqui vemos o análogo contínuo — somar um elemento arbitrário do espaço de ciclos
a um tour quase sempre o derruba para fora da variedade Hamiltoniana. A
hamiltonicidade (conectividade de componente única) é a obstrução global, não a
estrutura local de grau, exatamente como em [[project_residual_search_done]].

## Direção que poderia salvar a ideia (não implementada)

Usar uma **base de ciclos restrita a 4-ciclos tour-preservantes** (a vizinhança
2-opt explícita), ou somas de pares de ciclos fundamentais que se cancelam nos
vértices de grau-4 — efetivamente reconstruir os movimentos k-opt dentro do
espaço de ciclos. Isso recuperaria o 2-opt, mas não o superaria: seria o 2-opt
reescrito em álgebra GF(2).
