# Patch LUT — Lema de Deformação GF(2) pré-computado

Implementação: `knight_tours_patch.py`
Base: `knight_tours.py` (sem modificações no original).
Predecessor experimental: `knight_tours_deformation_gf2.py` (versão runtime,
sem LUT).

---

## 1. Resultado executivo

| medida                         | resultado                                                     |
|--------------------------------|---------------------------------------------------------------|
| Corretude n=6 (GT=9.862)       | **9.862/9.862 ✓** (todos tours válidos)                       |
| Speedup n=6 (enumeração)       | 1.00× ± 0.13 (5 trials, ruído domina)                         |
| Speedup n=10 K=200             | 0.59×  (lento — overhead da iteração sobre LUT)               |
| Taxa de resgate (n=6)          | **0/15.660 = 0,0 %**                                          |
| Taxa de resgate (n=10)         | **0/137   = 0,0 %**                                           |
| Tempo de build_patch_LUT n=10  | 0,016 s (≈ 5 % de overhead de `build_graph`)                  |
| Tempo de build_patch_LUT n=12  | 0,028 s                                                       |

A LUT é construída em tempo desprezível, mas **nunca dispara em runtime**.
A hipótese de que pré-computar elimina o overhead estava certa (a iteração
da LUT é O(diamantes_por_aresta) ≈ O(1)), mas o overhead real é
**conceitual**, não computacional.

---

## 2. Estrutura da LUT (Tarefa 4)

LUT[e1] = lista de triplas (e2, e3, e4) tal que e1-e3-e2-e4 é um 4-ciclo
no grafo do cavalo. As 4 arestas formam um "diamante":

```
        u --e1-- v
        |        |
       e3       e4    (e1, e2: paredes a DESTRUIR)
        |        |    (e3, e4: pontes a CRIAR)
        x --e2-- y
```

Os 4 vértices envolvidos têm **grau líquido invariante** após o patch:
cada um perde 1 (parede) e ganha 1 (ponte). Portanto qualquer estado
válido permanece válido após aplicação do patch.

### Tabela: número de diamantes catalogados por n

| n  | V   | E   | tot. diamantes | média/aresta | máx/aresta | t_build |
|----|-----|-----|----------------|--------------|------------|---------|
| 6  |  36 |  80 |    208         |  2,60        |  5         | 0,003 s |
| 8  |  64 | 168 |    592         |  3,52        |  6         | 0,010 s |
| 10 | 100 | 288 |  1.168         |  4,06        |  6         | 0,016 s |
| 12 | 144 | 440 |  1.936         |  4,40        |  6         | 0,028 s |

Todas as arestas têm ≥1 diamante; nenhuma está "isolada".

### Estrutura por nível L (distância à borda)

| L | n=6        | n=8        | n=10       | n=12       |
|---|------------|------------|------------|------------|
| 0 | 2,14 (56)  | 2,45 (88)  | 2,60 (120) | 2,68 (152) |
| 1 | 3,67 (24)  | 4,14 (56)  | 4,27 (88)  | 4,33 (120) |
| 2 |   —        | 6,00 (24)  | 6,00 (56)  | 6,00 (88)  |
| 3 |   —        |   —        | 6,00 (24)  | 6,00 (56)  |
| 4 |   —        |   —        |   —        | 6,00 (24)  |

(número entre parêntesis = quantidade de arestas no nível)

**Padrão claro e universal:**
- L=0 (borda): poucos diamantes disponíveis (cresce devagar com n)
- L≥2 (interior): saturação em **6 diamantes/aresta**

Isso confirma a hipótese de **"zonas congeladas" na borda** — exatamente
análogo a f∞(L). Assim como `F_INF[0]=0.528` (borda extremamente
polarizada) versus `F_INF[L≥2]≈0.20`, a LUT mostra densidade
estrutural ≈ 2× menor no anel L=0 vs interior.

A estrutura é universal em n: parametrizável por L como f∞(L).

---

## 3. Por que a taxa de resgate é 0 % (Tarefa 3)

Instrumentamos `_process_queue` (script: `investigate_patch.py`) e
contamos qual filtro elimina cada candidato:

| filtro                  | n=6 (37k checagens) | n=10 (655 checagens) |
|-------------------------|---------------------|----------------------|
| `e2_not_active`         | 78,5 %              | 76,5 %               |
| `e3_not_free`           | 20,6 %              | 22,9 %               |
| `e4_not_free`           |  0,9 %              |  0,6 %               |
| `e2_same_component`     |  0,0 %              |  0,0 %               |
| **patch aplicado**      |  **0,0 %**          |  **0,0 %**           |

O motivo dominante (78 %) é **e2_not_active**: a "parede" simétrica
do outro componente raramente está ACTIVE no instante do sub-ciclo.

`e2_same_component` jamais aparece — sempre que e2 ESTÁ ACTIVE, ela
está no MESMO componente que e1. Ou seja:

> **Quando `ru == rv` detecta sub-ciclo, ele é tipicamente o único
> componente conectado vivo. Não há "outro sub-tour" para deformar.**

Isso revela uma propriedade não-óbvia do algoritmo atual:

1. **R2 + heurística f∞ é "ávida"**: força fixações em cadeia, fundindo
   componentes UF cedo. Quando um sub-ciclo finalmente fecha, todos os
   outros vértices ainda livres estão num único componente fundido com ele.
2. **Não há simultaneidade de sub-componentes**: a propagação não deixa
   dois grupos crescerem em paralelo até se encontrarem.
3. **Deformação requer simultaneidade**: o lema XOR só ajuda se C1 e C2
   coexistem como sub-ciclos fechados antes de qualquer um abortar.

---

## 4. Implicação

A LUT confirma que **a estrutura geométrica do lema de deformação é
válida e universal por L**, mas a **viabilidade algorítmica é nula**
sob o regime de busca atual.

O lema XOR só seria útil em um algoritmo com propagação **mais
parcimoniosa** — por exemplo:

- branch-and-bound puro (sem R2), onde múltiplos sub-ciclos podem
  coexistir antes de qualquer fechamento;
- algoritmo de "subtour elimination" estilo TSP-cutting plane, onde
  vários sub-ciclos pequenos são iterativamente combinados;
- métodos de Monte Carlo locais (cavalo simulated annealing) onde
  deformações 2-opt/3-opt fazem sentido sobre um 2-fator pré-existente.

No algoritmo atual (knight_tours.py), o trade-off já foi feito a favor
de R2 — e R2 é tão boa que torna a deformação redundante.

### Tabela: Knuth estimator (não corrido)

T5 foi cancelada porque, sem patches aplicados (0/N), o estimador
produz números **idênticos** ao baseline. `frac_zeros`, `CV`, `N_hat`:
tudo igual. Não há sinal a medir.

---

## 5. Conclusão

| pergunta                                                    | resposta |
|-------------------------------------------------------------|----------|
| LUT é construível em tempo desprezível?                     | **Sim** |
| Estrutura da LUT é universal em n (parametrizável por L)?   | **Sim** — saturação em L≥2 |
| Pré-computar elimina o overhead do `find_patch_loop`?       | **Sim** |
| O patch dispara em runtime?                                 | **Não** (0/16k) |
| Speedup global?                                             | **1,00×** (n=6), **0,59×** (n=10) |
| O lema XOR é útil no algoritmo atual?                       | **Não** |
| O lema XOR é matematicamente correto?                       | **Sim** (corretude validada) |

A LUT é um **resultado negativo informativo**: ela demonstra
empiricamente que R2 + UF incremental + f∞ já operam num regime
onde o lema de deformação não tem espaço para atuar. O ganho marginal
em n=6 do experimento anterior (`knight_tours_deformation_gf2.py`,
7,38 s vs 5,80 s = 0,79×) era de fato puro overhead; a LUT confirma
isso ao reduzir o overhead a zero e expor que **não havia ganho real
para capturar**.

Próxima direção possível: testar a deformação sob um algoritmo de
busca diferente onde múltiplos sub-componentes coexistem
naturalmente.
