# transfer_matrix/: espectro do passeio do cavalo (2-fatores)

## Hipoteses a testar

  (i)  o autovetor dominante v_1 da transfer matrix coincide com
       f_inf(L) medido empiricamente em `incremental_subtour/`
  (ii) o gap espectral Delta(n) = lambda_1(n) - lambda_2(n) eh
       constante em n  (=> grafo de configuracoes eh expansor)
  (iii) o numero de 2-fatores 10x10 eh estimavel via potencia de lambda_1

## Resumo executivo (n=6)

  - T_bulk(6) construida: K=33.346 estados alcancaveis, nnz=1.966.428
  - count via T = 36.236 (= #2-fatores 6x6); validado por
    backtracking aresta-a-aresta independente
  - lambda_1 = 70.476614 (real, positivo)
  - lambda_2 = -39.086729 (real, NEGATIVO -- evidencia de
    estrutura bipartida no espaco de estados)
  - gap Delta = |lambda_1| - |lambda_2| = 31.39
  - razao rho = |lambda_2|/|lambda_1| = 0.555 (longe de 1: bom mixing)
  - marginal exata de 2-fatores (DP finito) bate com f_inf empirico
    de tours em ~5%; mas autovetor (bulk) diverge ~30%, sugerindo
    que n=6 eh estreito demais para o regime bulk

## Hipotese (i): autovetor vs f_inf

| comparacao                                    | L=0   | L=1   |
|-----------------------------------------------|-------|-------|
| f_inf empirico (tours, incremental_subtour)   | 0.584 | 0.136 |
| marginal exata 2-fatores (DP finito, n=6)     | 0.589 | 0.127 |
| marginal autovetor v_1 (todas arestas)        | 0.387 | 0.202 |
| marginal autovetor v_1 (so colunas bulk)      | 0.421 | 0.213 |

**Conclusoes:**
1. **2-fatores ~ tours na marginal de aresta** (5% de diferenca).
   Validacao tecnica: usar 2-fatores como veiculo eh adequado.
2. **Autovetor diverge ~30% de f_inf no n=6.** Provavel causa:
   n=6 eh estreito (6 colunas) demais para o regime bulk; as 4
   arestas "bulk" sao poucas para a media estabilizar. Para teste
   limpo seria necessario n>=10.

## Hipotese (ii): gap espectral

  Apenas n=6 disponivel ate agora. Razao |lambda_2|/|lambda_1| = 0.555.
  Para validar constancia precisariamos n=8 e n=10.

## Hipotese (iii): contagem por potencia

  3 calibracoes possiveis usando n=6 (count=36.236, lambda_1=70.477):

  | modelo                    | log C            | valor C          |
  |---------------------------|------------------|------------------|
  | count ~ C * lambda^(n-1)  | -10.78           | 2.08e-05         |
  | count ~ C * lambda^n      | -15.03           | 2.96e-07         |
  | count ~ C * lambda^(n^2)  | -142.69          | 1.07e-62         |

  O modelo `lambda^(n^2)` da especificacao original assume crescimento
  super-exponencial em n^2; mas como o tamanho da matriz T (e logo
  lambda_1) ja cresce exponencialmente com n, o modelo correto
  assintotico eh `count ~ C * lambda_1(n)^(n-1)` com lambda_1 vivendo
  num espaco que escala com n. Sem n>=8 nao da' para distinguir.

## Formulacao matematica

### Modelo

Estado: `s = (deg_c, deg_{c+1})` com `deg_*[r] in {0,1,2}`.
Passo c decide arestas com extremo esquerdo em col c:
  - (c, c+1) via knight (col_diff=1, row_diff=+-2)
  - (c, c+2) via knight (col_diff=2, row_diff=+-1)

Restricao: `deg_c[r] + #(E_c em (r,c)) = 2` para todo r
(col c fica fechada a' grau 2 apos o passo c).

Numero de 2-fatores:

  `count(n) = <s_0 | T_bulk^{n-2} . T_nocp2 | s_target>`

onde:
  - `s_0 = ((0,...,0), (0,...,0))`
  - `s_target = ((2,...,2), (0,...,0))`
  - T_nocp2 = mesma regra de T_bulk mas proibindo arestas a' col c+2
    (necessario no passo c=n-2 porque col n nao existe)

### Espaco bruto

  | n  | 3^n . 3^n   | nota                       |
  |----|-------------|----------------------------|
  | 6  | 530k        | factivel pura-Python (10s)|
  | 8  | 43M         | precisa otimizacao         |
  | 10 | 3.5B        | provavelmente intratavel ate em C |

## Verificacao independente

  count(4) = 1     (brute-force agree, transfer agree)
  count(5) = 0     (bipartido: 25 vertices, sem 2-fator)
  count(6) = 36236 (brute-force agree, transfer agree, DP agree)

## Status do escalonamento

  - **n=6** completo
  - **n=8**: build pura-Python travou em ~7GB RAM apos ~15min,
    nao finalizou. Bottleneck: BFS exploratoria com tuple-hashing
    em ~milhoes de estados. Otimizacoes pendentes:
    - codificacao int + dict-int (parcialmente feito em
      `transfer_build_fast.py`)
    - numba JIT na funcao de transicao
    - reescrita em Cython/C com hash table de capacidade fixa
    - alternativa: power iteration sem materializar T (eigenvalue
      dominante so precisa de T@v, que pode ser computado on the fly)
  - **n=10**: nao tentado; provavelmente requer >32GB RAM
    e implementacao em C/Rust

## Arquivos

  - `transfer_build.py`: versao recursiva (n<=6 testada)
  - `transfer_build_fast.py`: versao com precompute (~2x mais rapida em n=6)
  - `verify_2factor_count.py`: backtracking aresta independente
  - `spectral_analysis.py`: eigs direito e esquerdo
  - `compare_finf.py`: marginais autovetor vs f_inf
  - `exact_2factor_marginals.py`: marginais EXATAS via DP forward/backward
  - `tour_count_estimate.py`: 3 modelos de calibracao
  - `make_plots.py`: gera os 3 plots
  - `data/transfer_n6_*.npz`: T_bulk, T_nocp2
  - `data/eigenvalues_n6.json`
  - `data/eigenvector_dominant_n6.npz`: v (direita) e w (esquerda)
  - `data/eigenvector_marginals_n6.json`
  - `data/exact_2factor_marginals_n6.json`
  - `data/state_count.json`
  - `data/tour_count_estimates.json`
  - `data/plots/{spectral_gap,eigenvector_vs_finf,tour_count_scaling}.png`

## Conclusoes

1. **Hipotese (i) NAO confirmada em n=6.** O autovetor dominante
   da' marginais com erro relativo ~30% em relacao a f_inf. Mas
   isto pode ser efeito de borda (n=6 tem so 6 colunas; quase tudo
   eh boundary). As marginais EXATAS de 2-fatores (DP finito) batem
   com f_inf em 5%, indicando que a discrepancia eh do regime bulk
   vs finito, nao do conceito.

2. **Hipotese (ii) inconclusiva** (so um n).

3. **Hipotese (iii) inconclusiva** (so um n). Tres modelos compativeis
   com n=6 sozinho.

4. **Item bonus confirmado**: marginal de aresta nos 2-fatores aproxima
   muito bem a marginal de aresta nos tours (5% de diferenca). Isso
   justifica usar 2-fatores como veiculo computacional, ainda que
   tours sejam o objeto de interesse final.

## Proximos passos

Para validar (i), (ii), (iii) seria necessario completar n=8 (e idealmente
n=10). Isso requer reescrever transfer_build em linguagem compilada
ou usar numba JIT. Estimativa de esforco: alguns dias em C/Rust,
ou ~1 dia em numba se compatibilidade com Python 3.14 permitir.
