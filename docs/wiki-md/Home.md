# KnightMove: Wiki

Documentação de pesquisa sobre o **passeio do cavalo**: da otimização de
backtracking à reformulação no espaço de ciclos sobre GF(2), e daí aos
invariantes topológicos que essa reformulação expõe.

Esta wiki é a versão de leitura do projeto. O código está em
[KnightMove](https://github.com/matheus-fsc/KnightMove) e o visualizador
interativo em
[knight-tour-visualizer](https://github.com/matheus-fsc/knight-tour-visualizer).

## Por onde começar

Se você chegou aqui sem contexto, leia nesta ordem:

1. **[[Historia]]**: como as ideias apareceram, e por que o projeto mudou de rumo no meio.
2. **[[Teoria-GF2]]**: o que são os "loops" e por que eles são o espaço de ciclos.
3. **[[Invariante-Q]]**: o resultado central: `Q(n) = 3`.

## Índice

| página | assunto |
|---|---|
| [[Historia]] | a evolução do projeto, etapa por etapa |
| [[Teoria-GF2]] | espaço de ciclos, `beta_1`, tours como vetores binários |
| [[Invariante-Q]] | o deficit `Q(n) = 3`, sua localidade e sua prova |
| [[Algoritmo]] | backtracking + espaço de ciclos + union-find incremental |
| [[Topologia]] | cilindro, toro, Klein, toro cisalhado, e a hierarquia de `Q` |
| [[Benchmarks]] | comparações entre solvers, com baselines honestos |
| [[Estimativas]] | estimador de Knuth, `N(10)`, `N(12)`, divide-and-conquer |
| [[Resultados-Negativos]] | o que foi testado e não funcionou, e por quê |
| [[Paper-Status]] | o que o paper corrente reivindica, com os selos provado/verificado/aberto |
| [[Prior-Art]] | o que já existia na literatura |
| [[Formalizacao-Lean]] | os teoremas em Lean 4 |
| [[Guia-do-Codigo]] | mapa do repositório |

## Estado

Fase de revisão. Código e resultados estáveis. O que está em aberto é a
posição da contribuição em relação à literatura, comentários e referências são
bem-vindos via issues.
