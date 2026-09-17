# Agent Prompt — Viabilidade: fusão XOR de loops (GF(2)) aplicada a busca em proteínas

## 0. Regras de operação

- Trabalhe em `/home/math/Dev/knight_tour`. **Crie e use o diretório
  `protein_study/`**; não escreva nada fora dele (exceto ler os arquivos de
  contexto listados abaixo).
- Escreva em **português (pt-BR)**, inclusive nomes de arquivo, comentários e
  relatórios. É o padrão do repositório.
- Python: use `../venv/bin/python` (Python 3.14, já tem `networkx`, `numpy`).
- **Este é um estudo de viabilidade, não de implementação.** O produto final é
  um relatório com um veredito defensável. Só escreva código se for para
  produzir uma medição que sustente o veredito — e, nesse caso, mantenha o
  escopo mínimo (um protótipo de brinquedo, não uma ferramenta).
- **Honestidade acima de entusiasmo.** Se a conclusão for "não vale a pena",
  esse é um resultado bom e deve ser escrito com a mesma clareza que um
  positivo. Veja a §3 — o histórico deste repositório é majoritariamente de
  resultados negativos sobre exatamente essa ideia, e você precisa explicar por
  que o caso das proteínas seria diferente (ou admitir que não é).

---

## 1. De onde vem a ideia (contexto interno — leia antes de qualquer coisa)

Este repositório é uma pesquisa sobre a estrutura em **GF(2)** do passeio do
cavalo. O objeto central: dado um grafo `G`, o **espaço de ciclos**
`Z₁ = ker ∂₁ ⊆ F₂^E` tem dimensão `β₁ = |E| − |V| + 1`, e qualquer subgrafo par
é soma XOR de ciclos fundamentais. A ideia recorrente — a que se quer testar em
bioinformática — é:

> em vez de **buscar** cada solução do zero, tome uma solução conhecida e some
> (XOR) ciclos do espaço de ciclos para gerar outras soluções. A busca vira
> álgebra linear sobre F₂.

Leituras obrigatórias, nesta ordem:

| Arquivo | O que extrair |
|---|---|
| `paper/knight_tour_tightness.tex` | §1.1 (os objetos: `Z₁`, `β₁`, `𝒯`, `Ham`, deficit) e §1.3 (a cadeia "subgrafo par → 2-fator → tour"). **É a §1.3 que importa mais**: ela diz exatamente onde a álgebra linear para de funcionar. |
| `pathfinding_xor_experiment/RESULTADOS.md` + `README.md` | Tentativa de usar XOR para enumerar caminhos `s→t` em grades esparsas. |
| `path_decomposition_experiment/RESULTS.md` | Tentativa análoga por decomposição em caminhos. |
| `tsp_cycle_space/` | Tentativa anterior no TSP (grafo denso). |

---

## 2. O que se quer testar (a hipótese do usuário)

O professor do usuário é **proprietário do algoritmo GASS**, um algoritmo
genético para busca de padrões em proteínas. Ele comentou que existe um modelo
com **hierarquia de busca por grafo**, cujo limite prático seria da ordem de
**~5** (o usuário registrou "5 proteínas"; **muito provavelmente são ~5
resíduos/nós no padrão do sítio ativo, não 5 proteínas** — confirme isso na
literatura antes de construir qualquer coisa em cima, e registre a ambiguidade
no relatório).

A intuição do usuário, que você deve **testar e não assumir**:

> o cavalo tem 8 movimentos possíveis e cada passo depende do anterior; a busca
> por padrão em proteína também é uma cadeia de escolhas com grau de ramificação
> limitado e dependência local. Se a estrutura de ciclos em GF(2) organiza o
> espaço de soluções num caso, talvez organize no outro.

Trate isso como uma **analogia a ser falsificada**. As perguntas concretas:

1. **O problema-alvo tem espaço de ciclos não-trivial?** Busca de padrão em
   proteína é tipicamente **matching de subgrafo / clique** num grafo de
   produto, não busca de subgrafo gerador par. Se as soluções não são subgrafos
   pares, `Z₁` não age sobre elas e a ideia morre aqui — diga isso e pare.
2. **Se age: qual é o análogo de "2-fator"?** No cavalo, `∂₁c = 0` (condição
   linear) é barata e a **conexidade** é a barreira não-linear. Qual é a
   condição barata e qual é a barreira no caso proteico?
3. **A restrição é linear sobre F₂ ou é métrica?** Sítios ativos são casados com
   **tolerância geométrica** (RMSD, distâncias euclidianas entre resíduos). Uma
   restrição com limiar contínuo não é uma restrição linear sobre F₂. Este é o
   ponto mais provável de falha — investigue-o cedo, não no fim.
4. **De onde viria o ganho?** No 10×10 o ganho real medido não veio do XOR: veio
   de **pressão de vértice** (poda por grau) e de **detecção incremental de
   sub-ciclos com Union-Find**. Se o ganho candidato aqui também for "poda
   melhor", diga isso — pode ser uma contribuição real, mas então **não é a
   ideia do XOR** e deve ser reportada com o nome certo.
5. **O limite ~5 é de quê?** Combinatório (explosão do produto), de tempo, de
   memória, ou de qualidade estatística do casamento? A resposta muda
   completamente se GF(2) tem algo a oferecer.

---

## 3. Contra-evidência interna (não ignore — é o núcleo do prompt)

A mesma ideia já foi testada três vezes neste repositório e **falhou nas três**:

- **TSP (denso):** só ~3% dos XORs preservavam hamiltonicidade.
- **Pathfinding em grades esparsas:** compatibilidade ~18%, **decrescente com
  `n`**; o modo de falha dominante foi **desconexão em 90,8%** dos casos — não
  violação de grau.
- **Decomposição em caminhos (cavalo):** cobertura de **0,0086%** dos caminhos
  `s→t` (4 de 46.666 no 6×6); modo de falha dominante `bad_count` 73–84%. Um
  backtracking com poda de grau + conexidade venceu com folga.

O padrão é consistente e vale como lei operacional deste projeto:

> **A condição linear sobre F₂ é barata; a conexidade não é uma condição linear
> e é ela que domina o custo.** O XOR entrega objetos que satisfazem a parte
> barata e violam a cara.

Portanto o ônus da prova está do lado positivo. **A pergunta que organiza todo o
estudo é: existe alguma razão estrutural para que o caso proteico escape desse
padrão?** Um candidato honesto a "sim" seria: se o critério de aceitação
proteico for *local* (por resíduo/aresta) e não exigir um objeto globalmente
conexo, a barreira que matou os três casos anteriores simplesmente não existe.
Verifique se esse é o caso. Se não for, o veredito é negativo.

---

## 4. Pesquisa externa a fazer (use WebSearch/WebFetch)

Você **não** tem conhecimento confiável sobre GASS — pesquise, não recite.

- **GASS**: artigo original, autores, instituição, o que exatamente busca
  (sítio ativo? motivo estrutural?), representação do indivíduo, função de
  fitness, limites relatados. Identifique **qual autor é o professor do
  usuário** e cite os trabalhos dele com precisão.
- **Métodos por grafo em sítios ativos**: quais são, e onde está o limite de
  tamanho do padrão. Termos: *protein structure comparison graph*, *active site
  search*, *substructure matching*, *maximum common subgraph*, *product graph
  clique*, *geometric hashing*, *RMSD threshold matching*.
- **Álgebra sobre F₂ / homologia em bioinformática estrutural**: existe alguma
  literatura? Adjacente e relevante: *persistent homology* em estrutura
  proteica, *cycle basis* em redes bioquímicas. Note que homologia persistente
  é sobre **filtração métrica**, não sobre F₂-solving — não confunda as duas.
- **Cut & Count / soluções paramétricas mod 2**: onde o "mod 2" realmente ajuda
  em problemas de conectividade, e por que a técnica é de **contagem/decisão**,
  não de **enumeração** — isto é diretamente relevante ao veredito.

Anote toda citação com DOI/URL em `protein_study/referencias.md`. Não cite nada
que você não abriu.

---

## 5. Entregáveis (em `protein_study/`)

1. `README.md` — o que é o estudo, como rodar o que houver, estado atual.
2. `referencias.md` — bibliografia anotada, com o que cada trabalho decide.
3. `RELATORIO.md` — **o produto principal**, com esta estrutura:
   - Resumo do veredito em ≤5 linhas, logo no topo.
   - O que é GASS e qual é o gargalo real (com fonte).
   - O mapeamento formal proposto: o que é `V`, `E`, `Z₁`, o que é uma
     "solução", e se soluções são ou não subgrafos pares. Se o mapeamento não
     fecha, **mostre onde não fecha** — isso é o resultado.
   - Respostas às 5 perguntas da §2.
   - Confronto explícito com a §3: por que seria (ou não seria) diferente.
   - **Veredito**: viável / inviável / viável sob a condição X. Se inviável,
     aponte a alternativa mais próxima que faz sentido (frequentemente: a poda
     por pressão de vértice + Union-Find incremental, que é o que de fato
     ganhou no cavalo).
   - **O que falta para decidir**, caso reste incerteza — lista de perguntas
     concretas para o professor.
4. Só se necessário: um protótipo mínimo (`prototipo.py` + `resultados/`) sobre
   dados de brinquedo ou um punhado de PDBs, medindo *uma* quantidade que
   decida *uma* das perguntas. Registre no relatório o que ele mede e o que
   **não** mede.

---

## 6. Critérios de aceitação

- Nenhuma afirmação sobre GASS sem fonte verificada.
- A ambiguidade "5 proteínas vs 5 resíduos" explicitamente resolvida ou
  explicitamente marcada como não resolvida.
- A contra-evidência da §3 confrontada de frente, não omitida.
- O veredito é uma frase que se pode discordar. "Promissor, precisa de mais
  estudo" **não é** um veredito e será rejeitado.
- Nada fora de `protein_study/` foi modificado.
