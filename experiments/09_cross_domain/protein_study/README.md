# protein_study — viabilidade de GF(2)/XOR na busca de sítios ativos (GASS)

**Estado: concluído. Veredito: inviável.** Ver [`RELATORIO.md`](RELATORIO.md).

## O que é este estudo

Testa se a ideia recorrente deste repositório — *tomar uma solução conhecida e somar
(XOR) ciclos do espaço de ciclos `Z₁ = ker ∂₁ ⊆ F₂^E` para gerar outras soluções* —
se transfere para a busca de padrões de sítio ativo em proteínas, no formato do
algoritmo **GASS** (Izidoro, de Melo-Minardi & Pappa, 2015).

Não se transfere. Uma solução do problema proteico é uma **k-clique transversal** num
grafo de correspondência, aceita por **limiar métrico contínuo** (RMSD) — não um
subgrafo par. `Z₁` não age sobre o conjunto de soluções. Medido em estruturas reais:
**0 de 17.367** XORs de pares de soluções válidas é uma solução válida, embora
**58,8 %** deles tenham todos os graus pares.

## Arquivos

| arquivo | o que é |
|---|---|
| `RELATORIO.md` | **produto principal** — veredito, mapeamento formal, respostas às 5 perguntas, confronto com a contra-evidência interna, o que falta decidir |
| `referencias.md` | bibliografia anotada com DOI/URL; marca o que veio de texto integral e o que veio de resumo de busca |
| `prototipo.py` | constrói o grafo de correspondência `Γ` a partir de PDBs reais, enumera exaustivamente os casamentos e testa fechamento sob XOR |
| `mcsa_tamanho_sitio.py` | mede a distribuição do nº de resíduos por sítio no M-CSA — resolve a ambiguidade "5 proteínas vs 5 resíduos" |
| `dados/` | 1PPF, 1ACB, 1LYZ (RCSB) e `mcsa_curated_data.csv` (EBI/M-CSA) |
| `resultados/` | `medicoes.json`, `mcsa_tamanho_sitio.json` |

## Como rodar

```sh
cd protein_study
../venv/bin/python mcsa_tamanho_sitio.py   # ~2 s
../venv/bin/python prototipo.py            # ~2 s
```

Os dados já estão em `dados/`. Para rebaixá-los:

```sh
cd dados
for id in 1PPF 1ACB 1LYZ; do curl -sSL -o $id.pdb https://files.rcsb.org/download/$id.pdb; done
curl -sSL -o mcsa_curated_data.csv \
  https://www.ebi.ac.uk/thornton-srv/m-csa/media/flat_files/curated_data.csv
```

## O que o protótipo mede — e o que não mede

**Mede:** tamanho de `Γ` e do espaço bruto em função de `k`; quantos candidatos passam
no filtro par-a-par (local) e quantos sobrevivem ao RMSD global; e se o conjunto de
soluções é fechado sob XOR de conjuntos de arestas.

**Não mede:** desempenho do GASS; varredura de banco de estruturas (mede um alvo por
vez); substituições conservativas de aminoácido; sítios inter-cadeia; qualidade
biológica dos casamentos.

O argumento central do relatório (§2.2) é **estrutural** e não depende destes números;
as medições mostram que a estrutura se manifesta em dados reais.
