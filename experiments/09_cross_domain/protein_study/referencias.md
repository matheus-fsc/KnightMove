# Referências anotadas

Convenção: só entram trabalhos cuja página/PDF eu **abri** nesta sessão. Quando o
dado veio de resumo de busca e não do texto integral, está marcado
**[via resumo de busca]**. O que cada item **decide** para o estudo vem em itálico.

---

## GASS e família

**Izidoro, S. C.; de Melo-Minardi, R. C.; Pappa, G. L.**
*GASS: identifying enzyme active sites with genetic algorithms.*
Bioinformatics 31(6):864–870, 2015.
DOI: https://doi.org/10.1093/bioinformatics/btu746
Página aberta: https://academic.oup.com/bioinformatics/article/31/6/864/215207

Enuncia o problema como casamento de um sítio `A₁` de `pA` numa proteína `pB` de
função desconhecida. Indivíduo do AG = **vetor** de resíduos, cada posição com nome,
cadeia, posição na sequência, **último átomo pesado (LHA)** da cadeia lateral e
coordenadas 3D. Fitness `Fit(v,w) = Σ|vᵢ − wᵢ|²` (RMSD modificada; a escolha é
justificada por sítios distintos terem RMSD tradicional parecida). Discute PINTS
(DFS em grafo de pseudoátomos), ASSAM (isomorfismo de subgrafo) e CatSId (busca de
subgrafo com poda a 1,5 Å). Registra que o **ASSAM limita o template a 12
aminoácidos**, "enquanto o GASS não tem limite". Afiliações: UNIFEI Itabira e
DCC/UFMG. Corresponding: glpappa@dcc.ufmg.br.
*Decide: a representação (vetor, sem sequência/adjacência) e o fitness (métrico
contínuo). Ambos são premissas do veredito. Também decide que o gargalo declarado é
tamanho de template e rigidez do casamento — não tempo.*

---

**Moraes, J. P. A.; Pappa, G. L.; Pires, D. E. V.; Izidoro, S. C.**
*GASS-WEB: a web server for identifying enzyme active sites based on genetic
algorithms.* Nucleic Acids Research 45(W1):W315–W319, 2017.
DOI: https://doi.org/10.1093/nar/gkx337
Página aberta: https://academic.oup.com/nar/article/45/W1/W315/3782604

Fitness = RMSD modificada entre template e resíduos casados. Limiar operacional
**≤ 5 Å**: 97,78 % dos sítios encontrados conforme o CSA têm fitness ≤ 5 Å; 51,11 %
têm fitness < 1 Å. O usuário escolhe o **tamanho do template** (número de resíduos
do sítio).
*Decide: o valor numérico do limiar de aceitação, usado no protótipo (τ = 5 Å frouxo,
2 Å estrito). Confirma que a aceitação é um corte contínuo.*

---

**Paiva, V. A.; Mendonça, M. V.; Silveira, S. A.; Ascher, D. B.; Pires, D. E. V.;
Izidoro, S. C.**
*GASS-Metal: identifying metal-binding sites on protein structures using genetic
algorithms.* Briefings in Bioinformatics 23(5):bbac178, 2022.
DOI: https://doi.org/10.1093/bib/bbac178
Página aberta: https://academic.oup.com/bib/article/23/5/bbac178/6590153

AG paralelo; indivíduo = vetor de características dos resíduos do sítio; fitness
compara **distâncias entre átomos de referência**. Banco de 928 templates de sítio
de ligação a metal, 15 tipos de íon. Sem restrição de número de resíduos.
*Decide: a família GASS mantém a mesma representação e o mesmo tipo de fitness
métrico — o veredito não é específico do artigo de 2015.*

---

**Página do projeto GASS (DCC/UFMG).**
https://homepages.dcc.ufmg.br/~glpappa/gass/ — aberta.
Autores: Sandro Carvalho Izidoro, Raquel C. de Melo-Minardi, Gisele Lobo Pappa.
Enfatiza ausência de restrição no número de aminoácidos e no tamanho do sítio ativo.
Distribui código-fonte (`GASS.tar.gz`).
*Decide: identifica os autores; confirma que "sem limite de tamanho" é a alegação
central de venda do método.*

---

## Métodos por grafo em sítios ativos (o "modelo com hierarquia de busca por grafo")

**Artymiuk, P. J.; Poirrette, A. R.; Grindley, H. M.; Rice, D. W.; Willett, P.**
*A Graph-theoretic Approach to the Identification of Three-dimensional Patterns of
Amino Acid Side-chains in Protein Structures.*
Journal of Molecular Biology 243(2):327–344, 1994.
DOI: https://doi.org/10.1006/jmbi.1994.1657 (metadados via Crossref)

Cada cadeia lateral é representada por pseudoátomos; padrões são buscados por
**isomorfismo de subgrafo** (variante do algoritmo de Ullmann) — o programa ASSAM.
*Decide: é a formulação (b) da §2.1 do relatório — grafo de correspondência e clique.
É o alvo real do mapeamento testado. [via resumo de busca para o conteúdo; DOI e
metadados bibliográficos verificados no Crossref]*

---

**Barker, J. A.; Thornton, J. M.**
*An algorithm for constraint-based structural template matching: application to 3D
templates with statistical analysis.* Bioinformatics 19(13):1644–1649, 2003.
DOI: https://doi.org/10.1093/bioinformatics/btg226 (metadados via Crossref)

Algoritmo **Jess**, casamento de template por satisfação de restrições geométricas.
Templates parciais de 3 a 8 resíduos.
*Decide: confirma que a faixa de tamanho usada na prática (3–8) é a mesma medida no
M-CSA, e que o critério é geométrico/por restrição — não combinatório-linear.
[via resumo de busca; DOI verificado no Crossref]*

---

**Stark, A.; Russell, R. B.**
*Annotation in three dimensions. PINTS: Patterns in Non-homologous Tertiary
Structures.* Nucleic Acids Research 31(13):3341–3344, 2003.
Página de ajuda aberta: https://www.russelllab.org/pints/help.shtml
DOI: https://doi.org/10.1093/nar/gkg506

Padrões de até **10 resíduos** na publicação de 2003 (o servidor atual declara
limites maiores).
*Decide: junto com o teto de 12 do ASSAM, mostra que nenhum método por grafo
documentado tem corte em 5 — o "~5" do usuário não vem daqui. [via resumo de busca]*

---

## Estatística do objeto (tamanho de sítio ativo)

**Ribeiro, A. J. M. *et al.*** *Mechanism and Catalytic Site Atlas (M-CSA): a database
of enzyme reaction mechanisms and active sites.*
Nucleic Acids Research 46(D1):D618–D623, 2018.
DOI: https://doi.org/10.1093/nar/gkx1012
Página aberta: https://academic.oup.com/nar/article/46/D1/D618/4584620
Dados baixados: https://www.ebi.ac.uk/thornton-srv/m-csa/media/flat_files/curated_data.csv
(cópia em `dados/mcsa_curated_data.csv`, 961 entradas, 28 188 linhas)

"M-CSA contains 961 entries, 423 of these with detailed mechanism information, and
538 with information on the catalytic site residues only." O artigo **não** publica a
distribuição de tamanho de sítio; eu a computei do arquivo curado
(`mcsa_tamanho_sitio.py`): média **4,62** resíduos/sítio, mediana 4, p90 = 8,
máx = 23; 69,3 % dos sítios com ≤ 5 resíduos. Excluindo `role type = spectator`:
média **3,44**, mediana 3.
*Decide: a ambiguidade "5 proteínas vs 5 resíduos" — são resíduos, e ~5 é a
estatística do objeto biológico, não um limite de algoritmo.*

---

**Bartlett, G. J.; Porter, C. T.; Borkakoti, N.; Thornton, J. M.**
*Analysis of Catalytic Residues in Enzyme Active Sites.*
Journal of Molecular Biology 324(1):105–121, 2002.
DOI: https://doi.org/10.1016/S0022-2836(02)01036-7
Metadados verificados no Crossref e na API do Semantic Scholar; **texto integral não
aberto** (paywall). Valor citado (615 resíduos catalíticos em 178 enzimas ⇒ média
3,5) vem **[via resumo de busca]** e é usado apenas como corroboração secundária do
meu 3,44 medido no M-CSA.
*Decide: nada sozinho; serve de checagem cruzada.*

---

## "mod 2" e conectividade — onde a técnica realmente funciona

**Cygan, M.; Nederlof, J.; Pilipczuk, Ma.; Pilipczuk, Mi.; van Rooij, J. M. M.;
Wojtaszczyk, J. O.** *Solving Connectivity Problems Parameterized by Treewidth in
Single Exponential Time.* FOCS 2011; versão de revista em ACM Transactions on
Algorithms, 2022.
DOI (FOCS): https://doi.org/10.1109/FOCS.2011.23 —
DOI (TALG): https://doi.org/10.1145/3506707
Páginas abertas: https://dl.acm.org/doi/10.1145/3506707 ,
https://research.google.com/pubs/pub37373.html

Técnica **Cut & Count**: algoritmos Monte Carlo `cᵗʷ·|V|^{O(1)}` para problemas de
conectividade (Hamiltonian Path, Steiner Tree, Feedback Vertex Set, Connected
Dominating Set), baseados em contagem de **cortes consistentes** — soluções
desconexas cancelam aos pares mod 2.
*Decide: o "mod 2" resolve conectividade em regime de **decisão/contagem
randomizada**, parametrizado por treewidth. O GASS precisa **enumerar e ranquear**,
`Γ` é denso e o predicado é métrico. Fecha a última rota pela qual GF(2) poderia
entrar. [via resumo de busca para o conteúdo técnico]*

---

## Topologia em bioinformática estrutural — o que NÃO é isto

**Xia, K.; Wei, G.-W.** *Persistent homology analysis of protein structure,
flexibility and folding.* arXiv:1412.2779 (2014).
URL: https://arxiv.org/abs/1412.2779

Homologia persistente construída por **filtração métrica** (raio crescente dos átomos);
features em dimensões 0, 1, 2 = componentes, ciclos, cavidades. Aplicações a
flexibilidade, enovelamento e afinidade de ligação.
*Decide: existe topologia algébrica aplicada a proteínas, mas ela é **descritiva
sobre uma filtração métrica** — extrai invariantes de uma estrutura dada. Não é
resolver/enumerar soluções em F₂. Não confundir com a ideia testada aqui; a
literatura de "GF(2)-solving" em busca de sítio ativo, até onde procurei, **não
existe**.*

---

## Contexto interno (deste repositório)

- `paper/knight_tour_tightness.tex`, §1.3 — a cadeia
  `Tours ⊆ 2-fatores ⊆ ker ∂₁ = Z₁`, e a frase que organiza tudo: o termo do meio é
  definido por condições **lineares** sobre F₂, mas a passagem para "tour" exige
  **conexidade**, que não é condição linear alguma.
- `tsp_cycle_space/RESULTADOS.md` — TSP denso: 3/91 ciclos fundamentais preservam
  hamiltonicidade; 0/30 instâncias em que o método de ciclos bate o 2-opt.
- `pathfinding_xor_experiment/RESULTADOS.md` — grades esparsas: compatibilidade
  17,9 % média (26,9 % → 11,7 % com `n` crescendo); **90,8 %** das falhas são
  desconexão.
- `path_decomposition_experiment/RESULTS.md` — caminhos `s→t` no cavalo: cobertura
  0,0086 % (4 de 46 666 no 6×6); falha dominante `bad_count` 73–84 %.
- `ratio_analysis` (memória do projeto) — razão tours/2-fatores decai
  exponencialmente: 0,272 (n=6) → 0,10 (n=12).
