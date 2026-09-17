# Relatório — Viabilidade da fusão XOR de loops (GF(2)) na busca de sítios ativos em proteínas

Data: 2026-08-17. Escopo: estudo de viabilidade. Diretório: `protein_study/`.

---

## Veredito (resumo em 5 linhas)

**Inviável.** As soluções do problema-alvo (casamento de template de sítio ativo,
estilo GASS) são **k-cliques transversais num grafo de correspondência**, sujeitas a
um **limiar métrico contínuo** (RMSD) — **não são subgrafos pares**, logo
`Z₁ = ker ∂₁` não age sobre elas: para `k` par a própria solução nem pertence a `Z₁`,
e para `k` ímpar o conjunto de soluções não é coset de subespaço algum.
Medido em dados reais (1PPF→1ACB, 1PPF→1LYZ): **0 de 17.367** XORs de pares de
soluções válidas produz uma solução válida, embora **10.216 deles tenham grau par**
— a condição linear é satisfeita de graça e não compra nada.
Pior para a hipótese: no caso proteico **não existe a barreira global cara** que o XOR
prometia contornar (razão global/local medida entre 0,68 e 1,00, contra decaimento
exponencial no cavalo). Não há nem o mecanismo, nem o problema que ele resolveria.

---

## 1. O que é o GASS e qual é o gargalo real

**GASS (Genetic Active Site Search)** — Izidoro, de Melo-Minardi e Pappa,
*Bioinformatics* 31(6):864–870, 2015, DOI [10.1093/bioinformatics/btu746](https://doi.org/10.1093/bioinformatics/btu746).
Afiliações: UNIFEI (campus avançado de Itabira) e DCC/UFMG. Desdobramentos:
**GASS-WEB** (Moraes, Pappa, Pires, Izidoro, *NAR* 45(W1):W315–W319, 2017) e
**GASS-Metal** (Paiva *et al.*, *Brief. Bioinform.* 23(5):bbac178, 2022).
O professor do usuário é, com alta probabilidade, **Sandro C. Izidoro** (UNIFEI) —
autor principal do GASS e coautor das duas extensões.

**Problema, como o próprio artigo o enuncia:** dado um conjunto de *N* aminoácidos que
compõem o sítio ativo `A₁` de uma proteína `pA` de função conhecida, e uma segunda
proteína `pB` de função desconhecida e tamanho *M*, buscar um casamento de `A₁` em `pB`.

**Representação e fitness.** Um indivíduo do AG é um **vetor** de *N* posições; cada
posição carrega nome do aminoácido, cadeia, posição na sequência, o **último átomo
pesado (LHA)** da cadeia lateral e as coordenadas 3D. O fitness é
`Fit(v,w) = Σ |vᵢ − wᵢ|²` — uma RMSD modificada, escolhida porque sítios distintos
podem ter RMSD tradicional parecida. Na prática o limiar operacional é **≤ 5 Å**
(no GASS-WEB, 97,78% dos sítios recuperados conforme o CSA têm fitness ≤ 5 Å).

**Gargalo real.** O artigo do GASS posiciona o método contra métodos por grafo
(**PINTS**, **ASSAM**, **CatSId**) exatamente pelas restrições destes:
casamento *exato* de resíduos e **limite no tamanho do template** — o GASS registra
que o ASSAM limita o template de referência a **12 aminoácidos**, "enquanto o GASS
não tem limite". Ou seja: o gargalo que o GASS ataca é **combinatório-estrutural**
(tamanho do template e rigidez do casamento), não numérico.

### A ambiguidade "5 proteínas vs 5 resíduos" — **resolvida: são resíduos**

Baixei os dados curados do M-CSA (`dados/mcsa_curated_data.csv`, 961 entradas) e medi
diretamente o número de resíduos por sítio (`mcsa_tamanho_sitio.py`):

| conjunto | entradas | média | mediana | p90 | máx | ≤5 resíduos |
|---|---:|---:|---:|---:|---:|---:|
| todos os resíduos anotados | 961 | **4,62** | 4 | 8 | 23 | **69,3 %** |
| excluindo `role type = spectator` | 888 | **3,44** | 3 | 7 | 21 | **82,4 %** |

O valor 3,44 reproduz de forma independente o ~3,5 resíduos catalíticos por enzima
de Bartlett *et al.* (2002). **"~5" é o tamanho típico de um sítio ativo**, não um
número de proteínas — e não é um limite de algoritmo: é a estatística do objeto
biológico. Nenhuma fonte que abri atribui a métodos por grafo um teto de 5; o teto
citado é **12** (ASSAM) e **10** (PINTS, na publicação de 2003).

**Consequência para o estudo:** o "limite prático ~5" que motivou a pergunta
provavelmente **não é um limite computacional a ser vencido**. Se o professor tinha
em mente outra coisa (p. ex. um método específico com corte em 5), isso precisa ser
confirmado — ver §6.

---

## 2. O mapeamento formal — e onde ele não fecha

### 2.1 O que se poderia querer que fosse `V`, `E`, `Z₁`

Há duas escolhas naturais, e nenhuma funciona.

**(a) Grafo da proteína.** `V` = resíduos, `E` = pares dentro de um cutoff.
Aqui `Z₁` existe e é enorme, mas uma **solução não é um subconjunto de arestas**:
é uma **injeção** `φ: {1..k} → V` do template nos resíduos. Objetos de tipo diferente.
`Z₁ ⊆ F₂^E` simplesmente não atua sobre o conjunto de injeções.

**(b) Grafo de correspondência (*product graph*).** Esta é a formulação canônica dos
métodos por grafo (Artymiuk *et al.* 1994 / ASSAM; clique detection):

* `V(Γ) = { (i, a) : i` índice do template, `a` resíduo do alvo, tipos compatíveis `}`
* `{(i,a),(j,b)} ∈ E(Γ)` sse `i ≠ j`, `a ≠ b` e `| d_T(i,j) − d_P(a,b) | ≤ ε`
* **solução** = **k-clique transversal** (um nó por índice `i`) que ainda passe no
  teste global `RMSD(template, candidato) ≤ τ`.

É em (b) que a analogia com o cavalo é testável, e é onde ela quebra.

### 2.2 Onde não fecha — três obstruções independentes

**(i) A solução não é um subgrafo par (metade dos casos nem está em `Z₁`).**
Todo vértice de uma `k`-clique tem grau `k−1`. Logo
`1_clique ∈ Z₁(Γ) ⟺ k` é **ímpar**. Para `k` par — e `k` par cobre 51,8 % dos sítios
do M-CSA (498 de 961) — o vetor indicador da solução **não pertence sequer ao espaço de ciclos**,
e `∂₁ · 1_solução ≠ 0`. Não existe o análogo do "2-fator": no cavalo há uma cadeia
`Tours ⊆ 2-fatores ⊆ Z₁` em que o termo do meio é definido por condições **lineares**
sobre F₂; aqui **o termo do meio não existe**.

**(ii) Mesmo para `k` ímpar, o conjunto de soluções não é coset de subespaço.**
Toda solução tem exatamente `C(k,2)` arestas. Para `A ⊕ B` ter `C(k,2)` arestas é
preciso `|A ∩ B| = C(k,2)/2`, e além disso o resultado teria de ser de novo uma clique
transversal — duas condições que nada garante. Isto é uma afirmação sobre a
**cardinalidade e a forma** do suporte, não sobre paridade; nenhuma quantidade de
álgebra sobre F₂ a torna linear.

**(iii) O critério de aceitação é um limiar métrico contínuo.**
`|d_T − d_P| ≤ ε` e `RMSD ≤ τ` são desigualdades sobre reais. Não são o núcleo de
nenhuma aplicação linear sobre F₂. Este era, como o prompt antecipou, o ponto mais
provável de falha — e é o mais estrutural dos três: mesmo que (i) e (ii) fossem
contornados por uma reformulação esperta, o predicado "é um sítio" não é uma
condição mod 2.

### 2.3 A medição que fecha o caso

`prototipo.py` constrói `Γ` de verdade a partir de estruturas do PDB (LHA por resíduo,
como o GASS), enumera **exaustivamente** as soluções por backtracking, e testa o
fechamento sob XOR. Template: tríade catalítica His57/Asp102/Ser195 da elastase
leucocitária humana (1PPF, cadeia E), estendida pelos resíduos mais próximos do
centroide para `k = 4..7`. Alvos: quimotripsina (1ACB, cadeia E — homóloga) e
lisozima (1LYZ — controle negativo). `ε = 1,5 Å`, `τ = 2 Å` (e 5 Å).

Modo `geometrico` (sem restrição de tipo de aminoácido — o regime com muitas soluções,
que é o único em que o teste de XOR tem poder estatístico):

| alvo | k | \|V(Γ)\| | \|E(Γ)\| | β₁(Γ) | espaço bruto | cliques par-a-par | RMSD ≤ 2 Å | k-clique ∈ Z₁ |
|---|--:|--:|--:|--:|--:|--:|--:|:--:|
| 1ACB_E | 3 | 723 | 2 610 | 1 904 | 1,4·10⁷ | 344 | 344 | sim |
| 1ACB_E | 4 | 964 | 5 242 | 4 303 | 3,4·10⁹ | 190 | 178 | **não** |
| 1ACB_E | 5 | 1 205 | 9 772 | 8 580 | 8,1·10¹¹ | 31 | 21 | sim |
| 1ACB_E | 6 | 1 446 | 17 064 | 15 628 | 2,0·10¹⁴ | 5 | 4 | **não** |
| 1ACB_E | 7 | 1 687 | 27 278 | 25 600 | 4,7·10¹⁶ | 1 | 1 | sim |
| 1LYZ_A | 3 | 387 | 1 188 | 812 | 2,1·10⁶ | 155 | 155 | sim |
| 1LYZ_A | 4 | 516 | 2 348 | 1 846 | 2,8·10⁸ | 73 | 66 | **não** |

**Teste de fechamento sob XOR, agregado sobre todas as execuções:**

```
pares de soluções testados : 17 367
XOR que é solução válida   :      0     (0,000 %)
XOR com todos os graus pares: 10 216     (58,8 %)
```

Os 10.216 são exatamente os casos de `k` ímpar, onde ambas as soluções estão em `Z₁`
e portanto a soma também está — **e ainda assim nenhuma é solução**. Esta é a lei
operacional do repositório na sua forma mais nua: *a condição linear sobre F₂ é
barata; ela é satisfeita de graça; e sozinha não vale nada.*

---

## 3. Respostas às cinco perguntas

**1. O problema-alvo tem espaço de ciclos não-trivial?**
O grafo `Γ` tem `β₁` enorme (até 25.600 no teste). Mas **`Z₁` não age sobre o conjunto
de soluções**: as soluções são k-cliques transversais, não subgrafos pares — e para
`k` par não estão nem em `Z₁`. Conforme o critério da §2 do prompt, **a ideia morre
aqui.** As respostas seguintes são registradas porque respondem a perguntas
independentes do usuário, não porque a ideia sobreviva.

**2. Qual é o análogo de "2-fator"?**
**Não há.** No cavalo: `∂₁c = 0` (barato, linear) → conexidade (caro, não-linear).
No caso proteico: o filtro **par-a-par** `|d_T(i,j) − d_P(a,b)| ≤ ε` é o barato, e o
**RMSD global por superposição** é o caro. Mas o barato aqui **não é linear** — é
métrico —, então não existe o subespaço intermediário que dá sentido ao XOR.

**3. A restrição é linear sobre F₂ ou é métrica?**
**Métrica**, sem ambiguidade: `ε = 1,5 Å` par-a-par e fitness/RMSD ≤ 5 Å global,
ambos limiares contínuos. Investigado cedo, como o prompt pedia; é a obstrução mais
fundamental porque sobrevive a qualquer reformulação do grafo.

**4. De onde viria o ganho?**
De **poda**, não de XOR — exatamente como no 10×10. Medido: o espaço bruto do produto
vai de 1,4·10⁷ (k=3) a 4,7·10¹⁶ (k=7), e o filtro par-a-par o colapsa para
344 → 190 → 31 → 5 → 1 candidatos, em **menos de 0,4 s por proteína**.
E há um agravante para a hipótese: **a barreira global é fraca**. A razão
(soluções que passam no RMSD)/(cliques par-a-par) medida foi
1,00 / 0,94 / 0,68 / 0,80 / 1,00 (1ACB) e 1,00 / 0,90 / 1,00 / 0,50 (1LYZ) —
isto é, o filtro local **já quase determina** o resultado global. Compare com o
cavalo, onde a razão tours/2-fatores decai exponencialmente (0,272 em n=6 → 0,10 em
n=12; ver `ratio_analysis`). **Não há um gargalo global caro para o XOR atacar.**

**5. O limite ~5 é de quê?**
De **estatística do objeto biológico**, não de algoritmo: 4,62 resíduos por sítio em
média no M-CSA (mediana 4; 69,3 % dos sítios com ≤ 5). Nos meus testes, `k` maior
**barateia** a busca em vez de encarecê-la (mais restrições ⇒ menos candidatos:
344 → 1). O que degrada com `k` grande não é o tempo, é a **disponibilidade de
templates** e a **cobertura** — e, do lado do usuário, a chance de o casamento ser
espúrio quando `k` é pequeno (em k=3 obtive 155 casamentos "válidos" na lisozima,
que não tem sítio de serino-protease: k=3 puramente geométrico não discrimina nada).
**Ressalva de escopo:** meço o custo de **uma** proteína alvo. O GASS varre bancos de
estruturas; nada aqui mede esse custo.

---

## 4. Confronto com a contra-evidência interna (§3 do prompt)

O prompt propõe o único "sim" honesto possível: *se o critério de aceitação proteico
for local (por resíduo/aresta) e não exigir objeto globalmente conexo, a barreira que
matou os três casos anteriores não existe.*

**A premissa é verdadeira e a conclusão não segue.** O critério proteico é, de fato,
quase local (medido: razão global/local entre 0,68 e 1,00, §3 pergunta 4). Mas isso
**não salva a ideia — remove a motivação dela**:

| | TSP | Grades | Caminhos (cavalo) | **Proteínas** |
|---|---|---|---|---|
| solução ∈ `Z₁` (ou coset)? | sim | sim | sim | **não** (k par: nem em `Z₁`) |
| XOR preserva a parte barata? | sim | sim | sim | **irrelevante** |
| taxa de compatibilidade | ~3 % | 12–27 % | 0,1–0,7 % | **0 % (0/17 367)** |
| modo de falha dominante | desconexão | desconexão 90,8 % | `bad_count` 73–84 % | **cardinalidade/forma do suporte** |
| barreira global cara existe? | sim | sim | sim | **não** |

Nos três casos anteriores o XOR ao menos **gerava candidatos do tipo certo** e falhava
na conectividade. No caso proteico ele nem chega lá: **o XOR de duas soluções não é
um objeto do tipo "solução"**, e o modo de falha (58,8 % dos XORs têm grau par e
mesmo assim 0 % são solução) mostra que a paridade — a única coisa que o XOR preserva
— é ortogonal ao que define uma solução. É o pior dos dois mundos para a hipótese:
**não há o mecanismo, e não há o problema que ele resolveria.**

Sobre a analogia do usuário ("o cavalo tem 8 movimentos e cada passo depende do
anterior"): ela falha num ponto anterior a tudo isto. No cavalo, a **sequência** é
parte da solução — um tour é um caminho. Num template de sítio ativo **não há
sequência**: o indivíduo do GASS é um conjunto de `k` resíduos avaliado por uma soma
sobre **todos os pares**. Não é uma cadeia de escolhas com ramificação limitada; é um
casamento simultâneo. A ramificação existe no *algoritmo de busca*, não no *objeto*.

### Sobre "mod 2" na literatura de conectividade

A técnica **Cut & Count** (Cygan, Nederlof, Pilipczuk, Pilipczuk, van Rooij,
Wojtaszczyk, FOCS 2011 / *TALG* 2022) é o lugar onde "mod 2" realmente resolve
conectividade: ela conta cortes consistentes, e as soluções desconexas se cancelam
**aos pares** mod 2. Três razões pelas quais isso não socorre o caso proteico:
(1) é **decisão/contagem**, Monte Carlo — não **enumeração**, que é o que o GASS
precisa (ranquear candidatos); (2) o parâmetro é **treewidth**, e `Γ` é denso
(17 mil arestas para 1,4 mil nós em k=6); (3) o predicado a decidir teria de ser
combinatório, e o do GASS é métrico. Homologia persistente em estrutura proteica
(Xia & Wei 2014) é sobre **filtração métrica**, não sobre resolver sistemas em F₂ —
não confundir com a ideia aqui testada.

---

## 5. Veredito

> **Inviável.** A fusão XOR de loops em GF(2) não se aplica à busca de sítios ativos
> no formato do GASS, porque as soluções são k-cliques transversais num grafo de
> correspondência sob limiar métrico contínuo — não subgrafos pares — de modo que
> `Z₁ = ker ∂₁` não atua sobre o conjunto de soluções (0 de 17.367 XORs de soluções
> reais é solução, com 58,8 % deles satisfazendo a paridade). E, mesmo que atuasse,
> não haveria ganho: o caso proteico não tem a barreira global cara que o XOR
> prometia contornar.

Isto é uma frase da qual se pode discordar: para refutá-la basta exibir uma
formulação do problema de sítio ativo em que o conjunto de soluções seja um coset
de um subespaço de `F₂^E` (§6, pergunta 4).

### A alternativa mais próxima que faz sentido

A mesma que de fato ganhou no cavalo: **poda por pressão + detecção incremental**.
Concretamente, e no vocabulário do GASS:

1. **Busca exata em `Γ` com poda transversal**, no lugar (ou como semente) do AG.
   No protótipo isto enumerou **todas** as soluções em < 0,4 s por proteína, para
   k = 3..7, sem heurística. O AG do GASS é estocástico e sem garantia de completude;
   um enumerador exato por proteína pode ser mais rápido *e* completo — desde que
   `Γ` continue esparso, o que precisa ser medido em escala de banco.
2. **Ordenação por "pressão"**: escolher primeiro o índice `i` do template com menos
   nós candidatos em `Γ` (análogo direto da pressão de vértice), e propagar a poda
   par-a-par a cada escolha. É a única transferência técnica honesta deste projeto
   para o do professor.
3. **Detecção incremental** (Union-Find) **não** transfere: ela serve para detectar
   sub-ciclos, e aqui não há requisito de conectividade.

Nomeando corretamente: isso é **poda melhor**, não a ideia do XOR. Como o prompt
antecipou, pode ser uma contribuição real — mas deve ser apresentada com esse nome,
e checada contra ASSAM/Jess/CatSId, que já fazem busca exata em grafo (a novidade do
GASS é ausência de limite de tamanho e mutações conservativas, não a busca em si).

---

## 6. O que falta para decidir — perguntas para o professor

1. **O "~5" é o quê, exatamente?** Meu dado diz que 4,62 é o tamanho médio de um
   sítio no M-CSA e que ASSAM/PINTS declaram tetos de 12 e 10 resíduos. Existe um
   método (ou uma configuração do GASS) com corte real em 5? Se sim, é corte de
   tempo, de memória, ou de qualidade estatística do casamento?
2. **Qual "modelo com hierarquia de busca por grafo"** ele tinha em mente? Se for
   ASSAM/PINTS/CatSId, este relatório já o cobre; se for outro, o mapeamento da §2
   precisa ser refeito para ele.
3. **Onde dói de verdade no GASS hoje** — tempo por proteína, varredura do banco
   inteiro, taxa de falso-positivo, ou construção de templates? Meu protótipo mede
   só o primeiro, e ele não parece ser o gargalo.
4. **Existe alguma formulação em que a solução seja um subgrafo do alvo** (e não uma
   injeção), com aceitação por predicado combinatório em vez de limiar métrico?
   É a única porta pela qual a ideia GF(2) poderia voltar. Não conheço nenhuma.
5. **Mutações conservativas**: o protótipo testou só casamento de tipo exato e
   casamento puramente geométrico. O regime real do GASS fica entre os dois, e é
   nele que `|V(Γ)|` cresce. Vale medir com a matriz de substituição real do GASS.

---

## 7. Limites deste estudo (o que o protótipo NÃO mede)

- Um alvo por vez (241 e 129 resíduos). **Não** mede varredura de banco de estruturas,
  que é o caso de uso real do GASS/GASS-WEB.
- Sem substituições conservativas; só tipo exato e geométrico puro (§6.5).
- Um único template de origem (tríade de serino-protease de 1PPF), estendido
  geometricamente para k > 3 — os resíduos extras não são catalíticos reais.
- Fitness: usei RMSD de Kabsch, não a soma `Σ|vᵢ−wᵢ|²` do GASS. Para as conclusões
  sobre estrutura algébrica isso é irrelevante (ambas são métricas contínuas); para
  comparações de sensibilidade com o GASS, não.
- O argumento central da §2.2 é **estrutural, não estatístico**: ele não depende dos
  números. As medições servem para mostrar que a estrutura se manifesta em dados
  reais, não para estabelecê-la.
