# Agent Prompt — Wiki do Projeto Knight Tour

## CONTEXTO

Projeto: `knight_tour/`
O projeto tem um site existente com visualizador interativo.

**Objetivo:** transformar a documentação técnica dispersa em READMEs
em uma wiki navegável e didática, integrada ao site existente.

**Público-alvo:** estudante de graduação em computação ou matemática
que conhece teoria de grafos básica mas não conhece o projeto.

**Tom:** didático, honesto sobre o que está provado vs conjecturado,
com números concretos em toda página.

---

## REGRA ABSOLUTA

O agente deve executar as tarefas **em ordem**.
Cada tarefa marcada com `PARAR` exige confirmação antes de continuar.
Não criar arquivos de wiki antes da **Tarefa 0 ser aprovada**.

---

## TAREFA 0 — AUDITORIA COMPLETA

**Executar antes de qualquer outra coisa.**

### 0.1 — Mapear o site existente

Localizar e inspecionar:
- Diretório do site (provavelmente `knight_tour/site/` ou similar)
- Tecnologia usada (HTML puro? React? Jekyll? MkDocs? outro?)
- O visualizador: quais arquivos o compõem e o que ele faz
- Como o site é servido (gh-pages? local? outro?)

### 0.2 — Mapear todos os READMEs e `.md` do projeto

Ler **todos** os arquivos `.md` encontrados. Para cada um extrair:
- Resultado principal (número concreto quando existir)
- Status: `completo` / `em andamento` / `abandonado`
- Dependências de outros módulos

Caminhos esperados (verificar se existem):
```
knight_tour/README.md
knight_tour/board_10x10/README.md
knight_tour/forbidden_cycles_6x6/README.md
knight_tour/deficit_theorem/README.md
knight_tour/incremental_subtour/README.md
knight_tour/local_phase_heuristic/README.md
knight_tour/martingale_analysis/README.md
knight_tour/transfer_matrix/README.md
knight_tour/complex_orbit/README.md
knight_tour/6x6_higher_order/README.md
knight_tour/residual_search/README.md
knight_tour/residual_search_10x10/README.md
knight_tour/xor_clauses_benchmark/README.md
knight_tour/patch_lut_results.md
```

### 0.3 — Mapear plots e dados

Listar todos os `.png` em `data/plots/` de cada subdiretório.
Listar JSONs com resultados finais (não dados brutos).
Anotar o caminho relativo de cada um — serão usados na wiki.

### 0.4 — Mapear entry points de código

Identificar:
- `knight_tours.py` standalone (existe? onde?)
- `tour_count_estimator.py` (existe? onde?)
- Qualquer pacote Python instalável
- Como um usuário novo rodaria o código do zero

### 0.5 — Gerar `wiki_audit.md` e PARAR

Criar `knight_tour/wiki_audit.md` com:

```markdown
# Auditoria do Projeto Knight Tour

## Site existente
- Tecnologia: ...
- Arquivos: (listar todos)
- Visualizador: o que faz exatamente?
- Como servir: ...

## READMEs encontrados (N total)

| Módulo | Arquivo | Status | Resultado principal |
|--------|---------|--------|---------------------|
| ...    | ...     | ...    | ...                 |

## Plots disponíveis

| Arquivo | Localização | Assunto |
|---------|-------------|---------|
| ...     | ...         | ...     |

## Estrutura de wiki proposta

(baseada no que foi lido — não inventar)

## Pontos de atenção

(conflitos, tecnologias incompatíveis, dependências)
```

**→ PARAR. Reportar o conteúdo de `wiki_audit.md`.**
**→ NÃO criar nenhum arquivo de wiki antes da aprovação.**

---

## TAREFA 1 — ESTRUTURA E TECNOLOGIA

Após auditoria aprovada.

### 1.1 — Decidir tecnologia

Com base no que o site existente usa:

| Site usa | Wiki usa |
|----------|----------|
| HTML puro | `.html` com sidebar em `wiki/` |
| React/Next.js | componentes novos no mesmo projeto |
| Jekyll/Hugo | `.md` com front matter adequado |
| MkDocs | `mkdocs.yml` + `.md` |
| Outro | adaptar — reportar antes |

**Se HTML puro:** criar `wiki/` com arquivos `.html` independentes.
**Nunca modificar o visualizador existente.**

### 1.2 — Criar esqueleto de arquivos

Para HTML puro, criar os arquivos vazios:

```
wiki/
  index.html          ← página inicial da wiki
  style.css           ← CSS da wiki (prefixo .w- para não conflitar)
  nav.js              ← sidebar e navegação
  01-introducao.html
  02-teoria-gf2.html
  03-teorema-q3.html
  04-algoritmo.html
  05-fase-local.html
  06-benchmark.html
  07-estimativas.html
  08-topologia.html
  09-extensoes.html
  10-resultados.html
  11-codigo.html
```

### 1.3 — Template base

Cada `.html` deve incluir:

```html
<!DOCTYPE html>
<html lang="pt-BR">
<head>
  <meta charset="UTF-8">
  <title>[Título da Página] — Knight Tour Wiki</title>
  <link rel="stylesheet" href="style.css">
  <!-- MathJax para fórmulas LaTeX -->
  <script src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-chtml.js"></script>
  <!-- highlight.js para blocos de código -->
  <link rel="stylesheet"
    href="https://cdnjs.cloudflare.com/ajax/libs/highlight.js/11.9.0/styles/github.min.css">
  <script
    src="https://cdnjs.cloudflare.com/ajax/libs/highlight.js/11.9.0/highlight.min.js"></script>
  <script>hljs.highlightAll();</script>
</head>
<body>
  <nav class="w-sidebar" id="sidebar">
    <!-- preenchido por nav.js -->
  </nav>
  <main class="w-content">
    <header class="w-page-header">
      <a href="../index.html" class="w-back-link">← Visualizador</a>
      <span class="w-breadcrumb">Wiki / [Título]</span>
    </header>
    <!-- CONTEÚDO DA PÁGINA AQUI -->
    <footer class="w-nav-footer">
      <a class="w-prev" href="[anterior].html">← [Nome anterior]</a>
      <a class="w-next" href="[próximo].html">[Nome próximo] →</a>
    </footer>
  </main>
  <script src="nav.js"></script>
</body>
</html>
```

---

## TAREFA 2 — CSS E NAVEGAÇÃO

### 2.1 — `wiki/style.css`

Requisitos:
- Prefixo `.w-` em todas as classes para não conflitar com o site
- Sidebar fixa à esquerda, largura 240px
- Conteúdo com `max-width: 820px`, centralizado
- Fonte sans-serif, tamanho base 16px, line-height 1.7
- Blocos de código com fundo cinza claro, padding, borda
- Tabelas com zebra striping, bordas suaves
- Mobile: sidebar colapsável com botão hamburguer
- Cores: verificar paleta do visualizador existente e reusar

Estrutura CSS mínima:
```css
.w-sidebar { /* sidebar fixa */ }
.w-content { /* área de conteúdo */ }
.w-page-header { /* header com link ao visualizador */ }
.w-nav-footer { /* botões anterior/próximo */ }
.w-status-proven { color: #2a7a2a; } /* ✓ Provado */
.w-status-empirical { color: #7a5a00; } /* ~ Empírico */
.w-status-open { color: #7a0000; } /* ? Em aberto */
.w-status-negative { color: #555; } /* ✗ Negativo */
```

### 2.2 — `wiki/nav.js`

Gerar sidebar dinamicamente. O arquivo deve conter:

```javascript
const PAGES = [
  { id: 'index',          title: '🏠 Início',                  file: 'index.html' },
  { id: '01-introducao',  title: '1. Introdução',              file: '01-introducao.html' },
  { id: '02-teoria-gf2',  title: '2. GF(2) e Espaço de Ciclos',file: '02-teoria-gf2.html' },
  { id: '03-teorema-q3',  title: '3. Teorema Q(n) = 3',        file: '03-teorema-q3.html' },
  { id: '04-algoritmo',   title: '4. O Algoritmo',             file: '04-algoritmo.html' },
  { id: '05-fase-local',  title: '5. Fase Local f∞(L)',         file: '05-fase-local.html' },
  { id: '06-benchmark',   title: '6. Benchmark',               file: '06-benchmark.html' },
  { id: '07-estimativas', title: '7. Estimativas N(n)',         file: '07-estimativas.html' },
  { id: '08-topologia',   title: '8. Análise Topológica',      file: '08-topologia.html' },
  { id: '09-extensoes',   title: '9. Extensões',               file: '09-extensoes.html' },
  { id: '10-resultados',  title: '10. Todos os Resultados',    file: '10-resultados.html' },
  { id: '11-codigo',      title: '11. Como Usar o Código',     file: '11-codigo.html' },
];
// Detectar página atual, gerar sidebar, marcar link ativo
```

**Testar:** abrir `wiki/index.html` no browser e confirmar que
a sidebar aparece e os links funcionam antes de continuar.

---

## TAREFA 3 — PÁGINAS DE TEORIA (prioridade alta)

### `02-teoria-gf2.html` — GF(2) e Espaço de Ciclos

Estrutura da página:

**Seção 1: O que é GF(2)?**
- Aritmética com 0 e 1 onde 1+1=0
- É a operação XOR — familiar para programadores
- Por que aparece naturalmente em grafos

**Seção 2: Tours como vetores**
- Cada aresta do tabuleiro = 1 bit
- Um tour ativo/inativo em cada aresta → vetor binário
- $\mathbf{1}_t \in \mathbb{F}_2^{|E|}$

**Seção 3: A restrição de grau-2**
- Todo vértice num tour usa exatamente 2 arestas
- Em GF(2): $\partial_1 \cdot \mathbf{1}_t = \mathbf{0}$
- Logo: $\text{Ham}(n) \subseteq \ker(\partial_1) = Z_1(G_n; \mathbb{F}_2)$

**Seção 4: O número de Betti**
$$\beta_1(n) = |E_n| - |V_n| + 1$$

Tabela com valores concretos:

| $n$ | $\|V\|$ | $\|E\|$ | $\beta_1$ |
|-----|---------|---------|-----------|
| 6   | 36      | 80      | **45**    |
| 8   | 64      | 168     | **105**   |
| 10  | 100     | 288     | **189**   |
| 12  | 144     | 440     | **297**   |
| 14  | 196     | 624     | **429**   |

**Seção 5: O deficit**
- $\text{rank}(\text{Ham}(n)) < \beta_1(n)$
- Quanto menor, mais "especiais" são os tours
- Spoiler: sempre é $\beta_1(n) - 3$ → ver próxima página

Incluir diagrama ASCII mostrando a relação:
```
F₂^|E|  ⊃  Z₁ = ker(∂₁)  ⊃  Span(Ham)
          ↑                ↑
          β₁ dimensões     β₁ - 3 dimensões
```

### `03-teorema-q3.html` — O Teorema Q(n) = 3

Estrutura da página em ordem pedagógica:

**Seção 1: A observação empírica**

| $n$ | $\beta_1$ | $\text{rank(Ham)}$ | deficit |
|-----|-----------|---------------------|---------|
| 6   | 45        | 42                  | **3**   |
| 8   | 105       | 102                 | **3**   |
| 10  | 189       | 186                 | **3**   |
| 12  | 297       | 294                 | **3**   |

"Para todo tabuleiro testado, o deficit é sempre 3.
Por que sempre 3? O que são essas 3 dimensões?"

**Seção 2: Os cantos de grau 2**
- Diagrama mostrando as 8 arestas obrigatórias
- Canto $(0,0)$: seus únicos movimentos são $(1,2)$ e $(2,1)$
- Em qualquer tour: ambas DEVEM ser usadas
- Relação: $e_1(c) + e_2(c) \in \text{row}(\partial_1)$
- No quociente: $[e_1(c)] = [e_2(c)] =: r_c$

**Seção 3: A estrutura do quociente**
- 4 cantos → 4 representantes $r_{c_1}, r_{c_2}, r_{c_3}, r_{c_4}$
- São linearmente independentes no quociente (Lema 2.9)
- Geram $W \cong \mathbb{F}_2^4$

**Seção 4: O teorema (linguagem acessível)**
- Os pares XOR entre arestas obrigatórias geram um espaço de dimensão 3
- Não 4, porque os 4 cantos têm uma relação global
- A soma de todos os 4 representantes = 0 (paridade global)
- $\dim(\ker(\sigma)) = 4 - 1 = 3$

**Seção 5: Os 3 ciclos proibidos (6×6)**

Incluir plots dos 3 ciclos se disponíveis, ou diagrama ASCII:
```
f₀ (hexágono NW, 6 arestas):
  A6 → C5 → A4 → B6 → D5 → B4 → A6

f₁ (hexágono NE, reflexo de f₀):
  A6 → C5 → E4 → F6 → D5 → B4 → A6

f₂ (octógono diagonal, 8 arestas):
  A6 → C5 → E4 → D2 → F1 → E3 → D5 → B4 → A6
```

"Esses ciclos existem em $H_1$ mas NENHUM tour hamiltoniano
os contém. Eles são os 3 'espaços proibidos'."

**Seção 6: O detector mais simples**
$$\phi_1: x_{F6\text{-}D5} \oplus x_{B3\text{-}A1} \equiv 0 \pmod{2}$$

"Duas arestas obrigatórias. Para qualquer tour: $1 \oplus 1 = 0$ ✓.
Para $f_1$: $1 \oplus 0 = 1 \neq 0$ → detectado."

**Seção 7: O que falta (lacuna aberta)**

<div class="w-status-open">
? **Tightness:** provar que $\text{rank}(\text{Ham}(n)) = \beta_1(n) - 3$
para todo $n \geq 6$ par. A desigualdade $\leq$ está provada;
a igualdade $=$ permanece em aberto.
</div>

---

## TAREFA 4 — PÁGINAS DE ALGORITMO

### `04-algoritmo.html` — O Algoritmo

**Seção 1: O problema computacional**
"Queremos gerar K tours fechados num tabuleiro n×n.
Um tour fechado é um circuito hamiltoniano no grafo do cavalo."

**Seção 2: Três camadas de eficiência**

Apresentar como progressão:

```
Backtracking puro
  ↓ + Propagação R2
Backtracking + R2           (elimina arestas impossíveis)
  ↓ + Union-Find incremental
Backtracking + R2 + UF      (detecta sub-ciclos cedo)
  ↓ + Heurística f∞(L)
Backtracking + R2 + UF + f∞ (guia a busca para regiões férteis)
```

**Camada 1 — Propagação R2:**
- Se um vértice tem grau 2, suas duas arestas são obrigatórias
- Isso propaga: fixar uma aresta pode forçar outras
- Implementação: fila de propagação em O(1) por aresta

**Camada 2 — Union-Find incremental:**
- A cada aresta fixada em 1, atualizar componentes conexas
- Regra: componente com $k < n^2$ vértices e todos com grau 2 → sub-ciclo → podar
- Reduz razão 2-fatores/tour de **18.1×** para **1.00×** no 10×10

**Camada 3 — Heurística f∞(L):**
- Usar a tabela de frequências por nível (ver Seção 5)
- Ordenar arestas: priorizar as que mais se afastam de 0.5
- Guia o backtracking para regiões onde tours são mais prováveis

**Seção 3: Pseudocódigo**

```
função knight_tours(n, K):
  inicializar grafo, Union-Find, fila R2
  propagar R2 das arestas obrigatórias (cantos)
  
  tours_encontrados = []
  
  enquanto len(tours_encontrados) < K:
    aresta = escolher_melhor_aresta_por_f∞()
    
    para valor em [melhor_valor, outro_valor]:
      fixar aresta = valor
      propagar R2()
      
      se sub-ciclo detectado por UF:
        desfazer, continuar
      
      se tour completo:
        tours_encontrados.append(tour)
        desfazer
        continuar
      
      recursão()
      desfazer
  
  retornar tours_encontrados
```

**Seção 4: Como usar**

```python
# Instalar: pip install numpy (única dependência)
from knight_tours import knight_tours, verify_tour

# Gerar 10 tours fechados num tabuleiro 10×10
tours = knight_tours(n=10, K=10, seed=42)

# Verificar um tour
assert verify_tour(tours[0], n=10)

# Caminhos abertos: de A8 (0,0) até H1 (7,7)
from knight_tours import knight_path
path = knight_path(n=8, start=(0,0), end=(7,7), K=1)
```

### `05-fase-local.html` — Fase Local f∞(L)

**Seção 1: A observação**
"Ao longo de muitos tours em tabuleiros grandes, a frequência
de cada aresta converge para um valor que depende apenas
de sua distância à borda — não de sua posição exata."

**Seção 2: O experimento**
Incluir plot de convergência se disponível.

**Seção 3: A tabela universal**

| Nível L | $f_\infty(L)$ | Interpretação |
|---------|---------------|---------------|
| 0 | 0.528 | Borda: aresta ativa em ~53% dos tours |
| 1 | 0.193 | Anel interior 1: muito menos ativa |
| 2 | 0.198 | Anel interior 2 |
| 3 | 0.294 | Anel interior 3 |
| 4 | 0.247 | Anel interior 4 |
| 5 | 0.261 | Observado só em n=14 |

"6 números descrevem o comportamento assintótico de um grafo
infinito. Isso é o que torna a heurística universal."

**Seção 4: Por que funciona**
- Distância à borda captura a "dificuldade local"
- Arestas de borda (L=0) são frequentes porque os cantos as forçam
- Interior (L≥1) tem menos restrições → frequências menores
- O valor 0.5 seria "máxima incerteza" — desviar de 0.5 é informação

**Seção 5: Conexão com a matriz de transferência**
"Computacionalmente, $f_\infty$ é o autovetor principal da
matriz de transferência de 2-fatores no grafo infinito $G_\infty$
— hipótese em aberto mas sustentada pelos dados."

---

## TAREFA 5 — PÁGINAS DE RESULTADOS

### `06-benchmark.html` — Benchmark

**Seção 1: Progressão de melhorias (10×10, K=200)**

| Método | Tempo (s) | Speedup | O que adicionou |
|--------|-----------|---------|-----------------|
| Z3 puro | 33.11 | 1× | baseline |
| Z3 + 8 obrigatórias | 54.42 | 0.61× | *mais lento* |
| BT v1 (pressão + R2) | 7.18 | 4.6× | R2 + Warnsdorff |
| BT v2 (+ Union-Find) | 1.31 | 25× | detecção de sub-ciclos |
| BT theory (+ f∞) | 0.67 | **48×** | heurística universal |

Incluir gráfico se disponível.

**Seção 2: Por que Z3 + obrigatórias ficou mais lento?**

"Resultado contra-intuitivo: adicionar informação ao Z3 piorou.

Razão: o Z3 já infere as 8 arestas obrigatórias por propagação
unitária a partir das restrições de grau-2. Adicionar essas
cláusulas explicitamente cria trabalho redundante para o solver."

**Seção 3: Time-to-first-tour**

| Método | Tempo para 1º tour |
|--------|-------------------|
| Z3 puro | 0.85s |
| BT theory | **0.05s** |
| Diferença | **17×** |

### `07-estimativas.html` — Estimativas N(n)

**Seção 1: O estimador de Knuth**

"Imagine percorrer a árvore de busca aleatoriamente,
tomando um filho ao acaso em cada nó. No final, anote
o produto do número de filhos em cada nó do caminho.
A média desse produto sobre muitas amostras é exatamente N(n)."

Fórmula:
$$W = \prod_{v \in \text{caminho}} b(v), \qquad \mathbb{E}[W] = N(n)$$

"Não-viesado: cada tour é alcançado com probabilidade $1/W$."

**Seção 2: Validação**

| $n$ | Método | Resultado |
|-----|--------|-----------|
| 6 | Exaustivo | **9.862** (exato) |
| 6 | Estimador (10k amostras) | IC 95% contém 9.862 ✓ |
| 8 | Literatura | $1.33 \times 10^{13}$ |
| 8 | Estimador | IC 95% contém literatura ✓ |

**Seção 3: Resultados**

| $n$ | $\hat{N}$ | IC 95% |
|-----|-----------|--------|
| 6 | 9.862 | — |
| 8 | $7.4 \times 10^{12}$ | $[10^{12},\ 1.8 \times 10^{13}]$ |
| 10 | $2.4 \times 10^{22}$ | $[5.6 \times 10^{20},\ 6.8 \times 10^{22}]$ |
| 12 | $1.3 \times 10^{33}$ | $[7.8 \times 10^{31},\ 3.5 \times 10^{33}]$ |

**Seção 4: Em perspectiva**

$N(10 \times 10) \approx 2.4 \times 10^{22}$ — 24 sextilhões:
- **Mais** que grãos de areia na Terra ($\approx 7 \times 10^{18}$)
- **Comparável** ao número de estrelas no universo observável ($\sim 10^{23}$)
- **Menos** que moléculas de água num copo ($\sim 10^{25}$)

**Seção 5: Lei de crescimento**
$$N(n) \approx 1.82^{n^2}$$

"O número de tours cresce mais rápido que qualquer polinômio em n."

### `08-topologia.html` — Análise Topológica

**Seção 1: A razão tours/2-fatores**
$$r(n) := \frac{N_{\text{tours}}(n)}{N_{2\text{-fatores}}(n)}$$

"Quantos dos grafos com grau-2 em todo vértice são de fato tours
(conectados)? Essa razão mede quão difícil é a conectividade."

**Seção 2: Decaimento exponencial**

| $n$ | $r(n)$ | $\log_2(1/r)$ |
|-----|--------|----------------|
| 6 | 0.272 | 1.88 |
| 8 | 0.250 | 2.00 |
| 10 | 0.166 | 2.59 |
| 12 | 0.101 | 3.31 |

$$r(n) \approx 0.84 \cdot e^{-0.169n}$$

Incluir plot de decaimento se disponível.

**Seção 3: A transição de fase**

"O deficit linear é $Q(n) = 3$. O custo informacional da
conectividade é $\log_2(1/r(n))$. Quando esses dois números
se cruzam, ocorre uma transição:"

```
n ≤ 10:  log₂(1/r) < 3  →  cantos dominam a dificuldade
n ≥ 12:  log₂(1/r) > 3  →  conectividade adiciona restrições independentes
```

"Para tabuleiros pequenos, os 4 cantos explicam a maior parte
da diferença entre 2-fatores e tours. Para tabuleiros grandes,
a conectividade do interior domina — e o Q(n)=3 explica cada
vez menos."

---

## TAREFA 6 — EXTENSÕES, RESULTADOS E CÓDIGO

### `09-extensoes.html` — Extensões

**Seção 1: Simetria D₄**
- O tabuleiro tem 8 simetrias (4 rotações + 4 reflexões)
- Tours relacionados por D₄ são geometricamente idênticos
- No 6×6: 9.862 tours → 1.245 órbitas canônicas (~8× redução)
- 13 tours têm simetria própria ($|\text{Stab}| \in \{2, 4\}$)
- Aplicação prática: gerar só canônicas e expandir por D₄

**Seção 2: Caminhos Abertos**
- Tour fechado: grau 2 em todo vértice
- Caminho aberto: grau 1 em A e B, grau 2 nos demais
- Diferença na restrição: $\partial_1 \cdot \mathbf{1}_t = \mathbf{1}_A + \mathbf{1}_B$
- No 10×10: um caminho encontrado em ≈ 6ms
- Interface: `knight_path(n, start, end, K)`

**Seção 3: Divide-and-Conquer (em desenvolvimento)**
- Técnica de McKay (1997)
- Cortar o tabuleiro em duas metades
- Enumerar fragmentos hamiltonianos em cada metade
- Combinar pelos perfis de fronteira
- Para o 8×8: viabilizaria contagem exata (estimativa atual: $\sim 10^{13}$)

### `10-resultados.html` — Todos os Resultados

**Tabela principal:**

| Resultado | Tipo | Status | Valor |
|-----------|------|--------|-------|
| $Q(n) = 3$ | Teórico | <span class="w-status-proven">✓ Provado</span> | para todo $n \geq 4$ |
| deficit$(n) = 3$ | Empírico | <span class="w-status-empirical">~ Verificado</span> | $n \in \{6,8,10,12\}$ |
| 3 ciclos proibidos | Teórico | <span class="w-status-proven">✓ Construído</span> | $f_0, f_1, f_2$ no 6×6 |
| nós/tour = O(1) | Empírico | <span class="w-status-empirical">~ Verificado</span> | $[3.9, 5.3]$ para $n=6..14$ |
| $f_\infty(L)$ universal | Empírico | <span class="w-status-empirical">~ Verificado</span> | 6 valores, $n \leq 14$ |
| Speedup vs Z3 | Benchmark | <span class="w-status-proven">✓ Medido</span> | 48× no 10×10 |
| $N(10 \times 10)$ | Estimativa | <span class="w-status-empirical">~ IC 95%</span> | $[5.6 \times 10^{20},\ 6.8 \times 10^{22}]$ |
| Transição de fase | Topológico | <span class="w-status-empirical">~ Descoberta</span> | entre $n=10$ e $n=12$ |
| Tightness | Teórico | <span class="w-status-open">? Em aberto</span> | $\text{rank(Ham)} = \beta_1 - 3$? |
| $\mu_{\min} > 0$ | Teórico | <span class="w-status-open">? Em aberto</span> | nós/tour = O(1) provado? |

**Resultados negativos (informativos):**

| Resultado | O que foi tentado | Por que não funcionou |
|-----------|------------------|-----------------------|
| Patch LUT | Deformar sub-ciclos via XOR | R2 funde componentes cedo demais; taxa de resgate 0% |
| XOR no Z3 | Adicionar cláusulas Q(n)=3 | Z3 já as infere; sem ganho |
| Decomposição N = N₂ · r(n) | Fatorar contagem | r(n) decai exponencialmente; não simplifica |

### `11-codigo.html` — Como Usar o Código

**Seção 1: Instalação**
```bash
git clone https://github.com/[usuario]/knight_tour
cd knight_tour
pip install numpy  # única dependência
```

**Seção 2: Uso básico**
```python
from knight_tours import knight_tours, verify_tour

# Tours fechados
tours = knight_tours(n=10, K=100, seed=42)
print(f"Encontrados: {len(tours)} tours")

# Verificar um tour
ok = verify_tour(tours[0], n=10)
print(f"Tour válido: {ok}")

# Com simetria D₄ (só canônicas)
from knight_tours import knight_tours_canonical
canonical = knight_tours_canonical(n=6, K=500)

# Caminhos abertos
from knight_tours import knight_path
path = knight_path(n=10, start=(0,0), end=(4,3), K=1)

# Estimar N(n)
from tour_count_estimator import estimate_N
mean, ci_low, ci_high = estimate_N(n=10, M=1000)
print(f"N(10×10) ≈ {mean:.2e}  IC95%: [{ci_low:.1e}, {ci_high:.1e}]")
```

**Seção 3: Estrutura do projeto**
```
knight_tour/
├── knight_tours.py          ← motor principal (standalone)
├── tour_count_estimator.py  ← estimador de Knuth
├── knight_tours_optimized/  ← pacote com extensões
│   ├── __init__.py
│   ├── core.py
│   ├── d4_symmetry.py
│   └── open_paths.py
├── wiki/                    ← esta wiki
└── [subdiretórios de experimentos]
```

**Seção 4: Como rodar os experimentos**
```bash
# Benchmark completo no 10×10
python incremental_subtour/scaling_benchmark.py --n 10 --K 200

# Verificar Q(n)=3
python deficit_theorem/verify_small_cases.py

# Fase local: convergência de f∞(L)
python local_phase_heuristic/compute_f_inf.py --n 12 14
```

---

## TAREFA 7 — PÁGINA INICIAL

### `wiki/index.html`

**Header:**
```
Knight Tour — Wiki de Pesquisa
[link para o visualizador]
```

**Introdução (2 parágrafos):**

"O passeio do cavalo é um problema clássico de combinatória:
pode um cavalo de xadrez visitar todas as casas de um tabuleiro
n×n exatamente uma vez e retornar ao início? Esta wiki documenta
uma investigação matemática e computacional da estrutura interna
desse problema — não apenas como encontrar um tour, mas
*por que* os tours existem e como eles se organizam."

"Os resultados incluem uma prova de que certos ciclos são
topologicamente proibidos, um algoritmo 48× mais rápido
que abordagens SAT, e estimativas do número de tours
(24 sextilhões no tabuleiro 10×10)."

**Mapa de resultados (tabela resumida):**

| Pergunta | Resposta | Página |
|----------|----------|--------|
| Quantas dimensões os tours "perdem"? | Sempre 3 ($Q(n)=3$) | [Teorema Q(n)=3] |
| Qual o algoritmo mais eficiente? | BT theory: 48× vs Z3 | [Benchmark] |
| Quantos tours existem no 10×10? | ≈ $2.4 \times 10^{22}$ | [Estimativas] |

**3 destaques com links diretos:**
1. 📐 **Teorema Q(n) = 3** — a obstrução algébrica dos cantos
2. ⚡ **O algoritmo** — de 33s a 0.67s no 10×10
3. 🔢 **24 sextilhões** — estimativa de N(10×10)

**Linha do tempo cronológica:**
```
[Módulo → resultado → quando foi feito]
(preencher com base no wiki_audit.md)
```

---

## TAREFA 8 — INTEGRAÇÃO COM SITE EXISTENTE

Após auditoria, identificar o ponto de entrada do site e
adicionar um único link para a wiki.

**Exemplo para HTML puro:**
```html
<!-- Adicionar ao menu do site existente -->
<a href="wiki/index.html" class="nav-link">📖 Wiki</a>
```

**Nunca modificar o visualizador.**
**A wiki é adição, não substituição.**

---

## TAREFA 9 — REVISÃO FINAL

Verificar antes de reportar conclusão:

- [ ] Todos os links internos da sidebar funcionam
- [ ] MathJax renderiza fórmulas em todas as páginas com LaTeX
- [ ] highlight.js formata todos os blocos de código
- [ ] Imagens carregam com caminhos relativos corretos
- [ ] Mobile: sidebar colapsável funciona
- [ ] Link "← Visualizador" funciona em todas as páginas
- [ ] Botões anterior/próximo corretos em todas as páginas
- [ ] `wiki_audit.md` foi atualizado com o que foi efetivamente criado

---

## PRINCÍPIOS DE QUALIDADE

**1. Didático antes de completo**
Cada página deve ser entendível por um estudante de
computação que conhece teoria de grafos básica.
Não assumir que o leitor leu qualquer outro arquivo do projeto.

**2. Números concretos sempre**
Cada página deve ter pelo menos um resultado numérico
e uma tabela ou plot quando disponível.

**3. Honestidade sobre status**
Separar claramente com indicadores visuais:
- ✓ Provado formalmente
- ~ Verificado empiricamente  
- ? Conjectura em aberto
- ✗ Resultado negativo (informativo)

**4. Preservar o visualizador**
A wiki é complementar. Nunca modificar o visualizador.
O visualizador é o "wow factor" do site.

**5. Links verificados**
Não criar links para páginas que não existem ainda.
Criar primeiro o esqueleto de todos os arquivos, depois o conteúdo.
