# Resultados negativos

Esta página é deliberada. A maior parte do que foi testado não funcionou, e o
padrão comum entre as falhas é mais informativo que os sucessos.

## O padrão

**GF(2) captura a condição de grau par com perfeição. Não captura
conectividade.**

Grau par é linear sobre `F_2^E`; "ser um único ciclo" não é. Todo método que
gera candidatos por XOR produz 2-fatores de grau correto que se quebram em
componentes. Em todo domínio onde a solução precisa ser conexa, a taxa de
acerto desaba.

## No próprio passeio do cavalo

| hipótese | veredito | detalhe |
|---|---|---|
| cláusulas XOR aceleram o Z3 | **refutada** | speedup 1,00x em 6×6 e 10×10. DPLL(XOR) já é padrão. Subproduto útil: confirmou `deficit(10×10) = 3` |
| a tabela de fase local `f∞` acelera a busca | **refutada** | speedup constante ~1,0 em `n ∈ {6..16}`, com outliers catastróficos em `n >= 12`. O ganho atribuído a ela vinha da pressão de vértice |
| speedup `∝ (n-8)²/n²` | **falsa** | speedup é flat |
| LUT de patch XOR resgata sub-tours | **refutada** | corretude OK (9862/9862), mas taxa de resgate **0%** em 16k casos. O union-find funde componentes cedo demais para o patch disparar |
| `NOT(A ∧ B)` acelera o 10×10 | **não** | as bifurcações não são estritas: `P(A∧B) ≈ 0,10` no 10×10 contra 0 no 6×6. Transição qualitativa entre os dois tamanhos |
| exclusões de ordem 3 ajudam no 10×10 | **não** | 118 triplas provadas, mas **todas estruturais** (o grau 2 já bloqueia). Zero genuínas. `NOT(A∧B∧C)` dá slowdown de 0,71x |
| propagação de restrições resolve o 6×6 | **parcial** | R1..R6 fixa 0 variáveis livres (os cantos já estão isolados). O backtracking acha os 9862 tours em 53s, mas revela 3560 2-fatores multi-ciclo válidos localmente |
| `H_MART` forte (custo i.i.d. por nó) | **refutada** | spread 28%, chi-quadrado rejeita homogeneidade. A versão fraca se sustenta |
| conjectura `rank = número de órbitas` (toro) | **refutada** | o toro tem rank cheio universalmente |
| XOR de ciclos enumera caminhos `s -> t` | **refutada** | cobre 0,0086% (4 de 46666 no 6×6); satura em ~10 no 8×8. Falha dominante: 73-84% dos ciclos não alternam |

## Fora do passeio do cavalo

| domínio | veredito |
|---|---|
| **TSP** | XOR de ciclos fundamentais da MST fica ~15% pior que 2-opt |
| **Pathfinding em grades** | compatibilidade ~18% (contra 3% no TSP), cai com `n`. Falha dominante: **desconexão, 90,8%**, não grau |
| **Sítios ativos de proteínas (GASS)** | **inviável**. As soluções são k-cliques transversais sob limiar métrico, não subgrafos pares. 0 de 17.367 XORs válidos. O "~5" da literatura é o número médio de resíduos por sítio (4,62 no M-CSA), não um limite algorítmico |

## Erros metodológicos cometidos e corrigidos

Registrados porque custaram tempo:

- **Amostra pequena gera falsos contraexemplos.** A verificação de tightness na
  família híbrida precisa de `K >= 250k`. Com 20k, contraexemplos apareceram e
  eram artefato.
- **Uma única semente falseia rank.** A medição `rank ≈ 21` no toro era efeito
  de semente única; a medição correta dá rank cheio.
- **Um corte de `max_w = 20` contaminou medições de certificação.** O "limiar
  `n_0 = 12`" reportado era artefato do corte; o valor real, sem limiar, é 94.
- **Um argumento circular foi retratado** na prova de tightness e substituído
  por prova por construção.

## Por que isso está documentado

Um repositório que só mostra o que funcionou não permite avaliar o método. O
diagnóstico "conectividade é a obstrução global, e GF(2) não a vê" só é
defensável porque foi testado em cinco domínios independentes e falhou do
mesmo jeito nos cinco.

**Código:** `experiments/09_cross_domain/`, `experiments/08_heuristicas/`
