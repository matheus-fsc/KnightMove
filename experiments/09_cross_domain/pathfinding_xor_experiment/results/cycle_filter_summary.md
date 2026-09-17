# Cycle-length filter + multi-cycle XOR + connectivity predictor

Follow-up de `pathfinding_xor.py`. Mesmas 60 instancias (reconstruidas deterministicamente via `used_seed`); A* repetido NAO foi re-executado (diversidade recarregada do JSON anterior).

## EXPERIMENT A — filtro de comprimento

| Filtro | n_cycles | compat% | valid_paths | diversity | desconexao %falhas |
|---|---:|---:|---:|---:|---:|
| Filter_4 | 101±69 | 16.7±6.6% | 13.9±7.1 | 0.125±0.042 | 98.4% |
| Filter_6 | 149±102 | 16.8±7.1% | 20.3±10.4 | 0.138±0.046 | 97.5% |
| Filter_8 | 177±122 | 17.7±7.9% | 24.7±11.7 | 0.155±0.054 | 97.1% |
| Filter_10 | 197±135 | 18.0±8.1% | 27.7±12.9 | 0.167±0.056 | 96.7% |
| Filter_20 | 241±167 | 18.3±8.4% | 34.2±16.4 | 0.200±0.065 | 96.0% |
| Baseline | 326±238 | 17.9±7.8% | 45.5±24.2 | 0.307±0.077 | 90.8% |

**Melhor por compatibilidade:** Filter_20  
**Melhor por valid_paths×diversity:** Baseline

## EXPERIMENT B — diversidade multi-ciclo (Filter_8)

| Conjunto | diversidade |
|---|---:|
| single | 0.155±0.054 |
| pair | 0.239±0.065 |
| triple | 0.175±0.183 |
| combined | 0.198±0.082 |
| A* repetido (prior) | 0.477±0.142 |

Gap fechado (combined−single): +4.3 pp · gap restante p/ A*: +27.8 pp

## EXPERIMENT C — preditor de conectividade

- Candidatos: 19578 (taxa-base válida 13.9%)
- **AUC-ROC: 0.977**
- Feature mais preditiva: `shares_edges` (coef padronizado +5.273)
- Precisão@0.5: 0.615 · Recall@0.5: 0.987
- Coeficientes (padronizados): {'cycle_length': 0.133, 'is_local': 0.56, 'shares_edges': 5.273, 'n_shared_edges': -0.586, 'path_length_ratio': -0.518}
- Útil como pré-filtro (AUC>0.75): **SIM**

## Respostas às perguntas-chave (honestas)

1. **Filtro curto aumenta compatibilidade?** Filter_8 = 17.7% vs Baseline = 17.9% → **0.99×** — hipótese de 2× REFUTADA.
2. **Multi-ciclo fecha o gap de diversidade?** combined = 0.198; alvo 0.40+ → **REFUTADA**. A* = 0.477.
3. **Propriedades preveem validade?** AUC = 0.977 → CONFIRMADA (>0.75).
4. **Threshold ótimo (compat×n_cycles):** Baseline (≈58.4 caminhos válidos esperados).

## Achado central (contraintuitivo)

O filtro de comprimento **NÃO** aumenta a compatibilidade (plana em ~17–18% de Filter_4 a Baseline). A correlação prévia comprimento×compat (r=−0.665) era **entre instâncias** (confundida com o tamanho da grade), não **por ciclo dentro** de uma instância. O preditor real é `shares_edges`: o ciclo compartilha alguma aresta com o caminho-base? (coef padronizado **+5.27** vs `cycle_length` +0.13). Recall@0.5 = 0.99: quase todo ciclo válido toca o caminho-base. O pré-filtro útil é **'compartilha aresta com P'**, não 'é curto'.

## Conexão com o trabalho anterior

1. **A obstrução de conectividade muda com o filtro?** NÃO — piora. Coluna `desconexao %falhas`: 90.8% (Baseline) sobe para 98.4% (Filter_4). Ciclos curtos que não tocam P falham por pura desconexão (sem violação de grau), então filtrar por comprimento *concentra* as falhas em desconexão em vez de eliminá-las. Confirma o achado do Knight's Tour: conectividade é a obstrução GLOBAL, indiferente a filtros locais ([[project_residual_search_10x10_done]]).
2. **Relação com o teorema de suficiência local (D&C do Knight's Tour):** o que torna um ciclo 'seguro' não é a localidade geométrica (`is_local` tem coef só +0.56) e sim **compartilhar um trecho contíguo com o caminho-base** (`shares_edges` +5.27). É a mesma lógica do D&C: a compatibilidade exige sobreposição na fronteira (interface comum), não mera proximidade — blocos adjacentes colam quando casam na borda compartilhada, não quando são pequenos ([[project_dnc_blocks]]).
3. **Veredito final:** o XOR GF(2) é viável como ENUMERADOR de alternativas em grafos esparsos, mas o filtro de comprimento é a alavanca errada: não muda a compatibilidade nem fecha o gap de diversidade (combined 0.198 « A* 0.477). A conectividade global permanece a obstrução insuperável por filtragem local — só é **contornável** restringindo-se a ciclos que tocam P (pré-filtro `shares_edges`, AUC 0.977), o que reduz drasticamente o custo sem sacrificar caminhos válidos. Em diversidade, porém, o A* penalizado continua estruturalmente superior: XOR produz variações locais, A* produz rotas globalmente distintas.
