# Prompt para o agente — Fase 3: construtibilidade do parity-switcher

## Contexto

Repositório: `/home/math/Dev/knight_tour`. Trabalho anterior relevante:

- `tightness_witnesses/` — Fase 2 já executada. `knight_gf2.py` (infra GF(2) +
  enumerador de tours), `witness_map.py` (análise), `RESULTS.md` (resultados),
  `data/witness_map_n{6,8}.json`.
- `paper/notes/tightness_plan.md` — o plano geral.
- `paper/referencia/CNP2026_hamilton_space_pseudorandom.pdf` — **o paper-fonte**
  (Christoph–Nenadov–Petrova). Def. 2.2 e recipe (S1)–(S5) estão na §2.
- `paper/referencia/Conrad1994_knights_hamiltonian_path.pdf` — Teorema 3.1,
  caracteriza caminho hamiltoniano `s→t` no tabuleiro cheio.

Notação (fixada, use esta):
- `G_n` = grafo do cavalo `n×n`, bipartido, `|V|=n²`, `δ=2` (4 cantos), `Δ=8`.
- `C^⊥` = espaço de cortes = `row(∂₁)`, dim `|V|−1`.
- `Mand` = as 8 arestas incidentes a canto ("obrigatórias" — todo tour as contém).
- `X = C^⊥ + Span(XOR_pairs)`, onde `XOR_pairs` = `1_e + 1_f` para `e,f ∈ Mand`.
  `dim X = |V| + 2`.
- `C_n^⊥` = `{R : ⟨τ,R⟩ = 0 para todo tour τ}`. **Tightness ⟺ `C_n^⊥ = X`.**
- `Z_bulk` = espaço de ciclos de `G_n` menos os 4 cantos (ciclos que evitam cantos).

## O que já está estabelecido (não re-derive)

**Lema A (provado).** `perp(Z_bulk) = C^⊥ + ⟨vetores suportados em Mand⟩`, e
`perp(Z_bulk) = X ⊕ ⟨1_e⟩` para qualquer `e ∈ Mand` — codimensão exatamente 1.
Verificado para `n ∈ {6,8,10}`: dims `(35,38,39)`, `(63,66,67)`, `(99,102,103)`.

**Lema B (provado).** `C_n^⊥ ∩ perp(Z_bulk) = X`.
*Prova:* se `R ∈ perp(Z_bulk)` então `R = x + c·1_e`, `x ∈ X`. Se `c=1`, então
`⟨τ,R⟩ = 0 + 1 = 1` para todo tour `τ` (toda obrigatória está em todo tour),
logo `R ∉ C_n^⊥`. ∎

**Corolário (o que reabre a Fase 3).** Todo `R ∈ C_n^⊥ ∖ X` tem interseção
**ímpar** com algum ciclo que evita os cantos. Como `G_n` é bipartido, todo
ciclo é par. Logo o passo **(S2.a) do CNP está resolvido de graça** para
exatamente os `R` que uma prova por contradição precisaria tratar.

**Consequência de projeto — leia com atenção:** (S2.b) e (S3) **não dependem
de `R`**. São perguntas puramente estruturais sobre `G_n`. Portanto podem ser
testadas **sem** exibir um `R ∉ X` (que não existe no alcance computável, já
que tightness vale para `n ≤ 12`). É isso que torna esta fase testável.

**Restrição de cor (já verificada).** Os dois caminhos hamiltonianos de `W` vão
de `v₁` a `v_{k+1}`, que estão a distância `k` no ciclo `C`. Bipartido ⟹ cores
opostas ⟺ `k` ímpar. Conrad et al. Thm 3.1 exige cores opostas (n par). Logo
**`k` deve ser ímpar**. O menor caso não-degenerado é `k=3` (hexágono).

---

## Tarefa 1 — Estender a verificação dimensional (rápida, ~30 min)

Confirme os Lemas A e B para `n ∈ {12, 14, 16}` e para tabuleiros
**retangulares** `n×m` com `min(n,m) ≥ 6` (pelo menos `6×8`, `6×10`, `8×10`).

Isto é álgebra linear pura sobre GF(2) — **não precisa enumerar tours**. Basta
construir `cuts`, `cuts_bulk`, `mand`, `xors` e comparar ranks.

Predição falseável: `dim perp(Z_bulk) − dim X = 1` **sempre**.
Se der ≠ 1 em algum caso, pare e reporte — invalida o Corolário.

---

## Tarefa 2 — O teste principal: o switcher é construtível? (o grosso do trabalho)

Instancie a Definição 2.2 do CNP em `G_n` diretamente. **Leia a Def. 2.2 no PDF
original**, não em resumo — a Fase 2 registrou esta como a peça frágil.

Def. 2.2: `W ⊆ G` é um `R`-parity-switcher se consiste de um ciclo par
`C = (v₁,…,v_{2k})` com número ímpar de arestas em `R`, mais caminhos
**vértice-disjuntos** `P_i` ligando `v_i` a `v_{2k−i+2}`, para `2 ≤ i ≤ k`.

Para `k=3` (hexágono `C = v₁…v₆`): `P₂` liga `v₂↔v₆`, `P₃` liga `v₃↔v₅`.

### 2a. (S2.b) — os caminhos existem?

Para cada hexágono `C` no **bulk** de `G_n` (evitando os 4 cantos), `n ∈ {8,10,12}`:

- Existem `P₂` (`v₂↔v₆`) e `P₃` (`v₃↔v₅`) vértice-disjuntos entre si e
  internamente disjuntos de `C`?
- Enumere os `W` resultantes. Reporte: quantos hexágonos admitem `W`, e a
  distribuição de `|V(W)|`.

Cuidado com os graus: `v₂,v₃,v₅,v₆` ficam com `deg_W = 3`; `v₁,v₄` com
`deg_W = 2`. Nenhum vértice de `W` pode exceder o grau disponível em `G_n`.

### 2b. (S3) — o completamento hamiltoniano existe? ⭐ **este é o teste decisivo**

Para cada `W` válido da etapa 2a: existe caminho hamiltoniano de `v₁` a `v₄`
em `G_n ∖ (V(W) ∖ {v₁, v₄})`?

- **Filtro barato primeiro:** critério de cor. `v₁` e `v₄` têm cores opostas
  (garantido por `k=3`), mas o grafo com buraco precisa ter as classes de cor
  balanceadas corretamente. Cheque contagem de cores antes de buscar.
- **Busca:** use o motor BT v2 do repo (R2 + Union-Find incremental + pressão
  de vértice) adaptado para caminho `s→t` com vértices removidos. Já existe
  interface `knight_path` — verifique se aceita conjunto de vértices excluídos;
  se não, estenda.
- Reporte taxa de sucesso por `n` e por tipo de `W`.

### Interpretação dos resultados (defina antes de rodar)

| Resultado de 2b | Significado |
|---|---|
| Existe `W` com (S3) satisfeito | **O esquema CNP é instanciável em `G_n`.** A prova de tightness passa a ter todos os passos disponíveis; próximo passo é uniformizar em `n`. Resultado forte. |
| (S3) falha para todo `W`, com causa identificada | Obstrução estrutural real. Caracterize-a (paridade de cor? o buraco desconecta? grau?). Resultado negativo publicável. |
| (S3) falha sem causa clara | Provável limitação de busca, não teorema. Reporte como inconclusivo — **não** conclua impossibilidade. |

---

## Armadilhas conhecidas — leia antes de rodar

1. **⚠ Viés de amostragem (já queimou uma rodada).** Enumerador com Warnsdorff
   puro fecha tour na primeira tentativa e nunca explora saídas de níveis
   rasos, deixando arestas com frequência **zero** — que viram witnesses
   espúrios de peso 1 e "refutam" tightness falsamente. Aconteceu em `n=8`:
   deu `rank=96, deficit=9, 786 witnesses desconhecidos`; corrigido com ordem
   uniforme nos 4 primeiros níveis + jitter abaixo, `rank` foi a 102.
   **Se qualquer parte desta tarefa amostrar tours: certifique frequência > 0
   em TODAS as arestas antes de acreditar em qualquer resultado.**

2. **Não confunda `R ∈ X` com `R ∉ X`.** O `RESULTS.md` §2.2 atual contém um
   erro: conclui "nenhum switcher no bulk" a partir de "toda classe tem
   representante em Mand", mas isso **assume tightness** (é a conclusão) e só
   vale para `R ∈ X`. Para o `R ∉ X` hipotético, o Lema B garante o oposto.
   Não reproduza esse erro.

3. **Canto interno a `P_i` é permitido** (grau 2 em `W`), mas a paridade força
   o canto a estar em `C`, não em `P_i`. Se você encontrar `W` com canto,
   verifique qual dos dois casos é.

---

## Entregáveis

1. `tightness_witnesses/FASE3_RESULTS.md` — resultados, com a tabela de
   interpretação preenchida e veredicto explícito.
2. Código em `tightness_witnesses/`, reproduzível por linha de comando, com os
   comandos exatos documentados no topo do RESULTS.
3. Dados brutos em `tightness_witnesses/data/`.
4. **Correção do `RESULTS.md` §2.2** incorporando os Lemas A e B e o
   Corolário, substituindo a conclusão circular.
5. Se algo contradisser os Lemas A/B: **pare e reporte**, não contorne.

## Padrão de honestidade

- Reporte falhas e resultados negativos com o mesmo destaque dos positivos.
- Distinga sempre **"verificado computacionalmente para n ∈ {...}"** de
  **"provado para todo n"**.
- Se uma busca não achou algo, diga se é *"provado que não existe"* ou
  *"não encontrado com este orçamento"* — nunca confunda os dois.
- Não declare tightness provada ou refutada. Esta fase testa **um passo** de
  **uma** estratégia de prova.
