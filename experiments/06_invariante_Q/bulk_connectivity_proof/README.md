# Conexidade do bulk — prova + verificação

Prova (em papel) e verificação computacional do **Lema de Conexidade do
Bulk**, a única obstrução em aberto do paper:

> Para todo `n ≥ 6`, `Bulk(n) := Gₙ[Vₙ \ Corners(n)]` — o grafo do cavalo
> induzido em todas as casas exceto os 4 cantos — é **conexo**.

Equivale a: `Gₙ \ Mand(n)` é conexo nos vértices não-canto (remover as 8
arestas obrigatórias apenas isola os 4 cantos, que têm grau 2).

## Resultado

**PROVADO para todo `n ≥ 6`** (antes: só verificado por kernel em `n ∈ {6,8}`).

Estratégia: **indução com passo `n → n+2`**.
- Bases `n ∈ {6,7}` (uma por paridade), por BFS finito.
- Passo: imerge-se o `n`-tabuleiro como bloco interno `B` do `(n+2)`-tabuleiro;
  por invariância de translação `Gₙ₊₂[B] ≅ Gₙ`, logo `B` menos os cantos
  imersos é uma cópia de `Bulk(n)` (conexa por HI); cada canto imerso anexa-se
  ao bulk; e cada casa da moldura (anel externo) anexa-se a `B` por um
  movimento de cavalo explícito.

A prova é **completa e sem lacunas**. Threshold `n₀ = 6`.

## Arquivos

- `verify_bulk_connectivity.py` — (1) BFS direto até `n = 30` (incl. ímpares);
  (2) verificação explícita das 3 sub-afirmações do passo indutivo
  (Claims I, II, IV) para `n` até 30.
- `../paper/notes/bulk_connectivity_proof_draft.tex` — redação LaTeX da prova,
  pronta para revisão/integração (não editar `knight_tour_complete.tex` ainda).

## Rodar

```bash
python3 verify_bulk_connectivity.py
```

Saída esperada: todas as linhas `connected = YES`, todas as sub-afirmações
`PASS`, `VEREDICTO GERAL: PASS`.

## Impacto

Fecha o argumento de papel do `Lema indep(a)` (independência dos 4
representantes) para **todo `n ≥ 6`** → fecha o `Teorema Q(n)=3` geral no
nível de papel. Em Lean, cascateia para 3 dos 4 `sorry` estruturais
(`Q(n)=3` geral, Teorema U geral, `Q(n,m)=3` geral). Formalização é etapa
separada e posterior.
