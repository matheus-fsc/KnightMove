# 06 — O invariante de deficit Q

A pergunta que organizou o resto do projeto: **por que o rank do espaço gerado
pelos tours fica exatamente 3 abaixo de `beta_1`?**

| diretório / arquivo | resultado |
|---|---|
| `deficit_theorem/` | `Q(n) = 3` provado estruturalmente para `n >= 4` (4 representantes independentes + kernel da soma); verificado em `n=8` (rank 102) e `n=12` (rank 294) |
| `Q_locality_theorem.py` | **localidade**: `Q(G_T \ S) = 3 <=> W_corners ⊆ S`, e `Q = max(0, k_deg2 - 1)`. Verificado em `n ∈ {6,8,10}` |
| `bulk_connectivity_proof/` | prova de conexidade do bulk para `n >= 6` por indução `n -> n+2` com bases `{6,7}`; é o que fecha `Q(n)=3` no caso geral. Verificada até `n = 30` |
| `tightness_witnesses/` | redução: tightness `<=>` `Z_bulk ⊆ Span(Ham)`. Provada por construção em `n ∈ {8,10,12}` via pares de tours cujo XOR é um hexágono, com enumeração completa e auditoria |
| `tightness_intermediate.py` | família híbrida toro→plano: `rank = beta_1 - Q` em todos os 10 pontos. Atenção: precisa de `K >= 250k` amostras; com 20k aparecem falsos contraexemplos |
| `rank_ham_torus.py`, `characterize_torus_rank.py` | o toro tem deficit 0 universalmente (`n ∈ {4,6,8,10}`), rank cheio |
| `residual_search/`, `residual_search_10x10/` | propagação de restrições + backtracking sobre o espaço residual |
| `cornerless_plane_q.py` | `Q` no plano sem cantos |

**O que Q significa:** `Q = 8 - c(Punc)` (identidade fechada), com os 3 graus de
liberdade perdidos localizados nos 4 cantos do tabuleiro. A formalização em
Lean 4 está em `formalization/lean/`.

**Contraexemplo importante:** o tabuleiro 5×8 tem bulk conexo e `Q = 3`, mas
deficit 9. Bulk conexo e `Q=3` **não** implicam tightness.
