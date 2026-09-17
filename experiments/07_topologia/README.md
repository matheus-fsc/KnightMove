# 07 — Topologias alternativas

Se o deficit `Q = 3` vem dos cantos, então mudar a topologia do tabuleiro deve
mudar `Q`. Foi assim que o invariante deixou de ser uma curiosidade numérica e
virou algo com conteúdo geométrico.

## Hierarquia de Q

| superfície | `beta_1` | `Q` | por quê |
|---|---|---|---|
| plano `n×n` | `E - V + 1` | **3** | os 4 cantos têm grau 2 |
| cilindro | | **0** | X modular remove 2 cantos; paridade do winding 50/50 em 18k tours |
| toro | | **0** | grau uniforme 8, rank cheio (`n ∈ {4,6,8,10}`) |
| Klein | | 0 ou 1 | `H1 = Z`, não `Z²`: a paridade off-diagonal `(0,1)/(1,0)` é sempre proibida. `n+m` par implica `Q=1`; `n+m` ímpar implica `Q=0` |

| arquivo | conteúdo |
|---|---|
| `torus_*.py` | varreduras, cirurgias e experimentos no toro |
| `klein_*.py` | busca exaustiva, varredura de paridades (11 paridades, 76k tours), poda preditiva |
| `cylinder_winding*.py` | winding numbers no cilindro |
| `sheared_*.py`, `sheared_shared/` | toro cisalhado com parâmetro `s` |
| `run_research_experiments.py` | orquestrador das varreduras |

**Conclusão:** `Q` mede a obstrução introduzida pela borda. Numa superfície sem
borda, ela desaparece.
