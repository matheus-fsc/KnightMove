# Plano de decisões — fechamento da *tightness*

**Data:** 2026-08-07
**Substitui:** `paper/notes/tightness_plan.md` (cujo bloco de status ficou obsoleto após a Fase 3)
**Estado do paper:** `paper/knight_tour_complete_v3.tex` — *tightness* provada para `n ∈ {8,10,12}`

---

## 0. Onde estamos

### 0.1 Estabelecido (não re-derivar)

| # | Fato | Status |
|---|---|---|
| E1 | `Q(n) = 3` e `Q(n,m) = 3` para `n,m ≥ 6` | provado + Lean sem `sorry` |
| E2 | Conexidade do *bulk* `G_n ∖ Corners` | provado + Lean sem `sorry` |
| E3 | `dim Z_bulk(n) = β₁(n) − 4` | provado (usa E2) |
| E4 | **Redução:** `Z_bulk ⊆ Span(Ham) ⟹ rank(Ham) = β₁ − 3` | provado |
| E5 | **Tightness para `n ∈ {8,10,12}`** | **provado por construção**, sem amostragem |
| E6 | Hexágonos geram `Z₁(G_m)` para `m ∈ {5,…,12}` | verificado — **corner-free, inclui `m` ímpar** |
| E7 | `dim Z₁(m) = 3(m−1)(m−3)`, `|E| = 4(m−1)(m−2)` | identidade fechada |
| E8 | Salto `dim Z₁(n+2) − dim Z₁(n) = 12(n−1)` | identidade fechada |
| E9 | Localidade `r = 0`: todo quadrilátero do 10×10 é soma de hexágonos da própria caixa delimitadora | verificado (292/292) |
| E10 | `codim(X em perp(Z_bulk)) = 1`; `C_n^⊥ ∩ perp(Z_bulk) = X` | provado (12/12 tabuleiros) |
| E11 | Switcher ancorado em canto é impossível (bloqueio de grau em (S3)) | provado |
| E12 | Gadget tem `|V(W)| = 8` **constante** em `n` | verificado `n ∈ {6,8,10,12}` |

### 0.2 O alvo, decomposto

A formulação antiga — *"hexágonos certificados geram `Z_bulk`"* — esconde **duas** uniformidades independentes:

- **(U2) Geração.** Hexágonos geram `Z₁`. *Quase pronto:* enunciado corner-free (E6), localidade `r=0` (E9), alvo numérico exato `12(n−1)` (E8), subgrafo induzido interno idêntico ao tabuleiro isolado (invariância por translação do movimento).
- **(U1) Certificação.** Toda família geradora admite os pares `H_A, H_B` com `H_A ⊕ H_B = C`. *Gargalo real.* Hoje é **busca** (`ham_path` com orçamento), não construção. Exigiria teorema de caminho hamiltoniano em `n×n` **menos um patch limitado** com extremos prescritos — precisamente o que Conrad et al. e Cull–De Curtins **não** cobrem.

> **As Rotas A/B/C atacam apenas (U2).** Provar só (U2) deixa o resultado condicional:
> *"tightness ⟺ os `12(n−1)` hexágonos da moldura são certificáveis"*.
> Melhor que hoje, mas ainda condicional.

### 0.3 Correções registradas

- **Rota C, como originalmente formulada, era vacuosa.** "Quadrilátero é soma de hexágonos certificados" é trivialmente verdadeiro quando os certificados já geram `Z_bulk`. A versão informativa é **local** (E9).
- **A intuição de Cayley aponta para hexágonos, não quadriláteros.** "Comutadores geram" é falso por contagem: o núcleo de `ℤ⁸ → ℤ²` tem posto 6, então os comutadores apresentam `ℤ⁸` livre-abeliano; faltam 6 relatores independentes, e o menor tem comprimento 6 — são exatamente os hexágonos.
- **Rota B (cobertura por blocos) está descartada.** `r = 0` já é o resultado que a cobertura tentaria obter, sem o problema dos ciclos do nervo.

---

## 1. Os gates

### Gate 0 — Lean da redução ⟨começa já, paralelo, ~1 semana⟩

Formalizar E3 + E4: `dim Z_bulk = β₁ − 4` e `Z_bulk ⊆ Span(Ham) ⟹ rank = β₁ − 3`.

**Por que primeiro:** não é só "vale a pena", é **de-risking**. A redução é a espinha da cadeia inteira; falha sutil ali invalida tudo a jusante. É curto, é álgebra linear sobre `𝔽₂` que o Mathlib cobre bem, e senta sobre E2, já formalizado.

- **Saída:** teorema sem `sorry`.
- **Valor independente:** se (U1)/(U2) empacarem, este teorema sozinho já é contribuição publicável.

### Gate 1 — Identificar templates de translação ⟨2 dias⟩

Não é preciso certificar **todo** hexágono — basta que **alguma** família geradora seja certificada. Há folga: em `n=12`, 8.848/9.000 certificam e os certificados sozinhos geram.

**Pergunta:** qual a menor família geradora fechada por translação? Com `r=0` e alvo `12(n−1)`, a hipótese é ~12 formas-template repetidas por translação.

- **Se confirmar:** (U1) colapsa de *"teorema sobre todos os hexágonos"* para *"verificação finita de ~12 casos + argumento de translação"*. Deixa de ser gargalo.
- **Impacto colateral:** a base finita para Lean cai de ~17k pares de tours para ~12 — o que torna a formalização de (U1) viável, o que hoje não é.

⚠️ **Rodar antes do Gate 3.** É barato e pode reordenar tudo.

### Gate 2 — Teste do bloco 6×6 perfurado ⟨2 dias⟩

**A máquina para (U1) já existe e não está sendo usada: o catálogo D&C.**

O gadget tem `|V(W)| = 8` (E12) — cabe folgado num bloco 6×6. Então (S3) deixa de ser "caminho hamiltoniano em `n×n` menos um patch" e vira:

1. O patch está contido em **um** bloco 6×6
2. Enumerar no catálogo os caminhos do bloco **perfurado** com o perfil de fronteira certo
3. Demais blocos intactos — o catálogo de 25,6k já os cobre
4. Compatibilidade local ⟹ hamiltonicidade global (teorema próprio, 1639 trials, validado em 12×12 e 18×18)

Isso converte **busca com orçamento** em **consulta a catálogo finito** — e catálogo finito é o que se formaliza e o que se prova uniforme por translação. Custo: gerar o catálogo perfurado uma vez (finito, independente de `n`).

**Teste mínimo:** um template de hexágono, num bloco 6×6 interior; o bloco perfurado tem caminhos com os perfis necessários?

- **Sim:** (U1) vira construção; a indução carrega (U1) e (U2) juntas.
- **Não:** obstrução real identificada com 2 dias de custo.

### Gate 3 — (U2) por Rota A ⟨3 semanas, PARADA DURA⟩

Indução na moldura, com alvo numérico exato: exibir `12(n−1)` hexágonos que tocam a moldura, independentes módulo o span interno.

O subgrafo induzido na cópia interna `n×n` tem exatamente as mesmas arestas do tabuleiro `n×n` isolado (invariância por translação), então a hipótese de indução entra direto.

**Testar antes:** como E6 é corner-free e vale para `m` ímpar, tente indução **`m → m+1`** em vez de `m → m+2`. Passo mais fino, moldura em L em vez de anel completo, e mais bases disponíveis (`m = 5…12`). Verifique se o salto `m→m+1` tem estrutura mais simples que `12(n−1)`.

> 🛑 **Condição de parada: 3 semanas.** Se não sair, congelar (U2) como conjectura precisa e ir para o paper. Sem prazo isto vira item adiado indefinidamente — padrão já identificado neste projeto.

### Ordem de execução

```
Gate 0 (Lean redução) ──────────────────────────────►  [paralelo, 1 sem]

Gate 1 (templates, 2d) ─┐
                        ├─► Gate 3 ((U2) Rota A, 3 sem, PARADA DURA)
Gate 2 (bloco 6×6, 2d) ─┘        │
                                 └─► (U1) conforme resultado de 1–2
```

Gates 1 e 2 vêm antes do 3 porque são baratos e podem reordenar tudo.

---

## 2. Escopo do Lean

| Peça | Formalizar? | Razão |
|---|---|---|
| Redução (E3+E4) | ✅ **agora** | curto, Mathlib cobre, resultado durável, de-risking |
| (U2) geração | ⏸ depois de sair no papel | base finita = rank de matriz `𝔽₂`, viável por `native_decide` |
| (U1) certificação | ❌ **não** (hoje) | ~17k pares de tours como dado é insustentável |
| (U1) **se Gate 1 confirmar** | 🔄 reavaliar | base cai para ~12 templates — muda a conta radicalmente |

Certificados das bases `n ∈ {8,10,12}` vão para **apêndice computacional fora do Lean**, com hashes dos artefatos em `data/`.

---

## 3. Ações no paper

### 3.1 Imediata — ressalva na `tab:tightness`

A tabela da família híbrida usa rank de *ensembles* amostrais, que é **exatamente** o que a `sec:vies` da mesma versão manda desconfiar — e a armadilha do Warnsdorff já produziu deficit falso neste projeto (`rank 96 → 102`, 786 witnesses espúrios de peso 1).

Um parecerista que leia as duas seções em ordem vai perguntar. **Ação mínima:** nota de rodapé ligando `tab:tightness` à `sec:vies`, declarando os 5 pontos como evidência amostral não-certificada.

### 3.2 Quando houver tempo — refazer os 5 pontos com certificados

Generalizar `hexagon_in_tourspace.py` para tabuleiros híbridos. Trabalho pequeno — a máquina existe.

---

## 4. Armadilhas — não repetir

1. **⚠ Viés do amostrador.** Warnsdorff puro fecha tour na primeira tentativa; níveis rasos nunca explorados; arestas com frequência **zero** viram witnesses espúrios de peso 1. Produziu `rank=96, deficit=9` falso em `n=8`. **Antes de interpretar qualquer rank amostral, certifique frequência > 0 em TODAS as arestas.** A prova de E5 foi construída para não amostrar tours em ponto algum — mantenha essa disciplina.

2. **⚠ Circularidade `R ∈ X` vs `R ∉ X`.** "Toda classe não-trivial tem representante em `Mand`" descreve as classes de `X`; afirmar que são todas as de `C_n^⊥` **é** a tightness. Erro cometido e retratado na Fase 2 §2.2.

3. **⚠ Escopo do bloqueio de grau.** E11 é **verdadeiro**, mas restrito a `R ∈ X` (é a paridade que força o canto em `C`). Não é "o argumento caiu" — é "o argumento vale, com escopo".

4. **⚠ Compilação do paper.** Compilar da **raiz do repo**, não de dentro de `paper/`; caso contrário as 8 figuras somem **silenciosamente**:
   ```bash
   cd /home/math/Dev/knight_tour
   pdflatex -interaction=nonstopmode -output-directory=paper paper/knight_tour_complete_v3.tex
   ```

5. **⚠ Cortes de busca entram no predicado medido.** `build_switchers(max_w=20)`
   truncava candidatos, e a ordem de enumeração depende dos índices de vértice,
   logo de `n`. Os números "53 em n=10, 73 a partir de n=12, limiar n₀=12" são
   **artefato**; sem corte dá 94 para todo `n`. **Antes de reportar
   (in)dependência em `n`, sature todo corte e exiba a saturação.** Ver
   `GATES_RESULTS.md` §G1.12.

---

## 4-bis. As armadilhas, viradas artefato de código

Três protecções que antes eram nota mental e hoje são código. Estão em três
arquivos diferentes; esta é a lista canônica.

| protecção | onde vive | o que evita |
|---|---|---|
| **direção do viés no docstring** | `rect_exceptions.py`, `flip_paths.py`, `flip_realizable.py` | usar um método fora do sentido em que ele é confiável |
| **saturação de corte reportada** | `flip_graph_n8.py` (`ball_truncated`), `flip_paths` (`exhausted`) | repetir o `max_w = 20` |
| **fecho por simetria antes da busca** | `board_symmetry.py`, usado por `rect_exceptions.py` | gastar busca exaustiva no espelho de algo já certificado |

### A regra da direção do viés, em forma utilizável

Todo método aqui é **unilateral**. Antes de citar um resultado, pergunte em
que sentido o método pode errar:

| método | pode afirmar | **não** pode afirmar |
|---|---|---|
| bola de BFS | "C é flip-realizável" (certificado exibido) | "C não é" — a bola é truncada e assimétrica |
| `flip_paths` exaustivo | ambos, **se** `exhausted=True` | qualquer coisa com `exhausted=False` |
| enumeração completa de tours | ambos | — (é a verdade) |

O caso que motivou: a primeira `viable` de `flip_paths` podava demais, logo
**perdia flips existentes**, logo **fabricava exceções**. A conclusão em risco
("6×8 tem 4 exceções") era da mesma forma que o erro produzia. A proteção não
foi corrigir o bug — foi a validação **conter negativos verdadeiros** (os 12
do 6×6 e os 4 do 6×7, estes provados por enumeração de 1.067.638 tours).
**Validar só contra positivos deixaria o lado que importa sem cobertura.**

### Simetria: use antes e depois

O grupo é Klein (ordem 4) se `R ≠ C`, D₄ (ordem 8) se `R = C`, e
flip-realizabilidade é constante em cada órbita.

- **antes da busca** — teste só um representante por órbita
  (`orbit_representatives`): 828 → 223 em 6×7, 532 → 70 em 6×6;
- **depois da bola** — feche os certificados sob o grupo
  (`close_certificates`) antes de mandar o resíduo para a busca.

Custo medido de não fazer o segundo: no 7×8 o resíduo de 1 hexágono consumiu
**4.781 s de busca exaustiva e voltou inconclusivo**; o transporte por simetria
o resolveu em milissegundos, com certificado verificado. Todo transporte é
**verificado, nunca assumido** — os dois tours transportados têm de ser
hamiltonianos e o XOR tem de dar exatamente o hexágono imagem.

⚠️ Ao reportar, **expanda as órbitas de volta**: "1 representante provado" em
6×7 significa **4 hexágonos**. O sumário de `rect_exceptions.py` faz isso.

---

## 5. Critério de sucesso

| Cenário | Resultado |
|---|---|
| Gate 0 apenas | Teorema da redução formalizado — publicável isolado |
| Gate 0 + (U2) | *"Tightness ⟺ os `12(n−1)` hexágonos da moldura certificam"* — condicional, mas enunciado preciso |
| Gate 0 + (U1) + (U2) | **Tightness para todo `n` par** — fecha a lacuna central do trabalho |
| Nenhum em 3 semanas | Paper com E1–E12 e a Conjectura `hex-geram`; já é melhoria grande sobre a v2 |

**Em qualquer cenário o paper é escrito.** Os gates decidem quão forte ele é, não se ele existe.

---

## 6. Arquivos

```
tightness_witnesses/
  knight_gf2.py              infra GF(2) + enumerador (sem viés)
  witness_map.py             Fase 2 — mapa dual
  hexagon_in_tourspace.py    Fase 3 — certificados
  hexagon_locality.py        localidade r=0 / geração
  RESULTS.md                 Fase 2 (§2.2, §2.3 retratados)
  FASE3_RESULTS.md           Fase 3 (§3.45 redução, §3.6 auditoria)
  FASE3_PROMPT.md            prompt do agente
  PLAN.md                    este arquivo
  data/                      dados brutos + logs de auditoria

paper/
  knight_tour_complete_v3.tex   tightness provada n ∈ {8,10,12}
  notes/citation_audit.md       auditoria bibliográfica
  referencia/                   25 PDFs + BUSCA_CAPES.md
```
