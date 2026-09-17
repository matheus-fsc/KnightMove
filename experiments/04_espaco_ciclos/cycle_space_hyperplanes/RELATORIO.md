# Restrições GF(2) como hiperplanos em Z₁ — relatório

**Script:** `gf2_hyperplanes.py` · **Figura:** `results/gf2_hyperplanes.png`
**Tabuleiro:** cavalo 6×6. **Toda a aritmética das Partes 1–3 é sobre F₂** (mod 2),
nunca sobre ℝ.

## Achado central — a resposta à pergunta "qual é a direção/normal do plano?"

**Não existe direção/normal canônica.** A restrição não é um vetor que "aponta";
é uma **classe de equivalência de funcionais módulo `row(∂₁)`**. Concretamente:
o anulador dos tours tem dim 38, mas `row(∂₁)` (dim 35) zera em *todo* ciclo de
Z₁ — então `φ` e `φ ⊕ ρ` definem **exatamente a mesma restrição** para qualquer
`ρ ∈ row(∂₁)`. Verificação no script: a *mesma* obstrução aparece como `φ₀`,`φ₁`
de **suporte 2** e como `φ₂` de **suporte 20** — todos representantes legítimos
da mesma classe no quociente `anulador / row(∂₁)` (dim 3 = Q).

Sobre F₂ a forma bilinear `⟨φ,x⟩ = Σ φ_e x_e mod 2` é **degenerada**, e
`row(∂₁)` é precisamente seu radical sobre Z₁. Logo não há produto interno que
selecione uma "normal" única. A pergunta "para onde o plano aponta?" não tem
resposta geométrica — o dado invariante é a *classe* do funcional, não um vetor.
É isto que substitui a intuição da normal: **a restrição vive no quociente, não
no espaço.**

## Âncoras de fidelidade (todas confirmadas)

| quantidade | valor | papel |
|---|---|---|
| tours fechados | **9.862** | âncora-mestre; todos ∈ ker(∂₁), grau par em todo vértice |
| dim Z₁ = β₁ = \|E\|−\|V\|+1 | **45** = 80−36+1 | dimensão do espaço de ciclos |
| rank(Ham) | **42** | dim do span dos 9.862 tours |
| deficit = β₁ − rank(Ham) | **3** | = Q(6); o número de hiperplanos |

Se qualquer uma falhasse o script abortaria (`die`). Nenhuma falhou.

## O dicionário "plano que corta" → objeto combinatório real

A intuição vaga era: *"uma restrição é um plano cortando um volume"*. O objeto
real é mais simples e mais rígido.

| intuição contínua | objeto GF(2) correto | número |
|---|---|---|
| "espaço de soluções" | Z₁ = ker(∂₁) ⊂ F₂⁸⁰, **2⁴⁵ pontos** (não um contínuo) | dim 45 |
| "um plano" | hiperplano = **ker(φ)** de um funcional φ: F₂⁸⁰→F₂ | codim 1 |
| "o plano corta o volume pela metade" | cada φ independente divide a contagem por **exatamente 2**: 2⁴⁵→2⁴⁴ | ÷2 |
| "3 planos deixam o miolo" | 3 cortes ⇒ 2⁴⁵→2⁴⁴→2⁴³→**2⁴²** | Q=3 |
| "as soluções vivem no miolo" | os 9.862 tours ⊂ ker φ₀∩ker φ₁∩ker φ₂ (subespaço 2⁴²) | 9862 ≪ 2⁴² |

### Onde a intuição VALE
- **Restrição = hiperplano.** Uma paridade linear `φ(x)=0` é literalmente o
  núcleo de um funcional — um subespaço de codimensão 1. Isso é fiel.
- **Cortes sucessivos encolhem por fator 2.** Como cada φ é independente
  (mod row(∂₁)), cada interseção é metade da anterior. A imagem de
  "fatias cada vez menores" é correta — só que discreta: 2⁴⁵→2⁴²,
  não um volume diminuindo continuamente.
- **As soluções ficam na interseção.** Os tours estão dentro de
  ker φ₀∩ker φ₁∩ker φ₂; os 3 hiperplanos *contêm* o conjunto de soluções.

### Onde a intuição QUEBRA
- **Não há vetor normal.** Sobre F₂ não existe produto interno positivo-definido
  que dê uma "direção perpendicular" canônica; φ é um elemento do *dual*,
  não uma seta no espaço. (Há a forma bilinear `⟨φ,x⟩=Σφ_e x_e mod 2`, mas ela
  é degenerada — `row(∂₁)` é todo o radical sobre Z₁.)
- **Não há volume contínuo nem ângulo.** "Metade do espaço" significa
  `2⁴⁴` pontos, não metade de uma medida. Não existe ângulo entre hiperplanos:
  só a dimensão da interseção (45−k).
- **O funcional não é único.** φ vive no quociente `anulador / row(∂₁)`.
  Há muitas representações da *mesma* restrição (φ de suporte 2 e φ de suporte
  20 podem ser equivalentes mod row(∂₁)). A "posição do plano" não é um dado
  geométrico — é uma classe de equivalência.

## Parte 2 — a restrição φ de suporte 2 (o coração)

`φ(x) = x_[F6-D5] ⊕ x_[B3-A1]` (detector dual mínimo da Seção 2.3).

- **φ=0 em 100% dos 9.862 tours** ⇒ ker(φ) contém todas as soluções (restrição válida).
- **φ=0 em ~50,19% de 200 mil ciclos genéricos** de Z₁ (combinações lineares
  aleatórias da base, excluídos os que calham ser tour). O contraste
  100% vs ~50% é a prova de que φ **não é trivial**: um funcional trivial
  (∈ row(∂₁)) daria 100% em *qualquer* ciclo; um informativo separa tours
  de não-tours.
- **O anulador** {ψ : ψ·t=0 ∀ tour} tem dim 38 = 80 − rank(T). Módulo
  `row(∂₁)` (dim 35, os funcionais que zeram em *todo* ciclo), sobram
  **3 funcionais independentes** = os 3 hiperplanos Q. Encontrados:
  - φ₀ = x_[A6-C5] ⊕ x_[F6-D5]  (suporte 2)
  - φ₁ = x_[A6-C5] ⊕ x_[B3-A1]  (suporte 2)
  - φ₂ = suporte 20 (envolve o anel intermediário do tabuleiro)

  e `φ = φ₀ ⊕ φ₁` — confirmado não-trivial mod row(∂₁), logo é genuinamente
  um dos 3 hiperplanos Q.

## Parte 3 — o funcional σ (o análogo correto do "(1,1,1)")

A tentação é pensar Q=3 como "o vetor (1,1,1) apontando para algo". O objeto
correto é um **funcional de paridade global dos 4 cantos**:

`σ(x) = x_[A6-C5] ⊕ x_[F6-D5] ⊕ x_[B3-A1] ⊕ x_[E3-F1]` (um representante por canto).

- Cada gerador XOR de canto `r_ci ⊕ r_cj` tem **σ = 0** (soma de duas coordenadas
  marcadas = 0 mod 2) — verificado nos 6 pares.
- W = span(r_c) ≅ F₂⁴ (4 cantos independentes). `ker(σ|_W)` = span dos XOR-pairs
  tem dim **4 − 1 = 3**. Esse "−1" é exatamente o que produz **Q = 4 − 1 = 3**.
- Leitura: σ não "aponta" para lugar nenhum — ele *mede* a paridade somada
  dos 4 cantos. Q=3 é a dimensão do que sobra quando essa única paridade global
  é fixada (topologicamente, β₁ do K₄ dos cantos = 6−4+1 = 3).

## Parte 5 — pesca espectral (stub, NÃO executado)

`pesca_espectral(T_states, T_triples)` está definida mas não roda: precisa da
T_bulk real de `transfer_matrix/` com **33.346 estados**. Ela: (1) monta a
matriz esparsa; (2) confere fidelidade trace→**36.236** 2-fatores,
λ₁≈**70,48**, λ₂≈**−39,09**; (3) extrai o espectro via `scipy.sparse.linalg.eigs`;
(4) procura autovalores repetidos (degenerescência ⇒ simetria escondida);
(5) testa comutação `S·T = T·S` para S = paridade de (r+c) e reflexões D₄;
(6) para cada S comutante checa se o autovetor de Perron (λ₁) é S-invariante;
(7) reporta o fator de redução de estados. Documentado no docstring:
**destrava potencialmente o n=8** (>187M → cabível por bloco simétrico),
**mas NÃO o n=12**.
