# Plano de ataque à *tightness* via parity-switchers (CNP)

**Base:** M. Christoph, R. Nenadov, K. Petrova, *The Hamilton space of pseudorandom graphs*, JCTB **176** (2026), 254–267 (arXiv:2402.01447). Texto completo lido.
**Alvo:** $\rank_{\mathbb{F}_2}(\Ham(n)) = \beta_1(n) - 3$ para todo $n \geq 6$ par.
**Data:** 2026-08-06

> **STATUS (2026-08-06, após Fases 2 e 3):** ver
> `tightness_witnesses/RESULTS.md` e `tightness_witnesses/FASE3_RESULTS.md`.
> - **Fase 2 concluída** (n=6 exaustivo até peso 6; n=8 até peso 4).
>   Zero witnesses fora das 8 classes previstas em ambos.
> - **§1 (Lema 2.1 bipartido): confirmado e simplificado.** No caso bipartido
>   com `|V|` par, `1_E ∈ C^⊥` torna `R ↦ E∖R` uma involução que preserva a
>   classe; o `R` maximal é o complemento do detector mínimo, e (C1) sai de
>   graça. Não há argumento de extremalidade a fazer.
> - **Fase 3 executada e BEM-SUCEDIDA.** O esquema do CNP **é** instanciável
>   em `G_n`: switchers completos (C + P₂,P₃ + caminho de (S3)) verificados
>   para `n ∈ {6,8,10,12}`, com taxa 4,7% → 87,3% → 97,0% → **100%**.
> - **§3.1 do plano ("switcher no bulk") está certo, não impossível.** Uma
>   nota anterior afirmou o contrário; era erro circular, retratado no
>   §2.2 do RESULTS.
> - **§4.1 (controle de cor `k` ímpar): verificado** — e é satisfeito
>   automaticamente pela construção quando `k` é ímpar.
> - Extra: os **hexágonos do bulk geram `Z_bulk`** (n ∈ {6..14}), então `k=3`
>   é suficiente — não é escolha de conveniência.
> - **Fase 4 (uniformizar em `n`) é agora o passo seguinte**, com 3 lacunas
>   concretas listadas no §4 do FASE3_RESULTS.

---

## 0. A descoberta que muda o enquadramento

O Lema 2.1 do CNP **é a sua conjectura de tightness escrita em forma dual.** Traduzindo:

> **Lema 2.1 (CNP).** Seja $G$ hamiltoniano com $n$ ímpar e $\mathcal{C}_n(G) \neq \mathcal{C}(G)$. Então existe $R \subseteq G$ com:
> - **(C1)** $R \neq G$;
> - **(C2)** todo ciclo hamiltoniano contém um número **par** de arestas de $R$;
> - **(C3)** para toda partição $V = A \cup B$: $e_R(A,B) \geq e_G(A,B)/2$, e $R \neq G[A,B]$.

A prova é: $\mathcal{C}^\perp(G)$ é exatamente o **espaço de cortes**; tome $R$ como o elemento **de maior cardinalidade** em $\mathcal{C}_n^\perp(G) \setminus \mathcal{C}^\perp(G)$. A maximalidade força (C3), porque se algum corte tivesse $e_R(A,B) < e_G(A,B)/2$, então $R + G[A,B]$ seria um elemento maior do mesmo coset.

**Dicionário com o seu paper:**

| CNP | Você |
|---|---|
| $\mathcal{C}_n(G)$ | $\Span(\Ham(n))$ |
| $\mathcal{C}(G) = Z_1$ | $\ker \partial_1$ |
| $\mathcal{C}^\perp(G)$ = espaço de cortes | $R(n) = \mathrm{row}(\partial_1)$ ✅ **é a mesma coisa** |
| $R \in \mathcal{C}_n^\perp \setminus \mathcal{C}^\perp$ | um funcional anulador de tours fora de $\mathrm{row}(\partial_1)$ |
| $\dim(\mathcal{C}_n^\perp/\mathcal{C}^\perp)$ | **o seu `deficit`** |

E os seus **detectores duais $\varphi_0,\varphi_1,\varphi_2$ são exatamente witnesses $R$.** O seu $\varphi_1 = x_{F6\text{-}D5} \oplus x_{B3\text{-}A1}$ é um $R$ com $|R| = 2$ — duas arestas obrigatórias de cantos distintos. (C2) é trivial: todo tour contém as duas, logo $2 \equiv 0$.

**Reformulação da tightness:**

> $\rank(\Ham(n)) = \beta_1(n) - 3$ $\iff$ **todo** $R \subseteq G_n$ satisfazendo (C2) pertence a $\Span(\mathrm{XOR}_{\text{pairs}}) + \mathcal{C}^\perp(G_n)$.

Isto **já vale a pena escrever no paper mesmo que a prova falhe.** Transforma uma conjectura vaga ("existe um lema de deformação local?") numa pergunta precisa sobre certificados duais, conectada a uma linha de pesquisa ativa.

---

## 1. Lema 2.1 adapta para o seu caso? — SIM, e o argumento é seu

Obstáculo aparente: CNP exigem $n$ **ímpar**, e o autor destaca *"this is the only place where we use the fact that $n$ is odd"*. Você trabalha com $n$ par, logo $|V| = n^2$ par.

Onde a paridade é usada: para concluir **(C1)** $R \neq G$. O argumento deles é *"senão $R$ contém um ciclo hamiltoniano, que tem $n$ (ímpar) arestas em $R$, contradizendo (C2)"*.

Para $n^2$ par esse argumento morre — mas **(C1) continua valendo por outra razão**, e a razão é que **o grafo do cavalo é bipartido**:

> $R = E(G_n)$ satisfaz (C2) trivialmente (todo tour tem $n^2$ arestas, e $n^2$ é par). Mas $G_n$ é bipartido — o cavalo alterna cores a cada movimento. Logo $E(G_n) = G_n[A,B]$ com $A, B$ as duas classes de cor, ou seja $E(G_n) \in \mathcal{C}^\perp(G_n)$. Como witnesses vivem em $\mathcal{C}_n^\perp \setminus \mathcal{C}^\perp$, temos $R \neq E(G_n)$. ∎

**Resultado: Lema 2.1 vale para grafos bipartidos com $|V|$ par**, com (C3) saindo da maximalidade exatamente como no original. Isto é um lema pequeno, autocontido, provável em uma tarde — e **formalizável em Lean**, encaixando no seu pipeline existente.

⚠️ **Verificar antes de escrever:** que $\mathcal{C}_n^\perp \setminus \mathcal{C}^\perp \neq \emptyset$ requer $\mathcal{C}_n \neq \mathcal{C}$ — que você tem de graça (deficit $\geq 1$ pelos cantos, via Heinig). E confirmar que o argumento de maximalidade de (C3) não usa paridade em nenhum ponto. Pela leitura, não usa.

---

## 2. Onde o método deles NÃO transfere — e é melhor saber agora

Sendo direto: **a maquinaria quantitativa do CNP não serve para você.** O núcleo é o *Cycle Lemma* (2.3), cuja hipótese é:

> **(L1)** Para todo $S \subseteq V$ com $|S| \leq 2\ell$ e todos $x,y \notin S$, existe caminho $x \to y$ **dentro de $R$** de comprimento $\leq \ell - 1$.

Isso é uma condição de **expansão robusta de $R$**. Falha completamente no grafo do cavalo:

- $\Delta(G_n) = 8$, $\delta(G_n) = 2$ — grau limitado, sem expansão;
- diâmetro do grafo do cavalo é $\Theta(n)$, logo $\ell = \Theta(n)$;
- com $\ell = \Theta(n)$, "ciclo par de comprimento $\leq 2\ell$" deixa de ser uma restrição útil, e (L1) exigiria que $R$ tivesse diâmetro $\Theta(n)$ **após remover $\Theta(n)$ vértices** — o que num grafo de grau 8 é falso em geral.

Todos os Teoremas 1.1/1.3 e Corolários 1.4–1.6 do CNP são sobre $\delta(G) \geq n/2 + C$, jumbledness, $(n,d,\lambda)$-grafos. **Nada disso existe no seu regime.**

**A conclusão certa não é "abandonar", é inverter o uso:** você não vai usar os *teoremas* deles, vai usar o *esqueleto de prova* (S1)–(S5), substituindo os passos existenciais-probabilísticos por **construções explícitas** que você já tem. O grafo do cavalo é rígido demais para argumentos de expansão, mas rígido também significa **construtível à mão**.

---

## 3. A adaptação: switcher relativo, construído no *bulk*

Duas mudanças estruturais em relação ao CNP.

### 3.1 Versão *relativa* (você não quer deficit 0, quer deficit 3)

CNP provam $\mathcal{C}_n = \mathcal{C}$ (deficit 0). Você quer deficit **exatamente 3**. Então trabalhe no quociente:

$$X := \Span(\mathrm{XOR}_{\text{pairs}}) + \mathcal{C}^\perp(G_n), \qquad \text{Tightness} \iff \mathcal{C}_n^\perp = X.$$

Suponha $R \in \mathcal{C}_n^\perp \setminus X$, escolhido de cardinalidade máxima **no coset módulo $X$**. A maximalidade relativa deve dar uma versão de (C3) — e, crucialmente, deve dar informação sobre $R$ **restrito ao bulk**, porque módulo $X$ você já quocientou tudo que os cantos geram.

🎯 **Meta:** um ciclo par com número ímpar de arestas de $R$ que **evita os 4 cantos**. Se o switcher vive no bulk, o argumento não pode ser bloqueado pelas arestas obrigatórias — que são exatamente a fonte do deficit 3.

### 3.2 O switcher como *gadget finito fixo*

Definição 2.2 (parity-switcher): ciclo par $C = (v_1,\ldots,v_{2k})$ com número **ímpar** de arestas em $R$, mais caminhos vértice-disjuntos $P_i$ ligando $v_i$ a $v_{2k-i+2}$ para $2 \leq i \leq k$. A propriedade-chave: $W$ contém **dois caminhos hamiltonianos** de $v_1$ a $v_{k+1}$ com **paridades diferentes** de arestas de $R$ — cada $P_i$ aparece nos dois, e cada aresta de $C$ aparece em exatamente um.

Isso é o coração e é puramente combinatório — **não usa aleatoriedade nenhuma.** É um gadget.

**Sua vantagem:** você já tem a LUT de diamantes (`patch_LUT`, universal por nível $L$, satura em 6 diamantes/aresta para $L \geq 2$) e o catálogo de 25,6k paths dos blocos $6\times6$ do D&C. **Um switcher no grafo do cavalo é um objeto finito que você pode buscar por computador e depois fixar como constante da prova.** O menor $k$ que funcione (provavelmente $k = 2$ ou $3$, ou seja um 4-ciclo ou 6-ciclo no grafo do cavalo mais 1–2 caminhos) vira um lema com figura.

---

## 4. O passo (S3) já está resolvido na literatura — e você não sabia

Recipe do CNP, passo (S3): *"encontre um caminho hamiltoniano $H'$ de $v_1$ a $v_{k+1}$ em $G \setminus (V(W) \setminus \{v_1, v_{k+1}\})$."*

Para o grafo do cavalo isso é **exatamente** o teorema de Conrad–Hindrichs–Morsy–Wegener (que apareceu na bibliografia do KPW, §2.4 da auditoria):

> **[Conrad et al. 1994]** Para $n \geq 6$: existe caminho hamiltoniano do cavalo de $s$ a $t$ **sse** ($n$ par e $s,t$ têm cores diferentes) ou ($n$ ímpar e $s,t$ são ambos da cor dos cantos).

Caracterização **completa** de (S3) no tabuleiro cheio. Você precisa de duas coisas a mais:

1. **Controle de cor.** No switcher, $v_1$ e $v_{k+1}$ precisam ter cores opostas (para $n$ par). Como $C$ é um ciclo par e o grafo é bipartido, $v_1$ e $v_{k+1}$ estão a distância $k$ em $C$ — logo **cores opostas sse $k$ é ímpar**. ⚠️ Isso restringe $k$ a ímpar; verificar contra a construção do switcher, que exige $C$ par ($|C| = 2k$). $k=3$ ⟹ hexágono. Consistente com os seus $f_0, f_1$ serem hexágonos.
2. **Buraco.** Conrad et al. dão o tabuleiro cheio; você precisa de $G \setminus (\text{interior do switcher})$. Aqui entram: Miller & Farnsworth (tours com um quadrado removido, citado no Forrest–Teehan) e o seu próprio D&C por blocos, que já lida com fragmentos e perfis de fronteira.

Este é o ponto onde o **seu teorema de conexidade do bulk** entra: ele é exatamente a garantia de que remover um gadget finito do interior não desconecta o resto.

---

## 5. Bônus: por que o toro tem deficit 0 — há uma explicação na literatura

CNP citam:

> **[Alspach, Locke, Witte]** Se $\Gamma$ é grupo abeliano finito de ordem **ímpar**, então $\mathcal{C}_n(G) = \mathcal{C}(G)$ para todo Cayley graph conexo $G$ sobre $\Gamma$.

O grafo do cavalo **no toro** $T_{n,n}$ é um **Cayley graph sobre $\mathbb{Z}_n \times \mathbb{Z}_n$** (conjunto de conexão = os 8 vetores de movimento, simétrico). Isso é *exatamente* a estrutura que ALW tratam.

⚠️ ALW exigem ordem ímpar; $|\mathbb{Z}_n \times \mathbb{Z}_n| = n^2$ é par para os seus $n \in \{4,6,8,10\}$. Então **não se aplica diretamente** — mas o seu achado empírico `rank = β₁` (deficit 0, full rank) no toro é *moralmente* a mesma família de fenômeno: vértice-transitividade elimina o deficit.

**Ação de baixo custo e alto retorno:** buscar a referência ALW completa e testar computacionalmente o toro com $n$ **ímpar** (onde $n^2$ é ímpar e ALW **se aplica**). Se der deficit 0, você tem uma explicação *teórica citável* para a sua linha inteira de superfícies — o que resgata boa parte da §6 do dano do Forrest–Teehan. Isso também esclarece a estrutura real: **o deficit 3 é um fenômeno de fronteira, e a ausência de fronteira (toro) mais transitividade o elimina.**

---

## 6. Plano por fases

### Fase 0 — Reformulação dual *(1–2 dias, risco zero, valor garantido)*
- Reescrever a Conjectura de tightness na linguagem (C1)(C2)(C3).
- Identificar $\varphi_0,\varphi_1,\varphi_2$ explicitamente como witnesses $R$; exibir $|R|$ de cada.
- Escrever o dicionário da §0 numa observação do paper.
- ✅ **Entrega mesmo se todo o resto falhar.** Melhora o paper sozinho.

### Fase 1 — Lema 2.1 para bipartidos com $|V|$ par *(3–5 dias)*
- Provar a versão adaptada (§1). Cuidado com o argumento de maximalidade.
- **Formalizar em Lean** — encaixa no pipeline existente, é pequeno e autocontido.
- ✅ Lema publicável por si só; generaliza o CNP numa direção nova (bipartidos pares).

### Fase 2 — Reconhecimento computacional dos witnesses *(1–2 semanas)* ⭐ **fase crítica**
Antes de tentar provar qualquer coisa, **olhe para os $R$ de verdade**:
- Enumerar **todos** os $R$ maximais em $\mathcal{C}_n^\perp \setminus \mathcal{C}^\perp$ para $n=6$ (exaustivo — você tem os 9.862 tours) e amostrar para $n=8$.
- Medir: $|R|$; a densidade (C3) $e_R(A,B)/e_G(A,B)$ em cortes de vértice único; a distribuição de $R$ por nível $L$ (bulk vs. borda vs. canto).
- **Pergunta decisiva:** todo $R$ maximal é *equivalente módulo $X$* a um dos 3 conhecidos? Se sim para $n=6,8$, a tightness é verdadeira e você sabe o que provar. **Se aparecer um $R$ inesperado, a tightness pode ser FALSA para $n$ grande** — e descobrir isso agora vale mais que meses de tentativa de prova.
- Reaproveita direto: `forbidden_cycles_6x6/`, `deficit_theorem/`, `cycle_space_hyperplanes`.

### Fase 3 — Construção do switcher no bulk *(2–4 semanas)*
- Busca computacional pelo menor $R$-parity-switcher que evite os cantos, com $k$ ímpar (§4, controle de cor).
- Candidato natural: $k=3$ (hexágono) — note que os seus ciclos proibidos $f_0, f_1$ **são hexágonos**.
- Fixar como gadget constante; provar a propriedade dos dois caminhos hamiltonianos de paridades distintas (é verificação finita).
- Reaproveita: `patch_LUT`, catálogo D&C.

### Fase 4 — Completamento (S3) *(4–8 semanas, maior risco)*
- Conrad et al. dá o tabuleiro cheio; falta **evitar o buraco** do switcher.
- Combinar com conexidade do bulk (já provada e formalizada por você) + D&C por blocos.
- ⚠️ Aqui mora o risco real. É onde a rigidez do cavalo pode simplesmente travar.

### Fase 5 — Formalização *(depois, se 1–4 fecharem)*

---

## 7. Avaliação honesta de risco

**Chance de fechar tightness para todo $n$ par: baixa-a-moderada.** O gargalo é a Fase 4. CNP resolvem (S3) com maquinaria pseudoaleatória que você não tem, e a substituição construtiva em grafo de grau 8 com buraco é genuinamente difícil.

**Chance de resultado publicável: alta**, e não depende da Fase 4:

| Resultado | Depende de | Probabilidade |
|---|---|---|
| Reformulação dual + dicionário | Fase 0 | ~certa |
| Lema 2.1 para bipartidos pares (+ Lean) | Fase 1 | alta |
| Mapa computacional dos witnesses $n=6,8$ | Fase 2 | alta |
| Tightness para $n$ fixo, rigorosa | Fases 0–3 | média |
| Tightness geral $\forall n$ par | Fases 0–4 | baixa-média |

**Recomendação:** faça Fases 0–2 **antes de decidir** se vale investir nas 3–4. A Fase 2 é o teste decisivo, custa semanas e não meses, e usa código que você já tem. Se ela mostrar que todo witness maximal em $n=6,8$ colapsa nos 3 conhecidos, você tem evidência forte e sabe exatamente qual gadget construir. Se mostrar um witness exótico, você economizou meses e ganhou um resultado negativo interessante.

**Risco a monitorar:** a área está ativa (JCTB 2026, arXivs 2025–2026, Hou–Yin resolvendo o caso Dirac). O caso **esparso com deficit positivo exato** está em aberto e é visível. Vale alerta no arXiv (math.CO, termos `Hamilton space`, `Hamilton-generated`, `cycle space`).

---

## 8. Referências desta linha

```
✅[CNP2026]     M. Christoph, R. Nenadov, K. Petrova. The Hamilton space of pseudorandom
                  graphs. JCTB 176 (2026), 254-267.  [Lema 2.1, Def 2.2, Lema 2.3]
 [ALW]          B. Alspach, S. Locke, D. Witte. (Cayley graphs sobre grupos abelianos de
                  ordem impar: C_n(G)=C(G).)  <-- BUSCAR REF COMPLETA; ref [2] do CNP.
                  Relevante ao seu resultado de toro (§5).
✅[Heinig2013]  P. Heinig. When Hamilton circuits generate the cycle space of a random
                  graph. arXiv:1303.0026.   [necessidade de delta >= 3]
✅[HouYin2025]  X. Hou, Z. Yin. Dirac-type condition for Hamilton-generated graphs.
                  arXiv:2503.15950.
✅[Conrad1994]  A. Conrad, T. Hindrichs, H. Morsy, I. Wegener. Solution of the knight's
                  Hamiltonian path problem on a chessboard. DAM 50 (1994), 125-134.
                  <-- RESOLVE (S3) NO TABULEIRO CHEIO
 [MillerFarnsworth] Tours em cilindros/toros com um quadrado removido. <-- p/ o "buraco"
✅[Hartman1983] I. B.-A. Hartman. Long cycles generate the cycle space of a graph.
                  European J. Combin. 4 (1983), 237-246.
```

**Pendências desta nota:**
- [ ] Referência completa de Alspach–Locke–Witte (ref [2] do CNP — baixar a bibliografia).
- [ ] Confirmar que o argumento de maximalidade de (C3) não usa paridade (§1).
- [ ] Verificar o controle de cor $k$ ímpar contra a construção do switcher (§4.1).
