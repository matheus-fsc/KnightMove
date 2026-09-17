# Formalização em Lean 4

Os resultados que sobreviveram à verificação computacional foram formalizados.

## Estado

**Sem `sorry`:**

- `Q(n) = 3` para tabuleiros quadrados, `n >= 6`;
- `Q(n, m) = 3` para tabuleiros retangulares, `n, m >= 6`.

Módulos centrais: `QAbstract`, `RectQAbstract`, `RectBulkConnectivity`,
`EdgewiseConst`.

Outros módulos: `CycleSpaceDim`, `PuncGraph`, `PuncCount`, `Reduction`,
`CornerFiber`, `FaceCycles`, `GF2Path`, `ConjectureXOR`.

## O que a formalização NÃO estabelece

Esta seção importa mais que a anterior, e o paper a enuncia explicitamente:

> A cadeia formalizada é a **equivalência**, não a *tightness*.

O teorema `reduction` tem a forma

```
Z_bulk(n) <= H  ∧  tau ∈ H  ∧  dim H <= beta_1 - 3   ==>   dim H = beta_1 - 3
```

e a hipótese `Z_bulk(n) <= H` **é o problema inteiro**. O valor está em isolar
e verificar a implicação, não em descarregar a premissa.

A ponte com o objeto combinatório — definir `Ham(n)` como o span dos
indicadores de ciclos hamiltonianos de `G_n` e descarregar as hipóteses da
redução — permanece como comentário documentado, **não** como declaração.

### Uma escolha deliberada

O projeto optou por **não** escrever `theorem tightness := sorry`. Seria
trivial e daria aparência de progresso, mas colocaria uma afirmação não
provada no namespace, onde qualquer resultado construído sobre ela herdaria o
`sorryAx` sem aviso.

### Independência

A tightness por certificados (ver [[Invariante-Q]]) **não depende** da
formalização: sua prova é construtiva e cada certificado é auditado
independentemente. A formalização é contribuição paralela, não suporte lógico
daquele teorema.

## Outras ressalvas

- Os `sorry` que aparecem em `computeQ` são **sombras** de teoremas já
  superados pelas versões abstratas. Nenhum resultado citado depende deles.
- A formulação por `CornerFiber` foi escolhida justamente para evitar o
  quociente `F_2^E / row(∂_1)`, que é caro de formalizar.

## O que falta

A conexidade do bulk (`n >= 6`, indução `n -> n+2`) tem prova de papel completa
e verificação computacional até `n = 30`, mas ainda não foi formalizada. É a
próxima etapa natural.

**Código:** `formalization/lean/`
