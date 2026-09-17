# formalization

Formalização em **Lean 4** dos resultados que sobreviveram à verificação.

Estado: `Q(n) = 3` (quadrado) e `Q(n,m) = 3` (retangular), ambos para `n,m >=
6`, provados **sem `sorry`**. Os módulos centrais são `QAbstract`,
`RectQAbstract`, `RectBulkConnectivity` e `EdgewiseConst`.

Os `sorry` que aparecem em `computeQ` são sombras de teoremas já superados
pelas versões abstratas e não sustentam nenhum resultado citado.

O diretório de build `.lake/` (5,6 GB) não é versionado.

Ver `lean/README.md` para a lista de teoremas.
