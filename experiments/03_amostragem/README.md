# 03 — Amostragem massiva

A ideia que o autor chamou de "paridade": fixar um par (início, fim) e amostrar
o tabuleiro em vez de enumerá-lo. Isso compacta o espaço o suficiente para
procurar padrões em tabuleiros grandes.

| diretório | fase |
|---|---|
| `cavalo_8x8/` | **Fase A** (validação em 6×6) e **Fase B** (re-coleta 8×8 sem viés) |
| `cavalo_10x10/` | **Fase C**: 500k amostras × 8 simetrias D4 = 4M tours |

Resultados:

- **Fase A** — Z3 com `sym=False` validado contra ground truth (`rho = 0.96`,
  7 de 7 órbitas); com `sym=True` a correlação cai para 0.89, o que revelou o viés.
- **Fase B** — 100% dos top-13 preservados após a re-coleta sem viés; média de
  13,2 órbitas por pool. O par `B8-D7 ↔ B8-C6` é top-1 nos 18 pools.
- **Fase C** — 9 pools, média 18,6 órbitas, compactação 7,94x. O par
  `B10-D9 ↔ B10-C8` é top-1 universal, confirmando a hipótese em 10×10.

As amostras brutas (2,4 GB) estão em `data/archives/` via Git LFS.
