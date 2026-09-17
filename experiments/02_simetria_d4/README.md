# 02 — Simetria diedral D4

O grafo do cavalo no tabuleiro quadrado é invariante sob o grupo diedral D4
(4 rotações × 2 reflexões). Fixar um representante por órbita reduz o espaço
de busca por um fator de até 8.

Experimentos aqui testam ordenação de ramos por órbita (`d4_branch_order`,
`d4_orbital_branching`): a simetria é usada não só para deduplicar resultados
mas para escolher a ordem de exploração.

**Veredito:** a redução de 8x é real e é usada em toda a amostragem posterior
(Fase A/B/C), mas é um fator constante — não muda o crescimento. Foi, porém,
o que tornou a amostragem 10×10 viável na prática.
