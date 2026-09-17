# 09. Transferência para outros domínios

Se o espaço de ciclos GF(2) organiza o passeio do cavalo, ele deveria ajudar
em outros problemas de enumeração de estruturas conexas. Testamos quatro.
Todos negativos, e todos pelo mesmo motivo.

| diretório | domínio | veredito |
|---|---|---|
| `tsp_cycle_space/` | caixeiro viajante | XOR de ciclos fundamentais da MST fica ~15% pior que 2-opt |
| `pathfinding_xor_experiment/` | pathfinding em grades esparsas | compatibilidade ~18% (contra 3% no TSP), cai com `n`; falha dominante é **desconexão (90,8%)**, não grau |
| `path_decomposition_experiment/` | caminhos `s -> t` no próprio grafo do cavalo | XOR de ciclos fundamentais cobre 0,0086% dos caminhos (4 de 46666 no 6×6); satura em ~10 no 8×8; backtracking com poda de grau e conectividade vence |
| `protein_study/` | busca de sítio ativo em proteínas (GASS) | **inviável**: as soluções são k-cliques transversais sob limiar métrico, não subgrafos pares. 0 de 17.367 XORs válidos |
| `complex_orbit/` | invariantes sobre C (DFT, winding, CTQW) no grafo do cavalo | `k=3` domina a DFT; orbit-parity tem rank 4, não 7 |

**O padrão:** GF(2) captura perfeitamente a condição de grau par. Não captura
conectividade. Em todo domínio onde a solução precisa ser conexa, o XOR gera
candidatos que satisfazem grau mas se quebram em componentes, e por isso a taxa de
acerto desaba. Essa é a mesma obstrução que limita o método no passeio do
cavalo.
