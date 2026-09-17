
# --- bootstrap: torna core/ importavel a partir de qualquer subpasta ---
import sys as _sys, pathlib as _pathlib
for _anc in _pathlib.Path(__file__).resolve().parents:
    if (_anc / "core").is_dir() and (_anc / "pyproject.toml").exists():
        _sys.path.insert(0, str(_anc / "core"))
        break
# --- fim do bootstrap ---
import numpy as np
import knight_tours_hybrid as kt

def find_open_path(n, u, v):
    """Encontra um Caminho Hamiltoniano Aberto de u até v em um tabuleiro nxn."""
    ctx = kt.build_graph(n)
    V = ctx['V']
    E = ctx['E']
    
    # Verifica se (u, v) já é uma aresta válida de cavalo
    edge_idx = -1
    for ei in range(E):
        a, b = ctx['edge_endpoints'][ei]
        if (a == u and b == v) or (a == v and b == u):
            edge_idx = ei
            break
            
    # Se não for, adicionamos uma Aresta Virtual
    if edge_idx == -1:
        new_ei = E
        ctx['E'] = E + 1
        ctx['edge_endpoints'] = np.vstack([ctx['edge_endpoints'], np.array([[min(u, v), max(u, v)]], dtype=np.int32)])
        
        ctx['adj_edges'][u] = np.append(ctx['adj_edges'][u], new_ei)
        ctx['adj_edges'][v] = np.append(ctx['adj_edges'][v], new_ei)
        
        ctx['total_incident'][u] += 1
        ctx['total_incident'][v] += 1
        
        edge_idx = new_ei
        E = ctx['E']
        
    state = kt.State(V, E)
    
    # Preparar a fila de propagação R2
    queue = [(edge_idx, 1)] # Força a aresta U-V a ser ATIVA
    
    for w in range(V):
        if ctx['total_incident'][w] == 2:
            for ei in ctx['adj_edges'][w]:
                queue.append((int(ei), 1))
                
    status = kt._process_queue(state, ctx, queue)
    if status in (kt.CONTRADICTION, kt.SUBTOUR):
        return None
        
    rng = np.random.default_rng(42)
    tours = []
    
    if status == kt.COMPLETE_TOUR:
        if kt._is_complete_tour(state, ctx):
            tours.append(kt._extract_tour(state, ctx))
    else:
        kt._backtrack(state, ctx, rng, tours, 1)
        
    if not tours:
        return None
        
    tour = list(tours[0])
    idx_u = tour.index(u)
    
    idx_prev = (idx_u - 1) % V
    idx_next = (idx_u + 1) % V
    
    if tour[idx_prev] == v:
        # tour ... -> v -> u -> ... (desejado: começar em u, ir pra direita)
        open_path = tour[idx_u:] + tour[:idx_u]
    elif tour[idx_next] == v:
        # tour ... -> u -> v -> ... (inverter)
        tour_rev = tour[::-1]
        idx_u_rev = tour_rev.index(u)
        open_path = tour_rev[idx_u_rev:] + tour_rev[:idx_u_rev]
    else:
        raise RuntimeError("Aresta Virtual não foi usada!")
        
    return open_path

if __name__ == '__main__':
    n = 6
    # Coordenadas no formato v = row * n + col
    # Bloco 0,0: I=(5,1), E=(4,5)
    u1 = 5*n + 1; v1 = 4*n + 5
    # Bloco 0,1: I=(2,0), E=(5,4)
    u2 = 2*n + 0; v2 = 5*n + 4
    # Bloco 1,1: I=(0,2), E=(3,0)
    u3 = 0*n + 2; v3 = 3*n + 0
    # Bloco 1,0: I=(5,5), E=(0,3)
    u4 = 5*n + 5; v4 = 0*n + 3
    
    print("Testando busca de Caminho Aberto nos blocos 6x6...")
    for i, (u, v) in enumerate([(u1, v1), (u2, v2), (u3, v3), (u4, v4)]):
        path = find_open_path(n, u, v)
        if path:
            print(f"Bloco {i+1}: SUCESSO! ({u} -> {v}) | Nós: {len(path)}")
        else:
            print(f"Bloco {i+1}: FALHOU! ({u} -> {v}) não possui caminho hamiltoniano.")
