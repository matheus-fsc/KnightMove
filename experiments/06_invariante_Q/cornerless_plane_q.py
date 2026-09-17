import z3

def rank_gf2(matrix):
    """Calcula o posto (rank) de uma matriz sobre o corpo finito GF(2) via Eliminação Gaussiana."""
    if not matrix: return 0
    rows = []
    for r in matrix:
        val = 0
        for bit in r:
            val = (val << 1) | bit
        rows.append(val)
    rank = 0
    while rows:
        pivot_row = max(rows)
        if pivot_row == 0: break
        rank += 1
        pivot_bit = pivot_row.bit_length() - 1
        new_rows = []
        for r in rows:
            if (r >> pivot_bit) & 1:
                r ^= pivot_row
            if r > 0:
                new_rows.append(r)
        rows = new_rows
    return rank

def cornerless_plane_q(N, M, remove_corners=False):
    """
    Testa o colapso geométrico do Q.
    1. Constrói o grafo do cavalo (sem as quinas se solicitado).
    2. Usa o solver Z3 para resolver o sistema do 2-Fator e encontrar as arestas 
       completamente fixadas (0 graus de liberdade, por exemplo forçadas a 1 ou 0).
    3. Extrai as arestas "livres" e constrói a Matriz de Incidência sobre GF(2).
    4. Computa o Nullity sobre GF(2) do subespaço livre, revelando as dimensões reais perdidas.
    """
    corners = [(0,0), (0,M-1), (N-1,0), (N-1,M-1)]
    vertices = [(x,y) for x in range(N) for y in range(M)]
    if remove_corners:
        vertices = [v for v in vertices if v not in corners]
    
    v_idx = {v: i for i, v in enumerate(vertices)}
    V = len(vertices)
    
    edges = []
    for x, y in vertices:
        for dx, dy in [(1,2), (2,1), (-1,2), (-2,1), (1,-2), (2,-1), (-1,-2), (-2,-1)]:
            nx, ny = x + dx, y + dy
            if (nx, ny) in v_idx:
                if (x,y) < (nx, ny):
                    edges.append(((x,y), (nx,ny)))
    
    E = len(edges)
    
    # Montamos as equações do 2-Fator no Z3 para encontrar arestas rigidamente fixadas
    solver = z3.Solver()
    edge_vars = [z3.Int(f"e_{i}") for i in range(E)]
    for var in edge_vars:
        solver.add(var >= 0, var <= 1)
        
    for v in vertices:
        incident = [edge_vars[i] for i, e in enumerate(edges) if v in e]
        solver.add(z3.Sum(incident) == 2)
        
    free_edges = []
    for i, var in enumerate(edge_vars):
        # Verifica se pode assumir 0
        solver.push()
        solver.add(var == 0)
        can_be_0 = (solver.check() == z3.sat)
        solver.pop()
        
        # Verifica se pode assumir 1
        solver.push()
        solver.add(var == 1)
        can_be_1 = (solver.check() == z3.sat)
        solver.pop()
        
        if can_be_0 and can_be_1:
            free_edges.append(edges[i])
            
    # Matriz de incidência APENAS das arestas com graus de liberdade reais (não-fixadas)
    matrix = []
    for u, v in free_edges:
        row = [0] * V
        row[v_idx[u]] = 1
        row[v_idx[v]] = 1
        matrix.append(row)
        
    rank = rank_gf2(matrix)
    nullity = len(free_edges) - rank  # Nullity do subespaço ciclo-ativo
    beta_1 = E - V + 1
    Q = beta_1 - nullity
    
    print(f"[{'Plano Sem Quinas' if remove_corners else 'Plano Clássico'}] Tabuleiro {N}x{M}:")
    print(f"  V={V}, E={E}")
    print(f"  beta_1 (E - V + 1) = {beta_1}")
    print(f"  Arestas Totais = {E} | Arestas Livres = {len(free_edges)} (Arestas Fixadas pela Borda = {E - len(free_edges)})")
    print(f"  Rank GF(2) (sobre arestas livres) = {rank}")
    print(f"  Nullity (Graus de liberdade reais do 2-fator) = {nullity}")
    print(f"  --> Restrição Topológica/Geométrica Q = {Q}\n")

if __name__ == "__main__":
    print("=== TAREFA 1: O Plano Sem Cantos (Origens da Obstrução Q) ===")
    cornerless_plane_q(8, 8, remove_corners=False)
    cornerless_plane_q(8, 8, remove_corners=True)
