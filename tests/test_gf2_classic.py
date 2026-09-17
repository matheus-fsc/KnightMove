import numpy as np

def rank_gf2(matrix):
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

def calc_q(N, M, remove_corners=False):
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
    
    matrix = []
    for u, v in edges:
        row = [0] * V
        row[v_idx[u]] = 1
        row[v_idx[v]] = 1
        matrix.append(row)
        
    rank = rank_gf2(matrix)
    nullity = E - rank
    beta_1 = E - V + 1
    Q = beta_1 - nullity
    
    print(f"Board {N}x{M} (remove_corners={remove_corners}): V={V}, E={E}, beta_1={beta_1}, rank={rank}, nullity={nullity}, Q={Q}")

calc_q(8, 8, False)
calc_q(8, 8, True)
