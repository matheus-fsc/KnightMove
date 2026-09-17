#!/usr/bin/env python3
"""
Verification of Conjecture 6.1 — Rectangular Knight's Tour
==========================================================
Checks corner degrees and homological obstruction Q on rectangular boards.
"""

import numpy as np
import time

def gf2_rank(M):
    if M.size == 0 or M.shape[0] == 0 or M.shape[1] == 0:
        return 0
    M = M.copy().astype(np.int64) % 2
    rows, cols = M.shape
    pivot_row = 0
    for col in range(cols):
        found = -1
        for row in range(pivot_row, rows):
            if M[row, col] & 1:
                found = row
                break
        if found == -1:
            continue
        M[[pivot_row, found]] = M[[found, pivot_row]]
        for row in range(rows):
            if row != pivot_row and M[row, col] & 1:
                M[row] ^= M[pivot_row]
        pivot_row += 1
    return pivot_row

def build_planar_graph(n, m):
    edges = set()
    KNIGHT_MOVES = [(1,2),(1,-2),(-1,2),(-1,-2),(2,1),(2,-1),(-2,1),(-2,-1)]
    for r in range(n):
        for c in range(m):
            for dr, dc in KNIGHT_MOVES:
                r2, c2 = r + dr, c + dc
                if 0 <= r2 < n and 0 <= c2 < m:
                    v1 = r * m + c
                    v2 = r2 * m + c2
                    edges.add((min(v1, v2), max(v1, v2)))
    return sorted(list(edges))

def compute_Q(n, m):
    edges = build_planar_graph(n, m)
    E = len(edges)
    V = n * m
    
    # Calculate degree
    degree = np.zeros(V, dtype=int)
    adj = {v: [] for v in range(V)}
    for idx, (v1, v2) in enumerate(edges):
        degree[v1] += 1
        degree[v2] += 1
        adj[v1].append(idx)
        adj[v2].append(idx)
        
    corner_vids = [0, m-1, (n-1)*m, n*m-1]
    corner_degrees = [degree[vid] for vid in corner_vids]
    deg2_count = sum(1 for d in corner_degrees if d == 2)
    
    # Calculate Q
    inc = np.zeros((V, E), dtype=np.int64)
    for idx, (v1, v2) in enumerate(edges):
        inc[v1, idx] = 1
        inc[v2, idx] = 1
    rank_inc = gf2_rank(inc)
    
    deg2_corners = [vid for vid in corner_vids if degree[vid] == 2]
    
    if len(deg2_corners) < 2:
        return corner_degrees, deg2_count, 0
        
    forced = [adj[v][0] for v in deg2_corners]
    xvecs = []
    for i in range(len(forced)):
        for j in range(i+1, len(forced)):
            vec = np.zeros(E, dtype=np.int64)
            vec[forced[i]] = 1
            vec[forced[j]] = 1
            xvecs.append(vec)
            
    combined = np.vstack([inc, np.array(xvecs, dtype=np.int64)])
    Q = gf2_rank(combined) - rank_inc
    return corner_degrees, deg2_count, Q

def run_tests():
    start_time = time.time()
    
    conjecture_cases = [(4,6),(4,8),(4,10),(5,6),(5,8),(6,8),(6,10),(6,12),(8,10),(8,12)]
    boundary_cases = [(3,4),(3,6),(3,8),(3,10)]
    square_cases = [(4,4),(6,6),(8,8),(10,10)]
    extreme_cases = [(4,20),(4,50),(6,20),(6,50)]
    
    all_cases = [
        ("Conjecture (min >= 4)", conjecture_cases),
        ("Boundary (min = 3)", boundary_cases),
        ("Square (control)", square_cases),
        ("Extreme (large ratio)", extreme_cases)
    ]
    
    results = []
    
    print("| Grupo | n | m | grau(0,0) | grau(0,m-1) | grau(n-1,0) | grau(n-1,m-1) | deg2_count | Q |")
    print("|---|---|---|---|---|---|---|---|---|")
    
    for group_name, cases in all_cases:
        for n, m in cases:
            c_degs, d2c, Q = compute_Q(n, m)
            results.append((group_name, n, m, c_degs, d2c, Q))
            print(f"| {group_name} | {n} | {m} | {c_degs[0]} | {c_degs[1]} | {c_degs[2]} | {c_degs[3]} | {d2c} | {Q} |")
            
    total_time = time.time() - start_time
    print(f"\nTempo total de execução: {total_time:.4f} segundos\n")
    
    # Evaluate Hypotheses
    h_grau = True
    h_q = True
    h_min = False
    h_ratio = True
    
    q_vals_min_ge_4 = []
    
    for grp, n, m, c_degs, d2c, Q in results:
        min_dim = min(n, m)
        if min_dim >= 4:
            q_vals_min_ge_4.append(Q)
            if d2c != 4:
                h_grau = False
                print(f"Falha H_GRAU no caso ({n},{m}): deg2_count = {d2c}")
            if Q != 3:
                h_q = False
                print(f"Falha H_Q no caso ({n},{m}): Q = {Q}")
        elif min_dim == 3:
            if Q != 3:
                h_min = True
                
    # Check if Q is independent of n/m ratio for min>=4
    if len(set(q_vals_min_ge_4)) != 1 or list(set(q_vals_min_ge_4))[0] != 3:
        h_ratio = False
        
    print("=" * 50)
    print("STATUS DAS HIPÓTESES:")
    print("=" * 50)
    print(f"H_GRAU  : {'VERDADEIRO' if h_grau else 'FALSO'}")
    print(f"H_Q     : {'VERDADEIRO' if h_q else 'FALSO'}")
    print(f"H_MIN   : {'VERDADEIRO' if h_min else 'FALSO'}")
    print(f"H_RATIO : {'VERDADEIRO' if h_ratio else 'FALSO'}")
    print("=" * 50)

if __name__ == "__main__":
    run_tests()
