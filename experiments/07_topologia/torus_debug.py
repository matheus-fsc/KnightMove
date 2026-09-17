#!/usr/bin/env python3
"""
Torus 4×4 Knight Graph — Progressive Wrap Removal Sweep (v2)
=============================================================

Q measures independent XOR constraints from degree-2 vertices
projected into the CYCLE SPACE (kernel of ∂1 over GF(2)).

Key mathematical framework:
- A Hamiltonian cycle H ⊆ E has characteristic vector x ∈ F2^E
  satisfying ∂1·x = 0 (mod 2) and x lies in the cycle space.
- At a degree-2 vertex v with edges {e_a, e_b}, both must be used:
  x_{e_a} = 1 and x_{e_b} = 1.
- The constraint "e_a = 1" can be written as: x · δ_{e_a} = 1 (mod 2)
  where δ_{e_a} is the standard basis vector.
- Projecting δ_{e_a} into the cycle space quotient F2^E / rowspace(∂1),
  the constraint becomes a coset condition.
- Q counts independent such constraints.

Actually, the correct formulation for the "XOR clause" from vertex v of 
degree 2 is: the INDIVIDUAL edge indicator δ_{e_a} (or equivalently δ_{e_b},
since they differ by a row of ∂1).

So Q = rank( [∂1 ; δ_{e_a1} ; δ_{e_a2} ; ... ] ) - rank(∂1)
where {e_a_i} is one chosen forced edge per degree-2 vertex.

But wait — we should count constraints from vertices with degree ≤ 2 more
carefully. Let me reconsider from first principles.

For a Hamiltonian cycle candidate x ∈ {0,1}^E:
1) ∂1 · x ≡ 0 (mod 2) — even degree at every vertex in the subgraph
2) For each vertex v of degree d(v) in G: Σ_{e∋v} x_e = 2 
   (exactly 2 edges at each vertex — but this is NOT a GF(2) constraint,
   it's an integer constraint: Σ x_e = 2)

Over GF(2), constraint (2) reduces to: Σ_{e∋v} x_e ≡ 0 (mod 2),
which is exactly constraint (1). So the GF(2) system alone gives
us the cycle space, and the integer constraint "= 2" is additional.

For a degree-2 vertex v with edges {e_a, e_b}:
- Over GF(2): x_{e_a} + x_{e_b} ≡ 0  (from ∂1)
- Over Z: x_{e_a} + x_{e_b} = 2, so x_{e_a} = x_{e_b} = 1
- The constraint x_{e_a} = 1 over GF(2) is: e_a^T · x = 1

So the "obligation" from a degree-2 vertex is an AFFINE constraint
x_{e_a} = 1 (not a linear one). The question is whether the 
AFFINE SYSTEM {∂1·x = 0, x_{e_a} = 1 for each deg-2 vertex}
is consistent, and Q measures its deficiency.

But Q in the user's framework = dimension of independent XOR *clauses*.
Let me re-examine the definition more carefully.

The user says:
"Q = dim( π( Span( XORpairs das arestas obrigatórias ) ) )
 onde π é a projeção em F2^E / rowspace(∂1 restrita às active_edges)"
and "XOR-pairs só entre arestas de vértices DISTINTOS de grau 2"

XOR-pairs between edges of DISTINCT degree-2 vertices:
For degree-2 vertices u,v with forced edges e_u, e_v:
The XOR pair is δ_{e_u} + δ_{e_v}.

The constraint x_{e_u} = 1 and x_{e_v} = 1 implies
x_{e_u} + x_{e_v} = 0 over GF(2), i.e., (δ_{e_u} + δ_{e_v})^T · x = 0.

Wait, but x_{e_u} = 1 and x_{e_v} = 1 gives x_{e_u} + x_{e_v} = 0 mod 2.
That IS a linear constraint (not affine), and it's non-trivial only if
δ_{e_u} + δ_{e_v} is not in rowspace(∂1).

But for ANY pair of edges e_u, e_v of degree-2 vertices, we know both
must be 1 in a HC, so their sum must be 0 mod 2. The question is how
many of these pairwise XOR constraints are independent modulo ∂1.

If there are m degree-2 vertices with forced edges e_1,...,e_m,
the pairwise XOR constraints are:
  δ_{e_i} + δ_{e_j} for all i<j

The span of all pairwise sums of a set {v_1,...,v_m} equals the span 
of {v_1+v_2, v_1+v_3, ..., v_1+v_m} which has dimension ≤ m-1
(it's the "difference space").

Actually, Q should be computed as:
Q = rank([∂1; δ_{e_1}+δ_{e_2}; δ_{e_1}+δ_{e_3}; ... ; δ_{e_1}+δ_{e_m}]) - rank(∂1)

Let me think again about what happens for the plane 4×4...
In the plane 4×4, corners (0,0),(0,2),(1,3),(2,0),(3,1),(3,3) etc have degree 2.
Actually wait, let me check what the actual degree-2 vertices are for
the plane knight graph on 4×4.

Let me just verify: for 4×4 plane, the corner (0,0) has only moves to
(1,2) and (2,1) — degree 2. Similarly other corners.

Actually, the 4 corners of the 4×4 board are (0,0), (0,3), (3,0), (3,3)
and each has degree 2. But there may be other degree-2 vertices too.

Let me just fix the computation and run it.
"""

import numpy as np
from itertools import combinations, product
import random

# ============================================================
# 1. BUILD THE KNIGHT GRAPH ON THE 4×4 TORUS
# ============================================================

N = 4
KNIGHT_MOVES = [(1,2),(1,-2),(-1,2),(-1,-2),(2,1),(2,-1),(-2,1),(-2,-1)]

def vertex_id(r, c):
    return r * N + c

def build_torus_graph():
    """Build the knight graph on the 4×4 torus. Classify each edge."""
    edges = []
    edge_set = set()
    
    for r in range(N):
        for c in range(N):
            for dr, dc in KNIGHT_MOVES:
                r2_raw, c2_raw = r + dr, c + dc
                r2 = r2_raw % N
                c2 = c2_raw % N
                
                v1 = vertex_id(r, c)
                v2 = vertex_id(r2, c2)
                
                if v1 >= v2:
                    continue
                edge_key = (v1, v2)
                if edge_key in edge_set:
                    continue
                edge_set.add(edge_key)
                
                # Classify the edge
                r_wraps = (r2_raw != r2)
                c_wraps = (c2_raw != c2)
                
                if r_wraps and c_wraps:
                    cat = 'diag-wrap'
                elif r_wraps:
                    cat = 'y-wrap'
                elif c_wraps:
                    cat = 'x-wrap'
                else:
                    cat = 'internal'
                
                edges.append(((r, c), (r2, c2), v1, v2, cat))
    
    return edges

all_edges = build_torus_graph()

# Print category counts
categories = {}
for e in all_edges:
    cat = e[4]
    categories[cat] = categories.get(cat, 0) + 1

print("=" * 70)
print("TORUS 4×4 KNIGHT GRAPH — EDGE CLASSIFICATION")
print("=" * 70)
print(f"Total edges: {len(all_edges)}")
for cat in ['internal', 'x-wrap', 'y-wrap', 'diag-wrap']:
    print(f"  {cat:12s}: {categories.get(cat, 0)}")

# Build edge index
edge_to_idx = {}
for i, e in enumerate(all_edges):
    edge_to_idx[(e[2], e[3])] = i

# Identify wrap edges and their indices
wrap_edge_indices = []
wrap_edge_cats = []
for i, e in enumerate(all_edges):
    if e[4] != 'internal':
        wrap_edge_indices.append(i)
        wrap_edge_cats.append(e[4])

print(f"\nTotal wrap edges: {len(wrap_edge_indices)}")
for i, idx in enumerate(wrap_edge_indices):
    e = all_edges[idx]
    print(f"  [{i:2d}] edge#{idx:2d}  {e[0]}→{e[1]}  ({e[4]})")


# ============================================================
# 2. GF(2) RANK via Gaussian elimination
# ============================================================

def gf2_rank(M):
    """Compute rank of matrix M over GF(2)."""
    if M.size == 0 or M.shape[0] == 0 or M.shape[1] == 0:
        return 0
    M = M.copy().astype(np.int64) % 2
    rows, cols = M.shape
    pivot_row = 0
    for col in range(cols):
        found = -1
        for row in range(pivot_row, rows):
            if M[row, col] % 2 == 1:
                found = row
                break
        if found == -1:
            continue
        M[[pivot_row, found]] = M[[found, pivot_row]]
        for row in range(rows):
            if row != pivot_row and M[row, col] % 2 == 1:
                M[row] = (M[row] + M[pivot_row]) % 2
        pivot_row += 1
    return pivot_row


# ============================================================
# 3. COMPUTE Q
# ============================================================

def compute_Q(active_edge_mask):
    """
    Compute Q for the subgraph induced by active_edge_mask.
    
    For each degree-2 vertex v with edges {e_a, e_b}:
      - Both edges are forced in any Hamiltonian cycle: x_{e_a}=x_{e_b}=1
      - Pick one "representative" forced edge per vertex, say e_a.
    
    Method 1 (XOR-pairs between distinct deg-2 vertices):
      Generate δ_{e_i} + δ_{e_j} for all pairs (i,j) of deg-2 vertices.
      Q = rank([∂1; all XOR-pairs]) - rank(∂1)
    
    Method 2 (equivalent, simpler):
      Q = rank([∂1; δ_{e_1}; δ_{e_2}; ...; δ_{e_m}]) - rank(∂1) - (1 if m>0 else 0)
      No wait, let me think...
      
    Actually, the span of all pairwise XOR {δ_{e_i}+δ_{e_j} : i<j} equals
    the "difference subspace" of {δ_{e_1},...,δ_{e_m}}, which has dimension
    rank({δ_{e_i}}) - 1 when projected suitably. But we also need to account
    for ∂1.
    
    The correct formula:
      Let W = rowspace(∂1)
      Let S = Span(δ_{e_1},...,δ_{e_m}) (one forced edge per deg-2 vertex)
      Let S_pair = Span(δ_{e_i}+δ_{e_j} : i<j)
      
      Q = dim(π(S_pair)) where π: F2^E → F2^E / W
        = dim((S_pair + W) / W)
        = rank([∂1; all XOR-pairs]) - rank(∂1)
    
    Returns: (Q, deg2_vertex_coords, degree_dict)
    """
    num_active = int(np.sum(active_edge_mask))
    if num_active == 0:
        return 0, [], {}
    
    # Map active edges to contiguous indices
    active_indices = np.where(active_edge_mask)[0]
    active_map = {old: new for new, old in enumerate(active_indices)}
    E = num_active
    V = N * N
    
    # Build adjacency and compute degrees
    degree = np.zeros(V, dtype=int)
    adj = {v: [] for v in range(V)}
    
    for idx in active_indices:
        e = all_edges[idx]
        v1, v2 = e[2], e[3]
        degree[v1] += 1
        degree[v2] += 1
        adj[v1].append(idx)
        adj[v2].append(idx)
    
    # Find degree-2 vertices
    deg2_verts = [v for v in range(V) if degree[v] == 2]
    
    # Build incidence matrix ∂1 over GF(2): vertices × active_edges
    incidence = np.zeros((V, E), dtype=np.int64)
    for new_idx, old_idx in enumerate(active_indices):
        e = all_edges[old_idx]
        v1, v2 = e[2], e[3]
        incidence[v1, new_idx] = 1
        incidence[v2, new_idx] = 1
    
    rank_incidence = gf2_rank(incidence)
    
    if len(deg2_verts) < 2:
        # Need at least 2 deg-2 vertices to form XOR pairs
        deg2_coords = [(v // N, v % N) for v in deg2_verts]
        return 0, deg2_coords, {v: int(degree[v]) for v in range(V)}
    
    # For each degree-2 vertex, pick one forced edge (the first one)
    forced_edges = []
    for v in deg2_verts:
        e_a = adj[v][0]  # pick first edge
        forced_edges.append(e_a)
    
    # Build XOR-pair vectors: δ_{e_i} + δ_{e_j} for all pairs (i,j)
    m = len(forced_edges)
    xor_vectors = []
    for i in range(m):
        for j in range(i+1, m):
            vec = np.zeros(E, dtype=np.int64)
            vec[active_map[forced_edges[i]]] = 1
            vec[active_map[forced_edges[j]]] = 1
            xor_vectors.append(vec)
    
    xor_matrix = np.array(xor_vectors, dtype=np.int64)
    combined = np.vstack([incidence, xor_matrix])
    rank_combined = gf2_rank(combined)
    
    Q = rank_combined - rank_incidence
    
    deg2_coords = [(v // N, v % N) for v in deg2_verts]
    return Q, deg2_coords, {v: int(degree[v]) for v in range(V)}


# ============================================================
# Also implement a direct "affine deficiency" version for cross-check
# ============================================================

def compute_Q_affine(active_edge_mask):
    """
    Alternative computation:
    For each degree-2 vertex v with edges {e_a, e_b}, the HC constraint
    forces x_{e_a} = 1. Over GF(2), this is an affine constraint.
    
    We want: how many independent affine constraints do the forced edges
    impose on the cycle space?
    
    The cycle space is ker(∂1) over GF(2), dimension = E - rank(∂1).
    (assuming connected graph; for disconnected, adjust)
    
    Forced edge constraints: x_{e_i} = 1 for i=1..m
    
    These are affine hyperplanes in F2^E. Their intersection with ker(∂1)
    has codimension (in ker(∂1)) equal to the number of independent
    constraints among {δ_{e_i}} projected onto ker(∂1).
    
    This equals: rank([∂1; δ_{e_1}; ...; δ_{e_m}]) - rank(∂1)
    """
    num_active = int(np.sum(active_edge_mask))
    if num_active == 0:
        return 0, [], {}
    
    active_indices = np.where(active_edge_mask)[0]
    active_map = {old: new for new, old in enumerate(active_indices)}
    E = num_active
    V = N * N
    
    degree = np.zeros(V, dtype=int)
    adj = {v: [] for v in range(V)}
    
    for idx in active_indices:
        e = all_edges[idx]
        v1, v2 = e[2], e[3]
        degree[v1] += 1
        degree[v2] += 1
        adj[v1].append(idx)
        adj[v2].append(idx)
    
    deg2_verts = [v for v in range(V) if degree[v] == 2]
    
    incidence = np.zeros((V, E), dtype=np.int64)
    for new_idx, old_idx in enumerate(active_indices):
        e = all_edges[old_idx]
        v1, v2 = e[2], e[3]
        incidence[v1, new_idx] = 1
        incidence[v2, new_idx] = 1
    
    rank_inc = gf2_rank(incidence)
    
    if len(deg2_verts) == 0:
        deg2_coords = []
        return 0, deg2_coords, {v: int(degree[v]) for v in range(V)}
    
    # For each degree-2 vertex, pick one forced edge
    forced_edge_vecs = []
    for v in deg2_verts:
        e_a = adj[v][0]
        vec = np.zeros(E, dtype=np.int64)
        vec[active_map[e_a]] = 1
        forced_edge_vecs.append(vec)
    
    forced_matrix = np.array(forced_edge_vecs, dtype=np.int64)
    combined = np.vstack([incidence, forced_matrix])
    rank_combined = gf2_rank(combined)
    
    Q_affine = rank_combined - rank_inc
    
    deg2_coords = [(v // N, v % N) for v in deg2_verts]
    return Q_affine, deg2_coords, {v: int(degree[v]) for v in range(V)}


# ============================================================
# SANITY CHECK: Verify plane 4×4 first
# ============================================================

print("\n" + "=" * 70)
print("SANITY CHECK: Plane 4×4 Knight Graph")
print("=" * 70)

# Build plane graph directly
plane_edges = []
plane_edge_set = set()
for r in range(N):
    for c in range(N):
        for dr, dc in KNIGHT_MOVES:
            r2, c2 = r + dr, c + dc
            if 0 <= r2 < N and 0 <= c2 < N:
                v1 = vertex_id(r, c)
                v2 = vertex_id(r2, c2)
                if v1 < v2:
                    key = (v1, v2)
                    if key not in plane_edge_set:
                        plane_edge_set.add(key)
                        plane_edges.append((r, c, r2, c2, v1, v2))

print(f"Plane 4×4 edges: {len(plane_edges)}")

# Compute degrees
plane_deg = [0]*16
plane_adj = {v: [] for v in range(16)}
for i, (r1,c1,r2,c2,v1,v2) in enumerate(plane_edges):
    plane_deg[v1] += 1
    plane_deg[v2] += 1
    plane_adj[v1].append(i)
    plane_adj[v2].append(i)

print("Vertex degrees (plane 4×4):")
for r in range(N):
    for c in range(N):
        v = vertex_id(r, c)
        print(f"  ({r},{c}): degree = {plane_deg[v]}")

deg2_plane = [v for v in range(16) if plane_deg[v] == 2]
print(f"\nDegree-2 vertices: {[(v//N, v%N) for v in deg2_plane]}")

# Build incidence for plane graph
E_plane = len(plane_edges)
V = 16
inc_plane = np.zeros((V, E_plane), dtype=np.int64)
for i, (r1,c1,r2,c2,v1,v2) in enumerate(plane_edges):
    inc_plane[v1, i] = 1
    inc_plane[v2, i] = 1

rank_inc_plane = gf2_rank(inc_plane)
print(f"\nrank(∂1) for plane 4×4 = {rank_inc_plane}")
print(f"E = {E_plane}, so dim(cycle space) = E - rank(∂1) = {E_plane - rank_inc_plane}")

# Method: forced edge indicators  
if len(deg2_plane) > 0:
    forced_vecs = []
    for v in deg2_plane:
        e_idx = plane_adj[v][0]
        vec = np.zeros(E_plane, dtype=np.int64)
        vec[e_idx] = 1
        forced_vecs.append(vec)
        print(f"  Forced edge for vertex ({v//N},{v%N}): edge#{e_idx} = {plane_edges[e_idx][:4]}")
    
    forced_mat = np.array(forced_vecs, dtype=np.int64)
    combined = np.vstack([inc_plane, forced_mat])
    rank_combined = gf2_rank(combined)
    Q_plane_affine = rank_combined - rank_inc_plane
    print(f"\nQ (affine, individual forced edges) = {Q_plane_affine}")
    
    # XOR pairs method
    xor_vecs = []
    for i in range(len(deg2_plane)):
        for j in range(i+1, len(deg2_plane)):
            vec = np.zeros(E_plane, dtype=np.int64)
            e_i = plane_adj[deg2_plane[i]][0]
            e_j = plane_adj[deg2_plane[j]][0]
            vec[e_i] = 1
            vec[e_j] = 1
            xor_vecs.append(vec)
    
    xor_mat = np.array(xor_vecs, dtype=np.int64)
    combined2 = np.vstack([inc_plane, xor_mat])
    rank_combined2 = gf2_rank(combined2)
    Q_plane_xor = rank_combined2 - rank_inc_plane
    print(f"Q (XOR pairs) = {Q_plane_xor}")

# Now verify using the torus graph with all wraps removed
print(f"\n--- Cross-check via torus framework ---")
all_wrap_idx = [i for i, e in enumerate(all_edges) if e[4] != 'internal']
mask_plane = np.ones(len(all_edges), dtype=bool)
for idx in all_wrap_idx:
    mask_plane[idx] = False

Q_torus_xor, deg2_torus, _ = compute_Q(mask_plane)
Q_torus_aff, deg2_torus2, _ = compute_Q_affine(mask_plane)
print(f"Q (XOR pairs, via torus framework) = {Q_torus_xor}")
print(f"Q (affine, via torus framework) = {Q_torus_aff}")
print(f"Deg-2 vertices: {deg2_torus}")

# Check: are the plane edges the same as torus-minus-wraps?
torus_internal = [(e[2],e[3]) for e in all_edges if e[4] == 'internal']
plane_eset = [(e[4],e[5]) for e in plane_edges]
print(f"\nTorus internal edges: {len(torus_internal)}")
print(f"Plane edges: {len(plane_eset)}")
print(f"Same set? {set(torus_internal) == set(plane_eset)}")

# ============================================================
# DEEP DEBUG: check if forced edges are in rowspace of ∂1
# ============================================================
print("\n--- Debug: are forced edge vectors in rowspace(∂1)? ---")

for v in deg2_plane:
    e_idx = plane_adj[v][0]
    vec = np.zeros(E_plane, dtype=np.int64)
    vec[e_idx] = 1
    
    # Check if vec ∈ rowspace(inc_plane)
    test = np.vstack([inc_plane, vec.reshape(1,-1)])
    r_test = gf2_rank(test)
    in_rowspace = (r_test == rank_inc_plane)
    print(f"  δ_{{e_{e_idx}}} for vertex ({v//N},{v%N}): in rowspace(∂1)? {in_rowspace}")

# XOR pairs
print("\n--- Debug: are XOR pair vectors in rowspace(∂1)? ---")
for i in range(len(deg2_plane)):
    for j in range(i+1, len(deg2_plane)):
        e_i = plane_adj[deg2_plane[i]][0]
        e_j = plane_adj[deg2_plane[j]][0]
        vec = np.zeros(E_plane, dtype=np.int64)
        vec[e_i] = 1
        vec[e_j] = 1
        
        test = np.vstack([inc_plane, vec.reshape(1,-1)])
        r_test = gf2_rank(test)
        in_rs = (r_test == rank_inc_plane)
        vi, vj = deg2_plane[i], deg2_plane[j]
        print(f"  δ_{{e_{e_i}}}+δ_{{e_{e_j}}} for verts ({vi//N},{vi%N}),({vj//N},{vj%N}): in rowspace? {in_rs}")

print("\n✓ Sanity check complete.")
