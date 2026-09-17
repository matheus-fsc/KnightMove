#!/usr/bin/env python3
"""
Torus 4×4 Knight Graph — Progressive Wrap Removal Sweep
========================================================
Computes Q = dim(π(Span(XOR-pairs))) for the homological invariant.

Key insight: On the 4×4 torus, many "wrap" edges also exist as plane edges
(reached by non-wrapping moves). Only edges exclusive to the torus are
"removable wraps". It turns out ALL 8 removable edges are diagonal-wraps
(both row and column wrap simultaneously).
"""

import numpy as np
from itertools import combinations
import random

N = 4
KNIGHT_MOVES = [(1,2),(1,-2),(-1,2),(-1,-2),(2,1),(2,-1),(-2,1),(-2,-1)]

def vid(r, c):
    return r * N + c

# ============================================================
# 1. BUILD TORUS AND PLANE GRAPHS, CLASSIFY EDGES
# ============================================================

# Build torus edges
torus_edge_set = {}  # (v1,v2) → list of (wraps, dr, dc, r1, c1, r2, c2)
for r in range(N):
    for c in range(N):
        for dr, dc in KNIGHT_MOVES:
            r2 = (r + dr) % N
            c2 = (c + dc) % N
            wraps = ((r + dr) != r2) or ((c + dc) != c2)
            v1, v2 = vid(r, c), vid(r2, c2)
            if v1 == v2:
                continue
            key = (min(v1,v2), max(v1,v2))
            if key not in torus_edge_set:
                torus_edge_set[key] = []
            torus_edge_set[key].append((wraps, dr, dc, r, c, r2, c2))

# Build plane edges
plane_edge_set = set()
for r in range(N):
    for c in range(N):
        for dr, dc in KNIGHT_MOVES:
            r2, c2 = r + dr, c + dc
            if 0 <= r2 < N and 0 <= c2 < N:
                v1, v2 = vid(r, c), vid(r2, c2)
                plane_edge_set.add((min(v1,v2), max(v1,v2)))

# Classify: "removable wrap" = torus-only edge; "internal" = also exists in plane
edges_to_remove = set(torus_edge_set.keys()) - plane_edge_set

all_edges = []  # ((r1,c1),(r2,c2), v1, v2, category)
for (v1, v2), moves in sorted(torus_edge_set.items()):
    r1, c1, r2, c2 = v1//N, v1%N, v2//N, v2%N
    if (v1, v2) in edges_to_remove:
        # Classify wrap sub-type
        wrap_types = set()
        for w, dr, dc, *_ in moves:
            rr, cc = r1 + dr, c1 + dc  # raw coords
            r_w = (rr % N != rr)
            c_w = (cc % N != cc)
            if r_w and c_w: wrap_types.add('diag-wrap')
            elif r_w: wrap_types.add('y-wrap')
            elif c_w: wrap_types.add('x-wrap')
        cat = 'diag-wrap' if 'diag-wrap' in wrap_types else (
              'y-wrap' if 'y-wrap' in wrap_types else (
              'x-wrap' if 'x-wrap' in wrap_types else 'wrap-other'))
    else:
        cat = 'internal'
    all_edges.append(((r1,c1),(r2,c2), v1, v2, cat))

# Count by category
cats = {}
for e in all_edges:
    cats[e[4]] = cats.get(e[4], 0) + 1

print("=" * 70)
print("TORUS 4×4 KNIGHT GRAPH — EDGE CLASSIFICATION")
print("=" * 70)
print(f"Total torus edges: {len(all_edges)}  |  Plane edges: {len(plane_edge_set)}")
print(f"Removable wrap edges: {len(edges_to_remove)}")
for cat in ['internal', 'x-wrap', 'y-wrap', 'diag-wrap', 'wrap-other']:
    if cat in cats:
        print(f"  {cat:12s}: {cats[cat]}")

wrap_indices = []
wrap_cats = []
for i, e in enumerate(all_edges):
    if e[4] != 'internal':
        wrap_indices.append(i)
        wrap_cats.append(e[4])
        print(f"  wrap[{len(wrap_indices)-1}] edge#{i:2d}  {e[0]}→{e[1]}  ({e[4]})")


# ============================================================
# 2. GF(2) RANK
# ============================================================

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


# ============================================================
# 3. COMPUTE Q (XOR-pairs of distinct degree-2 vertices)
# ============================================================

def compute_Q(active_mask):
    """
    Q = rank([∂1; XOR-pairs]) - rank(∂1)
    
    XOR-pairs: for degree-2 vertices u,v with forced edges e_u, e_v,
    the vector δ_{e_u} + δ_{e_v} encodes the constraint
    x_{e_u} ≡ x_{e_v} (mod 2) that any Hamiltonian cycle must satisfy.
    
    Q counts how many such constraints are independent in the quotient
    F2^E / rowspace(∂1), i.e., independent modulo the cycle space relations.
    """
    active_idx = np.where(active_mask)[0]
    E = len(active_idx)
    V = N * N
    if E == 0:
        return 0, [], {}
    
    amap = {old: new for new, old in enumerate(active_idx)}
    
    degree = np.zeros(V, dtype=int)
    adj = {v: [] for v in range(V)}
    for idx in active_idx:
        e = all_edges[idx]
        v1, v2 = e[2], e[3]
        degree[v1] += 1
        degree[v2] += 1
        adj[v1].append(idx)
        adj[v2].append(idx)
    
    deg2 = [v for v in range(V) if degree[v] == 2]
    
    # Incidence matrix ∂1
    inc = np.zeros((V, E), dtype=np.int64)
    for new, old in enumerate(active_idx):
        e = all_edges[old]
        inc[e[2], new] = 1
        inc[e[3], new] = 1
    
    rank_inc = gf2_rank(inc)
    
    if len(deg2) < 2:
        return 0, [(v//N, v%N) for v in deg2], {v: int(degree[v]) for v in range(V)}
    
    # One forced edge per degree-2 vertex
    forced = [adj[v][0] for v in deg2]
    
    # XOR pairs: δ_{e_i} + δ_{e_j} for all i < j
    xvecs = []
    for i in range(len(forced)):
        for j in range(i+1, len(forced)):
            vec = np.zeros(E, dtype=np.int64)
            vec[amap[forced[i]]] = 1
            vec[amap[forced[j]]] = 1
            xvecs.append(vec)
    
    combined = np.vstack([inc, np.array(xvecs, dtype=np.int64)])
    Q = gf2_rank(combined) - rank_inc
    
    return Q, [(v//N, v%N) for v in deg2], {v: int(degree[v]) for v in range(V)}


# ============================================================
# 4. SWEEP
# ============================================================

def make_mask(removed):
    mask = np.ones(len(all_edges), dtype=bool)
    for idx in removed:
        mask[idx] = False
    return mask

print("\n" + "=" * 70)
print("MAIN SWEEP RESULTS")
print("=" * 70)

results = []

# (a) Toro puro
Q, d2, dg = compute_Q(make_mask([]))
deg_set = sorted(set(dg.values()))
results.append(('toro puro', 0, len(d2), Q, f'todos grau {deg_set}'))
print(f"\n(a) Toro puro: Q={Q}, deg2={len(d2)}, graus={deg_set}")

# Category-based removals
def cat_idx(cat):
    return [i for i, e in enumerate(all_edges) if e[4] == cat]

configs = [
    ('só x-wraps', ['x-wrap']),
    ('só y-wraps', ['y-wrap']),
    ('só diag-wraps', ['diag-wrap']),
    ('x+y wraps', ['x-wrap', 'y-wrap']),
    ('x+diag wraps', ['x-wrap', 'diag-wrap']),
    ('y+diag wraps', ['y-wrap', 'diag-wrap']),
    ('todas wraps (plano)', ['x-wrap', 'y-wrap', 'diag-wrap']),
]

for name, cat_list in configs:
    rem = []
    for c in cat_list:
        rem.extend(cat_idx(c))
    Q, d2, _ = compute_Q(make_mask(rem))
    results.append((name, len(rem), len(d2), Q, ''))
    deg2_str = ','.join(f"({r},{c})" for r,c in d2)
    print(f"  {name}: removed={len(rem)}, deg2={len(d2)} [{deg2_str}], Q={Q}")


# ============================================================
# 4f. FINE SWEEP — all subsets of wrap edges
# ============================================================

print("\n" + "=" * 70)
print("FINE SWEEP — All subsets of removable wrap edges")
print("=" * 70)

n_wraps = len(wrap_indices)
total_sub = 2 ** n_wraps
print(f"Removable wrap edges: {n_wraps}")
print(f"Total subsets: {total_sub}")

SAMPLE = total_sub > 100000

sweep = {}  # k → [(subset, Q, deg2_count, deg2_coords)]
size_Q = {}  # k → set of Q values

for k in range(n_wraps + 1):
    if SAMPLE and 0 < k < n_wraps:
        from math import comb
        tot = comb(n_wraps, k)
        if tot > 500:
            s = set()
            while len(s) < 500:
                s.add(tuple(sorted(random.sample(range(n_wraps), k))))
            subs = list(s)
        else:
            subs = list(combinations(range(n_wraps), k))
    else:
        subs = list(combinations(range(n_wraps), k))
    
    sweep[k] = []
    if k not in size_Q:
        size_Q[k] = set()
    
    for sub in subs:
        rem = [wrap_indices[i] for i in sub]
        Q, d2, _ = compute_Q(make_mask(rem))
        sweep[k].append((sub, Q, len(d2), d2))
        size_Q[k].add(Q)
    
    Qs = [x[1] for x in sweep[k]]
    qd = {}
    for q in Qs:
        qd[q] = qd.get(q, 0) + 1
    print(f"  k={k:2d}: {len(subs):6d} subsets | Q: {dict(sorted(qd.items()))}")


# ============================================================
# 5. FORMATTED OUTPUT
# ============================================================

# Summary table
print("\n" + "=" * 70)
print("SUMMARY TABLE")
print("=" * 70)
print(f"| {'config':<22} | {'removidas':>9} | {'deg2_count':>10} | {'Q':>3} | {'notas':<30} |")
print(f"|{'-'*24}|{'-'*11}|{'-'*12}|{'-'*5}|{'-'*32}|")
for name, rem, d2, Q, notes in results:
    print(f"| {name:<22} | {rem:>9d} | {d2:>10d} | {Q:>3d} | {notes:<30} |")

# Q thresholds
print(f"\n--- First Q threshold achievements ---")
for thr in [1, 2, 3]:
    best = None
    for k in sorted(sweep.keys()):
        for sub, Q, d2c, d2v in sweep[k]:
            if Q >= thr and (best is None or k < best[0]):
                best = (k, sub, Q, d2c, d2v)
        if best and best[0] == k:
            break
    if best:
        k, sub, Q, d2c, d2v = best
        es = [f"{all_edges[wrap_indices[i]][0]}→{all_edges[wrap_indices[i]][1]}" for i in sub]
        print(f"  Q≥{thr}: first at k={k} removals (Q={Q})")
        print(f"    edges removed: {es}")
        print(f"    {d2c} deg-2 vertices: {d2v}")
    else:
        print(f"  Q≥{thr}: never reached")

# Q=3 configurations
print(f"\n--- Q=3 configurations ---")
q3 = [(k, sub, d2c, d2v) for k in sweep for sub, Q, d2c, d2v in sweep[k] if Q == 3]
print(f"Total Q=3 configs: {len(q3)}")
for k, sub, d2c, d2v in q3:
    es = [f"{all_edges[wrap_indices[i]][0]}→{all_edges[wrap_indices[i]][1]}" for i in sub]
    print(f"  k={k}: {es}")
    print(f"    deg-2 verts ({d2c}): {d2v}")

# Small-k interesting configs
print(f"\n--- Interesting configs (Q>0, k≤5) ---")
for k in sorted(sweep.keys()):
    for sub, Q, d2c, d2v in sweep[k]:
        if Q > 0 and k <= 5:
            es = [f"{all_edges[wrap_indices[i]][0]}→{all_edges[wrap_indices[i]][1]}" for i in sub]
            print(f"  k={k} Q={Q}: removed {es}, deg2={d2c} {d2v}")


# ============================================================
# 6. HYPOTHESIS TESTING
# ============================================================

print("\n" + "=" * 70)
print("HYPOTHESIS TESTING")
print("=" * 70)

# H1: Q monotone with |removals|?
print("\nH1: Q cresce monotonamente com o número de remoções?")
h1_weak = True
h1_detail = []
for k in sorted(sweep.keys()):
    if k + 1 in sweep:
        maxQ_k = max(x[1] for x in sweep[k])
        minQ_k1 = min(x[1] for x in sweep[k+1])
        if maxQ_k > minQ_k1:
            h1_weak = False
            h1_detail.append((k, maxQ_k, k+1, minQ_k1))

if h1_weak:
    print("  VERDADEIRO (cardinalidade) — Q nunca decresce com mais remoções")
else:
    print("  FALSO (cardinalidade) — Q pode diminuir com mais remoções!")
    for k1, mq1, k2, mq2 in h1_detail:
        print(f"    max_Q(k={k1})={mq1} > min_Q(k={k2})={mq2}")

# Check subset-monotonicity (stronger)
h1_sub = True
h1_sub_ex = None
for k in sorted(sweep.keys()):
    for s1, Q1, _, _ in sweep[k]:
        if Q1 > 0:
            for k2 in range(k+1, min(k+4, n_wraps+1)):
                if k2 not in sweep:
                    continue
                for s2, Q2, _, _ in sweep[k2]:
                    if Q2 < Q1 and set(s1).issubset(set(s2)):
                        h1_sub = False
                        h1_sub_ex = (k, s1, Q1, k2, s2, Q2)
                        break
                if not h1_sub:
                    break
        if not h1_sub:
            break
    if not h1_sub:
        break

if h1_sub:
    print("  VERDADEIRO (subconjunto) — para S⊂T, Q(S)≤Q(T) sempre")
else:
    k1, s1, q1, k2, s2, q2 = h1_sub_ex
    e1 = [f"{all_edges[wrap_indices[i]][0]}→{all_edges[wrap_indices[i]][1]}" for i in s1]
    e2 = [f"{all_edges[wrap_indices[i]][0]}→{all_edges[wrap_indices[i]][1]}" for i in s2]
    print(f"  FALSO (subconjunto) — contraexemplo encontrado:")
    print(f"    S={s1} (k={k1}, Q={q1}): {e1}")
    print(f"    T={s2} (k={k2}, Q={q2}): {e2}")

# H2: Q=3 requires ≥4 degree-2 vertices
print("\nH2: Q=3 requer pelo menos 4 vértices de grau 2?")
h2_ok = True
min_d2_q3 = None
for k in sweep:
    for sub, Q, d2c, d2v in sweep[k]:
        if Q == 3:
            if min_d2_q3 is None or d2c < min_d2_q3:
                min_d2_q3 = d2c
            if d2c < 4:
                h2_ok = False

if min_d2_q3 is not None:
    print(f"  {'VERDADEIRO' if h2_ok else 'FALSO'} — min deg2 quando Q=3: {min_d2_q3}")
else:
    print(f"  VACUAMENTE VERDADEIRO — Q=3 {'atingido' if q3 else 'nunca atingido'}")

# H3: Q(só x-wraps) ≠ Q(só y-wraps)?
print("\nH3: Remover só x-wraps dá Q ≠ remover só y-wraps?")
Qx = [r for r in results if r[0] == 'só x-wraps'][0]
Qy = [r for r in results if r[0] == 'só y-wraps'][0]
if Qx[3] != Qy[3]:
    print(f"  VERDADEIRO — Q(x-wraps, {Qx[1]} rem)={Qx[3]} ≠ Q(y-wraps, {Qy[1]} rem)={Qy[3]}")
else:
    print(f"  FALSO — Q(x-wraps, {Qx[1]} rem)={Qx[3]} = Q(y-wraps, {Qy[1]} rem)={Qy[3]}")
    if Qx[1] == 0 and Qy[1] == 0:
        print(f"  NOTA: nenhuma aresta é removível nessas categorias!")
        print(f"  (Todas x-wraps e y-wraps do toro TAMBÉM existem no plano.)")
        print(f"  H3 é trivialmente falsa — a pergunta original pressupõe wrap-edges")
        print(f"  separáveis por eixo, mas no 4×4 todas as wrap-only edges são diag-wraps.")

# H4: same-size subsets with different Q?
print("\nH4: Existem subsets do mesmo tamanho com Q diferentes?")
h4 = False
for k in sorted(size_Q.keys()):
    if len(size_Q[k]) > 1:
        h4 = True
        exs = {}
        for sub, Q, d2c, d2v in sweep[k]:
            if Q not in exs:
                es = [f"{all_edges[wrap_indices[i]][0]}→{all_edges[wrap_indices[i]][1]}" for i in sub]
                exs[Q] = (sub, es, d2c, d2v)
            if len(exs) >= 2:
                break
        print(f"  VERDADEIRO — primeiro em k={k}: Q values = {sorted(size_Q[k])}")
        for qv, (sub, es, d2c, d2v) in sorted(exs.items()):
            print(f"    Q={qv}: subset={sub}, edges={es}, deg2={d2c}")
        print(f"  → Geometria (quais wraps) importa, não só cardinalidade!")
        break

if not h4:
    print(f"  FALSO — todos subsets do mesmo tamanho dão o mesmo Q")


# ============================================================
# Q vs removal count
# ============================================================
print("\n" + "=" * 70)
print("Q vs. REMOVAL COUNT — Summary")
print("=" * 70)
for k in sorted(sweep.keys()):
    Qs = [x[1] for x in sweep[k]]
    print(f"  k={k:2d}: Q ∈ [{min(Qs):d}, {max(Qs):d}]  "
          f"mean={np.mean(Qs):.3f}  unique={sorted(set(Qs))}")

print("\n✓ Experiment complete.")
