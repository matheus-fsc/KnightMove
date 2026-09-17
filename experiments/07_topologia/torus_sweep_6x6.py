#!/usr/bin/env python3
"""
Torus 6×6 Knight Graph — Progressive Wrap Removal Sweep
========================================================
Repeats the 4×4 experiment on 6×6 where ±2 ≢ ∓2 (mod 6),
so x-only and y-only wrap edges should exist as distinct categories.
"""

import numpy as np
from itertools import combinations
from math import comb
import random
import time

N = 6
KNIGHT_MOVES = [(1,2),(1,-2),(-1,2),(-1,-2),(2,1),(2,-1),(-2,1),(-2,-1)]

def vid(r, c):
    return r * N + c

t0_global = time.time()

# ============================================================
# 1. BUILD TORUS AND PLANE GRAPHS, CLASSIFY EDGES
# ============================================================

t0 = time.time()

# Build torus edges: (v1,v2) → list of (wraps_r, wraps_c, dr, dc, r1, c1, r2, c2)
torus_edge_set = {}
for r in range(N):
    for c in range(N):
        for dr, dc in KNIGHT_MOVES:
            r2_raw, c2_raw = r + dr, c + dc
            r2 = r2_raw % N
            c2 = c2_raw % N
            wraps_r = (r2_raw != r2)
            wraps_c = (c2_raw != c2)
            v1, v2 = vid(r, c), vid(r2, c2)
            if v1 == v2:
                continue
            key = (min(v1,v2), max(v1,v2))
            if key not in torus_edge_set:
                torus_edge_set[key] = []
            torus_edge_set[key].append((wraps_r, wraps_c, dr, dc, r, c, r2, c2))

# Build plane edges
plane_edge_set = set()
for r in range(N):
    for c in range(N):
        for dr, dc in KNIGHT_MOVES:
            r2, c2 = r + dr, c + dc
            if 0 <= r2 < N and 0 <= c2 < N:
                v1, v2 = vid(r, c), vid(r2, c2)
                plane_edge_set.add((min(v1,v2), max(v1,v2)))

# Removable edges = torus-only (not in plane)
edges_to_remove = set(torus_edge_set.keys()) - plane_edge_set

# Classify each edge
all_edges = []  # ((r1,c1),(r2,c2), v1, v2, category)
for (v1, v2), moves in sorted(torus_edge_set.items()):
    r1, c1, r2, c2 = v1//N, v1%N, v2//N, v2%N
    if (v1, v2) in edges_to_remove:
        # Classify wrap sub-type from the moves that realize this edge
        has_r_wrap = any(wr for wr, wc, *_ in moves)
        has_c_wrap = any(wc for wr, wc, *_ in moves)
        # Check if ALL realizations wrap — and which axes
        # We want: is there a move realization where only r wraps? only c? both?
        wrap_types = set()
        for wr, wc, *_ in moves:
            if wr and wc:
                wrap_types.add('diag')
            elif wr and not wc:
                wrap_types.add('y-only')
            elif wc and not wr:
                wrap_types.add('x-only')
            else:
                wrap_types.add('internal')  # shouldn't happen for removable
        
        # Classify by the "simplest" wrap type present
        if 'internal' in wrap_types:
            cat = 'internal'  # has a non-wrapping realization — NOT removable??
            # This shouldn't happen since we filtered by edges_to_remove
            print(f"WARNING: edge ({r1},{c1})-({r2},{c2}) is removable but has internal realization!")
        elif 'x-only' in wrap_types:
            cat = 'x-only'
        elif 'y-only' in wrap_types:
            cat = 'y-only'
        elif 'diag' in wrap_types:
            cat = 'diag'
        else:
            cat = 'unknown'
        
        all_edges.append(((r1,c1),(r2,c2), v1, v2, cat))
    else:
        all_edges.append(((r1,c1),(r2,c2), v1, v2, 'internal'))

# Count by category
cats = {}
for e in all_edges:
    cats[e[4]] = cats.get(e[4], 0) + 1

print("=" * 70)
print("TORUS 6×6 KNIGHT GRAPH — EDGE CLASSIFICATION")
print("=" * 70)
print(f"Total torus edges: {len(all_edges)}  |  Plane edges: {len(plane_edge_set)}")
print(f"Removable wrap edges: {len(edges_to_remove)}")
for cat in ['internal', 'x-only', 'y-only', 'diag', 'unknown']:
    if cat in cats:
        print(f"  {cat:12s}: {cats[cat]}")

# Verify: do x-only and y-only exist?
has_x = cats.get('x-only', 0) > 0
has_y = cats.get('y-only', 0) > 0
print(f"\n✓ x-only wraps exist: {has_x} ({cats.get('x-only',0)} edges)")
print(f"✓ y-only wraps exist: {has_y} ({cats.get('y-only',0)} edges)")
print(f"✓ diag wraps exist: {cats.get('diag',0) > 0} ({cats.get('diag',0)} edges)")
print(f"  Confirmation: ±2 mod 6 = {{2,4}} ≠ {{-2,-4}} — x/y wraps ARE separable")

# Build wrap index arrays
wrap_indices = []   # index into all_edges
wrap_cats = []
x_only_idx = []     # indices into wrap_indices
y_only_idx = []
diag_idx = []

for i, e in enumerate(all_edges):
    if e[4] != 'internal':
        wi = len(wrap_indices)
        wrap_indices.append(i)
        wrap_cats.append(e[4])
        if e[4] == 'x-only':
            x_only_idx.append(wi)
        elif e[4] == 'y-only':
            y_only_idx.append(wi)
        elif e[4] == 'diag':
            diag_idx.append(wi)

print(f"\nWrap edge breakdown (wrap_indices):")
print(f"  x-only: {len(x_only_idx)} edges")
print(f"  y-only: {len(y_only_idx)} edges")
print(f"  diag:   {len(diag_idx)} edges")
print(f"  total:  {len(wrap_indices)} edges")

# Print each wrap edge
print(f"\nAll removable wrap edges:")
for wi, ei in enumerate(wrap_indices):
    e = all_edges[ei]
    print(f"  wrap[{wi:2d}] edge#{ei:3d}  {e[0]}→{e[1]}  ({e[4]})")

t1 = time.time()
print(f"\n⏱ Graph construction: {t1-t0:.2f}s")


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
# 3. COMPUTE Q
# ============================================================

def compute_Q(active_mask):
    """
    Q = rank([∂1; XOR-pairs]) - rank(∂1)
    XOR-pairs = {δ_{e_i} + δ_{e_j}} for all pairs of degree-2 vertices.
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
    
    # Incidence matrix
    inc = np.zeros((V, E), dtype=np.int64)
    for new, old in enumerate(active_idx):
        e = all_edges[old]
        inc[e[2], new] = 1
        inc[e[3], new] = 1
    
    rank_inc = gf2_rank(inc)
    
    if len(deg2) < 2:
        return 0, [(v//N, v%N) for v in deg2], {v: int(degree[v]) for v in range(V) if degree[v] > 0}
    
    # One forced edge per degree-2 vertex
    forced = [adj[v][0] for v in deg2]
    
    # XOR pairs
    xvecs = []
    for i in range(len(forced)):
        for j in range(i+1, len(forced)):
            vec = np.zeros(E, dtype=np.int64)
            vec[amap[forced[i]]] = 1
            vec[amap[forced[j]]] = 1
            xvecs.append(vec)
    
    if len(xvecs) == 0:
        return 0, [(v//N, v%N) for v in deg2], {v: int(degree[v]) for v in range(V) if degree[v] > 0}
    
    combined = np.vstack([inc, np.array(xvecs, dtype=np.int64)])
    Q = gf2_rank(combined) - rank_inc
    
    return Q, [(v//N, v%N) for v in deg2], {v: int(degree[v]) for v in range(V) if degree[v] > 0}


# ============================================================
# Helper
# ============================================================

def make_mask(removed_edge_indices):
    """removed_edge_indices are indices into all_edges."""
    mask = np.ones(len(all_edges), dtype=bool)
    for idx in removed_edge_indices:
        mask[idx] = False
    return mask


# ============================================================
# 3a. ANCHOR CONFIGURATIONS
# ============================================================

print("\n" + "=" * 70)
print("3a. ANCHOR CONFIGURATIONS")
print("=" * 70)

t0 = time.time()

anchor_results = []

def run_anchor(name, wrap_idx_list):
    """Remove the wrap edges at these wrap_indices positions."""
    removed = [wrap_indices[i] for i in wrap_idx_list]
    # Count by category
    cx = sum(1 for i in wrap_idx_list if wrap_cats[i] == 'x-only')
    cy = sum(1 for i in wrap_idx_list if wrap_cats[i] == 'y-only')
    cd = sum(1 for i in wrap_idx_list if wrap_cats[i] == 'diag')
    
    Q, deg2, degs = compute_Q(make_mask(removed))
    anchor_results.append((name, len(wrap_idx_list), cx, cy, cd, len(deg2), Q, deg2))
    deg2_str = ', '.join(f"({r},{c})" for r,c in deg2) if len(deg2) <= 8 else f"{len(deg2)} vertices"
    print(f"  {name:<28s}: removed={len(wrap_idx_list):2d} (x={cx},y={cy},d={cd}), "
          f"deg2={len(deg2):2d}, Q={Q}  [{deg2_str}]")
    return Q

# Toro puro
run_anchor("toro puro", [])

# Remove by category
run_anchor("só x-only wraps", x_only_idx)
run_anchor("só y-only wraps", y_only_idx)
run_anchor("só diag wraps", diag_idx)
run_anchor("x-only + y-only", x_only_idx + y_only_idx)
run_anchor("x-only + diag", x_only_idx + diag_idx)
run_anchor("y-only + diag", y_only_idx + diag_idx)

Q_plane = run_anchor("TODAS wraps (plano)", list(range(len(wrap_indices))))

if Q_plane == 3:
    print(f"  ✓ Confirmado Q=3 para o plano 6×6")
else:
    print(f"  ✗ INESPERADO: Q={Q_plane} para o plano 6×6")

t1 = time.time()
print(f"\n⏱ Anchor configs: {t1-t0:.2f}s")


# ============================================================
# 3b. FINE SWEEP by cardinality
# ============================================================

print("\n" + "=" * 70)
print("3b. FINE SWEEP BY CARDINALITY")
print("=" * 70)

t0 = time.time()

n_wraps = len(wrap_indices)
print(f"Total removable wrap edges: {n_wraps}")
print(f"Total subsets: 2^{n_wraps} = {2**n_wraps}")

SAMPLE_SIZE = 200

sweep = {}  # k → [(subset_of_wrap_indices, Q, deg2_count)]
sweep_summary = {}  # k → (min_Q, max_Q, count_Q_gt0, count_Q_eq3, total, Q_dist)

for k in range(n_wraps + 1):
    total_k = comb(n_wraps, k)
    
    if total_k <= SAMPLE_SIZE:
        subs = list(combinations(range(n_wraps), k))
        sampled = False
    else:
        s = set()
        attempts = 0
        while len(s) < SAMPLE_SIZE and attempts < SAMPLE_SIZE * 20:
            s.add(tuple(sorted(random.sample(range(n_wraps), k))))
            attempts += 1
        subs = list(s)
        sampled = True
    
    sweep[k] = []
    q_dist = {}
    
    for sub in subs:
        removed = [wrap_indices[i] for i in sub]
        Q, deg2, _ = compute_Q(make_mask(removed))
        sweep[k].append((sub, Q, len(deg2)))
        q_dist[Q] = q_dist.get(Q, 0) + 1
    
    Qs = [x[1] for x in sweep[k]]
    minQ, maxQ = min(Qs), max(Qs)
    n_gt0 = sum(1 for q in Qs if q > 0)
    n_eq3 = sum(1 for q in Qs if q == 3)
    n_total = len(Qs)
    sweep_summary[k] = (minQ, maxQ, n_gt0, n_eq3, n_total, q_dist)
    
    samp_str = f" (sampled {n_total}/{total_k})" if sampled else f" (exhaustive {n_total})"
    print(f"  k={k:2d}: Q∈[{minQ},{maxQ}]  Q>0: {n_gt0}/{n_total} ({100*n_gt0/n_total:.0f}%)  "
          f"Q=3: {n_eq3}/{n_total}  dist={dict(sorted(q_dist.items()))}{samp_str}")

t1 = time.time()
print(f"\n⏱ Fine sweep: {t1-t0:.2f}s")

# Minimum k for Q=3
min_k_q3 = None
for k in sorted(sweep.keys()):
    for sub, Q, d2c in sweep[k]:
        if Q == 3:
            if min_k_q3 is None or k < min_k_q3[0]:
                cats_sub = [wrap_cats[i] for i in sub]
                min_k_q3 = (k, sub, d2c, cats_sub)
            break
    if min_k_q3 and min_k_q3[0] == k:
        pass  # keep looking at same k for potentially better example

print(f"\n--- Minimum k for Q=3 ---")
if min_k_q3:
    k, sub, d2c, cs = min_k_q3
    cx = cs.count('x-only')
    cy = cs.count('y-only')
    cd = cs.count('diag')
    edges = [f"{all_edges[wrap_indices[i]][0]}→{all_edges[wrap_indices[i]][1]}" for i in sub]
    print(f"  k={k}: x-only={cx}, y-only={cy}, diag={cd}, deg2={d2c}")
    print(f"  edges: {edges[:12]}{'...' if len(edges)>12 else ''}")
    _, d2_coords, _ = compute_Q(make_mask([wrap_indices[i] for i in sub]))
    print(f"  deg-2 vertices: {d2_coords}")
else:
    print(f"  Q=3 never reached in sweep")


# ============================================================
# 3c. TEST H3: x-only vs y-only subsets of same size
# ============================================================

print("\n" + "=" * 70)
print("3c. H3 TEST: x-only vs y-only wraps")
print("=" * 70)

t0 = time.time()

# Test at various sizes
max_test_size = min(len(x_only_idx), len(y_only_idx))
print(f"x-only count: {len(x_only_idx)}, y-only count: {len(y_only_idx)}")
print(f"Testing sizes k=1..{max_test_size}")

h3_results = []
h3_found = False

for k in range(1, max_test_size + 1):
    # Sample x-only subsets
    total_x = comb(len(x_only_idx), k)
    total_y = comb(len(y_only_idx), k)
    
    n_sample = min(100, total_x)
    if total_x <= 100:
        x_subs = list(combinations(x_only_idx, k))
    else:
        s = set()
        while len(s) < n_sample:
            s.add(tuple(sorted(random.sample(x_only_idx, k))))
        x_subs = list(s)
    
    n_sample_y = min(100, total_y)
    if total_y <= 100:
        y_subs = list(combinations(y_only_idx, k))
    else:
        s = set()
        while len(s) < n_sample_y:
            s.add(tuple(sorted(random.sample(y_only_idx, k))))
        y_subs = list(s)
    
    x_Qs = []
    for sub in x_subs:
        removed = [wrap_indices[i] for i in sub]
        Q, _, _ = compute_Q(make_mask(removed))
        x_Qs.append(Q)
    
    y_Qs = []
    for sub in y_subs:
        removed = [wrap_indices[i] for i in sub]
        Q, _, _ = compute_Q(make_mask(removed))
        y_Qs.append(Q)
    
    x_set = sorted(set(x_Qs))
    y_set = sorted(set(y_Qs))
    
    print(f"  k={k}: x-only Q={x_set} (mean={np.mean(x_Qs):.2f})  "
          f"y-only Q={y_set} (mean={np.mean(y_Qs):.2f})  "
          f"{'DIFFERENT!' if x_set != y_set else 'same'}")
    
    h3_results.append((k, x_set, y_set, np.mean(x_Qs), np.mean(y_Qs)))
    
    if x_set != y_set and not h3_found:
        h3_found = True
        # Find concrete examples
        x_ex = None
        y_ex = None
        for sub, Q in zip(x_subs, x_Qs):
            if Q not in [q for q in y_Qs]:
                x_ex = (sub, Q)
                break
        for sub, Q in zip(y_subs, y_Qs):
            if Q not in [q for q in x_Qs]:
                y_ex = (sub, Q)
                break
        if x_ex:
            print(f"    x-only example: subset={x_ex[0]}, Q={x_ex[1]}")
        if y_ex:
            print(f"    y-only example: subset={y_ex[0]}, Q={y_ex[1]}")

# Also test: remove ALL x-only vs ALL y-only
print(f"\n  Full removal comparison:")
Q_all_x, d2_x, _ = compute_Q(make_mask([wrap_indices[i] for i in x_only_idx]))
Q_all_y, d2_y, _ = compute_Q(make_mask([wrap_indices[i] for i in y_only_idx]))
print(f"    Remove ALL x-only ({len(x_only_idx)}): Q={Q_all_x}, deg2={len(d2_x)}")
print(f"    Remove ALL y-only ({len(y_only_idx)}): Q={Q_all_y}, deg2={len(d2_y)}")
print(f"    H3 (full removal): {'VERDADEIRO — Q differs!' if Q_all_x != Q_all_y else 'FALSO — same Q'}")

t1 = time.time()
print(f"\n⏱ H3 test: {t1-t0:.2f}s")


# ============================================================
# 3d. INCLUSION CHAIN MONOTONICITY TEST
# ============================================================

print("\n" + "=" * 70)
print("3d. INCLUSION CHAIN MONOTONICITY TEST")
print("=" * 70)

t0 = time.time()

# Build a chain S1 ⊂ S2 ⊂ ... ⊂ S_all
# Strategy: add wrap edges one at a time in various orders

def test_chain(order, label):
    """Test monotonicity along a chain defined by the order of wrap index additions."""
    chain_Qs = []
    for i in range(len(order) + 1):
        sub = order[:i]
        removed = [wrap_indices[j] for j in sub]
        Q, d2, _ = compute_Q(make_mask(removed))
        chain_Qs.append((i, Q, len(d2)))
    
    monotone = True
    violations = []
    for i in range(1, len(chain_Qs)):
        if chain_Qs[i][1] < chain_Qs[i-1][1]:
            monotone = False
            violations.append((i-1, chain_Qs[i-1], i, chain_Qs[i]))
    
    print(f"\n  Chain '{label}':")
    print(f"    Q sequence: {[q for _,q,_ in chain_Qs]}")
    print(f"    Monotone: {'YES ✓' if monotone else 'NO ✗'}")
    if violations:
        for i1, (k1,q1,d1), i2, (k2,q2,d2) in violations[:3]:
            print(f"      Violation at step {i1}→{i2}: Q={q1}→{q2}")
    return monotone, chain_Qs

# Chain 1: natural order (as enumerated)
natural_order = list(range(n_wraps))
mon1, _ = test_chain(natural_order, "natural order")

# Chain 2: x-only first, then y-only, then diag
typed_order = x_only_idx + y_only_idx + diag_idx
mon2, _ = test_chain(typed_order, "x-only → y-only → diag")

# Chain 3: diag first, then x, then y
typed_order2 = diag_idx + x_only_idx + y_only_idx
mon3, _ = test_chain(typed_order2, "diag → x-only → y-only")

# Chain 4-8: random orders
random_mons = []
for trial in range(5):
    rng_order = list(range(n_wraps))
    random.shuffle(rng_order)
    mon, _ = test_chain(rng_order, f"random #{trial+1}")
    random_mons.append(mon)

all_chains_monotone = mon1 and mon2 and mon3 and all(random_mons)
print(f"\n  Overall: {'ALL chains monotone ✓' if all_chains_monotone else 'Some chains NOT monotone ✗'}")

# Test subset monotonicity more carefully via sampling
print(f"\n  --- Subset monotonicity sampling ---")
h_mon_counter = None
n_tests = 0
for _ in range(2000):
    # Pick a random subset S and a random superset T
    k1 = random.randint(1, n_wraps - 1)
    S = tuple(sorted(random.sample(range(n_wraps), k1)))
    # Add 1-3 more elements to form T
    remaining = [i for i in range(n_wraps) if i not in S]
    if not remaining:
        continue
    extra = random.randint(1, min(3, len(remaining)))
    T_extra = random.sample(remaining, extra)
    T = tuple(sorted(list(S) + T_extra))
    
    Q_S, _, _ = compute_Q(make_mask([wrap_indices[i] for i in S]))
    Q_T, _, _ = compute_Q(make_mask([wrap_indices[i] for i in T]))
    n_tests += 1
    
    if Q_T < Q_S:
        h_mon_counter = (S, Q_S, T, Q_T)
        break

if h_mon_counter:
    S, qs, T, qt = h_mon_counter
    print(f"  FALSO — counterexample found in {n_tests} tests:")
    print(f"    S={S} (k={len(S)}, Q={qs})")
    print(f"    T={T} (k={len(T)}, Q={qt})")
else:
    print(f"  VERDADEIRO — no violations in {n_tests} subset pairs tested")

t1 = time.time()
print(f"\n⏱ Monotonicity tests: {t1-t0:.2f}s")


# ============================================================
# 3e. ANSWER CENTRAL QUESTION:
#     Can Q=3 be reached by removing only x-wraps or only y-wraps?
# ============================================================

print("\n" + "=" * 70)
print("CENTRAL QUESTION: Can Q=3 be reached with single-axis wraps only?")
print("=" * 70)

# Already computed above, but let's be explicit
print(f"  Remove ALL x-only wraps ({len(x_only_idx)}): Q = {Q_all_x}")
print(f"  Remove ALL y-only wraps ({len(y_only_idx)}): Q = {Q_all_y}")
print(f"  Remove ALL diag wraps   ({len(diag_idx)}):   Q = ", end="")
Q_all_d, d2_d, _ = compute_Q(make_mask([wrap_indices[i] for i in diag_idx]))
print(f"{Q_all_d}")

if Q_all_x == 3:
    print(f"  → YES: Q=3 reachable by removing only x-wraps!")
elif Q_all_y == 3:
    print(f"  → YES: Q=3 reachable by removing only y-wraps!")
else:
    print(f"  → NO: Q=3 requires removing wraps from multiple directions")
    print(f"         (x-only gives Q={Q_all_x}, y-only gives Q={Q_all_y})")

# Check minimum mixed requirement for Q=3
# Try: x-only + just enough of other types
print(f"\n  Searching for minimal Q=3 configurations by category mix...")
found_q3 = []
# Try adding diag to x-only
for k_extra in range(1, len(diag_idx) + len(y_only_idx) + 1):
    extras = diag_idx + y_only_idx
    if k_extra > len(extras):
        break
    for combo in (combinations(extras, k_extra) if comb(len(extras), k_extra) <= 200 
                  else [tuple(sorted(random.sample(extras, k_extra))) for _ in range(200)]):
        sub = list(set(x_only_idx + list(combo)))
        removed = [wrap_indices[i] for i in sub]
        Q, d2, _ = compute_Q(make_mask(removed))
        if Q == 3:
            cats_sub = [wrap_cats[i] for i in sub]
            cx = cats_sub.count('x-only')
            cy = cats_sub.count('y-only')
            cd = cats_sub.count('diag')
            found_q3.append((len(sub), cx, cy, cd, sub, d2))
    if found_q3:
        break

if found_q3:
    found_q3.sort()
    k, cx, cy, cd, sub, d2 = found_q3[0]
    print(f"  x-only + extras → Q=3 at k={k}: x={cx}, y={cy}, diag={cd}")

# Try: y-only + extras
found_q3_y = []
for k_extra in range(1, len(diag_idx) + len(x_only_idx) + 1):
    extras = diag_idx + x_only_idx
    if k_extra > len(extras):
        break
    for combo in (combinations(extras, k_extra) if comb(len(extras), k_extra) <= 200
                  else [tuple(sorted(random.sample(extras, k_extra))) for _ in range(200)]):
        sub = list(set(y_only_idx + list(combo)))
        removed = [wrap_indices[i] for i in sub]
        Q, d2, _ = compute_Q(make_mask(removed))
        if Q == 3:
            cats_sub = [wrap_cats[i] for i in sub]
            cx = cats_sub.count('x-only')
            cy = cats_sub.count('y-only')
            cd = cats_sub.count('diag')
            found_q3_y.append((len(sub), cx, cy, cd, sub, d2))
    if found_q3_y:
        break

if found_q3_y:
    found_q3_y.sort()
    k, cx, cy, cd, sub, d2 = found_q3_y[0]
    print(f"  y-only + extras → Q=3 at k={k}: x={cx}, y={cy}, diag={cd}")


# ============================================================
# 4. FINAL REPORT
# ============================================================

print("\n" + "=" * 70)
print("ANCHOR TABLE")
print("=" * 70)

print(f"| {'config':<28s} | {'|wraps|':>7s} | {'x-only':>6s} | {'y-only':>6s} | {'diag':>4s} | {'deg2':>4s} | {'Q':>3s} |")
print(f"|{'-'*30}|{'-'*9}|{'-'*8}|{'-'*8}|{'-'*6}|{'-'*6}|{'-'*5}|")
for name, nw, cx, cy, cd, d2c, Q, d2v in anchor_results:
    print(f"| {name:<28s} | {nw:>7d} | {cx:>6d} | {cy:>6d} | {cd:>4d} | {d2c:>4d} | {Q:>3d} |")

print(f"\nSWEEP TABLE")
print(f"| {'k':>3s} | {'min_Q':>5s} | {'max_Q':>5s} | {'frac_Q>0':>10s} | {'frac_Q=3':>10s} |")
print(f"|{'-'*5}|{'-'*7}|{'-'*7}|{'-'*12}|{'-'*12}|")
for k in sorted(sweep_summary.keys()):
    minQ, maxQ, ngt0, neq3, ntot, _ = sweep_summary[k]
    f1 = f"{ngt0}/{ntot}" if ntot <= 200 else f"{100*ngt0/ntot:.1f}%"
    f3 = f"{neq3}/{ntot}" if ntot <= 200 else f"{100*neq3/ntot:.1f}%"
    print(f"| {k:>3d} | {minQ:>5d} | {maxQ:>5d} | {f1:>10s} | {f3:>10s} |")


# ============================================================
# HYPOTHESIS SUMMARY
# ============================================================

print("\n" + "=" * 70)
print("HYPOTHESIS TESTING SUMMARY")
print("=" * 70)

# H1: Cardinality monotonicity
h1_card = True
h1_card_ex = None
for k in sorted(sweep_summary.keys()):
    if k + 1 in sweep_summary:
        maxQ_k = sweep_summary[k][1]
        minQ_k1 = sweep_summary[k+1][0]
        if maxQ_k > minQ_k1:
            h1_card = False
            h1_card_ex = (k, maxQ_k, k+1, minQ_k1)
            break

print(f"\nH1 (monotonia cardinalidade):")
if h1_card:
    print(f"  VERDADEIRO — max_Q(k) ≤ min_Q(k+1) para todo k")
else:
    k1, mq1, k2, mq2 = h1_card_ex
    print(f"  FALSO — max_Q(k={k1})={mq1} > min_Q(k={k2})={mq2}")

# H3
print(f"\nH3 (x-wraps vs y-wraps):")
if Q_all_x != Q_all_y:
    print(f"  VERDADEIRO — Q(all x-only)={Q_all_x} ≠ Q(all y-only)={Q_all_y}")
    print(f"  Simetria x↔y é QUEBRADA no 6×6!")
elif h3_found:
    print(f"  VERDADEIRO — encontrado em subsets de mesmo tamanho com Q diferentes")
else:
    print(f"  FALSO — Q(all x-only)={Q_all_x} = Q(all y-only)={Q_all_y}")
    # Check if per-k also matches
    any_diff = any(x_set != y_set for _, x_set, y_set, _, _ in h3_results)
    if any_diff:
        print(f"  MAS: existem tamanhos k onde a distribuição de Q difere!")
    else:
        print(f"  Simetria x↔y preservada em todos os tamanhos testados")

# H_MON
print(f"\nH_MON (monotonia inclusão):")
if h_mon_counter is None:
    print(f"  VERDADEIRO — nenhum contraexemplo em {n_tests} testes S⊂T")
    print(f"  + todas as {2+5} cadeias testadas são monotônicas: {all_chains_monotone}")
else:
    S, qs, T, qt = h_mon_counter
    print(f"  FALSO — S={S} (Q={qs}) ⊂ T={T} (Q={qt})")

# H_NEW: proportion comparison with 4×4
print(f"\nH_NEW (comparação proporcional com 4×4):")
if min_k_q3:
    k_6x6 = min_k_q3[0]
    prop_6x6 = k_6x6 / n_wraps
    # 4×4: Q=3 at k=8 out of 8 wraps = 100%
    prop_4x4 = 8 / 8
    print(f"  4×4: Q=3 em k=8/{8} = {prop_4x4:.0%} das wraps")
    print(f"  6×6: Q=3 em k={k_6x6}/{n_wraps} = {prop_6x6:.0%} das wraps")
    if prop_6x6 < prop_4x4:
        print(f"  → 6×6 atinge Q=3 proporcionalmente ANTES que o 4×4")
    elif prop_6x6 > prop_4x4:
        print(f"  → 6×6 requer proporcionalmente MAIS remoções que o 4×4")
    else:
        print(f"  → Mesma proporção")
else:
    print(f"  Q=3 não atingido no sweep — inconclusivo")

t_total = time.time() - t0_global
print(f"\n{'='*70}")
print(f"⏱ Tempo total de execução: {t_total:.1f}s")
print(f"✓ Experimento 6×6 completo.")
