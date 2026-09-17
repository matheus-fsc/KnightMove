#!/usr/bin/env python3
"""
Torus 6×6 — Surgical Analysis of Diag-Wrap Role
================================================
Central question: Do the 8 diag-wraps alone initiate the obstruction?
"""

import numpy as np
from itertools import combinations
import random
import time
import math

N = 6
KNIGHT_MOVES = [(1,2),(1,-2),(-1,2),(-1,-2),(2,1),(2,-1),(-2,1),(-2,-1)]

def vid(r, c):
    return r * N + c

# ============================================================
# GRAPH CONSTRUCTION (reused from previous experiment)
# ============================================================

torus_edge_set = {}
for r in range(N):
    for c in range(N):
        for dr, dc in KNIGHT_MOVES:
            r2_raw, c2_raw = r + dr, c + dc
            r2, c2 = r2_raw % N, c2_raw % N
            wraps_r = (r2_raw != r2)
            wraps_c = (c2_raw != c2)
            v1, v2 = vid(r, c), vid(r2, c2)
            if v1 == v2:
                continue
            key = (min(v1,v2), max(v1,v2))
            if key not in torus_edge_set:
                torus_edge_set[key] = []
            torus_edge_set[key].append((wraps_r, wraps_c, dr, dc, r, c, r2, c2))

plane_edge_set = set()
for r in range(N):
    for c in range(N):
        for dr, dc in KNIGHT_MOVES:
            r2, c2 = r + dr, c + dc
            if 0 <= r2 < N and 0 <= c2 < N:
                v1, v2 = vid(r, c), vid(r2, c2)
                plane_edge_set.add((min(v1,v2), max(v1,v2)))

edges_to_remove = set(torus_edge_set.keys()) - plane_edge_set

all_edges = []
for (v1, v2), moves in sorted(torus_edge_set.items()):
    r1, c1, r2, c2 = v1//N, v1%N, v2//N, v2%N
    if (v1, v2) in edges_to_remove:
        wrap_types = set()
        for wr, wc, *_ in moves:
            if wr and wc: wrap_types.add('diag')
            elif wr and not wc: wrap_types.add('y-only')
            elif wc and not wr: wrap_types.add('x-only')
        cat = 'x-only' if 'x-only' in wrap_types else ('y-only' if 'y-only' in wrap_types else 'diag')
    else:
        cat = 'internal'
    all_edges.append(((r1,c1),(r2,c2), v1, v2, cat))

# Build category indices
wrap_indices = []  # index into all_edges
wrap_cats = []
x_only_idx = []    # index into wrap_indices
y_only_idx = []
diag_idx = []

for i, e in enumerate(all_edges):
    if e[4] != 'internal':
        wi = len(wrap_indices)
        wrap_indices.append(i)
        wrap_cats.append(e[4])
        if e[4] == 'x-only': x_only_idx.append(wi)
        elif e[4] == 'y-only': y_only_idx.append(wi)
        elif e[4] == 'diag': diag_idx.append(wi)

print(f"Graph: {len(all_edges)} edges, {len(wrap_indices)} wraps "
      f"(x={len(x_only_idx)}, y={len(y_only_idx)}, d={len(diag_idx)})")

# Build x↔y mirror mapping
# Mirror: swap r↔c for each wrap edge
def mirror_wrap_index(wi):
    """Find the y-only wrap that mirrors x-only wrap wi (and vice versa)."""
    ei = wrap_indices[wi]
    e = all_edges[ei]
    (r1,c1), (r2,c2) = e[0], e[1]
    # Mirror: (r,c) → (c,r)
    mr1, mc1 = c1, r1
    mr2, mc2 = c2, r2
    mv1, mv2 = vid(mr1, mc1), vid(mr2, mc2)
    mkey = (min(mv1,mv2), max(mv1,mv2))
    # Find in wrap_indices
    for wj, ej in enumerate(wrap_indices):
        if (all_edges[ej][2], all_edges[ej][3]) == mkey:
            return wj
    return None

# Build mirror map
mirror_map = {}
for wi in x_only_idx:
    mwi = mirror_wrap_index(wi)
    if mwi is not None:
        mirror_map[wi] = mwi
        mirror_map[mwi] = wi


# ============================================================
# GF(2) RANK + COMPUTE Q
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

def compute_Q(active_mask):
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
    
    inc = np.zeros((V, E), dtype=np.int64)
    for new, old in enumerate(active_idx):
        e = all_edges[old]
        inc[e[2], new] = 1
        inc[e[3], new] = 1
    rank_inc = gf2_rank(inc)
    
    if len(deg2) < 2:
        return 0, [(v//N, v%N) for v in deg2], {v: int(degree[v]) for v in range(V) if degree[v] > 0}
    
    forced = [adj[v][0] for v in deg2]
    xvecs = []
    for i in range(len(forced)):
        for j in range(i+1, len(forced)):
            vec = np.zeros(E, dtype=np.int64)
            vec[amap[forced[i]]] = 1
            vec[amap[forced[j]]] = 1
            xvecs.append(vec)
    
    if not xvecs:
        return 0, [(v//N, v%N) for v in deg2], {v: int(degree[v]) for v in range(V) if degree[v] > 0}
    
    combined = np.vstack([inc, np.array(xvecs, dtype=np.int64)])
    Q = gf2_rank(combined) - rank_inc
    return Q, [(v//N, v%N) for v in deg2], {v: int(degree[v]) for v in range(V) if degree[v] > 0}

def make_mask(removed_wrap_idx_list):
    """removed_wrap_idx_list: indices into wrap_indices."""
    mask = np.ones(len(all_edges), dtype=bool)
    for wi in removed_wrap_idx_list:
        mask[wrap_indices[wi]] = False
    return mask


# ============================================================
# 1a. REMOVE ONLY THE 8 DIAG-WRAPS
# ============================================================

print("\n" + "=" * 70)
print("1a. REMOVE ONLY THE 8 DIAG-WRAPS")
print("=" * 70)

Q_diag, deg2_diag, degs_diag = compute_Q(make_mask(diag_idx))
print(f"Q = {Q_diag}")
print(f"deg2_count = {len(deg2_diag)}")
print(f"deg2 vertices: {deg2_diag}")
print(f"Degree distribution: {sorted(set(degs_diag.values()))}")

# Show degrees of the 4 plane corners
corners = [(0,0), (0,5), (5,0), (5,5)]
for r,c in corners:
    v = vid(r,c)
    d = degs_diag.get(v, 0)
    print(f"  Corner ({r},{c}): degree = {d}")


# ============================================================
# 1b. DIAG + 1 x-only (each of 28)
# ============================================================

print("\n" + "=" * 70)
print("1b. REMOVE 8 DIAG + 1 x-only (each of 28 x-only edges)")
print("=" * 70)

results_1b = []
for xi in x_only_idx:
    removed = diag_idx + [xi]
    Q, deg2, degs = compute_Q(make_mask(removed))
    ei = wrap_indices[xi]
    e = all_edges[ei]
    results_1b.append((xi, Q, len(deg2), deg2, e[0], e[1]))

q_dist_1b = {}
for _, Q, *_ in results_1b:
    q_dist_1b[Q] = q_dist_1b.get(Q, 0) + 1
print(f"Q distribution: {dict(sorted(q_dist_1b.items()))}")

for xi, Q, d2c, d2v, src, dst in results_1b:
    if Q > 0:
        print(f"  Q={Q}: +x-only {src}→{dst}, deg2={d2c} {d2v}")


# ============================================================
# 1c. DIAG + 1 y-only (each of 28)
# ============================================================

print("\n" + "=" * 70)
print("1c. REMOVE 8 DIAG + 1 y-only (each of 28 y-only edges)")
print("=" * 70)

results_1c = []
for yi in y_only_idx:
    removed = diag_idx + [yi]
    Q, deg2, degs = compute_Q(make_mask(removed))
    ei = wrap_indices[yi]
    e = all_edges[ei]
    results_1c.append((yi, Q, len(deg2), deg2, e[0], e[1]))

q_dist_1c = {}
for _, Q, *_ in results_1c:
    q_dist_1c[Q] = q_dist_1c.get(Q, 0) + 1
print(f"Q distribution: {dict(sorted(q_dist_1c.items()))}")

for yi, Q, d2c, d2v, src, dst in results_1c:
    if Q > 0:
        print(f"  Q={Q}: +y-only {src}→{dst}, deg2={d2c} {d2v}")


# ============================================================
# 1d. DIAG + k x-only + k y-only (balanced)
# ============================================================

print("\n" + "=" * 70)
print("1d. REMOVE 8 DIAG + k x-only + k y-only (balanced, k=1..14)")
print("=" * 70)

for k in range(1, 15):
    n_sample = min(200, math.comb(len(x_only_idx), k) * math.comb(len(y_only_idx), k))
    
    Qs = []
    d2s = []
    d2_verts_seen = set()
    
    for trial in range(min(200, n_sample)):
        sx = sorted(random.sample(x_only_idx, k))
        sy = sorted(random.sample(y_only_idx, k))
        removed = diag_idx + sx + sy
        Q, deg2, _ = compute_Q(make_mask(removed))
        Qs.append(Q)
        d2s.append(len(deg2))
        for v in deg2:
            d2_verts_seen.add(v)
    
    q_dist = {}
    for q in Qs:
        q_dist[q] = q_dist.get(q, 0) + 1
    
    print(f"  k={k:2d}: total_removed={8+2*k:2d}  Q∈[{min(Qs)},{max(Qs)}] "
          f"mean={np.mean(Qs):.2f}  dist={dict(sorted(q_dist.items()))}  "
          f"deg2∈[{min(d2s)},{max(d2s)}]  unique_deg2_verts={sorted(d2_verts_seen)[:8]}{'...' if len(d2_verts_seen)>8 else ''}")


# ============================================================
# 1e. GREEDY SEARCH: diag fixed, add x/y one at a time to maximize Q
# ============================================================

print("\n" + "=" * 70)
print("1e. GREEDY SEARCH: 8 diag fixed + add x/y to maximize Q")
print("=" * 70)

current_set = list(diag_idx)  # start with 8 diag wraps
remaining = list(x_only_idx) + list(y_only_idx)
greedy_log = []

Q0, d2_0, _ = compute_Q(make_mask(current_set))
greedy_log.append((0, 'start', Q0, len(d2_0), d2_0, len(current_set)))
print(f"  Step 0: Q={Q0}, deg2={len(d2_0)}, total_removed={len(current_set)}")

step = 0
while remaining and Q0 < 3:
    best_wi = None
    best_Q = -1
    best_d2 = None
    
    for wi in remaining:
        trial = current_set + [wi]
        Q, d2, _ = compute_Q(make_mask(trial))
        if Q > best_Q or (Q == best_Q and best_d2 is not None and len(d2) > len(best_d2)):
            best_Q = Q
            best_wi = wi
            best_d2 = d2
    
    current_set.append(best_wi)
    remaining.remove(best_wi)
    Q0 = best_Q
    step += 1
    
    ei = wrap_indices[best_wi]
    e = all_edges[ei]
    cat = e[4]
    greedy_log.append((step, f"+{cat} {e[0]}→{e[1]}", best_Q, len(best_d2), best_d2, len(current_set)))
    print(f"  Step {step}: +{cat} {e[0]}→{e[1]}  → Q={best_Q}, deg2={len(best_d2)} {best_d2}, total_removed={len(current_set)}")

print(f"\n  GREEDY RESULT: Q=3 reached after {step} additions beyond 8 diag-wraps")
print(f"  Total wraps removed: {len(current_set)} (8 diag + {step} x/y)")

# Category breakdown of greedy solution
cats_greedy = [wrap_cats[wi] for wi in current_set]
cx = cats_greedy.count('x-only')
cy = cats_greedy.count('y-only')
cd = cats_greedy.count('diag')
print(f"  Mix: x-only={cx}, y-only={cy}, diag={cd}")

# Also try greedy in reverse: start from plane (all removed), add back wraps to minimize Q
print("\n--- Reverse greedy: start from plane, add back wraps ---")
current_removed = set(range(len(wrap_indices)))  # all wraps removed = plane
remaining_add = list(range(len(wrap_indices)))

Q_plane, d2_plane, _ = compute_Q(make_mask(list(current_removed)))
print(f"  Start (plane): Q={Q_plane}, deg2={len(d2_plane)}")

reverse_log = []
while Q_plane > 0 and remaining_add:
    best_wi_add = None
    best_Q_add = 999
    best_d2_add = None
    
    for wi in remaining_add:
        trial = current_removed - {wi}  # add back = un-remove
        Q, d2, _ = compute_Q(make_mask(list(trial)))
        if Q < best_Q_add or (Q == best_Q_add and best_d2_add is not None and len(d2) < len(best_d2_add)):
            best_Q_add = Q
            best_wi_add = wi
            best_d2_add = d2
    
    current_removed.discard(best_wi_add)
    remaining_add.remove(best_wi_add)
    Q_plane = best_Q_add
    
    ei = wrap_indices[best_wi_add]
    e = all_edges[ei]
    cat = e[4]
    reverse_log.append((f"restore {cat} {e[0]}→{e[1]}", best_Q_add, len(best_d2_add), best_d2_add))
    print(f"  Restore {cat} {e[0]}→{e[1]}  → Q={best_Q_add}, deg2={len(best_d2_add)}, "
          f"remaining_removed={len(current_removed)}")

print(f"\n  REVERSE GREEDY: Q dropped to 0 after restoring {64 - len(current_removed)} wraps")
print(f"  Minimum wraps to keep removed for Q>0: >{len(current_removed)}")


# ============================================================
# 2. STRUCTURE OF DEGREE-2 VERTICES
# ============================================================

print("\n" + "=" * 70)
print("2. STRUCTURE OF DEGREE-2 VERTICES")
print("=" * 70)

# Collect deg2 vertices from all configs that have deg2 > 0
print("\n--- Survey: which vertices can become degree-2? ---")

# Sample many configs around the transition zone (k=50..64)
all_d2_verts = {}  # (r,c) → count of appearances

for trial in range(500):
    k_total = random.randint(45, 63)
    sub = random.sample(range(len(wrap_indices)), k_total)
    Q, d2, _ = compute_Q(make_mask(sub))
    if len(d2) > 0:
        for v in d2:
            all_d2_verts[v] = all_d2_verts.get(v, 0) + 1

print(f"Unique (r,c) positions that appeared as degree-2 (500 random configs, k=45..63):")
for v, count in sorted(all_d2_verts.items(), key=lambda x: -x[1]):
    print(f"  ({v[0]},{v[1]}): appeared {count} times")

# Check if they're always the same set
print(f"\nCorners: (0,0),(0,5),(5,0),(5,5)")
corner_set = {(0,0),(0,5),(5,0),(5,5)}
non_corner_d2 = {v for v in all_d2_verts if v not in corner_set}
print(f"Non-corner vertices that became degree-2: {sorted(non_corner_d2)}")

# For configs with Q=3, check deg2 vertices
print(f"\n--- Degree-2 vertices in Q=3 configs ---")
q3_d2_sets = []
for trial in range(500):
    k_total = random.randint(51, 64)
    sub = random.sample(range(len(wrap_indices)), k_total)
    Q, d2, _ = compute_Q(make_mask(sub))
    if Q == 3:
        q3_d2_sets.append(tuple(sorted(d2)))

q3_d2_unique = set(q3_d2_sets)
print(f"Q=3 configs found: {len(q3_d2_sets)}")
print(f"Unique deg2 vertex sets in Q=3 configs: {len(q3_d2_unique)}")
for vset in sorted(q3_d2_unique):
    count = q3_d2_sets.count(vset)
    print(f"  {list(vset)}: appeared {count}/{len(q3_d2_sets)} times")


# ============================================================
# 3. FINE SYMMETRY TEST: Q(diag + S_x) == Q(diag + mirror(S_x))
# ============================================================

print("\n" + "=" * 70)
print("3. FINE SYMMETRY TEST: Q(diag + S_x) == Q(diag + mirror(S_x))")
print("=" * 70)

sym_violations = 0
sym_tests = 0

for k in range(1, min(15, len(x_only_idx)+1)):
    n_test = min(100, math.comb(len(x_only_idx), k) if k <= 10 else 100)
    
    tested = set()
    local_violations = 0
    
    for _ in range(n_test):
        sx = tuple(sorted(random.sample(x_only_idx, k)))
        if sx in tested:
            continue
        tested.add(sx)
        
        # Mirror: map each x-only to its y-only counterpart
        sy_mirror = tuple(sorted([mirror_map[xi] for xi in sx if xi in mirror_map]))
        if len(sy_mirror) != k:
            continue  # skip if mirror not found for all
        
        Q_x, d2_x, _ = compute_Q(make_mask(diag_idx + list(sx)))
        Q_y, d2_y, _ = compute_Q(make_mask(diag_idx + list(sy_mirror)))
        
        sym_tests += 1
        if Q_x != Q_y:
            sym_violations += 1
            local_violations += 1
            if sym_violations <= 3:
                print(f"  VIOLATION at k={k}:")
                print(f"    S_x={sx} → Q={Q_x}, deg2={d2_x}")
                print(f"    mirror(S_x)={sy_mirror} → Q={Q_y}, deg2={d2_y}")
    
    if local_violations > 0:
        print(f"  k={k}: {local_violations}/{len(tested)} violations")
    else:
        print(f"  k={k}: {len(tested)} tests, all symmetric ✓")

print(f"\nTotal: {sym_violations}/{sym_tests} violations")
if sym_violations == 0:
    print("→ Symmetry r↔c perfectly preserved even with diag-wraps present")
else:
    print("→ Symmetry r↔c BROKEN in presence of diag-wraps!")


# ============================================================
# 4. COMPARISON TABLE 4×4 vs 6×6
# ============================================================

print("\n" + "=" * 70)
print("4. COMPARISON TABLE 4×4 vs 6×6")
print("=" * 70)

# 4×4 values (from previous experiment)
# Q(só diag removidas) = 3 (all 8 wraps are diag in 4×4, and Q=3 when all removed)
# Wait: in the 4×4, removing all 8 diag-wraps = removing all wraps = plane = Q=3

# 6×6 values
# Q(só diag removidas) = computed above
# wraps mínimas para Q=3: from greedy

# Find actual minimum via more thorough search
print("\nSearching for minimum wraps for Q=3 more carefully...")

# We know from previous experiment: first Q=3 at k=51 in random sampling
# Let's be more precise: try all diag + subsets of x+y of increasing size
min_q3_total = None
xy_pool = x_only_idx + y_only_idx

for k_extra in range(1, 57):
    found_q3 = False
    n_try = min(500, max(1, math.comb(len(xy_pool), k_extra)))  
    
    for _ in range(n_try):
        extra = random.sample(range(len(xy_pool)), k_extra)
        extra_wi = [xy_pool[i] for i in extra]
        removed = diag_idx + extra_wi
        Q, d2, _ = compute_Q(make_mask(removed))
        if Q == 3:
            found_q3 = True
            if min_q3_total is None or len(removed) < min_q3_total[0]:
                cats_rm = [wrap_cats[wi] for wi in removed]
                min_q3_total = (len(removed), cats_rm.count('x-only'), 
                               cats_rm.count('y-only'), cats_rm.count('diag'),
                               len(d2), d2, removed)
            break
    
    if found_q3:
        print(f"  Q=3 found at k_extra={k_extra} (total={8+k_extra})")
        break

if min_q3_total:
    total, cx, cy, cd, d2c, d2v, rm = min_q3_total
    print(f"  Minimum found: total={total}, x={cx}, y={cy}, diag={cd}")
    print(f"  deg2 vertices ({d2c}): {d2v}")

# Now the comparison table
print(f"\n{'='*70}")
print(f"| {'métrica':<40s} | {'4×4':>8s} | {'6×6':>8s} |")
print(f"|{'-'*42}|{'-'*10}|{'-'*10}|")
print(f"| {'total wraps':<40s} | {'8':>8s} | {'64':>8s} |")
print(f"| {'wraps diag':<40s} | {'8':>8s} | {'8':>8s} |")
print(f"| {'Q(só diag removidas)':<40s} | {'3':>8s} | {str(Q_diag):>8s} |")

min_wraps_q3_6x6 = len(current_set) if Q0 == 3 else '?'
print(f"| {'wraps mínimas para Q=3 (greedy)':<40s} | {'8':>8s} | {str(min_wraps_q3_6x6):>8s} |")

if isinstance(min_wraps_q3_6x6, int):
    pct = f"{100*min_wraps_q3_6x6/64:.0f}%"
else:
    pct = '?'
print(f"| {'wraps mínimas / total wraps':<40s} | {'100%':>8s} | {pct:>8s} |")

# deg2 when Q=3 in plane
Q_plane, d2_plane, _ = compute_Q(make_mask(list(range(len(wrap_indices)))))
print(f"| {'deg2_count quando Q=3 (plano)':<40s} | {'4':>8s} | {str(len(d2_plane)):>8s} |")

is_corners_4x4 = "sim"
is_corners_6x6 = "sim" if set(d2_plane) == corner_set else "não"
print(f"| {'vértices de grau 2 são os cantos?':<40s} | {is_corners_4x4:>8s} | {is_corners_6x6:>8s} |")


# ============================================================
# 5. CORNER INVARIANCE
# ============================================================

print("\n" + "=" * 70)
print("5. CORNER INVARIANCE ANALYSIS")
print("=" * 70)

print(f"\nPlane 6×6 degree-2 vertices: {d2_plane}")
print(f"Are they the 4 corners? {set(d2_plane) == corner_set}")

# Check degrees of corners and near-corners on the torus
print(f"\nTorus 6×6 vertex degrees (all vertices are degree 8 on full torus):")
print(f"After removing ONLY diag wraps:")
for r in range(N):
    for c in range(N):
        v = vid(r, c)
        d = degs_diag.get(v, 0)
        if d < 8:
            print(f"  ({r},{c}): degree = {d}  {'← CORNER' if (r,c) in corner_set else ''}")

# How many wraps touch each corner?
print(f"\nWrap edges incident to each corner:")
for r, c in corners:
    v = vid(r, c)
    incident_wraps = []
    for wi, ei in enumerate(wrap_indices):
        e = all_edges[ei]
        if e[2] == v or e[3] == v:
            incident_wraps.append((wi, e[4], e[0], e[1]))
    print(f"  ({r},{c}): {len(incident_wraps)} wrap edges")
    for wi, cat, src, dst in incident_wraps:
        print(f"    wrap[{wi:2d}] {cat:7s} {src}→{dst}")

# Check: after removing all wraps incident to each corner
print(f"\n--- What happens when we remove all wraps of a single corner? ---")
for r, c in corners:
    v = vid(r, c)
    corner_wraps = [wi for wi, ei in enumerate(wrap_indices) 
                    if all_edges[ei][2] == v or all_edges[ei][3] == v]
    Q, d2, degs = compute_Q(make_mask(corner_wraps))
    d_corner = degs.get(v, 0)
    print(f"  Remove all wraps of ({r},{c}): {len(corner_wraps)} wraps removed → "
          f"Q={Q}, deg({r},{c})={d_corner}, deg2={d2}")


# ============================================================
# FINAL ANSWER
# ============================================================

print("\n" + "=" * 70)
print("RESPOSTA FINAL")
print("=" * 70)

print(f"""
PERGUNTA: As 8 diag-wraps sozinhas produzem Q>0?
RESPOSTA: Q(só diag) = {Q_diag}
  → {'SIM, as diag-wraps sozinhas já iniciam a obstrução!' if Q_diag > 0 else 'NÃO, as diag-wraps sozinhas NÃO produzem Q>0.'}
  → {f'Q={Q_diag}, produz {len(deg2_diag)} vértices de grau 2' if Q_diag > 0 else f'deg2_count={len(deg2_diag)}, nenhum vértice isolado'}

{'O que falta para Q=3:' if Q_diag < 3 else ''}
{f'  Greedy encontrou Q=3 com {step} adições além das 8 diag (total={len(current_set)} wraps)' if Q_diag < 3 else ''}
{f'  Mix final: x-only={cx}, y-only={cy}, diag={cd}' if Q_diag < 3 else ''}
""")

print("✓ Análise cirúrgica completa.")
