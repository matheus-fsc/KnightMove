#!/usr/bin/env python3
r"""
ADVERSARIAL re-derivation of the bulk-connectivity inductive proof.

This script does NOT import verify_bulk_connectivity.py. Everything is
re-derived from scratch to avoid a shared-bug blind spot. The point is not
to re-confirm the STATEMENT (known true to n=30) but to attack the PROOF:
does the *specific* set of edges the proof explicitly invokes suffice to
connect Bulk(n+2), assuming only the induction hypothesis (B\EC connected)?

Key adversarial test (KILLER TEST):
  Treat B\EC as a single connected blob (the IH, given as a black box).
  Add ONLY the edges the proof explicitly names:
     - each embedded corner -> its proof target,
     - each frame cell      -> its proof target.
  Then check every vertex of Bulk(n+2) is connected to the blob using ONLY
  those proof-edges. If yes for all applicable n, the proof's logic is sound
  independent of whether the real graph happens to be connected.
"""

from collections import deque

KNIGHT = [(1, 2), (2, 1), (2, -1), (1, -2),
          (-1, -2), (-2, -1), (-2, 1), (-1, 2)]


def is_knight(u, v):
    di, dj = abs(u[0] - v[0]), abs(u[1] - v[1])
    return di * dj == 2  # {1,2} or {2,1}


def on_board(m, c):
    return 0 <= c[0] < m and 0 <= c[1] < m


def knight_nbrs(m, c):
    return [(c[0] + di, c[1] + dj) for di, dj in KNIGHT
            if on_board(m, (c[0] + di, c[1] + dj))]


def big_corners(m):
    return {(0, 0), (0, m - 1), (m - 1, 0), (m - 1, m - 1)}


def mand_edges(m):
    """The 8 mandatory edges = edges at the 4 board corners."""
    E = set()
    for c in big_corners(m):
        for w in knight_nbrs(m, c):
            E.add(frozenset((c, w)))
    return E


# ---- the proof's EXPLICIT moves, re-derived independently here ----

def embedded_corner_move(n, ec):
    """Per-corner reattachment move committed to by the proof (D4 images of
    (1,1)->(2,3)). Returns the target cell in the (n+2)-board."""
    (i, j) = ec
    table = {
        (1, 1): (2, 3),
        (1, n): (3, n - 1),
        (n, n): (n - 1, n - 2),
        (n, 1): (n - 2, 2),
    }
    return table[ec]


def frame_move(n, cell):
    """The proof's explicit frame->B move template (D4 of the top-edge rule
    (0,k)->(2,k+1) for k<=n-1 else (2,n-1)). Coords in (n+2)-board, B=1..n."""
    m = n + 2
    i, j = cell
    if i == 0:                       # top
        k = j
        return (2, k + 1) if k <= n - 1 else (2, n - 1)
    if i == m - 1:                   # bottom
        k = j
        return (n - 1, k + 1) if k <= n - 1 else (n - 1, n - 1)
    if j == 0:                       # left
        k = i
        return (k + 1, 2) if k <= n - 1 else (n - 1, 2)
    if j == m - 1:                   # right
        k = i
        return (k + 1, n - 1) if k <= n - 1 else (n - 1, n - 1)
    raise ValueError("not a frame cell")


def regions(n):
    m = n + 2
    B = {(i, j) for i in range(1, n + 1) for j in range(1, n + 1)}
    EC = {(1, 1), (1, n), (n, 1), (n, n)}
    R = big_corners(m)
    F = {(i, j) for i in range(m) for j in range(m)
         if i in (0, m - 1) or j in (0, m - 1)}
    return m, B, EC, R, F


def check2_embedded(n):
    r"""CHECK 2: each embedded corner's explicit move is a knight move landing
    in B\EC."""
    m, B, EC, R, F = regions(n)
    results = {}
    for ec in [(1, 1), (1, n), (n, 1), (n, n)]:
        tgt = embedded_corner_move(n, ec)
        ok = (is_knight(ec, tgt) and on_board(m, tgt)
              and tgt in B and tgt not in EC)
        results[ec] = (ok, tgt)
    return results


def check3_mand_internal_to_B(n):
    """CHECK 3: no Mand(n+2) edge is internal to B; removing Mand deletes 0
    B-internal edges. Returns (#mand_internal_to_B, #mand_touching_B)."""
    m, B, EC, R, F = regions(n)
    M = mand_edges(m)
    internal = 0
    touching = 0
    for e in M:
        a, b = tuple(e)
        if a in B and b in B:
            internal += 1
        if a in B or b in B:
            touching += 1
    return internal, touching


def check5_frame_enum(n):
    r"""CHECK 5: frame enumeration. Returns counts and verifies F\R partitions
    into the 4 edges, |F\R| = 4n, corners are exactly R and excluded."""
    m, B, EC, R, F = regions(n)
    FminusR = F - R
    # the 4 edges as the proof slices them
    top = {(0, k) for k in range(1, n + 1)}
    bot = {(m - 1, k) for k in range(1, n + 1)}
    left = {(k, 0) for k in range(1, n + 1)}
    right = {(k, m - 1) for k in range(1, n + 1)}
    union = top | bot | left | right
    disjoint = (len(top) + len(bot) + len(left) + len(right) == len(union))
    covers = (union == FminusR)
    size_ok = (len(FminusR) == 4 * n)
    corners_in_F = R.issubset(F)
    return covers and disjoint and size_ok and corners_in_F, len(FminusR), 4 * n


def check6_moves_valid(n):
    """CHECK 6: every proof-move (frame + embedded) is a valid knight move,
    lands in B, target is a real Bulk(n+2) vertex (in B, not in R), and the
    connecting edge is NOT a Mand(n+2) edge (so it survives)."""
    m, B, EC, R, F = regions(n)
    M = mand_edges(m)
    bulk_vertices = (B | F) - R
    failures = []
    # frame moves
    for cell in (F - R):
        tgt = frame_move(n, cell)
        if not is_knight(cell, tgt):
            failures.append(("frame-not-knight", cell, tgt)); continue
        if tgt not in B:
            failures.append(("frame-not-in-B", cell, tgt)); continue
        if tgt not in bulk_vertices or cell not in bulk_vertices:
            failures.append(("frame-not-bulk-vertex", cell, tgt)); continue
        if frozenset((cell, tgt)) in M:
            failures.append(("frame-edge-is-mand", cell, tgt)); continue
    # embedded corner moves
    for ec in EC:
        tgt = embedded_corner_move(n, ec)
        if not is_knight(ec, tgt):
            failures.append(("ec-not-knight", ec, tgt)); continue
        if tgt not in B or tgt in EC:
            failures.append(("ec-not-in-B-minus-EC", ec, tgt)); continue
        if frozenset((ec, tgt)) in M:
            failures.append(("ec-edge-is-mand", ec, tgt)); continue
    return failures


def killer_test(n):
    r"""KILLER TEST: assume ONLY that B\EC is connected (collapse it to one
    blob). Add ONLY the proof's explicit edges. Check every Bulk(n+2) vertex
    reaches the blob. This validates the proof's logic, not the graph."""
    m, B, EC, R, F = regions(n)
    BLOB = "BLOB"  # represents the connected B\EC from IH
    B_minus_EC = B - EC

    adj = {v: set() for v in ((B | F) - R)}
    adj[BLOB] = set()

    def link(u, v):
        # map any B\EC endpoint to BLOB
        uu = BLOB if u in B_minus_EC else u
        vv = BLOB if v in B_minus_EC else v
        if uu == vv:
            return
        adj.setdefault(uu, set()).add(vv)
        adj.setdefault(vv, set()).add(uu)

    # proof edge set: embedded corners -> target ; frame -> target
    for ec in EC:
        link(ec, embedded_corner_move(n, ec))
    for cell in (F - R):
        link(cell, frame_move(n, cell))

    # BFS from BLOB over ONLY proof edges
    seen = {BLOB}
    dq = deque([BLOB])
    while dq:
        x = dq.popleft()
        for y in adj.get(x, ()):
            if y not in seen:
                seen.add(y); dq.append(y)

    # every real bulk vertex must be reached (B\EC folded into BLOB)
    targets = (EC | (F - R))  # vertices that must attach via proof edges
    unreached = [v for v in targets if v not in seen]
    # also confirm vertex count bookkeeping
    expected_vertices = m * m - 4
    real_vertices = len((B | F) - R)
    return unreached, real_vertices, expected_vertices


def real_bfs(n):
    """Independent BFS on the actual Bulk(n) for cross-check."""
    m = n  # here build G_n minus its corners directly
    C = big_corners(m)
    start = next(v for v in ((i, j) for i in range(m) for j in range(m))
                 if v not in C)
    seen = {start}; dq = deque([start])
    while dq:
        v = dq.popleft()
        for w in knight_nbrs(m, v):
            if w not in C and w not in seen:
                seen.add(w); dq.append(w)
    return len(seen), m * m - 4


def main():
    print("=== ADVERSARIAL RE-DERIVATION ===\n")

    even = [6, 8, 10, 12, 14, 16, 18, 20]
    odd = [7, 9, 11, 13, 15, 17, 19, 21]
    step_ns = sorted(even + odd)

    print("CHECK 2 — embedded-corner explicit moves (knight + in B\\EC):")
    c2_all = True
    for n in step_ns:
        res = check2_embedded(n)
        allok = all(ok for ok, _ in res.values())
        c2_all &= allok
        if n in (6, 7, 8):  # show the thin/base instances in detail
            detail = "  ".join(f"{ec}->{t}:{'OK' if ok else 'FAIL'}"
                               for ec, (ok, t) in res.items())
            print(f"  n={n}->{n+2}: {detail}")
    print(f"  ALL n in {step_ns[0]}..{step_ns[-1]} (step): "
          f"{'4/4 PASS each' if c2_all else 'FAILURES'}\n")

    print("CHECK 3 — Mand(n+2) edges internal to B (must be 0):")
    c3_all = True
    for n in step_ns:
        internal, touching = check3_mand_internal_to_B(n)
        c3_all &= (internal == 0)
        if n in (6, 7, 8, 20):
            print(f"  n={n}->{n+2}: internal_to_B={internal}  "
                  f"touching_B={touching} (8 mand edges, all corner->B)")
    print(f"  ALL: {'SOUND (0 internal)' if c3_all else 'GAP'}\n")

    print("CHECK 5 — frame enumeration (F\\R = 4 disjoint edges, size 4n):")
    c5_all = True
    for n in step_ns:
        ok, got, exp = check5_frame_enum(n)
        c5_all &= ok
        if n in (6, 7, 20):
            print(f"  n={n}->{n+2}: |F\\R|={got} expected 4n={exp}  "
                  f"{'COMPLETE' if ok else 'INCOMPLETE'}")
    print(f"  ALL: {'COMPLETE' if c5_all else 'MISSING CASES'}\n")

    print("CHECK 6 — every proof-move valid (knight + in B + not a Mand edge):")
    c6_all = True
    for n in step_ns:
        fails = check6_moves_valid(n)
        c6_all &= (len(fails) == 0)
        if fails:
            print(f"  n={n}->{n+2}: {len(fails)} FAIL e.g. {fails[:3]}")
    print(f"  ALL n in step range: "
          f"{'EVERY PROOF-MOVE VALID' if c6_all else 'FAILURES'}\n")

    print("KILLER TEST — proof's explicit edges suffice (assume only IH):")
    kt_all = True
    for n in step_ns:
        unreached, realv, expv = killer_test(n)
        ok = (len(unreached) == 0 and realv == expv)
        kt_all &= ok
        if n in (6, 7, 8) or not ok:
            print(f"  n={n}->{n+2}: unreached={len(unreached)}  "
                  f"vertices={realv}/{expv}  "
                  f"{'PROOF-LOGIC SOUND' if ok else 'PROOF-LOGIC GAP: '+str(unreached[:5])}")
    print(f"  ALL: {'PROOF EDGES SUFFICE EVERYWHERE' if kt_all else 'GAP'}\n")

    print("CROSS-CHECK — independent BFS on real Bulk(n):")
    bfs_all = True
    for n in [6, 7, 8, 9, 10, 15, 20, 21]:
        got, exp = real_bfs(n)
        ok = got == exp
        bfs_all &= ok
        print(f"  n={n}: reached {got}/{exp}  {'YES' if ok else 'NO'}")

    print()
    overall = c2_all and c3_all and c5_all and c6_all and kt_all and bfs_all
    print(f"OVERALL: {'ALL ADVERSARIAL CHECKS PASS' if overall else 'GAP FOUND'}")
    return 0 if overall else 1


if __name__ == "__main__":
    raise SystemExit(main())
