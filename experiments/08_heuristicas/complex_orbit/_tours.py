"""_tours.py — internal helper.

Exhaustive enumeration of all closed knight tours (Hamiltonian cycles)
on a 6×6 board, plus a small cache to ``data/tours_closed_6x6.npy``.

Backtracking is adapted from ``../cavalo_loop_destruicao_6x6.py`` and
returns **canonical** directed cycles: of the two directions in which an
undirected Hamiltonian cycle can be traversed starting at vertex 0, we
keep the one where ``path[1] < path[-1]`` so each undirected tour is
counted exactly once.

Output shape: ``(9862, 36)`` int8 array of vertex indices ``0..35``.
"""
from __future__ import annotations

import time
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# Board + knight graph
# ---------------------------------------------------------------------------

BOARD: int = 6
TOTAL: int = BOARD * BOARD
N: int = TOTAL  # alias for clarity in DFT code

_KNIGHT_MOVES_RC: tuple[tuple[int, int], ...] = (
    (2, 1), (2, -1), (-2, 1), (-2, -1),
    (1, 2), (1, -2), (-1, 2), (-1, -2),
)


def _vid(r: int, c: int) -> int:
    return r * BOARD + c


def _vrc(v: int) -> tuple[int, int]:
    return divmod(v, BOARD)


def _build_adj() -> tuple[list[list[int]], list[tuple[int, int]]]:
    adj: list[list[int]] = [[] for _ in range(TOTAL)]
    edges_set: set[tuple[int, int]] = set()
    for r in range(BOARD):
        for c in range(BOARD):
            v = _vid(r, c)
            for dr, dc in _KNIGHT_MOVES_RC:
                nr, nc = r + dr, c + dc
                if 0 <= nr < BOARD and 0 <= nc < BOARD:
                    u = _vid(nr, nc)
                    adj[v].append(u)
                    if u > v:
                        edges_set.add((v, u))
    return adj, sorted(edges_set)


ADJ, EDGES = _build_adj()


def vertex_label(v: int) -> str:
    r, c = _vrc(v)
    return f"{chr(65 + c)}{BOARD - r}"


# ---------------------------------------------------------------------------
# Exhaustive search
# ---------------------------------------------------------------------------

def _enumerate_directed_cycles(start: int = 0, verbose: bool = True) -> list[tuple[int, ...]]:
    """All directed Hamiltonian cycles anchored at ``start`` (every undirected
    cycle appears twice — once per direction).  Warnsdorff-ordered backtracking.
    """
    solutions: list[tuple[int, ...]] = []
    path = [start]
    alive = bytearray([1] * TOTAL)
    alive[start] = 0
    deg = [sum(1 for u in ADJ[v] if alive[u]) for v in range(TOTAL)]
    counter = [0]
    t0 = time.time()

    def backtrack() -> None:
        counter[0] += 1
        if verbose and counter[0] % 1_000_000 == 0:
            print(f"    [enumeration] backtracks={counter[0]:>11,}  "
                  f"sols={len(solutions):>6,}  depth={len(path)}/{TOTAL}  "
                  f"elapsed={time.time()-t0:5.1f}s")
        current = path[-1]
        n = len(path)
        if n == TOTAL:
            if start in ADJ[current]:
                solutions.append(tuple(path))
            return
        candidates = [u for u in ADJ[current] if alive[u]]
        if not candidates:
            return
        candidates.sort(key=lambda u: deg[u])
        # Warnsdorff prune: deg-0 candidate is only viable on the closing step.
        if deg[candidates[0]] == 0 and n < TOTAL - 1:
            return
        for nxt in candidates:
            alive[nxt] = 0
            for w in ADJ[nxt]:
                deg[w] -= 1
            path.append(nxt)
            backtrack()
            path.pop()
            alive[nxt] = 1
            for w in ADJ[nxt]:
                deg[w] += 1

    backtrack()
    return solutions


def _canonicalize_undirected(cycles: list[tuple[int, ...]]) -> np.ndarray:
    """From the directed list, keep one direction per undirected cycle.

    Convention: keep the direction where ``path[1] < path[-1]``.
    """
    seen: set[tuple[int, ...]] = set()
    canonical: list[tuple[int, ...]] = []
    for cyc in cycles:
        # rotate so 0 is at index 0 — already true, since we anchor at vertex 0.
        forward = cyc
        # Reversed direction with same start:
        reverse = (cyc[0],) + tuple(reversed(cyc[1:]))
        # Take whichever has smaller second-vertex; on ties take min lex.
        key_f = (forward[1], forward)
        key_r = (reverse[1], reverse)
        canon = forward if key_f < key_r else reverse
        if canon not in seen:
            seen.add(canon)
            canonical.append(canon)
    canonical.sort()
    return np.asarray(canonical, dtype=np.int8)


# ---------------------------------------------------------------------------
# Public API with on-disk cache
# ---------------------------------------------------------------------------

_DATA_DIR = Path(__file__).parent / "data"
_CACHE_PATH = _DATA_DIR / "tours_closed_6x6.npy"


def load_or_compute_tours(verbose: bool = True) -> np.ndarray:
    """Return all 9,862 closed tours of the 6×6 board as an int8 array of
    shape ``(9862, 36)``.  Cached to ``data/tours_closed_6x6.npy``.
    """
    if _CACHE_PATH.exists():
        if verbose:
            print(f"[tours] loading cached tours from {_CACHE_PATH}")
        return np.load(_CACHE_PATH)
    _DATA_DIR.mkdir(parents=True, exist_ok=True)
    if verbose:
        print("[tours] no cache — running exhaustive backtracking on 6×6 board")
        print("        (expected: 19,724 directed cycles → 9,862 canonical)")
    directed = _enumerate_directed_cycles(start=0, verbose=verbose)
    if verbose:
        print(f"[tours] directed cycles found: {len(directed):,}")
    canonical = _canonicalize_undirected(directed)
    if verbose:
        print(f"[tours] canonical (undirected) tours: {len(canonical):,}")
    np.save(_CACHE_PATH, canonical)
    if verbose:
        print(f"[tours] cached to {_CACHE_PATH}")
    return canonical


def tour_to_complex(tour: np.ndarray) -> np.ndarray:
    """Convert a tour (int vertex sequence of length 36) to a complex array
    ``z(t) = col + i·row`` of length 36.
    """
    rows, cols = np.divmod(tour.astype(np.int64), BOARD)
    return cols.astype(np.complex128) + 1j * rows.astype(np.complex128)


def tours_to_complex(tours: np.ndarray) -> np.ndarray:
    """Vectorized version: shape ``(K, 36) → (K, 36)`` complex."""
    rows, cols = np.divmod(tours.astype(np.int64), BOARD)
    return cols.astype(np.complex128) + 1j * rows.astype(np.complex128)


# ---------------------------------------------------------------------------
# Edge-orbit parity bits (proxy for the 7 H₁ invariants under D₄)
# ---------------------------------------------------------------------------

def _d4_apply_rc(rc: tuple[int, int], t: int) -> tuple[int, int]:
    r, c = rc
    rot = t % 4
    flip = t // 4
    if flip:
        r = BOARD - 1 - r
    for _ in range(rot):
        r, c = c, BOARD - 1 - r
    return (r, c)


def _edge_canonical(u: int, v: int) -> tuple[int, int]:
    return (u, v) if u < v else (v, u)


def _edge_orbit_index() -> tuple[np.ndarray, int]:
    """Build the edge → orbit-id map.  Returns ``(idx[80], n_orbits)``."""
    edge_to_orbit = {e: -1 for e in EDGES}
    next_id = 0
    for e in EDGES:
        if edge_to_orbit[e] != -1:
            continue
        orbit_id = next_id
        next_id += 1
        u, v = e
        ur, uc = _vrc(u)
        vr, vc = _vrc(v)
        for t in range(8):
            ur2, uc2 = _d4_apply_rc((ur, uc), t)
            vr2, vc2 = _d4_apply_rc((vr, vc), t)
            ue = _vid(ur2, uc2)
            ve = _vid(vr2, vc2)
            edge_to_orbit[_edge_canonical(ue, ve)] = orbit_id
    idx = np.array([edge_to_orbit[e] for e in EDGES], dtype=np.int32)
    return idx, next_id


EDGE_ORBIT_ID, N_EDGE_ORBITS = _edge_orbit_index()


def tour_edge_vector(tour: np.ndarray) -> np.ndarray:
    """Return the 0/1 edge-incidence vector of a closed tour (length = 80)."""
    edge_pos = {e: i for i, e in enumerate(EDGES)}
    vec = np.zeros(len(EDGES), dtype=np.uint8)
    for t in range(TOTAL):
        u = int(tour[t])
        v = int(tour[(t + 1) % TOTAL])
        vec[edge_pos[_edge_canonical(u, v)]] = 1
    return vec


def tour_orbit_parity_bits(tour: np.ndarray) -> np.ndarray:
    """For each of the ``N_EDGE_ORBITS`` D₄ orbits of edges, return the
    parity (XOR / GF(2)-sum) of the tour's edge-incidence inside that orbit.
    These are D₄-invariant per tour and are the natural proxy for the
    project's "7 H₁ invariant bits" — the bit-vector lives in GF(2)^10
    with rank 7 across all 9862 tours (verified empirically in Task C).
    """
    ev = tour_edge_vector(tour)
    bits = np.zeros(N_EDGE_ORBITS, dtype=np.uint8)
    for orbit_id in range(N_EDGE_ORBITS):
        mask = EDGE_ORBIT_ID == orbit_id
        bits[orbit_id] = int(ev[mask].sum()) & 1
    return bits


if __name__ == "__main__":
    tours = load_or_compute_tours(verbose=True)
    print(f"\nshape: {tours.shape}, dtype: {tours.dtype}")
    print(f"first tour (first 10 verts): {tours[0][:10]}  → {[vertex_label(v) for v in tours[0][:10]]}")
    print(f"edge orbits: {N_EDGE_ORBITS}  (sizes via bincount: {np.bincount(EDGE_ORBIT_ID).tolist()})")
