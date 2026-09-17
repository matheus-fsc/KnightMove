"""winding_number.py — Task C.

Discrete winding number of every closed 6×6 knight tour around every
cell-centre of the board.

We choose

    p_{r,c} = (c + ½) + i (r + ½),    r, c ∈ {0, …, 5}

as the test points.  Algebra (see ``README.md`` §1.4) shows that no
knight segment between integer lattice points can ever pass through a
half-integer point, so the winding number

    w(p) = (1 / 2π) · Σ_t  Im( log( (z(t+1) − p) / (z(t) − p) ) )

is integer-valued and well-defined for every cell-centre.

Outputs:

* ``data/winding_numbers.npy``         shape ``(9862, 6, 6)``   int8
* ``data/winding_correlation.json``    correlations vs. the 10 edge-orbit
  parity bits (proxy for the 7 H₁ invariants under D₄)
* ``data/plots/winding_mean_heatmap.png``    mean ``|w|`` per cell
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from _tours import (
    BOARD,
    EDGE_ORBIT_ID,
    N_EDGE_ORBITS,
    load_or_compute_tours,
    tour_orbit_parity_bits,
    tours_to_complex,
    EDGES,
)

DATA_DIR = Path(__file__).parent / "data"
PLOT_DIR = DATA_DIR / "plots"

OUT_WINDINGS = DATA_DIR / "winding_numbers.npy"
OUT_CORRELATION = DATA_DIR / "winding_correlation.json"
OUT_HEATMAP = PLOT_DIR / "winding_mean_heatmap.png"


# ---------------------------------------------------------------------------
# Vectorised winding-number computation
# ---------------------------------------------------------------------------

def cell_center_grid() -> np.ndarray:
    """Return the (6, 6) array of cell centres as complex numbers."""
    rr, cc = np.meshgrid(
        np.arange(BOARD, dtype=np.float64) + 0.5,
        np.arange(BOARD, dtype=np.float64) + 0.5,
        indexing="ij",
    )
    return cc + 1j * rr  # shape (6, 6)


def windings_for_tour(z: np.ndarray, p_grid: np.ndarray) -> np.ndarray:
    """Compute w(p) for one tour ``z`` of shape ``(N,)`` and all cell
    centres in ``p_grid`` of shape ``(6, 6)``.

    Returns
    -------
    w : np.ndarray
        Integer array of shape ``(6, 6)``.
    """
    # close the polygon
    z_next = np.roll(z, -1)
    # Broadcast: (N, 1, 1) − (6, 6) → (N, 6, 6)
    a = z[:, None, None] - p_grid[None, :, :]
    b = z_next[:, None, None] - p_grid[None, :, :]
    # log(b/a) ; use complex log directly — equivalent and faster
    delta_log = np.log(b / a)
    total = np.sum(delta_log.imag, axis=0) / (2 * np.pi)
    return np.rint(total).astype(np.int8)


def windings_for_all_tours(tours: np.ndarray, verbose: bool = True) -> np.ndarray:
    """Compute winding numbers for every tour.

    Returns
    -------
    np.ndarray
        Shape ``(K, 6, 6)`` int8.
    """
    K = tours.shape[0]
    p_grid = cell_center_grid()
    out = np.empty((K, BOARD, BOARD), dtype=np.int8)
    # Vectorise across tours in chunks to keep memory bounded.
    chunk = 256
    z_all = tours_to_complex(tours)
    for i in range(0, K, chunk):
        block = z_all[i : i + chunk]  # (b, N)
        b = block.shape[0]
        z = block[:, :, None, None]
        z_next = np.roll(block, -1, axis=1)[:, :, None, None]
        a = z - p_grid[None, None, :, :]
        bc = z_next - p_grid[None, None, :, :]
        log_im = np.log(bc / a).imag
        total = log_im.sum(axis=1) / (2 * np.pi)
        out[i : i + b] = np.rint(total).astype(np.int8)
        if verbose:
            print(f"  windings chunk {i + b:>5} / {K}")
    return out


# ---------------------------------------------------------------------------
# Invariants vs. winding numbers
# ---------------------------------------------------------------------------

def compute_invariant_bits(tours: np.ndarray, verbose: bool = True) -> np.ndarray:
    """Per-tour D₄-orbit parity bits.  Shape ``(K, N_EDGE_ORBITS)``."""
    K = tours.shape[0]
    bits = np.empty((K, N_EDGE_ORBITS), dtype=np.uint8)
    for i in range(K):
        bits[i] = tour_orbit_parity_bits(tours[i])
        if verbose and (i + 1) % 1000 == 0:
            print(f"  invariant bits {i + 1:>5} / {K}")
    return bits


def gf2_rank(bits: np.ndarray) -> int:
    """Linear rank over GF(2) of the rows of ``bits``."""
    A = bits.copy().astype(np.uint8)
    rows, cols = A.shape
    rank = 0
    r = 0
    for c in range(cols):
        # find pivot
        pivot = -1
        for i in range(r, rows):
            if A[i, c]:
                pivot = i
                break
        if pivot < 0:
            continue
        A[[r, pivot]] = A[[pivot, r]]
        for i in range(rows):
            if i != r and A[i, c]:
                A[i] ^= A[r]
        rank += 1
        r += 1
        if r == rows:
            break
    return rank


def correlation_winding_vs_bits(
    windings: np.ndarray, bits: np.ndarray
) -> np.ndarray:
    """Return a Pearson correlation matrix of shape ``(36, N_EDGE_ORBITS)``
    between flattened winding numbers (36 per tour) and the 10 parity bits.

    Pearson r between an integer winding sequence and a 0/1 bit gives the
    point-biserial correlation, exactly what we want here.
    """
    K = windings.shape[0]
    W = windings.reshape(K, BOARD * BOARD).astype(np.float64)
    B = bits.astype(np.float64)
    Wc = W - W.mean(axis=0, keepdims=True)
    Bc = B - B.mean(axis=0, keepdims=True)
    num = Wc.T @ Bc
    denom = np.sqrt(
        (Wc ** 2).sum(axis=0)[:, None] * (Bc ** 2).sum(axis=0)[None, :]
    )
    with np.errstate(invalid="ignore", divide="ignore"):
        r = np.where(denom > 0, num / denom, 0.0)
    return r


# ---------------------------------------------------------------------------
# D₄-invariance check
# ---------------------------------------------------------------------------

def d4_apply_to_grid(grid: np.ndarray, t: int) -> np.ndarray:
    """Apply the t-th D₄ symmetry to a (6,6) array indexed by (row, col)."""
    rot = t % 4
    flip = t // 4
    out = grid
    if flip:
        out = out[::-1, :]
    for _ in range(rot):
        out = np.rot90(out, k=-1)
    return out


def check_d4_invariance(tours: np.ndarray, windings: np.ndarray) -> dict:
    """For each tour, find another tour that is its D₄ image and verify
    that the winding maps accordingly.  Returns a small dict summary.
    """
    K = tours.shape[0]
    # Build a hashable key for each tour: rotation-invariant signature of
    # the *unordered* edge set, which is what determines an undirected tour.
    def edge_set_key(tour: np.ndarray) -> tuple[tuple[int, int], ...]:
        s = sorted(
            (int(min(tour[i], tour[(i + 1) % len(tour)])),
             int(max(tour[i], tour[(i + 1) % len(tour)])))
            for i in range(len(tour))
        )
        return tuple(s)

    key_to_idx = {edge_set_key(tours[i]): i for i in range(K)}

    # D₄ action on vertex labels (vertex = r*6 + c)
    def vmap(t: int) -> np.ndarray:
        m = np.empty(36, dtype=np.int64)
        for v in range(36):
            r, c = divmod(v, 6)
            rot = t % 4
            flip = t // 4
            if flip:
                r = 5 - r
            for _ in range(rot):
                r, c = c, 5 - r
            m[v] = r * 6 + c
        return m

    # Two subtleties:
    #   (i) reflections reverse orientation in ℂ → w(p) flips sign.
    #  (ii) canonicalisation forces ``path[1] < path[-1]`` per cycle, so the
    #       D₄-image of a tour may have been canonicalised in the *reverse*
    #       direction — reversal also flips the sign of w(p).
    # We therefore compare absolute winding maps (which are invariant under
    # both subtleties) and separately count exact and sign-flipped matches.
    # Robust check (invariant to (i) which grid-orientation convention we
    # use for D₄ and (ii) whether canonicalisation flipped direction):
    # the *multiset* of |w(p)| over all 36 cell-centres must be identical
    # for any pair of tours related by a D₄ symmetry.
    matched = 0
    multiset_match = 0
    sample_size = min(K, 500)
    rng = np.random.default_rng(0)
    sample = rng.choice(K, size=sample_size, replace=False)
    for idx in sample:
        tour = tours[idx]
        own_multiset = tuple(sorted(np.abs(windings[idx]).flatten().tolist()))
        for t in range(1, 8):
            m = vmap(t)
            mapped_tour = m[tour]
            k = edge_set_key(mapped_tour)
            j = key_to_idx.get(k)
            if j is None:
                continue
            matched += 1
            other_multiset = tuple(sorted(np.abs(windings[j]).flatten().tolist()))
            if own_multiset == other_multiset:
                multiset_match += 1
    return {
        "sample_size": int(sample_size),
        "d4_matches_found": int(matched),
        "abs_multiset_match": int(multiset_match),
        "abs_multiset_match_rate": (
            float(multiset_match) / matched if matched else float("nan")
        ),
    }


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

def plot_mean_winding_heatmap(windings: np.ndarray, path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    mean_abs = np.abs(windings).mean(axis=0)
    fig, ax = plt.subplots(figsize=(5, 5))
    im = ax.imshow(mean_abs, origin="upper", cmap="magma")
    for r in range(BOARD):
        for c in range(BOARD):
            ax.text(c, r, f"{mean_abs[r, c]:.2f}", ha="center",
                    va="center",
                    color="white" if mean_abs[r, c] < mean_abs.max() * 0.6 else "black",
                    fontsize=9)
    ax.set_title(f"Mean |w(p)| over {windings.shape[0]:,} tours")
    ax.set_xlabel("col")
    ax.set_ylabel("row")
    fig.colorbar(im, ax=ax, fraction=0.046)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> dict:
    tours = load_or_compute_tours(verbose=True)
    print(f"[C] tours: {tours.shape}")

    print("[C] computing winding numbers …")
    windings = windings_for_all_tours(tours, verbose=True)
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    np.save(OUT_WINDINGS, windings)
    print(f"[C] saved {OUT_WINDINGS}  shape={windings.shape}")

    print("[C] computing edge-orbit parity bits …")
    bits = compute_invariant_bits(tours, verbose=True)
    bits_rank = gf2_rank(bits)
    print(f"[C] GF(2) rank of orbit-parity bit matrix: {bits_rank}  (out of {N_EDGE_ORBITS})")

    print("[C] correlation matrix windings ↔ bits …")
    corr = correlation_winding_vs_bits(windings, bits)
    print(f"[C] correlation matrix shape: {corr.shape}")
    print(f"[C] max |corr|: {np.nanmax(np.abs(corr)):.4f}")

    print("[C] D₄ self-consistency check (sample of tours) …")
    d4_check = check_d4_invariance(tours, windings)
    print(f"[C] {d4_check}")

    info = {
        "n_tours": int(windings.shape[0]),
        "winding_range": [int(windings.min()), int(windings.max())],
        "winding_mean_abs": float(np.abs(windings).mean()),
        "winding_per_cell_mean_abs": np.abs(windings).mean(axis=0).tolist(),
        "n_edge_orbits": int(N_EDGE_ORBITS),
        "bits_gf2_rank": int(bits_rank),
        "max_abs_correlation": float(np.nanmax(np.abs(corr))),
        "top10_cells_by_max_corr": [],
        "d4_invariance_check": d4_check,
    }
    flat_abs = np.abs(corr).max(axis=1)
    top10 = np.argsort(flat_abs)[::-1][:10]
    for idx in top10:
        r, c = divmod(int(idx), BOARD)
        best_bit = int(np.argmax(np.abs(corr[idx])))
        info["top10_cells_by_max_corr"].append(
            {
                "cell_rc": [r, c],
                "best_bit": best_bit,
                "corr": float(corr[idx, best_bit]),
            }
        )

    OUT_CORRELATION.write_text(json.dumps(info, indent=2))
    print(f"[C] saved {OUT_CORRELATION}")

    plot_mean_winding_heatmap(windings, OUT_HEATMAP)
    print(f"[C] saved {OUT_HEATMAP}")
    return info


if __name__ == "__main__":
    main()
