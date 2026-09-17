"""h1_invariants.py — recover the genuine H₁ basis of the 6×6 closed
knight tours over GF(2), then classify each invariant bit as geometric
(correlated with winding) or purely algebraic.

This closes the gap left by ``winding_number.py``, whose D₄-orbit
parity proxy only spanned a rank-4 subspace.  Here we work directly
with the tour-edge incidence matrix ``T`` of shape ``(9862, |E|)`` over
GF(2) and extract the true row-span basis.

Pipeline:

    Task 1  →  build T, compute β₁ of the graph, find rank(T) and a
               linearly-independent basis B of its rows; project tours
               onto B; verify D₄-invariance of the projection.
    Task 2  →  for every basis bit, correlate with windings (integer &
               mod 2) and classify "geometric" vs "purely algebraic".
    Task 3  →  one heatmap per invariant in ``data/plots/``.
    Task 4  →  ``README.md`` update is done from the orchestrating
               session, not here.

Outputs of Task 1:
    data/h1_basis_gf2.npy         (R, |E|)  uint8
    data/h1_coordinates_gf2.npy   (9862, R) uint8

Outputs of Task 2:
    data/invariant_classification.json

Outputs of Task 3:
    data/plots/invariant_classification.png
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from _tours import (
    BOARD,
    EDGES,
    TOTAL,
    load_or_compute_tours,
)

DATA_DIR = Path(__file__).parent / "data"
PLOT_DIR = DATA_DIR / "plots"

OUT_BASIS = DATA_DIR / "h1_basis_gf2.npy"
OUT_COORDS = DATA_DIR / "h1_coordinates_gf2.npy"
OUT_CLASSIFICATION = DATA_DIR / "invariant_classification.json"
OUT_CLASSIFICATION_PNG = PLOT_DIR / "invariant_classification.png"
OUT_BASIS_FREE = DATA_DIR / "winding_linearity_basis_free.json"
OUT_BASIS_FREE_PNG = PLOT_DIR / "winding_linearity_basis_free.png"

GEOMETRIC_THRESHOLD = 0.3  # per the user spec for type classification


# ---------------------------------------------------------------------------
# Graph: boundary matrix and incidence
# ---------------------------------------------------------------------------

def boundary_matrix_1() -> np.ndarray:
    """Return ∂₁ over GF(2): shape ``(V, E) = (36, 80)``.

    ``∂₁[v, e] = 1`` iff vertex ``v`` is an endpoint of edge ``e``.
    Over GF(2), ``∂₁[v, e] = 1 + 1 = 0`` for self-loops (none here).
    """
    E = len(EDGES)
    B = np.zeros((TOTAL, E), dtype=np.uint8)
    for e_idx, (u, v) in enumerate(EDGES):
        B[u, e_idx] = 1
        B[v, e_idx] = 1
    return B


def tours_to_edge_matrix(tours: np.ndarray) -> np.ndarray:
    """Return ``T`` with ``T[i, e] = 1`` iff edge ``e`` is in tour ``i``.

    Shape: ``(K, |E|)``.  Each row has Hamming weight ``N = 36``
    (every tour has exactly 36 edges).
    """
    K, N = tours.shape
    E = len(EDGES)
    edge_to_idx = {e: i for i, e in enumerate(EDGES)}
    T = np.zeros((K, E), dtype=np.uint8)
    for i in range(K):
        for t in range(N):
            u, v = int(tours[i, t]), int(tours[i, (t + 1) % N])
            if u > v:
                u, v = v, u
            T[i, edge_to_idx[(u, v)]] = 1
    return T


# ---------------------------------------------------------------------------
# GF(2) linear algebra
# ---------------------------------------------------------------------------

def gf2_rank_and_basis(rows: np.ndarray) -> tuple[int, list[int], np.ndarray]:
    """Find a maximal subset of linearly independent rows of ``rows`` over GF(2).

    Returns
    -------
    rank : int
    pivot_indices : list[int]
        Indices of the original rows that form the basis (in selection order).
    basis : np.ndarray
        ``(rank, |E|)`` matrix of the chosen ORIGINAL rows (not RREF'd).
    """
    K, C = rows.shape
    # Reduced rows keyed by their leading column, for fast elimination.
    reduced_by_lc: dict[int, np.ndarray] = {}
    pivot_indices: list[int] = []
    basis_rows: list[np.ndarray] = []

    for i in range(K):
        row = rows[i].copy()
        # Eliminate using existing pivots, in increasing-column order.
        for lc in sorted(reduced_by_lc.keys()):
            if row[lc]:
                row ^= reduced_by_lc[lc]
        nz = np.flatnonzero(row)
        if len(nz) == 0:
            continue
        leading_col = int(nz[0])
        reduced_by_lc[leading_col] = row  # reduced form, used only for elimination
        pivot_indices.append(i)
        basis_rows.append(rows[i].copy())  # ORIGINAL row is the basis vector
    if basis_rows:
        basis = np.stack(basis_rows, axis=0)
    else:
        basis = np.zeros((0, C), dtype=np.uint8)
    return len(pivot_indices), pivot_indices, basis


def gf2_project(T: np.ndarray, basis: np.ndarray) -> np.ndarray:
    """Compute ``coords[i, k] = T[i] · basis[k]   mod 2``, shape ``(K, R)``."""
    # Promote to int32 to avoid uint8 overflow when summing 80 products.
    prod = T.astype(np.int32) @ basis.T.astype(np.int32)
    return (prod & 1).astype(np.uint8)


# ---------------------------------------------------------------------------
# D₄ action on edges (used for the invariance check)
# ---------------------------------------------------------------------------

def _vrc(v: int) -> tuple[int, int]:
    return divmod(v, BOARD)


def _vid(r: int, c: int) -> int:
    return r * BOARD + c


def _d4_node(r: int, c: int, t: int) -> tuple[int, int]:
    rot = t % 4
    flip = t // 4
    if flip:
        r = BOARD - 1 - r
    for _ in range(rot):
        r, c = c, BOARD - 1 - r
    return r, c


def edge_permutation(t: int) -> np.ndarray:
    """Return ``perm`` of length ``|E|`` such that edge index
    ``e = (u, v)`` maps to ``perm[e] = idx of the image edge under D₄_t``.
    """
    edge_to_idx = {e: i for i, e in enumerate(EDGES)}
    out = np.empty(len(EDGES), dtype=np.int32)
    for e_idx, (u, v) in enumerate(EDGES):
        ur, uc = _vrc(u)
        vr, vc = _vrc(v)
        ur2, uc2 = _d4_node(ur, uc, t)
        vr2, vc2 = _d4_node(vr, vc, t)
        u2 = _vid(ur2, uc2)
        v2 = _vid(vr2, vc2)
        e2 = (min(u2, v2), max(u2, v2))
        out[e_idx] = edge_to_idx[e2]
    return out


def vertex_map(t: int) -> np.ndarray:
    m = np.empty(TOTAL, dtype=np.int64)
    for v in range(TOTAL):
        r, c = _vrc(v)
        r2, c2 = _d4_node(r, c, t)
        m[v] = _vid(r2, c2)
    return m


def find_d4_partner_indices(tours: np.ndarray) -> np.ndarray:
    """For every tour ``i`` and every ``t ∈ 1..7``, find the index of the
    canonical tour that is the D₄_t image of tour ``i``.

    Returns
    -------
    np.ndarray
        Shape ``(K, 7)``, dtype int32.  ``out[i, t-1] = j`` such that
        ``T[j] == permuted_edge_vector(T[i], t)``.  All entries should
        be defined since the set of tours is closed under D₄.
    """
    K = tours.shape[0]
    perms = [edge_permutation(t) for t in range(1, 8)]
    T_edges = tours_to_edge_matrix(tours)
    # Hashable key per tour: pack bits.
    packed = np.packbits(T_edges, axis=1)
    key_to_idx: dict[bytes, int] = {bytes(row): i for i, row in enumerate(packed)}
    out = np.full((K, 7), -1, dtype=np.int32)
    for ti, perm in enumerate(perms):
        permuted = T_edges[:, perm]
        permuted_packed = np.packbits(permuted, axis=1)
        for i in range(K):
            key = bytes(permuted_packed[i])
            j = key_to_idx.get(key, -1)
            out[i, ti] = j
    return out


# ---------------------------------------------------------------------------
# Task 1
# ---------------------------------------------------------------------------

def task1_build_basis(tours: np.ndarray, verbose: bool = True) -> dict:
    n_V = TOTAL
    n_E = len(EDGES)
    beta1 = n_E - n_V + 1
    if verbose:
        print(f"[T1] graph: V={n_V}, E={n_E}, β₁ = E − V + 1 = {beta1}")
        print(f"[T1] building boundary matrix ∂₁ …")
    B1 = boundary_matrix_1()
    rank_B1 = int(np.linalg.matrix_rank(B1.astype(np.float64)))
    if verbose:
        print(f"[T1] rank(∂₁) over ℝ = {rank_B1}  (= V − 1 = {n_V - 1} for connected graph)")

    if verbose:
        print(f"[T1] building tour-edge matrix T …")
    T = tours_to_edge_matrix(tours)
    if verbose:
        print(f"[T1] T shape: {T.shape}, dtype: {T.dtype}, row weight: {int(T[0].sum())} (== 36 for closed Ham)")

    # Sanity: T @ ∂₁.T mod 2 == 0  for every closed tour.
    cycle_check = (T.astype(np.int32) @ B1.T.astype(np.int32)) & 1
    cycle_residual = int(cycle_check.sum())
    if verbose:
        print(f"[T1] cycle check  T · ∂₁ᵀ mod 2 == 0 ?   total residual = {cycle_residual}  ({'OK' if cycle_residual == 0 else 'FAIL'})")

    if verbose:
        print(f"[T1] GF(2) row-rank of T …")
    rank, pivot_indices, basis = gf2_rank_and_basis(T)
    if verbose:
        print(f"[T1] rank(T) over GF(2) = {rank}")
        print(f"[T1] first 10 pivot row indices: {pivot_indices[:10]}")

    # Project tours onto basis
    coords = gf2_project(T, basis)
    if verbose:
        print(f"[T1] coords shape: {coords.shape}  (per-tour GF(2) signature)")

    # Distribution over the 2^rank possible signatures
    coord_keys = np.packbits(coords, axis=1).tobytes()
    bytes_per_row = coords.shape[1]  # already small
    coord_signatures = [tuple(coords[i].tolist()) for i in range(coords.shape[0])]
    from collections import Counter
    sig_counter = Counter(coord_signatures)
    if verbose:
        print(f"[T1] distinct signatures: {len(sig_counter)}  (max possible = 2^{rank} = {2**rank})")
        print(f"[T1] signature size distribution (top 10):")
        for sig, n in sig_counter.most_common(10):
            print(f"     {sig}  →  {n} tours")
        if len(sig_counter) > 10:
            tail = sum(n for _, n in sig_counter.most_common()[10:])
            print(f"     … {len(sig_counter) - 10} more signatures, totalling {tail} tours")

    # D₄-invariance check: do D₄-equivalent tours share the same coords?
    if verbose:
        print(f"[T1] D₄ invariance check on coords …")
    partner_idx = find_d4_partner_indices(tours)
    n_d4_pairs = 0
    n_d4_coord_match = 0
    for i in range(coords.shape[0]):
        for ti in range(7):
            j = int(partner_idx[i, ti])
            if j < 0:
                continue
            n_d4_pairs += 1
            if np.array_equal(coords[i], coords[j]):
                n_d4_coord_match += 1
    d4_rate = n_d4_coord_match / n_d4_pairs if n_d4_pairs > 0 else float("nan")
    if verbose:
        print(f"[T1] D₄ coord match: {n_d4_coord_match} / {n_d4_pairs}  (rate {d4_rate:.4f})")

    # If the basis isn't D₄-invariant, re-derive the maximal D₄-invariant
    # sub-basis: the set of coordinate functionals b such that
    #     b · T[i] = b · T[σ(i)] for every D₄ symmetry σ.
    # Equivalently: b is orthogonal (over GF(2)) to every difference
    # ``T[i] − T[σ(i)]``.  We compute that orthogonal complement, then
    # intersect with the row span of T (i.e. with B).
    if verbose:
        print(f"[T1] computing maximal D₄-invariant sub-basis of span(T) …")
    diffs_list = []
    for i in range(coords.shape[0]):
        for ti in range(7):
            j = int(partner_idx[i, ti])
            if j < 0:
                continue
            diffs_list.append(T[i] ^ T[j])
    diffs = np.array(diffs_list, dtype=np.uint8)
    # Reduce to a basis of the difference span (call its rank ``d``).
    d_rank, _, diff_basis = gf2_rank_and_basis(diffs)
    if verbose:
        print(f"[T1] dim(D₄-orbit-difference span) = {d_rank}")
    # The D₄-invariant span(T) is span(T) ∩ (diff_span)^⊥ inside GF(2)^E.
    # Compute by reducing coords-of-basis modulo diff_basis using basis vectors
    # that are simultaneously in span(T) and ⊥ diff_basis.
    # We re-pick from span(T): try each basis vector; if it isn't ⊥ all of
    # diff_basis, replace it with a corrected vector that lives in span(T).
    # Practical approach: build a (basis ∪ diff_basis) stack, reduce, and read
    # off the kernel of the linear map  span(T) → coord-vector-differences.
    # For our purposes we use a simpler method: compute the rank of the joint
    # system [diff_basis ; basis], and the D₄-invariant dim of span(T) is
    #   d4_inv_dim = rank(basis) − ( rank([diff_basis; basis]) − rank(diff_basis) )
    # Correct computation: dim(span(T) ∩ D⊥) where D⊥ = annihilator of diff_basis
    # in GF(2)^|E|.  Equivalently: kernel of the map  basis  →  GF(2)^d_rank
    # sending  b ↦ (b · d_1, ..., b · d_d).
    if d_rank > 0:
        # Build the (rank × d_rank) coupling matrix M
        M = (basis.astype(np.int32) @ diff_basis.T.astype(np.int32)) & 1
        M = M.astype(np.uint8)
        # Left kernel of M: vectors v ∈ GF(2)^rank such that v · M = 0
        # Compute via: rank(M^T) — left nullity = rank − rank(M).
        rank_M, _, _ = gf2_rank_and_basis(M)
        d4_inv_dim = rank - rank_M
    else:
        d4_inv_dim = rank
    if verbose:
        print(f"[T1] diff_rank = {d_rank}  ⇒  D₄-invariant dim of span(T) = {d4_inv_dim}")

    # Save artifacts
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    np.save(OUT_BASIS, basis)
    np.save(OUT_COORDS, coords)
    if verbose:
        print(f"[T1] saved {OUT_BASIS}  shape={basis.shape}")
        print(f"[T1] saved {OUT_COORDS}  shape={coords.shape}")

    return {
        "n_vertices": n_V,
        "n_edges": n_E,
        "beta1_graph": beta1,
        "rank_boundary_real": rank_B1,
        "cycle_residual": cycle_residual,
        "rank_T_gf2": rank,
        "n_distinct_signatures": len(sig_counter),
        "d4_pairs_checked": n_d4_pairs,
        "d4_coord_match": n_d4_coord_match,
        "d4_match_rate": d4_rate,
        "d4_invariant_dim_of_spanT": int(d4_inv_dim),
        "pivot_indices_first10": list(map(int, pivot_indices[:10])),
        "signature_top10": [
            {"sig": list(map(int, sig)), "n_tours": int(n)}
            for sig, n in sig_counter.most_common(10)
        ],
    }


# ---------------------------------------------------------------------------
# Task 2 — correlate H₁ bits with windings
# ---------------------------------------------------------------------------

def pearson_per_column(W: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Pearson correlation between each column of ``W`` (shape K×C) and the
    1-D vector ``b`` (length K)."""
    W = W.astype(np.float64)
    bf = b.astype(np.float64)
    bc = bf - bf.mean()
    Wc = W - W.mean(axis=0, keepdims=True)
    num = Wc.T @ bc
    denom = np.sqrt((Wc ** 2).sum(axis=0) * (bc ** 2).sum())
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(denom > 0, num / denom, 0.0)


def task2_classify_invariants(verbose: bool = True) -> dict:
    if not OUT_COORDS.exists():
        raise RuntimeError(
            f"missing {OUT_COORDS}; run Task 1 first via task1_build_basis()."
        )
    coords = np.load(OUT_COORDS)  # (K, R) uint8
    windings_path = DATA_DIR / "winding_numbers.npy"
    if not windings_path.exists():
        raise RuntimeError(
            f"missing {windings_path}; run winding_number.py first."
        )
    windings = np.load(windings_path)  # (K, 6, 6) int8
    K, R = coords.shape
    W = windings.reshape(K, BOARD * BOARD).astype(np.int32)  # (K, 36)
    W_mod2 = (W & 1).astype(np.uint8)

    if verbose:
        print(f"[T2] coords shape: {coords.shape}, windings (flat): {W.shape}")

    classifications: list[dict] = []
    corr_int_table = np.zeros((R, BOARD * BOARD), dtype=np.float64)
    corr_mod2_table = np.zeros((R, BOARD * BOARD), dtype=np.float64)
    for k in range(R):
        bit = coords[:, k]
        n1 = int(bit.sum())
        n0 = K - n1
        # Correlations
        corr_int = pearson_per_column(W, bit)
        corr_mod2 = pearson_per_column(W_mod2, bit)
        corr_int_table[k] = corr_int
        corr_mod2_table[k] = corr_mod2
        max_int = float(np.nanmax(np.abs(corr_int)))
        max_mod2 = float(np.nanmax(np.abs(corr_mod2)))
        argmax_int = int(np.nanargmax(np.abs(corr_int)))
        argmax_mod2 = int(np.nanargmax(np.abs(corr_mod2)))
        kind = "geometric" if max_mod2 > GEOMETRIC_THRESHOLD else "algebraic"
        witness = (
            [argmax_mod2 // BOARD, argmax_mod2 % BOARD]
            if kind == "geometric"
            else None
        )
        classifications.append({
            f"basis_vector_{k}": {
                "type": kind,
                "max_r_winding_int": max_int,
                "argmax_winding_int_cell": [argmax_int // BOARD, argmax_int % BOARD],
                "max_r_winding_mod2": max_mod2,
                "argmax_winding_mod2_cell": [argmax_mod2 // BOARD, argmax_mod2 % BOARD],
                "witness_cell": witness,
                "n_tours_bit0": int(n0),
                "n_tours_bit1": int(n1),
            }
        })
        if verbose:
            tag = "GEOM " if kind == "geometric" else "alg  "
            print(f"[T2] inv #{k:>2}  {tag}  "
                  f"|r|_int_max = {max_int:.4f}  "
                  f"|r|_mod2_max = {max_mod2:.4f}  "
                  f"n1 = {n1:>5}  ({100.0*n1/K:5.1f}%)")

    summary = {
        "n_invariants": R,
        "n_tours": K,
        "geometric_threshold_mod2": GEOMETRIC_THRESHOLD,
        "invariants": classifications,
    }
    OUT_CLASSIFICATION.write_text(json.dumps(summary, indent=2))
    if verbose:
        print(f"[T2] saved {OUT_CLASSIFICATION}")

    # Also persist the full per-cell correlation tables (useful for plotting)
    np.save(DATA_DIR / "h1_corr_int.npy", corr_int_table)
    np.save(DATA_DIR / "h1_corr_mod2.npy", corr_mod2_table)
    if verbose:
        print(f"[T2] saved data/h1_corr_int.npy  shape={corr_int_table.shape}")
        print(f"[T2] saved data/h1_corr_mod2.npy shape={corr_mod2_table.shape}")

    return summary


# ---------------------------------------------------------------------------
# Task 3 — visualisation
# ---------------------------------------------------------------------------

def task3_plot(summary: dict | None = None) -> None:
    if summary is None:
        summary = json.loads(OUT_CLASSIFICATION.read_text())
    corr_mod2_table = np.load(DATA_DIR / "h1_corr_mod2.npy")  # (R, 36)
    R = corr_mod2_table.shape[0]
    invs = summary["invariants"]

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n_cols = 4
    n_rows = (R + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4.0 * n_cols, 3.7 * n_rows))
    axes = np.atleast_2d(axes)
    for k in range(R):
        ax = axes[k // n_cols, k % n_cols]
        info = invs[k][f"basis_vector_{k}"]
        grid = np.abs(corr_mod2_table[k]).reshape(BOARD, BOARD)
        im = ax.imshow(grid, origin="upper", cmap="Reds", vmin=0, vmax=1)
        for r in range(BOARD):
            for c in range(BOARD):
                ax.text(c, r, f"{grid[r, c]:.2f}", ha="center", va="center",
                        color="black" if grid[r, c] < 0.5 else "white",
                        fontsize=7)
        kind = info["type"].upper()
        max_r = info["max_r_winding_mod2"]
        title = (
            f"Inv {k} — GEOMÉTRICO  (max|r|={max_r:.2f})"
            if info["type"] == "geometric"
            else f"Inv {k} — ALGÉBRICO  (max|r|={max_r:.2f})"
        )
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("col")
        ax.set_ylabel("row")
    # Hide leftover axes
    for k in range(R, n_rows * n_cols):
        axes[k // n_cols, k % n_cols].axis("off")
    fig.colorbar(im, ax=axes.ravel().tolist(), fraction=0.025, pad=0.02)
    fig.suptitle(
        f"Correlação |r| entre bits de invariante H₁ e winding mod 2 "
        f"(threshold geométrico = {GEOMETRIC_THRESHOLD})",
        fontsize=12,
        y=1.00,
    )
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_CLASSIFICATION_PNG, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[T3] saved {OUT_CLASSIFICATION_PNG}")


# ---------------------------------------------------------------------------
# Task 2b — basis-independent classification
# ---------------------------------------------------------------------------
#
# The Pearson-r classification of Task 2 depends on the choice of basis for
# span(T): with my row-pivot basis the bits have no semantic alignment with
# the board geometry, so *every* bit ends up algebraic.  A more meaningful
# question — independent of basis — is:
#
#     Q1.  Is the integer winding W[:, p] an ℝ-linear function of T's
#          columns?  I.e. is there w_p ∈ ℝ^|E| with T @ w_p = W[:, p]
#          exactly?  Equivalently: W[:, p] ∈ col(T) over ℝ.
#
#     Q2.  Is W[:, p] mod 2 a GF(2)-linear function of T's columns?
#          I.e. is W[:, p] mod 2 ∈ col(T) over GF(2)?
#
# A "yes" to Q1 means that the winding number around cell p is a sum of
# fixed per-edge weights — a 1-cochain.  A "yes" to Q2 is the GF(2)
# analogue and is exactly what "geometric H₁ invariant" should mean.
#
# Subtlety: my edge incidence matrix ``T`` is *undirected*.  The actual
# discrete winding formula uses *directed* edges.  Linearity over an
# undirected basis would be a stronger property (it would mean: the
# winding's value doesn't depend on traversal direction).  We test both
# the undirected (T) and the directed (T_dir) settings.

def tours_to_directed_edge_matrix(tours: np.ndarray) -> tuple[np.ndarray, list[tuple[int, int]]]:
    """Return ``T_dir`` of shape ``(K, 2|E|)`` over GF(2) and the column key.

    ``T_dir[i, k] = 1`` iff tour ``i`` traverses the directed edge
    indexed by ``k`` (i.e. ``key[k] = (u, v)`` with ``u → v``).
    """
    K, N = tours.shape
    directed_keys: list[tuple[int, int]] = []
    key_to_idx: dict[tuple[int, int], int] = {}
    for u, v in EDGES:
        for a, b in ((u, v), (v, u)):
            key_to_idx[(a, b)] = len(directed_keys)
            directed_keys.append((a, b))
    T_dir = np.zeros((K, len(directed_keys)), dtype=np.uint8)
    for i in range(K):
        for t in range(N):
            a, b = int(tours[i, t]), int(tours[i, (t + 1) % N])
            T_dir[i, key_to_idx[(a, b)]] = 1
    return T_dir, directed_keys


def lstsq_residual(A: np.ndarray, y: np.ndarray) -> tuple[float, float, np.ndarray]:
    """Solve A b ≈ y over ℝ.  Return (max_abs_residual, r2, b)."""
    A_f = A.astype(np.float64)
    y_f = y.astype(np.float64)
    b, *_ = np.linalg.lstsq(A_f, y_f, rcond=None)
    pred = A_f @ b
    res = pred - y_f
    ss_res = float((res ** 2).sum())
    var = float(((y_f - y_f.mean()) ** 2).sum())
    r2 = 1.0 - (ss_res / var) if var > 0 else float("nan")
    return float(np.max(np.abs(res))), r2, b


def gf2_in_col_space(T: np.ndarray, y: np.ndarray, rank_T: int) -> bool:
    """Return True iff ``y ∈ col(T)`` over GF(2).

    Implementation: rank of ``[T | y]`` equals rank of ``T``.
    """
    augmented = np.concatenate([T, y.reshape(-1, 1)], axis=1)
    rank_aug, _, _ = gf2_rank_and_basis(augmented)
    return rank_aug == rank_T


def task2b_basis_free(verbose: bool = True) -> dict:
    """Test winding-number linearity in T (undirected) and T_dir (directed)."""
    tours = load_or_compute_tours(verbose=False)
    T = tours_to_edge_matrix(tours)
    windings = np.load(DATA_DIR / "winding_numbers.npy")
    K = T.shape[0]
    W = windings.reshape(K, BOARD * BOARD).astype(np.int32)
    W_mod2 = (W & 1).astype(np.uint8)
    rank_T = gf2_rank_and_basis(T)[0]
    if verbose:
        print(f"[T2b] rank(T) over GF(2) = {rank_T}")
        print(f"[T2b] building directed-edge matrix T_dir …")
    T_dir, _ = tours_to_directed_edge_matrix(tours)
    if verbose:
        print(f"[T2b] T_dir shape: {T_dir.shape}")
        rank_T_dir = gf2_rank_and_basis(T_dir)[0]
        print(f"[T2b] rank(T_dir) over GF(2) = {rank_T_dir}")

    cells = []
    n_R_exact = 0  # cells with integer winding exactly ℝ-linear in T
    n_R_exact_dir = 0
    n_GF2_in_T = 0
    n_GF2_in_T_dir = 0
    for p in range(BOARD * BOARD):
        r, c = divmod(p, BOARD)
        y_int = W[:, p]
        y_mod2 = W_mod2[:, p]
        # Trivial case: cells outside convex hull have y_int = 0
        if int(np.abs(y_int).sum()) == 0:
            cells.append({
                "cell_rc": [r, c],
                "trivial_zero": True,
                "int_R_linear_in_T": True,
                "int_R_linear_in_Tdir": True,
                "mod2_GF2_linear_in_T": True,
                "mod2_GF2_linear_in_Tdir": True,
                "R_max_residual_T": 0.0,
                "R2_T": float("nan"),
                "R_max_residual_Tdir": 0.0,
                "R2_Tdir": float("nan"),
            })
            n_R_exact += 1
            n_R_exact_dir += 1
            n_GF2_in_T += 1
            n_GF2_in_T_dir += 1
            continue
        # ℝ-linearity in undirected T
        max_res_T, r2_T, _ = lstsq_residual(T, y_int)
        is_R_lin_T = max_res_T < 1e-6
        # ℝ-linearity in directed T_dir
        max_res_Td, r2_Td, _ = lstsq_residual(T_dir, y_int)
        is_R_lin_Td = max_res_Td < 1e-6
        # GF(2)-linearity in T and T_dir
        is_GF2_T = gf2_in_col_space(T, y_mod2, rank_T)
        is_GF2_Td = gf2_in_col_space(T_dir, y_mod2,
                                     gf2_rank_and_basis(T_dir)[0]
                                     if not hasattr(task2b_basis_free, "_cached_rank_Tdir")
                                     else task2b_basis_free._cached_rank_Tdir)
        cells.append({
            "cell_rc": [r, c],
            "trivial_zero": False,
            "int_R_linear_in_T": bool(is_R_lin_T),
            "int_R_linear_in_Tdir": bool(is_R_lin_Td),
            "mod2_GF2_linear_in_T": bool(is_GF2_T),
            "mod2_GF2_linear_in_Tdir": bool(is_GF2_Td),
            "R_max_residual_T": max_res_T,
            "R2_T": r2_T,
            "R_max_residual_Tdir": max_res_Td,
            "R2_Tdir": r2_Td,
        })
        n_R_exact += int(is_R_lin_T)
        n_R_exact_dir += int(is_R_lin_Td)
        n_GF2_in_T += int(is_GF2_T)
        n_GF2_in_T_dir += int(is_GF2_Td)
        if verbose:
            print(f"[T2b] cell ({r},{c})  "
                  f"int-lin/T:{'Y' if is_R_lin_T else 'N'}  R²={r2_T:.4f}   "
                  f"int-lin/T_dir:{'Y' if is_R_lin_Td else 'N'}  R²={r2_Td:.4f}   "
                  f"mod2/GF2/T:{'Y' if is_GF2_T else 'N'}  "
                  f"mod2/GF2/T_dir:{'Y' if is_GF2_Td else 'N'}")

    # The decisive count: rank of the 36 winding-mod-2 vectors when viewed
    # as functionals on tours.  This is the number of LINEARLY-INDEPENDENT
    # geometric H₁ invariants extractable from windings.
    W_mod2_matrix = W_mod2  # shape (K, 36)
    rank_W_mod2 = gf2_rank_and_basis(W_mod2_matrix.T)[0]  # rank over GF(2)
    rank_W_int = int(np.linalg.matrix_rank(W.astype(np.float64)))
    if verbose:
        print(f"[T2b] rank of {{winding mod 2 per cell}} as GF(2) functionals: {rank_W_mod2}")
        print(f"[T2b] rank of {{integer winding per cell}}  over ℝ:           {rank_W_int}")

    summary = {
        "rank_T_gf2": int(rank_T),
        "rank_Tdir_gf2": int(gf2_rank_and_basis(T_dir)[0]),
        "n_cells_total": BOARD * BOARD,
        "n_cells_int_R_linear_in_T": int(n_R_exact),
        "n_cells_int_R_linear_in_Tdir": int(n_R_exact_dir),
        "n_cells_mod2_GF2_linear_in_T": int(n_GF2_in_T),
        "n_cells_mod2_GF2_linear_in_Tdir": int(n_GF2_in_T_dir),
        "rank_winding_mod2_functionals_gf2": int(rank_W_mod2),
        "rank_winding_int_functionals_real": int(rank_W_int),
        "cells": cells,
    }
    OUT_BASIS_FREE.write_text(json.dumps(summary, indent=2))
    if verbose:
        print(f"[T2b] saved {OUT_BASIS_FREE}")
        print()
        print(f"[T2b] === SUMMARY (basis-independent) ===")
        print(f"[T2b] cells with int winding ℝ-linear in T       (undirected): {n_R_exact} / 36")
        print(f"[T2b] cells with int winding ℝ-linear in T_dir   (directed):   {n_R_exact_dir} / 36")
        print(f"[T2b] cells with mod-2 winding GF(2)-linear in T (undirected): {n_GF2_in_T} / 36")
        print(f"[T2b] cells with mod-2 winding GF(2)-linear in T_dir (directed): {n_GF2_in_T_dir} / 36")
    return summary


def task3b_plot(summary: dict | None = None) -> None:
    """Plot the basis-independent classification grid."""
    if summary is None:
        summary = json.loads(OUT_BASIS_FREE.read_text())
    cells = summary["cells"]
    # Encode each cell as one of 4 colours (4 = perfect both, 3 = only directed,
    # 2 = only mod2/dir, 1 = trivial zero, 0 = nothing).
    grid_int = np.zeros((BOARD, BOARD), dtype=np.int8)
    grid_mod2 = np.zeros((BOARD, BOARD), dtype=np.int8)
    grid_r2 = np.zeros((BOARD, BOARD), dtype=np.float64)
    for cell in cells:
        r, c = cell["cell_rc"]
        if cell["trivial_zero"]:
            grid_int[r, c] = -1
            grid_mod2[r, c] = -1
            grid_r2[r, c] = np.nan
        else:
            grid_int[r, c] = (
                2 if cell["int_R_linear_in_T"]
                else 1 if cell["int_R_linear_in_Tdir"]
                else 0
            )
            grid_mod2[r, c] = (
                2 if cell["mod2_GF2_linear_in_T"]
                else 1 if cell["mod2_GF2_linear_in_Tdir"]
                else 0
            )
            grid_r2[r, c] = cell["R2_Tdir"]

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    labels_int = {-1: "trivial 0", 0: "alg", 1: "lin in T_dir", 2: "lin in T"}
    cmap_int = plt.cm.get_cmap("RdYlGn", 4)
    ax = axes[0]
    im = ax.imshow(grid_int, origin="upper", cmap=cmap_int, vmin=-1, vmax=2)
    for r in range(BOARD):
        for c in range(BOARD):
            ax.text(c, r, labels_int[int(grid_int[r, c])],
                    ha="center", va="center", fontsize=8,
                    color="black")
    ax.set_title("Integer winding — ℝ-linearity\nin edge incidence")
    ax.set_xlabel("col"); ax.set_ylabel("row")

    ax = axes[1]
    im = ax.imshow(grid_mod2, origin="upper", cmap=cmap_int, vmin=-1, vmax=2)
    for r in range(BOARD):
        for c in range(BOARD):
            ax.text(c, r, labels_int[int(grid_mod2[r, c])],
                    ha="center", va="center", fontsize=8,
                    color="black")
    ax.set_title("Mod-2 winding — GF(2)-linearity\nin edge incidence")
    ax.set_xlabel("col"); ax.set_ylabel("row")

    ax = axes[2]
    im = ax.imshow(grid_r2, origin="upper", cmap="viridis", vmin=0, vmax=1)
    for r in range(BOARD):
        for c in range(BOARD):
            txt = "—" if np.isnan(grid_r2[r, c]) else f"{grid_r2[r, c]:.3f}"
            ax.text(c, r, txt, ha="center", va="center", fontsize=8,
                    color="white" if (not np.isnan(grid_r2[r, c]) and grid_r2[r, c] < 0.6) else "black")
    ax.set_title("R² of OLS regression\n(integer winding on T_dir)")
    ax.set_xlabel("col"); ax.set_ylabel("row")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    fig.suptitle(
        "Basis-independent: is winding number linearly recoverable\n"
        "from the tour's edge incidence?",
        fontsize=12,
    )
    fig.tight_layout()
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_BASIS_FREE_PNG, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[T3b] saved {OUT_BASIS_FREE_PNG}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    tours = load_or_compute_tours(verbose=True)
    print(f"[main] tours loaded: {tours.shape}")
    print()
    print("=" * 72)
    print("TASK 1 — Reconstruir invariantes H₁ reais sobre GF(2)")
    print("=" * 72)
    t1 = task1_build_basis(tours, verbose=True)
    print()
    print("=" * 72)
    print("TASK 2 — Correlação H₁ ↔ winding")
    print("=" * 72)
    t2 = task2_classify_invariants(verbose=True)
    print()
    print("=" * 72)
    print("TASK 3 — Visualização")
    print("=" * 72)
    task3_plot(t2)
    print()
    print("=" * 72)
    print("TASK 2b — Classificação basis-independent (winding ∈ col(T)?)")
    print("=" * 72)
    t2b = task2b_basis_free(verbose=True)
    print()
    print("=" * 72)
    print("TASK 3b — Visualização basis-independent")
    print("=" * 72)
    task3b_plot(t2b)
    print()
    print("Done.")


if __name__ == "__main__":
    main()
