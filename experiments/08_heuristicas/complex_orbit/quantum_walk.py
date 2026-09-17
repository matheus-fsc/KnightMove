"""quantum_walk.py — Task D.

Continuous-time quantum walk on the 6×6 knight graph.

The user's spec mentions "Grover coin OR adjacency normalised to be
unitary".  On an irregular graph the row-normalised adjacency is not
Hermitian and exp(...) of it is not the most natural object, so we
implement the standard **continuous-time quantum walk** (CTQW):

    |ψ(t+1)⟩ = U_τ |ψ(t)⟩,        U_τ = exp(−i · A · τ).

``A`` is the 0/1 adjacency of the 6×6 knight graph (Hermitian → ``U_τ``
unitary by spectral calculus).  We pick a step ``τ = 1``, since the
spectral radius of ``A`` on the 6×6 knight graph is finite (max degree
on a 6×6 board is 8, so ‖A‖₂ ≤ 8 and ``τ = 1`` gives non-trivial
rotation without making each step a near-2π loop).

We start from the corner cell A6 (vertex 0) and evolve for 100 steps.
At step ``t = N = 36`` the spec asks for a comparison to the uniform
distribution over the 36 cells (the "endpoint" of a classical closed
tour).  We report total variation distance from uniform.

Outputs:

* ``data/quantum_walk_probabilities.npy``     shape ``(101, 36)``  float64
* ``data/quantum_walk_phases.npy``            shape ``(101, 36)``  float64 in (−π, π]
* ``data/plots/quantum_walk_snapshots.png``   4 snapshots t=0, 10, 36, 100
* ``data/plots/quantum_walk_uniform_distance.png``   ‖ψ(t)|² − uniform‖_TV vs t
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from _tours import ADJ, BOARD, TOTAL, vertex_label

DATA_DIR = Path(__file__).parent / "data"
PLOT_DIR = DATA_DIR / "plots"

OUT_PROB = DATA_DIR / "quantum_walk_probabilities.npy"
OUT_PHASE = DATA_DIR / "quantum_walk_phases.npy"
OUT_SNAPSHOTS = PLOT_DIR / "quantum_walk_snapshots.png"
OUT_TVDIST = PLOT_DIR / "quantum_walk_uniform_distance.png"
OUT_SUMMARY = DATA_DIR / "quantum_walk_summary.json"

TAU = 1.0
N_STEPS = 100
START_VERTEX = 0  # A6 = (row 0, col 0)


def build_adjacency() -> np.ndarray:
    """Symmetric 0/1 adjacency of the 6×6 knight graph."""
    A = np.zeros((TOTAL, TOTAL), dtype=np.float64)
    for v in range(TOTAL):
        for u in ADJ[v]:
            A[v, u] = 1.0
    assert np.allclose(A, A.T), "knight graph adjacency must be symmetric"
    return A


def evolve(A: np.ndarray, psi0: np.ndarray, n_steps: int, tau: float) -> np.ndarray:
    """Return ``ψ(0..n_steps)`` of shape ``(n_steps+1, TOTAL)`` complex.

    Uses one diagonalisation of ``A`` (Hermitian) followed by step-wise
    multiplication in the eigenbasis — much cheaper than re-evaluating
    ``expm`` at every step.
    """
    eigvals, eigvecs = np.linalg.eigh(A)  # A = V diag(λ) V†
    coeffs = eigvecs.conj().T @ psi0  # length TOTAL
    out = np.empty((n_steps + 1, TOTAL), dtype=np.complex128)
    out[0] = psi0
    for t in range(1, n_steps + 1):
        phase = np.exp(-1j * eigvals * t * tau)
        out[t] = eigvecs @ (phase * coeffs)
    return out


def total_variation_to_uniform(prob: np.ndarray) -> float:
    """½ Σ |p_i − 1/N|."""
    return float(0.5 * np.sum(np.abs(prob - 1.0 / TOTAL)))


def plot_snapshots(probs: np.ndarray, times: list[int], path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, len(times), figsize=(3.6 * len(times), 4.0))
    if len(times) == 1:
        axes = [axes]
    vmax = max(probs[t].max() for t in times)
    for ax, t in zip(axes, times):
        grid = probs[t].reshape(BOARD, BOARD)
        im = ax.imshow(grid, origin="upper", cmap="viridis", vmin=0, vmax=vmax)
        for r in range(BOARD):
            for c in range(BOARD):
                ax.text(c, r, f"{grid[r, c]:.2f}", ha="center", va="center",
                        color="white" if grid[r, c] < vmax * 0.6 else "black",
                        fontsize=8)
        ax.set_title(f"t = {t}    TV(p, uniform) = {total_variation_to_uniform(probs[t]):.3f}")
        ax.set_xlabel("col")
        ax.set_ylabel("row")
    fig.colorbar(im, ax=axes, fraction=0.025, pad=0.04)
    fig.suptitle(
        f"CTQW on the 6×6 knight graph, start = {vertex_label(START_VERTEX)}, τ = {TAU}",
        y=1.02,
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def plot_tv_distance(probs: np.ndarray, path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    tv = np.array([total_variation_to_uniform(probs[t]) for t in range(probs.shape[0])])
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(tv, lw=1.6, color="#c44")
    ax.axvline(36, ls=":", color="grey", label="t = N = 36 (tour length)")
    ax.set_xlabel("step t")
    ax.set_ylabel("TV distance to uniform distribution")
    ax.set_title("CTQW probability mass vs. uniform on 6×6 knight graph")
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def main() -> dict:
    A = build_adjacency()
    spectral_radius = float(np.max(np.abs(np.linalg.eigvalsh(A))))
    print(f"[D] 6×6 knight graph: |V|={TOTAL}, spectral radius ‖A‖₂ = {spectral_radius:.4f}")
    psi0 = np.zeros(TOTAL, dtype=np.complex128)
    psi0[START_VERTEX] = 1.0
    print(f"[D] starting from vertex {START_VERTEX} ({vertex_label(START_VERTEX)})")
    print(f"[D] evolving for {N_STEPS} steps with τ = {TAU} …")
    psi = evolve(A, psi0, N_STEPS, TAU)
    probs = np.abs(psi) ** 2
    phases = np.angle(psi)

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    np.save(OUT_PROB, probs)
    np.save(OUT_PHASE, phases)
    print(f"[D] saved {OUT_PROB}")
    print(f"[D] saved {OUT_PHASE}")

    # Snapshots
    snapshot_times = [0, 10, 36, N_STEPS]
    plot_snapshots(probs, snapshot_times, OUT_SNAPSHOTS)
    print(f"[D] saved {OUT_SNAPSHOTS}")

    plot_tv_distance(probs, OUT_TVDIST)
    print(f"[D] saved {OUT_TVDIST}")

    # Numerical summary
    tv36 = total_variation_to_uniform(probs[36])
    tv100 = total_variation_to_uniform(probs[100])
    norm_err = float(np.abs(probs.sum(axis=1) - 1.0).max())
    info = {
        "n_vertices": TOTAL,
        "start_vertex_label": vertex_label(START_VERTEX),
        "spectral_radius": spectral_radius,
        "tau": TAU,
        "n_steps": N_STEPS,
        "unitarity_max_norm_error": norm_err,
        "tv_to_uniform_at_t36": tv36,
        "tv_to_uniform_at_t100": tv100,
        "min_tv_to_uniform": float(min(total_variation_to_uniform(probs[t]) for t in range(N_STEPS + 1))),
        "argmin_tv_step": int(np.argmin([total_variation_to_uniform(probs[t]) for t in range(N_STEPS + 1)])),
        "max_prob_at_t36": float(probs[36].max()),
        "argmax_prob_at_t36_label": vertex_label(int(np.argmax(probs[36]))),
    }
    OUT_SUMMARY.write_text(json.dumps(info, indent=2))
    print(f"[D] saved {OUT_SUMMARY}")
    print("\n=== Task D — quantum walk summary ===")
    for k, v in info.items():
        print(f"  {k}: {v}")
    return info


if __name__ == "__main__":
    main()
