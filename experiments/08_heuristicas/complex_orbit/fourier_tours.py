"""fourier_tours.py — Task B.

DFT decomposition of every closed knight tour on the 6×6 board.

For each of the 9,862 canonical closed tours (one per undirected
Hamiltonian cycle), the position sequence

    z(t) = col(t) + i·row(t),   t = 0, 1, …, 35

is fed to the unitary-normalised Discrete Fourier Transform

    c_k = (1/N) · Σ_t  z(t) · exp(−2π i · k · t / N),     k = 0, …, N−1.

Outputs:

* ``data/fourier_coefficients.npy``    shape ``(9862, 36, 2)``   float64
  (last axis: ``[Re c_k, Im c_k]``)
* ``data/plots/fourier_mean_spectrum.png``    mean ``|c_k|`` ± 1σ vs k

A short stdout report summarises the spectral profile.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from _tours import N, load_or_compute_tours, tours_to_complex

DATA_DIR = Path(__file__).parent / "data"
PLOT_DIR = DATA_DIR / "plots"

OUT_NPY = DATA_DIR / "fourier_coefficients.npy"
OUT_PLOT = PLOT_DIR / "fourier_mean_spectrum.png"


def compute_dft(z_tours: np.ndarray) -> np.ndarray:
    """Compute c_k for every tour.

    Parameters
    ----------
    z_tours
        Complex array of shape ``(K, N)`` with ``K`` tours, ``N`` steps each.

    Returns
    -------
    coeffs
        Complex array of shape ``(K, N)`` with ``coeffs[t, k] = c_k`` of
        tour ``t``.  Uses the same normalisation as in the docstring of
        this module (i.e. the result of :func:`numpy.fft.fft` divided by
        ``N``).
    """
    return np.fft.fft(z_tours, axis=1) / z_tours.shape[1]


def coefficients_as_real_pairs(coeffs: np.ndarray) -> np.ndarray:
    """Convert ``(K, N)`` complex → ``(K, N, 2)`` float [Re, Im]."""
    out = np.empty(coeffs.shape + (2,), dtype=np.float64)
    out[..., 0] = coeffs.real
    out[..., 1] = coeffs.imag
    return out


def plot_mean_spectrum(coeffs: np.ndarray, path: Path) -> None:
    """Plot the cross-tour mean of ``|c_k|`` with ±1σ shading."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    mag = np.abs(coeffs)  # (K, N)
    mean = mag.mean(axis=0)
    std = mag.std(axis=0)
    k = np.arange(N)

    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.fill_between(k, mean - std, mean + std, alpha=0.25,
                    label="±1σ across tours", color="#3a8")
    ax.plot(k, mean, color="#085", lw=1.6, label="mean |c_k|")
    ax.axvline(N / 2, ls=":", color="grey", alpha=0.6,
               label=f"k = N/2 = {N // 2}")
    ax.set_xlabel("DFT index k")
    ax.set_ylabel("|c_k|")
    ax.set_title(
        f"DFT spectrum of 6×6 closed knight tours\n"
        f"({coeffs.shape[0]:,} tours, N = {N})"
    )
    ax.set_xticks(np.arange(0, N + 1, 4))
    ax.grid(alpha=0.25)
    ax.legend(loc="upper right")
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140)
    plt.close(fig)


def report(coeffs: np.ndarray) -> dict:
    """Numerical summary printed to stdout and returned for the run-log."""
    mag = np.abs(coeffs)
    mean = mag.mean(axis=0)
    std = mag.std(axis=0)

    # Sanity: c_0 is the centroid; for a Hamiltonian tour on a 6×6 board
    # the centroid must be the board centroid (2.5, 2.5) → 2.5 + 2.5i.
    expected_centroid = 2.5 + 2.5j
    c0_err = np.max(np.abs(coeffs[:, 0] - expected_centroid))

    info = {
        "n_tours": int(coeffs.shape[0]),
        "n_steps": int(coeffs.shape[1]),
        "centroid_max_error": float(c0_err),
        "mean_abs_c1": float(mean[1]),
        "mean_abs_c2": float(mean[2]),
        "top3_k_by_mean_mag": [int(x) for x in np.argsort(mean)[::-1][:5]],
        "spectral_floor_max_k": int(np.argmin(mean[1:]) + 1),
    }
    print("\n=== Task B — DFT spectrum summary ===")
    print(f"  tours: {info['n_tours']:,}    N: {info['n_steps']}")
    print(f"  centroid check |c_0 − 2.5(1+i)|_max = {info['centroid_max_error']:.3e}")
    print(f"  k with largest mean |c_k| (top-5): {info['top3_k_by_mean_mag']}")
    print(f"  smallest mean |c_k| outside k=0 is at k = {info['spectral_floor_max_k']}")
    print(f"  mean |c_1| = {info['mean_abs_c1']:.4f}  "
          f"(fundamental rotation amplitude)")
    return info


def main() -> dict:
    tours = load_or_compute_tours(verbose=True)
    print(f"[B] loaded tours: shape={tours.shape}")
    z = tours_to_complex(tours)  # (K, N) complex
    print(f"[B] complex sequences: shape={z.shape}")
    coeffs = compute_dft(z)
    print(f"[B] DFT done — coeffs shape={coeffs.shape}")
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    out = coefficients_as_real_pairs(coeffs)
    np.save(OUT_NPY, out)
    print(f"[B] saved {OUT_NPY}  (shape={out.shape}, dtype={out.dtype})")
    plot_mean_spectrum(coeffs, OUT_PLOT)
    print(f"[B] saved {OUT_PLOT}")
    return report(coeffs)


if __name__ == "__main__":
    main()
