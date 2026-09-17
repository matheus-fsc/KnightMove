"""
Analise espectral de T_bulk(n): autovalores dominantes e autovetor.

Entrada: data/transfer_n{n}_bulk.npz e data/transfer_n{n}_idx.json
Saida:
  - data/eigenvalues_n{n}.json
  - data/eigenvector_dominant_n{n}.npz
"""

import json
import sys
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

DATA = Path(__file__).parent / "data"


def load(n):
    T = sp.load_npz(DATA / f"transfer_n{n}_bulk.npz").astype(np.float64)
    with open(DATA / f"transfer_n{n}_idx.json") as f:
        idx_raw = json.load(f)
    # idx_raw chave eh "deg_c_tuple|deg_cp1_tuple" string
    idx = {}
    for k, v in idx_raw.items():
        a, b = k.split("|")
        a_t = tuple(int(x) for x in a.strip("()").split(",") if x.strip())
        b_t = tuple(int(x) for x in b.strip("()").split(",") if x.strip())
        idx[(a_t, b_t)] = v
    return T, idx


def top_eigs(T, k=10, side="right"):
    """Retorna (eigvals, eigvecs) ordenados por |lambda| decrescente.

    side="right": T @ v = lam * v
    side="left":  T.T @ w = lam * w  (autovetores a esquerda de T)
    """
    K = T.shape[0]
    k = min(k, K - 2)

    M = T.astype(np.float64) if side == "right" else T.T.astype(np.float64)
    t0 = time.time()
    vals, vecs = spla.eigs(M, k=k, which="LM", maxiter=5000, tol=1e-9)
    t1 = time.time()
    order = np.argsort(-np.abs(vals))
    vals = vals[order]
    vecs = vecs[:, order]
    return vals, vecs, t1 - t0


def main(ns):
    summary = []
    for n in ns:
        print(f"\n=== n = {n} ===")
        T, idx = load(n)
        K = T.shape[0]
        print(f"  K = {K}, nnz = {T.nnz}")

        # Autovalores e autovetores a direita
        vals, vecs, dt = top_eigs(T, k=10, side="right")
        print(f"  eigs (direita) em {dt:.2f}s")

        # Autovetor a esquerda (T.T @ w = lam * w)
        vals_L, vecs_L, dt_L = top_eigs(T, k=10, side="left")
        print(f"  eigs (esquerda) em {dt_L:.2f}s")

        print(f"  10 autovalores top (|lambda| decrescente):")
        for i, v in enumerate(vals):
            re = v.real
            im = v.imag
            mag = abs(v)
            tag = ""
            if abs(im) < 1e-8:
                tag = " (real)"
            print(f"    [{i}] lambda = {re:+.6f} {im:+.6f}i  |lambda|={mag:.6f}{tag}")

        lam1 = vals[0]
        lam2 = vals[1]
        gap = abs(lam1) - abs(lam2)
        ratio = abs(lam2) / abs(lam1) if abs(lam1) > 0 else None
        print(f"  gap = |lambda1|-|lambda2| = {gap:.6f}")
        print(f"  ratio = |lambda2|/|lambda1| = {ratio:.6f}")

        # Salvar
        with open(DATA / f"eigenvalues_n{n}.json", "w") as f:
            json.dump(
                {
                    "n": n,
                    "K": K,
                    "nnz": int(T.nnz),
                    "eigenvalues_real": [float(v.real) for v in vals],
                    "eigenvalues_imag": [float(v.imag) for v in vals],
                    "eigenvalues_abs": [float(abs(v)) for v in vals],
                    "gap_abs": float(gap),
                    "ratio_abs": float(ratio) if ratio else None,
                    "time_eigs_s": dt,
                },
                f,
                indent=2,
            )

        # Salvar autovetores dominantes direito e esquerdo
        v1 = vecs[:, 0]
        w1 = vecs_L[:, 0]
        np.savez(
            DATA / f"eigenvector_dominant_n{n}.npz",
            v_real=v1.real,
            v_imag=v1.imag,
            w_real=w1.real,
            w_imag=w1.imag,
            lambda1_real=float(lam1.real),
            lambda1_imag=float(lam1.imag),
        )
        print(f"  autovetores salvos em data/eigenvector_dominant_n{n}.npz")

        summary.append({
            "n": n,
            "K": K,
            "lambda1": complex(lam1),
            "lambda2": complex(lam2),
            "gap": float(gap),
            "ratio": float(ratio) if ratio else None,
        })

    return summary


if __name__ == "__main__":
    ns = [int(x) for x in sys.argv[1:]] if len(sys.argv) > 1 else [6]
    main(ns)
