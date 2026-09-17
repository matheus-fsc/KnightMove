"""
Estudo de calibração de f∞(L).

Mede sensibilidade de nós/tour a perturbações em f∞(L) no 10×10:
  - 2.1 perturbação uniforme δ ∈ {-0.10, -0.05, 0, +0.05, +0.10}
  - 2.2 perturbação por nível L ∈ {0..4} com δ ∈ {-0.10, +0.10}
  - 2.3 baselines: H_uniforme (f∞=0.25), H_invertida, H_aleatoria

Resposta: quanto da eficiência vem da pressão de vértice vs da fase local?
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from theory_heuristic import (  # noqa: E402
    TheoryHeuristicV2, F_INF_DEFAULT, F_INF_DEFAULT_FALLBACK,
    construir_freq_teorica,
)
from theory_heuristic import build_knight_graph as _bkg  # noqa: F401


N_BENCH = 10
K_BENCH = 200
TIMEOUT_BENCH = 600.0


def perturb_uniforme(delta: float) -> dict[int, float]:
    return {L: float(np.clip(v + delta, 0.01, 0.99))
            for L, v in F_INF_DEFAULT.items()}


def perturb_nivel(L_alvo: int, delta: float) -> dict[int, float]:
    out = dict(F_INF_DEFAULT)
    out[L_alvo] = float(np.clip(out[L_alvo] + delta, 0.01, 0.99))
    return out


def rodar(n: int, K: int, timeout: float,
          f_inf: dict[int, float] | None,
          fallback: float = F_INF_DEFAULT_FALLBACK,
          custom_freq: np.ndarray | None = None,
          label: str = "") -> dict:
    """Roda backtracking teórico e retorna métricas resumidas.
    Se custom_freq é fornecido, usa-o em vez da tabela f_inf (útil para
    H_invertida / H_aleatoria que não são derivadas de f∞ por nível)."""
    bt = TheoryHeuristicV2(
        n=n, alvo=K, timeout=timeout, f_inf=f_inf, fallback=fallback,
    )
    if custom_freq is not None:
        assert custom_freq.shape == bt.freq.shape
        bt.freq = custom_freq.astype(np.float64)
    t0 = time.perf_counter()
    r = bt.executar()
    dt = time.perf_counter() - t0
    out = {
        "label": label,
        "n_tours": r["n_tours"], "n_2fatores": r["n_2fatores"],
        "n_nos": r["n_nos"], "nos_por_tour": r["nos_por_tour"],
        "razao_2fat_tour": r["razao_2fat_tour"],
        "tempo_s": r["tempo_s"], "tempo_real_s": dt,
        "parou_por_timeout": r.get("parou_por_timeout", False),
    }
    return out


def passo_21() -> list[dict]:
    print("\n--- 2.1 Perturbação uniforme em todos os níveis ---")
    out = []
    for delta in [-0.10, -0.05, 0.0, +0.05, +0.10]:
        f = perturb_uniforme(delta)
        label = f"uniform_delta={delta:+.2f}"
        r = rodar(N_BENCH, K_BENCH, TIMEOUT_BENCH, f_inf=f, label=label)
        r["delta"] = delta
        print(f"  δ={delta:+.2f}  nós/tour={r['nos_por_tour']:6.2f}  "
              f"t={r['tempo_s']:6.2f}s  razão={r['razao_2fat_tour']:.3f}×")
        out.append(r)
    return out


def passo_22() -> list[dict]:
    print("\n--- 2.2 Perturbação por nível (δ ∈ {-0.10, +0.10}) ---")
    out = []
    for L in [0, 1, 2, 3, 4]:
        for delta in [-0.10, +0.10]:
            f = perturb_nivel(L, delta)
            label = f"L={L}_delta={delta:+.2f}"
            r = rodar(N_BENCH, K_BENCH, TIMEOUT_BENCH, f_inf=f, label=label)
            r["L_alvo"] = L
            r["delta"] = delta
            r["f_inf_perturbado"] = f[L]
            print(f"  L={L}  δ={delta:+.2f}  nós/tour={r['nos_por_tour']:6.2f}  "
                  f"t={r['tempo_s']:6.2f}s")
            out.append(r)
    return out


def passo_23() -> list[dict]:
    print("\n--- 2.3 Baselines (H_uniforme, H_invertida, H_aleatoria) ---")
    out = []

    # H_teoria (referência)
    r = rodar(N_BENCH, K_BENCH, TIMEOUT_BENCH, f_inf=None, label="H_teoria")
    print(f"  H_teoria      → nós/tour={r['nos_por_tour']:6.2f}  "
          f"t={r['tempo_s']:6.2f}s")
    out.append(r)

    # H_uniforme: f∞(L) = 0.25 para todo L
    f_uniforme = {L: 0.25 for L in range(0, 10)}
    r = rodar(N_BENCH, K_BENCH, TIMEOUT_BENCH, f_inf=f_uniforme,
              label="H_uniforme")
    print(f"  H_uniforme    → nós/tour={r['nos_por_tour']:6.2f}  "
          f"t={r['tempo_s']:6.2f}s")
    out.append(r)

    # H_invertida: f∞ → 1 - f∞
    f_inv = {L: float(np.clip(1.0 - v, 0.01, 0.99))
             for L, v in F_INF_DEFAULT.items()}
    r = rodar(N_BENCH, K_BENCH, TIMEOUT_BENCH, f_inf=f_inv, label="H_invertida")
    print(f"  H_invertida   → nós/tour={r['nos_por_tour']:6.2f}  "
          f"t={r['tempo_s']:6.2f}s")
    out.append(r)

    # H_aleatoria: freq[e] ~ U(0,1) sorteado por aresta
    bt_tmp = TheoryHeuristicV2(n=N_BENCH, alvo=1, timeout=10.0)
    rng = np.random.default_rng(2024)
    freq_aleatoria = rng.uniform(0.01, 0.99, size=bt_tmp.E)
    del bt_tmp
    r = rodar(N_BENCH, K_BENCH, TIMEOUT_BENCH, f_inf=None,
              custom_freq=freq_aleatoria, label="H_aleatoria")
    print(f"  H_aleatoria   → nós/tour={r['nos_por_tour']:6.2f}  "
          f"t={r['tempo_s']:6.2f}s")
    out.append(r)

    return out


def plot_heatmap(res_22: list[dict], baselines: list[dict],
                 out_path: Path) -> None:
    """Heatmap: linhas = δ (-0.10, baseline, +0.10), colunas = L."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    Ls = sorted({r["L_alvo"] for r in res_22})
    deltas = [-0.10, 0.0, +0.10]
    Z = np.zeros((len(deltas), len(Ls)))

    baseline = next(b for b in baselines if b["label"] == "H_teoria")
    bl_nos = baseline["nos_por_tour"]

    # δ=0 linha (referência) usa o valor de H_teoria para todos L
    for j, L in enumerate(Ls):
        for i, d in enumerate(deltas):
            if d == 0.0:
                Z[i, j] = bl_nos
            else:
                cand = [r for r in res_22 if r["L_alvo"] == L and r["delta"] == d]
                if cand:
                    Z[i, j] = cand[0]["nos_por_tour"]
                else:
                    Z[i, j] = np.nan

    fig, ax = plt.subplots(figsize=(8, 4))
    im = ax.imshow(Z, aspect="auto", cmap="RdYlBu_r",
                   origin="lower",
                   vmin=min(Z.min(), bl_nos - 0.5),
                   vmax=max(Z.max(), bl_nos + 0.5))
    ax.set_xticks(range(len(Ls)))
    ax.set_xticklabels([f"L={L}" for L in Ls])
    ax.set_yticks(range(len(deltas)))
    ax.set_yticklabels([f"δ={d:+.2f}" for d in deltas])
    ax.set_xlabel("nível L perturbado (apenas)")
    ax.set_ylabel("δ aplicado a f∞(L)")
    ax.set_title(
        f"Sensibilidade de nós/tour a perturbações em f∞(L)  "
        f"[10×10, K={K_BENCH}]\n"
        f"baseline H_teoria: nós/tour = {bl_nos:.2f}"
    )
    for i in range(len(deltas)):
        for j in range(len(Ls)):
            val = Z[i, j]
            color = "white" if abs(val - bl_nos) > 1.5 else "black"
            ax.text(j, i, f"{val:.2f}", ha="center", va="center",
                    color=color, fontsize=10)
    plt.colorbar(im, label="nós/tour")
    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"plot salvo: {out_path}")


def main():
    print(f"=== Calibração — 10×10, K={K_BENCH} ===")

    res_21 = passo_21()
    res_22 = passo_22()
    res_23 = passo_23()

    # análise: qual L é mais sensível?
    max_sens = {}
    baseline = next(r for r in res_23 if r["label"] == "H_teoria")
    bl_nos = baseline["nos_por_tour"]
    for r in res_22:
        L = r["L_alvo"]
        diff = abs(r["nos_por_tour"] - bl_nos)
        max_sens[L] = max(max_sens.get(L, 0.0), diff)
    print(f"\n--- Sensibilidade por nível (|Δ nós/tour| vs baseline) ---")
    for L in sorted(max_sens):
        print(f"  L={L}: max |Δ| = {max_sens[L]:.2f}  "
              f"(f∞(L)={F_INF_DEFAULT.get(L, '-')})")
    L_max = max(max_sens, key=max_sens.get)
    print(f"\n  → Nível mais sensível: L={L_max} "
          f"(|Δ|={max_sens[L_max]:.2f})")

    out = {
        "config": {
            "n": N_BENCH, "K": K_BENCH, "timeout": TIMEOUT_BENCH,
            "f_inf_default": F_INF_DEFAULT,
        },
        "passo_21_uniforme": res_21,
        "passo_22_por_nivel": res_22,
        "passo_23_baselines": res_23,
        "sensibilidade_por_nivel": max_sens,
        "nivel_mais_sensivel": int(L_max),
    }
    out_path = ROOT / "data" / "calibration_results.json"
    out_path.write_text(json.dumps(out, indent=2))
    print(f"\nSalvo em {out_path}")

    plot_heatmap(res_22, res_23, ROOT / "data" / "plots" / "sensitivity_heatmap.png")


if __name__ == "__main__":
    main()
