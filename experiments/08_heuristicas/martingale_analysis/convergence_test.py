"""
T3 — testes formais de convergência (versão refinada pós-T2).

T2 mostrou que H_MART forte (p(b,n) constante em n) está REFUTADA — o processo
não é i.i.d. exato. Mas O(1) ainda pode valer se:

  (1) μ(n) não DECAI com n     — não rejeita inclinação > 0 no limite
  (2) Var[M_t] ≈ t·σ²          — variância de incrementos uniformemente limitada
                                  (TCL para martingales com dependência fraca)
  (3) μ_min(n) ≥ μ_min* > 0    — bound inferior empírico

Salvas em data/convergence_tests.json.
Plot: data/plots/convergence_evidence.png e martingale_process.png
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
PLOTS = DATA / "plots"
PLOTS.mkdir(parents=True, exist_ok=True)

NS = [6, 8, 10, 12, 14]


# ----------------------------------------------------------------------
# Carregar
# ----------------------------------------------------------------------
def load_X_seq(n: int) -> np.ndarray:
    ev = json.loads((DATA / f"event_log_n{n}.json").read_text())
    return np.array([1 if e["outcome"] == "TOUR" else 0 for e in ev],
                    dtype=np.int32)


# ----------------------------------------------------------------------
# Wilson lower / upper bound
# ----------------------------------------------------------------------
def wilson_lower(k: int, n: int, z: float = 1.96) -> float:
    if n == 0:
        return 0.0
    p = k / n
    denom = 1 + z * z / n
    center = (p + z * z / (2 * n)) / denom
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return max(0.0, center - half)


def wilson_upper(k: int, n: int, z: float = 1.96) -> float:
    if n == 0:
        return 1.0
    p = k / n
    denom = 1 + z * z / n
    center = (p + z * z / (2 * n)) / denom
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return min(1.0, center + half)


# ----------------------------------------------------------------------
# Teste 1 — Regressão μ(n) = a·n + b
# ----------------------------------------------------------------------
def test_no_decay(mus: dict[int, float], counts: dict[int, int]) -> dict:
    """H0: a = 0 em μ(n) = a·n + b. Usa regressão linear com pesos = √n_amostras."""
    ns = sorted(mus.keys())
    x = np.array(ns, dtype=float)
    y = np.array([mus[n] for n in ns])
    # ponderação por √N (variância ∝ 1/N)
    w = np.array([math.sqrt(counts[n]) for n in ns])

    # WLS via numpy
    X = np.vstack([x, np.ones_like(x)]).T
    W = np.diag(w * w)
    XtWX = X.T @ W @ X
    XtWy = X.T @ W @ y
    beta = np.linalg.solve(XtWX, XtWy)
    a, b = float(beta[0]), float(beta[1])

    yhat = X @ beta
    resid = y - yhat
    # SE de a via scipy.stats.linregress (não-ponderado, mas N pequeno)
    lin = stats.linregress(x, y)
    return {
        "slope_wls": a,
        "intercept_wls": b,
        "slope_ols": float(lin.slope),
        "intercept_ols": float(lin.intercept),
        "stderr_ols": float(lin.stderr),
        "pvalue_slope_eq_0": float(lin.pvalue),
        "r_squared": float(lin.rvalue ** 2),
        "ns": ns,
        "mus": [mus[n] for n in ns],
    }


# ----------------------------------------------------------------------
# Teste 2 — Var[M_t] linear em t
# ----------------------------------------------------------------------
def test_variance_linearity(X: np.ndarray, mu: float,
                            n_block_sizes: int = 10) -> dict:
    """Estima Var[M_t] usando BLOCOS NÃO-SOBREPOSTOS.

    Para vários tamanhos t, divide a sequência em K = floor(N/t) blocos
    disjuntos. Cada bloco gera S_i = sum(X_block - μ). Var[M_t] ≈ var(S_1..S_K).

    Para um martingale i.i.d.: Var[M_t] = t·μ(1-μ) (linear).
    Se Var ∝ t^α com α<1: incrementos negativamente correlacionados
    (concentração sub-linear, mais forte que martingale)."""
    N = len(X)
    incr = X.astype(float) - mu
    # tamanhos de bloco geometricamente espaçados — garantir K ≥ 8 blocos
    t_min = max(10, N // 200)
    t_max = N // 8
    ts = np.unique(np.geomspace(t_min, t_max, n_block_sizes).astype(int))
    var_mt = []
    n_blocks_per_t = []
    for t in ts:
        K = N // int(t)
        if K < 4:
            continue
        S = np.array([incr[i*int(t):(i+1)*int(t)].sum() for i in range(K)])
        var_mt.append(float(S.var(ddof=1)))
        n_blocks_per_t.append(int(K))

    # Ajuste 1: linear Var[M_t] = σ² · t (intercepto 0)
    ts = np.array(ts[:len(var_mt)], dtype=float)
    var_mt = np.array(var_mt)
    if len(ts) < 3:
        return {
            "N": N, "mu": mu, "sigma2_estimated": float("nan"),
            "sigma2_iid": float(mu*(1-mu)), "ratio_to_iid": float("nan"),
            "r2_linear": float("nan"), "power_alpha": float("nan"),
            "power_c": float("nan"), "r2_power": float("nan"),
            "ts": [], "var_mt": [], "n_blocks_per_t": [],
            "note": "blocos insuficientes",
        }
    sigma2 = float((ts * var_mt).sum() / (ts * ts).sum())
    ss_res_lin = float(((var_mt - sigma2 * ts) ** 2).sum())
    ss_tot = float(((var_mt - var_mt.mean()) ** 2).sum())
    r2_lin = 1 - ss_res_lin / max(ss_tot, 1e-12)

    # Ajuste 2: lei de potência Var = c · t^α (regressão em log)
    log_t = np.log(ts)
    log_v = np.log(np.clip(var_mt, 1e-12, None))
    lin_pow = stats.linregress(log_t, log_v)
    alpha = float(lin_pow.slope)
    log_c = float(lin_pow.intercept)
    r2_pow = float(lin_pow.rvalue ** 2)

    # comparação com o predito i.i.d.: σ²_iid = μ(1-μ)
    sigma2_iid = mu * (1 - mu)
    return {
        "N": N,
        "mu": mu,
        "sigma2_estimated": sigma2,
        "sigma2_iid": float(sigma2_iid),
        "ratio_to_iid": sigma2 / max(sigma2_iid, 1e-12),
        "r2_linear": r2_lin,
        "power_alpha": alpha,
        "power_c": math.exp(log_c),
        "r2_power": r2_pow,
        "ts": ts.tolist(),
        "var_mt": var_mt.tolist(),
        "n_blocks_per_t": n_blocks_per_t,
    }


# ----------------------------------------------------------------------
# Teste 3 — Bound inferior empírico de μ_min
# ----------------------------------------------------------------------
def test_mu_lower_bound(mus: dict[int, float],
                       counts_tour: dict[int, int],
                       counts_tot: dict[int, int]) -> dict:
    """IC de Wilson 95% individual para cada μ(n); LB = min sobre n."""
    rows = []
    for n in sorted(mus.keys()):
        lo = wilson_lower(counts_tour[n], counts_tot[n])
        hi = wilson_upper(counts_tour[n], counts_tot[n])
        rows.append({"n": n, "mu": mus[n],
                     "tours": counts_tour[n],
                     "events": counts_tot[n],
                     "wilson_lo": lo, "wilson_hi": hi})
    lb_min = min(r["wilson_lo"] for r in rows)
    return {
        "rows": rows,
        "min_lower_bound_95pct": lb_min,
        "implied_E_nodes_per_tour_upper": 1.0 / max(lb_min, 1e-12),
    }


# ----------------------------------------------------------------------
# Teste 4 — TCL para martingale (M_T/√T normal?)
# ----------------------------------------------------------------------
def test_martingale_clt(X: np.ndarray, mu: float,
                        block_size: int = 50) -> dict:
    """Divide a sequência em blocos de tamanho `block_size`. Para cada bloco,
    calcula S = sum(X_block - mu). Sob TCL, S/√block_size deve ser N(0, σ²)."""
    incr = X.astype(float) - mu
    n_blk = len(X) // block_size
    if n_blk < 5:
        return {"ok": False, "reason": "blocos insuficientes",
                "n_blocks": n_blk, "block_size": block_size}
    S = np.array([incr[i * block_size:(i + 1) * block_size].sum()
                  for i in range(n_blk)])
    S_norm = S / math.sqrt(block_size)
    # Shapiro-Wilk
    try:
        w_stat, w_p = stats.shapiro(S_norm)
    except Exception as exc:
        w_stat, w_p = float("nan"), float("nan")
    # KS contra N(0, σ²) estimado
    sigma = float(S_norm.std(ddof=1))
    ks_stat, ks_p = stats.kstest(
        (S_norm - S_norm.mean()) / max(sigma, 1e-9),
        "norm",
    )
    return {
        "ok": True,
        "n_blocks": int(n_blk),
        "block_size": int(block_size),
        "S_norm_mean": float(S_norm.mean()),
        "S_norm_sigma": sigma,
        "shapiro_stat": float(w_stat),
        "shapiro_pvalue": float(w_p),
        "ks_stat": float(ks_stat),
        "ks_pvalue": float(ks_p),
        "samples": S_norm.tolist(),
    }


# ----------------------------------------------------------------------
# Plot principal: convergence_evidence
# ----------------------------------------------------------------------
def plot_convergence_evidence(decay: dict, var_lin: dict, mu_lb: dict,
                              clt: dict, out_path: Path) -> None:
    fig = plt.figure(figsize=(15, 9))
    gs = fig.add_gridspec(2, 3, hspace=0.4, wspace=0.3)

    # painel 1: regressão μ(n) vs n
    ax = fig.add_subplot(gs[0, 0])
    ns = decay["ns"]
    mus = decay["mus"]
    ax.scatter(ns, mus, s=80, color="black", zorder=3)
    x_line = np.linspace(min(ns) - 1, max(ns) + 1, 50)
    y_line = decay["slope_ols"] * x_line + decay["intercept_ols"]
    ax.plot(x_line, y_line, "--", color="red", lw=2,
            label=f"μ ≈ {decay['slope_ols']:+.4f}·n + {decay['intercept_ols']:.3f}")
    pv = decay["pvalue_slope_eq_0"]
    tag = "NÃO rejeita H0 (sem tendência)" if pv > 0.05 else "REJEITA (há tendência)"
    ax.set_title(f"Teste 1: regressão μ(n)\np(slope=0) = {pv:.3f} — {tag}",
                 fontsize=10)
    ax.set_xlabel("n")
    ax.set_ylabel("μ(n)")
    ax.legend(fontsize=8, loc="lower right")
    ax.grid(alpha=0.3)
    ax.set_ylim(0, 0.4)

    # painel 2: μ_min com IC de Wilson
    ax = fig.add_subplot(gs[0, 1])
    rows = mu_lb["rows"]
    ns_lb = [r["n"] for r in rows]
    mus_lb = [r["mu"] for r in rows]
    los = [r["wilson_lo"] for r in rows]
    his = [r["wilson_hi"] for r in rows]
    ax.errorbar(ns_lb, mus_lb,
                yerr=[np.array(mus_lb) - np.array(los),
                      np.array(his) - np.array(mus_lb)],
                fmt="o", capsize=4, color="black", lw=1.5)
    lb_min = mu_lb["min_lower_bound_95pct"]
    ax.axhline(lb_min, color="green", ls="--", lw=2,
               label=f"min LB 95% = {lb_min:.3f}")
    ax.axhline(0.10, color="orange", ls=":", lw=1,
               label="threshold 0.10")
    ax.set_title(
        f"Teste 3: μ_min ≥ {lb_min:.3f} (95%)\n"
        f"⇒ E[nós/tour] ≤ {mu_lb['implied_E_nodes_per_tour_upper']:.2f}",
        fontsize=10,
    )
    ax.set_xlabel("n")
    ax.set_ylabel("μ(n) com IC 95%")
    ax.legend(fontsize=8, loc="lower right")
    ax.grid(alpha=0.3)
    ax.set_ylim(0, 0.4)

    # painel 3: tabela de p-values dos testes principais
    ax = fig.add_subplot(gs[0, 2])
    ax.axis("off")
    rows_tbl = [
        ["Teste 1: μ(n) sem tendência (p)", decay["pvalue_slope_eq_0"], "p"],
        ["Teste 2a: Var ~ t^α  (α médio)", var_lin["mean_alpha"], "alpha"],
        ["Teste 2b: R² do ajuste t^α", var_lin["mean_r2_power"], "r2"],
        ["Teste 4 (TCL): Shapiro avg", clt["mean_shapiro_p"], "p"],
        ["Teste 4 (TCL): KS avg",      clt["mean_ks_p"], "p"],
    ]
    cell_text = []
    cell_colors = []
    for label, val, kind in rows_tbl:
        txt = f"{val:.3f}"
        if kind == "p":
            color = "#a8e6a3" if val > 0.05 else "#f4a3a3"
        elif kind == "alpha":
            color = "#a8e6a3" if 0.7 < val < 1.3 else "#f4a3a3"
        elif kind == "r2":
            color = "#a8e6a3" if val > 0.9 else ("#ffe28a" if val > 0.7 else "#f4a3a3")
        else:
            color = "white"
        cell_text.append([label, txt])
        cell_colors.append(["white", color])
    table = ax.table(cellText=cell_text, cellColours=cell_colors,
                     cellLoc="center", colLabels=["Teste", "estatística"],
                     loc="center")
    table.scale(1.0, 1.6)
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    ax.set_title("Resumo dos testes formais", fontsize=10)

    # painel 4-5: Var[M_t] vs t para cada n
    for idx, n in enumerate(NS):
        if idx >= 3:
            break
        ax = fig.add_subplot(gs[1, idx])
        block = var_lin["per_n"][n]
        ts = np.array(block["ts"])
        vm = np.array(block["var_mt"])
        ax.scatter(ts, vm, s=20, color="#1f77b4", zorder=3)
        s2 = block["sigma2_estimated"]
        ax.plot(ts, s2 * ts, "--", color="red", lw=1.5,
                label=f"σ²·t  (σ²={s2:.4f})")
        s2_iid = block["sigma2_iid"]
        ax.plot(ts, s2_iid * ts, ":", color="gray", lw=1,
                label=f"σ²_iid·t ({s2_iid:.4f})")
        ax.set_title(f"n={n}  Var[M_t] vs t  R²={block['r2_linear']:.3f}",
                     fontsize=10)
        ax.set_xlabel("t (eventos)")
        ax.set_ylabel("Var[M_t]")
        ax.legend(fontsize=7)
        ax.grid(alpha=0.3)

    fig.suptitle(
        "T3 — convergência: μ não decai, Var[M_t] linear, μ_min com IC 95%",
        fontsize=12,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out_path, dpi=130)
    plt.close(fig)


# ----------------------------------------------------------------------
# Plot martingale_process.png — n=10 detalhado
# ----------------------------------------------------------------------
def plot_martingale_process(X: np.ndarray, mu: float, n_label: int,
                            clt_block: dict, out_path: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(13, 8))

    # (a) X_t binário
    ax = axes[0, 0]
    idx = np.arange(len(X))
    ax.scatter(idx, X, s=4, color="black", alpha=0.5)
    ax.axhline(mu, color="red", ls="--", lw=1, label=f"μ={mu:.3f}")
    ax.set_xlabel("t")
    ax.set_ylabel("X_t")
    ax.set_title(f"(a) X_t para n={n_label}")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    # (b) M_t = Σ(X - μ)
    ax = axes[0, 1]
    Mt = np.cumsum(X.astype(float) - mu)
    ax.plot(idx, Mt, color="#1f77b4", lw=1)
    ax.axhline(0, color="black", ls="--", lw=0.8)
    sigma2_iid = mu * (1 - mu)
    bound = 2 * np.sqrt(sigma2_iid * idx)
    ax.fill_between(idx, -bound, bound, color="gray", alpha=0.2,
                    label="±2σ_iid·√t")
    ax.set_xlabel("t")
    ax.set_ylabel("M_t")
    ax.set_title("(b) M_t = Σ(X_k − μ)  (deve oscilar em torno de 0)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    # (c) Var[M_t] vs t
    ax = axes[1, 0]
    # bootstrap pequeno: pega blocos de tamanho t
    N = len(X)
    rng = np.random.default_rng(0)
    ts = np.linspace(N // 20, N - 1, 20).astype(int)
    var_ts = []
    incr = X.astype(float) - mu
    for t in ts:
        offsets = rng.integers(0, N - t + 1, size=200)
        m_t = np.array([incr[o:o + t].sum() for o in offsets])
        var_ts.append(m_t.var())
    ax.scatter(ts, var_ts, s=25, color="#1f77b4")
    # ajuste linear
    s2 = float((ts * var_ts).sum() / (ts * ts).sum())
    ax.plot(ts, s2 * ts, "--", color="red", lw=1.5,
            label=f"Var ≈ {s2:.4f}·t")
    ax.plot(ts, sigma2_iid * ts, ":", color="gray", lw=1,
            label=f"σ²_iid·t = {sigma2_iid:.4f}·t")
    ax.set_xlabel("t")
    ax.set_ylabel("Var[M_t]")
    ax.set_title("(c) Var[M_t] vs t (linear → comportamento martingale)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    # (d) QQ-plot de S_norm
    ax = axes[1, 1]
    samples = np.array(clt_block["samples"])
    samples_std = (samples - samples.mean()) / max(samples.std(ddof=1), 1e-9)
    stats.probplot(samples_std, dist="norm", plot=ax)
    ax.set_title(
        f"(d) QQ-plot M_T/√T  (Shapiro p={clt_block['shapiro_pvalue']:.3f})",
        fontsize=10,
    )
    ax.grid(alpha=0.3)

    fig.suptitle(
        f"Processo de martingale empírico  (n={n_label}, N={len(X)} eventos)",
        fontsize=12,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out_path, dpi=130)
    plt.close(fig)


# ----------------------------------------------------------------------
# main
# ----------------------------------------------------------------------
def main() -> None:
    print("Carregando sequências X^(n)...", flush=True)
    Xs = {n: load_X_seq(n) for n in NS}
    mus = {n: float(X.mean()) for n, X in Xs.items()}
    counts_tot = {n: int(len(X)) for n, X in Xs.items()}
    counts_tour = {n: int(X.sum()) for n, X in Xs.items()}

    print(f"  μ(n) = {mus}")

    print("\n[1] Teste de não-decaimento de μ(n)...", flush=True)
    decay = test_no_decay(mus, counts_tot)
    print(f"  slope = {decay['slope_ols']:+.5f}  "
          f"(SE = {decay['stderr_ols']:.5f})")
    print(f"  p-value (slope=0) = {decay['pvalue_slope_eq_0']:.4f}  "
          f"R² = {decay['r_squared']:.3f}")
    tag = "NÃO rejeita H0" if decay["pvalue_slope_eq_0"] > 0.05 else "REJEITA H0"
    print(f"  → {tag}")

    print("\n[2] Teste de variância de M_t (linear vs lei de potência)...",
          flush=True)
    var_per_n = {}
    r2_lin_list = []
    alphas = []
    r2_pow_list = []
    for n in NS:
        v = test_variance_linearity(Xs[n], mus[n])
        var_per_n[n] = v
        r2_lin_list.append(v["r2_linear"])
        alphas.append(v["power_alpha"])
        r2_pow_list.append(v["r2_power"])
        print(f"  n={n:>2}: σ²={v['sigma2_estimated']:.4f}  "
              f"σ²_iid={v['sigma2_iid']:.4f}  ratio={v['ratio_to_iid']:.2f}  "
              f"R²_lin={v['r2_linear']:.3f}  "
              f"α={v['power_alpha']:.2f}  R²_pow={v['r2_power']:.3f}")
    mean_r2_lin = float(np.mean(r2_lin_list))
    mean_alpha = float(np.mean(alphas))
    mean_r2_pow = float(np.mean(r2_pow_list))
    var_summary = {"per_n": var_per_n,
                   "mean_r2": mean_r2_lin,
                   "mean_alpha": mean_alpha,
                   "mean_r2_power": mean_r2_pow}
    print(f"  R²_lin médio = {mean_r2_lin:.3f}")
    print(f"  α médio       = {mean_alpha:.3f}  "
          f"(α=1: martingale i.i.d.; α<1: concentração sub-linear)")
    print(f"  R²_pow médio  = {mean_r2_pow:.3f}")

    print("\n[3] Bound inferior empírico de μ_min...", flush=True)
    mu_lb = test_mu_lower_bound(mus, counts_tour, counts_tot)
    for r in mu_lb["rows"]:
        print(f"  n={r['n']:>2}: μ={r['mu']:.4f}  "
              f"IC95=[{r['wilson_lo']:.4f}, {r['wilson_hi']:.4f}]")
    print(f"  min LB 95% = {mu_lb['min_lower_bound_95pct']:.4f}")
    print(f"  ⇒ E[nós/tour] ≤ {mu_lb['implied_E_nodes_per_tour_upper']:.2f}")
    if mu_lb["min_lower_bound_95pct"] >= 0.10:
        print("  → μ_min ≥ 0.10 CONFIRMADO (95% conf.) — bound O(1) ok")
    else:
        print("  → μ_min ABAIXO de 0.10 — teorema mais frouxo")

    print("\n[4] TCL para martingale (Shapiro+KS por n)...", flush=True)
    clt_per_n = {}
    shap_ps, ks_ps = [], []
    for n in NS:
        c = test_martingale_clt(Xs[n], mus[n], block_size=50)
        clt_per_n[n] = c
        if c["ok"]:
            shap_ps.append(c["shapiro_pvalue"])
            ks_ps.append(c["ks_pvalue"])
            print(f"  n={n:>2}: blocos={c['n_blocks']:>3}  "
                  f"Shapiro p={c['shapiro_pvalue']:.3f}  "
                  f"KS p={c['ks_pvalue']:.3f}")
        else:
            print(f"  n={n:>2}: {c.get('reason')}")
    mean_shap = float(np.mean(shap_ps)) if shap_ps else float("nan")
    mean_ks = float(np.mean(ks_ps)) if ks_ps else float("nan")
    clt_summary = {
        "per_n": clt_per_n,
        "mean_shapiro_p": mean_shap,
        "mean_ks_p": mean_ks,
    }
    print(f"  Shapiro p médio = {mean_shap:.3f}  KS p médio = {mean_ks:.3f}")

    # ------------------------------------------------------------------
    # Salvar
    # ------------------------------------------------------------------
    def to_json_safe(x):
        if isinstance(x, dict):
            return {str(k): to_json_safe(v) for k, v in x.items()}
        if isinstance(x, (list, tuple)):
            return [to_json_safe(v) for v in x]
        if isinstance(x, (np.integer,)):
            return int(x)
        if isinstance(x, (np.floating,)):
            return float(x)
        if isinstance(x, np.ndarray):
            return x.tolist()
        return x

    payload = {
        "test1_no_decay": decay,
        "test2_variance_linear": var_summary,
        "test3_mu_lower_bound": mu_lb,
        "test4_clt": clt_summary,
        "verdict": {
            "no_decay": decay["pvalue_slope_eq_0"] > 0.05,
            "variance_linear": mean_r2_lin > 0.9,
            "variance_power_alpha_near_1": 0.7 < mean_alpha < 1.3,
            "mu_min_geq_010": mu_lb["min_lower_bound_95pct"] >= 0.10,
            "clt_supported": mean_shap > 0.05 and mean_ks > 0.05,
        },
    }
    out_json = DATA / "convergence_tests.json"
    out_json.write_text(json.dumps(to_json_safe(payload), indent=2))
    print(f"\nSalvo {out_json}")

    # ------------------------------------------------------------------
    # Plots
    # ------------------------------------------------------------------
    print("Gerando plots...", flush=True)
    plot_convergence_evidence(
        decay, var_summary, mu_lb, clt_summary,
        PLOTS / "convergence_evidence.png",
    )
    plot_martingale_process(
        Xs[10], mus[10], n_label=10, clt_block=clt_per_n[10],
        out_path=PLOTS / "martingale_process.png",
    )
    print("  plots gerados.")

    # ------------------------------------------------------------------
    # Veredito final
    # ------------------------------------------------------------------
    print("\n" + "=" * 60)
    print("VEREDITO T3")
    print("=" * 60)
    v = payload["verdict"]
    print(f"  μ não decai com n ................ {'SIM' if v['no_decay'] else 'NÃO'}")
    print(f"  Var[M_t] linear (R²>0.9) ......... {'SIM' if v['variance_linear'] else 'NÃO'}")
    print(f"  α de Var ~ t^α ∈ [0.7,1.3] ....... {'SIM' if v['variance_power_alpha_near_1'] else 'NÃO'}  (α≈{mean_alpha:.2f})")
    print(f"  μ_min ≥ 0.10 (95% conf.) ......... {'SIM' if v['mu_min_geq_010'] else 'NÃO'}")
    print(f"  TCL martingale suportado ......... {'SIM' if v['clt_supported'] else 'NÃO'}")
    all_ok = all(v.values())
    print(f"\n  Teorema O(1) viável?  ......... "
          f"{'SIM (todos os 4 critérios passam)' if all_ok else 'PARCIAL — ver lacunas'}")


if __name__ == "__main__":
    main()
