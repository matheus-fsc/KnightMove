"""
T2 — análise condicional dupla.

Pergunta refinada (pós-T1):
  Em quais bins de estado o n=6 diverge dos n grandes?
  A divergência é efeito de borda (L baixo) ou distribuída?

Quatro análises:
  A) p(b, n) por variável de estado (depth, n_components, L, pressure)
  B) μ_interior (L≥2) vs μ_borda (L<2) — efeito de borda
  C) μ_padronizado (composição vs p genuíno)
  D) χ² de homogeneidade por bin (todos n e só n≥8)

Saídas:
  data/p_conditional_analysis.json
  data/plots/p_distribution_by_n.png
  data/plots/mu_decomposition.png
  data/plots/chi2_heatmap.png
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
MIN_EVENTS_FOR_CHI2 = 5    # mínimo por célula; regra de Cochran ≥5
WILSON_Z = 1.96            # 95% CI


# ----------------------------------------------------------------------
# Carregar logs e construir arrays planos
# ----------------------------------------------------------------------
def load_logs() -> dict[int, list[dict]]:
    out = {}
    for n in NS:
        path = DATA / f"event_log_n{n}.json"
        if not path.exists():
            print(f"!! falta {path}")
            continue
        out[n] = json.loads(path.read_text())
    return out


def to_arrays(events: list[dict]) -> dict[str, np.ndarray]:
    """Constrói os arrays do log e DESLOCA is_tour por +1: para o evento t,
    'target' = X_{t+1} (próximo evento é TOUR?). O último evento é descartado
    (não tem sucessor)."""
    n_ev = len(events)
    out = {
        "depth": np.empty(n_ev, dtype=np.float32),
        "n_components": np.empty(n_ev, dtype=np.int32),
        "max_degree_free": np.empty(n_ev, dtype=np.int32),
        "L_current": np.empty(n_ev, dtype=np.int32),
        "L_parent": np.empty(n_ev, dtype=np.int32),
        "score_pressure": np.empty(n_ev, dtype=np.float32),
        "is_tour": np.empty(n_ev, dtype=np.int32),
        "outcome": np.empty(n_ev, dtype=object),
    }
    for i, ev in enumerate(events):
        out["depth"][i] = ev["depth"]
        out["n_components"][i] = ev["n_components"]
        out["max_degree_free"][i] = ev["max_degree_free"]
        out["L_current"][i] = ev["L_current"]
        out["L_parent"][i] = ev.get("L_parent", -1)
        sp = ev["score_pressure"]
        out["score_pressure"][i] = float("nan") if sp is None else sp
        out["is_tour"][i] = int(ev["outcome"] == "TOUR")
        out["outcome"][i] = ev["outcome"]

    # truncar últimos para alinhar com target = X_{t+1}
    target = out["is_tour"][1:].copy()
    for k in out:
        out[k] = out[k][:-1]
    out["target"] = target
    # is_tour mantido como X_t (do próprio evento) p/ análise de μ_total
    return out


# ----------------------------------------------------------------------
# Bins
# ----------------------------------------------------------------------
DEPTH_BINS = [(0.0, 0.80), (0.80, 0.90), (0.90, 0.95), (0.95, 1.001)]
DEPTH_LABELS = ["[0,0.80)", "[0.80,0.90)", "[0.90,0.95)", "[0.95,1.0]"]

NC_BINS = [(1, 1), (2, 2), (3, 3), (4, 5), (6, 10), (11, 10**6)]
NC_LABELS = ["nc=1", "nc=2", "nc=3", "nc=4-5", "nc=6-10", "nc≥11"]

L_BINS = [(0, 0), (1, 1), (2, 2), (3, 3), (4, 4), (5, 10**6)]
L_LABELS = ["L=0", "L=1", "L=2", "L=3", "L=4", "L≥5"]


def bin_depth(x: float) -> int:
    for i, (lo, hi) in enumerate(DEPTH_BINS):
        if lo <= x < hi:
            return i
    return len(DEPTH_BINS) - 1


def bin_range(value: int, ranges: list[tuple[int, int]]) -> int:
    for i, (lo, hi) in enumerate(ranges):
        if lo <= value <= hi:
            return i
    return -1


# ----------------------------------------------------------------------
# Wilson interval
# ----------------------------------------------------------------------
def wilson(k: int, n: int, z: float = WILSON_Z) -> tuple[float, float, float]:
    if n == 0:
        return 0.0, 0.0, 1.0
    p = k / n
    denom = 1 + z * z / n
    center = (p + z * z / (2 * n)) / denom
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return p, max(0.0, center - half), min(1.0, center + half)


# ----------------------------------------------------------------------
# A) p(b, n) por variável de estado
# ----------------------------------------------------------------------
def conditional_by_variable(arrs_by_n: dict[int, dict]) -> dict:
    """Para cada variável de estado e cada bin, calcula p(b,n) com Wilson."""
    out: dict = {}

    # depth
    out["depth"] = {"labels": DEPTH_LABELS, "p_by_n": {}}
    for n, A in arrs_by_n.items():
        bins_idx = np.array([bin_depth(d) for d in A["depth"]])
        rows = []
        for i, lbl in enumerate(DEPTH_LABELS):
            mask = bins_idx == i
            cnt = int(mask.sum())
            tours = int(A["target"][mask].sum())
            p, lo, hi = wilson(tours, cnt)
            rows.append({"bin": lbl, "n_events": cnt, "n_tours": tours,
                         "p": p, "ci_lo": lo, "ci_hi": hi})
        out["depth"]["p_by_n"][n] = rows

    # n_components
    out["n_components"] = {"labels": NC_LABELS, "p_by_n": {}}
    for n, A in arrs_by_n.items():
        bins_idx = np.array([bin_range(v, NC_BINS) for v in A["n_components"]])
        rows = []
        for i, lbl in enumerate(NC_LABELS):
            mask = bins_idx == i
            cnt = int(mask.sum())
            tours = int(A["target"][mask].sum())
            p, lo, hi = wilson(tours, cnt)
            rows.append({"bin": lbl, "n_events": cnt, "n_tours": tours,
                         "p": p, "ci_lo": lo, "ci_hi": hi})
        out["n_components"]["p_by_n"][n] = rows

    # L_current (apenas eventos com L>=0; L=-1 são podas pré-escolha ou folhas)
    out["L_current"] = {"labels": L_LABELS, "p_by_n": {}}
    for n, A in arrs_by_n.items():
        mask_valid = A["L_current"] >= 0
        L_valid = A["L_current"][mask_valid]
        tour_valid = A["target"][mask_valid]
        bins_idx = np.array([bin_range(v, L_BINS) for v in L_valid])
        rows = []
        for i, lbl in enumerate(L_LABELS):
            mask = bins_idx == i
            cnt = int(mask.sum())
            tours = int(tour_valid[mask].sum())
            p, lo, hi = wilson(tours, cnt)
            rows.append({"bin": lbl, "n_events": cnt, "n_tours": tours,
                         "p": p, "ci_lo": lo, "ci_hi": hi})
        out["L_current"]["p_by_n"][n] = rows

    # L_parent (L do vértice escolhido no pai — definido para todos exceto raiz)
    out["L_parent"] = {"labels": L_LABELS, "p_by_n": {}}
    for n, A in arrs_by_n.items():
        mask_valid = A["L_parent"] >= 0
        L_valid = A["L_parent"][mask_valid]
        tour_valid = A["target"][mask_valid]
        bins_idx = np.array([bin_range(v, L_BINS) for v in L_valid])
        rows = []
        for i, lbl in enumerate(L_LABELS):
            mask = bins_idx == i
            cnt = int(mask.sum())
            tours = int(tour_valid[mask].sum())
            p, lo, hi = wilson(tours, cnt)
            rows.append({"bin": lbl, "n_events": cnt, "n_tours": tours,
                         "p": p, "ci_lo": lo, "ci_hi": hi})
        out["L_parent"]["p_by_n"][n] = rows

    # score_pressure — quintis por n (bins relativos)
    out["pressure"] = {"labels": ["Q1", "Q2", "Q3", "Q4", "Q5"], "p_by_n": {}}
    for n, A in arrs_by_n.items():
        mask = ~np.isnan(A["score_pressure"])
        sp = A["score_pressure"][mask]
        tr = A["target"][mask]
        if len(sp) < 5:
            out["pressure"]["p_by_n"][n] = []
            continue
        qs = np.quantile(sp, [0.2, 0.4, 0.6, 0.8])
        bins_idx = np.digitize(sp, qs)  # 0..4
        rows = []
        for i, lbl in enumerate(["Q1", "Q2", "Q3", "Q4", "Q5"]):
            mm = bins_idx == i
            cnt = int(mm.sum())
            tours = int(tr[mm].sum())
            p, lo, hi = wilson(tours, cnt)
            rows.append({"bin": lbl, "n_events": cnt, "n_tours": tours,
                         "p": p, "ci_lo": lo, "ci_hi": hi,
                         "qrange": [float(qs[i-1]) if i > 0 else None,
                                    float(qs[i]) if i < 4 else None]})
        out["pressure"]["p_by_n"][n] = rows

    return out


# ----------------------------------------------------------------------
# B) μ_interior (L_parent≥2) vs μ_borda (L_parent<2)
# ----------------------------------------------------------------------
def mu_interior_borda(arrs_by_n: dict[int, dict]) -> dict:
    """Para cada evento, classifica pelo L do vértice escolhido no PAI
    (L_parent ≥ 0). Eventos sem pai (raiz, L_parent=-1) ficam em 'root'.

    Usar L_parent — e não L_current — torna a tabela válida para eventos
    terminais (TOURs, podas) que não têm vértice escolhido próprio."""
    out: dict = {}
    for n, A in arrs_by_n.items():
        Lp = A["L_parent"]
        tr = A["target"]
        m_int = Lp >= 2
        m_bor = (Lp >= 0) & (Lp < 2)
        m_root = Lp < 0
        out[n] = {
            "interior": {
                "n_events": int(m_int.sum()),
                "n_tours": int(tr[m_int].sum()),
                "mu": float(tr[m_int].mean()) if m_int.sum() else float("nan"),
            },
            "borda": {
                "n_events": int(m_bor.sum()),
                "n_tours": int(tr[m_bor].sum()),
                "mu": float(tr[m_bor].mean()) if m_bor.sum() else float("nan"),
            },
            "root": {
                "n_events": int(m_root.sum()),
                "n_tours": int(tr[m_root].sum()),
                "mu": float(tr[m_root].mean()) if m_root.sum() else float("nan"),
            },
            "mu_total": float(tr.mean()),
        }
    return out


# ----------------------------------------------------------------------
# C) μ_padronizado — composição vs p genuíno
# ----------------------------------------------------------------------
def mu_padronizado(arrs_by_n: dict[int, dict],
                   p_by_var: dict,
                   ref_n: int = 10) -> dict:
    """Calcula μ padronizado usando a composição de bins de n=ref.
    Faz para depth, L_current e n_components. Pressure não faz sentido
    (quintis são relativos)."""
    out: dict = {}

    for var_name in ["depth", "L_current", "L_parent", "n_components"]:
        block = p_by_var[var_name]
        ref_rows = block["p_by_n"].get(ref_n, [])
        if not ref_rows:
            continue
        total_ref = sum(r["n_events"] for r in ref_rows)
        if total_ref == 0:
            continue
        weights_ref = [r["n_events"] / total_ref for r in ref_rows]

        out[var_name] = {"ref_n": ref_n,
                         "weights_ref": weights_ref,
                         "labels": block["labels"],
                         "mu_pad_by_n": {}}
        for n, rows in block["p_by_n"].items():
            if not rows or len(rows) != len(weights_ref):
                continue
            mu_pad = sum(r["p"] * w for r, w in zip(rows, weights_ref))
            out[var_name]["mu_pad_by_n"][n] = float(mu_pad)
    return out


# ----------------------------------------------------------------------
# D) χ² de homogeneidade por bin
# ----------------------------------------------------------------------
def chi2_homogeneity(p_by_var: dict, ns: list[int]) -> dict:
    """Para cada variável de estado e cada bin, χ² 2xN (TOUR vs não-TOUR
    por n). Retorna p-value por bin. Bins com menos de MIN_EVENTS_FOR_CHI2
    em algum n são marcados como inconclusivos."""
    out: dict = {}
    for var_name, block in p_by_var.items():
        labels = block["labels"]
        rows_per_n = [block["p_by_n"].get(n, []) for n in ns]
        if not all(rows_per_n):
            continue
        out[var_name] = {}
        for i, lbl in enumerate(labels):
            table = []  # [[TOUR, NOT_TOUR] per n]
            ok = True
            for j, n in enumerate(ns):
                rows = rows_per_n[j]
                if not rows or i >= len(rows):
                    ok = False
                    break
                r = rows[i]
                tour = r["n_tours"]
                not_t = r["n_events"] - r["n_tours"]
                if r["n_events"] < MIN_EVENTS_FOR_CHI2:
                    ok = False
                table.append([tour, not_t])
            if not ok:
                out[var_name][lbl] = {"pvalue": None, "ok": False,
                                      "table": table}
                continue
            arr = np.array(table)
            if arr.sum() == 0:
                out[var_name][lbl] = {"pvalue": None, "ok": False,
                                      "table": table}
                continue
            try:
                chi2, pval, dof, exp = stats.chi2_contingency(
                    arr, correction=False
                )
            except ValueError:
                out[var_name][lbl] = {"pvalue": None, "ok": False,
                                      "table": table}
                continue
            out[var_name][lbl] = {
                "pvalue": float(pval),
                "chi2": float(chi2),
                "dof": int(dof),
                "ok": True,
                "table": table,
            }
    return out


# ----------------------------------------------------------------------
# Plots
# ----------------------------------------------------------------------
def plot_p_distribution(p_by_var: dict, mu_total: dict,
                        chi2_all: dict, out_path: Path) -> None:
    vars_to_plot = ["depth", "n_components", "L_current", "L_parent",
                    "pressure"]
    fig, axes = plt.subplots(
        len(vars_to_plot), len(NS),
        figsize=(3.0 * len(NS), 2.6 * len(vars_to_plot)),
        sharey="row",
    )
    for ri, var in enumerate(vars_to_plot):
        block = p_by_var[var]
        labels = block["labels"]
        chi2_block = chi2_all.get(var, {})
        for ci, n in enumerate(NS):
            ax = axes[ri, ci]
            rows = block["p_by_n"].get(n, [])
            if not rows:
                ax.axis("off")
                continue
            xs = np.arange(len(rows))
            ps = [r["p"] for r in rows]
            lo = [max(0.0, r["p"] - r["ci_lo"]) for r in rows]
            hi = [max(0.0, r["ci_hi"] - r["p"]) for r in rows]
            ax.bar(xs, ps, color="#4a90d9", edgecolor="black",
                   linewidth=0.5, alpha=0.85)
            ax.errorbar(xs, ps, yerr=[lo, hi], fmt="none",
                        ecolor="black", capsize=2, linewidth=0.8)
            ax.axhline(mu_total[n], color="red", lw=1, ls="--",
                       label=f"μ(n)={mu_total[n]:.2f}")
            ax.set_xticks(xs)
            ax.set_xticklabels(labels, rotation=30, fontsize=7)
            ax.set_ylim(0, 1.0)
            if ri == 0:
                ax.set_title(f"n={n}", fontsize=10)
            if ci == 0:
                ax.set_ylabel(var, fontsize=9)
            ax.tick_params(labelsize=7)
            ax.legend(fontsize=6, loc="upper left")

            # estrelas para bins onde χ² rejeita
            for i, lbl in enumerate(labels):
                info = chi2_block.get(lbl)
                if info and info.get("ok") and info["pvalue"] is not None \
                        and info["pvalue"] < 0.05:
                    ax.text(xs[i], ps[i] + 0.04, "*", ha="center",
                            fontsize=9, color="darkred", fontweight="bold")

    fig.suptitle(
        "p(TOUR | bin) por variável de estado e n  "
        "(* = χ²(homog) p<0.05 entre n)",
        fontsize=11,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(out_path, dpi=130)
    plt.close(fig)


def plot_mu_decomposition(mu_total: dict,
                          mu_ib: dict,
                          mu_pad: dict,
                          out_path: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.0))

    # painel 1
    ax = axes[0]
    ns = sorted(mu_total)
    ys = [mu_total[n] for n in ns]
    ax.plot(ns, ys, "o-", color="black", lw=2)
    ax.set_title("μ_total(n)")
    ax.set_xlabel("n")
    ax.set_ylabel("μ")
    ax.set_ylim(0, 0.6)
    ax.grid(alpha=0.3)
    for n, y in zip(ns, ys):
        ax.annotate(f"{y:.3f}", (n, y), textcoords="offset points",
                    xytext=(0, 8), fontsize=8, ha="center")

    # painel 2 — interior vs borda
    ax = axes[1]
    ys_i = [mu_ib[n]["interior"]["mu"] for n in ns]
    ys_b = [mu_ib[n]["borda"]["mu"] for n in ns]
    ax.plot(ns, ys_i, "o-", color="#1f77b4", lw=2, label="μ_interior (L≥2)")
    ax.plot(ns, ys_b, "s-", color="#d62728", lw=2, label="μ_borda (L<2)")
    ax.set_title("μ por região")
    ax.set_xlabel("n")
    ax.set_ylim(0, 0.8)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=9)
    for n, y in zip(ns, ys_i):
        ax.annotate(f"{y:.3f}", (n, y), textcoords="offset points",
                    xytext=(0, 8), fontsize=7, ha="center", color="#1f77b4")
    for n, y in zip(ns, ys_b):
        if not math.isnan(y):
            ax.annotate(f"{y:.3f}", (n, y), textcoords="offset points",
                        xytext=(0, -12), fontsize=7, ha="center",
                        color="#d62728")

    # painel 3 — μ_padronizado
    ax = axes[2]
    ax.plot(ns, ys, "--", color="gray", lw=1.5, label="μ_total")
    for var, color in [("depth", "#1f77b4"),
                       ("L_parent", "#2ca02c"),
                       ("n_components", "#ff7f0e")]:
        if var not in mu_pad:
            continue
        d = mu_pad[var]["mu_pad_by_n"]
        ys_pad = [d.get(n, float("nan")) for n in ns]
        ax.plot(ns, ys_pad, "o-", color=color, lw=2,
                label=f"μ_pad (ref n=10, var={var})")
    ax.set_title("μ_padronizado (composição n=10)")
    ax.set_xlabel("n")
    ax.set_ylim(0, 0.6)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)

    fig.suptitle("Decomposição da variação de μ(n)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out_path, dpi=130)
    plt.close(fig)


def plot_chi2_heatmap(chi2_all: dict, chi2_n8: dict,
                      out_path: Path) -> None:
    vars_order = ["depth", "n_components", "L_current", "L_parent",
                  "pressure"]
    # construir grade: linhas = todos os bins concatenados; colunas = 2
    fig, axes = plt.subplots(1, 2, figsize=(11, 6))
    for col, (title, src) in enumerate([
        ("Todos n (6..14)", chi2_all),
        ("Apenas n ≥ 8", chi2_n8),
    ]):
        ax = axes[col]
        labels = []
        vals = []
        for var in vars_order:
            block = src.get(var, {})
            for lbl, info in block.items():
                labels.append(f"{var}: {lbl}")
                if info["pvalue"] is None or not info["ok"]:
                    vals.append(np.nan)
                else:
                    # -log10(p) clipped
                    p = max(info["pvalue"], 1e-50)
                    vals.append(-math.log10(p))
        arr = np.array(vals).reshape(-1, 1)
        im = ax.imshow(arr, aspect="auto", cmap="RdYlGn_r",
                       vmin=0, vmax=10)
        ax.set_yticks(range(len(labels)))
        ax.set_yticklabels(labels, fontsize=7)
        ax.set_xticks([])
        ax.set_title(title, fontsize=11)
        ax.axvline(-0.5, color="black", lw=0.5)
        # anotar valores
        for i, v in enumerate(vals):
            if math.isnan(v):
                ax.text(0, i, "n/a", ha="center", va="center",
                        fontsize=7, color="gray")
            else:
                color = "white" if v > 5 else "black"
                ax.text(0, i, f"{v:.2f}", ha="center", va="center",
                        fontsize=7, color=color)
        cbar = fig.colorbar(im, ax=ax, shrink=0.8)
        cbar.set_label("-log10(p)", fontsize=8)
        cbar.ax.axhline(1.3, color="black", lw=1)  # threshold p=0.05

    fig.suptitle(
        "χ² de homogeneidade por bin entre tabuleiros  "
        "(threshold -log10(0.05) ≈ 1.30 = preto)",
        fontsize=11,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out_path, dpi=130)
    plt.close(fig)


# ----------------------------------------------------------------------
# Tabela de spread
# ----------------------------------------------------------------------
def spread(values: list[float]) -> float:
    vs = [v for v in values if not (isinstance(v, float) and math.isnan(v))]
    if not vs:
        return float("nan")
    mx = max(vs)
    if mx <= 0:
        return float("nan")
    return (mx - min(vs)) / mx


# ----------------------------------------------------------------------
# main
# ----------------------------------------------------------------------
def main() -> None:
    print("Carregando event logs...", flush=True)
    logs = load_logs()
    print(f"  carregados n={list(logs.keys())}, "
          f"eventos totais = {sum(len(v) for v in logs.values())}",
          flush=True)

    arrs = {n: to_arrays(ev) for n, ev in logs.items()}
    # μ usado em toda análise condicional: target = X_{t+1}
    mu_total = {n: float(a["target"].mean()) for n, a in arrs.items()}

    print("\n[A] Distribuição condicional p(b,n)...", flush=True)
    p_by_var = conditional_by_variable(arrs)

    print("[B] μ_interior vs μ_borda...", flush=True)
    mu_ib = mu_interior_borda(arrs)

    print("[C] μ_padronizado (composição n=10)...", flush=True)
    mu_pad = mu_padronizado(arrs, p_by_var, ref_n=10)

    print("[D] χ² de homogeneidade...", flush=True)
    chi2_all = chi2_homogeneity(p_by_var, NS)
    chi2_n8 = chi2_homogeneity(p_by_var, [n for n in NS if n >= 8])

    # ------------------------------------------------------------------
    # Reportar em texto
    # ------------------------------------------------------------------
    print("\n" + "=" * 72)
    print("RESULTADO T2 — análise condicional dupla")
    print("=" * 72)

    print("\n--- μ_total(n) por tabuleiro ---")
    for n in NS:
        print(f"  n={n:>2}: μ = {mu_total[n]:.4f}")
    sp_all = spread(list(mu_total.values()))
    sp_n8 = spread([mu_total[n] for n in NS if n >= 8])
    print(f"  spread(todos)  = {sp_all*100:.1f}%")
    print(f"  spread(n≥8)    = {sp_n8*100:.1f}%")

    # B) interior / borda (classificado por L_parent)
    print("\n--- (B) μ por região do PAI (interior L_p≥2 / borda L_p<2) ---")
    print(f"  {'n':>3} {'μ_interior':>11} {'n_ev_int':>10} "
          f"{'μ_borda':>9} {'n_ev_bor':>10} {'μ_root':>9} {'n_root':>7}")
    mu_int_list, mu_bor_list = [], []
    for n in NS:
        d = mu_ib[n]
        mu_int_list.append(d["interior"]["mu"])
        mu_bor_list.append(d["borda"]["mu"])
        print(f"  {n:>3} {d['interior']['mu']:>11.4f} "
              f"{d['interior']['n_events']:>10} "
              f"{d['borda']['mu']:>9.4f} {d['borda']['n_events']:>10} "
              f"{d['root']['mu']:>9.4f} {d['root']['n_events']:>7}")
    sp_int = spread(mu_int_list)
    sp_bor = spread(mu_bor_list)
    print(f"  spread μ_interior = {sp_int*100:.1f}%  "
          f"(vs μ_total {sp_all*100:.1f}%)")
    print(f"  spread μ_borda    = {sp_bor*100:.1f}%")
    if not math.isnan(sp_int):
        print(f"  → interior {'MAIS' if sp_int < sp_all else 'menos'} "
              f"constante que total")

    # C) μ padronizado
    print("\n--- (C) μ_padronizado (composição n=10) ---")
    print(f"  {'n':>3} {'μ_total':>9} {'μ_pad/depth':>13} "
          f"{'μ_pad/L_par':>13} {'μ_pad/nc':>11}")
    for n in NS:
        row = [f"{n:>3}", f"{mu_total[n]:>9.4f}"]
        for var in ["depth", "L_parent", "n_components"]:
            v = mu_pad.get(var, {}).get("mu_pad_by_n", {}).get(n,
                                                              float("nan"))
            row.append(f"{v:>13.4f}" if not math.isnan(v) else f"{'-':>13}")
        print("  " + " ".join(row))
    for var in ["depth", "L_parent", "n_components"]:
        d = mu_pad.get(var, {}).get("mu_pad_by_n", {})
        if not d:
            continue
        sp_pad = spread(list(d.values()))
        print(f"  spread μ_pad/{var:<12} = {sp_pad*100:.1f}%  "
              f"(vs μ_total {sp_all*100:.1f}%)")

    # D) χ²
    print("\n--- (D) χ² de homogeneidade por bin (todos n) ---")
    for var in ["depth", "n_components", "L_current", "L_parent", "pressure"]:
        if var not in chi2_all:
            continue
        print(f"  [{var}]")
        for lbl, info in chi2_all[var].items():
            if not info.get("ok"):
                print(f"    {lbl:>14}: poucos dados / inconclusivo")
                continue
            pv = info["pvalue"]
            tag = "REJ" if pv < 0.05 else "ok"
            print(f"    {lbl:>14}: p={pv:.2e}  χ²={info['chi2']:.2f}  "
                  f"dof={info['dof']}  [{tag}]")

    print("\n--- (D') χ² apenas para n ≥ 8 ---")
    for var in ["depth", "n_components", "L_current", "L_parent", "pressure"]:
        if var not in chi2_n8:
            continue
        print(f"  [{var}]")
        for lbl, info in chi2_n8[var].items():
            if not info.get("ok"):
                print(f"    {lbl:>14}: poucos dados / inconclusivo")
                continue
            pv = info["pvalue"]
            tag = "REJ" if pv < 0.05 else "ok"
            print(f"    {lbl:>14}: p={pv:.2e}  [{tag}]")

    # ------------------------------------------------------------------
    # salvar
    # ------------------------------------------------------------------
    # converter chaves int->str para JSON
    def stringify_keys(d):
        if isinstance(d, dict):
            return {str(k): stringify_keys(v) for k, v in d.items()}
        if isinstance(d, list):
            return [stringify_keys(x) for x in d]
        if isinstance(d, (np.integer,)):
            return int(d)
        if isinstance(d, (np.floating,)):
            return float(d)
        return d

    payload = {
        "mu_total": mu_total,
        "by_variable": p_by_var,
        "mu_interior_borda": mu_ib,
        "mu_padronizado": mu_pad,
        "chi2_all": chi2_all,
        "chi2_n8plus": chi2_n8,
        "spread_summary": {
            "mu_total_all": sp_all,
            "mu_total_n8": sp_n8,
            "mu_interior": spread(mu_int_list),
            "mu_borda": spread(mu_bor_list),
        },
    }
    out_json = DATA / "p_conditional_analysis.json"
    out_json.write_text(json.dumps(stringify_keys(payload), indent=2,
                                   default=str))
    print(f"\nSalvo {out_json}")

    # ------------------------------------------------------------------
    # plots
    # ------------------------------------------------------------------
    print("\nGerando plots...", flush=True)
    plot_p_distribution(p_by_var, mu_total, chi2_all,
                        PLOTS / "p_distribution_by_n.png")
    plot_mu_decomposition(mu_total, mu_ib, mu_pad,
                          PLOTS / "mu_decomposition.png")
    plot_chi2_heatmap(chi2_all, chi2_n8,
                      PLOTS / "chi2_heatmap.png")
    print("  plots gerados em data/plots/")


if __name__ == "__main__":
    main()
