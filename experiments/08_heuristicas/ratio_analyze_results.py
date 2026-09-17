"""Análise dos resultados de ratio_analysis.

Lê data/ratio_main_results.json + data/residual_enumeration.json (def B p/ n=6).
Gera:
  - relatório textual (T3..T7)
  - data/plots/ratio_analysis.png

Uso: python ratio_analyze_results.py
"""

import json
import os
import sys
import numpy as np


def load_results():
    with open('data/ratio_main_results.json') as f:
        raw = json.load(f)
    # Normalizar para dict[int → registro]
    out = {}
    for k, v in raw.items():
        try:
            n = int(k)
            out[n] = v
        except ValueError:
            pass
    # Definição B (n=6, R1..R6)
    try:
        with open('residual_search/data/residual_enumeration.json') as f:
            b6 = json.load(f)
    except FileNotFoundError:
        b6 = None
    return out, b6, raw


def report(results, b6, raw):
    print("=" * 70)
    print("RELATÓRIO FINAL — razão tours/2-fatores")
    print("=" * 70)

    print()
    print("=== Tabela resumo (Definição A: 2-fator padrão grau-2) ===")
    print(f"{'n':>3} | {'modo':>8} | {'N(2-fat)':>10} | {'N(tours)':>8} | "
          f"{'ratio':>8} | {'log2(1/r)':>10}")
    print("-" * 64)
    for n in sorted(results.keys()):
        r = results[n]
        if r['mode'] == 'exact':
            ratio = r['ratio']
            mode_s = 'EXATO'
            n2f, nt = r['n_2fatores'], r['n_tours']
        else:
            ratio = r['ratio_mean']
            mode_s = f'RST×{len(r["seeds"])}'
            n2f = int(np.mean(r['n_2fatores_per_seed']))
            nt = int(np.mean(r['n_tours_per_seed']))
        log2_inv = np.log2(1 / ratio) if ratio > 0 else float('inf')
        std_s = f'±{r["ratio_std"]:.4f}' if r['mode'] == 'restart' else ''
        print(f"{n:>3} | {mode_s:>8} | {n2f:>10d} | {nt:>8d} | "
              f"{ratio:>5.4f}{std_s:>0} | {log2_inv:>10.3f}")
    print()

    # === Definição B (n=6, R1..R6) ===
    if b6 is not None:
        print("=== Definição B (filtro R1..R6, somente n=6) ===")
        r_b6 = b6['n_tours_conexos'] / b6['n_2fatores']
        print(f"  N(2-fat | R1..R6) = {b6['n_2fatores']}")
        print(f"  N(tours)          = {b6['n_tours_conexos']}")
        print(f"  ratio_B(6)        = {r_b6:.4f}")
        print(f"  log2(1/r)         = {np.log2(1/r_b6):.3f}")
        print(f"  (não-extensível a n>6 sem precomputar exclusões)")
        print()

    # === T3 — distribuição de componentes ===
    print("=== T3 — Distribuição de componentes ===")
    print(f"{'n':>3} | {'comp=1':>8} | {'comp=2':>8} | {'comp=3':>8} | "
          f"{'comp=4':>8} | {'comp≥5':>8}")
    print("-" * 56)
    for n in sorted(results.keys()):
        r = results[n]
        if r['mode'] == 'exact':
            dist = r['component_distribution']
            total = r['n_2fatores']
        else:
            dist = r.get('comp_dist', r.get('comp_dist_seed42', {}))
            total = sum(int(v) for v in dist.values())
        dist = {int(k): int(v) for k, v in dist.items()}
        if total == 0:
            continue
        c = lambda k: dist.get(k, 0) / total
        c5p = sum(v for k, v in dist.items() if k >= 5) / total
        print(f"{n:>3} | {c(1):>8.4f} | {c(2):>8.4f} | {c(3):>8.4f} | "
              f"{c(4):>8.4f} | {c5p:>8.4f}")
    print()
    print("Observação: comp=2 dominante → falhas de conectividade")
    print("são predominantemente fusões de 2 ciclos.")
    print()

    # === T4 — padrão da razão ===
    print("=== T4 — Padrão de razão(n) ===")
    ns = sorted([n for n in results.keys()])
    ratios = []
    for n in ns:
        r = results[n]
        v = r['ratio'] if r['mode'] == 'exact' else r['ratio_mean']
        ratios.append(v)
    ns_arr = np.array(ns)
    rs_arr = np.array(ratios)

    print(f"n         : {ns_arr}")
    print(f"ratio     : {[round(r, 4) for r in rs_arr]}")
    log2_strs = [round(float(np.log2(1/r)), 3) if r > 0 else float('inf')
                  for r in rs_arr]
    print(f"log2(1/r) : {log2_strs}")
    print()

    valid = rs_arr > 0
    if valid.sum() >= 2:
        log_r = np.log(rs_arr[valid])
        ns_v = ns_arr[valid]

        # Fit linear em n
        s, b = np.polyfit(ns_v, log_r, 1)
        print(f"Fit em n:  ratio(n) ≈ {np.exp(b):.4f} × exp({s:.4f}·n)")
        for n_p in [6, 8, 10, 12, 14, 16]:
            pred = np.exp(b + s * n_p)
            mk = ' ✓' if n_p in ns else ''
            print(f"  n={n_p:2d}: pred={pred:.4f}{mk}")

        # Fit em n²
        s2, b2 = np.polyfit(ns_v ** 2, log_r, 1)
        print()
        print(f"Fit em n²: ratio(n) ≈ {np.exp(b2):.4f} × exp({s2:.6f}·n²)")

        # Fit constante (apenas média)
        log_const = np.mean(log_r)
        mean_ratio = np.exp(log_const)
        residuals_const = log_r - log_const
        residuals_linear = log_r - (b + s * ns_v)
        rss_const = float(np.sum(residuals_const ** 2))
        rss_linear = float(np.sum(residuals_linear ** 2))
        print()
        print(f"Modelo constante: ratio ≈ {mean_ratio:.4f}, RSS={rss_const:.4f}")
        print(f"Modelo linear-em-n:                          RSS={rss_linear:.4f}")

        print()
        print("Classificação automática:")
        if abs(s) < 0.05 and rss_const < 0.1:
            print("  → CASO 1: ratio aproximadamente CONSTANTE em n")
        elif s < -0.1:
            print(f"  → CASO 3: ratio DECAI exponencialmente (slope={s:.3f})")
        else:
            print(f"  → CASO 2: ratio varia moderadamente (slope={s:.3f})")

        fit = {'slope_n': float(s), 'intercept_n': float(b),
               'slope_n2': float(s2), 'intercept_n2': float(b2),
               'mean_log': float(log_const), 'rss_const': rss_const,
               'rss_linear': rss_linear}
    else:
        fit = {}
    print()

    # === T5 — deficit=3 ===
    print("=== T5 — Conexão com deficit=3 ===")
    print(f"2^(-3) = {2**(-3):.4f} = 12.5%")
    print()
    print(f"{'n':>3} | {'ratio':>8} | {'log2(1/r)':>10} | {'extra além 3':>12}")
    for n, r in zip(ns_arr, rs_arr):
        if r <= 0:
            continue
        li = np.log2(1 / r)
        print(f"{n:>3} | {r:>8.4f} | {li:>10.3f} | {li - 3:>+12.3f}")
    print()
    print("log2(1/r) ≈ 3 e estável → deficit=3 explica r(n)")
    print("log2(1/r) cresce com n → conectividade adiciona dimensões")
    print()

    # === T6 — N implícito ===
    print("=== T6 — N(2-fatores) implícito via N(tours) e ratio ===")
    lit = {6: 9862, 8: 1.33e13, 10: None, 12: None}
    for n in ns:
        r = results[n]
        ratio = r['ratio'] if r['mode'] == 'exact' else r['ratio_mean']
        if ratio <= 0:
            continue
        if lit.get(n) is None:
            print(f"  n={n:2d}: ratio={ratio:.4f}  (sem N(tours) lit.)")
            continue
        n_imp = lit[n] / ratio
        log_imp = np.log10(n_imp)
        print(f"  n={n:2d}: ratio={ratio:.4f}  N(tours)_lit={lit[n]:.3e}  "
              f"N(2-fat)_imp={n_imp:.3e}  log10={log_imp:.3f}")
    print()

    return fit


def make_plots(results, b6, fit):
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib indisponível — pulando plots")
        return

    ns = sorted(results.keys())
    ratios = []
    stds = []
    for n in ns:
        r = results[n]
        if r['mode'] == 'exact':
            ratios.append(r['ratio'])
            stds.append(0.0)
        else:
            ratios.append(r['ratio_mean'])
            stds.append(r['ratio_std'])
    ns_arr = np.array(ns)
    rs_arr = np.array(ratios)
    err_arr = np.array(stds)

    os.makedirs('data/plots', exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(17, 7.5))

    # Plot 1: ratio(n) com barras de erro
    axes[0].errorbar(ns_arr, rs_arr, yerr=err_arr, fmt='o-',
                      color='steelblue', linewidth=2, markersize=10,
                      capsize=5, label='ratio (def A)')
    if b6 is not None:
        r_b6 = b6['n_tours_conexos'] / b6['n_2fatores']
        axes[0].scatter([6], [r_b6], color='purple', s=120, marker='s',
                         label=f'def B (n=6): {r_b6:.3f}', zorder=5)
    axes[0].axhline(y=0.125, color='red', linestyle='--', alpha=0.6,
                     label='1/8 = 12.5% (deficit=3)')
    if fit:
        ns_smooth = np.linspace(min(ns), max(ns) + 4, 100)
        pred = np.exp(fit['intercept_n'] + fit['slope_n'] * ns_smooth)
        axes[0].plot(ns_smooth, pred, ':', color='orange', alpha=0.7,
                     label=f"fit exp({fit['slope_n']:.3f}·n)")
    axes[0].set_xlabel('n', fontsize=12)
    axes[0].set_ylabel('ratio = tours / 2-fatores', fontsize=12)
    axes[0].set_title('ratio(n) — Definição A (formal)', fontsize=12)
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    # Plot 2: log2(1/r)
    log2_inv = np.array([np.log2(1/r) if r > 0 else np.nan for r in rs_arr])
    axes[1].plot(ns_arr, log2_inv, 'o-', color='darkorange',
                  linewidth=2, markersize=10)
    axes[1].axhline(y=3, color='black', linestyle='--', alpha=0.6,
                     label='deficit = 3')
    axes[1].set_xlabel('n', fontsize=12)
    axes[1].set_ylabel('log2(1 / ratio)', fontsize=12)
    axes[1].set_title('dimensões "proibidas" por conectividade', fontsize=12)
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)

    # Plot 3: distribuição de componentes — todos os n
    width = 0.18
    x_base = np.array([1, 2, 3, 4, 5])  # comp=1,2,3,4,≥5
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
    for i, n in enumerate(sorted(results.keys())):
        r = results[n]
        if r['mode'] == 'exact':
            dist = r['component_distribution']
            total = r['n_2fatores']
        else:
            dist = r.get('comp_dist', r.get('comp_dist_seed42', {}))
            total = sum(int(v) for v in dist.values())
        if total == 0:
            continue
        dist = {int(k): int(v) for k, v in dist.items()}
        vals = [dist.get(1, 0)/total, dist.get(2, 0)/total,
                dist.get(3, 0)/total, dist.get(4, 0)/total,
                sum(v for k, v in dist.items() if k >= 5)/total]
        axes[2].bar(x_base + (i - 1.5) * width, vals, width,
                     label=f'n={n}', color=colors[i % 4])
    axes[2].set_xticks(x_base)
    axes[2].set_xticklabels(['1', '2', '3', '4', '≥5'])
    axes[2].set_xlabel('número de componentes', fontsize=12)
    axes[2].set_ylabel('fração de 2-fatores', fontsize=12)
    axes[2].set_title('distribuição de componentes por n', fontsize=12)
    axes[2].legend()
    axes[2].grid(True, axis='y', alpha=0.3)

    plt.tight_layout()
    out_path = 'data/plots/ratio_analysis.png'
    plt.savefig(out_path, dpi=150)
    print(f"Plot salvo: {out_path}")


def main():
    results, b6, raw = load_results()
    fit = report(results, b6, raw)
    make_plots(results, b6, fit)

    summary = {
        'definicao_A': {str(n): {
            'mode': r['mode'],
            'ratio': r['ratio'] if r['mode'] == 'exact' else r['ratio_mean'],
            'ratio_std': r['ratio_std'] if r['mode'] == 'restart' else 0.0,
        } for n, r in results.items()},
        'definicao_B_n6': (b6 if b6 else None),
        'fit': fit,
    }
    with open('data/ratio_summary.json', 'w') as f:
        json.dump(summary, f, indent=2, default=str)
    print("Salvo: data/ratio_summary.json")


if __name__ == '__main__':
    main()
