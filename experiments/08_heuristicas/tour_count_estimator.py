"""tour_count_estimator.py — Estimador de Knuth para N(n).

Estima o número total de tours fechados do cavalo no tabuleiro n×n
percorrendo caminhos aleatórios na árvore de busca do knight_tours.py
e tomando a média do produto dos branching factors.

E[W] = N(n)  onde W = Π_{nó no caminho} b(nó).

Dependências: numpy, matplotlib (opcional para plot), knight_tours.py.
"""

# --- bootstrap: torna core/ importavel a partir de qualquer subpasta ---
import sys as _sys, pathlib as _pathlib
for _anc in _pathlib.Path(__file__).resolve().parents:
    if (_anc / "core").is_dir() and (_anc / "pyproject.toml").exists():
        _sys.path.insert(0, str(_anc / "core"))
        break
# --- fim do bootstrap ---

import json
import os
import time
import numpy as np

from knight_tours import (
    build_graph, State,
    _propagate_initial, _choose_next_edge,
    fix_and_propagate, _is_complete_tour,
    OK, CONTRADICTION, SUBTOUR, COMPLETE_TOUR,
)


# ---------------------------------------------------------------------------
# Núcleo: caminho aleatório na árvore de busca
# ---------------------------------------------------------------------------

def random_path_weight(n, ctx, rng):
    """Percorre um caminho aleatório raiz→folha na árvore de busca.

    Retorna W = Π_{nó} b(nó), onde b(nó) ∈ {0, 1, 2} é o número de
    filhos válidos do nó (valores 0/1 cuja propagação R2 não produz
    CONTRADICTION ou SUBTOUR dead-end).

    E[W] = N(n) (estimador de Knuth, não-viesado).
    """
    state = State(ctx['V'], ctx['E'])
    status = _propagate_initial(state, ctx)
    if status == CONTRADICTION or status == SUBTOUR:
        return 0.0
    if status == COMPLETE_TOUR:
        return 1.0 if _is_complete_tour(state, ctx) else 0.0

    weight = 1.0
    while True:
        e, _ = _choose_next_edge(state, ctx, rng)
        if e == -1:
            return weight if _is_complete_tour(state, ctx) else 0.0

        valid = []
        for val in (0, 1):
            snap = state.snapshot()
            s = fix_and_propagate(state, ctx, e, val)
            if s == OK:
                valid.append((val, False))
            elif s == COMPLETE_TOUR and _is_complete_tour(state, ctx):
                valid.append((val, True))
            state.restore(snap)

        b = len(valid)
        if b == 0:
            return 0.0
        weight *= b

        idx = 0 if b == 1 else int(rng.integers(b))
        chosen_val, is_complete = valid[idx]
        fix_and_propagate(state, ctx, e, chosen_val)
        if is_complete:
            return weight


# ---------------------------------------------------------------------------
# Estimador com IC
# ---------------------------------------------------------------------------

def estimate_tour_count(n, M=10000, seed=None, verbose=True,
                        bootstrap_B=2000):
    """Estima N(n) com IC 95% via normal e bootstrap.

    Parâmetros:
      n:           tamanho do tabuleiro
      M:           amostras do estimador
      seed:        semente aleatória
      verbose:     imprime progresso a cada 10% das amostras
      bootstrap_B: réplicas do bootstrap para IC robusto (0 = desliga)

    Retorna dict com estatísticas.
    """
    rng = np.random.default_rng(seed)
    ctx = build_graph(n)

    weights = np.zeros(M)
    t0 = time.time()
    report_every = max(1, M // 10)
    for i in range(M):
        weights[i] = random_path_weight(n, ctx, rng)
        if verbose and (i + 1) % report_every == 0:
            seen = weights[:i + 1]
            nz = seen[seen > 0]
            mu = float(np.mean(seen))
            elapsed = time.time() - t0
            frac0 = 100.0 * (i + 1 - len(nz)) / (i + 1)
            print(f"    [{i+1:6d}/{M}] E[W] ≈ {mu:.3e}  "
                  f"zeros={frac0:5.1f}%  t={elapsed:6.1f}s")
    elapsed = time.time() - t0

    nonzero = weights[weights > 0]
    mu = float(np.mean(weights))
    sigma = float(np.std(weights, ddof=1)) if M > 1 else 0.0
    se = sigma / np.sqrt(M) if M > 0 else 0.0

    ci_low = max(0.0, mu - 1.96 * se)
    ci_high = mu + 1.96 * se

    # IC log-normal (sobre a geom. mean dos não-zeros)
    if len(nonzero) > 10:
        log_w = np.log(nonzero)
        log_mu = float(np.mean(log_w))
        log_std = float(np.std(log_w, ddof=1))
        log_se = log_std / np.sqrt(len(nonzero))
        ci_log_low = float(np.exp(log_mu - 1.96 * log_se))
        ci_log_high = float(np.exp(log_mu + 1.96 * log_se))
    else:
        ci_log_low = ci_log_high = mu

    # IC bootstrap-percentil da média (robusto para cauda pesada)
    if bootstrap_B and M >= 30:
        bs_rng = np.random.default_rng(
            None if seed is None else int(seed) + 7919)
        boot_means = np.empty(bootstrap_B)
        for b in range(bootstrap_B):
            idx = bs_rng.integers(0, M, size=M)
            boot_means[b] = float(np.mean(weights[idx]))
        ci_boot_low = float(np.percentile(boot_means, 2.5))
        ci_boot_high = float(np.percentile(boot_means, 97.5))
    else:
        ci_boot_low = ci_boot_high = mu

    def safe_log10(x):
        return float(np.log10(x)) if x > 0 else 0.0

    return {
        "n": n,
        "M": M,
        "estimate": mu,
        "std": sigma,
        "se": float(se),
        "cv": float(sigma / mu) if mu > 0 else float('inf'),
        "ci_95_low": float(ci_low),
        "ci_95_high": float(ci_high),
        "ci_log_low": float(ci_log_low),
        "ci_log_high": float(ci_log_high),
        "ci_log10": [safe_log10(ci_log_low), safe_log10(ci_log_high)],
        "ci_boot_low": ci_boot_low,
        "ci_boot_high": ci_boot_high,
        "ci_boot_log10": [safe_log10(ci_boot_low), safe_log10(ci_boot_high)],
        "n_zeros": int(np.sum(weights == 0)),
        "frac_zeros": float(np.mean(weights == 0)),
        "n_nonzero": int(len(nonzero)),
        "log10_estimate": safe_log10(mu),
        "elapsed_s": float(elapsed),
        "time_per_sample_ms": float(1000.0 * elapsed / M) if M > 0 else 0.0,
        "max_weight": float(np.max(weights)) if M > 0 else 0.0,
    }


# ---------------------------------------------------------------------------
# Validação n=6 (ground truth = 9862)
# ---------------------------------------------------------------------------

GROUND_TRUTH_N6 = 9862  # número exato de tours fechados em 6×6 (orientados)


def validate_estimator(verbose=True):
    """Valida o estimador em n=6 onde N(6) é conhecido."""
    print("=" * 64)
    print("Validação n=6 — ground truth N(6) = 9.862 tours fechados")
    print("=" * 64)
    results = {}
    for M in (100, 1000, 5000, 10000):
        r = estimate_tour_count(6, M=M, seed=0, verbose=False,
                                bootstrap_B=1000 if M >= 1000 else 0)
        results[M] = r
        contains = r['ci_95_low'] <= GROUND_TRUTH_N6 <= r['ci_95_high']
        contains_boot = r['ci_boot_low'] <= GROUND_TRUTH_N6 <= r['ci_boot_high']
        flag = "OK" if (contains or contains_boot) else "FAIL"
        print(f"  M={M:6d}: Ê = {r['estimate']:9.1f}  "
              f"IC95 = [{r['ci_95_low']:.0f}, {r['ci_95_high']:.0f}]  "
              f"boot = [{r['ci_boot_low']:.0f}, {r['ci_boot_high']:.0f}]  "
              f"zeros={r['frac_zeros']*100:5.1f}%  "
              f"CV={r['cv']:5.2f}  [{flag}]")
    print(f"  ground truth: {GROUND_TRUTH_N6}")
    print()
    return results


# ---------------------------------------------------------------------------
# Execução completa: n ∈ {6, 8, 10, 12}
# ---------------------------------------------------------------------------

def run_all_estimates(M6=30000, M8=30000, M10=30000, M12=10000, seed=0):
    results = {}

    # n=6: confirmação contra ground truth
    print(f"=== n=6 (ground truth: 9.862)  M={M6} ===")
    results[6] = estimate_tour_count(6, M=M6, seed=seed, verbose=True)
    _report(results[6], extra=f"GT=9.862  log10(GT)={np.log10(9862):.3f}")

    # n=8: comparação com literatura
    print(f"\n=== n=8 (literatura: 13.267.364.410.532 ≈ 1.33e13)  M={M8} ===")
    results[8] = estimate_tour_count(8, M=M8, seed=seed, verbose=True)
    _report(results[8], extra=f"literatura log10 ≈ {np.log10(1.33e13):.3f}")

    # n=10: resultado principal
    print(f"\n=== n=10 (desconhecido)  M={M10} ===")
    results[10] = estimate_tour_count(10, M=M10, seed=seed, verbose=True)
    _report(results[10])

    # n=12: bonus
    print(f"\n=== n=12  M={M12} ===")
    results[12] = estimate_tour_count(12, M=M12, seed=seed, verbose=True)
    _report(results[12])

    return results


def _report(r, extra=""):
    n = r['n']
    print(f"  Ê               = {r['estimate']:.3e}")
    print(f"  log10(Ê)        = {r['log10_estimate']:.3f}")
    print(f"  CV              = {r['cv']:.2f}")
    print(f"  frac zeros      = {r['frac_zeros']*100:.1f}%")
    print(f"  IC95 normal     = [{r['ci_95_low']:.3e}, {r['ci_95_high']:.3e}]")
    print(f"  IC95 bootstrap  = [{r['ci_boot_low']:.3e}, {r['ci_boot_high']:.3e}]")
    print(f"  IC95 log10 boot = [{r['ci_boot_log10'][0]:.3f}, "
          f"{r['ci_boot_log10'][1]:.3f}]  "
          f"largura={r['ci_boot_log10'][1]-r['ci_boot_log10'][0]:.3f} dex")
    print(f"  tempo           = {r['elapsed_s']:.1f}s "
          f"({r['time_per_sample_ms']:.2f} ms/amostra)")
    if extra:
        print(f"  obs             = {extra}")


# ---------------------------------------------------------------------------
# Persistência e plot
# ---------------------------------------------------------------------------

def save_results(results, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    out = {f"n{n}": r for n, r in results.items()}
    with open(path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"  → resultados salvos em {path}")


def plot_results(results, path):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:
        print(f"  (matplotlib indisponível, plot pulado: {exc})")
        return

    os.makedirs(os.path.dirname(path), exist_ok=True)
    ns = sorted(results.keys())
    fig, ax = plt.subplots(figsize=(9, 6))

    for n in ns:
        r = results[n]
        log_est = r['log10_estimate']
        lo, hi = r['ci_boot_log10']
        if n == 6:
            ax.scatter([n], [np.log10(GROUND_TRUTH_N6)], s=110, marker='*',
                       color='black', zorder=5, label='n=6 (exato)')
            ax.errorbar([n], [log_est],
                        yerr=[[log_est - lo], [hi - log_est]],
                        fmt='o', color='tab:blue', alpha=0.6,
                        capsize=4, label='estimativa Knuth')
        elif n == 8:
            lit = np.log10(1.33e13)
            ax.scatter([n], [lit], s=80, marker='s',
                       color='tab:green', zorder=5, label='n=8 (literatura)')
            ax.errorbar([n], [log_est],
                        yerr=[[log_est - lo], [hi - log_est]],
                        fmt='o', color='tab:blue', alpha=0.6, capsize=4)
        elif n == 10:
            ax.errorbar([n], [log_est],
                        yerr=[[log_est - lo], [hi - log_est]],
                        fmt='D', color='tab:red', markersize=10,
                        capsize=5, lw=2, zorder=6,
                        label='n=10 (resultado principal)')
        else:
            ax.errorbar([n], [log_est],
                        yerr=[[log_est - lo], [hi - log_est]],
                        fmt='o', color='tab:blue', alpha=0.6, capsize=4)

    ax.set_xlabel('tamanho n do tabuleiro')
    ax.set_ylabel(r'$\log_{10} N(n)$')
    ax.set_title('Estimador de Knuth — tours fechados do cavalo em '
                 'tabuleiros n×n\n(barras = IC 95% bootstrap)')
    ax.grid(True, alpha=0.3)
    ax.set_xticks(ns)
    ax.legend(loc='lower right')

    plt.tight_layout()
    plt.savefig(path, dpi=130)
    plt.close(fig)
    print(f"  → plot salvo em {path}")


def print_summary(results):
    print()
    print("=" * 72)
    print("Estimativa do número de tours fechados do cavalo")
    print("=" * 72)
    print(f"  n=6 :  9.862                                        "
          f"(exato, ground truth)")
    for n in (8, 10, 12):
        if n not in results:
            continue
        r = results[n]
        lo, hi = r['ci_boot_low'], r['ci_boot_high']
        print(f"  n={n:<2d}:  {r['estimate']:.3e}   "
              f"IC95 = [{lo:.3e}, {hi:.3e}]   "
              f"CV={r['cv']:.2f}  "
              f"log10 ∈ [{r['ci_boot_log10'][0]:.2f}, "
              f"{r['ci_boot_log10'][1]:.2f}]")
    print()


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    # T3 — validação obrigatória
    val = validate_estimator()

    # critério: IC95 (normal ou bootstrap) deve conter 9862 com M=10000
    r10k = val[10000]
    inside = (r10k['ci_95_low'] <= GROUND_TRUTH_N6 <= r10k['ci_95_high']) or \
             (r10k['ci_boot_low'] <= GROUND_TRUTH_N6 <= r10k['ci_boot_high'])
    if not inside:
        print("VALIDAÇÃO FALHOU: IC95 (normal e bootstrap) com M=10000 "
              "não contém 9862. Abortando.")
        return

    # T4 — estimativas (M ajustado por n; n=10 é o resultado principal)
    results = run_all_estimates(M6=30000, M8=30000, M10=30000,
                                M12=10000, seed=0)

    # T6 — persistência e plot
    out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
    save_results(results, os.path.join(out_dir, "tour_count_estimates.json"))
    plot_results(results,
                 os.path.join(out_dir, "plots", "tour_count_estimates.png"))

    print_summary(results)


if __name__ == '__main__':
    main()
