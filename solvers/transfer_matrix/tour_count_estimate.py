"""
Estimativa do numero de 2-fatores via potencia dominante: N(n) ~ C * lam_1(n)^(n^2).

Calibra C usando n=6 (conhecido por DP: 36236). Estimar para n=8, 10.

Reporta intervalo de confianca baseado em correcoes ρ^(n^2).
"""

import json
import sys
from pathlib import Path

import numpy as np

THIS = Path(__file__).parent
DATA = THIS / "data"


def load_eig(n):
    with open(DATA / f"eigenvalues_n{n}.json") as f:
        return json.load(f)


def get_count(n):
    """Conta exata de 2-fatores (do state_count.json)."""
    with open(DATA / "state_count.json") as f:
        sc = json.load(f)
    for r in sc:
        if r["n"] == n:
            return r["count_dp"] or r["count_matrix"]
    return None


def main(ns):
    """A contagem de 2-fatores escala como
        count(n) = <s_0 | T_bulk^{n-2} . T_nocp2 | s_target>
              ~ C(n) * lambda_1(n)^{n-1}    [efeitos de borda em C(n)]
    A potencia n^2 da especifica do usuario seria correta apenas se
    lambda_1(n) crescesse exponencialmente em n (o que de fato ocorre,
    pois a matriz tem ~9^n entries). Vamos reportar AMBAS as calibracoes."""
    rows = []
    for n in ns:
        eig = load_eig(n)
        lam1 = eig["eigenvalues_abs"][0]
        lam2 = eig["eigenvalues_abs"][1]
        rho = lam2 / lam1
        cnt = get_count(n)
        rows.append({
            "n": n,
            "lam1": lam1,
            "lam2": lam2,
            "rho": rho,
            "log_lam1": np.log(lam1),
            "count": cnt,
        })

    if not any(r["count"] for r in rows):
        print("Nenhum count exato disponivel para calibracao.")
        return

    # Calibracao A: count ~ C * lam1^(n-1)
    # Calibracao B: count ~ C * lam1^n
    # Calibracao C: count ~ C * lam1^(n^2)
    print("\n=== Calibracoes ===")
    print(f"{'n':>3} {'lam1':>10} {'count':>12} {'log_count':>12} "
          f"{'log_lam^(n-1)':>14} {'log_lam^(n)':>14} {'log_lam^(n^2)':>14}")
    for r in rows:
        if r["count"] is None:
            continue
        n = r["n"]
        lc = np.log(r["count"])
        e1 = (n - 1) * r["log_lam1"]
        e2 = n * r["log_lam1"]
        e3 = n * n * r["log_lam1"]
        print(f"{n:>3} {r['lam1']:>10.4f} {r['count']:>12} {lc:>12.4f} "
              f"{e1:>14.4f} {e2:>14.4f} {e3:>14.4f}")

    # Calibrar usando o n maior com count exato
    cal_rows = [r for r in rows if r["count"] is not None]
    cal = max(cal_rows, key=lambda r: r["n"])
    n_cal = cal["n"]
    log_count_cal = np.log(cal["count"])

    C_A_log = log_count_cal - (n_cal - 1) * cal["log_lam1"]
    C_B_log = log_count_cal - n_cal * cal["log_lam1"]
    C_C_log = log_count_cal - n_cal * n_cal * cal["log_lam1"]

    print(f"\nCalibracao com n={n_cal}:")
    print(f"  A: log C = log count - (n-1) log lam1 = {C_A_log:.4f}   C = {np.exp(C_A_log):.4e}")
    print(f"  B: log C = log count -  n  log lam1 = {C_B_log:.4f}   C = {np.exp(C_B_log):.4e}")
    print(f"  C: log C = log count - n^2 log lam1 = {C_C_log:.4f}   C = {np.exp(C_C_log):.4e}")

    print("\n=== Extrapolacao ===")
    print(f"{'n':>3} {'log_pred_A':>12} {'log_pred_B':>12} {'log_pred_C':>12} {'log_exato':>12}")
    for r in rows:
        n = r["n"]
        lp_A = C_A_log + (n - 1) * r["log_lam1"]
        lp_B = C_B_log + n * r["log_lam1"]
        lp_C = C_C_log + n * n * r["log_lam1"]
        lc = np.log(r["count"]) if r["count"] else None
        print(f"{n:>3} {lp_A:>12.4f} {lp_B:>12.4f} {lp_C:>12.4f} "
              f"{lc if lc is not None else 'TBD':>12}")
        r["log_pred_A"] = lp_A
        r["log_pred_B"] = lp_B
        r["log_pred_C"] = lp_C
        r["log_count_exato"] = float(lc) if lc is not None else None

    out = {
        "calibration_n": n_cal,
        "log_C_A": float(C_A_log),  # for count ~ C * lam1^(n-1)
        "log_C_B": float(C_B_log),  # for count ~ C * lam1^n
        "log_C_C": float(C_C_log),  # for count ~ C * lam1^(n^2)
        "rows": rows,
    }
    with open(DATA / "tour_count_estimates.json", "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nSalvo em data/tour_count_estimates.json")


if __name__ == "__main__":
    ns = [int(x) for x in sys.argv[1:]] if len(sys.argv) > 1 else [6]
    main(ns)
