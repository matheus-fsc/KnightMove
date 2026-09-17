#!/usr/bin/env python3
"""
compare_validation.py
=====================
Compara lado-a-lado dois relatórios de validate_sampler.py.
Uso típico: Fase A4 — quantificar viés do break_symmetry.

Uso:
  python compare_validation.py data_6x6_sampled data_6x6_sym
"""

import argparse
import json
import os
import sys


def load(d):
    p = os.path.join(d, "validation.json")
    with open(p, "r", encoding="utf-8") as f:
        return json.load(f)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dir_a", help="diretório A (ex: data_6x6_sampled — sym=False)")
    ap.add_argument("dir_b", help="diretório B (ex: data_6x6_sym — sym=True)")
    ap.add_argument("--label-a", default=None)
    ap.add_argument("--label-b", default=None)
    args = ap.parse_args()

    a = load(args.dir_a)
    b = load(args.dir_b)
    la = args.label_a or os.path.basename(os.path.normpath(args.dir_a))
    lb = args.label_b or os.path.basename(os.path.normpath(args.dir_b))

    print("=" * 84)
    print(f"COMPARAÇÃO  A = {la}   B = {lb}")
    print("=" * 84)

    rows = [
        ("orbits_recovered",      f"{a['orbits_recovered']}/{a['gt_n_orbits']}", f"{b['orbits_recovered']}/{b['gt_n_orbits']}"),
        ("ranking_correlation",   f"{a['ranking_correlation']:.4f}",              f"{b['ranking_correlation']:.4f}"),
        ("r_mae",                 f"{a['r_mae']:.4f}",                            f"{b['r_mae']:.4f}"),
        ("z3_n_orbits_aggregated",str(a['z3_n_orbits_aggregated']),               str(b['z3_n_orbits_aggregated'])),
        ("mandatory_subset_ok",   str(a['mandatory']['subset_match']),            str(b['mandatory']['subset_match'])),
        ("verdict",               a['verdict'],                                    b['verdict']),
    ]
    print(f"  {'métrica':<24} {'A':>18} {'B':>18}   delta")
    print(f"  {'-'*24} {'-'*18} {'-'*18}   {'-'*12}")
    for name, va, vb in rows:
        try:
            fa = float(va); fb = float(vb)
            delta = f"{fb - fa:+.4f}"
        except ValueError:
            delta = ""
        print(f"  {name:<24} {va:>18} {vb:>18}   {delta}")

    # Tabela órbita a órbita
    print()
    print("Órbitas matched, lado-a-lado:")
    print(f"  {'representante':<32} {'r_gt':>8} {'r_z3(A)':>9} {'Δ(A)':>8} {'r_z3(B)':>9} {'Δ(B)':>8}")
    print(f"  {'-'*32} {'-'*8} {'-'*9} {'-'*8} {'-'*9} {'-'*8}")
    # Indexar matched por representativo
    def key_of(m):
        return tuple(m["representative"])
    a_by = {key_of(m): m for m in a["matched_orbits"]}
    b_by = {key_of(m): m for m in b["matched_orbits"]}
    all_keys = sorted(set(a_by) | set(b_by), key=lambda k: a_by.get(k, b_by[k])["r_gt"])
    for k in all_keys:
        ma = a_by.get(k); mb = b_by.get(k)
        rep = " ↔ ".join(k)
        r_gt = (ma or mb)["r_gt"]
        ra = f"{ma['r_z3']:.4f}" if ma else "  —  "
        da = f"{ma['delta_r']:+.4f}" if ma else "  —  "
        rb = f"{mb['r_z3']:.4f}" if mb else "  —  "
        db = f"{mb['delta_r']:+.4f}" if mb else "  —  "
        print(f"  {rep:<32} {r_gt:>8.4f} {ra:>9} {da:>8} {rb:>9} {db:>8}")

    print()
    print("Interpretação:")
    a_pass = a["verdict"] == "VALID"
    b_pass = b["verdict"] == "VALID"
    if a_pass and not b_pass:
        print(f"  → {la}: método APROVADO; {lb}: método REPROVADO/PARCIAL.")
        print(f"  → A flag que difere os dois introduz viés mensurável.")
    elif b_pass and not a_pass:
        print(f"  → {lb}: método APROVADO; {la}: método REPROVADO/PARCIAL.")
    elif a_pass and b_pass:
        print(f"  → Ambos aprovados. Diferenças são marginais.")
    else:
        print(f"  → Ambos reprovados/parciais. Reavaliar metodologia.")

    # Delta nas métricas-chave
    drho = b["ranking_correlation"] - a["ranking_correlation"]
    dmae = b["r_mae"] - a["r_mae"]
    print(f"  → Δρ = {drho:+.4f}  (negativo = B degrada ranking vs A)")
    print(f"  → Δmae = {dmae:+.4f}  (positivo = B aumenta erro absoluto vs A)")
    print()


if __name__ == "__main__":
    sys.exit(main() or 0)
