#!/usr/bin/env python3
"""
report.py
=========
Tabela comparativa Z3 vs Backtracking 6×6 a partir dos JSONs em results/.

Lê:
  baseline_bt_6x6.json          — BT instrumentado
  benchmark_z3_aggregated.json  — Z3 agregado (80k amostras existentes)
  benchmark_z3_fresh.json       — Z3 fresco com K controlado
  benchmark_z3_fresh_invariants.json — H4
  kl_divergence_6x6.json        — viés de amostragem

Escreve:
  benchmark/results/report.txt
"""

import json
import sys
from pathlib import Path
from statistics import mean, stdev

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "benchmark" / "results"
OUT = RESULTS / "report.txt"


def load(name):
    p = RESULTS / name
    if not p.exists():
        return None
    return json.loads(p.read_text())


def fmt(v, w=10, prec=3):
    if v is None:
        return "—".rjust(w)
    if isinstance(v, float):
        return f"{v:>{w}.{prec}f}"
    if isinstance(v, int):
        return f"{v:>{w},}"
    return str(v).rjust(w)


def main():
    lines = []
    w = lines.append

    bt = load("baseline_bt_6x6.json")
    z3_agg = load("benchmark_z3_aggregated.json")
    z3_fresh = load("benchmark_z3_fresh.json")
    z3_inv = load("benchmark_z3_fresh_invariants.json")
    kl = load("kl_divergence_6x6.json")

    w("=" * 78)
    w("BENCHMARK 6×6  —  Z3 GF(2)  vs  BACKTRACKING")
    w("=" * 78)
    w("")

    # ── Baseline BT ───────────────────────────────────────────────────────
    if bt is not None:
        w("── BASELINE BACKTRACKING (busca a partir de v=0 = A6) ─────────────")
        w(f"  Ciclos dirigidos em v=0      : {bt['total_directed_cycles_at_v0']:,}")
        w(f"  Total c/ rotação (×36)       : {bt['total_solutions_with_rotation']:,}")
        w("")
        w(f"  {'método':22s} {'target':>8s} {'sols':>8s} {'tempo(s)':>10s} "
          f"{'sols/s':>10s} {'de/sol':>10s} {'nodes/sol':>12s} "
          f"{'t_1ª(s)':>10s} {'mem(MB)':>10s}")
        w("  " + "─" * 110)
        for r in bt["runs"]:
            w(f"  {r['method']:22s} {fmt(r['target_k'], 8, 0)} "
              f"{fmt(r['n_solutions'], 8, 0)} {fmt(r['elapsed_s'], 10, 2)} "
              f"{fmt(r['solutions_per_second'], 10, 2)} "
              f"{fmt(r['dead_ends_per_solution'], 10, 3)} "
              f"{fmt(r['nodes_per_solution'], 12, 3)} "
              f"{fmt(r['time_to_first_s'], 10, 4)} "
              f"{fmt(r['peak_memory_mb'], 10, 1)}")
        # curva time-to-K se disponível
        for r in bt["runs"]:
            if r.get("time_to_k_s"):
                w("")
                w(f"  curva time-to-K  ({r['method']}, target={r['target_k']}):")
                for k, t in r["time_to_k_s"].items():
                    w(f"    K={k:>6s} → {fmt(t, 10, 4)}s")
        w("")

    # ── Z3 agregado ───────────────────────────────────────────────────────
    if z3_agg is not None:
        t = z3_agg["totals"]
        w("── Z3 GF(2)  AGREGADO (dados existentes em data_6x6_sampled) ───────")
        w(f"  fonte                : {z3_agg['source']}")
        w(f"  amostras totais      : {t['n_samples']:,}")
        w(f"  attempts totais      : {t['n_attempts']:,}")
        w(f"  attempts/sol         : {t['attempts_per_solution']}")
        w(f"  failed/sol           : {t['failed_per_solution']}")
        w(f"  tempo agregado (s)   : {t['elapsed_s']:.1f}")
        w(f"  sols/s (1 worker)    : {t['solutions_per_second']}")
        w("")
        w("  Por par canônico:")
        w(f"    {'(start)':>10s} {'(end)':>10s} {'amostras':>10s} "
          f"{'att/sol':>10s} {'fail/sol':>10s} {'s/sol':>10s}")
        for key, e in z3_agg["by_pair"].items():
            w(f"    {str(tuple(e['start'])):>10s} {str(tuple(e['end'])):>10s} "
              f"{fmt(e['n_samples'], 10, 0)} "
              f"{fmt(e['attempts_per_solution'], 10, 3)} "
              f"{fmt(e['failed_per_solution'], 10, 3)} "
              f"{fmt(e['s_per_solution'], 10, 4)}")
        w("")

    # ── Z3 fresco vs Z3 com invariantes ───────────────────────────────────
    def emit_fresh(data, header):
        if data is None:
            return
        w(f"── {header} ─────────────────────────────────────────")
        w(f"  {'par':>20s} {'K':>6s} {'found':>6s} "
          f"{'tempo(s)':>10s} {'att/sol':>10s} {'fail/sol':>10s} "
          f"{'t_1ª(s)':>10s} {'t_K(s)':>10s} {'mem(MB)':>10s}")
        w("  " + "─" * 110)
        rows = []
        for pair_block in data["pairs"]:
            for r in pair_block["runs"]:
                pair_str = f"{tuple(pair_block['start'])}→{tuple(pair_block['end'])}"
                t_k = r.get("time_to_k_s", {}).get(str(r["target_k"]))
                w(f"  {pair_str:>20s} {fmt(r['target_k'], 6, 0)} "
                  f"{fmt(r['found_k'], 6, 0)} "
                  f"{fmt(r['elapsed_s'], 10, 2)} "
                  f"{fmt(r['attempts_per_solution'], 10, 3)} "
                  f"{fmt(r['failed_per_solution'], 10, 3)} "
                  f"{fmt(r['time_to_first_s'], 10, 4)} "
                  f"{fmt(t_k, 10, 4)} "
                  f"{fmt(r['peak_memory_mb'], 10, 1)}")
                rows.append(r)
        # médias por K
        if rows:
            w("")
            w("  Médias por K:")
            ks = sorted({r["target_k"] for r in rows})
            for k in ks:
                rs = [r for r in rows if r["target_k"] == k and r["found_k"] >= 1]
                if not rs:
                    continue
                mean_att = mean(r["attempts_per_solution"] for r in rs)
                mean_t1 = mean(r["time_to_first_s"] for r in rs if r["time_to_first_s"])
                mean_total = mean(r["elapsed_s"] for r in rs)
                w(f"    K={k:>5}  att/sol={mean_att:.3f}  "
                  f"t_1ª={mean_t1:.4f}s  tempo_total={mean_total:.2f}s "
                  f"({len(rs)} pares)")
        w("")
        return rows

    rows_plain = emit_fresh(z3_fresh, "Z3 GF(2)  FRESCO  (sym=False, sem invariantes)")
    rows_inv = emit_fresh(z3_inv, "Z3 + INVARIANTES  (H4: NOT(A AND B) para r<-0.5)")

    # H4 comparison
    if rows_plain and rows_inv:
        w("── H4: comparação direta (mesmos pares, mesmo K) ───────────────────")
        plain_idx = {(tuple(r['start']), tuple(r['end']), r['target_k']): r for r in rows_plain}
        for r in rows_inv:
            key = (tuple(r['start']), tuple(r['end']), r['target_k'])
            base = plain_idx.get(key)
            if base is None or not base.get("found_k"):
                continue
            d_att = (r["attempts_per_solution"] - base["attempts_per_solution"])
            d_time = r["elapsed_s"] - base["elapsed_s"]
            pct_att = (d_att / base["attempts_per_solution"]) * 100
            pct_time = (d_time / base["elapsed_s"]) * 100
            w(f"  par {key[0]}→{key[1]}  K={key[2]}: "
              f"att/sol {base['attempts_per_solution']:.3f}→{r['attempts_per_solution']:.3f} "
              f"({pct_att:+.1f}%)   "
              f"tempo {base['elapsed_s']:.2f}→{r['elapsed_s']:.2f}s ({pct_time:+.1f}%)")
        w("")

    # ── KL-divergência ────────────────────────────────────────────────────
    if kl is not None:
        w("── KL-DIVERGÊNCIA  Z3 vs uniforme sobre GT (6×6) ───────────────────")
        w(f"  {'método':25s} {'par':>20s} {'K':>6s} "
          f"{'KL(bits)':>10s} {'cobertura':>10s} {'fora_sup':>10s}")
        w("  " + "─" * 88)
        for r in kl["results"]:
            pair = f"{tuple(r['start'])}→{tuple(r['end'])}"
            w(f"  {r['method']:25s} {pair:>20s} {fmt(r['target_k'], 6, 0)} "
              f"{fmt(r['kl_partial_bits'], 10, 4)} "
              f"{fmt(r['coverage_gt'] * 100 if r['coverage_gt'] else None, 10, 2)} "
              f"{fmt(r['out_support'], 10, 0)}")
        w("")
        w("  KL=0  : amostragem uniforme. KL>0 : viés (KL alto = pesos não-uniformes).")
        w("  cobertura: fração das assinaturas GT distintas que o Z3 visitou.")
        w("")

    w("=" * 78)
    OUT.write_text("\n".join(lines))
    print("\n".join(lines))
    print(f"\nRelatório salvo em {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
