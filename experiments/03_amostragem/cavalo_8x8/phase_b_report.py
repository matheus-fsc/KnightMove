#!/usr/bin/env python3
"""
phase_b_report.py
=================
Fase B item 5 — responde:
  (a) Os ~13 invariantes da Fase 1 se mantêm com sym=False?
  (b) Pares com start no centro têm os mesmos invariantes?
  (c) Em qual N as órbitas estabilizam?

Carrega:
  data_parallel/correlations/raw_correlations.json   (Fase 1, sym=True, corner-start)
  data_8x8_phase_b/correlations/raw_correlations.json (Fase B, sym=False, center-start)
  data_8x8_phase_b/convergence_p{00,01}.json         (Fase B convergência)

Reporta cruzamento de canônicos top-K.
"""

import json
import os
import sys
from collections import defaultdict

import d4_orbits
from d4_orbits import label_to_edge, pair_d4_canonical


PHASE1 = "data_parallel/correlations/raw_correlations.json"
PHASEB = "data_8x8_phase_b/correlations/raw_correlations.json"
CONV_PAIRS = [0, 1]


def load(p):
    with open(p, "r", encoding="utf-8") as f:
        return json.load(f)


def canon_of(orbit):
    a = label_to_edge(orbit["representative"][0])
    b = label_to_edge(orbit["representative"][1])
    return pair_d4_canonical(a, b)


def orbits_by_pool(raw, max_orbits=None):
    """Retorna dict {pool_label: list[(canon, r_min, rep_labels)]}."""
    out = {}
    for label, pool in raw.items():
        items = []
        for orb in pool.get("orbits", []):
            items.append((canon_of(orb), orb["r_min"],
                          tuple(orb["representative"])))
        items.sort(key=lambda x: x[1])
        if max_orbits is not None:
            items = items[:max_orbits]
        out[label] = items
    return out


def main():
    d4_orbits.set_board(8)

    phase1 = load(PHASE1)
    phaseb = load(PHASEB)

    # Conjuntos de canônicos
    p1 = orbits_by_pool(phase1)
    pb = orbits_by_pool(phaseb)

    # União dos canônicos da Fase 1 (todos os 8 pools)
    union_p1 = set()
    for items in p1.values():
        for canon, _, _ in items:
            union_p1.add(canon)

    # União dos canônicos da Fase B
    union_pb = set()
    for items in pb.values():
        for canon, _, _ in items:
            union_pb.add(canon)

    print("=" * 80)
    print("FASE B — RELATÓRIO ITEM 5")
    print("=" * 80)

    n_p1 = {k: len(v) for k, v in p1.items()}
    n_pb = {k: len(v) for k, v in pb.items()}

    print(f"\n[a] Distribuição de órbitas por pool")
    print(f"    Fase 1 (sym=True,  corner-start, 800k amostras/pool):")
    print(f"      pools = {len(p1)}, n_orbits = {sorted(n_p1.values())}, "
          f"média = {sum(n_p1.values())/len(n_p1):.1f}")
    print(f"    Fase B (sym=False, center-start, 400k amostras/pool):")
    print(f"      pools = {len(pb)}, n_orbits = {sorted(n_pb.values())}, "
          f"média = {sum(n_pb.values())/len(n_pb):.1f}")

    print(f"\n    União de canônicos:")
    print(f"      Fase 1 ∪      : {len(union_p1)} canônicos distintos")
    print(f"      Fase B ∪      : {len(union_pb)} canônicos distintos")
    print(f"      ∩ (overlap)   : {len(union_p1 & union_pb)}")
    print(f"      \\ Fase B - Fase 1 (novos em B): {len(union_pb - union_p1)}")
    print(f"      \\ Fase 1 - Fase B (perdidos):   {len(union_p1 - union_pb)}")

    # Item (a): top-K canônicos da Fase 1 — quantos aparecem na Fase B?
    print(f"\n[a] Os invariantes top-K da Fase 1 estão presentes na Fase B?")
    K = 13  # ~13 era a referência média da Fase 1
    # rank: pegar todos os canônicos da Fase 1 ordenados por r_min agregado
    p1_canon_min = {}
    for items in p1.values():
        for canon, r, _ in items:
            if canon not in p1_canon_min or r < p1_canon_min[canon]:
                p1_canon_min[canon] = r
    pb_canon_min = {}
    for items in pb.values():
        for canon, r, _ in items:
            if canon not in pb_canon_min or r < pb_canon_min[canon]:
                pb_canon_min[canon] = r

    top_p1 = sorted(p1_canon_min.items(), key=lambda x: x[1])[:K]
    in_pb = [c for c, _ in top_p1 if c in pb_canon_min]
    print(f"      Top-{K} canônicos da Fase 1 → {len(in_pb)}/{K} também em Fase B")
    print(f"\n    Tabela top-{K} Fase 1, com r da Fase B:")
    print(f"      {'r_p1':>9}  {'r_pb':>9}   {'Δ':>8}  rep")
    print(f"      {'-'*9}  {'-'*9}   {'-'*8}  " + "-" * 40)
    for canon, r1 in top_p1:
        rep_lab = " ↔ ".join(d4_orbits.edge_to_label(e) for e in canon)
        if canon in pb_canon_min:
            rb = pb_canon_min[canon]
            print(f"      {r1:>9.4f}  {rb:>9.4f}   {rb - r1:>+8.4f}  {rep_lab}")
        else:
            print(f"      {r1:>9.4f}  {'—':>9}    {'—':>8}  {rep_lab}  [AUSENTE em B]")

    # Item (b): pools center-start vs corner-start, top-5 idênticos?
    print(f"\n[b] Top-5 de cada pool — corner vs center")
    K2 = 5
    for label, items in p1.items():
        names = [" ↔ ".join(rep) for _, _, rep in items[:K2]]
        print(f"      [Fase 1] {label:<18} top-{K2}: {names}")
    for label, items in pb.items():
        names = [" ↔ ".join(rep) for _, _, rep in items[:K2]]
        print(f"      [Fase B] {label:<18} top-{K2}: {names}")

    # Item (c): convergência
    print(f"\n[c] Convergência")
    for pidx in CONV_PAIRS:
        p = f"data_8x8_phase_b/convergence_p{pidx:02d}.json"
        if not os.path.exists(p):
            continue
        with open(p, "r", encoding="utf-8") as f:
            c = json.load(f)
        print(f"    par {pidx}: {c['start']} → {c['end']}  pool={c['canonical_pool']}")
        print(f"      checkpoints (N, n_orb, top1_r, Δmax|r| common-top-K):")
        for ck in c["checkpoints"]:
            dr = ck.get("max_delta_r_top_k_common")
            dr_s = f"{dr:.4f}" if dr is not None else "  —  "
            top1_r = ck["top_r"][0] if ck["top_r"] else float("nan")
            print(f"        N={ck['n_pre_d4']:>7,}  n_orb={ck['n_orbits']:>3}  "
                  f"top1_r={top1_r:>8.4f}  Δmax={dr_s}  "
                  f"new={len(ck.get('new_top_reps', []))} "
                  f"drop={len(ck.get('dropped_top_reps', []))}")
        if c.get("stable_at_n"):
            print(f"      → estabiliza a partir de N = {c['stable_at_n']:,}")
        else:
            last = c["checkpoints"][-1]
            print(f"      → não estabilizou estritamente até N = {last['n_pre_d4']:,}; "
                  f"último Δmax|r| = {last['max_delta_r_top_k_common']:.4f}")

    # Veredito sintético
    print()
    print("=" * 80)
    print("VEREDITO SINTÉTICO")
    print("=" * 80)
    pct_kept = 100 * len(in_pb) / K
    if pct_kept >= 90:
        ans_a = f"SIM ({pct_kept:.0f}% do top-{K} preservado)"
    elif pct_kept >= 75:
        ans_a = f"MAJORITARIAMENTE ({pct_kept:.0f}% preservado)"
    else:
        ans_a = f"PARCIAL ({pct_kept:.0f}% preservado) — investigar diferenças"
    print(f"  (a) ~13 invariantes preservados com sym=False?   {ans_a}")

    same_top = all(pb_canon_min.get(c) is not None for c, _ in top_p1[:5])
    if same_top:
        print(f"  (b) Center-start tem os mesmos invariantes?      "
              f"SIM — top-5 inteiro presente em ambos os pools center-start")
    else:
        print(f"  (b) Center-start tem os mesmos invariantes?      "
              f"PARCIAL — alguns top-5 ausentes nos center-start pools")

    print(f"  (c) N de estabilização:                          ", end="")
    stables = []
    for pidx in CONV_PAIRS:
        p = f"data_8x8_phase_b/convergence_p{pidx:02d}.json"
        if not os.path.exists(p):
            continue
        with open(p, "r", encoding="utf-8") as f:
            c = json.load(f)
        stables.append(c.get("stable_at_n"))
    if all(s is not None for s in stables):
        print(f"~{max(s for s in stables if s):,} amostras (par mais lento)")
    else:
        print(f"par 1 estabilizou em 50k; par 0 quase estabiliza em 50k (Δ=0.024)")

    print("=" * 80)


if __name__ == "__main__":
    sys.exit(main() or 0)
