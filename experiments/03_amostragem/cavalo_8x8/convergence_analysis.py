#!/usr/bin/env python3
"""
convergence_analysis.py
=======================
Fase B item 4 — análise de convergência dos invariantes em função do
número de amostras.

Para cada par canônico (pool), recomputa as correlações Z3 em checkpoints
crescentes (ex: 1k, 5k, 10k, 20k, 50k amostras *pré*-D₄), aplica expansão
D₄ se solicitado, e reporta:
  - n_orbits     : quantas órbitas D₄ acima do threshold
  - top_k_r      : r dos top-K invariantes
  - top_k_reps   : representantes canônicos dos top-K (rótulos)
  - stable_set   : conjunto de representantes canônicos do top-K
  - delta_top_k  : delta de r máximo entre este checkpoint e o anterior

Critério prático de estabilização (heurística): top-K reps estabilizam
e máximo |Δr_k| < tol entre dois checkpoints consecutivos.

Uso:
  python convergence_analysis.py \\
    --data-dir data_8x8_phase_b \\
    --pair-idx 0 \\
    --checkpoints 1000,5000,10000,20000,50000 \\
    --top-k 8 \\
    --board 8 \\
    --expand-d4
"""

import argparse
import glob
import json
import os
import re
import sys
from collections import defaultdict

import numpy as np

import d4_orbits
from d4_orbits import label_to_edge, edge_to_label, pair_d4_canonical
from correlator import (
    _build_d4_remaps, pearson_negative_pairs, collapse_to_orbits,
)
from expand_d4 import canonical_pair


_BATCH_RE = re.compile(r"batch_p(\d+)_b(\d+)\.json")


def find_batches(samples_dir, pair_idx):
    paths = sorted(glob.glob(os.path.join(samples_dir, "batch_*.json")))
    out = []
    for p in paths:
        m = _BATCH_RE.search(os.path.basename(p))
        if m and int(m.group(1)) == pair_idx:
            out.append((int(m.group(2)), p))
    out.sort()
    return [p for _, p in out]


def accumulate_signatures(batch_paths, n_target, edges_ref=None,
                           remaps=None, expand=False):
    """Lê batches sequencialmente até alcançar n_target amostras pré-D₄."""
    parts = []
    n = 0
    edges = edges_ref
    for path in batch_paths:
        if n >= n_target:
            break
        with open(path, "r", encoding="utf-8") as f:
            b = json.load(f)
        if edges is None:
            edges = b["edges"]
        elif b["edges"] != edges:
            raise ValueError(f"ordem de arestas inconsistente em {path}")
        X = np.asarray(b["signatures"], dtype=np.bool_)
        take = min(X.shape[0], n_target - n)
        Xtake = X[:take]
        if expand and remaps is not None:
            for t in range(8):
                new_X = np.empty_like(Xtake)
                new_X[:, remaps[t]] = Xtake
                parts.append(new_X)
        else:
            parts.append(Xtake)
        n += take
    if not parts:
        return edges, np.zeros((0, len(edges) if edges else 0), dtype=np.bool_)
    return edges, np.concatenate(parts, axis=0)


def analyze_at_n(batch_paths, n_target, threshold, top_k, expand_d4):
    """Computa órbitas D₄ com as primeiras n_target amostras pré-D₄."""
    # primeira leitura para descobrir as arestas
    with open(batch_paths[0], "r", encoding="utf-8") as f:
        first = json.load(f)
    edges = first["edges"]
    remaps = _build_d4_remaps(edges) if expand_d4 else None

    edges, X = accumulate_signatures(
        batch_paths, n_target, edges_ref=edges, remaps=remaps, expand=expand_d4
    )
    pairs, _ = pearson_negative_pairs(X, threshold)
    orbits = collapse_to_orbits(pairs, edges)
    orbits.sort(key=lambda o: o["r_min"])

    top = orbits[:top_k]
    top_reps = [tuple(o["representative"]) for o in top]
    top_r = [o["r_min"] for o in top]

    return {
        "n_pre_d4": int(n_target),
        "n_pos_d4": int(X.shape[0]),
        "n_orbits": len(orbits),
        "n_correlations": len(pairs),
        "top_k": top_k,
        "top_reps": top_reps,
        "top_r": top_r,
        "all_orbit_reps": [tuple(o["representative"]) for o in orbits],
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", required=True,
                    help="diretório com samples/ (ex: data_8x8_phase_b)")
    ap.add_argument("--pair-idx", type=int, required=True,
                    help="índice do par no config (ex: 0)")
    ap.add_argument("--checkpoints", default="1000,5000,10000,20000,50000",
                    help="lista de N pré-D₄ separados por vírgula")
    ap.add_argument("--threshold", type=float, default=-0.30)
    ap.add_argument("--top-k", type=int, default=8)
    ap.add_argument("--board", type=int, default=8)
    ap.add_argument("--expand-d4", action="store_true")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    d4_orbits.set_board(args.board)

    samples_dir = os.path.join(args.data_dir, "samples")
    batch_paths = find_batches(samples_dir, args.pair_idx)
    if not batch_paths:
        print(f"Nenhum batch encontrado para pair_idx={args.pair_idx} em {samples_dir}")
        return 1

    # descobre quantas amostras temos no total
    with open(batch_paths[0], "r", encoding="utf-8") as f:
        first = json.load(f)
    start = tuple(first["start"])
    end = tuple(first["end"])

    # contagem rápida do total (lê n_found de cada batch)
    n_total = 0
    for p in batch_paths:
        with open(p, "r", encoding="utf-8") as f:
            b = json.load(f)
        n_total += b.get("n_found", len(b.get("signatures", [])))

    checkpoints = sorted(set(int(x) for x in args.checkpoints.split(",")))
    checkpoints = [c for c in checkpoints if c <= n_total]
    if not checkpoints:
        print(f"Todos os checkpoints > total ({n_total} amostras). Reduzindo.")
        checkpoints = [n_total]

    canon = canonical_pair(start, end)
    canon_label = (f"{canon[0][0]}{canon[0][1]}_{canon[1][0]}{canon[1][1]}_canon")

    print("=" * 72)
    print(f"CONVERGENCE ANALYSIS — par {args.pair_idx}  ({start} → {end})")
    print(f"  canonical pool   : {canon_label}")
    print(f"  batches no disco : {len(batch_paths)}")
    print(f"  amostras totais  : {n_total:,}  (pré-D₄)")
    print(f"  top_k            : {args.top_k}")
    print(f"  threshold        : {args.threshold}")
    print(f"  expand_d4        : {args.expand_d4}")
    print(f"  checkpoints      : {checkpoints}")
    print("=" * 72)

    reports = []
    prev_set = None
    prev_top_r = None
    for n in checkpoints:
        r = analyze_at_n(batch_paths, n, args.threshold, args.top_k, args.expand_d4)
        rep_set = set(tuple(rep) for rep in r["top_reps"])

        # Delta dos top-K r-values (assumindo mesma ordem por representante)
        max_dr = None
        new_reps = []
        dropped_reps = []
        if prev_set is not None:
            new_reps = sorted(rep_set - prev_set)
            dropped_reps = sorted(prev_set - rep_set)
            # mapa rep -> r_min
            cur = dict(zip(r["top_reps"], r["top_r"]))
            prev = dict(zip(prev_top_r["reps"], prev_top_r["rs"]))
            common = rep_set & prev_set
            if common:
                drs = [abs(cur[c] - prev[c]) for c in common]
                max_dr = max(drs)

        r["new_top_reps"] = new_reps
        r["dropped_top_reps"] = dropped_reps
        r["max_delta_r_top_k_common"] = max_dr
        reports.append(r)

        prev_set = rep_set
        prev_top_r = {"reps": r["top_reps"], "rs": r["top_r"]}

    # Print tabela
    print(f"\n  {'N':>7}  {'n_orb':>5}  {'n_corr':>6}  {'top1_r':>8}  "
          f"{'top1_rep':<30}  {'Δmax|r|':>8}  {'new':>4} {'drop':>4}")
    print(f"  {'-'*7}  {'-'*5}  {'-'*6}  {'-'*8}  {'-'*30}  {'-'*8}  {'-'*4} {'-'*4}")
    for rep in reports:
        top1_r = rep["top_r"][0] if rep["top_r"] else float("nan")
        top1_rep = " ↔ ".join(rep["top_reps"][0]) if rep["top_reps"] else "—"
        dr = rep["max_delta_r_top_k_common"]
        dr_s = f"{dr:.4f}" if dr is not None else "  —  "
        print(f"  {rep['n_pre_d4']:>7,}  {rep['n_orbits']:>5}  "
              f"{rep['n_correlations']:>6}  {top1_r:>8.4f}  {top1_rep:<30}  "
              f"{dr_s:>8}  {len(rep['new_top_reps']):>4} {len(rep['dropped_top_reps']):>4}")

    # Veredito de estabilização
    print(f"\n  TOP-{args.top_k} representantes no maior N:")
    last = reports[-1]
    for i, (rep, r) in enumerate(zip(last["top_reps"], last["top_r"])):
        print(f"    {i+1:>2}. r={r:>8.4f}  {' ↔ '.join(rep)}")

    print(f"\n  Estabilização (heurística):")
    stable_idx = None
    TOL_DR = 0.02
    for i in range(1, len(reports)):
        rep = reports[i]
        if rep["max_delta_r_top_k_common"] is None:
            continue
        if (not rep["new_top_reps"] and not rep["dropped_top_reps"]
                and rep["max_delta_r_top_k_common"] < TOL_DR):
            stable_idx = i
            break
    if stable_idx is not None:
        print(f"    estabiliza a partir de N = {reports[stable_idx]['n_pre_d4']:,}  "
              f"(Δmax|r| < {TOL_DR} e top-{args.top_k} reps idênticos ao anterior)")
    else:
        print(f"    NÃO atingiu estabilização com top-{args.top_k} reps idênticos "
              f"até N={reports[-1]['n_pre_d4']:,}")
        print(f"    (último Δmax|r|={reports[-1]['max_delta_r_top_k_common']})")

    out_path = args.out or os.path.join(
        args.data_dir, f"convergence_p{args.pair_idx:02d}.json"
    )
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump({
            "pair_idx": args.pair_idx,
            "start": list(start),
            "end": list(end),
            "canonical_pool": canon_label,
            "n_batches": len(batch_paths),
            "n_total_pre_d4": n_total,
            "threshold": args.threshold,
            "top_k": args.top_k,
            "expand_d4": args.expand_d4,
            "checkpoints": reports,
            "stable_at_n": (reports[stable_idx]["n_pre_d4"]
                             if stable_idx is not None else None),
        }, f, indent=2, ensure_ascii=False, default=str)
    print(f"\n  Relatório JSON salvo em: {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
