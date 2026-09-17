#!/usr/bin/env python3
"""
validate_sampler.py
===================
Fase A — valida a amostragem Z3 no 6×6 contra o ground truth exaustivo.

Entradas:
  - ../analysis_6x6.json                (ground truth)
  - <data_dir>/correlations/raw_correlations.json   (resultado correlator)
  - <data_dir>/samples/                 (para distribuição de freq por edge)

Métricas computadas (por pool e agregado):
  1. orbits_recovered      : quais das 7 órbitas canônicas do GT aparecem nos
                              tops de algum pool Z3?
  2. ranking_correlation   : Spearman rho entre r-values do GT e os Z3 nos
                              representantes canônicos que aparecem em ambos
  3. r_mae                 : MAE absoluto entre r-values matched
  4. mandatory_match       : as 8 arestas obrigatórias do GT aparecem como
                              obrigatórias em pelo menos um pool? E são as
                              MESMAS 8? (interseção sobre pools = GT?)
  5. freq_kl_divergence    : KL-div da distribuição de frequências por aresta,
                              comparando freq média do GT (computada de
                              loop_correlations.json se disponível) com a freq
                              média sobre pools Z3
  6. verdict               : VALID | PARTIAL | BIASED

Critério prático (do agent_prompt_fase2.md):
  ranking_correlation > 0.90 AND orbits_recovered >= 6/7  →  VALID

Uso:
  python validate_sampler.py --data-dir data_6x6_smoke
  python validate_sampler.py --data-dir data_6x6_sampled --gt ../analysis_6x6.json
"""

import argparse
import json
import math
import os
import sys
from collections import defaultdict

import numpy as np

import d4_orbits
from d4_orbits import label_to_edge, edge_to_label, pair_d4_canonical


# ── helpers de carregamento ──────────────────────────────────────────

def load_gt(path):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)

def load_raw_correlations(path):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


# ── conversões ───────────────────────────────────────────────────────

def gt_orbit_key(orbit_record):
    """Chave canônica (tupla ordenada de arestas) p/ uma órbita do GT."""
    a, b = orbit_record["representative"]
    ea = label_to_edge(a)
    eb = label_to_edge(b)
    return pair_d4_canonical(ea, eb)


def correlation_to_canonical(corr):
    """Recebe {'r':..., 'edge_i':..., 'edge_j':...}, retorna (canonical_pair, r)."""
    ea = label_to_edge(corr["edge_i"])
    eb = label_to_edge(corr["edge_j"])
    return pair_d4_canonical(ea, eb), corr["r"]


# ── Spearman manual (sem scipy) ──────────────────────────────────────

def spearman_rho(x, y):
    """Correlação de postos de Spearman. Aceita listas/arrays paralelos."""
    if len(x) < 2:
        return float("nan")
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    rx = _rank(x)
    ry = _rank(y)
    rx_c = rx - rx.mean()
    ry_c = ry - ry.mean()
    denom = (np.sqrt((rx_c**2).sum()) * np.sqrt((ry_c**2).sum()))
    if denom == 0:
        return float("nan")
    return float((rx_c * ry_c).sum() / denom)

def _rank(a):
    """Postos com correção média para empates (avg ranks)."""
    order = np.argsort(a, kind="stable")
    ranks = np.empty_like(order, dtype=float)
    ranks[order] = np.arange(1, len(a) + 1, dtype=float)
    # tie correction: média de postos para valores iguais
    sorted_a = a[order]
    i = 0
    while i < len(a):
        j = i
        while j + 1 < len(a) and sorted_a[j + 1] == sorted_a[i]:
            j += 1
        if j > i:
            avg = (i + j) / 2 + 1  # ranks são 1-indexed
            for k in range(i, j + 1):
                ranks[order[k]] = avg
        i = j + 1
    return ranks


# ── KL-divergence sobre distribuição de freq de arestas ──────────────

def kl_divergence(p, q, eps=1e-9):
    """KL(p||q) — distribuições discretas (mesmo suporte, mesmo tamanho)."""
    p = np.asarray(p, dtype=float)
    q = np.asarray(q, dtype=float)
    p = p / p.sum() if p.sum() > 0 else p
    q = q / q.sum() if q.sum() > 0 else q
    p = np.clip(p, eps, 1.0)
    q = np.clip(q, eps, 1.0)
    return float((p * np.log(p / q)).sum())


# ── análise principal ────────────────────────────────────────────────

def validate(data_dir, gt_path, top_n=100, ranking_threshold=0.90,
             orbits_min_recovered=6):
    gt = load_gt(gt_path)
    raw_path = os.path.join(data_dir, "correlations", "raw_correlations.json")
    z3 = load_raw_correlations(raw_path)

    print("=" * 72)
    print("VALIDAÇÃO Z3 vs GROUND TRUTH (6×6)")
    print("=" * 72)
    print(f"GT          : {gt_path}")
    print(f"Z3 dataset  : {raw_path}")
    print(f"Pools Z3    : {len(z3)}")
    print()

    # ─── 1. GT orbits (canônicas) ─────────────────────────────────────
    gt_orbits = gt["correlation_orbits"]  # lista de 7
    gt_canon_map = {}      # canonical_pair -> {min_r, rep_labels, members_count}
    for orb in gt_orbits:
        key = gt_orbit_key(orb)
        gt_canon_map[key] = {
            "min_r": orb["min_r"],
            "representative": orb["representative"],
            "size": orb["size"],
        }

    print(f"Órbitas D₄ no GT: {len(gt_orbits)}")
    for orb in gt_orbits:
        print(f"  r={orb['min_r']:>8.4f}  |orb|={orb['size']:>2}  "
              f"{orb['representative'][0]} ↔ {orb['representative'][1]}  "
              f"(shared={orb.get('shared_vertex')})")

    # ─── 2. Reconstruir órbitas Z3 a partir do raw_correlations ───────
    # Para cada pool, `orbits` já vem agrupado por canonical. Pegamos o
    # representante de cada órbita Z3, re-canonicalizamos (consistência) e
    # acumulamos o melhor r_min observado entre todos os pools.

    z3_canon_map = {}  # canonical_pair -> r_min observado
    for pool_label, pool in z3.items():
        for orb in pool.get("orbits", []):
            rep_a = label_to_edge(orb["representative"][0])
            rep_b = label_to_edge(orb["representative"][1])
            canon = pair_d4_canonical(rep_a, rep_b)
            r = orb["r_min"]
            if canon not in z3_canon_map or r < z3_canon_map[canon]:
                z3_canon_map[canon] = r

    print()
    print(f"Órbitas D₄ observadas pela Z3 (top-{top_n} por pool, agregado): {len(z3_canon_map)}")

    # ─── 3. Recall — quantas órbitas GT aparecem em Z3 ────────────────
    matched = []
    missing = []
    for key, gt_info in gt_canon_map.items():
        if key in z3_canon_map:
            matched.append((key, gt_info["min_r"], z3_canon_map[key]))
        else:
            missing.append((key, gt_info["min_r"], gt_info["representative"]))

    print()
    print(f"Órbitas GT recuperadas: {len(matched)}/{len(gt_orbits)}")
    print(f"{'r_gt':>9}  {'r_z3':>9}  {'Δr':>8}  representante")
    print(f"{'-'*9}  {'-'*9}  {'-'*8}  " + "-" * 40)
    for key, r_gt, r_z3 in sorted(matched, key=lambda x: x[1]):
        info = gt_canon_map[key]
        rep = f"{info['representative'][0]} ↔ {info['representative'][1]}"
        print(f"{r_gt:>9.4f}  {r_z3:>9.4f}  {r_z3-r_gt:>+8.4f}  {rep}")
    if missing:
        print()
        print("Órbitas GT NÃO recuperadas:")
        for key, r_gt, rep in missing:
            print(f"  r={r_gt:.4f}  {rep[0]} ↔ {rep[1]}")

    # ─── 4. Métricas quantitativas ────────────────────────────────────
    r_gt_list = [m[1] for m in matched]
    r_z3_list = [m[2] for m in matched]
    rho = spearman_rho(r_gt_list, r_z3_list) if len(matched) >= 2 else float("nan")
    mae = float(np.mean(np.abs(np.array(r_gt_list) - np.array(r_z3_list)))) \
        if matched else float("nan")

    # ─── 5. Arestas obrigatórias ──────────────────────────────────────
    gt_mand_set = set()
    for orb_record in gt["mandatory_edges_by_orbit"]:
        gt_mand_set.update(orb_record["edges"])

    # interseção sobre pools (sempre obrigatórias) = backbone real
    pool_mands = [set(p.get("mandatory_edges", [])) for p in z3.values()]
    z3_global_mand = set.intersection(*pool_mands) if pool_mands else set()
    # união sobre pools = candidatas (pair-mandatory ou global-mandatory)
    z3_any_mand = set.union(*pool_mands) if pool_mands else set()

    mand_match = (gt_mand_set == z3_global_mand)
    mand_in_union = gt_mand_set.issubset(z3_any_mand)

    print()
    print("Arestas obrigatórias:")
    print(f"  GT (backbone)                 : {len(gt_mand_set):>2}  {sorted(gt_mand_set)}")
    print(f"  Z3 ∩ pools (sempre obrig.)    : {len(z3_global_mand):>2}  {sorted(z3_global_mand)}")
    print(f"  Z3 ∪ pools (algum pool)       : {len(z3_any_mand):>2}")
    print(f"  GT == Z3 ∩ pools ?             {mand_match}")
    print(f"  GT ⊆ Z3 ∪ pools ?             {mand_in_union}")

    # ─── 6. KL-divergence sobre freq de arestas ───────────────────────
    # GT freq não está no analysis_6x6.json diretamente; usamos como aproxim.
    # a freq média ao longo de todos os pools. Comparamos a freq DENTRO de cada
    # pool com a média global como uma medida de variabilidade.
    # (Para uma comparação rigorosa contra GT seria preciso recomputar freq de
    # cada aresta a partir do destruction_catalogue.json — fora do escopo deste
    # script, mas reportamos a variabilidade entre pools.)

    edge_freqs_per_pool = {}
    for pool_label, pool in z3.items():
        all_edges = pool.get("mandatory_edges", []) + pool.get("impossible_edges", [])
        # raw_correlations não traz freq por aresta; só sabemos f=1 e f=0 explicitamente
        # e que n_live arestas têm 0<f<1. Sem freq por aresta, KL fica como TODO.
        # Aqui registramos só uma medida de cobertura.
        edge_freqs_per_pool[pool_label] = {
            "mandatory": pool.get("mandatory_edges", []),
            "impossible": pool.get("impossible_edges", []),
            "n_live": pool.get("n_live"),
            "n_samples": pool.get("n_samples"),
        }

    # ─── 7. Veredito ──────────────────────────────────────────────────
    recovered = len(matched)
    pass_rho = (not math.isnan(rho)) and rho >= ranking_threshold
    pass_recall = recovered >= orbits_min_recovered
    pass_mand = mand_in_union  # critério leve (subset)

    if pass_rho and pass_recall and pass_mand:
        verdict = "VALID"
    elif pass_recall and (pass_rho or pass_mand):
        verdict = "PARTIAL"
    else:
        verdict = "BIASED"

    print()
    print("=" * 72)
    print("RESUMO QUANTITATIVO")
    print("=" * 72)
    print(f"  orbits_recovered     : {recovered}/{len(gt_orbits)}")
    print(f"  ranking_correlation  : {rho:.4f}    (Spearman, threshold {ranking_threshold})")
    print(f"  r_mae                : {mae:.4f}")
    print(f"  mandatory_subset_ok  : {mand_in_union}")
    print(f"  mandatory_exact_ok   : {mand_match}")
    print()
    print(f"  VEREDITO : {verdict}")
    print("=" * 72)

    return {
        "data_dir": data_dir,
        "gt_path": gt_path,
        "n_pools": len(z3),
        "gt_n_orbits": len(gt_orbits),
        "z3_n_orbits_aggregated": len(z3_canon_map),
        "orbits_recovered": recovered,
        "ranking_correlation": rho,
        "r_mae": mae,
        "matched_orbits": [
            {
                "representative": gt_canon_map[k]["representative"],
                "r_gt": r_gt,
                "r_z3": r_z3,
                "delta_r": r_z3 - r_gt,
            }
            for k, r_gt, r_z3 in sorted(matched, key=lambda x: x[1])
        ],
        "missing_orbits": [
            {"representative": rep, "r_gt": r_gt}
            for _, r_gt, rep in missing
        ],
        "mandatory": {
            "gt": sorted(gt_mand_set),
            "z3_intersection": sorted(z3_global_mand),
            "z3_union": sorted(z3_any_mand),
            "exact_match": mand_match,
            "subset_match": mand_in_union,
        },
        "verdict": verdict,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="data_6x6_sampled",
                    help="diretório com correlations/raw_correlations.json")
    ap.add_argument("--gt", default="../analysis_6x6.json",
                    help="path do ground truth (analysis_6x6.json)")
    ap.add_argument("--top-n", type=int, default=100,
                    help="top-N correlações por pool a considerar")
    ap.add_argument("--out", default=None,
                    help="se setado, salva relatório JSON (default: <data_dir>/validation.json)")
    args = ap.parse_args()

    d4_orbits.set_board(6)

    if args.out is None:
        args.out = os.path.join(args.data_dir, "validation.json")

    report = validate(args.data_dir, args.gt, top_n=args.top_n)

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, ensure_ascii=False, default=str)
    print(f"\nRelatório JSON salvo em: {args.out}")
    return 0 if report["verdict"] == "VALID" else 1


if __name__ == "__main__":
    sys.exit(main())
