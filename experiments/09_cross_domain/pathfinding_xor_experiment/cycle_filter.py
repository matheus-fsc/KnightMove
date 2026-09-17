#!/usr/bin/env python3
"""
Follow-up: filtro de comprimento de ciclo + XOR multi-ciclo + preditor de
conectividade.
============================================================================

Estudo anterior (pathfinding_xor.py): comprimento do ciclo e o preditor mais
forte de compatibilidade do XOR (r = -0.665); modo de falha dominante e
desconexao (90.8%). Hipotese aqui: filtrar a base de ciclos para ciclos
curtos/locais ANTES do XOR aumenta a taxa de compatibilidade e pode melhorar a
diversidade.

NAO regenera mapas nem reroda o A* repetido. As instancias sao reconstruidas
de forma DETERMINISTICA a partir de `used_seed` (mesmos 60 mapas), e a
diversidade do A* repetido e RECARREGADA do JSON anterior.

Uso:
    ../venv/bin/python cycle_filter.py
"""

import json
import os
import time
from collections import defaultdict

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, precision_score, recall_score

# reusa exatamente o pipeline do experimento anterior
import pathfinding_xor as base

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(HERE, "results")
PRIOR_JSON = os.path.join(RESULTS_DIR, "pathfinding_xor_experiment.json")
OUT_JSON = os.path.join(RESULTS_DIR, "cycle_filter_experiment.json")
OUT_MD = os.path.join(RESULTS_DIR, "cycle_filter_summary.md")
OUT_PRED = os.path.join(RESULTS_DIR, "connectivity_predictor.json")

FILTERS = [4, 6, 8, 10, 20]          # comprimento maximo de ciclo
BASELINE_KEY = "Baseline"
N_TRIPLES = 500
N_PAIRS = 500
MAX_SINGLE_FOR_B = 50
BBOX_LOCAL = 5                        # ciclo local = caixa <= 5x5 celulas


# --------------------------------------------------------------------------
# Reconstroi uma instancia identica a partir do used_seed (deterministico).
# Retorna grafo, S, T, caminho-base (edge set), e lista de ciclos.
# --------------------------------------------------------------------------
def rebuild_instance(n, used_seed):
    G, S, T, seed_ok = base.gen_graph(n, used_seed)
    assert seed_ok == used_seed, f"seed drift {seed_ok} != {used_seed}"
    base_seq = base.astar_path(G, S, T)
    P = base.seq_to_edges(base_seq)
    cycles = base.fundamental_cycles(G, S)
    return G, S, T, P, cycles


def cycle_bbox_local(cyc):
    """True se todos os vertices do ciclo cabem numa caixa <= 5x5."""
    rs = [v[0] for e in cyc for v in tuple(e)]
    cs = [v[1] for e in cyc for v in tuple(e)]
    return (max(rs) - min(rs) <= BBOX_LOCAL - 1) and \
           (max(cs) - min(cs) <= BBOX_LOCAL - 1)


# --------------------------------------------------------------------------
# EXPERIMENT A — filtro de comprimento
# --------------------------------------------------------------------------
def run_filter(P, cycles, S, T, thr):
    """thr=None -> baseline (todos). Retorna dict de metricas + caminhos validos."""
    pool = cycles if thr is None else [c for c in cycles if len(c) <= thr]
    n_avail = len(pool)
    if n_avail == 0:
        return {
            "n_cycles": 0, "compat": float("nan"), "valid_paths": 0,
            "diversity": 0.0, "time_ms": 0.0,
            "fail_degree": 0, "fail_disconnected": 0, "fail_both": 0,
        }, []
    t0 = time.perf_counter()
    valid = []
    fails = {"degree": 0, "disconnected": 0, "both": 0}
    for c in pool:
        cand = base.xor(P, set(c))
        ok, reason = base.check_path(cand, S, T)
        if ok:
            valid.append(cand)
        else:
            fails[reason] += 1
    t_ms = (time.perf_counter() - t0) * 1000.0
    return {
        "n_cycles": n_avail,
        "compat": len(valid) / n_avail,
        "valid_paths": len(valid),
        "diversity": base.diversity(valid),
        "time_ms": t_ms,
        "fail_degree": fails["degree"],
        "fail_disconnected": fails["disconnected"],
        "fail_both": fails["both"],
    }, valid


# --------------------------------------------------------------------------
# EXPERIMENT B — XOR multi-ciclo (Filter_8)
# --------------------------------------------------------------------------
def run_multi(P, cycles, S, T, rng):
    """Usa apenas ciclos com len<=8. single/pair/triple + combinado."""
    f8 = [c for c in cycles if len(c) <= 8]
    # single
    single_valid = []
    for c in f8:
        cand = base.xor(P, set(c))
        ok, _ = base.check_path(cand, S, T)
        if ok:
            single_valid.append(cand)
            if len(single_valid) >= MAX_SINGLE_FOR_B:
                # ainda continuamos so para diversidade? spec: ate 50 single
                pass
    single_capped = single_valid[:MAX_SINGLE_FOR_B]

    k = len(f8)
    # pares (amostra ate 500)
    pair_valid = []
    if k >= 2:
        pairs = _sample_combos(k, 2, N_PAIRS, rng)
        for (i, j) in pairs:
            cand = base.xor(base.xor(P, set(f8[i])), set(f8[j]))
            ok, _ = base.check_path(cand, S, T)
            if ok:
                pair_valid.append(cand)

    # triplas (amostra ate 500; com reposicao de indices se k<3 -- documentado)
    triple_valid = []
    replacement_used = False
    if k >= 1:
        triples, replacement_used = _sample_triples(k, N_TRIPLES, rng)
        for (i, j, l) in triples:
            cand = base.xor(base.xor(base.xor(P, set(f8[i])), set(f8[j])), set(f8[l]))
            ok, _ = base.check_path(cand, S, T)
            if ok:
                triple_valid.append(cand)

    combined = single_capped + pair_valid + triple_valid
    return {
        "n_f8": k,
        "div_single": base.diversity(single_capped),
        "div_pair": base.diversity(pair_valid),
        "div_triple": base.diversity(triple_valid),
        "div_combined": base.diversity(combined),
        "n_single": len(single_capped),
        "n_pair": len(pair_valid),
        "n_triple": len(triple_valid),
        "triple_replacement_used": replacement_used,
    }


def _sample_combos(k, r, limit, rng):
    """Amostra ate `limit` combinacoes distintas (i<j) sem reposicao."""
    import math
    total = math.comb(k, r)
    if total <= limit:
        if r == 2:
            return [(i, j) for i in range(k) for j in range(i + 1, k)]
    seen = set()
    out = []
    tries = 0
    while len(out) < min(limit, total) and tries < limit * 50:
        tries += 1
        idx = tuple(sorted(rng.choice(k, size=r, replace=False).tolist()))
        if idx in seen:
            continue
        seen.add(idx)
        out.append(idx)
    return out


def _sample_triples(k, limit, rng):
    """Amostra ate `limit` triplas (i<j<l). Se k<3, amostra COM reposicao
    de indices (documentado via flag). Retorna (lista, replacement_used)."""
    import math
    if k >= 3:
        total = math.comb(k, 3)
        if total <= limit:
            return ([(i, j, l) for i in range(k) for j in range(i + 1, k)
                     for l in range(j + 1, k)], False)
        seen = set()
        out = []
        tries = 0
        while len(out) < min(limit, total) and tries < limit * 50:
            tries += 1
            idx = tuple(sorted(rng.choice(k, size=3, replace=False).tolist()))
            if idx in seen:
                continue
            seen.add(idx)
            out.append(idx)
        return out, False
    # k < 3: precisa de reposicao
    out = []
    for _ in range(min(limit, k ** 3)):
        idx = tuple(int(x) for x in rng.integers(0, k, size=3))
        out.append(idx)
    return out, True


# --------------------------------------------------------------------------
# EXPERIMENT C — features por candidato (todos os ciclos, valido ou nao)
# --------------------------------------------------------------------------
def collect_features(P, cycles, S, T):
    base_len = len(P)
    rows = []
    for c in cycles:
        cand = base.xor(P, set(c))
        ok, reason = base.check_path(cand, S, T)
        n_shared = len(c & P)
        feat = {
            "cycle_length": len(c),
            "is_local": 1 if cycle_bbox_local(c) else 0,
            "shares_edges": 1 if n_shared > 0 else 0,
            "n_shared_edges": n_shared,
            "path_length_ratio": len(c) / base_len if base_len else 0.0,
            "valid": 1 if ok else 0,
            "result": "valid" if ok else (
                "degree_fail" if reason in ("degree", "both") else "disconnected"),
        }
        rows.append(feat)
    return rows


# --------------------------------------------------------------------------
def nanmean(a):
    a = np.array(a, dtype=float)
    return float(np.nanmean(a)) if np.any(~np.isnan(a)) else float("nan")


def nanstd(a):
    a = np.array(a, dtype=float)
    return float(np.nanstd(a)) if np.any(~np.isnan(a)) else float("nan")


def main():
    prior = json.load(open(PRIOR_JSON))
    prior_by = {(r["grid_size"], r["seed"]): r for r in prior}

    per_instance = []
    feat_rows = []

    t_start = time.perf_counter()
    for r in prior:
        n = r["grid_size"]
        seed = r["seed"]
        used = r["used_seed"]
        rng = np.random.default_rng(used)

        G, S, T, P, cycles = rebuild_instance(n, used)

        # sanidade: bate com o numero de ciclos do estudo anterior
        assert len(cycles) == r["cycles_collected"], \
            f"cycle count drift n={n} seed={seed}: {len(cycles)} vs {r['cycles_collected']}"

        # EXP A
        filt = {}
        for thr in FILTERS:
            m, _ = run_filter(P, cycles, S, T, thr)
            filt[f"Filter_{thr}"] = m
        mb, _ = run_filter(P, cycles, S, T, None)
        filt[BASELINE_KEY] = mb

        # EXP B
        multi = run_multi(P, cycles, S, T, rng)
        multi["diversity_astar_prior"] = r["diversity_astar"]

        # EXP C features
        feat_rows.extend(collect_features(P, cycles, S, T))

        per_instance.append({
            "grid_size": n, "seed": seed, "used_seed": used,
            "filters": filt, "multi": multi,
        })

    elapsed = time.perf_counter() - t_start

    # ---------------- EXP A agregacao ----------------
    filter_names = [f"Filter_{t}" for t in FILTERS] + [BASELINE_KEY]
    aggA = {}
    for fn in filter_names:
        ncyc = [p["filters"][fn]["n_cycles"] for p in per_instance]
        comp = [p["filters"][fn]["compat"] for p in per_instance]
        vp = [p["filters"][fn]["valid_paths"] for p in per_instance]
        dv = [p["filters"][fn]["diversity"] for p in per_instance]
        tm = [p["filters"][fn]["time_ms"] for p in per_instance]
        fdis = sum(p["filters"][fn]["fail_disconnected"] for p in per_instance)
        fdeg = sum(p["filters"][fn]["fail_degree"] for p in per_instance)
        fboth = sum(p["filters"][fn]["fail_both"] for p in per_instance)
        ftot = fdis + fdeg + fboth
        aggA[fn] = {
            "n_cycles_mean": float(np.mean(ncyc)), "n_cycles_std": float(np.std(ncyc)),
            "compat_mean": nanmean(comp), "compat_std": nanstd(comp),
            "valid_paths_mean": float(np.mean(vp)), "valid_paths_std": float(np.std(vp)),
            "diversity_mean": float(np.mean(dv)), "diversity_std": float(np.std(dv)),
            "time_ms_mean": float(np.mean(tm)),
            "n_instances_no_cycles": int(sum(1 for x in ncyc if x == 0)),
            "disconnection_rate_of_fails": (fdis / ftot) if ftot else float("nan"),
            "degree_rate_of_fails": (fdeg / ftot) if ftot else float("nan"),
        }

    best_compat = max(filter_names, key=lambda fn: (aggA[fn]["compat_mean"]
                      if not np.isnan(aggA[fn]["compat_mean"]) else -1))
    # tradeoff: valid_paths * diversity
    def tradeoff(fn):
        return aggA[fn]["valid_paths_mean"] * aggA[fn]["diversity_mean"]
    best_tradeoff = max(filter_names, key=tradeoff)

    # ---------------- EXP B agregacao ----------------
    ds = [p["multi"]["div_single"] for p in per_instance]
    dp = [p["multi"]["div_pair"] for p in per_instance]
    dt = [p["multi"]["div_triple"] for p in per_instance]
    dc = [p["multi"]["div_combined"] for p in per_instance]
    da = [p["multi"]["diversity_astar_prior"] for p in per_instance]
    repl_count = sum(1 for p in per_instance if p["multi"]["triple_replacement_used"])
    aggB = {
        "div_single_mean": float(np.mean(ds)), "div_single_std": float(np.std(ds)),
        "div_pair_mean": float(np.mean(dp)), "div_pair_std": float(np.std(dp)),
        "div_triple_mean": float(np.mean(dt)), "div_triple_std": float(np.std(dt)),
        "div_combined_mean": float(np.mean(dc)), "div_combined_std": float(np.std(dc)),
        "div_astar_mean": float(np.mean(da)), "div_astar_std": float(np.std(da)),
        "gap_closed_pp": (np.mean(dc) - np.mean(ds)) * 100.0,
        "gap_to_astar_pp": (np.mean(da) - np.mean(dc)) * 100.0,
        "instances_triple_replacement": repl_count,
    }

    # ---------------- EXP C preditor ----------------
    X = np.array([[f["cycle_length"], f["is_local"], f["shares_edges"],
                   f["n_shared_edges"], f["path_length_ratio"]] for f in feat_rows],
                 dtype=float)
    y = np.array([f["valid"] for f in feat_rows], dtype=int)
    feat_labels = ["cycle_length", "is_local", "shares_edges",
                   "n_shared_edges", "path_length_ratio"]

    # padroniza para coeficientes comparaveis
    mu, sd = X.mean(axis=0), X.std(axis=0)
    sd[sd == 0] = 1.0
    Xs = (X - mu) / sd
    clf = LogisticRegression(max_iter=1000, class_weight="balanced")
    clf.fit(Xs, y)
    proba = clf.predict_proba(Xs)[:, 1]
    pred = (proba >= 0.5).astype(int)
    auc = float(roc_auc_score(y, proba))
    prec = float(precision_score(y, pred, zero_division=0))
    rec = float(recall_score(y, pred, zero_division=0))
    coefs = {lbl: float(c) for lbl, c in zip(feat_labels, clf.coef_[0])}
    top_feat = max(coefs, key=lambda k: abs(coefs[k]))

    predictor = {
        "n_candidates": int(len(y)),
        "n_valid": int(y.sum()),
        "base_rate_valid": float(y.mean()),
        "auc_roc": auc,
        "precision_at_0.5": prec,
        "recall_at_0.5": rec,
        "coefficients_standardized": coefs,
        "most_predictive_feature": top_feat,
        "most_predictive_coef": coefs[top_feat],
        "useful_pre_filter_auc_gt_0.75": auc > 0.75,
        "feature_order": feat_labels,
    }

    # ---------------- salvar ----------------
    os.makedirs(RESULTS_DIR, exist_ok=True)
    json.dump({
        "per_instance": per_instance,
        "aggregate_A": aggA,
        "aggregate_B": aggB,
        "best_filter_compat": best_compat,
        "best_filter_tradeoff": best_tradeoff,
        "elapsed_s": elapsed,
    }, open(OUT_JSON, "w"), indent=2)
    json.dump(predictor, open(OUT_PRED, "w"), indent=2)

    print_report(aggA, aggB, predictor, best_compat, best_tradeoff,
                 filter_names, per_instance)
    write_md(aggA, aggB, predictor, best_compat, best_tradeoff, filter_names)
    print(f"\nSaved:\n  {OUT_JSON}\n  {OUT_MD}\n  {OUT_PRED}")
    print(f"(new computation time: {elapsed:.1f}s)")


# --------------------------------------------------------------------------
def print_report(aggA, aggB, pred, best_compat, best_tradeoff, fnames, per):
    print("\n=== EXPERIMENT A: CYCLE LENGTH FILTER ===\n")
    print("          n_cycles   compat%      valid_paths   diversity")
    print("          --------   -------      -----------   ---------")
    for fn in fnames:
        a = aggA[fn]
        cm = a["compat_mean"] * 100
        cs = a["compat_std"] * 100
        print(f"{fn:9s} {a['n_cycles_mean']:5.0f}+-{a['n_cycles_std']:<4.0f} "
              f"{cm:5.1f}+-{cs:<4.1f}% "
              f"{a['valid_paths_mean']:5.1f}+-{a['valid_paths_std']:<4.1f} "
              f"{a['diversity_mean']:.3f}+-{a['diversity_std']:.3f}")
    print(f"\nBest filter by compatibility rate: {best_compat}")
    print(f"Best filter by valid_paths x diversity tradeoff: {best_tradeoff}")

    print("\n=== EXPERIMENT B: MULTI-CYCLE DIVERSITY (Filter_8) ===\n")
    print(f"Filter_8 single XOR diversity:  {aggB['div_single_mean']:.3f}+-{aggB['div_single_std']:.3f}")
    print(f"Filter_8 pair   XOR diversity:  {aggB['div_pair_mean']:.3f}+-{aggB['div_pair_std']:.3f}")
    print(f"Filter_8 triple XOR diversity:  {aggB['div_triple_mean']:.3f}+-{aggB['div_triple_std']:.3f}")
    print(f"Combined (single+pair+triple):  {aggB['div_combined_mean']:.3f}+-{aggB['div_combined_std']:.3f}")
    print(f"Repeated A* diversity (prior):  {aggB['div_astar_mean']:.3f}+-{aggB['div_astar_std']:.3f}  [reloaded]")
    print(f"\nGap closed (combined vs single): {aggB['gap_closed_pp']:+.1f} pp")
    print(f"Remaining gap to A*:             {aggB['gap_to_astar_pp']:+.1f} pp")
    if aggB["instances_triple_replacement"]:
        print(f"(triplas com reposicao em {aggB['instances_triple_replacement']} instancias com <3 ciclos f8)")

    print("\n=== EXPERIMENT C: CONNECTIVITY PREDICTOR ===\n")
    print(f"Candidates: {pred['n_candidates']} (valid base rate {pred['base_rate_valid']*100:.1f}%)")
    print(f"AUC-ROC: {pred['auc_roc']:.3f}")
    print(f"Most predictive feature: {pred['most_predictive_feature']} "
          f"(coef={pred['most_predictive_coef']:+.3f}, standardized)")
    print(f"Precision at 0.5: {pred['precision_at_0.5']:.3f}  |  "
          f"Recall at 0.5: {pred['recall_at_0.5']:.3f}")
    print(f"Standardized coefficients: "
          f"{ {k: round(v,3) for k,v in pred['coefficients_standardized'].items()} }")
    print(f"\nPredictor useful as pre-filter (AUC > 0.75): "
          f"{'YES' if pred['useful_pre_filter_auc_gt_0.75'] else 'NO'}")

    # ---- key questions ----
    base_compat = aggA["Baseline"]["compat_mean"]
    f8_compat = aggA["Filter_8"]["compat_mean"]
    print("\n=== KEY QUESTIONS ===")
    ratio = f8_compat / base_compat if base_compat else float("nan")
    print(f"1. Filtro curto aumenta compatibilidade? Filter_8={f8_compat*100:.1f}% "
          f"vs Baseline={base_compat*100:.1f}% = {ratio:.2f}x "
          f"({'CONFIRMADA (>2x)' if ratio > 2 else 'REFUTADA (<2x)'})")
    print(f"2. Multi-ciclo fecha gap de diversidade? combined="
          f"{aggB['div_combined_mean']:.3f} "
          f"({'aproxima 0.40+' if aggB['div_combined_mean'] >= 0.40 else 'NAO atinge 0.40 -> REFUTADA'}); "
          f"A*={aggB['div_astar_mean']:.3f}")
    print(f"3. Propriedades preveem validade? AUC={pred['auc_roc']:.3f} "
          f"({'CONFIRMADA (>0.75)' if pred['auc_roc']>0.75 else 'REFUTADA (<0.75)'})")
    # optimal tradeoff compat x n_cycles
    def cxn(fn):
        c = aggA[fn]["compat_mean"]
        return (c if not np.isnan(c) else 0) * aggA[fn]["n_cycles_mean"]
    best_cxn = max(fnames, key=cxn)
    print(f"4. Threshold otimo (compat x n_cycles): {best_cxn} "
          f"(= valid_paths esperados {cxn(best_cxn):.1f})")

    # ---- connectivity vs filter ----
    print("\n=== CONNECTIVITY OBSTRUCTION vs FILTER ===")
    for fn in fnames:
        dr = aggA[fn]["disconnection_rate_of_fails"]
        print(f"  {fn:9s}: desconexao = {dr*100:5.1f}% das falhas" if not np.isnan(dr)
              else f"  {fn:9s}: (sem falhas)")


def write_md(aggA, aggB, pred, best_compat, best_tradeoff, fnames):
    L = []
    L.append("# Cycle-length filter + multi-cycle XOR + connectivity predictor\n")
    L.append("Follow-up de `pathfinding_xor.py`. Mesmas 60 instancias "
             "(reconstruidas deterministicamente via `used_seed`); A* repetido "
             "NAO foi re-executado (diversidade recarregada do JSON anterior).\n")
    L.append("## EXPERIMENT A — filtro de comprimento\n")
    L.append("| Filtro | n_cycles | compat% | valid_paths | diversity | desconexao %falhas |")
    L.append("|---|---:|---:|---:|---:|---:|")
    for fn in fnames:
        a = aggA[fn]
        dr = a["disconnection_rate_of_fails"]
        drs = f"{dr*100:.1f}%" if not np.isnan(dr) else "—"
        cm = a["compat_mean"] * 100
        L.append(f"| {fn} | {a['n_cycles_mean']:.0f}±{a['n_cycles_std']:.0f} | "
                 f"{cm:.1f}±{a['compat_std']*100:.1f}% | "
                 f"{a['valid_paths_mean']:.1f}±{a['valid_paths_std']:.1f} | "
                 f"{a['diversity_mean']:.3f}±{a['diversity_std']:.3f} | {drs} |")
    L.append(f"\n**Melhor por compatibilidade:** {best_compat}  ")
    L.append(f"**Melhor por valid_paths×diversity:** {best_tradeoff}\n")

    base_compat = aggA["Baseline"]["compat_mean"]
    f8 = aggA["Filter_8"]["compat_mean"]
    ratio = f8 / base_compat if base_compat else float("nan")

    L.append("## EXPERIMENT B — diversidade multi-ciclo (Filter_8)\n")
    L.append("| Conjunto | diversidade |")
    L.append("|---|---:|")
    L.append(f"| single | {aggB['div_single_mean']:.3f}±{aggB['div_single_std']:.3f} |")
    L.append(f"| pair | {aggB['div_pair_mean']:.3f}±{aggB['div_pair_std']:.3f} |")
    L.append(f"| triple | {aggB['div_triple_mean']:.3f}±{aggB['div_triple_std']:.3f} |")
    L.append(f"| combined | {aggB['div_combined_mean']:.3f}±{aggB['div_combined_std']:.3f} |")
    L.append(f"| A* repetido (prior) | {aggB['div_astar_mean']:.3f}±{aggB['div_astar_std']:.3f} |")
    L.append(f"\nGap fechado (combined−single): {aggB['gap_closed_pp']:+.1f} pp · "
             f"gap restante p/ A*: {aggB['gap_to_astar_pp']:+.1f} pp\n")

    L.append("## EXPERIMENT C — preditor de conectividade\n")
    L.append(f"- Candidatos: {pred['n_candidates']} (taxa-base válida {pred['base_rate_valid']*100:.1f}%)")
    L.append(f"- **AUC-ROC: {pred['auc_roc']:.3f}**")
    L.append(f"- Feature mais preditiva: `{pred['most_predictive_feature']}` "
             f"(coef padronizado {pred['most_predictive_coef']:+.3f})")
    L.append(f"- Precisão@0.5: {pred['precision_at_0.5']:.3f} · Recall@0.5: {pred['recall_at_0.5']:.3f}")
    L.append(f"- Coeficientes (padronizados): "
             f"{ {k: round(v,3) for k,v in pred['coefficients_standardized'].items()} }")
    L.append(f"- Útil como pré-filtro (AUC>0.75): "
             f"**{'SIM' if pred['useful_pre_filter_auc_gt_0.75'] else 'NÃO'}**\n")

    L.append("## Respostas às perguntas-chave (honestas)\n")
    L.append(f"1. **Filtro curto aumenta compatibilidade?** Filter_8 = {f8*100:.1f}% vs "
             f"Baseline = {base_compat*100:.1f}% → **{ratio:.2f}×** — "
             f"{'CONFIRMADA (>2×)' if ratio>2 else 'hipótese de 2× REFUTADA' if ratio<2 else 'no limite'}.")
    L.append(f"2. **Multi-ciclo fecha o gap de diversidade?** combined = "
             f"{aggB['div_combined_mean']:.3f}; alvo 0.40+ → "
             f"{'atingido' if aggB['div_combined_mean']>=0.40 else '**REFUTADA**'}. "
             f"A* = {aggB['div_astar_mean']:.3f}.")
    L.append(f"3. **Propriedades preveem validade?** AUC = {pred['auc_roc']:.3f} → "
             f"{'CONFIRMADA (>0.75)' if pred['auc_roc']>0.75 else '**REFUTADA (<0.75)**'}.")
    def cxn(fn):
        c = aggA[fn]["compat_mean"]
        return (c if not np.isnan(c) else 0) * aggA[fn]["n_cycles_mean"]
    best_cxn = max(fnames, key=cxn)
    L.append(f"4. **Threshold ótimo (compat×n_cycles):** {best_cxn} "
             f"(≈{cxn(best_cxn):.1f} caminhos válidos esperados).\n")

    L.append("## Achado central (contraintuitivo)\n")
    L.append("O filtro de comprimento **NÃO** aumenta a compatibilidade "
             "(plana em ~17–18% de Filter_4 a Baseline). A correlação prévia "
             "comprimento×compat (r=−0.665) era **entre instâncias** "
             "(confundida com o tamanho da grade), não **por ciclo dentro** de "
             "uma instância. O preditor real é `shares_edges`: o ciclo "
             "compartilha alguma aresta com o caminho-base? (coef padronizado "
             "**+5.27** vs `cycle_length` +0.13). Recall@0.5 = 0.99: quase todo "
             "ciclo válido toca o caminho-base. O pré-filtro útil é "
             "**'compartilha aresta com P'**, não 'é curto'.\n")

    L.append("## Conexão com o trabalho anterior\n")
    L.append("1. **A obstrução de conectividade muda com o filtro?** NÃO — "
             "piora. Coluna `desconexao %falhas`: 90.8% (Baseline) sobe para "
             "98.4% (Filter_4). Ciclos curtos que não tocam P falham por pura "
             "desconexão (sem violação de grau), então filtrar por comprimento "
             "*concentra* as falhas em desconexão em vez de eliminá-las. "
             "Confirma o achado do Knight's Tour: conectividade é a obstrução "
             "GLOBAL, indiferente a filtros locais ([[project_residual_search_10x10_done]]).")
    L.append("2. **Relação com o teorema de suficiência local (D&C do Knight's "
             "Tour):** o que torna um ciclo 'seguro' não é a localidade "
             "geométrica (`is_local` tem coef só +0.56) e sim **compartilhar um "
             "trecho contíguo com o caminho-base** (`shares_edges` +5.27). É a "
             "mesma lógica do D&C: a compatibilidade exige sobreposição na "
             "fronteira (interface comum), não mera proximidade — blocos "
             "adjacentes colam quando casam na borda compartilhada, não quando "
             "são pequenos ([[project_dnc_blocks]]).")
    L.append("3. **Veredito final:** o XOR GF(2) é viável como ENUMERADOR de "
             "alternativas em grafos esparsos, mas o filtro de comprimento é a "
             "alavanca errada: não muda a compatibilidade nem fecha o gap de "
             "diversidade (combined 0.198 « A* 0.477). A conectividade global "
             "permanece a obstrução insuperável por filtragem local — só é "
             "**contornável** restringindo-se a ciclos que tocam P (pré-filtro "
             "`shares_edges`, AUC 0.977), o que reduz drasticamente o custo sem "
             "sacrificar caminhos válidos. Em diversidade, porém, o A* "
             "penalizado continua estruturalmente superior: XOR produz "
             "variações locais, A* produz rotas globalmente distintas.")
    open(OUT_MD, "w").write("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
