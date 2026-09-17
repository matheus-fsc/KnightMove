#!/usr/bin/env python3
"""
Structural characterization of η (crossing efficiency).
========================================================

η(G) = compat / coverage = #valid_XOR / #cycles_touching_path
     = P(contiguous crossing | overlap≥1, G)

Goal: predict η from measurable graph features and, if possible, give a simple
interpretable formula. Honest: if the top predictor is NOT frac_degree_2, we
lead with whatever actually wins.

Per-instance η was NOT stored by benchmarks 1/2 (only cell aggregates), but the
generators are deterministic, so we REGENERATE the identical graphs and compute
η + features freshly. path_selection (20) and cycle_filter (60) instances are
reconstructed via used_seed and used as genuinely HELD-OUT data.

Uso:
    ../venv/bin/python eta_characterization.py build     # build dataset cache
    ../venv/bin/python eta_characterization.py analyze   # regression + outputs
    ../venv/bin/python eta_characterization.py all
"""

import json
import os
import signal
import sys

import numpy as np
import networkx as nx
from scipy import stats as sstats
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import PolynomialFeatures, StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import cross_val_score, train_test_split
from sklearn.metrics import r2_score, mean_squared_error
from sklearn.impute import SimpleImputer

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import pathfinding_xor as base                       # noqa: E402
import validation_benchmark as vb                    # noqa: E402

RES = os.path.join(HERE, "results")
CACHE = os.path.join(RES, "eta_dataset_cache.json")
OUT_MAIN = os.path.join(RES, "eta_characterization.json")
OUT_PLOTS = os.path.join(RES, "eta_characterization_plots.json")
OUT_VALID = os.path.join(RES, "eta_formula_validation.json")
OUT_MD = os.path.join(RES, "eta_summary.md")

FEATURES = [
    "mean_degree", "std_degree", "max_degree", "frac_degree_2",
    "frac_degree_3", "frac_degree_4plus",
    "mean_cycle_length", "frac_local_cycles", "cycle_space_dim", "cycle_density",
    "diameter", "avg_clustering", "avg_path_length", "algebraic_connectivity",
    "min_vertex_cut", "avg_vertex_cut", "treewidth_width",
    "mean_path_coverage", "path_length_ratio",
]


# --------------------------------------------------------------------------
# timeout helper (main thread only)
# --------------------------------------------------------------------------
class _TO(Exception):
    pass


def _alarm(signum, frame):
    raise _TO()


def with_timeout(fn, secs=15, default=float("nan")):
    old = signal.signal(signal.SIGALRM, _alarm)
    signal.alarm(secs)
    try:
        return fn()
    except Exception:
        return default
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old)


# --------------------------------------------------------------------------
# η + features for one instance
# --------------------------------------------------------------------------
def compute_record(G, S, T, gtype, meta=None):
    cycles = base.fundamental_cycles(G, S)
    if not cycles:
        return None
    try:
        P = base.seq_to_edges(nx.astar_path(G, S, T, heuristic=lambda a, b: 0))
    except Exception:
        return None
    # η = valid / touched (touched = cycles with overlap>=1)
    valid = 0
    touched = 0
    for c in cycles:
        ov = bool(c & P)
        if ov:
            touched += 1
        cand = base.xor(P, set(c))
        ok, _ = base.check_path(cand, S, T)
        if ok:
            valid += 1
    if touched == 0:
        return None
    eta = valid / touched
    coverage = touched / len(cycles)

    # GROUP A — degrees
    degs = np.array([d for _, d in G.degree()])
    V = G.number_of_nodes()
    Ecnt = G.number_of_edges()
    rec = {
        "gtype": gtype, "V": V, "E": Ecnt,
        "eta": eta, "coverage": coverage, "compat": valid / len(cycles),
        "n_cycles": len(cycles), "touched": touched, "valid": valid,
        "mean_degree": float(degs.mean()),
        "std_degree": float(degs.std()),
        "max_degree": int(degs.max()),
        "frac_degree_2": float(np.mean(degs == 2)),
        "frac_degree_3": float(np.mean(degs == 3)),
        "frac_degree_4plus": float(np.mean(degs >= 4)),
    }
    # GROUP B — cycle structure
    clen = np.array([len(c) for c in cycles])
    rec["mean_cycle_length"] = float(clen.mean())
    rec["frac_local_cycles"] = float(np.mean(clen <= 8))
    rec["cycle_space_dim"] = Ecnt - V + 1
    rec["cycle_density"] = (Ecnt - V + 1) / V
    # GROUP C — global topology (timeouts)
    rec["diameter"] = float(with_timeout(lambda: nx.diameter(G), 15))
    rec["avg_clustering"] = float(with_timeout(lambda: nx.average_clustering(G), 15))
    rec["avg_path_length"] = float(with_timeout(
        lambda: nx.average_shortest_path_length(G), 15))
    rec["algebraic_connectivity"] = float(with_timeout(
        lambda: nx.algebraic_connectivity(G), 15))
    # GROUP D — bottleneck geometry
    rec["min_vertex_cut"] = float(with_timeout(
        lambda: nx.node_connectivity(G, S, T), 15))
    def _avg_cut():
        nodes = list(G.nodes())
        rng = np.random.default_rng(0)
        vals = []
        for _ in range(10):
            a, b = nodes[int(rng.integers(0, len(nodes)))], nodes[int(rng.integers(0, len(nodes)))]
            if a != b:
                vals.append(nx.node_connectivity(G, a, b))
        return float(np.mean(vals)) if vals else float("nan")
    rec["avg_vertex_cut"] = float(with_timeout(_avg_cut, 15))
    rec["treewidth_width"] = float(with_timeout(
        lambda: nx.algorithms.approximation.treewidth_min_degree(G)[0], 15))
    # GROUP E — XOR specific
    rec["mean_path_coverage"] = coverage
    diam = rec["diameter"]
    rec["path_length_ratio"] = (len(P) / diam) if diam and not np.isnan(diam) and diam > 0 else float("nan")
    if meta:
        rec.update(meta)
    return rec


# --------------------------------------------------------------------------
# dataset construction (regenerate everything deterministically)
# --------------------------------------------------------------------------
def build_dataset():
    records = []
    done = 0

    def log():
        nonlocal done
        done += 1
        if done % 50 == 0:
            print(f"  ... {done} instances computed", flush=True)
            json.dump(records, open(CACHE, "w"))

    # TRAIN: grid size × density
    for size in [10, 20, 30, 50]:
        for dens in [0.10, 0.20, 0.30, 0.40, 0.50]:
            for s in range(30):
                inst = vb.gen_grid(size, dens, s)
                if inst is None:
                    continue
                G, S, T, _ = inst
                r = compute_record(G, S, T, "grid",
                                   {"split": "train", "size": size, "density": dens})
                if r:
                    records.append(r); log()
    # TRAIN: non-grid planar types
    for s in range(20):
        for tname, gen in [("delaunay", vb.make_delaunay),
                           ("gabriel", vb.make_gabriel),
                           ("random_planar", vb.make_random_planar)]:
            try:
                res = gen(s)
                G, S, T = res[0], res[2], res[3]
            except Exception:
                continue
            if G is None or S not in G or T not in G or not nx.has_path(G, S, T):
                continue
            r = compute_record(G, S, T, tname, {"split": "train", "seed": s})
            if r:
                records.append(r); log()

    # HELD-OUT: path_selection (20×20) and cycle_filter (20/30/50) grids
    ps = os.path.join(RES, "path_selection_experiment.json")
    if os.path.exists(ps):
        for p in json.load(open(ps))["per_instance"]:
            used = p["_meta"]["used_seed"]
            G, S, T, _ = base.gen_graph(20, used)
            r = compute_record(G, S, T, "grid",
                               {"split": "heldout", "source": "path_selection",
                                "size": 20, "density": 0.30})
            if r:
                records.append(r); log()
    cf = os.path.join(RES, "cycle_filter_experiment.json")
    if os.path.exists(cf):
        for p in json.load(open(cf))["per_instance"]:
            G, S, T, _ = base.gen_graph(p["grid_size"], p["used_seed"])
            r = compute_record(G, S, T, "grid",
                               {"split": "heldout", "source": "cycle_filter",
                                "size": p["grid_size"], "density": 0.30})
            if r:
                records.append(r); log()

    json.dump(records, open(CACHE, "w"))
    print(f"Built {len(records)} records → {CACHE}", flush=True)
    return records


# --------------------------------------------------------------------------
# cleaning
# --------------------------------------------------------------------------
def clean(records):
    out = []
    for r in records:
        e = r.get("eta")
        if e is None or np.isnan(e) or e > 1.05:
            continue
        out.append(r)
    return out


def matrix(records, feats):
    X = np.array([[r.get(f, np.nan) for f in feats] for r in records], dtype=float)
    y = np.array([r["eta"] for r in records], dtype=float)
    return X, y


# --------------------------------------------------------------------------
# manual VIF
# --------------------------------------------------------------------------
def vif(X, names):
    out = {}
    Xi = SimpleImputer(strategy="median").fit_transform(X)
    for i in range(Xi.shape[1]):
        others = np.delete(Xi, i, axis=1)
        if np.std(Xi[:, i]) == 0:
            out[names[i]] = float("inf"); continue
        lr = LinearRegression().fit(others, Xi[:, i])
        r2 = lr.score(others, Xi[:, i])
        out[names[i]] = float(1.0 / (1.0 - r2)) if r2 < 1 - 1e-9 else float("inf")
    return out


# --------------------------------------------------------------------------
# analysis
# --------------------------------------------------------------------------
def analyze():
    records = clean(json.load(open(CACHE)))
    train = [r for r in records if r.get("split") == "train"]
    held = [r for r in records if r.get("split") == "heldout"]
    print(f"Loaded {len(records)} valid η records "
          f"(train={len(train)}, heldout={len(held)}).")
    if len(train) < 100:
        print("WARNING: <100 training instances — results less reliable.")

    y_all = np.array([r["eta"] for r in train])

    # ---- STEP 2: correlations ----
    corr = []
    for f in FEATURES:
        x = np.array([r.get(f, np.nan) for r in train], dtype=float)
        mask = ~np.isnan(x)
        if mask.sum() < 10 or np.std(x[mask]) == 0:
            corr.append((f, float("nan"), float("nan"), float("nan"), int(mask.sum())))
            continue
        pr, pp = sstats.pearsonr(x[mask], y_all[mask])
        sr, _ = sstats.spearmanr(x[mask], y_all[mask])
        corr.append((f, float(pr), float(pp), float(sr), int(mask.sum())))
    corr_sorted = sorted(corr, key=lambda c: -abs(c[1]) if not np.isnan(c[1]) else 0)

    top5 = [c[0] for c in corr_sorted[:5]]
    strong = [c[0] for c in corr_sorted if not np.isnan(c[1]) and abs(c[1]) > 0.5]
    useless = [c[0] for c in corr_sorted if not np.isnan(c[1]) and abs(c[1]) < 0.1]

    Xtop, _ = matrix(train, top5)
    vif_top = vif(Xtop, top5)

    # ---- STEP 3: models (5-fold CV R² + 80/20 test RMSE) ----
    def cv_and_test(feats, poly=False, rf=False):
        X, y = matrix(train, feats)
        imp = SimpleImputer(strategy="median")
        if rf:
            model = make_pipeline(imp, RandomForestRegressor(
                n_estimators=100, random_state=0))
        elif poly:
            model = make_pipeline(imp, StandardScaler(),
                                  PolynomialFeatures(degree=2, include_bias=False),
                                  LinearRegression())
        else:
            model = make_pipeline(imp, StandardScaler(), LinearRegression())
        cv = cross_val_score(model, X, y, cv=5, scoring="r2")
        Xtr, Xte, ytr, yte = train_test_split(X, y, test_size=0.2, random_state=0)
        model.fit(Xtr, ytr)
        rmse = float(np.sqrt(mean_squared_error(yte, model.predict(Xte))))
        return float(cv.mean()), float(cv.std()), rmse, model

    # MODEL 1 (spec hypothesis): single feature frac_degree_2
    m1 = cv_and_test(["frac_degree_2"])
    # Also: single-feature CV for every top correlate, to pick the feature that
    # actually predicts best (CRITICAL: do not assume frac_degree_2 wins).
    single_cv = {}
    for f in dict.fromkeys(top5 + ["frac_degree_2", "mean_cycle_length"]):
        single_cv[f] = cv_and_test([f])[0]
    best_feat = max(single_cv, key=single_cv.get)   # best by CV, not by |r|
    m_interp = cv_and_test([best_feat])
    # best 2-feature pair among top5
    best_pair, best_pair_score = None, -np.inf
    import itertools
    for a, b in itertools.combinations(top5, 2):
        sc = cv_and_test([a, b])[0]
        if sc > best_pair_score:
            best_pair_score, best_pair = sc, (a, b)
    m2 = cv_and_test(list(best_pair))
    m3 = cv_and_test(top5)
    top3 = top5[:3]
    m4 = cv_and_test(top3, poly=True)
    m5 = cv_and_test(FEATURES, rf=True)

    # RF importances
    Xrf, yrf = matrix(train, FEATURES)
    rf = make_pipeline(SimpleImputer(strategy="median"),
                       RandomForestRegressor(n_estimators=100, random_state=0)).fit(Xrf, yrf)
    importances = dict(sorted(zip(FEATURES, rf.named_steps["randomforestregressor"]
                                  .feature_importances_), key=lambda z: -z[1]))

    # ---- MODEL 6: interpretable formula ----
    # fit several functional forms on the single best feature; pick best CV R²
    xb = np.array([r.get(best_feat, np.nan) for r in train], dtype=float)
    mask = ~np.isnan(xb)
    xb_, y_ = xb[mask], y_all[mask]
    forms = {}
    # linear
    lr = LinearRegression().fit(xb_.reshape(-1, 1), y_)
    forms["linear"] = {"r2": float(r2_score(y_, lr.predict(xb_.reshape(-1, 1)))),
                       "expr": f"{lr.intercept_:.4f} + {lr.coef_[0]:.4f}*{best_feat}",
                       "params": [float(lr.intercept_), float(lr.coef_[0])]}
    # saturating: 1 - exp(-c*x)  (only if feature >=0)
    if np.all(xb_ >= 0):
        try:
            from scipy.optimize import curve_fit
            def sat(x, a, c):
                return a * (1 - np.exp(-c * x))
            popt, _ = curve_fit(sat, xb_, y_, p0=[1.0, 1.0], maxfev=5000)
            forms["saturating"] = {
                "r2": float(r2_score(y_, sat(xb_, *popt))),
                "expr": f"{popt[0]:.4f}*(1 - exp(-{popt[1]:.4f}*{best_feat}))",
                "params": [float(popt[0]), float(popt[1])]}
        except Exception:
            pass
    best_form_name = max(forms, key=lambda k: forms[k]["r2"])
    # interpretable formula uses the BEST-CV single feature (= best_feat)
    interp = {"feature": best_feat, "forms": forms, "best_form": best_form_name,
              "cv_r2": m_interp[0], "cv_r2_std": m_interp[1], "test_rmse": m_interp[2]}

    # ---- assemble & decide claim strength ----
    candidates = {
        "frac_degree_2_hypothesis": m1, "best_single_feature": m_interp,
        "best_2feature": m2, "five_feature_linear": m3,
        "poly_deg2_top3": m4, "random_forest": m5,
    }
    best_model_name = max(candidates, key=lambda k: candidates[k][0])
    best_cv = candidates[best_model_name][0]
    # claim strength = best achievable INTERPRETABLE model (single or 2-feature)
    claim_cv = max(m_interp[0], m2[0])
    strength = ("STRONG" if claim_cv > 0.70 else
                "MODERATE" if claim_cv >= 0.50 else "WEAK")

    results = {
        "n_train": len(train), "n_heldout": len(held),
        "correlations": [{"feature": f, "pearson_r": r, "p_value": p,
                          "spearman_r": s, "n": n} for (f, r, p, s, n) in corr_sorted],
        "top5": top5, "strong_features": strong, "useless_features": useless,
        "vif_top5": vif_top,
        "single_feature_cv": single_cv,
        "models": {
            "frac_degree_2_hypothesis": {"features": ["frac_degree_2"], "cv_r2": m1[0],
                                         "cv_r2_std": m1[1], "test_rmse": m1[2]},
            "best_single_feature": {"features": [best_feat], "cv_r2": m_interp[0],
                                    "cv_r2_std": m_interp[1], "test_rmse": m_interp[2]},
            "best_2feature": {"features": list(best_pair), "cv_r2": m2[0],
                              "cv_r2_std": m2[1], "test_rmse": m2[2]},
            "five_feature_linear": {"features": top5, "cv_r2": m3[0],
                                    "cv_r2_std": m3[1], "test_rmse": m3[2]},
            "poly_deg2_top3": {"features": top3, "cv_r2": m4[0],
                               "cv_r2_std": m4[1], "test_rmse": m4[2]},
            "random_forest": {"features": FEATURES, "cv_r2": m5[0],
                              "cv_r2_std": m5[1], "test_rmse": m5[2]},
        },
        "rf_importances": importances,
        "interpretable_formula": interp,
        "best_model": best_model_name, "best_cv_r2": best_cv,
        "claim_cv_r2": claim_cv, "claim_strength": strength,
    }
    json.dump(results, open(OUT_MAIN, "w"), indent=2)

    # ---- STEP 5: held-out validation with interpretable single-feature model
    held_valid = validate_heldout(train, held, best_feat, m1[3] if False else None)
    json.dump(held_valid, open(OUT_VALID, "w"), indent=2)

    # ---- STEP 4: plots
    save_plots(train, best_feat, rf, importances)

    write_md(results, held_valid)
    report(results, held_valid)
    return results


def validate_heldout(train, held, best_feat, _unused):
    # refit interpretable single-feature linear on ALL train, predict held
    Xtr, ytr = matrix(train, [best_feat])
    model = make_pipeline(SimpleImputer(strategy="median"),
                          StandardScaler(), LinearRegression()).fit(Xtr, ytr)
    if not held:
        return {"note": "no held-out instances"}
    Xh, yh = matrix(held, [best_feat])
    pred = model.predict(Xh)
    err = np.abs(pred - yh)
    return {
        "feature": best_feat, "n_heldout": len(held),
        "mae": float(err.mean()), "rmse": float(np.sqrt(((pred - yh) ** 2).mean())),
        "max_error": float(err.max()),
        "r2_heldout": float(r2_score(yh, pred)),
        "mean_eta_heldout": float(yh.mean()),
    }


def save_plots(train, best_feat, rf_pipe, importances):
    def col(f):
        return [r.get(f, None) for r in train]
    eta = [r["eta"] for r in train]
    gt = [r["gtype"] for r in train]
    # model line for PLOT A (single-feature linear)
    x = np.array([r.get(best_feat, np.nan) for r in train], dtype=float)
    m = ~np.isnan(x)
    lr = LinearRegression().fit(x[m].reshape(-1, 1), np.array(eta)[m])
    xs = np.linspace(np.nanmin(x), np.nanmax(x), 50)
    plots = {
        "A": {"plot": "A", "x": col(best_feat), "y": eta, "color": gt,
              "x_label": best_feat, "y_label": "eta",
              "model_line": {"x": xs.tolist(),
                             "y_pred": lr.predict(xs.reshape(-1, 1)).tolist()}},
        "B": {"plot": "B", "x": col("mean_degree"), "y": eta, "color": gt,
              "x_label": "mean_degree", "y_label": "eta"},
        "D": {"plot": "D", "feature": list(importances.keys()),
              "importance": list(importances.values()),
              "x_label": "feature", "y_label": "RF importance"},
        "E": {"plot": "E", "actual": eta,
              "predicted": rf_pipe.predict(matrix(train, FEATURES)[0]).tolist(),
              "x_label": "actual eta", "y_label": "predicted eta (RF)"},
    }
    json.dump(plots, open(OUT_PLOTS, "w"), indent=2)


def report(res, hv):
    print("\n=== FEATURE CORRELATIONS WITH η ===")
    print(f"{'rank':>4}  {'feature':<24}{'pearson':>9}{'spearman':>10}{'p':>10}")
    for i, c in enumerate(res["correlations"], 1):
        if np.isnan(c["pearson_r"]):
            continue
        print(f"{i:>4}  {c['feature']:<24}{c['pearson_r']:+9.3f}"
              f"{c['spearman_r']:+10.3f}{c['p_value']:10.1e}")
    print(f"\nTop 5: {res['top5']}")
    print(f"|r|>0.5: {res['strong_features']}")
    print(f"|r|<0.1 (useless): {res['useless_features']}")
    print(f"VIF (top5): { {k: round(v,1) for k,v in res['vif_top5'].items()} }")

    print("\n=== η CHARACTERIZATION RESULTS ===")
    top = res["correlations"][0]
    sec = res["correlations"][1]
    print(f"Top predictor:    {top['feature']} (r={top['pearson_r']:+.3f})")
    print(f"Second predictor: {sec['feature']} (r={sec['pearson_r']:+.3f})")
    print(f"\n[hypothesis check] frac_degree_2 single-feature CV R² = "
          f"{res['models']['frac_degree_2_hypothesis']['cv_r2']:+.3f} "
          f"→ {'CONFIRMED as strong predictor' if res['models']['frac_degree_2_hypothesis']['cv_r2'] > 0.5 else 'NOT the dominant predictor (hypothesis falsified)'}")
    print("\nModel comparison (CV R²):")
    M = res["models"]
    print(f"  frac_degree_2 (hypothesis): R²={M['frac_degree_2_hypothesis']['cv_r2']:.3f}, "
          f"RMSE={M['frac_degree_2_hypothesis']['test_rmse']:.3f}")
    print(f"  best single ({M['best_single_feature']['features'][0]}): "
          f"R²={M['best_single_feature']['cv_r2']:.3f}, RMSE={M['best_single_feature']['test_rmse']:.3f}")
    print(f"  Best 2-feature {M['best_2feature']['features']}: "
          f"R²={M['best_2feature']['cv_r2']:.3f}, RMSE={M['best_2feature']['test_rmse']:.3f}")
    print(f"  5-feature linear: R²={M['five_feature_linear']['cv_r2']:.3f}, "
          f"RMSE={M['five_feature_linear']['test_rmse']:.3f}")
    print(f"  Polynomial deg-2: R²={M['poly_deg2_top3']['cv_r2']:.3f}, "
          f"RMSE={M['poly_deg2_top3']['test_rmse']:.3f}")
    print(f"  Random forest:    R²={M['random_forest']['cv_r2']:.3f}, "
          f"RMSE={M['random_forest']['test_rmse']:.3f}")
    int = res["interpretable_formula"]
    bf = int["forms"][int["best_form"]]
    print(f"\nBest interpretable formula ({int['best_form']}): "
          f"η ≈ {bf['expr']}  [single-feat CV R²={int['cv_r2']:.3f}]")
    print(f"\nHeld-out validation: MAE={hv.get('mae', float('nan')):.3f}, "
          f"RMSE={hv.get('rmse', float('nan')):.3f}, "
          f"max_error={hv.get('max_error', float('nan')):.3f}, "
          f"R²_heldout={hv.get('r2_heldout', float('nan')):.3f}")
    print(f"\nClaim strength: {res['claim_strength']} "
          f"(claim CV R²={res['claim_cv_r2']:.3f}; best model {res['best_model']} "
          f"R²={res['best_cv_r2']:.3f})")


def write_md(res, hv):
    L = ["# η characterization — crossing efficiency vs graph structure\n",
         f"Dataset: {res['n_train']} training + {res['n_heldout']} held-out "
         "(graphs regenerated deterministically; η = valid_XOR / cycles_touching_path).\n"]
    top = res["correlations"][0]
    L.append(f"## Top predictors\n")
    L.append("| rank | feature | Pearson r | Spearman | p |")
    L.append("|---:|---|---:|---:|---:|")
    for i, c in enumerate(res["correlations"][:10], 1):
        if np.isnan(c["pearson_r"]):
            continue
        L.append(f"| {i} | {c['feature']} | {c['pearson_r']:+.3f} "
                 f"| {c['spearman_r']:+.3f} | {c['p_value']:.1e} |")
    L.append(f"\n- |r|>0.5: {res['strong_features']}")
    L.append(f"- |r|<0.1 (useless): {res['useless_features']}\n")

    M = res["models"]
    L.append("## Model comparison (5-fold CV R², 20% test RMSE)\n")
    L.append("| model | features | CV R² | test RMSE |")
    L.append("|---|---|---:|---:|")
    L.append(f"| frac_degree_2 (hypothesis) | frac_degree_2 | {M['frac_degree_2_hypothesis']['cv_r2']:.3f} | {M['frac_degree_2_hypothesis']['test_rmse']:.3f} |")
    L.append(f"| best single | {M['best_single_feature']['features'][0]} | {M['best_single_feature']['cv_r2']:.3f} | {M['best_single_feature']['test_rmse']:.3f} |")
    L.append(f"| best-2 | {','.join(M['best_2feature']['features'])} | {M['best_2feature']['cv_r2']:.3f} | {M['best_2feature']['test_rmse']:.3f} |")
    L.append(f"| 5-feature | top5 | {M['five_feature_linear']['cv_r2']:.3f} | {M['five_feature_linear']['test_rmse']:.3f} |")
    L.append(f"| poly deg-2 | top3 | {M['poly_deg2_top3']['cv_r2']:.3f} | {M['poly_deg2_top3']['test_rmse']:.3f} |")
    L.append(f"| random forest | all | {M['random_forest']['cv_r2']:.3f} | {M['random_forest']['test_rmse']:.3f} |")

    int = res["interpretable_formula"]
    bf = int["forms"][int["best_form"]]
    L.append(f"\n## Best interpretable formula\n")
    L.append(f"η(G) ≈ **{bf['expr']}**  [single-feature CV R²={int['cv_r2']:.3f}, "
             f"RMSE={int['test_rmse']:.3f}]\n")
    L.append(f"Held-out validation ({hv.get('n_heldout','?')} grids): "
             f"MAE={hv.get('mae', float('nan')):.3f}, RMSE={hv.get('rmse', float('nan')):.3f}, "
             f"max_error={hv.get('max_error', float('nan')):.3f}.\n")
    L.append("## RF feature importances (top 6)\n")
    for k, v in list(res["rf_importances"].items())[:6]:
        L.append(f"- {k}: {v:.3f}")
    L.append("")

    st = res["claim_strength"]
    L.append(f"## Claim strength: **{st}** (CV R²={res['claim_cv_r2']:.3f})\n")
    L.append("### PAPER CLAIM (exact text)\n")
    tp = top["feature"]; r = top["pearson_r"]
    if st == "STRONG":
        L.append(f"> Crossing efficiency is well predicted by graph structure: "
                 f"`{bf['expr']}` (CV R²={int['cv_r2']:.2f}, held-out RMSE={hv.get('rmse',float('nan')):.2f}). "
                 f"The dominant predictor is `{tp}` (r={r:+.2f}).")
    elif st == "MODERATE":
        L.append(f"> Crossing efficiency η is *partially* predictable from graph "
                 f"structure. The dominant correlate is `{tp}` (Pearson r={r:+.2f}); "
                 f"a single-feature model gives η ≈ {bf['expr']} with CV R²={int['cv_r2']:.2f} "
                 f"(held-out RMSE={hv.get('rmse',float('nan')):.2f}). This is a rough estimate "
                 f"(±{hv.get('rmse',float('nan')):.2f}) suitable for order-of-magnitude "
                 f"reasoning, not precise prediction.")
    else:
        L.append(f"> η is **not** well predicted by simple structural features "
                 f"(best single-feature CV R²={int['cv_r2']:.2f}). The qualitative "
                 f"driver is real — the strongest correlate is `{tp}` (r={r:+.2f}) — "
                 f"but no simple formula captures η reliably; we report the "
                 f"qualitative trend only.")
    open(OUT_MD, "w").write("\n".join(L) + "\n")


def main():
    os.makedirs(RES, exist_ok=True)
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("build", "all"):
        build_dataset()
    if which in ("analyze", "all"):
        analyze()
    print(f"\nSaved: {OUT_MAIN}, {OUT_PLOTS}, {OUT_VALID}, {OUT_MD}")


if __name__ == "__main__":
    main()
