#!/usr/bin/env python3
"""
Face-cycle basis vs DFS basis for GF(2) XOR pathfinding.
========================================================

Tests the Mac Lane prediction: switching from the DFS spanning-tree cycle
basis (long, non-local back-edge cycles) to the PLANAR FACE basis (bounded
faces = geometrically local cycles, Mac Lane 2-basis) raises the single-XOR
compatibility rate, because a face cycle meets a path in a contiguous segment
(combinatorial Jordan curve).

Reuses the SAME 20×20 grid instances from pathfinding_xor_experiment
(reconstructed deterministically via used_seed). Does not regenerate maps.

Uso:
    ../../venv/bin/python face_basis_experiment.py
"""

import json
import os
import sys

import numpy as np
import networkx as nx

# localiza o projeto e reaproveita o pipeline existente
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PXE = os.path.join(ROOT, "pathfinding_xor_experiment")
sys.path.insert(0, PXE)
import pathfinding_xor as base   # noqa: E402

PRIOR_JSON = os.path.join(PXE, "results", "pathfinding_xor_experiment.json")
OUT_DIR = os.path.join(HERE, "face_basis_experiment", "results")
OUT_JSON = os.path.join(OUT_DIR, "face_vs_dfs.json")

GRID = 20   # apenas as instâncias 20×20 (conforme enunciado)


# --------------------------------------------------------------------------
# Base de faces de um grafo planar (via embedding do networkx)
# --------------------------------------------------------------------------
def face_basis(G):
    """Retorna (lista de ciclos-face como set de arestas, info). Exclui a face
    externa (não-limitada), identificada pela maior área absoluta."""
    ok, emb = nx.check_planarity(G)
    if not ok:
        raise RuntimeError("grafo não-planar (inesperado para grade)")
    seen = set()
    faces = []          # cada face = lista ordenada de nós
    for (v, w) in list(emb.edges()):
        if (v, w) in seen:
            continue
        nodes = emb.traverse_face(v, w, mark_half_edges=seen)
        faces.append(nodes)

    def face_edges(nodes):
        m = len(nodes)
        return {base.E(nodes[i], nodes[(i + 1) % m]) for i in range(m)}

    def shoelace(nodes):
        s = 0.0
        m = len(nodes)
        for i in range(m):
            x1, y1 = nodes[i]
            x2, y2 = nodes[(i + 1) % m]
            s += x1 * y2 - x2 * y1
        return abs(s) / 2.0

    if not faces:
        return [], {"n_faces_total": 0, "outer_face_len": 0}
    areas = [shoelace(f) for f in faces]
    outer = int(np.argmax(areas))          # face externa = maior área
    bounded = [face_edges(faces[i]) for i in range(len(faces)) if i != outer]
    # remove faces vazias/degeneradas e deduplica
    uniq = {}
    for fe in bounded:
        if fe:
            uniq[frozenset(fe)] = fe
    cycles = list(uniq.values())
    info = {
        "n_faces_total": len(faces),
        "outer_face_len": len(faces[outer]),
        "n_bounded_faces": len(cycles),
    }
    return cycles, info


# --------------------------------------------------------------------------
# Avalia uma base de ciclos: compat (raw e condicional) + diversidade
# --------------------------------------------------------------------------
def eval_basis(P, cycles, S, T):
    valid = []
    n_overlap_pos = 0          # ciclos com overlap >= 1
    n_valid_given_overlap = 0
    cyc_lens = []
    for c in cycles:
        cyc_lens.append(len(c))
        ov = len(c & P)
        cand = base.xor(P, set(c))
        ok, _ = base.check_path(cand, S, T)
        if ov >= 1:
            n_overlap_pos += 1
            if ok:
                n_valid_given_overlap += 1
        if ok:
            valid.append(cand)
    n = len(cycles)
    return {
        "n_cycles": n,
        "compat_raw": len(valid) / n if n else 0.0,
        "valid_paths": len(valid),
        "n_overlap_pos": n_overlap_pos,
        "compat_given_overlap": (n_valid_given_overlap / n_overlap_pos
                                 if n_overlap_pos else float("nan")),
        "diversity": base.diversity(valid),
        "mean_cycle_len": float(np.mean(cyc_lens)) if cyc_lens else 0.0,
    }


def main():
    prior = [r for r in json.load(open(PRIOR_JSON)) if r["grid_size"] == GRID]
    per = []
    for r in prior:
        used = r["used_seed"]
        G, S, T, seed_ok = base.gen_graph(GRID, used)
        assert seed_ok == used
        P = base.seq_to_edges(base.astar_path(G, S, T))

        dfs_cycles = base.fundamental_cycles(G, S)
        face_cycles, finfo = face_basis(G)

        dfs = eval_basis(P, dfs_cycles, S, T)
        face = eval_basis(P, face_cycles, S, T)
        per.append({"used_seed": used, "dfs": dfs, "face": face, "face_info": finfo})

    def agg(key, sub):
        vals = [p[key][sub] for p in per if not (
            isinstance(p[key][sub], float) and np.isnan(p[key][sub]))]
        return float(np.mean(vals)) if vals else float("nan")

    summary = {
        "grid_size": GRID, "n_instances": len(per),
        "dfs_compat_raw": agg("dfs", "compat_raw"),
        "face_compat_raw": agg("face", "compat_raw"),
        "dfs_compat_given_overlap": agg("dfs", "compat_given_overlap"),
        "face_compat_given_overlap": agg("face", "compat_given_overlap"),
        "dfs_diversity": agg("dfs", "diversity"),
        "face_diversity": agg("face", "diversity"),
        "dfs_mean_cycle_len": agg("dfs", "mean_cycle_len"),
        "face_mean_cycle_len": agg("face", "mean_cycle_len"),
        "dfs_n_cycles": agg("dfs", "n_cycles"),
        "face_n_cycles": agg("face", "n_cycles"),
    }
    improvement_raw = (summary["face_compat_raw"] / summary["dfs_compat_raw"]
                       if summary["dfs_compat_raw"] else float("inf"))
    improvement_cond = (summary["face_compat_given_overlap"]
                        / summary["dfs_compat_given_overlap"]
                        if summary["dfs_compat_given_overlap"] else float("inf"))
    summary["improvement_raw"] = improvement_raw
    summary["improvement_given_overlap"] = improvement_cond

    os.makedirs(OUT_DIR, exist_ok=True)
    json.dump({"summary": summary, "per_instance": per}, open(OUT_JSON, "w"), indent=2)

    print("=" * 50)
    print(f"Grid {GRID}×{GRID}, {len(per)} instances")
    print(f"DFS basis:   compat={summary['dfs_compat_raw']*100:5.1f}%  "
          f"diversity={summary['dfs_diversity']:.3f}  "
          f"(mean|C|={summary['dfs_mean_cycle_len']:.1f}, "
          f"n={summary['dfs_n_cycles']:.0f})")
    print(f"Face basis:  compat={summary['face_compat_raw']*100:5.1f}%  "
          f"diversity={summary['face_diversity']:.3f}  "
          f"(mean|C|={summary['face_mean_cycle_len']:.1f}, "
          f"n={summary['face_n_cycles']:.0f})")
    print(f"Improvement (raw compat):           {improvement_raw:.2f}×")
    print()
    print("Conditional on overlap >= 1 (where the Mac Lane / Jordan prediction")
    print("actually applies — does a touching cycle cross contiguously?):")
    print(f"DFS  compat|overlap>=1:  {summary['dfs_compat_given_overlap']*100:5.1f}%")
    print(f"Face compat|overlap>=1:  {summary['face_compat_given_overlap']*100:5.1f}%")
    print(f"Improvement (conditional):          {improvement_cond:.2f}×")
    print("=" * 50)
    hyp = improvement_raw > 2.0
    print(f"Hypothesis (raw >> 2×): {'CONFIRMED' if hyp else 'REFUTED'}")
    print(f"\nSaved to: {OUT_JSON}")


if __name__ == "__main__":
    main()
