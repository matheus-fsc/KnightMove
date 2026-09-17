#!/usr/bin/env python3
"""
correlator.py
=============
Carrega batches gravados pelo runner, computa correlações de Pearson
entre arestas (vetorizado, numpy) e colapsa pares por órbita D₄.

Saídas:
  data/correlations/raw_correlations.json    # por par (start,end)
  data/correlations/orbit_correlations.json  # consolidado em órbitas D₄

Uso:
  python3 correlator.py                       # roda sobre data/
  python3 correlator.py --threshold -0.30     # só correlações abaixo
  python3 correlator.py --top 100             # top-N por par
"""

import argparse
import glob
import json
import os
import re
from collections import defaultdict

import numpy as np

import d4_orbits
from d4_orbits import (
    label_to_edge, edge_to_label, pair_d4_canonical, edge_orbits,
    d4_apply_edge,
)
from expand_d4 import expand_and_pool_by_canonical, canonical_pair


# ── carregamento ──────────────────────────────────────────────────────

def load_batches(samples_dir):
    paths = sorted(glob.glob(os.path.join(samples_dir, "batch_*.json")))
    out = []
    for p in paths:
        with open(p, "r", encoding="utf-8") as f:
            out.append(json.load(f))
    return out


# ── streaming: agrupa file paths por par canônico sem carregar tudo ───

_BATCH_RE = re.compile(r"batch_p(\d+)_b\d+\.json")

def group_paths_by_canonical(samples_dir, config_pairs):
    """
    Agrupa file paths por par canônico SEM ler signatures (rápido).
    Usa pair_idx do nome do arquivo + config_pairs para mapear → canonical.
    Retorna dict: { canonical_pair_tuple: list[file_path] }.
    """
    pair_canon = {}
    for i, (s, e) in enumerate(config_pairs):
        pair_canon[i] = canonical_pair(tuple(s), tuple(e))

    paths = sorted(glob.glob(os.path.join(samples_dir, "batch_*.json")))
    groups = defaultdict(list)
    for p in paths:
        m = _BATCH_RE.search(os.path.basename(p))
        if not m:
            continue
        pidx = int(m.group(1))
        if pidx in pair_canon:
            groups[pair_canon[pidx]].append(p)
    return groups


def _build_d4_remaps(edge_labels):
    """Pré-calcula 8 permutações de índices para expansão D₄ vetorizada."""
    edge_tuples = [label_to_edge(s) for s in edge_labels]
    label_idx = {edge_to_label(e): k for k, e in enumerate(edge_tuples)}
    n = len(edge_tuples)
    remaps = np.zeros((8, n), dtype=np.int32)
    for t in range(8):
        for i, e in enumerate(edge_tuples):
            remaps[t, i] = label_idx[edge_to_label(d4_apply_edge(e, t))]
    return remaps


def streaming_consolidate(file_paths, expand_d4=False, log=None):
    """
    Lê batches de disco em sequência, converte signatures direto para
    np.bool_, opcionalmente aplica as 8 transformações D₄ via permutações
    de coluna vetorizadas, e concatena tudo num único array numpy.

    Retorna (edges_labels: list[str], X: np.ndarray[bool] de shape
    (n_total_samples, n_edges)).

    Memória O(tamanho_do_pool); jamais segura todos os batches como dicts
    Python ao mesmo tempo (descarta após concatenar).
    """
    if not file_paths:
        return [], np.zeros((0, 0), dtype=np.bool_)

    with open(file_paths[0], "r", encoding="utf-8") as f:
        first = json.load(f)
    edges = first["edges"]
    n_edges = len(edges)
    remaps = _build_d4_remaps(edges) if expand_d4 else None

    parts = []
    n_orig = 0
    for k, path in enumerate(file_paths):
        with open(path, "r", encoding="utf-8") as f:
            b = json.load(f)
        if b["edges"] != edges:
            raise ValueError(f"ordem de arestas inconsistente em {path}")
        X = np.asarray(b["signatures"], dtype=np.bool_)
        n_orig += X.shape[0]
        if expand_d4:
            for t in range(8):
                new_X = np.empty_like(X)
                new_X[:, remaps[t]] = X
                parts.append(new_X)
        else:
            parts.append(X)
        del b
        if log and (k + 1) % 200 == 0:
            log(f"    streamed {k+1}/{len(file_paths)} batches "
                f"(orig={n_orig:,}, parts_in_mem={len(parts)})")

    pool = np.concatenate(parts, axis=0)
    parts.clear()
    return edges, pool


def consolidate(batches):
    """Concatena assinaturas de batches que compartilham a mesma ordem de arestas."""
    if not batches:
        return [], np.zeros((0, 0))
    edges = batches[0]["edges"]
    for b in batches[1:]:
        if b["edges"] != edges:
            raise ValueError("ordem de arestas inconsistente entre lotes")
    sigs = []
    for b in batches:
        sigs.extend(b["signatures"])
    X = np.array(sigs, dtype=np.float64)
    return edges, X


# ── correlação vetorizada ────────────────────────────────────────────

def pearson_negative_pairs(X, threshold):
    """
    Retorna lista [(r, i, j)] com r < threshold, ordenada do mais negativo.
    Considera apenas arestas com 0 < freq < 1 (i.e., não-constantes).

    Aceita X como np.bool_, np.uint8, ou np.float64. Promove para float64
    apenas para a operação aritmética (cria 1 cópia do tamanho de X).
    """
    n, d = X.shape
    if n < 2:
        return [], np.zeros(d)

    Xf = X.astype(np.float64) if X.dtype != np.float64 else X
    freq = Xf.mean(axis=0)
    sd = Xf.std(axis=0)
    live = sd > 0

    Z = np.zeros_like(Xf)
    Z[:, live] = (Xf[:, live] - freq[live]) / sd[live]
    C = (Z.T @ Z) / n

    iu, ju = np.triu_indices(d, k=1)
    rs = C[iu, ju]
    mask = (rs < threshold) & live[iu] & live[ju]
    pairs = sorted(zip(rs[mask].tolist(), iu[mask].tolist(), ju[mask].tolist()))
    return pairs, freq


# ── colapso em órbitas D₄ ────────────────────────────────────────────

def collapse_to_orbits(pairs, edge_labels):
    buckets = defaultdict(list)
    for r, i, j in pairs:
        ea = label_to_edge(edge_labels[i])
        eb = label_to_edge(edge_labels[j])
        canon = pair_d4_canonical(ea, eb)
        buckets[canon].append((r, ea, eb))

    orbits = []
    for canon, members in buckets.items():
        rs = [r for r, _, _ in members]
        ca, cb = canon
        shared = set(ca) & set(cb)
        orbits.append({
            "representative": [edge_to_label(ca), edge_to_label(cb)],
            "shared_vertex": (
                chr(ord('A') + next(iter(shared))[1])
                + str(d4_orbits.BOARD - next(iter(shared))[0])
            ) if shared else None,
            "size": len(members),
            "r_min": round(min(rs), 4),
            "r_max": round(max(rs), 4),
            "r_mean": round(sum(rs) / len(rs), 4),
            "members": [
                {"r": round(r, 4),
                 "edges": [edge_to_label(a), edge_to_label(b)]}
                for r, a, b in members
            ],
        })
    orbits.sort(key=lambda o: o["r_min"])
    return orbits


# ── relatório por par ────────────────────────────────────────────────

def analyze_X(label, edges, X, threshold, top, start, end):
    """Variante que aceita X numpy pré-construído (rota streaming)."""
    n, d = X.shape
    freq = X.mean(axis=0)  # bool → float auto-promotion
    n_mand = int((freq == 1.0).sum())
    n_imp = int((freq == 0.0).sum())
    n_live = int(((freq > 0) & (freq < 1)).sum())

    pairs, _ = pearson_negative_pairs(X, threshold)
    pairs = pairs[:top]
    orbits = collapse_to_orbits(pairs, edges)

    print(f"\n── {label} " + "─" * (66 - len(label)))
    print(f"  amostras            : {n:,}")
    print(f"  arestas analisadas  : {d}")
    print(f"  obrigatórias (f=1)  : {n_mand}")
    print(f"  impossíveis  (f=0)  : {n_imp}")
    print(f"  livres              : {n_live}")
    print(f"  correlações < {threshold:>5.2f} : {len(pairs)}")
    print(f"  órbitas D₄ distintas: {len(orbits)}")
    if pairs:
        compact = len(pairs) / max(len(orbits), 1)
        print(f"  compactação D₄      : {compact:.1f}×")
        print(f"  r mais negativo     : {pairs[0][0]:.4f}  "
              f"({edges[pairs[0][1]]} ↔ {edges[pairs[0][2]]})")

    return {
        "label": label,
        "start": list(start),
        "end": list(end),
        "n_samples": int(n),
        "n_edges": int(d),
        "n_mandatory": n_mand,
        "n_impossible": n_imp,
        "n_live": n_live,
        "mandatory_edges": [edges[i] for i, f in enumerate(freq) if f == 1.0],
        "impossible_edges": [edges[i] for i, f in enumerate(freq) if f == 0.0],
        "n_correlations": len(pairs),
        "n_orbits": len(orbits),
        "compaction_d4": (
            round(len(pairs) / len(orbits), 2) if orbits else None
        ),
        "top_correlations": [
            {"r": round(r, 4), "edge_i": edges[i], "edge_j": edges[j]}
            for r, i, j in pairs[:20]
        ],
        "orbits": orbits,
    }


def analyze_pair(label, batches, threshold, top):
    edges, X = consolidate(batches)
    n, d = X.shape
    freq = X.mean(axis=0)
    n_mand = int((freq == 1.0).sum())
    n_imp = int((freq == 0.0).sum())
    n_live = int(((freq > 0) & (freq < 1)).sum())

    pairs, _ = pearson_negative_pairs(X, threshold)
    pairs = pairs[:top]
    orbits = collapse_to_orbits(pairs, edges)

    print(f"\n── {label} " + "─" * (66 - len(label)))
    print(f"  amostras            : {n}")
    print(f"  arestas analisadas  : {d}")
    print(f"  obrigatórias (f=1)  : {n_mand}")
    print(f"  impossíveis  (f=0)  : {n_imp}")
    print(f"  livres              : {n_live}")
    print(f"  correlações < {threshold:>5.2f} : {len(pairs)}")
    print(f"  órbitas D₄ distintas: {len(orbits)}")
    if pairs:
        compact = len(pairs) / max(len(orbits), 1)
        print(f"  compactação D₄      : {compact:.1f}×")
        print(f"  r mais negativo     : {pairs[0][0]:.4f}  "
              f"({edges[pairs[0][1]]} ↔ {edges[pairs[0][2]]})")

    return {
        "label": label,
        "start": list(batches[0]["start"]),
        "end": list(batches[0]["end"]),
        "n_samples": n,
        "n_edges": d,
        "n_mandatory": n_mand,
        "n_impossible": n_imp,
        "n_live": n_live,
        "mandatory_edges": [edges[i] for i, f in enumerate(freq) if f == 1.0],
        "impossible_edges": [edges[i] for i, f in enumerate(freq) if f == 0.0],
        "n_correlations": len(pairs),
        "n_orbits": len(orbits),
        "compaction_d4": (
            round(len(pairs) / len(orbits), 2) if orbits else None
        ),
        "top_correlations": [
            {"r": round(r, 4), "edge_i": edges[i], "edge_j": edges[j]}
            for r, i, j in pairs[:20]
        ],
        "orbits": orbits,
    }


# ── main ──────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="data")
    ap.add_argument("--board", type=int, default=8,
                    help="tamanho do tabuleiro (default 8; use 6 para validação)")
    ap.add_argument("--threshold", type=float, default=-0.30)
    ap.add_argument("--top", type=int, default=100,
                    help="máx. correlações negativas a manter por par")
    ap.add_argument("--expand-d4", action="store_true",
                    help="expandir cada amostra pelas 8 transformações D₄ "
                         "antes de correlacionar (pool canônico)")
    args = ap.parse_args()

    d4_orbits.set_board(args.board)

    samples_dir = os.path.join(args.data_dir, "samples")
    out_dir = os.path.join(args.data_dir, "correlations")
    os.makedirs(out_dir, exist_ok=True)

    print("=" * 70)
    print(f"CORRELATOR {args.board}×{args.board}")
    print("=" * 70)

    paths = sorted(glob.glob(os.path.join(samples_dir, "batch_*.json")))
    print(f"\nLotes em disco: {len(paths)}  (de {samples_dir})")
    if not paths:
        print("Nenhum lote encontrado. Execute o runner primeiro.")
        return 1

    results = {}

    if args.expand_d4:
        # Rota STREAMING: agrupa file paths por par canônico via config.json
        # (sem ler signatures), depois processa um pool por vez.
        cfg_path = os.path.join(args.data_dir, "config.json")
        if not os.path.exists(cfg_path):
            print(f"ERRO: precisa de {cfg_path} para mapear pair_idx → canonical")
            return 1
        with open(cfg_path) as f:
            config = json.load(f)
        groups = group_paths_by_canonical(samples_dir, config["pairs"])
        print(f"Pools canônicos: {len(groups)}")
        print(f"Modo streaming: cada pool consolidado em numpy bool e "
              f"liberado antes do próximo (memória O(maior_pool)).")
        for canon, plist in sorted(groups.items()):
            s, e = canon
            label = f"{s[0]}{s[1]}_{e[0]}{e[1]}_canon"
            print(f"\n[pool {label}] {len(plist):,} batches no disco — streaming...")
            edges, X = streaming_consolidate(plist, expand_d4=True, log=print)
            mb_X = X.nbytes / (1024 * 1024)
            print(f"[pool {label}] X.shape={X.shape}  dtype={X.dtype}  ({mb_X:.0f} MB)")
            results[label] = analyze_X(label, edges, X, args.threshold, args.top,
                                        start=s, end=e)
            del X  # libera antes do próximo pool
    else:
        # Rota original (não-streaming): só viável para datasets pequenos
        batches = load_batches(samples_dir)
        by_pair = defaultdict(list)
        for b in batches:
            key = (tuple(b["start"]), tuple(b["end"]))
            by_pair[key].append(b)
        print(f"Pares (start, end) distintos: {len(by_pair)}")
        for (s, e), pair_batches in sorted(by_pair.items()):
            label = f"{s[0]}{s[1]}_{e[0]}{e[1]}"
            results[label] = analyze_pair(label, pair_batches,
                                           args.threshold, args.top)

    raw_path = os.path.join(out_dir, "raw_correlations.json")
    with open(raw_path, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)

    # consolidado por órbitas (cross-pair) — top-level summary
    summary = {
        "n_pairs_analyzed": len(results),
        "by_pair": {
            label: {
                "n_samples": r["n_samples"],
                "n_correlations": r["n_correlations"],
                "n_orbits": r["n_orbits"],
                "compaction_d4": r["compaction_d4"],
                "top_5_orbits": r["orbits"][:5],
            }
            for label, r in results.items()
        },
    }
    orbit_path = os.path.join(out_dir, "orbit_correlations.json")
    with open(orbit_path, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)

    print("\n" + "=" * 70)
    print("Saídas:")
    print(f"  {os.path.relpath(raw_path)}")
    print(f"  {os.path.relpath(orbit_path)}")
    print("=" * 70)
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
