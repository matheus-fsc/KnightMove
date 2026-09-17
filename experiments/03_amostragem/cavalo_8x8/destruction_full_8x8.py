#!/usr/bin/env python3
"""
destruction_full_8x8.py
========================
Análise de destruição de loops sobre TODAS as amostras (não subsample),
paralelizada via multiprocessing.Pool.

Algoritmo (por caminho, O(V × max_deg × α(V))):
  - reconstrói o caminho hamiltoniano da assinatura (walk O(V))
  - simula a remoção em ORDEM REVERSA: começa com grafo vazio e adiciona
    path[n-1], path[n-2], ..., path[0]. Mantém união-find sobre as
    posições do caminho (índices 0..n-1) com path-halving + union-by-rank
  - rastreia (V, E, C) ao longo dos n+1 estados; H₁ = E - V + C
  - Δ H₁ no passo forward t = H₁_history[n-t] - H₁_history[n-t-1]

Por que não GPU: união-find e BFS em grafos pequenos (V=64, E≤168) são
serial-friendly. Throughput obtido (~80k caminhos/s/core) processa todo
o dataset em ~5s/core × 16 cores ≈ alguns segundos. CUDA daria ~10× a
mais código pra ganho marginal.

Uso:
  ../.venv/bin/python destruction_full_8x8.py
  ../.venv/bin/python destruction_full_8x8.py --expand-d4   # × 8 paths
  ../.venv/bin/python destruction_full_8x8.py --workers 16
"""

import argparse
import glob
import json
import os
import time
from collections import defaultdict, Counter
from multiprocessing import Pool

BOARD = 8
KNIGHT_DELTAS = [(1,2),(2,1),(-1,2),(-2,1),(1,-2),(2,-1),(-1,-2),(-2,-1)]


# ── grafo (precomputado, herdado pelos workers via fork) ─────────────

def _build_adj():
    adj = {(r, c): [] for r in range(BOARD) for c in range(BOARD)}
    for r in range(BOARD):
        for c in range(BOARD):
            for dr, dc in KNIGHT_DELTAS:
                nr, nc = r + dr, c + dc
                if 0 <= nr < BOARD and 0 <= nc < BOARD:
                    adj[(r, c)].append((nr, nc))
    return adj

ADJ = _build_adj()
NODES = list(ADJ.keys())
DEGREE = {v: len(ADJ[v]) for v in NODES}


def _build_canonical_edges():
    s = set()
    for v in NODES:
        for u in ADJ[v]:
            e = (v, u) if v <= u else (u, v)
            s.add(e)
    return sorted(s)

EDGES = _build_canonical_edges()
N_EDGES = len(EDGES)


def _node_label(rc):
    r, c = rc
    return chr(ord('A') + c) + str(BOARD - r)


# ── D₄ ────────────────────────────────────────────────────────────────

def _d4_node(rc, t):
    r, c = rc
    N = BOARD - 1
    if t == 0: return (r, c)
    if t == 1: return (c, N - r)
    if t == 2: return (N - r, N - c)
    if t == 3: return (N - c, r)
    if t == 4: return (r, N - c)
    if t == 5: return (N - r, c)
    if t == 6: return (c, r)
    if t == 7: return (N - c, N - r)
    raise ValueError(t)

def _d4_edge(e, t):
    u, v = e
    u2, v2 = _d4_node(u, t), _d4_node(v, t)
    return (u2, v2) if u2 <= v2 else (v2, u2)

_EDGE_TO_IDX = {e: i for i, e in enumerate(EDGES)}
D4_REMAPS = [[_EDGE_TO_IDX[_d4_edge(e, t)] for e in EDGES] for t in range(8)]


# ── reconstrução + destruição ────────────────────────────────────────

def reconstruct_path(signature, start, end):
    """Walk caminho de start até end usando arestas ativas em signature."""
    adj_local = {v: [] for v in NODES}
    for i, active in enumerate(signature):
        if active:
            u, v = EDGES[i]
            adj_local[u].append(v)
            adj_local[v].append(u)

    path = [start]
    prev = None
    curr = start
    while curr != end:
        nbrs = adj_local[curr]
        nxt = nbrs[0] if (prev is None or nbrs[0] != prev) else nbrs[1]
        path.append(nxt)
        prev, curr = curr, nxt
    return path


def destruction_sequence_uf(path):
    """
    Δ H₁ por passo via união-find em ordem reversa.
    Custo: O(n × max_deg × α(n)) com path-halving + union-by-rank.
    """
    n = len(path)
    pos = {v: i for i, v in enumerate(path)}
    parent = list(range(n))
    rank = [0] * n
    is_alive = [False] * n

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    V = E = C = 0
    h1_history = [0] * (n + 1)

    for k in range(n):
        i = n - 1 - k
        v = path[i]
        is_alive[i] = True
        V += 1
        C += 1
        for u in ADJ[v]:
            j = pos.get(u, -1)
            if j == -1 or not is_alive[j]:
                continue
            E += 1
            ri = find(i)
            rj = find(j)
            if ri != rj:
                if rank[ri] < rank[rj]:
                    ri, rj = rj, ri
                parent[rj] = ri
                if rank[ri] == rank[rj]:
                    rank[ri] += 1
                C -= 1
        h1_history[k + 1] = max(0, E - V + C)

    return [h1_history[n - t] - h1_history[n - t - 1] for t in range(n)]


def transform_signature(sig, t):
    remap = D4_REMAPS[t]
    new_sig = [False] * N_EDGES
    for i, val in enumerate(sig):
        new_sig[remap[i]] = val
    return new_sig


# ── worker ───────────────────────────────────────────────────────────

def process_batch_file(args):
    path, expand_d4 = args
    with open(path, "r", encoding="utf-8") as f:
        b = json.load(f)
    s, e = tuple(b["start"]), tuple(b["end"])

    kill_n = defaultdict(int)
    kill_sum = defaultdict(int)
    kill_sumsq = defaultdict(int)
    kill_max = defaultdict(int)
    kill_dist = defaultdict(Counter)
    peak = Counter()
    n_processed = 0

    for sig in b["signatures"]:
        # Original
        path_seq = reconstruct_path(sig, s, e)
        seq = destruction_sequence_uf(path_seq)
        for v, k in zip(path_seq, seq):
            kill_n[v] += 1
            kill_sum[v] += k
            kill_sumsq[v] += k * k
            if k > kill_max[v]:
                kill_max[v] = k
            kill_dist[v][k] += 1
        if seq:
            peak[seq.index(max(seq))] += 1
        n_processed += 1

        # 7 imagens D₄ (opcional)
        if expand_d4:
            for t in range(1, 8):
                sig_t = transform_signature(sig, t)
                s_t = _d4_node(s, t)
                e_t = _d4_node(e, t)
                path_t = reconstruct_path(sig_t, s_t, e_t)
                seq_t = destruction_sequence_uf(path_t)
                for v, k in zip(path_t, seq_t):
                    kill_n[v] += 1
                    kill_sum[v] += k
                    kill_sumsq[v] += k * k
                    if k > kill_max[v]:
                        kill_max[v] = k
                    kill_dist[v][k] += 1
                if seq_t:
                    peak[seq_t.index(max(seq_t))] += 1
                n_processed += 1

    return {
        "n_processed": n_processed,
        "kill_n": dict(kill_n),
        "kill_sum": dict(kill_sum),
        "kill_sumsq": dict(kill_sumsq),
        "kill_max": dict(kill_max),
        "kill_dist": {v: dict(c) for v, c in kill_dist.items()},
        "peak": dict(peak),
    }


def merge_stats(g, b):
    g["n_processed"] += b["n_processed"]
    for v, n in b["kill_n"].items():
        g["kill_n"][v] = g["kill_n"].get(v, 0) + n
        g["kill_sum"][v] = g["kill_sum"].get(v, 0) + b["kill_sum"][v]
        g["kill_sumsq"][v] = g["kill_sumsq"].get(v, 0) + b["kill_sumsq"][v]
        g["kill_max"][v] = max(g["kill_max"].get(v, 0), b["kill_max"][v])
        gd = g["kill_dist"].setdefault(v, {})
        for k, c in b["kill_dist"][v].items():
            gd[k] = gd.get(k, 0) + c
    for k, c in b["peak"].items():
        g["peak"][k] = g["peak"].get(k, 0) + c


def heatmap_grid(values_dict, fmt):
    print("            " + "  ".join(chr(ord('A') + c).center(7) for c in range(BOARD)))
    for r in range(BOARD):
        rank_lbl = BOARD - r
        row = []
        for c in range(BOARD):
            v = values_dict.get((r, c), 0)
            row.append(fmt.format(v))
        print(f"  rank {rank_lbl}  " + "  ".join(row))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="data_parallel")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--expand-d4", action="store_true",
                    help="também processar as 7 imagens D₄ de cada caminho")
    ap.add_argument("--chunksize", type=int, default=4)
    args = ap.parse_args()

    samples_dir = os.path.join(args.data_dir, "samples")
    out_dir = os.path.join(args.data_dir, "destruction")
    os.makedirs(out_dir, exist_ok=True)

    paths = sorted(glob.glob(os.path.join(samples_dir, "batch_*.json")))
    print("=" * 70)
    print("DESTRUIÇÃO DE LOOPS — 8×8 (TODAS AS AMOSTRAS)")
    print("=" * 70)
    print(f"  data_dir     : {args.data_dir}")
    print(f"  lotes        : {len(paths)}")
    print(f"  workers      : {args.workers}")
    print(f"  expand_d4    : {args.expand_d4}")
    print(f"  chunksize    : {args.chunksize}")
    print()
    if not paths:
        print("Nenhum lote encontrado.")
        return 1

    tasks = [(p, args.expand_d4) for p in paths]

    global_stats = {
        "n_processed": 0,
        "kill_n": {}, "kill_sum": {}, "kill_sumsq": {},
        "kill_max": {}, "kill_dist": {}, "peak": {},
    }

    t0 = time.time()
    progress = 0
    with Pool(processes=args.workers) as pool:
        for result in pool.imap_unordered(process_batch_file, tasks,
                                           chunksize=args.chunksize):
            merge_stats(global_stats, result)
            progress += 1
            if progress % 200 == 0 or progress == len(tasks):
                elapsed = time.time() - t0
                rate_paths = global_stats["n_processed"] / max(elapsed, 1e-9)
                eta = (len(tasks) - progress) / max(progress / max(elapsed, 1e-9), 1e-9)
                print(f"  [{progress}/{len(tasks)}] "
                      f"paths={global_stats['n_processed']:>10,}  "
                      f"rate={rate_paths:>9,.0f}/s  "
                      f"eta={eta:>5.0f}s")

    elapsed = time.time() - t0
    rate = global_stats["n_processed"] / elapsed
    print()
    print(f"FIM: {global_stats['n_processed']:,} caminhos em {elapsed:.1f}s  "
          f"({rate:,.0f} paths/s, {rate/args.workers:,.0f}/s/core)")
    print()

    # ── relatório ────────────────────────────────────────────────────
    means = {v: global_stats["kill_sum"][v] / global_stats["kill_n"][v]
             for v in global_stats["kill_n"]}
    maxes = global_stats["kill_max"]

    print("── HEATMAP — DESTRUIÇÃO MÉDIA POR CASA ───────────────────────")
    print()
    heatmap_grid(means, "{:6.3f}")
    print()

    print("── HEATMAP — DESTRUIÇÃO MÁXIMA POR CASA ──────────────────────")
    print()
    heatmap_grid(maxes, "{:>5d}  ")
    print()

    by_deg = defaultdict(list)
    for v, m in means.items():
        by_deg[DEGREE[v]].append(m)

    print("── AGREGADO POR GRAU (sanity D₄: desvio→0 quando expand_d4) ──")
    print(f"  {'grau':>4}  {'casas':>5}  {'média':>8}  {'desvio':>8}  {'mín':>8}  {'máx':>8}")
    for deg in sorted(by_deg):
        vlist = by_deg[deg]
        m = sum(vlist) / len(vlist)
        sd = (sum((x - m) ** 2 for x in vlist) / len(vlist)) ** 0.5
        print(f"  {deg:>4}  {len(vlist):>5}  {m:>8.4f}  {sd:>8.4f}  "
              f"{min(vlist):>8.4f}  {max(vlist):>8.4f}")
    print()

    rows = sorted(((m, maxes[v], DEGREE[v], _node_label(v))
                   for v, m in means.items()), reverse=True)
    print("── RANKING POR DESTRUIÇÃO MÉDIA ──────────────────────────────")
    print(f"  {'#':>3}  {'casa':>4}  {'grau':>4}  {'média':>8}  {'max':>4}")
    for i, (m, mx, d, lbl) in enumerate(rows[:10]):
        print(f"  {i+1:>3}  {lbl:>4}  {d:>4}  {m:>8.4f}  {mx:>4}")
    print(f"  ── últimas 10 ──")
    for i, (m, mx, d, lbl) in enumerate(rows[-10:]):
        rk = len(rows) - 10 + i + 1
        print(f"  {rk:>3}  {lbl:>4}  {d:>4}  {m:>8.4f}  {mx:>4}")
    print()

    peak = global_stats["peak"]
    n_peaks = sum(peak.values())
    mean_peak = sum(int(k) * v for k, v in peak.items()) / max(n_peaks, 1)
    cum = 0
    median_step = None
    for k in sorted(peak):
        cum += peak[k]
        if cum >= n_peaks / 2 and median_step is None:
            median_step = int(k)

    print("── PASSO DE PICO DE DESTRUIÇÃO ───────────────────────────────")
    print(f"  N caminhos: {n_peaks:,}  passo médio: {mean_peak:.2f}  "
          f"mediano: {median_step}")
    print()
    print(f"  {'passo':>5}  {'N':>10}  {'%':>6}  {'cum%':>6}  histograma")
    cum = 0
    for k in sorted(peak):
        c = peak[k]
        cum += c
        pct = 100 * c / n_peaks
        cum_pct = 100 * cum / n_peaks
        bar = "█" * max(1, int(round(pct))) if pct >= 0.5 else ""
        print(f"  {int(k):>5d}  {c:>10,}  {pct:>5.1f}%  {cum_pct:>5.1f}%  {bar}")
    print()

    # ── salva JSON ────────────────────────────────────────────────────
    out_payload = {
        "data_dir": args.data_dir,
        "expand_d4": args.expand_d4,
        "workers": args.workers,
        "n_paths_processed": global_stats["n_processed"],
        "elapsed_s": round(elapsed, 1),
        "throughput_paths_per_s": round(rate, 0),
        "destruction_map": {
            _node_label(v): {
                "row_col": list(v),
                "degree": DEGREE[v],
                "n_appearances": global_stats["kill_n"][v],
                "mean": round(means[v], 6),
                "stddev": round(
                    max(0,
                        (global_stats["kill_sumsq"][v] / global_stats["kill_n"][v])
                        - means[v] ** 2
                    ) ** 0.5, 6),
                "max": maxes[v],
                "min": min(global_stats["kill_dist"][v].keys())
                       if global_stats["kill_dist"][v] else 0,
                "distribution": {str(k): c for k, c in
                                  sorted(global_stats["kill_dist"][v].items())},
            }
            for v in global_stats["kill_n"]
        },
        "by_degree": {
            str(deg): {
                "n_squares": len(vlist),
                "mean": round(sum(vlist) / len(vlist), 6),
                "stddev": round((sum((x - sum(vlist) / len(vlist)) ** 2
                                     for x in vlist) / len(vlist)) ** 0.5, 6),
                "min": round(min(vlist), 6),
                "max": round(max(vlist), 6),
            }
            for deg, vlist in by_deg.items()
        },
        "peak_step": {
            "distribution": {str(int(k)): int(v) for k, v in sorted(peak.items())},
            "mean": round(mean_peak, 4),
            "median": median_step,
        },
    }

    out_name = "destruction_full_d4.json" if args.expand_d4 else "destruction_full.json"
    out_path = os.path.join(out_dir, out_name)
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(out_payload, f, indent=2, ensure_ascii=False)
    print(f"Salvo: {out_path}")


if __name__ == "__main__":
    main()
