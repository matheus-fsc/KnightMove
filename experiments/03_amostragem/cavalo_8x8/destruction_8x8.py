#!/usr/bin/env python3
"""
destruction_8x8.py
==================
Análise de destruição de loops fundamentais por casa no 8×8 — análoga à
análise H₁ do `cavalo_loop_destruicao_6x6.py` (closed tours), mas aqui
sobre os caminhos hamiltonianos abertos amostrados pelo Z3.

Para cada amostra:
  1. reconstrói a sequência de vértices a partir da assinatura de arestas
     (caminho start → end, grau 2 nos interiores, grau 1 nos endpoints)
  2. remove vértices na ordem de visita do caminho
  3. computa Δ H₁ = H₁(G_alive_anterior) − H₁(G_alive_atual) em cada passo
     (usa a fórmula H₁ = E − V + C sobre o subgrafo induzido)

Agrega por casa:
  - distribuição de Δ H₁ por casa
  - média, max, mínimo
  - heatmap ASCII 8×8 (média e max)
  - desvio padrão dentro de cada classe de grau (sanity da D₄-invariância:
    no 6×6 fechado deu 0.0000)
  - distribuição do passo de pico de destruição

Subsample padrão: 2.000 caminhos por par canônico (16.000 totais).
Suficiente porque, como vimos no 6×6, dentro de uma órbita D₄ os números
são quase idênticos. `--samples-per-pool N` escala se quiser.

Uso:
  ../.venv/bin/python destruction_8x8.py --data-dir data_parallel
  ../.venv/bin/python destruction_8x8.py --samples-per-pool 5000
"""

import argparse
import glob
import json
import os
import random
import re
import time
from collections import defaultdict, Counter

from d4_orbits import label_to_edge, BOARD

KNIGHT_DELTAS = [(1,2),(2,1),(-1,2),(-2,1),(1,-2),(2,-1),(-1,-2),(-2,-1)]


# ── grafo do cavalo 8×8 (estável, fora dos hot loops) ────────────────

def build_adj():
    adj = {(r, c): [] for r in range(BOARD) for c in range(BOARD)}
    for r in range(BOARD):
        for c in range(BOARD):
            for dr, dc in KNIGHT_DELTAS:
                nr, nc = r + dr, c + dc
                if 0 <= nr < BOARD and 0 <= nc < BOARD:
                    adj[(r, c)].append((nr, nc))
    return adj

ADJ = build_adj()
NODES = list(ADJ.keys())
DEGREE = {v: len(ADJ[v]) for v in NODES}


def label(v):
    r, c = v
    return chr(ord('A') + c) + str(BOARD - r)


# ── H₁ residual ───────────────────────────────────────────────────────

def count_components(alive):
    if not alive:
        return 0
    visited = set()
    comps = 0
    for s in alive:
        if s in visited:
            continue
        comps += 1
        stack = [s]
        while stack:
            v = stack.pop()
            if v in visited:
                continue
            visited.add(v)
            for u in ADJ[v]:
                if u in alive and u not in visited:
                    stack.append(u)
    return comps


def cyclomatic_dim(alive):
    """H₁(G_alive) = E − V + C."""
    V = len(alive)
    if V == 0:
        return 0
    E = 0
    for v in alive:
        for u in ADJ[v]:
            if u in alive and u > v:
                E += 1
    C = count_components(alive)
    return max(0, E - V + C)


# ── reconstrução do caminho a partir da assinatura ───────────────────

def reconstruct_path(signature, edge_tuples, start, end):
    """
    Dado o vetor booleano `signature` (uma entrada por aresta em
    `edge_tuples`), constrói a adjacência do caminho ativo e caminha
    de start até end.
    """
    adj_local = defaultdict(list)
    for i, active in enumerate(signature):
        if active:
            u, v = edge_tuples[i]
            adj_local[u].append(v)
            adj_local[v].append(u)

    path = [start]
    prev = None
    curr = start
    while curr != end:
        nbrs = adj_local[curr]
        if prev is None:
            nxt = nbrs[0]
        else:
            nxt = nbrs[0] if nbrs[0] != prev else nbrs[1]
        path.append(nxt)
        prev = curr
        curr = nxt
    return path


def destruction_sequence(path):
    """Retorna lista de Δ H₁ em cada passo, removendo vértices na ordem."""
    alive = set(NODES)
    h1_prev = cyclomatic_dim(alive)
    seq = []
    for v in path:
        alive.discard(v)
        h1_curr = cyclomatic_dim(alive)
        seq.append(h1_prev - h1_curr)
        h1_prev = h1_curr
    return seq


# ── carregamento de assinaturas (subsample por pool canônico) ────────

_BATCH_RE = re.compile(r"batch_p(\d+)_b\d+\.json")

def collect_subsample(samples_dir, config_pairs, samples_per_pool, seed=42):
    """
    Para cada par canônico (lex-min sobre D₄ do par original), retorna
    uma lista [(signature, edge_tuples, start, end), ...] com até
    `samples_per_pool` amostras escolhidas aleatoriamente.

    Lemos UM batch de cada vez, sorteamos índices, descartamos.
    """
    from d4_orbits import d4_apply_node

    def canonical_pair(s, e):
        images = set()
        for t in range(8):
            st = d4_apply_node(tuple(s), t)
            et = d4_apply_node(tuple(e), t)
            images.add(tuple(sorted([st, et])))
        return min(images)

    pair_canon = {i: canonical_pair(s, e) for i, (s, e) in enumerate(config_pairs)}
    paths = sorted(glob.glob(os.path.join(samples_dir, "batch_*.json")))

    paths_by_canon = defaultdict(list)
    for p in paths:
        m = _BATCH_RE.search(os.path.basename(p))
        if not m:
            continue
        pidx = int(m.group(1))
        if pidx in pair_canon:
            paths_by_canon[pair_canon[pidx]].append(p)

    rng = random.Random(seed)
    subsamples = {}
    for canon, plist in paths_by_canon.items():
        # quantas amostras tirar de cada batch
        per_batch = max(1, samples_per_pool // len(plist))
        bucket = []
        for path in plist:
            with open(path, "r", encoding="utf-8") as f:
                b = json.load(f)
            sigs = b["signatures"]
            s, e = tuple(b["start"]), tuple(b["end"])
            edge_tuples = [label_to_edge(s_) for s_ in b["edges"]]
            n = len(sigs)
            k = min(per_batch, n)
            for idx in rng.sample(range(n), k):
                bucket.append((sigs[idx], edge_tuples, s, e))
        rng.shuffle(bucket)
        subsamples[canon] = bucket[:samples_per_pool]
    return subsamples


# ── relatório ─────────────────────────────────────────────────────────

def heatmap_grid(values_dict, fmt):
    """Imprime grade 8×8 com valores formatados."""
    print("            " + "  ".join(chr(ord('A') + c).center(7) for c in range(BOARD)))
    for r in range(BOARD):
        rank = BOARD - r
        row = []
        for c in range(BOARD):
            v = values_dict.get((r, c), 0.0)
            row.append(fmt.format(v))
        print(f"  rank {rank}  " + "  ".join(row))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="data_parallel")
    ap.add_argument("--samples-per-pool", type=int, default=2000,
                    help="caminhos por par canônico (default: 2000)")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    samples_dir = os.path.join(args.data_dir, "samples")
    cfg_path = os.path.join(args.data_dir, "config.json")
    out_dir = os.path.join(args.data_dir, "destruction")
    os.makedirs(out_dir, exist_ok=True)

    with open(cfg_path, "r", encoding="utf-8") as f:
        config = json.load(f)

    print("=" * 70)
    print("DESTRUIÇÃO DE LOOPS — CAVALO 8×8 (caminhos abertos)")
    print("=" * 70)
    print(f"  data_dir              : {args.data_dir}")
    print(f"  samples_per_pool      : {args.samples_per_pool}")
    print(f"  pares canônicos       : até {len(config['pairs'])}")
    print(f"  total nominal de paths: {args.samples_per_pool * len(config['pairs']):,}")
    print()

    t0 = time.time()
    subsamples = collect_subsample(samples_dir, config["pairs"],
                                    args.samples_per_pool, args.seed)
    n_paths_total = sum(len(v) for v in subsamples.values())
    print(f"Coleta concluída em {time.time()-t0:.1f}s. "
          f"Pools: {len(subsamples)}, paths totais: {n_paths_total:,}.")
    print()

    # acumuladores globais
    global_kills_by_square = defaultdict(list)
    peak_steps = Counter()
    paths_processed = 0

    # acumuladores por pool (pra ver invariância)
    per_pool_summary = []

    for canon, paths_list in sorted(subsamples.items()):
        s_can, e_can = canon
        pool_label = f"{s_can[0]}{s_can[1]}_{e_can[0]}{e_can[1]}"
        t_pool = time.time()
        pool_kills = defaultdict(list)
        n_pool = len(paths_list)

        for sig, edges, s, e in paths_list:
            path = reconstruct_path(sig, edges, s, e)
            seq = destruction_sequence(path)
            for v, k in zip(path, seq):
                pool_kills[v].append(k)
                global_kills_by_square[v].append(k)
            if seq:
                peak_steps[seq.index(max(seq))] += 1
        paths_processed += n_pool

        means = {v: sum(ks)/len(ks) for v, ks in pool_kills.items()}
        maxes = {v: max(ks) for v, ks in pool_kills.items()}
        per_pool_summary.append({
            "pool": pool_label,
            "n_paths": n_pool,
            "elapsed_s": round(time.time() - t_pool, 1),
            "max_destroyed_ever": max(maxes.values()) if maxes else 0,
        })
        print(f"  pool {pool_label}: {n_pool} paths em {time.time()-t_pool:.1f}s, "
              f"max kill global={max(maxes.values()) if maxes else 0}")

    print()
    print(f"Total de caminhos processados: {paths_processed:,}  "
          f"em {time.time()-t0:.1f}s")
    print()

    # ── 1. Heatmap de destruição média ────────────────────────────────
    means_global = {v: sum(ks)/len(ks) for v, ks in global_kills_by_square.items()}
    maxes_global = {v: max(ks) for v, ks in global_kills_by_square.items()}

    print("── 1. HEATMAP — DESTRUIÇÃO MÉDIA POR CASA ────────────────────")
    print()
    heatmap_grid(means_global, "{:6.2f}")
    print()

    print("── 2. HEATMAP — DESTRUIÇÃO MÁXIMA POR CASA ───────────────────")
    print()
    heatmap_grid({k: v for k, v in maxes_global.items()}, "{:>5d}  ")
    print()

    # ── 3. Aggregate by degree (sanity D₄) ────────────────────────────
    by_deg = defaultdict(list)
    for v, ks in global_kills_by_square.items():
        m = sum(ks) / len(ks)
        by_deg[DEGREE[v]].append(m)

    print("── 3. AGREGADO POR GRAU DO VÉRTICE ───────────────────────────")
    print(f"  {'grau':>4}  {'casas':>5}  {'média':>8}  {'desvio':>8}  {'mín':>8}  {'máx':>8}")
    print(f"  {'-'*4}  {'-'*5}  {'-'*8}  {'-'*8}  {'-'*8}  {'-'*8}")
    for deg in sorted(by_deg):
        v = by_deg[deg]
        m = sum(v) / len(v)
        sd = (sum((x - m) ** 2 for x in v) / len(v)) ** 0.5
        print(f"  {deg:>4}  {len(v):>5}  {m:>8.4f}  {sd:>8.4f}  {min(v):>8.4f}  {max(v):>8.4f}")
    print()

    # ── 4. Top / bottom squares ───────────────────────────────────────
    rows = sorted(((m, maxes_global[v], DEGREE[v], label(v))
                   for v, m in means_global.items()), reverse=True)
    print("── 4. RANKING DAS CASAS POR DESTRUIÇÃO MÉDIA ─────────────────")
    print(f"  {'#':>3}  {'casa':>4}  {'grau':>4}  {'média':>8}  {'max':>4}")
    print(f"  {'-'*3}  {'-'*4}  {'-'*4}  {'-'*8}  {'-'*4}")
    for i, (m, mx, d, lbl) in enumerate(rows[:10]):
        print(f"  {i+1:>3}  {lbl:>4}  {d:>4}  {m:>8.4f}  {mx:>4}")
    print(f"  ── últimas 10 ──")
    for i, (m, mx, d, lbl) in enumerate(rows[-10:]):
        rk = len(rows) - 10 + i + 1
        print(f"  {rk:>3}  {lbl:>4}  {d:>4}  {m:>8.4f}  {mx:>4}")
    print()

    # ── 5. Peak step distribution ─────────────────────────────────────
    n_peaks = sum(peak_steps.values())
    print("── 5. PASSO DE PICO DE DESTRUIÇÃO ────────────────────────────")
    print("    (em qual passo o caminho atinge seu kill máximo)")
    print(f"  {'passo':>5}  {'N':>7}  {'%':>6}  {'cum%':>6}  histograma")
    cum = 0
    for k in sorted(peak_steps):
        c = peak_steps[k]
        cum += c
        pct = 100 * c / n_peaks
        cum_pct = 100 * cum / n_peaks
        bar = "█" * max(1, int(round(pct))) if pct >= 0.5 else ""
        print(f"  {k:>5d}  {c:>7,}  {pct:>5.1f}%  {cum_pct:>5.1f}%  {bar}")

    mean_peak = sum(k * c for k, c in peak_steps.items()) / max(n_peaks, 1)
    median_step = None
    cum = 0
    for k in sorted(peak_steps):
        cum += peak_steps[k]
        if cum >= n_peaks / 2:
            median_step = k
            break
    print()
    print(f"  Passo de pico médio    : {mean_peak:.2f}")
    print(f"  Passo de pico mediano  : {median_step}")
    print()

    # ── 6. Salvar JSON ────────────────────────────────────────────────
    out = {
        "data_dir": args.data_dir,
        "samples_per_pool": args.samples_per_pool,
        "n_paths_processed": paths_processed,
        "per_pool_summary": per_pool_summary,
        "destruction_map": {
            label(v): {
                "row_col": list(v),
                "degree": DEGREE[v],
                "n_appearances": len(ks),
                "mean": round(sum(ks)/len(ks), 6),
                "max": max(ks),
                "min": min(ks),
                "distribution": dict(sorted(Counter(ks).items())),
            }
            for v, ks in global_kills_by_square.items()
        },
        "by_degree": {
            str(deg): {
                "n_squares": len(v),
                "mean": round(sum(v)/len(v), 6),
                "stddev": round((sum((x - sum(v)/len(v))**2 for x in v)/len(v))**0.5, 6),
                "min": round(min(v), 6),
                "max": round(max(v), 6),
            }
            for deg, v in by_deg.items()
        },
        "peak_step_distribution": dict(sorted(peak_steps.items())),
        "peak_step_mean": round(mean_peak, 4),
        "peak_step_median": median_step,
    }
    out_path = os.path.join(out_dir, "destruction_8x8.json")
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)

    print("=" * 70)
    print(f"Salvo: {out_path}")
    print("=" * 70)


if __name__ == "__main__":
    main()
