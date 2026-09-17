#!/usr/bin/env python3
"""
triples_search.py
=================
Busca dirigida de triplas (e_a, e_b, e_c) com P(x_a ∧ x_b ∧ x_c) ≈ 0
nas 5.000 amostras 10×10 já coletadas.

Como C(288, 3) ≈ 3.9M é grande para força bruta com 5k amostras + RAM
limitada, usamos duas heurísticas dirigidas:

  1.1 — Estender pares correlatos (r<-0.5): para cada par e cada
        terceira aresta com freq>0.05, calcular P_triple.
  1.2 — Vértices de baixo grau (grau ≤ 4): triplas entre arestas
        incidentes ao mesmo vértice (já naturalmente excluídas pelo
        grau-2 se #incidentes>2, mas mantemos para sanidade).

Threshold:
  P_triple < 0.02 (em 5000 amostras = no máximo ~100 ocorrências)

Minimalidade:
  freq_pairs[ea,eb] > 0.02 e idem para os outros dois subpares
  (subpar ainda "vivo" → tripla não é redundância de par excluído)

Saídas:
  data/T_10x10.npy        — matriz 5000×288 reconstruída
  data/freq_singles.npy
  data/freq_pairs.npy
  data/triples_candidates.json
"""

import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
PARENT = ROOT.parent
SAMPLES = PARENT / "data" / "samples"
INV = PARENT / "data" / "invariants"

sys.path.insert(0, str(PARENT))
from graph_10x10 import build_graph, label  # noqa: E402

ADJ, EDGES = build_graph()


def edge_label(idx):
    u, v = EDGES[idx]
    return f"{label(u)}-{label(v)}"


def load_samples():
    files = sorted(SAMPLES.glob("tours_10x10_batch_*.npy"))
    if not files:
        raise FileNotFoundError(f"sem batches em {SAMPLES}")
    arrs = [np.load(f) for f in files]
    return np.concatenate(arrs, axis=0).astype(np.uint8)


def vertex_to_incident_edges():
    """vid -> lista de índices de arestas incidentes."""
    out = defaultdict(list)
    for i, (u, v) in enumerate(EDGES):
        out[u].append(i)
        out[v].append(i)
    return out


def main():
    DATA.mkdir(parents=True, exist_ok=True)

    print("Tarefa 0 — Reconstruir T e validar ...")
    T = load_samples()
    N, E = T.shape
    print(f"  T shape={T.shape}  (esperado N=5000, E=288)")
    assert N == 5000 and E == 288

    # cada tour usa 100 arestas (100 vértices, grau 2)
    sums = T.sum(axis=1)
    print(f"  arestas/tour: min={sums.min()} max={sums.max()} avg={sums.mean():.2f} "
          f"(esperado 100)")
    assert (sums == 100).all()

    freq = T.mean(axis=0).astype(np.float64)
    print(f"  freq.min={freq.min():.4f} freq.max={freq.max():.4f}")
    print(f"  freq média = {freq.mean():.4f} (esperado ~0.694)")
    # E * 2 / V = total incidências / V = grau médio. freq média de aresta:
    # 2*100/(2*E) = 100/E = 100/288 = 0.347 — note esperado do brief é
    # 2×100/288 = 0.694, mas isso só vale se contássemos cada aresta duas vezes.
    # O cálculo correto: 100 arestas/tour ÷ E = 100/288 ≈ 0.347.

    # rank GF(2)
    print("  computando rank GF(2) de T ...")

    def gf2_row_reduce(M):
        A = M.copy().astype(np.uint8)
        rows, cols = A.shape
        r = 0
        for c in range(cols):
            if r >= rows:
                break
            pivot = None
            for rr in range(r, rows):
                if A[rr, c]:
                    pivot = rr
                    break
            if pivot is None:
                continue
            if pivot != r:
                A[[r, pivot]] = A[[pivot, r]]
            for rr in range(rows):
                if rr != r and A[rr, c]:
                    A[rr] ^= A[r]
            r += 1
        return r

    rank = gf2_row_reduce(T)
    print(f"  rank(T) GF(2) = {rank} (esperado 186)")

    # ── freq pares ───────────────────────────────────────────────────
    print("\nCalculando freq_pairs ...")
    T32 = T.astype(np.int32)
    P2_count = T32.T @ T32  # contagem inteira
    P2 = P2_count.astype(np.float64) / N

    np.save(DATA / "T_10x10.npy", T)
    np.save(DATA / "freq_singles.npy", freq.astype(np.float32))
    np.save(DATA / "freq_pairs.npy", P2.astype(np.float32))

    # ── carrega pares correlatos ─────────────────────────────────────
    cand_path = INV / "candidate_clauses.json"
    cand_data = json.loads(cand_path.read_text())
    seed_pairs = []
    for c in cand_data["candidates"]:
        seed_pairs.append((int(c["edge_a_idx"]), int(c["edge_b_idx"]), c["r"]))
    print(f"\nPares semente (de candidate_clauses.json): {len(seed_pairs)}")

    THRESHOLD = 0.02
    PAIR_LIVE = 0.02
    FREQ_C_MIN = 0.05

    candidates = []
    seen = set()

    # ── 1.1 estender pares correlatos ────────────────────────────────
    print(f"\n[1.1] Estendendo pares com r<-0.5 (P_triple<{THRESHOLD}, freq_c>{FREQ_C_MIN}) ...")
    n_pair_ext = 0
    for ea, eb, r in seed_pairs:
        # iterar todas as outras arestas como ec
        T_ab = T32[:, ea] * T32[:, eb]   # (N,) — tours onde a e b ativos
        # P_triple = sum(T_ab * T[:, ec]) / N
        triple_counts = T_ab @ T32  # (E,) — contagem para cada ec
        for ec in range(E):
            if ec == ea or ec == eb:
                continue
            if freq[ec] < FREQ_C_MIN:
                continue
            p_triple = triple_counts[ec] / N
            if p_triple >= THRESHOLD:
                continue
            # minimal? subpares ainda vivos
            if (P2[ea, eb] < PAIR_LIVE or
                P2[ea, ec] < PAIR_LIVE or
                P2[eb, ec] < PAIR_LIVE):
                continue
            tup = tuple(sorted((ea, eb, ec)))
            if tup in seen:
                continue
            seen.add(tup)
            candidates.append({
                "edges": list(tup),
                "edge_names": [edge_label(e) for e in tup],
                "P_triple": float(p_triple),
                "P_pairs": [round(float(P2[tup[0], tup[1]]), 4),
                            round(float(P2[tup[0], tup[2]]), 4),
                            round(float(P2[tup[1], tup[2]]), 4)],
                "freq_individual": [round(float(freq[e]), 4) for e in tup],
                "source": "pair_extension",
            })
            n_pair_ext += 1

    print(f"  candidatas via pair-extension: {n_pair_ext}")

    # ── 1.2 vértices de baixo grau ──────────────────────────────────
    print(f"\n[1.2] Triplas em vértices de grau ≤ 4 ...")
    inc = vertex_to_incident_edges()
    from itertools import combinations
    n_lowdeg = 0
    for v, lst in inc.items():
        d = len(lst)
        if d > 4:
            continue
        for trip in combinations(lst, 3):
            ea, eb, ec = trip
            T_abc = T32[:, ea] * T32[:, eb] * T32[:, ec]
            p_triple = float(T_abc.sum() / N)
            if p_triple >= THRESHOLD:
                continue
            if (P2[ea, eb] < PAIR_LIVE or
                P2[ea, ec] < PAIR_LIVE or
                P2[eb, ec] < PAIR_LIVE):
                continue
            tup = tuple(sorted(trip))
            if tup in seen:
                continue
            seen.add(tup)
            candidates.append({
                "edges": list(tup),
                "edge_names": [edge_label(e) for e in tup],
                "P_triple": p_triple,
                "P_pairs": [round(float(P2[tup[0], tup[1]]), 4),
                            round(float(P2[tup[0], tup[2]]), 4),
                            round(float(P2[tup[1], tup[2]]), 4)],
                "freq_individual": [round(float(freq[e]), 4) for e in tup],
                "source": "low_degree_vertex",
                "vertex": label(v),
                "vertex_degree": d,
            })
            n_lowdeg += 1
    print(f"  candidatas via low-degree vertex: {n_lowdeg}")

    # ── relatório ─────────────────────────────────────────────────────
    print(f"\nTotal candidatas (dedup): {len(candidates)}")

    if candidates:
        # histograma de P_triple
        ps = [c["P_triple"] for c in candidates]
        bins = [0, 0.0001, 0.001, 0.005, 0.01, 0.02]
        hist = [0] * (len(bins) - 1)
        for p in ps:
            for i in range(len(bins) - 1):
                if bins[i] <= p < bins[i + 1]:
                    hist[i] += 1
                    break
        print("\nDistribuição P_triple:")
        for i in range(len(bins) - 1):
            print(f"  [{bins[i]}, {bins[i+1]}): {hist[i]}")

        # top 10 (menor P)
        candidates.sort(key=lambda c: c["P_triple"])
        print("\nTop 10 candidatas (menor P_triple):")
        for c in candidates[:10]:
            print(f"  {c['edge_names']}  "
                  f"P_triple={c['P_triple']:.5f}  "
                  f"P_pairs={c['P_pairs']}  "
                  f"src={c['source']}")

    with open(DATA / "triples_candidates.json", "w") as f:
        json.dump({
            "n_tours": int(N),
            "E": int(E),
            "threshold_p_triple": THRESHOLD,
            "threshold_pair_live": PAIR_LIVE,
            "threshold_freq_c": FREQ_C_MIN,
            "n_candidates": len(candidates),
            "candidates": candidates,
        }, f, indent=2)
    print(f"\nSalvo: data/triples_candidates.json")


if __name__ == "__main__":
    main()
