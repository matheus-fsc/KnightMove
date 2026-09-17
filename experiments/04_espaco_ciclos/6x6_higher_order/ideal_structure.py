#!/usr/bin/env python3
"""
ideal_structure.py
==================
Análise estrutural do ideal de infactibilidade do tour 6×6:
hierarquia, cobertura de vértices, geração e conexão com H₁.

Entradas:
  data/exclusions_pairs.json
  data/exclusions_triples.json
  data/exclusions_quads.json
  data/incidence_matrix_6x6.npy
  data/edges_6x6.json
  data/freq_singles.npy

Saídas:
  data/ideal_summary.json
  data/plots/vertex_coverage.png
  data/plots/ideal_structure.png
"""

import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
PLOTS = DATA / "plots"

sys.path.insert(0, str(ROOT.parent))
from cavalo_loop_destruicao_6x6 import (  # noqa: E402
    EDGES_LIST, ADJ, TOTAL, BOARD, vid, vrc, label
)


def edge_vertices(e_idx):
    u, v = EDGES_LIST[e_idx]
    return u, v


def vertex_class(v):
    r, c = vrc(v)
    on_h = r in (0, BOARD - 1)
    on_v = c in (0, BOARD - 1)
    if on_h and on_v:
        return "corner"
    if on_h or on_v:
        return "edge"
    if r in (1, BOARD - 2) or c in (1, BOARD - 2):
        return "near_edge"
    return "interior"


def load_minimals():
    pairs = json.loads((DATA / "exclusions_pairs.json").read_text())["exclusions"]
    triples = json.loads((DATA / "exclusions_triples.json").read_text())["minimal_triples"]
    quads = json.loads((DATA / "exclusions_quads.json").read_text())["minimal_quads"]
    return pairs, triples, quads


def main():
    PLOTS.mkdir(parents=True, exist_ok=True)

    pairs, triples, quads = load_minimals()
    edges = list(EDGES_LIST)
    E = len(edges)
    T = np.load(DATA / "incidence_matrix_6x6.npy")

    n_p, n_t, n_q = len(pairs), len(triples), len(quads)
    print(f"Minimais: pairs={n_p}  triples={n_t}  quads={n_q}")
    print(f"Total geradores do ideal de infactibilidade: {n_p + n_t + n_q}")

    # ── 4.1 Hierarquia (verificar que minimais são realmente independentes) ──
    print("\n[4.1] Verificando independência (nenhum minimal-superior contém minimal-inferior)...")
    pair_set = {tuple(p["edges"]) for p in pairs}
    triple_set = {tuple(t["edges"]) for t in triples}

    def has_subpair(edges_tup):
        from itertools import combinations
        for a, b in combinations(edges_tup, 2):
            if tuple(sorted((a, b))) in pair_set:
                return True
        return False

    def has_subtriple(edges_tup):
        from itertools import combinations
        for combo in combinations(edges_tup, 3):
            if tuple(sorted(combo)) in triple_set:
                return True
        return False

    bad_triples = sum(1 for t in triples if has_subpair(tuple(t["edges"])))
    bad_quads_p = sum(1 for q in quads if has_subpair(tuple(q["edges"])))
    bad_quads_t = sum(1 for q in quads if has_subtriple(tuple(q["edges"])))
    print(f"  triplas com sub-par excluído  : {bad_triples}  (esperado 0)")
    print(f"  quadras com sub-par excluído : {bad_quads_p}  (esperado 0)")
    print(f"  quadras com sub-tripla excl. : {bad_quads_t}  (esperado 0)")

    # ── 4.2 Cobertura de vértices ─────────────────────────────────────
    print("\n[4.2] Cobertura de vértices (quantos minimais envolvem cada vértice)...")
    vertex_count = Counter()
    vertex_count_by_order = {"pair": Counter(), "triple": Counter(),
                             "quad": Counter()}
    for p in pairs:
        verts = set()
        for ei in p["edges"]:
            u, v = edge_vertices(ei)
            verts.update([u, v])
        for v in verts:
            vertex_count[v] += 1
            vertex_count_by_order["pair"][v] += 1
    for t in triples:
        verts = set()
        for ei in t["edges"]:
            u, v = edge_vertices(ei)
            verts.update([u, v])
        for v in verts:
            vertex_count[v] += 1
            vertex_count_by_order["triple"][v] += 1
    for q in quads:
        verts = set()
        for ei in q["edges"]:
            u, v = edge_vertices(ei)
            verts.update([u, v])
        for v in verts:
            vertex_count[v] += 1
            vertex_count_by_order["quad"][v] += 1

    # agregar por classe posicional
    class_counts = defaultdict(list)
    for v in range(TOTAL):
        cls = vertex_class(v)
        class_counts[cls].append(vertex_count[v])
    print("\n  Por classe posicional (média de minimais por vértice):")
    for cls, counts in class_counts.items():
        avg = sum(counts) / len(counts) if counts else 0
        mn = min(counts) if counts else 0
        mx = max(counts) if counts else 0
        print(f"    {cls:>10s}  n={len(counts):>2d}  "
              f"min={mn:>4d}  avg={avg:>7.1f}  max={mx:>4d}")

    # top vértices
    top_vertices = sorted(vertex_count.items(), key=lambda x: -x[1])[:10]
    print("\n  Top 10 vértices mais cobertos:")
    for v, n in top_vertices:
        print(f"    {label(v):>3s} (deg={len(ADJ[v])}, {vertex_class(v)}): {n} minimais "
              f"[p={vertex_count_by_order['pair'][v]}, "
              f"t={vertex_count_by_order['triple'][v]}, "
              f"q={vertex_count_by_order['quad'][v]}]")

    # ── 4.3 Geração do ideal ──────────────────────────────────────────
    print("\n[4.3] Geração do ideal:")
    print(f"  ordem 2 (pares) é insuficiente? "
          f"{'SIM' if n_t > 0 else 'não — só pares geram tudo'}")
    print(f"  ordem 3 (triplas) é insuficiente? "
          f"{'SIM' if n_q > 0 else 'não — pares + triplas geram tudo'}")
    print(f"  => ideal NÃO é finitamente gerado por uma única ordem k ≤ 4")

    # ── 4.4 Conexão com H₁ ─────────────────────────────────────────────
    print("\n[4.4] Assinatura H₁ das exclusões minimais...")
    # base GF(2) via row-reduction local
    def gf2_row_reduce(M):
        A = M.copy().astype(np.uint8)
        rows, cols = A.shape
        r = 0
        pivots = []
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
            pivots.append(c)
            r += 1
        return A[:r], pivots

    # rank de T
    print("  computando base GF(2) ...")
    B, pivots = gf2_row_reduce(T)
    rank = B.shape[0]
    print(f"  rank(T) 6×6 = {rank} (esperado 45 = β₁)")

    # Coordenadas dos minimais: cada exclusão tem uma assinatura "se ativada"
    # — mas como xor dos vetores indicadores das arestas envolvidas é
    # apenas o vetor indicador (porque GF(2)). Mais interessante: para cada
    # par/tripla/quadra, projeta o vetor indicador 1_S no espaço dos ciclos
    # e mede quantas coordenadas H₁ ele atinge.

    def signature_h1(edge_set):
        """Vetor indicador 1_S projetado em B → coordenadas H₁."""
        s = np.zeros(E, dtype=np.uint8)
        for e in edge_set:
            s[e] = 1
        coords = (B @ s.astype(np.uint8)) % 2
        return coords

    def h1_weight(edge_set):
        return int(signature_h1(edge_set).sum())

    weights_p = [h1_weight(p["edges"]) for p in pairs]
    weights_t = [h1_weight(t["edges"]) for t in triples]
    weights_q = [h1_weight(q["edges"]) for q in quads]
    print(f"  peso H₁ (#coords não-zero) das assinaturas indicadoras:")
    print(f"    pares    : avg={np.mean(weights_p):.2f}  min={min(weights_p)}  max={max(weights_p)}")
    print(f"    triplas  : avg={np.mean(weights_t):.2f}  min={min(weights_t)}  max={max(weights_t)}")
    print(f"    quadras  : avg={np.mean(weights_q):.2f}  min={min(weights_q)}  max={max(weights_q)}")

    # ── 4.5 Plots ──────────────────────────────────────────────────────
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        # heatmap de cobertura por casa
        cov = np.zeros((BOARD, BOARD), dtype=np.int32)
        for v in range(TOTAL):
            r, c = vrc(v)
            cov[r, c] = vertex_count[v]

        fig, ax = plt.subplots(figsize=(7, 6))
        im = ax.imshow(cov, cmap="Reds", aspect="equal")
        for r in range(BOARD):
            for c in range(BOARD):
                ax.text(c, r, str(cov[r, c]),
                        ha="center", va="center",
                        color="white" if cov[r, c] > cov.max() / 2 else "black",
                        fontsize=11)
        ax.set_xticks(range(BOARD))
        ax.set_xticklabels(list("ABCDEF"))
        ax.set_yticks(range(BOARD))
        ax.set_yticklabels(list(range(BOARD, 0, -1)))
        ax.set_title(f"Cobertura por casa\n(N minimais envolvendo cada vértice)\n"
                     f"pares={n_p}  triplas={n_t}  quadras={n_q}")
        fig.colorbar(im, ax=ax)
        fig.tight_layout()
        plt.savefig(PLOTS / "vertex_coverage.png", dpi=120)
        plt.close(fig)
        print(f"\nSalvo: data/plots/vertex_coverage.png")

        # plot de peso H₁
        fig, axes = plt.subplots(1, 3, figsize=(13, 4))
        for ax, (data, name) in zip(axes, [
                (weights_p, f"Pares ({n_p})"),
                (weights_t, f"Triplas ({n_t})"),
                (weights_q, f"Quadras ({n_q})")]):
            ax.hist(data, bins=range(0, max(data) + 2), align='left',
                    edgecolor="black")
            ax.set_xlabel("peso H₁ (# coordenadas)")
            ax.set_ylabel("contagem")
            ax.set_title(name)
        fig.suptitle("Assinatura H₁ das exclusões minimais (β₁ = 45)")
        fig.tight_layout()
        plt.savefig(PLOTS / "ideal_structure.png", dpi=120)
        plt.close(fig)
        print(f"Salvo: data/plots/ideal_structure.png")
    except Exception as ex:
        print(f"AVISO: erro ao plotar: {ex}")

    # ── salvar sumário ─────────────────────────────────────────────────
    summary = {
        "n_pairs_minimal": n_p,
        "n_triples_minimal": n_t,
        "n_quads_minimal": n_q,
        "total_ideal_generators_to_order_4": n_p + n_t + n_q,
        "consistency_checks": {
            "triples_with_subpair_excluded": bad_triples,
            "quads_with_subpair_excluded": bad_quads_p,
            "quads_with_subtriple_excluded": bad_quads_t,
        },
        "vertex_coverage_by_class": {
            cls: {
                "n_vertices": len(counts),
                "min": min(counts), "max": max(counts),
                "avg": round(sum(counts) / len(counts), 2),
            } for cls, counts in class_counts.items()
        },
        "top_vertices": [
            {"label": label(v), "v": int(v),
             "degree": len(ADJ[v]),
             "class": vertex_class(v),
             "total": int(n),
             "pairs": int(vertex_count_by_order["pair"][v]),
             "triples": int(vertex_count_by_order["triple"][v]),
             "quads": int(vertex_count_by_order["quad"][v])}
            for v, n in top_vertices
        ],
        "h1_weight": {
            "rank_sampled": int(rank),
            "beta1_theoretical": int(E - 36 + 1),
            "pairs": {"min": int(min(weights_p)), "max": int(max(weights_p)),
                      "avg": round(float(np.mean(weights_p)), 3)},
            "triples": {"min": int(min(weights_t)), "max": int(max(weights_t)),
                        "avg": round(float(np.mean(weights_t)), 3)},
            "quads": {"min": int(min(weights_q)), "max": int(max(weights_q)),
                      "avg": round(float(np.mean(weights_q)), 3)},
        },
    }
    with open(DATA / "ideal_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\nSalvo: data/ideal_summary.json")


if __name__ == "__main__":
    main()
