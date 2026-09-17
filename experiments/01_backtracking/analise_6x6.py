#!/usr/bin/env python3
"""
analise_6x6.py
==============
Análise estatística dos resultados do catalogador 6×6.

Lê:
  summary.json
  destruction_map.json
  loop_correlations.json

Produz:
  - Sanity check das contagens (710064 = 9862 × 36 × 2)
  - Decomposição das arestas em órbitas D₄
  - Arestas obrigatórias agrupadas por órbita
  - Correlações negativas agrupadas por órbita D₄ de pares
    (é aqui que os "8 pares com r=-0.7714" colapsam em UMA invariante)
  - Heatmaps ASCII de destruição média e máxima por casa
  - Ranking das casas por periculosidade + agregação por grau
  - Distribuição do passo de pico de destruição
  - analysis_6x6.json com tudo serializado
"""

import json, collections, math, os
from itertools import combinations

BOARD = 6
TOTAL = BOARD * BOARD
MOVES = [(2,1),(2,-1),(-2,1),(-2,-1),(1,2),(1,-2),(-1,2),(-1,-2)]


# ── identidade de vértice ─────────────────────────────────────────────

def vid(r, c): return r * BOARD + c
def vrc(v):    return divmod(v, BOARD)

def label(v):
    r, c = vrc(v)
    return chr(65 + c) + str(BOARD - r)


# ── grupo D₄ sobre o tabuleiro ────────────────────────────────────────
# Convenção (com N = BOARD - 1):
#   t=0 identidade           (r, c)
#   t=1 rotação 90° (CW)     (c, N-r)
#   t=2 rotação 180°         (N-r, N-c)
#   t=3 rotação 270° (CW)    (N-c, r)
#   t=4 reflexão horizontal  (r, N-c)
#   t=5 reflexão vertical    (N-r, c)
#   t=6 transposição         (c, r)
#   t=7 anti-transposição    (N-c, N-r)

def d4_apply(rc, t):
    r, c = rc
    N = BOARD - 1
    if   t == 0: return (r,     c    )
    elif t == 1: return (c,     N - r)
    elif t == 2: return (N - r, N - c)
    elif t == 3: return (N - c, r    )
    elif t == 4: return (r,     N - c)
    elif t == 5: return (N - r, c    )
    elif t == 6: return (c,     r    )
    elif t == 7: return (N - c, N - r)
    raise ValueError(t)

def d4_vertex(v, t):
    return vid(*d4_apply(vrc(v), t))

def d4_edge(edge, t):
    a, b = edge
    a2, b2 = d4_vertex(a, t), d4_vertex(b, t)
    return (min(a2, b2), max(a2, b2))


# ── grafo do cavalo ───────────────────────────────────────────────────

def build_edges():
    EDGES = set()
    for r in range(BOARD):
        for c in range(BOARD):
            v = vid(r, c)
            for dr, dc in MOVES:
                nr, nc = r + dr, c + dc
                if 0 <= nr < BOARD and 0 <= nc < BOARD:
                    u = vid(nr, nc)
                    e = (min(u, v), max(u, v))
                    EDGES.add(e)
    return sorted(EDGES)


# ── órbitas D₄ ────────────────────────────────────────────────────────

def edge_orbits(edges):
    edge_set = set(edges)
    seen = set()
    orbits = []
    for e in edges:
        if e in seen:
            continue
        orb = {d4_edge(e, t) for t in range(8)}
        orbits.append(sorted(orb & edge_set))
        seen |= orb
    return orbits

def edge_to_orbit_id(orbits):
    m = {}
    for i, orb in enumerate(orbits):
        for e in orb:
            m[e] = i
    return m

def pair_d4_canonical(pair):
    """Representante canônico (lex-min) da órbita D₄ de um par de arestas."""
    a, b = pair
    images = set()
    for t in range(8):
        ea = d4_edge(a, t)
        eb = d4_edge(b, t)
        images.add(tuple(sorted([ea, eb])))
    return min(images)


# ── carregamento ──────────────────────────────────────────────────────

def load_json(path):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


# ── relatório ─────────────────────────────────────────────────────────

def section(title):
    print()
    print("── " + title + " " + "─" * max(0, 64 - len(title)))

def main():
    here = os.path.dirname(os.path.abspath(__file__))
    summary = load_json(os.path.join(here, "summary.json"))
    dmap    = load_json(os.path.join(here, "destruction_map.json"))
    corrs   = load_json(os.path.join(here, "loop_correlations.json"))

    EDGES = build_edges()
    label_to_vid = {label(v): v for v in range(TOTAL)}
    label_to_edge = {}
    for e in EDGES:
        label_to_edge[f"{label(e[0])}-{label(e[1])}"] = e

    print("=" * 70)
    print("ANÁLISE ESTATÍSTICA — CAVALO 6×6 FECHADO")
    print("=" * 70)

    # ── 1. SANITY CHECK ───────────────────────────────────────────────
    section("1. Sanity check das contagens")
    total = summary["total_solutions"]
    print(f"  Sequências totais:                    {total:>10,}")
    print(f"  H₁ inicial (dim. espaço de loops):    {summary['H1_initial']:>10}")
    print(f"  Arestas do grafo:                     {summary['total_edges']:>10}")
    print()
    print("  Decomposição esperada das sequências:")
    print(f"    total / 36 (rotação)         = {total // 36:>8,}  ciclos direcionados (anchored)")
    print(f"    total / (36×2)               = {total // 72:>8,}  ciclos não-direcionados  (esperado 9.862)")
    print(f"    total / (36×2×8)             = {total // 576:>8,}  classes D₄ aproximadas   (paper diz 1.232)")
    print()
    counts = list(summary["solutions_by_start"].values())
    if all(c == counts[0] for c in counts):
        print(f"  Soluções por casa de partida: {counts[0]:,} (uniforme — D₄ confirmada)")

    # ── 2. ÓRBITAS D₄ DAS ARESTAS ─────────────────────────────────────
    section("2. Decomposição das arestas em órbitas D₄")
    orbits = edge_orbits(EDGES)
    e2o = edge_to_orbit_id(orbits)
    sizes = collections.Counter(len(o) for o in orbits)
    print(f"  Total de arestas:  {len(EDGES)}")
    print(f"  Órbitas D₄:        {len(orbits)}")
    for sz in sorted(sizes):
        print(f"    {sizes[sz]:>2} órbita(s) de tamanho {sz}  →  {sz * sizes[sz]} arestas")
    print()
    print("  Representante de cada órbita (lex-min):")
    for i, orb in enumerate(orbits):
        rep = orb[0]
        print(f"    O{i:02d} (|{len(orb)}|): {label(rep[0])}-{label(rep[1])}")

    # ── 3. ARESTAS OBRIGATÓRIAS POR ÓRBITA ────────────────────────────
    section("3. Arestas obrigatórias (freq=100%) por órbita D₄")
    mand_edges = []
    for m in summary["mandatory_edges"]:
        u = label_to_vid[m["edge"][0]]
        v = label_to_vid[m["edge"][1]]
        mand_edges.append((min(u, v), max(u, v)))

    by_orbit = collections.defaultdict(list)
    for e in mand_edges:
        by_orbit[e2o[e]].append(e)

    print(f"  {len(mand_edges)} arestas obrigatórias  →  {len(by_orbit)} órbita(s)")
    for oid, eds in sorted(by_orbit.items()):
        orb = orbits[oid]
        print(f"    Órbita O{oid:02d} (tamanho {len(orb)}):  {len(eds)}/{len(orb)} arestas obrigatórias")
        for e in eds:
            print(f"        {label(e[0])}-{label(e[1])}")
    if all(len(eds) == len(orbits[oid]) for oid, eds in by_orbit.items()):
        print()
        print("  ✓ Todas as órbitas obrigatórias estão completas (consistência D₄).")

    # ── 4. CORRELAÇÕES NEGATIVAS POR ÓRBITA D₄ DE PARES ──────────────
    section("4. Correlações negativas — colapso por órbita D₄")
    pair_buckets = collections.defaultdict(list)
    for c in corrs["top_negative_correlations"]:
        ea = label_to_edge[c["edge_i"]]
        eb = label_to_edge[c["edge_j"]]
        pair = tuple(sorted([ea, eb]))
        canon = pair_d4_canonical(pair)
        pair_buckets[canon].append((c["r"], pair))

    bucket_summary = []
    for canon, members in pair_buckets.items():
        bucket_summary.append((min(r for r, _ in members), len(members), canon, members))
    bucket_summary.sort()

    n_pairs = len(corrs["top_negative_correlations"])
    n_orbits = len(bucket_summary)
    print(f"  {n_pairs} pares listados  →  {n_orbits} órbitas D₄ distintas")
    print(f"  (compactação: {n_pairs/n_orbits:.1f}× — confirma simetria do grafo)")
    print()
    print(f"  {'r':>9}  {'|orb|':>5}  representante                          tipo")
    print(f"  {'-'*9}  {'-'*5}  {'-'*38}  {'-'*16}")
    for r, n, canon, members in bucket_summary[:15]:
        ea, eb = canon
        shared = set(ea) & set(eb)
        if shared:
            tipo = f"compart. v={label(next(iter(shared)))}"
        else:
            tipo = "arestas disjuntas"
        rep = f"{label(ea[0])}-{label(ea[1])} ↔ {label(eb[0])}-{label(eb[1])}"
        print(f"  {r:>9.4f}  {n:>5}  {rep:<38}  {tipo}")
    if n_orbits > 15:
        print(f"  ... ({n_orbits - 15} órbitas adicionais omitidas)")

    # ── 5. HEATMAP DESTRUIÇÃO MÉDIA ───────────────────────────────────
    section("5. Heatmap de destruição média por casa")
    print("    (média de loops H₁ destruídos quando a casa é visitada)")
    print()
    header = "         " + "  ".join(chr(65 + c).center(6) for c in range(BOARD))
    print(header)
    for r in range(BOARD):
        rank = BOARD - r
        cells = []
        for c in range(BOARD):
            md = dmap[label(vid(r, c))]["mean_destruction"]
            cells.append(f"{md:.3f}".rjust(6))
        print(f"  rank {rank}  " + "  ".join(cells))
    print()
    print("    (degradê: corner→2 vizinhos, borda→3-4, miolo→6-8)")

    # ── 6. HEATMAP MAX DESTRUCTION ────────────────────────────────────
    section("6. Heatmap do max-destruction por casa")
    print()
    print(header)
    for r in range(BOARD):
        rank = BOARD - r
        cells = []
        for c in range(BOARD):
            mx = dmap[label(vid(r, c))]["max_destruction"]
            cells.append(f"{mx:>4d}  ")
        print(f"  rank {rank}  " + "  ".join(cells))

    # ── 7. RANKING POR PERICULOSIDADE ─────────────────────────────────
    section("7. Ranking de casas por destruição média")
    rows = []
    for lbl, data in dmap.items():
        rows.append((data["mean_destruction"], data["max_destruction"], data["degree"], lbl))
    rows.sort(reverse=True)

    print(f"  {'#':>3}  {'casa':>4}  {'grau':>4}  {'média':>7}  {'max':>4}  distribuição (kills:freq)")
    print(f"  {'-'*3}  {'-'*4}  {'-'*4}  {'-'*7}  {'-'*4}  {'-'*40}")
    for i, (mean, mx, deg, lbl) in enumerate(rows[:12]):
        dist = dmap[lbl]["destruction_distribution"]
        dist_str = " ".join(f"{k}:{int(v):,}" for k, v in sorted(dist.items(), key=lambda x: int(x[0])))
        print(f"  {i+1:>3}  {lbl:>4}  {deg:>4}  {mean:>7.4f}  {mx:>4}  {dist_str}")

    print()
    print("  ── menos perigosas (8 últimas) ──")
    for i, (mean, mx, deg, lbl) in enumerate(rows[-8:]):
        rk = len(rows) - 8 + i + 1
        print(f"  {rk:>3}  {lbl:>4}  {deg:>4}  {mean:>7.4f}  {mx:>4}")

    # ── 8. AGREGADO POR GRAU ──────────────────────────────────────────
    section("8. Destruição agregada por grau do vértice")
    by_deg = collections.defaultdict(list)
    for lbl, data in dmap.items():
        by_deg[data["degree"]].append(data["mean_destruction"])

    print(f"  {'grau':>4}  {'casas':>5}  {'média':>8}  {'desvio':>7}  {'mín':>7}  {'máx':>7}")
    print(f"  {'-'*4}  {'-'*5}  {'-'*8}  {'-'*7}  {'-'*7}  {'-'*7}")
    for deg in sorted(by_deg.keys()):
        v = by_deg[deg]
        m = sum(v) / len(v)
        sd = (sum((x - m)**2 for x in v) / len(v)) ** 0.5
        print(f"  {deg:>4}  {len(v):>5}  {m:>8.4f}  {sd:>7.4f}  {min(v):>7.4f}  {max(v):>7.4f}")

    # ── 9. PASSO DE PICO ──────────────────────────────────────────────
    section("9. Distribuição do passo de pico de destruição")
    peak = summary["peak_step_distribution"]
    peak_pairs = sorted(peak.items(), key=lambda kv: int(kv[0]))
    n = sum(int(v) for v in peak.values())
    print(f"  {'passo':>5}  {'N':>9}  {'%':>6}  {'cum%':>6}  histograma")
    cum = 0
    for k, v in peak_pairs:
        cum += v
        pct = 100 * v / n
        cum_pct = 100 * cum / n
        bar = "█" * max(1, int(round(pct)))
        print(f"  {int(k):>5d}  {v:>9,}  {pct:>5.1f}%  {cum_pct:>5.1f}%  {bar}")

    mean_peak = sum(int(k) * v for k, v in peak.items()) / n
    median_step = None
    cum = 0
    for k, v in peak_pairs:
        cum += v
        if cum >= n / 2:
            median_step = int(k)
            break
    print()
    print(f"  Passo de pico médio:    {mean_peak:.2f}")
    print(f"  Passo de pico mediano:  {median_step}")
    print(f"  Interpretação: a destruição de loops é precoce —")
    print(f"  metade dos casos atinge o pico antes do passo {median_step}.")

    # ── 10. SAÍDA JSON ────────────────────────────────────────────────
    analysis = {
        "counts": {
            "sequences_total": total,
            "directed_cycles_anchored": total // 36,
            "undirected_cycles": total // 72,
            "d4_classes_estimated": total // 576,
            "h1_initial": summary["H1_initial"],
            "edges": summary["total_edges"],
        },
        "edge_orbits": [
            {
                "id": i,
                "size": len(o),
                "representative": f"{label(o[0][0])}-{label(o[0][1])}",
                "members": [f"{label(e[0])}-{label(e[1])}" for e in o],
            }
            for i, o in enumerate(orbits)
        ],
        "mandatory_edges_by_orbit": [
            {
                "orbit_id": oid,
                "orbit_size": len(orbits[oid]),
                "edges": [f"{label(e[0])}-{label(e[1])}" for e in eds],
            }
            for oid, eds in sorted(by_orbit.items())
        ],
        "correlation_orbits": [
            {
                "min_r": r,
                "size": n_,
                "representative": [
                    f"{label(canon[0][0])}-{label(canon[0][1])}",
                    f"{label(canon[1][0])}-{label(canon[1][1])}",
                ],
                "shared_vertex": (
                    label(next(iter(set(canon[0]) & set(canon[1]))))
                    if (set(canon[0]) & set(canon[1])) else None
                ),
                "members": [
                    {
                        "r": rr,
                        "edges": [
                            f"{label(p[0][0])}-{label(p[0][1])}",
                            f"{label(p[1][0])}-{label(p[1][1])}",
                        ],
                    }
                    for rr, p in mems
                ],
            }
            for r, n_, canon, mems in bucket_summary
        ],
        "destruction_ranking": [
            {
                "rank": i + 1, "label": lbl, "degree": deg,
                "mean": mean, "max": mx,
                "distribution": dmap[lbl]["destruction_distribution"],
            }
            for i, (mean, mx, deg, lbl) in enumerate(rows)
        ],
        "destruction_by_degree": {
            str(deg): {
                "n_squares": len(v),
                "mean": round(sum(v) / len(v), 6),
                "stddev": round((sum((x - sum(v)/len(v))**2 for x in v) / len(v))**0.5, 6),
                "min": round(min(v), 6),
                "max": round(max(v), 6),
            }
            for deg, v in by_deg.items()
        },
        "peak_step": {
            "distribution": {str(int(k)): int(v) for k, v in peak_pairs},
            "mean": round(mean_peak, 4),
            "median": median_step,
        },
    }

    out = os.path.join(here, "analysis_6x6.json")
    with open(out, "w", encoding="utf-8") as f:
        json.dump(analysis, f, indent=2, ensure_ascii=False)

    print()
    print("=" * 70)
    print(f"Análise salva em: {os.path.basename(out)}")
    print("=" * 70)


if __name__ == "__main__":
    main()
