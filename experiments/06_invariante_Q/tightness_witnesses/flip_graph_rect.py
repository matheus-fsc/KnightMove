#!/usr/bin/env python3
"""
flip_graph_rect.py
==================
Segundo ponto EXAUSTIVO da conectividade do grafo de flips, em geometria
diferente: tabuleiros retangulares R x C.

Q(n,m)=3 ja' esta' provado para retangulos (n,m >= 6), entao a reducao
`Z_bulk ⊆ Span(Ham) => rank = beta1 - 3` se aplica igual.

Enumera TODOS os tours fechados por backtracking com poda de conectividade
(cada ciclo uma unica vez: comeca no vertice 0 e fixa a orientacao), monta o
grafo de flips com o criterio O(1) de `flip_graph.py` e conta componentes.

CLI:
  python flip_graph_rect.py --rows 6 --cols 7
  python flip_graph_rect.py --rows 6 --cols 8
"""
from __future__ import annotations

import argparse
import json
import time
from collections import Counter
from pathlib import Path
from typing import Dict, List, Tuple

from flip_graph import (
    Basis, DSU, alternating_matchings, cycle_mask, is_hamiltonian,
)

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
MOVES = [(1, 2), (2, 1), (-1, 2), (-2, 1), (1, -2), (2, -1), (-1, -2), (-2, -1)]


def build_rect(R: int, C: int):
    adj: List[List[int]] = [[] for _ in range(R * C)]
    edges: List[Tuple[int, int]] = []
    for r in range(R):
        for c in range(C):
            i = r * C + c
            for dr, dc in MOVES:
                rr, cc = r + dr, c + dc
                if 0 <= rr < R and 0 <= cc < C:
                    j = rr * C + cc
                    adj[i].append(j)
                    if i < j:
                        edges.append((i, j))
    for a in adj:
        a.sort()
    edges.sort()
    return adj, edges, {e: k for k, e in enumerate(edges)}


def is_corner_rect(R: int, C: int, v: int) -> bool:
    r, c = divmod(v, C)
    return r in (0, R - 1) and c in (0, C - 1)


def hexagons_bulk_rect(R: int, C: int, adj) -> List[Tuple[int, ...]]:
    total = R * C
    ok = [not is_corner_rect(R, C, v) for v in range(total)]
    out: List[Tuple[int, ...]] = []
    for v1 in range(total):
        if not ok[v1]:
            continue
        path = [v1]
        seen = {v1}

        def rec():
            cur = path[-1]
            if len(path) == 6:
                if v1 in adj[cur] and path[1] < path[-1]:
                    out.append(tuple(path))
                return
            for u in adj[cur]:
                if u <= v1 or u in seen or not ok[u]:
                    continue
                path.append(u)
                seen.add(u)
                rec()
                seen.discard(u)
                path.pop()
        rec()
    return out


def enumerate_tours(R: int, C: int, adj, eidx) -> List[int]:
    """Todos os tours fechados, como mascaras de arestas, cada um uma vez."""
    total = R * C
    used = [False] * total
    path = [0]
    used[0] = True
    out: List[int] = []

    def viable() -> bool:
        cur = path[-1]
        free = [v for v in range(total) if not used[v]]
        if not free:
            return True
        for v in free:
            k = 0
            for u in adj[v]:
                if (not used[u]) or u == cur or u == 0:
                    k += 1
                    if k >= 2:
                        break
            if k < 2:
                return False
        seen = {cur}
        stack = [cur]
        cnt = 0
        while stack:
            x = stack.pop()
            for u in adj[x]:
                if not used[u] and u not in seen:
                    seen.add(u)
                    stack.append(u)
                    cnt += 1
        return cnt == len(free)

    def rec():
        cur = path[-1]
        if len(path) == total:
            if 0 in adj[cur]:
                m = 0
                for a, b in zip(path, path[1:] + path[:1]):
                    m ^= 1 << eidx[(min(a, b), max(a, b))]
                out.append(m)
            return
        for u in adj[cur]:
            if used[u]:
                continue
            # fixa a orientacao: o segundo vertice do ciclo e' o menor vizinho
            if len(path) == 1 and u > min(adj[0]):
                continue
            used[u] = True
            path.append(u)
            if viable():
                rec()
            path.pop()
            used[u] = False

    rec()
    return out


def run(R: int, C: int) -> dict:
    t0 = time.time()
    adj, edges, eidx = build_rect(R, C)
    V, E = R * C, len(edges)
    beta1 = E - V + 1
    dimZb = (E - 8) - V + 5
    print(f"[{R}x{C}] V={V} E={E} beta1={beta1} dim_ZBulk={dimZb}")

    hexes = hexagons_bulk_rect(R, C, adj)
    hx = [(cycle_mask(c, eidx), *alternating_matchings(c, eidx)) for c in hexes]
    print(f"[{R}x{C}] hexagonos do bulk: {len(hexes)}")

    masks = enumerate_tours(R, C, adj, eidx)
    t_enum = time.time() - t0
    print(f"[{R}x{C}] tours fechados (exaustivo): {len(masks)}  "
          f"({t_enum:.0f}s)")
    if not masks:
        return {"R": R, "C": C, "n_tours": 0}

    idx = {m: i for i, m in enumerate(masks)}
    assert len(idx) == len(masks), "tours repetidos na enumeracao"

    dsu = DSU(len(masks))
    deg = [0] * len(masks)
    realized: Dict[int, int] = {}
    inside = disconn = external = 0

    class _B:
        pass

    def ham(mask: int) -> bool:
        nb: List[List[int]] = [[] for _ in range(V)]
        m = mask
        while m:
            b = m & -m
            u, v = edges[b.bit_length() - 1]
            nb[u].append(v)
            nb[v].append(u)
            m ^= b
        prev, cur, cnt = -1, 0, 0
        while True:
            a, b2 = nb[cur]
            nxt = a if a != prev else b2
            prev, cur = cur, nxt
            cnt += 1
            if cur == 0:
                break
            if cnt > V:
                return False
        return cnt == V

    for i, tm in enumerate(masks):
        for cm, m0, m1 in hx:
            it = tm & cm
            if it != m0 and it != m1:
                continue
            inside += 1
            nm = tm ^ cm
            if not ham(nm):
                disconn += 1
                continue
            deg[i] += 1
            realized[cm] = realized.get(cm, 0) + 1
            j = idx.get(nm)
            if j is None:
                external += 1
                continue
            dsu.union(i, j)

    B = Basis()
    for cm in realized:
        B.add(cm)
    Ball = Basis()
    for cm, _, _ in hx:
        Ball.add(cm)

    sizes = Counter(dsu.find(i) for i in range(len(masks)))
    comp = sorted(sizes.values(), reverse=True)

    res = {
        "R": R, "C": C, "V": V, "E": E, "beta1": beta1, "dim_ZBulk": dimZb,
        "n_hexagons_bulk": len(hexes),
        "n_tours": len(masks),
        "enum_seconds": round(t_enum, 1),
        "pairs_matching_ok": inside,
        "pairs_disconnected": disconn,
        "pairs_external": external,
        "flip_degree_min": min(deg), "flip_degree_max": max(deg),
        "flip_degree_mean": round(sum(deg) / len(deg), 3),
        "flip_degree_hist": dict(sorted(Counter(deg).items())),
        "n_isolated": sum(1 for d in deg if d == 0),
        "n_components": len(comp),
        "component_sizes": comp[:20],
        "n_hexagons_realized": len(realized),
        "rank_realized": B.rank,
        "rank_all_bulk_hexagons": Ball.rank,
        "realized_span_ZBulk": B.rank == dimZb,
        "all_hex_span_ZBulk": Ball.rank == dimZb,
        "seconds": round(time.time() - t0, 1),
    }
    return res


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rows", type=int, default=6)
    ap.add_argument("--cols", type=int, default=7)
    args = ap.parse_args()
    res = run(args.rows, args.cols)
    out = DATA / f"flip_graph_rect_{args.rows}x{args.cols}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(res, indent=2))
    print()
    print("=" * 66)
    print(f" GRAFO DE FLIPS  {args.rows}x{args.cols}  (EXAUSTIVO)")
    print("=" * 66)
    for k in ("n_tours", "n_hexagons_bulk", "pairs_matching_ok",
              "pairs_disconnected", "pairs_external", "flip_degree_min",
              "flip_degree_mean", "flip_degree_max", "n_isolated",
              "n_components", "component_sizes", "dim_ZBulk", "rank_realized",
              "realized_span_ZBulk", "rank_all_bulk_hexagons",
              "all_hex_span_ZBulk", "n_hexagons_realized", "seconds"):
        print(f" {k:26s}: {res.get(k)}")
    print(f" -> {out}")


if __name__ == "__main__":
    main()
