#!/usr/bin/env python3
"""
width6_threshold.py
===================
O limiar do enunciado forte e' em `n` ou em `min(n,m)`?

Dados:
  6x6 -> 12 excecoes ; 6x7 -> 4 (todas ancoradas em canto) ; 8x8 -> 0

Isso e' compativel com "limiar em min(n,m) >= 8" tanto quanto com "n >= 8".
As duas hipoteses so' se separam em tabuleiros com um lado 6 e o outro grande.

Teste: as 4 formas da orbita alongada (as que SOBREVIVEM em 6x7), colocadas
ancoradas em canto em 6x8, 6x10, 6x12.
  - realizadas   -> limiar NAO e' em min(n,m)  (conclusivo, e' certificado)
  - persistem    -> limiar e' em min(n,m)      (precisa da busca exaustiva)

Estagio 1 e' bola de BFS -- que, como se aprendeu em n=8, e' a ferramenta
certa para certificados positivos. Estagio 2 e' busca restrita exaustiva, que
aqui e' barata porque os alvos estao encostados no canto (ha' grau baixo para
podar).
"""
from __future__ import annotations

import argparse
import json
import random
import time
from collections import deque
from pathlib import Path
from typing import Dict, List, Tuple

from flip_graph import Basis, alternating_matchings, cycle_mask, is_hamiltonian
from flip_graph_rect import build_rect, hexagons_bulk_rect, is_corner_rect

DATA = Path(__file__).resolve().parent / "data"
MOVES = [(1, 2), (2, 1), (-1, 2), (-2, 1), (1, -2), (2, -1), (-1, -2), (-2, -1)]


def tour_mask_rect(seq, eidx) -> int:
    x = 0
    for a, b in zip(seq, tuple(seq[1:]) + (seq[0],)):
        x ^= 1 << eidx[(min(a, b), max(a, b))]
    return x


def warnsdorff_rect(R: int, C: int, adj, rng, tries: int = 20000):
    total = R * C
    for _ in range(tries):
        start = rng.randrange(total)
        path = [start]
        used = [False] * total
        used[start] = True
        ok = True
        for _ in range(total - 1):
            cur = path[-1]
            cand = [u for u in adj[cur] if not used[u]]
            if not cand:
                ok = False
                break
            best = min(sum(1 for w in adj[u] if not used[w]) for u in cand)
            cand = [u for u in cand
                    if sum(1 for w in adj[u] if not used[w]) == best]
            nxt = rng.choice(cand)
            path.append(nxt)
            used[nxt] = True
        if ok and start in adj[path[-1]]:
            return path
    return None


def ham_rect(mask: int, V: int, edges) -> bool:
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
        if cur == 0 or cnt > V:
            break
    return cnt == V


def cells_to_cycle(cells, C: int, adj):
    vs = [r * C + c for r, c in cells]
    if len(set(vs)) != 6:
        return None
    s = set(vs)
    start = vs[0]
    order = [start]
    seen = {start}

    def rec():
        cur = order[-1]
        if len(order) == 6:
            return start in adj[cur]
        for u in adj[cur]:
            if u in s and u not in seen:
                order.append(u)
                seen.add(u)
                if rec():
                    return True
                seen.discard(u)
                order.pop()
        return False
    return tuple(order) if rec() else None


def elongated_shapes() -> List[Tuple[Tuple[int, int], ...]]:
    """As 4 formas nao realizadas em 6x7 (a orbita alongada de Klein)."""
    d = json.load(open(DATA / "rect_unrealized_6x7.json"))
    out = set()
    for cells in d["unrealized_cells"]:
        cells = [tuple(x) for x in cells]
        r0 = min(r for r, _ in cells)
        c0 = min(c for _, c in cells)
        out.add(tuple(sorted((r - r0, c - c0) for r, c in cells)))
    return sorted(out)


def run(R: int, C: int, ball: int, seeds: int) -> dict:
    t0 = time.time()
    adj, edges, eidx = build_rect(R, C)
    V = R * C
    hexes = hexagons_bulk_rect(R, C, adj)
    hx = [(cycle_mask(c, eidx), *alternating_matchings(c, eidx)) for c in hexes]

    shapes = elongated_shapes()
    targets: Dict[int, Tuple[int, int, bool]] = {}
    for si, sh in enumerate(shapes):
        H = max(r for r, _ in sh) + 1
        W = max(c for _, c in sh) + 1
        for r0 in range(R - H + 1):
            for c0 in range(C - W + 1):
                cells = [(r + r0, c + c0) for r, c in sh]
                if any(is_corner_rect(R, C, r * C + c) for r, c in cells):
                    continue
                cyc = cells_to_cycle(cells, C, adj)
                if cyc is None:
                    continue
                anchored = (any(r in (0, R - 1) for r, _ in cells)
                            and any(c in (0, C - 1) for _, c in cells))
                targets[cycle_mask(cyc, eidx)] = (si, r0 * 100 + c0, anchored)

    n_anch = sum(1 for v in targets.values() if v[2])
    print(f"[{R}x{C}] hex bulk={len(hexes)}  alvos={len(targets)} "
          f"(ancorados em canto: {n_anch})", flush=True)

    realized = set()
    total_seen = 0
    for s in range(seeds):
        rng = random.Random(1000 + s)
        p = warnsdorff_rect(R, C, adj, rng)
        if p is None:
            continue
        root = tour_mask_rect(p, eidx)
        seen = {root}
        q = deque([root])
        while q and len(seen) < ball:
            m = q.popleft()
            for cm, m0, m1 in hx:
                it = m & cm
                if it != m0 and it != m1:
                    continue
                nm = m ^ cm
                if not ham_rect(nm, V, edges):
                    continue
                realized.add(cm)
                if nm not in seen:
                    seen.add(nm)
                    q.append(nm)
        total_seen += len(seen)
        hit = len(set(targets) & realized)
        print(f"  semente {s}: bola={len(seen)}  alvos certificados="
              f"{hit}/{len(targets)}  hex realizados={len(realized)}/{len(hexes)}"
              f"  ({time.time()-t0:.0f}s)", flush=True)
        if hit == len(targets):
            break

    hit = set(targets) & realized
    miss = {cm: targets[cm] for cm in targets if cm not in realized}
    miss_anch = sum(1 for v in miss.values() if v[2])

    dimZb = (len(edges) - 8) - V + 5
    B = Basis()
    for cm in realized:
        B.add(cm)

    res = {
        "R": R, "C": C, "min_side": min(R, C),
        "n_hexagons_bulk": len(hexes),
        "n_targets": len(targets), "n_targets_corner_anchored": n_anch,
        "targets_certified": len(hit),
        "targets_missing": len(miss),
        "targets_missing_corner_anchored": miss_anch,
        "hexagons_realized": len(realized),
        "rank_realized": B.rank, "dim_ZBulk": dimZb,
        "ball_cap": ball, "seeds_used": seeds, "tours_seen_total": total_seen,
        "missing_cells": [sorted(divmod(v, C) for v in
                                 next(c for c in hexes
                                      if cycle_mask(c, eidx) == cm))
                          for cm in miss],
        "seconds": round(time.time() - t0, 1),
    }
    return res


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--boards", type=str, default="6x8,6x10,6x12")
    ap.add_argument("--ball", type=int, default=400000)
    ap.add_argument("--seeds", type=int, default=3)
    args = ap.parse_args()

    out_all = {}
    for b in args.boards.split(","):
        R, C = (int(x) for x in b.split("x"))
        out_all[b] = run(R, C, args.ball, args.seeds)

    (DATA / "width6_threshold.json").write_text(json.dumps(out_all, indent=2))
    print()
    print("=" * 72)
    print(" A ORBITA ALONGADA (as 4 que sobrevivem em 6x7) EM TABULEIROS 6xM")
    print("=" * 72)
    for b, r in out_all.items():
        print(f" {b:>6}: alvos {r['targets_certified']}/{r['n_targets']} "
              f"certificados   faltando {r['targets_missing']} "
              f"(em canto: {r['targets_missing_corner_anchored']})"
              f"   rank {r['rank_realized']}/{r['dim_ZBulk']}")
    print("\n -> data/width6_threshold.json")


if __name__ == "__main__":
    main()
