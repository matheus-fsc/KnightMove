#!/usr/bin/env python3
"""
flip_graph.py
=============
Conectividade do grafo de flips hexagonais sobre os tours fechados.

Aresta do grafo de flips:  tau -- tau XOR C,  C hexagono do bulk.

Criterio barato (sem construir gadget W, sem busca hamiltoniana):

    flip(tau, C) valido  <=>  tau ∩ C e' um dos DOIS emparelhamentos
                              alternados de C  (grau 2 preservado)
                         AND  tau XOR C e' conexo

O primeiro teste e' O(1) com mascaras; o segundo e' um passeio de |V| passos.

FASE 1 (n=6): exaustivo sobre os 9862 tours -- NAO amostrar.

Medidas reportadas:
  1. numero de componentes conexas do grafo de flips (union-find)
  2. distribuicao do grau de flip (grau 0 => resposta e' NAO, imediatamente)
  3. os hexagonos REALIZADOS (os que aparecem em algum flip) geram Z_bulk?

O item 3 e' exatamente (U1)+(U2) em n=6: um flip valido com C = H_A XOR H_B
E' um certificado de C.

CLI:
  python flip_graph.py            # n=6 exaustivo
  python flip_graph.py --n 8 --tours <arquivo>   # fase 2, amostral
"""
from __future__ import annotations

import argparse
import json
import time
from collections import Counter
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"

MOVES = [(1, 2), (2, 1), (-1, 2), (-2, 1), (1, -2), (2, -1), (-1, -2), (-2, -1)]


# ------------------------------------------------------------------ grafo

def build_board(n: int):
    """Grafo do cavalo em n x n. Vertices indexados r*n + c."""
    adj: List[List[int]] = [[] for _ in range(n * n)]
    edges: List[Tuple[int, int]] = []
    for r in range(n):
        for c in range(n):
            i = r * n + c
            for dr, dc in MOVES:
                rr, cc = r + dr, c + dc
                if 0 <= rr < n and 0 <= cc < n:
                    j = rr * n + cc
                    adj[i].append(j)
                    if i < j:
                        edges.append((i, j))
    for a in adj:
        a.sort()
    edges.sort()
    eidx = {e: k for k, e in enumerate(edges)}
    return adj, edges, eidx


def is_corner(n: int, v: int) -> bool:
    r, c = divmod(v, n)
    return r in (0, n - 1) and c in (0, n - 1)


def hexagons_bulk(n: int, adj) -> List[Tuple[int, ...]]:
    """6-ciclos simples que nao tocam canto (= hexagonos de Punc n)."""
    total = n * n
    ok = [not is_corner(n, v) for v in range(total)]
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


# ------------------------------------------------------------- mascaras

def cycle_mask(cyc: Sequence[int], eidx) -> int:
    x = 0
    for a, b in zip(cyc, tuple(cyc[1:]) + (cyc[0],)):
        x ^= 1 << eidx[(min(a, b), max(a, b))]
    return x


def alternating_matchings(cyc: Sequence[int], eidx) -> Tuple[int, int]:
    """Os dois emparelhamentos perfeitos alternados do ciclo par `cyc`."""
    ebits = []
    for a, b in zip(cyc, tuple(cyc[1:]) + (cyc[0],)):
        ebits.append(1 << eidx[(min(a, b), max(a, b))])
    m0 = 0
    m1 = 0
    for k, bit in enumerate(ebits):
        if k % 2 == 0:
            m0 |= bit
        else:
            m1 |= bit
    return m0, m1


def tour_mask(seq: Sequence[int], eidx) -> int:
    x = 0
    for a, b in zip(seq, tuple(seq[1:]) + (seq[0],)):
        x ^= 1 << eidx[(min(a, b), max(a, b))]
    return x


def is_hamiltonian(mask: int, n: int, edges) -> bool:
    """`mask` ja' e' 2-regular por construcao; testa se e' UM unico ciclo."""
    total = n * n
    nb: List[List[int]] = [[] for _ in range(total)]
    m = mask
    while m:
        b = m & -m
        k = b.bit_length() - 1
        u, v = edges[k]
        nb[u].append(v)
        nb[v].append(u)
        m ^= b
    prev = -1
    cur = 0
    cnt = 0
    while True:
        a, b = nb[cur]
        nxt = a if a != prev else b
        prev, cur = cur, nxt
        cnt += 1
        if cur == 0:
            break
        if cnt > total:
            return False
    return cnt == total


# ------------------------------------------------------------ union-find

class DSU:
    def __init__(self, n: int):
        self.p = list(range(n))
        self.r = [0] * n
        self.comps = n

    def find(self, x: int) -> int:
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x

    def union(self, a: int, b: int) -> bool:
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return False
        if self.r[ra] < self.r[rb]:
            ra, rb = rb, ra
        self.p[rb] = ra
        if self.r[ra] == self.r[rb]:
            self.r[ra] += 1
        self.comps -= 1
        return True


# ---------------------------------------------------------------- GF(2)

class Basis:
    def __init__(self):
        self.piv: Dict[int, int] = {}

    def reduce(self, x: int) -> int:
        while x:
            h = x.bit_length() - 1
            p = self.piv.get(h)
            if p is None:
                return x
            x ^= p
        return 0

    def add(self, x: int) -> bool:
        r = self.reduce(x)
        if r:
            self.piv[r.bit_length() - 1] = r
            return True
        return False

    @property
    def rank(self) -> int:
        return len(self.piv)


def zbulk_rank(n: int, edges, eidx, adj) -> int:
    """dim Z_bulk = |E(Punc)| - |V| + #comps(Punc) = |E|-8 - n^2 + 5."""
    return (len(edges) - 8) - n * n + 5


# ------------------------------------------------------------------ main

def run(n: int, tours: np.ndarray, out_path: Path, label: str) -> dict:
    t0 = time.time()
    adj, edges, eidx = build_board(n)
    print(f"[{label}] V={n*n}  E={len(edges)}  beta1={len(edges)-n*n+1}")

    hexes = hexagons_bulk(n, adj)
    hx = []
    for cyc in hexes:
        cm = cycle_mask(cyc, eidx)
        m0, m1 = alternating_matchings(cyc, eidx)
        hx.append((cm, m0, m1))
    print(f"[{label}] hexagonos do bulk: {len(hexes)}  ({time.time()-t0:.1f}s)")

    # tours -> mascaras, dedup
    masks: List[int] = []
    seen: Dict[int, int] = {}
    dup = 0
    for seq in tours:
        m = tour_mask([int(x) for x in seq], eidx)
        if m in seen:
            dup += 1
            continue
        seen[m] = len(masks)
        masks.append(m)
    print(f"[{label}] tours: {len(tours)} lidos, {len(masks)} distintos "
          f"(como conjunto de arestas), {dup} duplicados")

    dsu = DSU(len(masks))
    deg = [0] * len(masks)
    realized: Dict[int, int] = {}   # mascara do hexagono -> contagem de flips
    inside = 0    # pares (tau, C) que passam no teste de emparelhamento
    disconn = 0   # ... e falham conexidade
    external = 0  # ... e caem fora do conjunto de tours conhecido (fase 2)

    for i, tm in enumerate(masks):
        for cm, m0, m1 in hx:
            inter = tm & cm
            if inter != m0 and inter != m1:
                continue
            inside += 1
            nm = tm ^ cm
            if not is_hamiltonian(nm, n, edges):
                disconn += 1
                continue
            deg[i] += 1
            realized[cm] = realized.get(cm, 0) + 1
            j = seen.get(nm)
            if j is None:
                external += 1
                continue
            dsu.union(i, j)

    dist = Counter(deg)
    iso = [i for i, d in enumerate(deg) if d == 0]

    # componentes
    sizes = Counter(dsu.find(i) for i in range(len(masks)))
    comp_sizes = sorted(sizes.values(), reverse=True)

    # (U1)+(U2): os hexagonos realizados geram Z_bulk?
    B = Basis()
    for cm in realized:
        B.add(cm)
    dimZb = zbulk_rank(n, edges, eidx, adj)

    # controle: TODOS os hexagonos do bulk geram Z_bulk?
    Ball = Basis()
    for cm, _, _ in hx:
        Ball.add(cm)

    res = {
        "n": n,
        "V": n * n,
        "E": len(edges),
        "beta1": len(edges) - n * n + 1,
        "dim_ZBulk": dimZb,
        "n_hexagons_bulk": len(hexes),
        "n_tours_read": int(len(tours)),
        "n_tours_distinct": len(masks),
        "pairs_matching_ok": inside,
        "pairs_disconnected": disconn,
        "pairs_external": external,
        "flip_degree_min": min(deg) if deg else None,
        "flip_degree_max": max(deg) if deg else None,
        "flip_degree_mean": (sum(deg) / len(deg)) if deg else None,
        "flip_degree_hist": dict(sorted(dist.items())),
        "n_isolated": len(iso),
        "isolated_examples": iso[:10],
        "n_components": len(comp_sizes),
        "component_sizes": comp_sizes[:20],
        "n_hexagons_realized": len(realized),
        "rank_realized_hexagons": B.rank,
        "rank_all_bulk_hexagons": Ball.rank,
        "realized_span_ZBulk": B.rank == dimZb,
        "all_hex_span_ZBulk": Ball.rank == dimZb,
        "seconds": round(time.time() - t0, 2),
    }
    out_path.write_text(json.dumps(res, indent=2))
    return res


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=6)
    ap.add_argument("--tours", type=str,
                    default=str(ROOT.parent / "complex_orbit" / "data"
                               / "tours_closed_6x6.npy"))
    ap.add_argument("--out", type=str, default=None)
    args = ap.parse_args()

    tours = np.load(args.tours)
    out = Path(args.out) if args.out else DATA / f"flip_graph_n{args.n}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    res = run(args.n, tours, out, f"n={args.n}")

    print()
    print("=" * 66)
    print(f" GRAFO DE FLIPS HEXAGONAIS  n={args.n}")
    print("=" * 66)
    print(f" tours distintos            : {res['n_tours_distinct']}")
    print(f" hexagonos do bulk          : {res['n_hexagons_bulk']}")
    print(f" pares com matching ok      : {res['pairs_matching_ok']}")
    print(f"   -> desconexos (rejeitados): {res['pairs_disconnected']}")
    print(f"   -> fora do conjunto       : {res['pairs_external']}")
    print()
    print(f" grau de flip  min/med/max  : {res['flip_degree_min']} / "
          f"{res['flip_degree_mean']:.2f} / {res['flip_degree_max']}")
    print(f" tours isolados (grau 0)    : {res['n_isolated']}")
    print()
    print(f" COMPONENTES                : {res['n_components']}")
    print(f"   maiores                  : {res['component_sizes']}")
    print()
    print(f" dim Z_bulk                 : {res['dim_ZBulk']}")
    print(f" rank(hexagonos realizados) : {res['rank_realized_hexagons']}"
          f"   -> gera Z_bulk? {res['realized_span_ZBulk']}")
    print(f" rank(todos hex. do bulk)   : {res['rank_all_bulk_hexagons']}"
          f"   -> gera Z_bulk? {res['all_hex_span_ZBulk']}")
    print(f" hexagonos realizados       : {res['n_hexagons_realized']}"
          f" / {res['n_hexagons_bulk']}")
    print()
    print(f" tempo: {res['seconds']}s   -> {out}")


if __name__ == "__main__":
    main()
