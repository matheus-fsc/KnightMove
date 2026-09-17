#!/usr/bin/env python3
"""
hexagon_locality.py
===================
Decide entre as rotas A/B/C para "hexagonos geram Z_bulk para todo n".

Tres testes, todos baratos:

(A) UNIFORMIDADE DO ENUNCIADO. Os hexagonos geram Z_1 do grafo do cavalo
    induzido num subretangulo m x m ARBITRARIO (sem remover cantos)?
    Se sim, o enunciado a provar deixa de mencionar bulk/cantos e vira uma
    propriedade local do grafo do cavalo -- muito mais facil de induzir.

(B) LOCALIDADE. Todo 4-ciclo (e todo ciclo curto) e' soma de hexagonos
    contidos numa JANELA de raio r em volta dele? O r minimo uniforme em n
    e' o que da' uniformidade de graca; sem isso, inducao nao fecha.

(C) PASSO INDUTIVO n -> n+2. Decompondo o tabuleiro N=n+2 em copia interna
    n x n + moldura, quanto do rank vem de cada parte? Se
    rank(interna) = dim Z_1(interna) e o resto fecha por contagem, o passo
    da Rota A e' numericamente viavel.

CLI:
  python hexagon_locality.py
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Dict, List, Set, Tuple

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"

MOVES = [(1, 2), (2, 1), (-1, 2), (-2, 1), (1, -2), (2, -1), (-1, -2), (-2, -1)]


def build_sub(cells: List[Tuple[int, int]]):
    """Grafo do cavalo induzido no conjunto de celulas `cells`."""
    idx = {c: i for i, c in enumerate(sorted(cells))}
    adj: List[List[int]] = [[] for _ in idx]
    edges: List[Tuple[int, int]] = []
    for (r, c), i in idx.items():
        for dr, dc in MOVES:
            j = idx.get((r + dr, c + dc))
            if j is None:
                continue
            adj[i].append(j)
            if i < j:
                edges.append((i, j))
    for a in adj:
        a.sort()
    return idx, adj, sorted(edges)


def cycles_len(adj, total: int, length: int) -> List[Tuple[int, ...]]:
    """Ciclos simples de comprimento `length`, cada um uma vez."""
    out = []
    for v1 in range(total):
        def rec(path: List[int]):
            cur = path[-1]
            if len(path) == length:
                if v1 in adj[cur] and path[1] < path[-1]:
                    out.append(tuple(path))
                return
            for u in adj[cur]:
                if u <= v1 or u in path:
                    continue
                path.append(u)
                rec(path)
                path.pop()
        rec([v1])
    return out


class Basis:
    """Base escalonada GF(2) sobre inteiros-mascara."""

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

    def contains(self, x: int) -> bool:
        return self.reduce(x) == 0

    @property
    def rank(self) -> int:
        return len(self.piv)


def cyc_vec(cyc, eidx) -> int:
    x = 0
    for a, b in zip(cyc, cyc[1:] + cyc[:1]):
        x ^= 1 << eidx[(min(a, b), max(a, b))]
    return x


# ---------------------------------------------------------------- teste (A)

def test_uniform(ms: List[int]) -> List[dict]:
    print("=" * 70)
    print(" (A) Hexagonos geram Z_1 do grafo do cavalo em m x m (sem cantos)?")
    print("=" * 70)
    out = []
    for m in ms:
        cells = [(r, c) for r in range(m) for c in range(m)]
        idx, adj, edges = build_sub(cells)
        eidx = {e: i for i, e in enumerate(edges)}
        V, E = len(idx), len(edges)
        # componentes
        seen, comps = set(), 0
        for s in range(V):
            if s in seen:
                continue
            comps += 1
            st = [s]
            seen.add(s)
            while st:
                u = st.pop()
                for w in adj[u]:
                    if w not in seen:
                        seen.add(w)
                        st.append(w)
        dimZ = E - V + comps
        b = Basis()
        for h in cycles_len(adj, V, 6):
            b.add(cyc_vec(h, eidx))
        ok = b.rank == dimZ
        print(f"   m={m:>3}: V={V:>4} E={E:>4} comp={comps}  "
              f"dim Z_1={dimZ:>4}  rank(hex)={b.rank:>4}  "
              f"{'GERA' if ok else 'FALTA ' + str(dimZ - b.rank)}")
        out.append({"m": m, "V": V, "E": E, "components": comps,
                    "dim_Z1": dimZ, "rank_hex": b.rank, "spans": ok})
    return out


# ---------------------------------------------------------------- teste (B)

def test_locality(n: int, radii: List[int]) -> dict:
    """4-ciclos do tabuleiro n x n sao soma de hexagonos numa janela de raio r?"""
    print("\n" + "=" * 70)
    print(f" (B) Localidade em {n}x{n}: 4-ciclos vs hexagonos de janela raio r")
    print("=" * 70)
    cells = [(r, c) for r in range(n) for c in range(n)]
    idx, adj, edges = build_sub(cells)
    rev = {i: c for c, i in idx.items()}
    eidx = {e: i for i, e in enumerate(edges)}
    V = len(idx)

    quads = cycles_len(adj, V, 4)
    hexs = cycles_len(adj, V, 6)
    print(f"   4-ciclos={len(quads)}  hexagonos={len(hexs)}")

    hex_cells = [set(rev[v] for v in h) for h in hexs]
    hex_vec = [cyc_vec(h, eidx) for h in hexs]

    res = {"n": n, "n_quads": len(quads), "n_hex": len(hexs), "by_radius": {}}
    unresolved = list(range(len(quads)))
    for r in radii:
        still = []
        for qi in unresolved:
            q = quads[qi]
            qc = [rev[v] for v in q]
            r0 = min(x for x, _ in qc) - r
            r1 = max(x for x, _ in qc) + r
            c0 = min(y for _, y in qc) - r
            c1 = max(y for _, y in qc) + r
            b = Basis()
            for k, hc in enumerate(hex_cells):
                if all(r0 <= x <= r1 and c0 <= y <= c1 for x, y in hc):
                    b.add(hex_vec[k])
            if not b.contains(cyc_vec(q, eidx)):
                still.append(qi)
        frac = 1 - len(still) / len(quads)
        print(f"   r={r}: resolvidos {len(quads) - len(still)}/{len(quads)} "
              f"({frac:.1%})")
        res["by_radius"][r] = {"resolved": len(quads) - len(still),
                               "total": len(quads)}
        unresolved = still
        if not unresolved:
            break
    res["min_radius_all"] = (None if unresolved
                             else max(int(k) for k in res["by_radius"]))
    return res


# ---------------------------------------------------------------- teste (C)

def test_induction(n: int) -> dict:
    """Tabuleiro N=n+2: quanto do rank vem dos hexagonos da copia interna?"""
    N = n + 2
    cells = [(r, c) for r in range(N) for c in range(N)]
    idx, adj, edges = build_sub(cells)
    rev = {i: c for c, i in idx.items()}
    eidx = {e: i for i, e in enumerate(edges)}
    V, E = len(idx), len(edges)

    inner = {(r, c) for r in range(1, N - 1) for c in range(1, N - 1)}
    _, adj_in, edges_in = build_sub(sorted(inner))
    V_in, E_in = len(inner), len(edges_in)
    dimZ_in = E_in - V_in + 1
    dimZ_N = E - V + 1

    hexs = cycles_len(adj, V, 6)
    b_in, b_all = Basis(), Basis()
    n_in = 0
    for h in hexs:
        x = cyc_vec(h, eidx)
        if all(rev[v] in inner for v in h):
            n_in += 1
            b_in.add(x)
    for h in hexs:
        b_all.add(cyc_vec(h, eidx))

    print(f"   n={n:>3} -> N={N:>3}: dim Z_1(N)={dimZ_N:>4}  "
          f"dim Z_1(interna)={dimZ_in:>4}  "
          f"rank(hex internos)={b_in.rank:>4} "
          f"{'=' if b_in.rank == dimZ_in else '!='}  "
          f"rank(todos)={b_all.rank:>4} "
          f"{'=' if b_all.rank == dimZ_N else '!='}  "
          f"salto={dimZ_N - dimZ_in}")
    return {"n": n, "N": N, "dim_Z1_N": dimZ_N, "dim_Z1_inner": dimZ_in,
            "n_hex_inner": n_in, "rank_hex_inner": b_in.rank,
            "rank_hex_all": b_all.rank,
            "inner_spans": b_in.rank == dimZ_in,
            "all_spans": b_all.rank == dimZ_N,
            "jump": dimZ_N - dimZ_in}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ms", type=str, default="5,6,7,8,9,10,11,12")
    ap.add_argument("--loc-n", type=int, default=10)
    ap.add_argument("--radii", type=str, default="0,1,2,3")
    ap.add_argument("--ind-ns", type=str, default="6,8,10,12")
    a = ap.parse_args()

    t0 = time.time()
    res = {}
    res["uniform"] = test_uniform([int(x) for x in a.ms.split(",")])
    res["locality"] = test_locality(a.loc_n, [int(x) for x in a.radii.split(",")])
    print("\n" + "=" * 70)
    print(" (C) Passo indutivo n -> n+2: copia interna + moldura")
    print("=" * 70)
    res["induction"] = [test_induction(int(x)) for x in a.ind_ns.split(",")]
    res["elapsed_s"] = round(time.time() - t0, 1)

    DATA.mkdir(parents=True, exist_ok=True)
    p = DATA / "hexagon_locality.json"
    p.write_text(json.dumps(res, indent=2))
    print(f"\n-> {p}   ({res['elapsed_s']}s)")


if __name__ == "__main__":
    main()
