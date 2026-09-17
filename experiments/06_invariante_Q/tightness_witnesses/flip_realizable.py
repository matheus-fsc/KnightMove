#!/usr/bin/env python3
"""
flip_realizable.py
==================
Os 12 hexagonos excepcionais de n=6 sao efeito de tamanho finito?

Teste: para um hexagono C (num tabuleiro n x n qualquer) e um dos seus dois
emparelhamentos alternados m, existe um tour que CONTEM m e EVITA C\\m?

Isso e' uma busca hamiltoniana com arestas forcadas e proibidas. Ela e'
exaustiva por padrao (sem orcamento), logo:
  - achar tour  => C e' matching-realizavel  (conclusivo)
  - esgotar     => C NAO e' matching-realizavel (conclusivo)
  - estourar o orcamento => INCONCLUSIVO, e e' reportado como tal.

Depois, para os que passam, testa tambem se algum desses tours da' um flip
completo (tau XOR C conexo).

CLI:
  python flip_realizable.py --check6          # valida contra a verdade exaustiva
  python flip_realizable.py --n 8 --budget 3000000
"""
from __future__ import annotations

import argparse
import json
import time
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple

import numpy as np

from flip_graph import (
    alternating_matchings, build_board, cycle_mask, hexagons_bulk,
    is_hamiltonian, tour_mask,
)

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"


class Budget(Exception):
    pass


def constrained_tours(n: int, adj, forced: Dict[int, int],
                      forbidden: Set[Tuple[int, int]],
                      budget: int, want: int):
    """Tours hamiltonianos fechados com arestas forcadas/proibidas.

    `forced` mapeia v -> w para cada uma das 6 pontas do emparelhamento
    (relacao simetrica). `forbidden` sao pares (min,max).

    Retorna (lista_de_tours, exhausted). `exhausted=False` significa que a
    busca PAROU CEDO -- por orcamento de nos OU por ter atingido `want`
    solucoes -- e portanto nao viu todas as solucoes.

    ATENCAO: parar em `want` tambem quebra a exaustividade. Um `flip=False`
    depois de `want` tours encontrados NAO e' prova de que C nao seja
    flip-realizavel; so' diz que nenhum dos `want` primeiros serviu.
    """
    total = n * n
    allowed = [[u for u in adj[v] if (min(v, u), max(v, u)) not in forbidden]
               for v in range(total)]
    used = [False] * total
    path: List[int] = []
    out: List[List[int]] = []
    nodes = 0
    start = 0

    def viable() -> bool:
        """Conectividade + grau residual das casas nao usadas."""
        cur = path[-1]
        free = [v for v in range(total) if not used[v]]
        if not free:
            return True
        # grau: cada casa livre precisa de >=2 ligacoes entre {livres, cur, start}
        for v in free:
            k = 0
            for u in allowed[v]:
                if not used[u] or u == cur or u == start:
                    k += 1
                    if k >= 2:
                        break
            if k < 2:
                return False
        # conectividade a partir de cur pelas casas livres
        seen = {cur}
        stack = [cur]
        cnt = 0
        while stack:
            x = stack.pop()
            for u in allowed[x]:
                if not used[u] and u not in seen:
                    seen.add(u)
                    stack.append(u)
                    cnt += 1
        return cnt == len(free)

    def rec() -> None:
        nonlocal nodes
        nodes += 1
        if nodes > budget:
            raise Budget
        cur = path[-1]
        if len(path) == total:
            if start not in allowed[cur]:
                return
            # verificacao final direta sobre o conjunto de arestas do ciclo
            es = {(min(a, b), max(a, b))
                  for a, b in zip(path, path[1:] + path[:1])}
            if es & forbidden:
                return
            for v, w2 in forced.items():
                if (min(v, w2), max(v, w2)) not in es:
                    return
            out.append(list(path))
            return
        prev = path[-2] if len(path) >= 2 else None
        w = forced.get(cur)
        if w is not None and w != prev:
            cand = [w] if not used[w] else []
        else:
            cand = [u for u in allowed[cur] if not used[u]]
            # Warnsdorff
            cand.sort(key=lambda u: sum(1 for z in allowed[u] if not used[z]))
        for u in cand:
            # se `u` tem forcada e nao e' `cur`, so' podemos entrar em `u` se
            # a forcada de `u` ainda estiver disponivel para a saida
            fu = forced.get(u)
            if fu is not None and fu != cur and used[fu]:
                continue
            used[u] = True
            path.append(u)
            if viable():
                rec()
                if len(out) >= want:
                    used[u] = False
                    path.pop()
                    return
            path.pop()
            used[u] = False

    used[start] = True
    path.append(start)
    try:
        rec()
        exhausted = True
    except Budget:
        exhausted = False
    # atingir `want` tambem interrompe a varredura: nao e' exaustivo
    if len(out) >= want:
        exhausted = False
    return out, exhausted


def ham_mask(mask: int, V: int, edges) -> bool:
    """`mask` 2-regular: e' um unico ciclo cobrindo os V vertices?"""
    nb: List[List[int]] = [[] for _ in range(V)]
    m = mask
    while m:
        b = m & -m
        u, v = edges[b.bit_length() - 1]
        nb[u].append(v)
        nb[v].append(u)
        m ^= b
    if any(len(x) != 2 for x in nb):
        return False
    prev, cur, cnt = -1, 0, 0
    while True:
        a, b2 = nb[cur]
        nxt = a if a != prev else b2
        prev, cur = cur, nxt
        cnt += 1
        if cur == 0 or cnt > V:
            break
    return cnt == V


def search_flip(n: int, adj, edges, eidx, cyc, budget: int,
                rng=None) -> Tuple[bool, bool]:
    """Versao quadrada: V = n*n. Ver `search_flip_V`."""
    return search_flip_V(n * n, adj, edges, eidx, cyc, budget, rng)


def search_flip_V(V: int, adj, edges, eidx, cyc, budget: int,
                  rng=None) -> Tuple[bool, bool]:
    """Procura tour tau com tau∩C alternado E tau XOR C conexo.

    Testa o flip INLINE em cada solucao, sem acumular lista, e para na
    primeira que serve. Sem `want`: ou acha, ou varre tudo, ou estoura o
    orcamento.

    Se `rng` for dado, a ordem dos candidatos e' embaralhada -- necessario
    porque o DFS em ordem fixa entrega solucoes que compartilham quase todo o
    prefixo, e portanto falham todas pelo mesmo motivo.

    Retorna (achou_flip, exhausted).
    """
    cm = cycle_mask(cyc, eidx)
    cyc_edges = [(min(a, b), max(a, b))
                 for a, b in zip(cyc, tuple(cyc[1:]) + (cyc[0],))]
    total = V
    exhausted_all = True
    for m in alternating_matchings(cyc, eidx):
        forced: Dict[int, int] = {}
        forb: Set[Tuple[int, int]] = set()
        for e in cyc_edges:
            if m >> eidx[e] & 1:
                forced[e[0]] = e[1]
                forced[e[1]] = e[0]
            else:
                forb.add(e)
        allowed = [[u for u in adj[v] if (min(v, u), max(v, u)) not in forb]
                   for v in range(total)]
        used = [False] * total
        path: List[int] = [0]
        used[0] = True
        nodes = 0
        found = False

        def viable() -> bool:
            cur = path[-1]
            free = [v for v in range(total) if not used[v]]
            if not free:
                return True
            for v in free:
                k = 0
                for u in allowed[v]:
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
                for u in allowed[x]:
                    if not used[u] and u not in seen:
                        seen.add(u)
                        stack.append(u)
                        cnt += 1
            return cnt == len(free)

        def rec() -> None:
            nonlocal nodes, found
            if found:
                return
            nodes += 1
            if nodes > budget:
                raise Budget
            cur = path[-1]
            if len(path) == total:
                if 0 not in allowed[cur]:
                    return
                tm = tour_mask(path, eidx)
                if (tm & cm) != m:
                    return
                if ham_mask(tm ^ cm, V, edges):
                    found = True
                return
            prev = path[-2] if len(path) >= 2 else None
            w = forced.get(cur)
            if w is not None and w != prev:
                cand = [w] if not used[w] else []
            else:
                cand = [u for u in allowed[cur] if not used[u]]
                if rng is not None:
                    rng.shuffle(cand)
                else:
                    cand.sort(key=lambda u:
                              sum(1 for z in allowed[u] if not used[z]))
            for u in cand:
                fu = forced.get(u)
                if fu is not None and fu != cur and used[fu]:
                    continue
                used[u] = True
                path.append(u)
                if viable():
                    rec()
                path.pop()
                used[u] = False
                if found:
                    return

        try:
            rec()
            ex = True
        except Budget:
            ex = False
        if found:
            return True, True
        exhausted_all = exhausted_all and ex
    return False, exhausted_all


def test_hexagon(n: int, adj, edges, eidx, cyc, budget: int, want: int) -> dict:
    cm = cycle_mask(cyc, eidx)
    ms = alternating_matchings(cyc, eidx)
    cyc_edges = [(min(a, b), max(a, b))
                 for a, b in zip(cyc, tuple(cyc[1:]) + (cyc[0],))]
    res = {"matching": False, "flip": False, "exhausted": True, "n_tours": 0}
    for m in ms:
        fmask = m
        forced: Dict[int, int] = {}
        forb: Set[Tuple[int, int]] = set()
        for e in cyc_edges:
            k = eidx[e]
            if fmask >> k & 1:
                forced[e[0]] = e[1]
                forced[e[1]] = e[0]
            else:
                forb.add(e)
        tours, ex = constrained_tours(n, adj, forced, forb, budget, want)
        res["exhausted"] = res["exhausted"] and ex
        res["n_tours"] += len(tours)
        if tours:
            res["matching"] = True
            for t in tours:
                if is_hamiltonian(tour_mask(t, eidx) ^ cm, n, edges):
                    res["flip"] = True
                    break
        if res["flip"]:
            break
    return res


# ------------------------------------------------------------ formas dos 12

def exceptional_shapes_n6():
    """As formas de translacao dos 12 excepcionais de n=6, e a orbita D4."""
    n = 6
    adj, edges, eidx = build_board(n)
    hexes = hexagons_bulk(n, adj)
    tours = np.load(ROOT.parent / "complex_orbit" / "data"
                    / "tours_closed_6x6.npy")
    masks = [tour_mask([int(x) for x in s], eidx) for s in tours]
    bad = []
    for cyc in hexes:
        cm = cycle_mask(cyc, eidx)
        m0, m1 = alternating_matchings(cyc, eidx)
        if not any((tm & cm) in (m0, m1) for tm in masks):
            bad.append(cyc)
    shapes = set()
    for cyc in bad:
        cells = [divmod(v, n) for v in cyc]
        r0 = min(r for r, _ in cells)
        c0 = min(c for _, c in cells)
        shapes.add(tuple(sorted((r - r0, c - c0) for r, c in cells)))
    return sorted(shapes)


def cells_to_cycle(cells: Sequence[Tuple[int, int]], n: int, adj):
    """Ordena as 6 casas num 6-ciclo do cavalo, se possivel."""
    vs = [r * n + c for r, c in cells]
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


def run_scale(shapes, ns: List[int], budget: int, want: int) -> dict:
    out = {}
    for n in ns:
        adj, edges, eidx = build_board(n)
        rows = []
        t0 = time.time()
        for si, sh in enumerate(shapes):
            R = max(r for r, _ in sh) + 1
            C = max(c for _, c in sh) + 1
            placed = 0
            m_ok = 0
            f_ok = 0
            incon = 0
            for r0 in range(n - R + 1):
                for c0 in range(n - C + 1):
                    cells = [(r + r0, c + c0) for r, c in sh]
                    # so' no bulk: nenhuma casa pode ser canto
                    if any((r in (0, n - 1)) and (c in (0, n - 1))
                           for r, c in cells):
                        continue
                    cyc = cells_to_cycle(cells, n, adj)
                    if cyc is None:
                        continue
                    placed += 1
                    res = test_hexagon(n, adj, edges, eidx, cyc, budget, want)
                    m_ok += res["matching"]
                    f_ok += res["flip"]
                    incon += (not res["exhausted"]) and (not res["matching"])
            rows.append({"shape": sh, "bbox": [R, C], "placements": placed,
                         "matching_realizable": m_ok, "flip_realizable": f_ok,
                         "inconclusive": incon})
            print(f"  [n={n}] forma {si+1}/{len(shapes)} bbox {R}x{C}: "
                  f"{placed} posicoes -> matching {m_ok}, flip {f_ok}, "
                  f"inconclusivo {incon}  ({time.time()-t0:.0f}s)", flush=True)
        out[str(n)] = rows
    return out


def check6(budget: int) -> None:
    """Valida o buscador restrito contra a verdade exaustiva de n=6."""
    n = 6
    adj, edges, eidx = build_board(n)
    hexes = hexagons_bulk(n, adj)
    tours = np.load(ROOT.parent / "complex_orbit" / "data"
                    / "tours_closed_6x6.npy")
    masks = [tour_mask([int(x) for x in s], eidx) for s in tours]
    truth = {}
    for cyc in hexes:
        cm = cycle_mask(cyc, eidx)
        m0, m1 = alternating_matchings(cyc, eidx)
        truth[cyc] = any((tm & cm) in (m0, m1) for tm in masks)
    disagree = 0
    incon = 0
    t0 = time.time()
    for k, cyc in enumerate(hexes):
        r = test_hexagon(n, adj, edges, eidx, cyc, budget, 1)
        if not r["exhausted"] and not r["matching"]:
            incon += 1
            continue
        if r["matching"] != truth[cyc]:
            disagree += 1
            print("  DIVERGENCIA em", [divmod(v, 6) for v in cyc])
    print(f"[check6] {len(hexes)} hexagonos: divergencias={disagree} "
          f"inconclusivos={incon}  ({time.time()-t0:.0f}s)")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check6", action="store_true")
    ap.add_argument("--ns", type=int, nargs="*", default=[8, 10, 12])
    ap.add_argument("--budget", type=int, default=2_000_000)
    ap.add_argument("--want", type=int, default=40)
    args = ap.parse_args()

    if args.check6:
        check6(args.budget)
        return

    shapes = exceptional_shapes_n6()
    print(f"formas de translacao dos 12 excepcionais: {len(shapes)}")
    res = run_scale(shapes, args.ns, args.budget, args.want)
    out = DATA / "flip_exceptional_scaling.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(res, indent=2, default=str))

    print()
    print("=" * 70)
    print(" AS 10 FORMAS DOS 12 EXCEPCIONAIS, COLOCADAS EM n MAIOR")
    print("=" * 70)
    for n, rows in res.items():
        P = sum(r["placements"] for r in rows)
        M = sum(r["matching_realizable"] for r in rows)
        F = sum(r["flip_realizable"] for r in rows)
        I = sum(r["inconclusive"] for r in rows)
        print(f" n={n:>2}: {P:4d} posicoes -> matching-realizavel {M:4d}"
              f"   flip-realizavel {F:4d}   inconclusivo {I:4d}")
    print(f"\n -> {out}")


if __name__ == "__main__":
    main()
