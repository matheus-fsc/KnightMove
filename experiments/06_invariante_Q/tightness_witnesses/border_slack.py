#!/usr/bin/env python3
"""
border_slack.py
===============
As formas que geram as 6 direcoes perdidas no 5x8 tornam-se realizaveis
quando se acrescenta UMA linha (5x8 -> 6x8)?

A filtracao de `missing_directions.py` mostrou que a hipotese "nao cabe" e'
FALSA: hexagonos de altura 3 (num tabuleiro de 5 linhas) ja' geram as 6
direcoes. O que discrimina e' o contato com as bordas do lado CURTO --
excluir hexagonos que tocam linha 0 ou linha R-1 derruba 6 -> 2.

Este script testa a hipotese revisada: nao e' a forma que nao cabe, e' o
*switcher* que nao tem folga do lado de fora quando o hexagono esta' encostado
na borda curta.

Metodo: busca DIRIGIDA por certificado positivo, nao amostragem. Para um
hexagono C com emparelhamentos alternados m0/m1, procura-se um tour que
contenha m0 e evite m1 (entao tau ∩ C = m0 exatamente); se achar, testa-se
tau XOR C hamiltoniano. Isso e' um CERTIFICADO quando encontra.

⚠️ Quando NAO encontra dentro do orcamento de nos, o resultado e' "nao
encontrado no orcamento", NAO "provado inexistente" -- salvo em `--exhaustive`,
onde a busca satura (e a saturacao e' reportada explicitamente).

Auto-validacao: rodar com --rows 5 --cols 8 --exhaustive deve reproduzir
exatamente a classificacao de `missing_directions.py` (435/151/2).

CLI:
  python border_slack.py --rows 5 --cols 8 --exhaustive     # validacao
  python border_slack.py --embed 5x8 --into 6x8             # o teste
"""
from __future__ import annotations

import argparse
import json
import time
from collections import Counter
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

from flip_graph import alternating_matchings, cycle_mask
from flip_graph_rect import build_rect, hexagons_bulk_rect

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"


# ------------------------------------------------------- busca com restricoes

class Searcher:
    """Busca de ciclos hamiltonianos com arestas forcadas e proibidas."""

    def __init__(self, R: int, C: int, adj, edges, eidx):
        self.R, self.C, self.V = R, C, R * C
        self.adj, self.edges, self.eidx = adj, edges, eidx

    def find(self, forced: Sequence[Tuple[int, int]],
             banned: Sequence[Tuple[int, int]],
             node_budget: int, accept) -> Tuple[int | None, int, bool, int]:
        """Procura um tour compativel que satisfaca `accept`.

        Nao basta parar no primeiro tour compativel: o XOR dele com o hexagono
        pode desconectar enquanto o de outro tour nao desconecta. A busca so'
        para quando `accept` aceita, ou quando esgota/estoura o orcamento.

        Retorna (mascara aceita | None, nos, saturou, n_compativeis_vistos).
        """
        V, adj = self.V, self.adj
        ban = {frozenset(e) for e in banned}
        forced_nbr: Dict[int, List[int]] = {}
        for u, v in forced:
            forced_nbr.setdefault(u, []).append(v)
            forced_nbr.setdefault(v, []).append(u)
        # forced e' um emparelhamento (m0 = 3 arestas disjuntas)
        assert all(len(w) <= 1 for w in forced_nbr.values())

        nbrs = [[u for u in adj[v] if frozenset((v, u)) not in ban]
                for v in range(V)]
        used = [False] * V
        used[0] = True
        path = [0]
        nodes = [0]
        found: List[int] = []
        n_compat = [0]

        def viable(cur: int) -> bool:
            free = [v for v in range(V) if not used[v]]
            if not free:
                return True
            for v in free:
                k = 0
                for u in nbrs[v]:
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
                for u in nbrs[x]:
                    if not used[u] and u not in seen:
                        seen.add(u)
                        stack.append(u)
                        cnt += 1
            return cnt == len(free)

        def rec() -> bool:
            """True se achou (ou estourou o orcamento)."""
            nodes[0] += 1
            if nodes[0] > node_budget:
                return True
            cur = path[-1]
            prev = path[-2] if len(path) > 1 else -1
            if len(path) == V:
                if 0 in nbrs[cur]:
                    # a aresta de fechamento respeita o forcado em cur e em 0?
                    fc = forced_nbr.get(cur)
                    f0 = forced_nbr.get(0)
                    ok = True
                    if fc and fc[0] != prev and fc[0] != 0:
                        ok = False
                    if f0 and f0[0] != path[1] and f0[0] != cur:
                        ok = False
                    if ok:
                        m = 0
                        for a, b in zip(path, path[1:] + path[:1]):
                            m ^= 1 << self.eidx[(min(a, b), max(a, b))]
                        n_compat[0] += 1
                        if accept(m):
                            found.append(m)
                            return True
                return False

            # regra do forcado: se cur tem aresta forcada e nao chegamos por
            # ela, temos de sair por ela
            fc = forced_nbr.get(cur)
            if fc and fc[0] != prev:
                cands = [fc[0]] if not used[fc[0]] else []
            else:
                cands = [u for u in nbrs[cur] if not used[u]]

            for u in cands:
                # nao entre num vertice por aresta nao-forcada se ele ja' tem
                # de usar a forcada dos dois lados -- coberto pela regra acima
                used[u] = True
                path.append(u)
                if viable(u) and rec():
                    return True
                path.pop()
                used[u] = False
            return False

        rec()
        saturated = nodes[0] <= node_budget
        return (found[0] if found else None), nodes[0], saturated, n_compat[0]


def hamiltonian(mask: int, V: int, edges) -> bool:
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
        if cur == 0:
            return cnt == V
        if cnt > V:
            return False


def try_realize(S: Searcher, cyc, budget: int) -> dict:
    """Tenta certificar que o hexagono `cyc` e' flip-realizavel."""
    eidx, edges, V = S.eidx, S.edges, S.V
    cm = cycle_mask(cyc, eidx)
    hedges = [(min(a, b), max(a, b))
              for a, b in zip(cyc, tuple(cyc[1:]) + (cyc[0],))]
    e0 = [hedges[k] for k in range(6) if k % 2 == 0]
    e1 = [hedges[k] for k in range(6) if k % 2 == 1]

    total_nodes = 0
    total_compat = 0
    sat = True
    accept = lambda m: hamiltonian(m ^ cm, V, edges)
    for forced, banned in ((e0, e1), (e1, e0)):
        tm, nodes, s, nc = S.find(forced, banned, budget, accept)
        total_nodes += nodes
        total_compat += nc
        sat = sat and s
        if tm is not None:
            return {"mode": "realized", "nodes": total_nodes,
                    "saturated": sat, "n_compatible": total_compat,
                    "tour": tm}
    if not sat:
        return {"mode": "not_found_in_budget", "nodes": total_nodes,
                "saturated": False, "n_compatible": total_compat}
    # busca saturou sem aceitar: a distincao e' entao PROVADA
    return {"mode": "no_match" if total_compat == 0 else "disconnect",
            "nodes": total_nodes, "saturated": True,
            "n_compatible": total_compat}


# ------------------------------------------------------------------ mapeamento

def embed(cellset, dr: int) -> List[Tuple[int, int]]:
    return [(r + dr, c) for r, c in cellset]


def cyc_of_cells(cl, C: int, adj) -> Tuple[int, ...] | None:
    """Reordena as 6 celulas num ciclo do grafo do cavalo, se existir."""
    vs = [r * C + c for r, c in cl]
    vset = set(vs)
    start = vs[0]
    path = [start]
    seen = {start}

    def rec():
        cur = path[-1]
        if len(path) == 6:
            return start in adj[cur]
        for u in adj[cur]:
            if u in vset and u not in seen:
                path.append(u)
                seen.add(u)
                if rec():
                    return True
                seen.discard(u)
                path.pop()
        return False

    return tuple(path) if rec() else None


# ----------------------------------------------------------------------- main

def validate(R: int, C: int, budget: int) -> dict:
    adj, edges, eidx = build_rect(R, C)
    S = Searcher(R, C, adj, edges, eidx)
    hexes = hexagons_bulk_rect(R, C, adj)
    t0 = time.time()
    modes = Counter()
    unsat = 0
    for i, cyc in enumerate(hexes):
        r = try_realize(S, cyc, budget)
        modes[r["mode"]] += 1
        if not r["saturated"]:
            unsat += 1
        if (i + 1) % 100 == 0:
            print(f"  {i+1}/{len(hexes)}  {dict(modes)} "
                  f"({time.time()-t0:.0f}s)", flush=True)
    return {"R": R, "C": C, "n_hexagons": len(hexes),
            "modes": dict(modes), "n_unsaturated": unsat,
            "budget": budget, "seconds": round(time.time() - t0, 1)}


def embed_test(src: Tuple[int, int], dst: Tuple[int, int],
               budget: int) -> dict:
    R0, C0 = src
    R1, C1 = dst
    src_json = DATA / f"missing_directions_{R0}x{C0}.json"
    data = json.loads(src_json.read_text())

    adj1, edges1, eidx1 = build_rect(R1, C1)
    S = Searcher(R1, C1, adj1, edges1, eidx1)

    # os hexagonos nao-realizados e nao-nulos no quociente sao os interessantes;
    # o JSON guarda os 6 geradores explicitamente, e usamos TODOS os
    # nao-realizados para nao depender da escolha gulosa.
    adj0, edges0, eidx0 = build_rect(R0, C0)
    hexes0 = hexagons_bulk_rect(R0, C0, adj0)

    # reconstitui a classificacao a partir do JSON (por celulas)
    gens = [[tuple(x) for x in g["cells"]] for g in data["generators"]]

    meta = [(g.get("interior", False), g.get("touches", {}))
            for g in data["generators"]]

    out = []
    for gi, cl in enumerate(gens):
        rowspan = max(r for r, _ in cl) - min(r for r, _ in cl) + 1
        entry = {"generator": gi, "cells_src": [list(x) for x in cl],
                 "bbox_rows": rowspan, "interior_in_src": meta[gi][0],
                 "placements": []}
        for dr in range(0, R1 - R0 + 1):
            for dc in range(0, C1 - C0 + 1):
                cl2 = [(r + dr, c + dc) for r, c in cl]
                cyc = cyc_of_cells(cl2, C1, adj1)
                if cyc is None:
                    entry["placements"].append(
                        {"dr": dr, "dc": dc, "mode": "not_a_cycle"})
                    continue
                r = try_realize(S, cyc, budget)
                r.pop("tour", None)
                r["dr"], r["dc"] = dr, dc
                r["cells"] = [list(x) for x in cl2]
                entry["placements"].append(r)
                # imprime na hora: se um timeout externo matar o processo, os
                # resultados parciais sobrevivem (com `| tail` perdiam-se todos)
                sat = "" if r.get("saturated", True) else " (NAO saturou)"
                print(f"   g[{gi}]{'*' if meta[gi][0] else ' '} dr={dr} dc={dc}"
                      f": {r['mode']} nos={r.get('nodes','-')}{sat}", flush=True)
        out.append(entry)
    return {"src": f"{R0}x{C0}", "dst": f"{R1}x{C1}",
            "budget": budget, "generators": out}


def span_test(R: int, C: int, budget: int) -> dict:
    """Os hexagonos CERTIFICADOS de R x C ja' geram Z_bulk?

    Unilateral na direcao segura: `span{tau+tau'} ⊇ span(realizados)` sempre,
    logo se os certificados ja' geram Z_bulk entao `deficit = 3` -- e isso vale
    como PROVA, sem enumerar os tours do tabuleiro. O contrario nao vale: se
    nao geram, pode ser so' orcamento.
    """
    from coset_tests import Basis as IntBasis, Ctx, zbulk_basis

    adj, edges, eidx = build_rect(R, C)
    S = Searcher(R, C, adj, edges, eidx)
    hexes = hexagons_bulk_rect(R, C, adj)
    Zb = IntBasis(zbulk_basis(Ctx(R, C)))
    dim_zb = Zb.rank
    print(f"[{R}x{C}] hex={len(hexes)}  dim Z_bulk={dim_zb}", flush=True)

    B = IntBasis()
    modes = Counter()
    t0 = time.time()
    for i, cyc in enumerate(hexes):
        cm = cycle_mask(cyc, eidx)
        assert Zb.contains(cm), "hexagono do bulk fora de Z_bulk"
        # se ja' esta' no span, nao gasta busca
        if B.contains(cm):
            modes["skipped_in_span"] += 1
            continue
        r = try_realize(S, cyc, budget)
        modes[r["mode"]] += 1
        if r["mode"] == "realized":
            B.add(cm)
        if (i + 1) % 100 == 0:
            print(f"  {i+1}/{len(hexes)}  rank={B.rank}/{dim_zb}  "
                  f"{dict(modes)} ({time.time()-t0:.0f}s)", flush=True)
        if B.rank == dim_zb:
            print(f"  ATINGIU dim Z_bulk em i={i+1}", flush=True)
            break
    assert B.rank <= dim_zb
    return {
        "R": R, "C": C, "n_hexagons": len(hexes), "dim_ZBulk": dim_zb,
        "rank_certified": B.rank,
        "certified_span_ZBulk": B.rank == dim_zb,
        "deficit_proved_3": B.rank == dim_zb,
        "modes": dict(modes), "budget": budget,
        "seconds": round(time.time() - t0, 1),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rows", type=int, default=5)
    ap.add_argument("--cols", type=int, default=8)
    ap.add_argument("--exhaustive", action="store_true")
    ap.add_argument("--budget", type=int, default=2_000_000)
    ap.add_argument("--embed", default="")
    ap.add_argument("--span", default="")
    ap.add_argument("--into", default="")
    args = ap.parse_args()

    if args.embed:
        R0, C0 = map(int, args.embed.lower().split("x"))
        R1, C1 = map(int, args.into.lower().split("x"))
        res = embed_test((R0, C0), (R1, C1), args.budget)
        out = DATA / f"border_slack_{R0}x{C0}_into_{R1}x{C1}.json"
        out.write_text(json.dumps(res, indent=2))
        print("=" * 70)
        print(f" EMBEDDING  {R0}x{C0} -> {R1}x{C1}   (orcamento "
              f"{args.budget:,} nos)")
        print("=" * 70)
        for g in res["generators"]:
            tag = "INTERIOR" if g.get("interior_in_src") else "borda"
            print(f" gerador [{g['generator']}] {tag}  altura bbox="
                  f"{g['bbox_rows']}  cells={g['cells_src']}")
            for p in g["placements"]:
                sat = "" if p.get("saturated", True) else "  (NAO saturou)"
                print(f"     dr={p['dr']} dc={p['dc']}: {p['mode']}"
                      f"  nos={p.get('nodes','-')}{sat}")
        print(f" -> {out}")
        return

    if args.span:
        R1, C1 = map(int, args.span.lower().split("x"))
        res = span_test(R1, C1, args.budget)
        out = DATA / f"border_slack_span_{R1}x{C1}.json"
        out.write_text(json.dumps(res, indent=2))
        print("=" * 70)
        print(f" SPAN DOS CERTIFICADOS  {R1}x{C1}")
        print("=" * 70)
        for k, v in res.items():
            print(f" {k:24s}: {v}")
        print(f" -> {out}")
        return

    budget = 10**18 if args.exhaustive else args.budget
    res = validate(args.rows, args.cols, budget)
    out = DATA / f"border_slack_validate_{args.rows}x{args.cols}.json"
    out.write_text(json.dumps(res, indent=2))
    print("=" * 70)
    print(f" VALIDACAO {args.rows}x{args.cols}  "
          f"{'(EXAUSTIVO)' if args.exhaustive else f'(orcamento {budget:,})'}")
    print("=" * 70)
    for k, v in res.items():
        print(f" {k:20s}: {v}")
    print(f" -> {out}")


if __name__ == "__main__":
    main()
