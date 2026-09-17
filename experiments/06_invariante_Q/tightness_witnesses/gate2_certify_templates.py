#!/usr/bin/env python3
"""
gate2_certify_templates.py
==========================
Gate 1 (2a metade) + Gate 2 do PLAN.md.

Gate 1 mostrou que 5 FORMAS de hexagono, fechadas por translacao, geram
Z_1(G_m) para m = 8,10,12 -- as MESMAS 5 formas em todo m, todas cabendo numa
caixa 4x4.

Aqui responde-se a pergunta que decide (U1):

    a certificacao (existencia dos dois tours H_A, H_B com H_A ^ H_B = C)
    depende SO' DA FORMA, ou depende da posicao no tabuleiro?

Se depender so' da forma (a menos de um efeito de borda que some no interior),
entao (U1) colapsa para "verificar 5 templates + argumento de translacao".

Para cada translacao de cada template mede-se:
  - certificado (sim/nao)
  - `inset` = distancia do bounding box a' borda do tabuleiro
e reporta-se a taxa de certificacao por inset. A previsao a testar e':
certificacao = 100% para inset >= d0, com d0 constante em n.

CLI:
  python gate2_certify_templates.py --ns 10,12
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Dict, List, Set, Tuple

from switcher_search import (build, build_switchers, color, corners, ham_path,
                             verify_switcher)

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"

# as 5 formas do Gate 1 (identicas em m = 8, 10, 12)
TEMPLATES = [
    ((0, 0), (0, 1), (1, 2), (1, 3), (2, 0), (2, 1)),
    ((0, 0), (0, 2), (1, 0), (1, 2), (2, 1), (3, 1)),
    ((0, 1), (0, 2), (1, 0), (1, 3), (2, 1), (2, 2)),
    ((0, 0), (0, 2), (1, 2), (2, 1), (2, 3), (3, 1)),
    ((0, 2), (0, 3), (1, 0), (1, 1), (2, 2), (2, 3)),
]


def cycle_from_cells(adj, cells: List[int]) -> Tuple[int, ...] | None:
    """Ordena as 6 celulas no ciclo hamiltoniano do subgrafo induzido."""
    S = set(cells)
    start = min(S)
    order = [start]
    seen = {start}

    def rec() -> bool:
        if len(order) == 6:
            return start in adj[order[-1]]
        for u in adj[order[-1]]:
            if u in S and u not in seen:
                seen.add(u)
                order.append(u)
                if rec():
                    return True
                order.pop()
                seen.discard(u)
        return False

    return tuple(order) if rec() else None


def certify(adj, n: int, hexa: Tuple[int, ...], cs: Set[int],
            max_internal: int, max_w: int, budget: int) -> bool:
    total = n * n
    cands: List[dict] = []
    for rot in range(3):
        cands += build_switchers(adj, hexa, rot, max_internal, max_w)
    cands.sort(key=lambda w: w["VW"])
    for w in cands:
        v1, v4 = w["poles"]
        removed = (set(w["cycle"]) | set(w["P2_internal"])
                   | set(w["P3_internal"])) - {v1, v4}
        allowed = set(range(total)) - removed
        ca = sum(1 for x in allowed if color(x, n) == 0)
        if not (color(v1, n) != color(v4, n) and ca * 2 == len(allowed)):
            continue
        if any(sum(1 for u in adj[c] if u in allowed) < 2 for c in cs):
            continue
        p, _ex, _nd = ham_path(adj, allowed, v1, v4, budget)
        if p is not None and not verify_switcher(adj, w, p, allowed, total):
            return True
    return False


def run(n: int, max_internal: int, max_w: int, budget: int) -> dict:
    t0 = time.time()
    adj = build(n)
    cs = set(corners(n))
    res = {"n": n, "templates": []}

    print(f"\n{'=' * 72}\n n = {n}\n{'=' * 72}")
    for ti, shape in enumerate(TEMPLATES):
        H = max(r for r, _ in shape)
        W = max(c for _, c in shape)
        by_inset: Dict[int, List[int]] = {}
        n_tr = n_ok = 0
        fails: List[dict] = []
        for r0 in range(n - H):
            for c0 in range(n - W):
                cells = [(r0 + r) * n + (c0 + c) for r, c in shape]
                if any(v in cs for v in cells):
                    continue
                hexa = cycle_from_cells(adj, cells)
                if hexa is None:          # forma nao realiza ciclo aqui
                    continue
                top, left = r0, c0
                bot, right = n - 1 - (r0 + H), n - 1 - (c0 + W)
                inset = min(top, left, bot, right)
                ok = certify(adj, n, hexa, cs, max_internal, max_w, budget)
                n_tr += 1
                n_ok += ok
                by_inset.setdefault(inset, [0, 0])
                by_inset[inset][1] += 1
                by_inset[inset][0] += ok
                if not ok:
                    # assinatura relativa ao canto MAIS PROXIMO: distancia a'
                    # cada par de lados + de qual lado se mede (a forma nao e'
                    # D4-simetrica, entao o lado importa)
                    fails.append({
                        "r0": r0, "c0": c0,
                        "top": top, "left": left, "bot": bot, "right": right,
                        "sig": [min(top, bot), 0 if top <= bot else 1,
                                min(left, right), 0 if left <= right else 1],
                    })

        rates = {k: (v[0], v[1]) for k, v in sorted(by_inset.items())}
        d0 = None
        for k in sorted(rates):
            if all(rates[j][0] == rates[j][1] for j in sorted(rates) if j >= k):
                d0 = k
                break
        print(f"  T{ti} {shape}")
        print(f"     translacoes={n_tr}  certificadas={n_ok} "
              f"({n_ok / n_tr:.1%})   inset -> ok/total: "
              + "  ".join(f"{k}:{v[0]}/{v[1]}" for k, v in rates.items())
              + f"   d0={d0}")
        print(f"     falhas (assinatura canto-relativa): "
              + "  ".join(sorted(set(str(tuple(f["sig"])) for f in fails))))
        res["templates"].append({
            "index": ti, "shape": [list(x) for x in shape],
            "translates": n_tr, "certified": n_ok,
            "by_inset": {str(k): list(v) for k, v in rates.items()},
            "d0_all_certified_from": d0,
            "failures": fails,
            "failure_signatures": sorted(tuple(f["sig"]) for f in fails),
        })
    res["elapsed_s"] = round(time.time() - t0, 1)
    return res


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ns", type=str, default="10,12")
    ap.add_argument("--max-internal", type=int, default=4)
    ap.add_argument("--max-w", type=int, default=20)
    ap.add_argument("--budget", type=int, default=400_000)
    a = ap.parse_args()

    print("=" * 72)
    print(" Gate 2 — a certificacao depende so' da FORMA?")
    print("=" * 72)
    out = [run(n, a.max_internal, a.max_w, a.budget)
           for n in [int(x) for x in a.ns.split(",")]]

    DATA.mkdir(parents=True, exist_ok=True)
    p = DATA / f"gate2_certify_templates_n{a.ns.replace(',', '_')}.json"
    p.write_text(json.dumps(out, indent=2))
    print(f"\n-> {p}")


if __name__ == "__main__":
    main()
