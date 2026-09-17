#!/usr/bin/env python3
"""
flip_graph_n8.py
================
Fase 2 do teste de conectividade do grafo de flips: n=8, SEM enumeracao
exaustiva (N(8) ~ 1.3e13).

Duas medidas que sao testaveis sem o conjunto completo:

  (2a) CONECTIVIDADE LOCAL. BFS por flips a partir de UM tour semente.
       Mede |bola(k)| em funcao de k. Crescimento exponencial e' evidencia
       de conectividade; platou precoce e' evidencia contra.
       A bola tambem da' um certificado construtivo: todo hexagono usado
       nela e' C = H_A XOR H_B com ambos hamiltonianos.

  (2b) GRAU DE FLIP. Distribuicao sobre duas amostras com vieses DIFERENTES:
         - tours da bola de BFS       (viesada para grau alto, por construcao)
         - tours de Warnsdorff aleatorizado (viesada, mas de outro jeito)
       Se as duas distribuicoes coincidirem, a estatistica e' robusta ao
       viés; se divergirem, nenhuma das duas vale como estimativa.
       Grau 0 em QUALQUER das duas e' refutacao valida (achar e' conclusivo,
       nao achar nao e').

Nenhum corte de busca e' usado: todos os hexagonos do bulk sao testados em
todo tour. O unico corte e' o tamanho da bola, e ele e' reportado.

CLI:
  python flip_graph_n8.py --n 8 --ball 20000 --sample 2000
"""
from __future__ import annotations

import argparse
import json
import random
import time
from collections import Counter, deque
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from flip_graph import (
    Basis, DSU, alternating_matchings, build_board, cycle_mask,
    hexagons_bulk, is_hamiltonian, tour_mask,
)

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"


# ------------------------------------------------- amostrador independente

def warnsdorff_tour(n: int, adj, rng: random.Random,
                    tries: int = 4000) -> Optional[List[int]]:
    """Tour FECHADO por Warnsdorff com desempate aleatorio e reinicio.

    ATENCAO: amostra NAO uniforme. Serve para procurar contraexemplos
    (grau 0) e para comparar com o viés oposto da bola de BFS; nao serve
    como estimador de frequencia.
    """
    total = n * n
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
            cand = [u for u in cand if
                    sum(1 for w in adj[u] if not used[w]) == best]
            nxt = rng.choice(cand)
            path.append(nxt)
            used[nxt] = True
        if ok and start in adj[path[-1]]:
            return path
    return None


# ------------------------------------------------------------------ nucleo

def flip_neighbours(tm: int, hx, n: int, edges) -> List[Tuple[int, int]]:
    """Vizinhos de `tm` no grafo de flips. Retorna (mascara_nova, hexagono)."""
    out = []
    for cm, m0, m1 in hx:
        inter = tm & cm
        if inter != m0 and inter != m1:
            continue
        nm = tm ^ cm
        if is_hamiltonian(nm, n, edges):
            out.append((nm, cm))
    return out


def run(n: int, ball: int, sample: int, seed: int) -> dict:
    t0 = time.time()
    rng = random.Random(seed)
    adj, edges, eidx = build_board(n)
    beta1 = len(edges) - n * n + 1
    dimZb = (len(edges) - 8) - n * n + 5
    print(f"[n={n}] V={n*n} E={len(edges)} beta1={beta1} dim_ZBulk={dimZb}")

    hexes = hexagons_bulk(n, adj)
    hx = [(cycle_mask(c, eidx), *alternating_matchings(c, eidx)) for c in hexes]
    print(f"[n={n}] hexagonos do bulk: {len(hexes)}  ({time.time()-t0:.1f}s)")

    seed_path = warnsdorff_tour(n, adj, rng)
    if seed_path is None:
        raise RuntimeError("nao achei tour semente")
    root = tour_mask(seed_path, eidx)

    # ---------------------------------------------------------- (2a) BFS
    seen: Dict[int, int] = {root: 0}
    order: List[int] = [root]
    q = deque([root])
    cur_layer = 0
    realized: Dict[int, int] = {}
    deg_ball: List[int] = []
    truncated = False
    t_bfs = time.time()
    while q:
        if len(seen) >= ball:
            truncated = True
            break
        m = q.popleft()
        d = seen[m]
        if d > cur_layer:
            cur_layer = d
            print(f"  [bfs] camada {d}: |bola|={len(seen)} "
                  f"({time.time()-t_bfs:.0f}s)", flush=True)
        nb = flip_neighbours(m, hx, n, edges)
        deg_ball.append(len(nb))
        for nm, cm in nb:
            realized[cm] = realized.get(cm, 0) + 1
            if nm not in seen:
                seen[nm] = d + 1
                order.append(nm)
                q.append(nm)
    # tamanho por camada, recalculado de forma limpa
    lay = Counter(seen.values())
    layers = [lay[k] for k in range(max(lay) + 1)]

    B = Basis()
    for cm in realized:
        B.add(cm)

    # ------------------------------------------------- (2b) amostra externa
    deg_ext: List[int] = []
    ext_in_ball = 0
    t_s = time.time()
    for k in range(sample):
        p = warnsdorff_tour(n, adj, rng)
        if p is None:
            continue
        tm = tour_mask(p, eidx)
        if tm in seen:
            ext_in_ball += 1
        deg_ext.append(len(flip_neighbours(tm, hx, n, edges)))
        if (k + 1) % 500 == 0:
            print(f"  [amostra] {k+1}/{sample} ({time.time()-t_s:.0f}s)",
                  flush=True)

    def stats(v: List[int]) -> dict:
        if not v:
            return {}
        v2 = sorted(v)
        return {
            "n": len(v), "min": v2[0], "max": v2[-1],
            "mean": round(sum(v2) / len(v2), 3),
            "median": v2[len(v2) // 2],
            "zeros": sum(1 for x in v2 if x == 0),
            "hist": dict(sorted(Counter(v2).items())),
        }

    res = {
        "n": n, "V": n * n, "E": len(edges), "beta1": beta1,
        "dim_ZBulk": dimZb,
        "n_hexagons_bulk": len(hexes),
        "ball_cap": ball,
        "ball_truncated": truncated,
        "ball_size": len(seen),
        "ball_layers": layers,
        "ball_depth_reached": max(seen.values()),
        "n_hexagons_realized_in_ball": len(realized),
        "rank_realized_in_ball": B.rank,
        "realized_span_ZBulk": B.rank == dimZb,
        "degree_ball": stats(deg_ball),
        "degree_external_sample": stats(deg_ext),
        "external_hits_inside_ball": ext_in_ball,
        "external_sample_requested": sample,
        "seed": seed,
        "seconds": round(time.time() - t0, 2),
    }
    return res


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=8)
    ap.add_argument("--ball", type=int, default=20000)
    ap.add_argument("--sample", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=20260808)
    ap.add_argument("--out", type=str, default=None)
    args = ap.parse_args()

    res = run(args.n, args.ball, args.sample, args.seed)
    out = Path(args.out) if args.out else DATA / f"flip_graph_ball_n{args.n}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(res, indent=2))

    print()
    print("=" * 66)
    print(f" BOLA DE FLIPS  n={args.n}")
    print("=" * 66)
    print(f" |bola| = {res['ball_size']}"
          f"   (corte {res['ball_cap']}, truncada={res['ball_truncated']})")
    print(f" camadas: {res['ball_layers'][:15]}")
    db, de = res["degree_ball"], res["degree_external_sample"]
    print(f" grau na bola      : min={db.get('min')} med={db.get('mean')} "
          f"max={db.get('max')}  zeros={db.get('zeros')}  (N={db.get('n')})")
    print(f" grau amostra ext. : min={de.get('min')} med={de.get('mean')} "
          f"max={de.get('max')}  zeros={de.get('zeros')}  (N={de.get('n')})")
    print(f" amostras externas que cairam dentro da bola: "
          f"{res['external_hits_inside_ball']}/{de.get('n')}")
    print()
    print(f" dim Z_bulk               : {res['dim_ZBulk']}")
    print(f" rank(hex. realizados)    : {res['rank_realized_in_ball']}"
          f"  -> gera Z_bulk? {res['realized_span_ZBulk']}")
    print(f" hexagonos realizados     : {res['n_hexagons_realized_in_ball']}"
          f" / {res['n_hexagons_bulk']}")
    print(f" tempo: {res['seconds']}s  -> {out}")


if __name__ == "__main__":
    main()
