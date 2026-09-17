#!/usr/bin/env python3
"""
switcher_search.py
==================
Fase 3, Tarefa 2: o parity-switcher do CNP e' construtivel em G_n?

Def. 2.2 (CNP, lida no PDF original): W e' um R-parity-switcher se consiste
de um ciclo par C = (v_1,...,v_2k) com numero impar de arestas em R, mais
caminhos VERTICE-DISJUNTOS P_i ligando v_i a v_{2k-i+2}, para 2 <= i <= k.
(Nota: o passo (S2.b) do artigo diz "edge-disjoint"; a Def. 2.2 diz
"vertex-disjoint". Usamos a definicao.)

Para k=3 (hexagono C = v1..v6): P2 liga v2<->v6, P3 liga v3<->v5.
Polos (grau 2 em W): v1 e v4. Todos os outros v_i tem grau 3 em W.

k deve ser IMPAR: v_1 e v_{k+1} estao a distancia k em C, e Conrad et al.
Thm 3.1 exige cores opostas (n par). k=3 e' o menor caso nao-degenerado.

Este script NAO depende de R: os passos (S2.b) e (S3) sao perguntas
puramente estruturais sobre G_n. E' isso que torna a fase testavel sem
exibir um R fora de X (que nao existe no alcance computavel).

Etapas:
  2a  para cada hexagono no bulk, existem P2 e P3 vertice-disjuntos?
  2b  para cada W valido, existe caminho hamiltoniano v1 -> v4 em
      G_n \\ (V(W) \\ {v1, v4})?   <-- teste decisivo

CLI:
  python switcher_search.py --n 8  --limit 400
  python switcher_search.py --n 10 --limit 200
  python switcher_search.py --n 12 --limit 100
"""
from __future__ import annotations

import argparse
import itertools
import json
import random
import time
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"

MOVES = [(2, 1), (2, -1), (-2, 1), (-2, -1),
         (1, 2), (1, -2), (-1, 2), (-1, -2)]


def build(n: int):
    total = n * n
    adj: List[List[int]] = [[] for _ in range(total)]
    for r in range(n):
        for c in range(n):
            v = r * n + c
            for dr, dc in MOVES:
                nr, nc = r + dr, c + dc
                if 0 <= nr < n and 0 <= nc < n:
                    adj[v].append(nr * n + nc)
    return [sorted(set(a)) for a in adj]


def label(v: int, n: int) -> str:
    r, c = divmod(v, n)
    return f"{chr(ord('A') + c)}{n - r}"


def color(v: int, n: int) -> int:
    r, c = divmod(v, n)
    return (r + c) & 1


def corners(n: int) -> List[int]:
    return [0, n - 1, n * (n - 1), n * n - 1]


# ── hexagonos no bulk ───────────────────────────────────────────────

def hexagons(adj, forbidden: Set[int], total: int) -> List[Tuple[int, ...]]:
    """
    Todos os 6-ciclos induzidos-ou-nao que evitam `forbidden`.
    Canonico: v1 = menor vertice do ciclo, e path[1] < path[-1].
    """
    out = []
    for v1 in range(total):
        if v1 in forbidden:
            continue
        # DFS de profundidade 6 voltando a v1
        def rec(path: List[int]):
            cur = path[-1]
            if len(path) == 6:
                if v1 in adj[cur] and path[1] < path[-1]:
                    out.append(tuple(path))
                return
            for u in adj[cur]:
                if u <= v1 or u in forbidden or u in path:
                    continue
                path.append(u)
                rec(path)
                path.pop()
        rec([v1])
    return out


# ── caminhos P_i vertice-disjuntos ──────────────────────────────────

def simple_paths(adj, s: int, t: int, blocked: Set[int],
                 max_internal: int) -> List[List[int]]:
    """Caminhos s->t cujos vertices internos evitam `blocked`, ordenados
    por comprimento crescente. Retorna a lista de internos."""
    res: List[List[int]] = []

    def rec(cur: int, internal: List[int]):
        if len(internal) > max_internal:
            return
        if cur == t:
            if internal:
                res.append(list(internal))
            return
        if len(internal) == max_internal:
            return
        for u in adj[cur]:
            if u == t:
                res.append(list(internal))
                continue
            if u in blocked or u in internal or u == s:
                continue
            internal.append(u)
            rec(u, internal)
            internal.pop()

    rec(s, [])
    # dedupe preservando ordem por tamanho
    seen = set()
    uniq = []
    for p in sorted(res, key=len):
        k = tuple(p)
        if k in seen:
            continue
        seen.add(k)
        uniq.append(p)
    return uniq


def build_switchers(adj, hexa: Tuple[int, ...], rot: int,
                    max_internal: int, max_w: int) -> List[dict]:
    """
    Para um hexagono e uma escolha de polos (rot = indice de v1 em 0..2),
    devolve ate max_w switchers W, dos menores |V(W)| para cima.
    """
    v = [hexa[(rot + i) % 6] for i in range(6)]      # v1..v6
    vc = set(v)
    out: List[dict] = []
    # P2: v2 <-> v6 ; P3: v3 <-> v5
    for int2 in simple_paths(adj, v[1], v[5], vc, max_internal):
        blocked = vc | set(int2)
        for int3 in simple_paths(adj, v[2], v[4], blocked, max_internal):
            out.append({
                "poles": [v[0], v[3]],
                "cycle": v,
                "P2_internal": int2,
                "P3_internal": int3,
                "VW": len(vc) + len(int2) + len(int3),
            })
            if len(out) >= max_w:
                return out
    return out


# ── (S3): caminho hamiltoniano com buraco ───────────────────────────

def ham_path(adj, allowed: Set[int], s: int, t: int,
             budget: int = 400_000) -> Tuple[Optional[List[int]], bool, int]:
    """
    Caminho hamiltoniano de s a t no subgrafo induzido por `allowed`.

    Retorna (caminho ou None, exhausted, nodes).
      exhausted=True  -> a busca TERMINOU: None significa PROVADO que nao existe
      exhausted=False -> estourou o orcamento: None significa "nao encontrado"
    """
    nbr = {v: [u for u in adj[v] if u in allowed] for v in allowed}
    need = len(allowed)
    visited: Set[int] = {s}
    path = [s]
    nodes = [0]
    hit_budget = [False]

    def prune(cur: int) -> bool:
        # (a) grau residual
        for x in allowed:
            if x in visited:
                continue
            avail = 0
            for u in nbr[x]:
                if u not in visited or u == cur:
                    avail += 1
                    if avail >= 2:
                        break
            floor = 1 if x == t else 2
            if avail < floor:
                return True
        # (b) conexidade do residual a partir de cur
        stack = [cur]
        seen = {cur}
        reach = 0
        while stack:
            x = stack.pop()
            for u in nbr[x]:
                if u in seen or u in visited:
                    continue
                seen.add(u)
                reach += 1
                stack.append(u)
        return reach != need - len(visited)

    def dfs(cur: int) -> bool:
        nodes[0] += 1
        if nodes[0] > budget:
            hit_budget[0] = True
            return False
        if len(path) == need:
            return cur == t
        if cur == t:                       # chegou cedo demais em t
            return False
        if prune(cur):
            return False
        cands = [u for u in nbr[cur] if u not in visited]
        cands.sort(key=lambda u: sum(1 for w in nbr[u] if w not in visited))
        for u in cands:
            visited.add(u)
            path.append(u)
            if dfs(u):
                return True
            path.pop()
            visited.discard(u)
        return False

    ok = dfs(s)
    return (list(path) if ok else None), (not hit_budget[0]), nodes[0]


# ── verificacao independente do objeto encontrado ───────────────────

def verify_switcher(adj, w: dict, hp: List[int], allowed: Set[int],
                    total: int) -> List[str]:
    """
    Confere, do zero, que (W, caminho hamiltoniano) satisfaz a Def. 2.2 + (S3).
    Retorna lista de problemas (vazia = tudo certo).
    """
    bad: List[str] = []
    v = w["cycle"]
    i2, i3 = w["P2_internal"], w["P3_internal"]

    # C e' um 6-ciclo de G
    if len(set(v)) != 6:
        bad.append("C tem vertices repetidos")
    for a, b in zip(v, v[1:] + v[:1]):
        if b not in adj[a]:
            bad.append(f"C: {a}-{b} nao e' aresta")

    # P2 liga v2..v6 ; P3 liga v3..v5 ; disjuncao
    for name, s, t, internal in (("P2", v[1], v[5], i2),
                                 ("P3", v[2], v[4], i3)):
        chain = [s] + internal + [t]
        for a, b in zip(chain, chain[1:]):
            if b not in adj[a]:
                bad.append(f"{name}: {a}-{b} nao e' aresta")
        if len(set(internal)) != len(internal):
            bad.append(f"{name}: internos repetidos")
        if set(internal) & set(v):
            bad.append(f"{name}: interno toca C")
    if set(i2) & set(i3):
        bad.append("P2 e P3 nao sao vertice-disjuntos")

    # os dois caminhos hamiltonianos de W (Def. 2.2, k=3)
    VW = set(v) | set(i2) | set(i3)
    pathA = [v[0], v[1]] + i2 + [v[5], v[4]] + i3[::-1] + [v[2], v[3]]
    pathB = [v[0], v[5]] + i2[::-1] + [v[1], v[2]] + i3 + [v[4], v[3]]
    for name, p in (("hamA(W)", pathA), ("hamB(W)", pathB)):
        if set(p) != VW or len(p) != len(VW):
            bad.append(f"{name} nao e' hamiltoniano em V(W)")
        for a, b in zip(p, p[1:]):
            if b not in adj[a]:
                bad.append(f"{name}: {a}-{b} nao e' aresta")
        if p[0] != v[0] or p[-1] != v[3]:
            bad.append(f"{name} nao vai de v1 a v4")
    # cada aresta de C em exatamente um dos dois caminhos
    def eset(p):
        return {frozenset((a, b)) for a, b in zip(p, p[1:])}
    cyc_edges = {frozenset((a, b)) for a, b in zip(v, v[1:] + v[:1])}
    ea, eb = eset(pathA), eset(pathB)
    for e in cyc_edges:
        if (e in ea) == (e in eb):
            bad.append("aresta de C nao esta em exatamente um dos caminhos")

    # (S3): hp e' caminho hamiltoniano de v1 a v4 no grafo com buraco
    if hp is not None:
        if hp[0] != v[0] or hp[-1] != v[3]:
            bad.append("(S3): extremos errados")
        if set(hp) != allowed or len(hp) != len(allowed):
            bad.append("(S3): nao cobre exatamente os vertices permitidos")
        for a, b in zip(hp, hp[1:]):
            if b not in adj[a]:
                bad.append(f"(S3): {a}-{b} nao e' aresta")
        if set(hp) & (VW - {v[0], v[3]}):
            bad.append("(S3): caminho invade o interior de W")
    return bad


# ── experimento ─────────────────────────────────────────────────────

def run(n: int, limit: int, max_internal: int, max_w: int,
        budget: int, seed: int) -> dict:
    t0 = time.time()
    adj = build(n)
    total = n * n
    cs = set(corners(n))
    rng = random.Random(seed)

    print(f"\n{'=' * 70}\n G_{n}x{n}  (V={total})\n{'=' * 70}")

    hexa = hexagons(adj, cs, total)
    print(f"[2a] hexagonos no bulk (evitando os 4 cantos): {len(hexa)}")

    sample = hexa if len(hexa) <= limit else rng.sample(hexa, limit)
    print(f"     testando {len(sample)} (amostra aleatoria, seed={seed})")

    n_with_W = 0
    per_hex: List[List[dict]] = []
    size_dist: Dict[int, int] = {}
    for h in sample:
        got = []
        for rot in range(3):                 # 3 pares de polos por hexagono
            got += build_switchers(adj, h, rot, max_internal, max_w)
        if got:
            n_with_W += 1
            got.sort(key=lambda w: w["VW"])
            per_hex.append(got)
            size_dist[got[0]["VW"]] = size_dist.get(got[0]["VW"], 0) + 1
    print(f"     hexagonos que admitem W: {n_with_W}/{len(sample)} "
          f"({n_with_W / max(1, len(sample)):.1%})")
    print(f"     |V(W)| minimo por hexagono: {dict(sorted(size_dist.items()))}")

    # ── (S3) ────────────────────────────────────────────────────────
    # Para cada hexagono tentamos os seus W em ordem de tamanho ate um
    # satisfazer (S3): o que interessa e' se o hexagono admite switcher
    # COMPLETO, nao se o menor W especifico funciona.
    print(f"\n[2b] (S3) caminho hamiltoniano v1 -> v4 em G \\ (V(W)\\{{v1,v4}}):")
    ok = hex_fail = 0
    n_color = n_corner = n_proved = n_budget = 0
    examples: List[dict] = []
    verify_problems: List[str] = []
    for cands in per_hex:
        solved = False
        for w in cands:
            v1, v4 = w["poles"]
            removed = (set(w["cycle"]) | set(w["P2_internal"])
                       | set(w["P3_internal"])) - {v1, v4}
            allowed = set(range(total)) - removed

            # filtro de cor
            ca = sum(1 for x in allowed if color(x, n) == 0)
            cb = len(allowed) - ca
            if not (color(v1, n) != color(v4, n) and ca == cb):
                n_color += 1
                continue
            # canto com < 2 vizinhos restantes e' fatal (canto tem grau 2)
            if any(sum(1 for u in adj[c] if u in allowed) < 2 for c in cs):
                n_corner += 1
                continue

            p, exhausted, nodes = ham_path(adj, allowed, v1, v4, budget)
            if p is not None:
                problems = verify_switcher(adj, w, p, allowed, total)
                if problems:
                    verify_problems += problems[:3]
                    continue
                solved = True
                if len(examples) < 3:
                    examples.append({
                        "cycle": [label(x, n) for x in w["cycle"]],
                        "P2_internal": [label(x, n) for x in w["P2_internal"]],
                        "P3_internal": [label(x, n) for x in w["P3_internal"]],
                        "VW": w["VW"], "nodes": nodes,
                        "ham_path_len": len(p),
                    })
                break
            elif exhausted:
                n_proved += 1
            else:
                n_budget += 1
        if solved:
            ok += 1
        else:
            hex_fail += 1

    print(f"     hexagonos com switcher COMPLETO : {ok}/{len(per_hex)} "
          f"({ok / max(1, len(per_hex)):.1%})")
    print(f"     hexagonos sem nenhum W completo : {hex_fail}")
    print(f"     (contagem por W: cor={n_color}, canto={n_corner}, "
          f"(S3) provado ausente={n_proved}, sem orcamento={n_budget})")
    if verify_problems:
        print(f"     !!! VERIFICACAO FALHOU: {verify_problems[:5]}")
    else:
        print(f"     verificacao independente de todos os W aceitos: OK")
    tested = len(per_hex)
    color_fail, corner_fail = n_color, n_corner
    fail_proved, fail_budget = n_proved, n_budget

    if examples:
        print(f"\n     exemplo de switcher completo:")
        e = examples[0]
        print(f"       C  = {e['cycle']}")
        print(f"       P2 = {e['P2_internal']}   P3 = {e['P3_internal']}")
        print(f"       |V(W)| = {e['VW']}, caminho ham. cobre "
              f"{e['ham_path_len']} vertices ({e['nodes']} nos)")

    res = {
        "n": n, "V": total,
        "n_hexagons_bulk": len(hexa), "n_sampled": len(sample),
        "hexagons_with_W": n_with_W,
        "hexagons_with_complete_switcher": ok,
        "hexagons_without_complete_switcher": hex_fail,
        "W_tested": tested,
        "s3_ok": ok, "s3_fail_proved": fail_proved,
        "s3_fail_budget": fail_budget,
        "rejected_color": color_fail, "rejected_corner": corner_fail,
        "examples": examples,
        "params": {"limit": limit, "max_internal": max_internal,
                   "max_w": max_w, "budget": budget, "seed": seed},
        "elapsed_s": round(time.time() - t0, 1),
    }
    DATA.mkdir(parents=True, exist_ok=True)
    out = DATA / f"switcher_n{n}.json"
    out.write_text(json.dumps(res, indent=2))
    print(f"\n-> {out}  ({res['elapsed_s']}s)")
    return res


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=8)
    ap.add_argument("--limit", type=int, default=400)
    ap.add_argument("--max-internal", type=int, default=3)
    ap.add_argument("--max-w", type=int, default=6)
    ap.add_argument("--budget", type=int, default=400_000)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    run(a.n, a.limit, a.max_internal, a.max_w, a.budget, a.seed)


if __name__ == "__main__":
    main()
