#!/usr/bin/env python3
"""
witness_map.py
==============
Fase 2 do plano de tightness (paper/notes/tightness_plan.md):
reconhecimento computacional dos *witnesses duais* R.

Definicoes (dicionario CNP <-> paper):
  C          = Z_1 = ker d_1              (espaco de ciclos, dim beta_1)
  C^perp     = row(d_1)                   (espaco de cortes, dim V-1)
  C_n        = Span(Ham(n))               (espaco gerado pelos tours)
  C_n^perp   = { R : |R ^ T| par para todo tour T }
  deficit    = dim(C_n^perp / C^perp) = beta_1 - rank(Ham)

Um *witness* e' um R em C_n^perp \\ C^perp: um funcional que anula todos os
tours mas nao e' um corte. Este script:

  1. calcula rank(Ham), beta_1, deficit;
  2. verifica que o quociente C_n^perp / C^perp tem exatamente as 8 classes
     previstas pelo teorema Q(n)=3 (0, seis pares de cantos, e a soma dos 4);
  3. ENUMERA TODOS os witnesses de peso <= W (meet-in-the-middle sobre
     sindromes), classifica cada um por classe do quociente e por
     localizacao (nivel de anel L, incidencia em canto);
  4. procura, por classe, um representante de peso MAXIMO (busca local) e
     mede a condicao (C3) do Lema 2.1 do CNP em cortes de vertice unico.

CLI:
  python witness_map.py --n 6 --max-weight 6
  python witness_map.py --n 8 --tours 4000 --max-weight 4
"""
from __future__ import annotations

import argparse
import itertools
import json
import random
import time
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np

from knight_gf2 import (Reducer, boundary_matrix, build_graph, corners,
                        edge_index, edge_label, edge_level, enumerate_tours,
                        gf2_rref, int_to_support, label, mandatory_edges,
                        ring_level, rows_to_ints, sample_tours, vec_to_int)

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"


# ── busca de baixo peso: meet-in-the-middle sobre sindromes ─────────

def low_weight_vectors(syndromes: Sequence[int], max_weight: int,
                       n_edges: int) -> Dict[int, List[Tuple[int, ...]]]:
    """
    Todos os R (como tuplas de indices de aresta) com |R| <= max_weight e
    sindrome nula, i.e. R em C_n^perp.

    Estrategia: todo subconjunto de tamanho w <= 2h parte-se em duas metades
    de tamanho <= h. Indexamos as metades pela sindrome e casamos colisoes.
    """
    h = max_weight // 2
    table: Dict[int, List[Tuple[int, ...]]] = {}
    for size in range(h + 1):
        for combo in itertools.combinations(range(n_edges), size):
            s = 0
            for e in combo:
                s ^= syndromes[e]
            table.setdefault(s, []).append(combo)

    found: Dict[int, List[Tuple[int, ...]]] = {}
    seen = set()
    for size in range(h + 1):
        for combo in itertools.combinations(range(n_edges), size):
            s = 0
            for e in combo:
                s ^= syndromes[e]
            partners = table.get(s)
            if not partners:
                continue
            aset = set(combo)
            for other in partners:
                if len(other) + size > max_weight:
                    continue
                if aset & set(other):
                    continue
                r = tuple(sorted(aset | set(other)))
                if not r or r in seen:
                    continue
                seen.add(r)
                found.setdefault(len(r), []).append(r)
    return found


# ── busca local por representante de peso maximo no coset ──────────

def max_weight_in_coset(R0: int, n_vertices: int, edges, adj,
                        restarts: int = 60, seed: int = 0) -> Tuple[int, int]:
    """
    Maximiza |R0 ^ delta(S)| sobre S subset V (delta(S) percorre C^perp).
    Busca local com flips de vertice; retorna (melhor R como int, peso).
    NAO e' certificado como otimo global (isto e' MAX-CUT com pesos +-1).
    """
    rng = random.Random(seed)
    inc: List[List[int]] = [[] for _ in range(n_vertices)]
    for i, (u, v) in enumerate(edges):
        inc[u].append(i)
        inc[v].append(i)

    best, best_w = R0, bin(R0).count("1")
    for t in range(restarts):
        cur = R0
        if t:
            for v in range(n_vertices):
                if rng.random() < 0.5:
                    for i in inc[v]:
                        cur ^= 1 << i
        improved = True
        while improved:
            improved = False
            order = list(range(n_vertices))
            rng.shuffle(order)
            for v in order:
                # flip de v: arestas incidentes trocam de lado
                delta = sum(1 if not ((cur >> i) & 1) else -1 for i in inc[v])
                if delta > 0:
                    for i in inc[v]:
                        cur ^= 1 << i
                    improved = True
        w = bin(cur).count("1")
        if w > best_w:
            best, best_w = cur, w
    return best, best_w


def c3_single_vertex(R: int, n_vertices: int, edges) -> Tuple[float, int, int]:
    """
    Condicao (C3) restrita a cortes de vertice unico S={v}:
        e_R(v) >= deg(v)/2.
    Retorna (razao minima, #vertices violando, V).
    """
    inc: List[List[int]] = [[] for _ in range(n_vertices)]
    for i, (u, v) in enumerate(edges):
        inc[u].append(i)
        inc[v].append(i)
    worst, viol = 1.0, 0
    for v in range(n_vertices):
        deg = len(inc[v])
        if not deg:
            continue
        hits = sum(1 for i in inc[v] if (R >> i) & 1)
        ratio = hits / deg
        worst = min(worst, ratio)
        if hits * 2 < deg:
            viol += 1
    return worst, viol, n_vertices


# ── analise principal ───────────────────────────────────────────────

def analyse(n: int, n_tours: int | None, max_weight: int,
            seed: int = 0) -> dict:
    t0 = time.time()
    adj, edges = build_graph(n)
    eidx = edge_index(edges)
    V, E = n * n, len(edges)
    beta1 = E - V + 1

    print(f"\n{'=' * 68}\n G_{n}x{n}:  V={V}  E={E}  beta_1={beta1}\n{'=' * 68}")

    # 1. tours -------------------------------------------------------
    if n_tours is None:
        tours = enumerate_tours(n)
        mode = "exaustivo"
    else:
        tours = sample_tours(n, n_tours, seed=seed)
        mode = f"amostrado (K={len(tours)})"
    print(f"[1] tours: {len(tours)}  ({mode}, {time.time() - t0:.1f}s)")

    Ham = np.zeros((len(tours), E), dtype=np.uint8)
    for i, t in enumerate(tours):
        Ham[i, list(t)] = 1
    ham_red, ham_piv, rank_ham = gf2_rref(Ham)
    deficit = beta1 - rank_ham
    print(f"    rank(Ham) = {rank_ham}   deficit = beta_1 - rank = {deficit}")

    # 2. espacos duais ----------------------------------------------
    B = boundary_matrix(edges, V)
    cut = Reducer(B)                      # C^perp = row(d_1), dim V-1
    dim_cnperp = E - rank_ham
    print(f"[2] dim C^perp = {cut.rank} (= V-1 = {V - 1})   "
          f"dim C_n^perp = E - rank = {dim_cnperp}")
    print(f"    dim(C_n^perp / C^perp) = {dim_cnperp - cut.rank}  "
          f"(deve ser = deficit = {deficit})")

    # 3. as 8 classes previstas por Q(n)=3 ---------------------------
    mand = mandatory_edges(n, adj, eidx)
    cs = corners(n)
    reps: Dict[str, int] = {}
    for k in range(0, 5, 2):                       # |S| par: 0, 2, 4
        for S in itertools.combinations(range(4), k):
            vec = vec_to_int([mand[cs[i]][0] for i in S])
            name = "0" if not S else "+".join(f"r{i + 1}" for i in S)
            reps[name] = vec
    residual_to_class = {}
    for name, vec in reps.items():
        res = cut.reduce(vec)
        if res in residual_to_class:
            raise SystemExit(f"classes {residual_to_class[res]} e {name} colidem")
        residual_to_class[res] = name
    print(f"[3] 8 classes previstas (Q(n)=3) sao distintas mod C^perp: OK")

    # cada rep previsto anula todos os tours?
    ham_rows = rows_to_ints(ham_red)
    def in_cnperp(x: int) -> bool:
        return all(bin(x & r).count("1") % 2 == 0 for r in ham_rows)
    assert all(in_cnperp(v) for v in reps.values()), "rep previsto nao e' witness"
    print(f"    todos os 8 reps estao em C_n^perp: OK")

    # 4. enumeracao exaustiva de witnesses de baixo peso -------------
    syndromes = [0] * E
    for j, row in enumerate(ham_red):
        for e in np.flatnonzero(row):
            syndromes[int(e)] |= 1 << j

    t1 = time.time()
    low = low_weight_vectors(syndromes, max_weight, E)
    n_low = sum(len(v) for v in low.values())
    print(f"[4] vetores de C_n^perp com peso <= {max_weight}: {n_low}"
          f"  ({time.time() - t1:.1f}s)")

    by_class: Dict[str, List[dict]] = {}
    unknown: List[dict] = []
    weight_hist: Dict[int, Dict[str, int]] = {}
    for w in sorted(low):
        for combo in low[w]:
            x = vec_to_int(combo)
            res = cut.reduce(x)
            cls = residual_to_class.get(res)
            corner_edges = {i for m in mand.values() for i in m}
            info = {
                "weight": w,
                "idx": list(combo),
                "edges": [edge_label(edges[i], n) for i in combo],
                "levels": [edge_level(edges[i], n) for i in combo],
                "touches_corner": bool(set(combo) & corner_edges),
            }
            if cls is None:
                info["class"] = "DESCONHECIDA"
                unknown.append(info)
                cls = "DESCONHECIDA"
            else:
                info["class"] = cls
            by_class.setdefault(cls, []).append(info)
            weight_hist.setdefault(w, {}).setdefault(cls, 0)
            weight_hist[w][cls] += 1

    print(f"\n    peso | classe -> contagem")
    for w in sorted(weight_hist):
        parts = ", ".join(f"{c}:{k}" for c, k in sorted(weight_hist[w].items()))
        print(f"    {w:>4} | {parts}")

    # 5. witnesses MINIMOS por classe -------------------------------
    print(f"\n[5] witness de peso minimo por classe do quociente:")
    minimal: Dict[str, dict] = {}
    for cls, lst in sorted(by_class.items()):
        if cls == "0":
            continue
        wmin = min(i["weight"] for i in lst)
        cands = [i for i in lst if i["weight"] == wmin]
        bulk = [i for i in cands if not i["touches_corner"]]
        minimal[cls] = {
            "min_weight": wmin,
            "n_minimal": len(cands),
            "n_minimal_bulk": len(bulk),
            "example": cands[0]["edges"],
            "example_levels": cands[0]["levels"],
            "bulk_example": bulk[0]["edges"] if bulk else None,
        }
        flag = "" if bulk else "   <-- TODOS tocam canto"
        print(f"    {cls:<12} |R|min={wmin}  #min={len(cands):<4} "
              f"#bulk={len(bulk):<4} ex={cands[0]['edges']}{flag}")

    # 5-bis. localizacao: existe ALGUM witness nao-trivial que evita os cantos?
    # (a Fase 3 do plano quer um parity-switcher vivendo no bulk)
    print(f"\n[5b] witnesses que EVITAM toda aresta incidente a canto, "
          f"por classe e peso:")
    bulk_stats: Dict[str, Dict[str, int]] = {}
    for cls, lst in sorted(by_class.items()):
        per_w: Dict[str, int] = {}
        for i in lst:
            if not i["touches_corner"]:
                per_w[str(i["weight"])] = per_w.get(str(i["weight"]), 0) + 1
        bulk_stats[cls] = per_w
        tot = sum(per_w.values())
        detail = ", ".join(f"w{w}:{c}" for w, c in sorted(per_w.items())) or "-"
        print(f"    {cls:<12} total={tot:<4} {detail}")
    nonzero_bulk = sum(sum(v.values()) for c, v in bulk_stats.items() if c != "0")
    if nonzero_bulk == 0:
        print(f"    => NENHUM witness nao-trivial de peso <= {max_weight} "
              f"evita os cantos.")

    if unknown:
        print(f"\n    !!! {len(unknown)} witnesses fora das 8 classes previstas "
              f"-- tightness FALHA neste n")
    else:
        print(f"\n    nenhum witness fora das 8 classes: consistente com "
              f"deficit={deficit}")

    # 6. representantes de peso maximo e condicao (C3) ---------------
    # Lema da complementacao: G_n e' bipartido (o cavalo alterna de cor) e
    # |V| e' par, logo o vetor de arestas cheio 1_E = delta(A) esta' em
    # C^perp.  Portanto R -> E \ R e' uma involucao de C_n^perp que PRESERVA
    # a classe mod C^perp, e em cada coset  max|R| = E - min|R|.
    # Consequencia: o representante maximal do Lema 2.1 do CNP e' exatamente
    # o COMPLEMENTO do detector minimo -- nao ha' nada a otimizar.
    full = (1 << E) - 1
    assert cut.contains(full), "1_E nao esta' em C^perp (grafo nao bipartido?)"
    print(f"\n[6] lema da complementacao: 1_E em C^perp (bipartido, |V| par): OK")
    print(f"    => em cada coset, max|R| = E - min|R|; o R maximal do CNP e'")
    print(f"       o complemento do detector minimo. Condicao (C3) medida nele:")

    maximal: Dict[str, dict] = {}
    for cls in sorted(minimal):
        base = min(by_class[cls], key=lambda i: i["weight"])
        Rmin = vec_to_int(base["idx"])
        Rmax = full ^ Rmin
        wmax = E - base["weight"]
        assert in_cnperp(Rmax) and cut.reduce(Rmax) == cut.reduce(Rmin)
        # confirmacao independente por busca local (deve empatar, nao superar)
        _, w_local = max_weight_in_coset(Rmin, V, edges, adj, seed=seed)
        assert w_local <= wmax, "busca local superou o limite da complementacao"
        ratio, viol, _ = c3_single_vertex(Rmax, V, edges)
        maximal[cls] = {
            "max_weight": wmax,
            "equals_complement_of_min": True,
            "local_search_best": w_local,
            "fraction_of_E": round(wmax / E, 4),
            "c3_min_ratio_single_vertex": round(ratio, 4),
            "c3_violating_vertices": viol,
            "certified_optimal": True,
        }
        print(f"    {cls:<12} |R|max={wmax:<4} ({wmax / E:.1%} de E)  "
              f"(C3) min e_R(v)/deg(v) = {ratio:.3f}  violacoes = {viol}"
              f"  [busca local: {w_local}]")

    result = {
        "n": n, "V": V, "E": E, "beta1": beta1,
        "n_tours": len(tours), "tour_mode": mode,
        "rank_ham": rank_ham, "deficit": deficit,
        "dim_cut_space": cut.rank, "dim_cnperp": dim_cnperp,
        "quotient_dim": dim_cnperp - cut.rank,
        "max_weight_searched": max_weight,
        "n_low_weight_vectors": n_low,
        "weight_histogram": {str(k): v for k, v in weight_hist.items()},
        "minimal_witnesses": minimal,
        "maximal_representatives": maximal,
        "unknown_witnesses": unknown,
        "elapsed_s": round(time.time() - t0, 1),
    }
    DATA.mkdir(parents=True, exist_ok=True)
    out = DATA / f"witness_map_n{n}.json"
    out.write_text(json.dumps(result, indent=2))
    print(f"\n-> {out}  ({result['elapsed_s']}s)")
    return result


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=6)
    ap.add_argument("--tours", type=int, default=None,
                    help="K amostras; omitido = exaustivo (so' viavel em n=6)")
    ap.add_argument("--max-weight", type=int, default=6)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    analyse(a.n, a.tours, a.max_weight, a.seed)


if __name__ == "__main__":
    main()
