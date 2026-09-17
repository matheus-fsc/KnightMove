#!/usr/bin/env python3
"""
gt_paths_6x6.py
===============
Enumera por backtracking exaustivo todos os caminhos hamiltonianos
start→end no tabuleiro 6×6 (não é ciclo — endpoints fixos).

Necessário para KL-divergência porque o sampler Z3 produz caminhos
start→end (via aresta virtual), não ciclos fechados. O catálogo
destruction_catalogue.json só tem ciclos.

Saída:
  benchmark/results/gt_paths_6x6_<start>_<end>.json
  - assinaturas (set de arestas) e contagens
"""

import argparse
import json
import sys
import time
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from cavalo_loop_destruicao_6x6 import (  # noqa: E402
    ADJ, EDGES_LIST, TOTAL, BOARD, vid, label
)


def enumerate_paths(start, end, *, progress_every=5_000_000):
    """
    Backtracking exaustivo de todos os caminhos hamiltonianos start→end.
    Não fecha em ciclo — apenas exige terminar em `end` cobrindo todos
    os vértices.

    Retorna lista de tuplas (sequência de vértices).
    """
    paths = []
    path = [start]
    alive = bytearray([1] * TOTAL)
    alive[start] = 0
    deg = [sum(1 for u in ADJ[v] if alive[u]) for v in range(TOTAL)]

    nodes_visited = [0]
    dead_ends = [0]
    t0 = time.perf_counter()

    def backtrack():
        nodes_visited[0] += 1
        if nodes_visited[0] % progress_every == 0:
            elapsed = time.perf_counter() - t0
            print(f"    [progresso] nós={nodes_visited[0]:>11,}  "
                  f"paths={len(paths):>7,}  "
                  f"dead={dead_ends[0]:>10,}  "
                  f"taxa={nodes_visited[0]/elapsed:>9,.0f}/s",
                  flush=True)

        n = len(path)
        current = path[-1]

        if n == TOTAL:
            if current == end:
                paths.append(tuple(path))
            return

        candidates = [u for u in ADJ[current] if alive[u]]
        if not candidates:
            dead_ends[0] += 1
            return

        # Warnsdorff: ordena por grau ascendente
        candidates.sort(key=lambda u: deg[u])

        # Poda básica
        if deg[candidates[0]] == 0 and n < TOTAL - 1:
            # vértice com grau 0 (isolado no subgrafo vivo) e não estamos
            # no penúltimo passo: dead-end garantido
            dead_ends[0] += 1
            return

        for nxt in candidates:
            alive[nxt] = 0
            for w in ADJ[nxt]:
                deg[w] -= 1
            path.append(nxt)
            backtrack()
            path.pop()
            alive[nxt] = 1
            for w in ADJ[nxt]:
                deg[w] += 1

    backtrack()
    elapsed = time.perf_counter() - t0
    return paths, {
        "elapsed_s": round(elapsed, 4),
        "nodes_visited": nodes_visited[0],
        "dead_ends": dead_ends[0],
    }


def path_to_signature(path_seq, edge_idx):
    """Vetor booleano sobre EDGES_LIST das arestas do caminho."""
    sig = [False] * len(edge_idx)
    for i in range(len(path_seq) - 1):
        u, v = path_seq[i], path_seq[i + 1]
        e = (min(u, v), max(u, v))
        if e in edge_idx:
            sig[edge_idx[e]] = True
    return tuple(sig)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--start", type=int, nargs=2, required=True,
                   help="row col (ex: 0 0)")
    p.add_argument("--end", type=int, nargs=2, required=True,
                   help="row col (ex: 0 5)")
    p.add_argument("--out", default=None)
    args = p.parse_args()

    start = vid(*args.start)
    end = vid(*args.end)
    print(f"Enumerando caminhos hamiltonianos {label(start)} → {label(end)}")
    print(f"  start={tuple(args.start)} (v={start})  end={tuple(args.end)} (v={end})")

    paths, stats = enumerate_paths(start, end)

    edge_idx = {e: i for i, e in enumerate(EDGES_LIST)}
    sig_counter = Counter()
    for pth in paths:
        sig_counter[path_to_signature(pth, edge_idx)] += 1

    print(f"  Total de caminhos     : {len(paths):,}")
    print(f"  Assinaturas distintas : {len(sig_counter):,}")
    print(f"  Tempo                 : {stats['elapsed_s']}s")
    print(f"  Nós visitados         : {stats['nodes_visited']:,}")
    print(f"  Dead-ends             : {stats['dead_ends']:,}")

    if args.out is None:
        out_path = (ROOT / "benchmark" / "results"
                    / f"gt_paths_6x6_{label(start)}_{label(end)}.json")
    else:
        out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # serializa: assinaturas como listas, contagens como inteiros
    payload = {
        "board": BOARD,
        "start": list(args.start),
        "end": list(args.end),
        "start_label": label(start),
        "end_label": label(end),
        "n_paths_total": len(paths),
        "n_unique_signatures": len(sig_counter),
        "stats": stats,
        "signature_counts": [
            {"sig": list(int(b) for b in sig), "count": cnt}
            for sig, cnt in sig_counter.items()
        ],
    }
    with open(out_path, "w") as f:
        json.dump(payload, f)
    print(f"Salvo: {out_path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
