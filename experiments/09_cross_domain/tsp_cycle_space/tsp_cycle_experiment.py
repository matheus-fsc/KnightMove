#!/usr/bin/env python3
"""
TSP Cycle-Space vs 2-opt
========================

Experimento: usar a estrutura do espaco de ciclos GF(2) (H_1(G; F_2)) para
guiar busca local no TSP Euclidiano, comparando com 2-opt cego.

Tour Hamiltoniano = vetor em F_2^|E| (conjunto de arestas, i<j).
Ciclo fundamental C_e (aresta nao-arvore da MST) = outro vetor em F_2^|E|.
Movimento = diferenca simetrica (XOR / soma em GF(2)) tour XOR C_e.
Aplica-se se o resultado e um ciclo Hamiltoniano valido E tem custo menor.

Pergunta-chave: ordenar os ciclos por delta estimado (estrategia d) torna
a busca guiada por ciclos competitiva com 2-opt? Relatorio honesto.

Uso:
    ./venv/bin/python tsp_cycle_space/tsp_cycle_experiment.py
"""

import json
import os
import time
import itertools
from collections import defaultdict

import numpy as np
import networkx as nx

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(HERE, "results")
RESULTS_PATH = os.path.join(RESULTS_DIR, "tsp_cycle_experiment.json")

MASTER_SEED = 42
SIZES = [15, 20, 25]
N_SEEDS = 10
STRATEGIES = ["random", "length", "min_edge", "delta"]


# --------------------------------------------------------------------------
# STEP 1 - geracao de instancia
# --------------------------------------------------------------------------
def gen_instance(n, seed):
    """Pontos uniformes em [0,100]^2 -> matriz de distancias Euclidianas."""
    rng = np.random.default_rng(seed)
    pts = rng.uniform(0.0, 100.0, size=(n, 2))
    diff = pts[:, None, :] - pts[None, :, :]
    D = np.sqrt((diff ** 2).sum(axis=2))
    return pts, D


def edge(i, j):
    """Aresta canonica (i<j)."""
    return (i, j) if i < j else (j, i)


def tour_edges(seq):
    """Sequencia de vertices (ciclo) -> conjunto de arestas canonicas."""
    n = len(seq)
    return {edge(seq[k], seq[(k + 1) % n]) for k in range(n)}


def edge_set_cost(edges, D):
    return sum(D[i, j] for (i, j) in edges)


def seq_cost(seq, D):
    n = len(seq)
    return sum(D[seq[k], seq[(k + 1) % n]] for k in range(n))


# --------------------------------------------------------------------------
# STEP 2 - tour inicial (nearest neighbor a partir do vertice 0)
# --------------------------------------------------------------------------
def nearest_neighbor_tour(D, start=0):
    n = D.shape[0]
    visited = [False] * n
    seq = [start]
    visited[start] = True
    cur = start
    for _ in range(n - 1):
        best, best_d = -1, np.inf
        for j in range(n):
            if not visited[j] and D[cur, j] < best_d:
                best, best_d = j, D[cur, j]
        seq.append(best)
        visited[best] = True
        cur = best
    return seq


# --------------------------------------------------------------------------
# STEP 3 - base do espaco de ciclos (ciclos fundamentais via MST)
# --------------------------------------------------------------------------
def fundamental_cycles(n, D):
    """
    Retorna lista de ciclos fundamentais como conjuntos (frozenset) de arestas.
    Um por aresta nao-arvore: C_e = caminho(u->v na MST) + (u,v).
    Total = |E| - |V| + 1.
    """
    G = nx.Graph()
    for i in range(n):
        for j in range(i + 1, n):
            G.add_edge(i, j, weight=float(D[i, j]))

    T = nx.minimum_spanning_tree(G, weight="weight")
    tree_edges = {edge(u, v) for u, v in T.edges()}

    cycles = []
    for i in range(n):
        for j in range(i + 1, n):
            e = (i, j)
            if e in tree_edges:
                continue
            # caminho i->j na arvore (unico)
            path = nx.shortest_path(T, source=i, target=j)
            cyc = {edge(path[k], path[k + 1]) for k in range(len(path) - 1)}
            cyc.add(e)
            cycles.append(frozenset(cyc))
    return cycles


# --------------------------------------------------------------------------
# STEP 4 - validade Hamiltoniana e XOR
# --------------------------------------------------------------------------
def is_hamiltonian(edges, n):
    """Grau 2 em todo vertice E uma unica componente conexa (n arestas)."""
    if len(edges) != n:
        return False
    deg = defaultdict(int)
    adj = defaultdict(list)
    for (i, j) in edges:
        deg[i] += 1
        deg[j] += 1
        adj[i].append(j)
        adj[j].append(i)
    if len(deg) != n or any(d != 2 for d in deg.values()):
        return False
    # conectividade: caminhada a partir de um vertice qualquer
    start = next(iter(deg))
    seen = {start}
    stack = [start]
    while stack:
        u = stack.pop()
        for w in adj[u]:
            if w not in seen:
                seen.add(w)
                stack.append(w)
    return len(seen) == n


def xor_delta(tour, cyc, D):
    """delta = peso(arestas adicionadas) - peso(arestas removidas).
    adicionadas = cyc \\ tour ; removidas = cyc & tour."""
    delta = 0.0
    for e in cyc:
        if e in tour:
            delta -= D[e[0], e[1]]
        else:
            delta += D[e[0], e[1]]
    return delta


def apply_xor(tour, cyc):
    """Diferenca simetrica (novo conjunto de arestas)."""
    return tour ^ set(cyc)


# --------------------------------------------------------------------------
# STEP 4 - busca guiada por ciclos (com ordenacao)
# --------------------------------------------------------------------------
def cycle_guided_search(seq0, cycles, D, strategy, rng):
    """
    Retorna (final_cost, iterations, time_ms).
    iterations = numero de movimentos aplicados (melhorias).
    """
    n = len(seq0)
    t0 = time.perf_counter()
    tour = tour_edges(seq0)
    cost = edge_set_cost(tour, D)
    iters = 0

    # ordenacoes estaticas (a,b,c) sao pre-computadas; (d) e dinamica
    def static_order():
        if strategy == "random":
            order = list(cycles)
            rng.shuffle(order)
            return order
        if strategy == "length":
            return sorted(cycles, key=len)
        if strategy == "min_edge":
            return sorted(cycles, key=lambda c: min(D[i, j] for (i, j) in c))
        return None  # delta -> dinamico

    improved = True
    while improved:
        improved = False
        if strategy == "delta":
            # escolhe o ciclo de menor delta estimado que produza tour valido
            scored = sorted(
                ((xor_delta(tour, c, D), c) for c in cycles),
                key=lambda x: x[0],
            )
            for delta, c in scored:
                if delta >= -1e-12:
                    break  # nenhum movimento melhora
                cand = apply_xor(tour, c)
                if is_hamiltonian(cand, n):
                    tour = cand
                    cost += delta
                    iters += 1
                    improved = True
                    break
        else:
            order = static_order()
            for c in order:
                delta = xor_delta(tour, c, D)
                if delta < -1e-12:
                    cand = apply_xor(tour, c)
                    if is_hamiltonian(cand, n):
                        tour = cand
                        cost += delta
                        iters += 1
                        improved = True
                        break  # reinicia varredura apos aplicar
    t_ms = (time.perf_counter() - t0) * 1000.0
    return cost, iters, t_ms


# --------------------------------------------------------------------------
# STEP 5 - baseline 2-opt
# --------------------------------------------------------------------------
def two_opt(seq0, D):
    """2-opt classico. Retorna (final_cost, iterations, time_ms)."""
    t0 = time.perf_counter()
    seq = list(seq0)
    n = len(seq)
    iters = 0
    improved = True
    while improved:
        improved = False
        for i in range(n - 1):
            a, b = seq[i], seq[i + 1]
            for k in range(i + 2, n):
                c = seq[k]
                d = seq[(k + 1) % n]
                if i == 0 and (k + 1) % n == 0:
                    continue  # mesma aresta
                # ganho ao reverter o segmento seq[i+1..k]
                delta = (D[a, c] + D[b, d]) - (D[a, b] + D[c, d])
                if delta < -1e-12:
                    seq[i + 1:k + 1] = seq[i + 1:k + 1][::-1]
                    iters += 1
                    improved = True
                    a, b = seq[i], seq[i + 1]
    cost = seq_cost(seq, D)
    t_ms = (time.perf_counter() - t0) * 1000.0
    return cost, iters, t_ms


# --------------------------------------------------------------------------
# STEP 6 - benchmark
# --------------------------------------------------------------------------
def save_partial(records):
    os.makedirs(RESULTS_DIR, exist_ok=True)
    with open(RESULTS_PATH, "w") as f:
        json.dump(records, f, indent=2)


def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)
    master_rng = np.random.default_rng(MASTER_SEED)

    records = []
    done = 0
    combos = list(itertools.product(SIZES, range(N_SEEDS)))

    for (n, s) in combos:
        # seed derivada do mestre 42 -> reprodutivel e distinta por instancia
        inst_seed = int(master_rng.integers(0, 2**31 - 1))
        pts, D = gen_instance(n, inst_seed)
        seq0 = nearest_neighbor_tour(D, start=0)
        init_cost = seq_cost(seq0, D)
        cycles = fundamental_cycles(n, D)

        # busca por ciclos em todas as estrategias (mesmo tour inicial)
        strat_results = {}
        for strat in STRATEGIES:
            strat_rng = np.random.default_rng(inst_seed ^ (hash(strat) & 0xFFFF))
            c_cost, c_it, c_ms = cycle_guided_search(seq0, cycles, D, strat, strat_rng)
            strat_results[strat] = {
                "final_cost": c_cost,
                "iterations": c_it,
                "time_ms": c_ms,
            }

        # melhor estrategia de ciclos nesta instancia
        best_strat = min(strat_results, key=lambda k: strat_results[k]["final_cost"])
        best = strat_results[best_strat]

        # baseline 2-opt do mesmo tour inicial
        o_cost, o_it, o_ms = two_opt(seq0, D)

        gap = (best["final_cost"] - o_cost) / o_cost * 100.0

        rec = {
            "n": n,
            "seed_index": s,
            "instance_seed": inst_seed,
            "n_cycles": len(cycles),
            "initial_cost": init_cost,
            "cycle_strategies": strat_results,
            "best_cycle_strategy": best_strat,
            "final_cost_cycle_method": best["final_cost"],
            "iterations_cycle": best["iterations"],
            "time_cycle_ms": best["time_ms"],
            "final_cost_2opt": o_cost,
            "iterations_2opt": o_it,
            "time_2opt_ms": o_ms,
            "gap_to_2opt": gap,
        }
        records.append(rec)
        done += 1

        print(f"  [{done:2d}/{len(combos)}] n={n} seed={s}: "
              f"init={init_cost:7.1f}  cycle({best_strat})={best['final_cost']:7.1f}  "
              f"2opt={o_cost:7.1f}  gap={gap:+6.2f}%")

        if done % 5 == 0:
            save_partial(records)

    save_partial(records)
    report(records)


# --------------------------------------------------------------------------
# Relatorio agregado
# --------------------------------------------------------------------------
def report(records):
    print("\n=== TSP CYCLE-SPACE vs 2-OPT ===")
    for n in SIZES:
        rs = [r for r in records if r["n"] == n]
        cyc_avg = np.mean([r["final_cost_cycle_method"] for r in rs])
        opt_avg = np.mean([r["final_cost_2opt"] for r in rs])
        gap = np.mean([r["gap_to_2opt"] for r in rs])
        tc = np.mean([r["time_cycle_ms"] for r in rs])
        to = np.mean([r["time_2opt_ms"] for r in rs])
        tratio = tc / to if to > 0 else float("inf")
        print(f"n={n}: cycle_avg={cyc_avg:7.2f}  2opt_avg={opt_avg:7.2f}  "
              f"gap={gap:+6.2f}%  time_ratio={tratio:6.2f}x")

    # melhor estrategia global (qual venceu como best_strategy com mais frequencia)
    strat_wins = defaultdict(int)
    for r in records:
        strat_wins[r["best_cycle_strategy"]] += 1
    # tambem: custo medio por estrategia (sobre todas as instancias)
    strat_avg = {}
    for strat in STRATEGIES:
        strat_avg[strat] = np.mean(
            [r["cycle_strategies"][strat]["final_cost"] for r in records]
        )
    best_overall = min(strat_avg, key=strat_avg.get)

    beats = sum(1 for r in records if r["final_cost_cycle_method"] < r["final_cost_2opt"])

    print(f"\nBest ordering strategy: {best_overall} "
          f"(custo medio={strat_avg[best_overall]:.2f})")
    print("  custo medio por estrategia:",
          {k: round(v, 2) for k, v in strat_avg.items()})
    print("  vezes escolhida como melhor:", dict(strat_wins))
    print(f"Cases where cycle-method beats 2-opt: {beats} / {len(records)}")

    # resposta a pergunta-chave
    delta_avg = strat_avg["delta"]
    overall_gap = np.mean([r["gap_to_2opt"] for r in records])
    print("\n--- PERGUNTA-CHAVE ---")
    print(f"Estrategia 'delta' custo medio = {delta_avg:.2f}")
    print(f"Gap medio (melhor ciclo vs 2-opt) = {overall_gap:+.2f}%")
    if overall_gap <= 0:
        print("=> busca guiada por ciclos EMPATA/VENCE o 2-opt em media.")
    else:
        print(f"=> 2-opt VENCE em media por {overall_gap:.2f}%. "
              "Honestamente, o 2-opt e mais forte.")
    print(f"\nResultados completos salvos em: {RESULTS_PATH}")


if __name__ == "__main__":
    main()
