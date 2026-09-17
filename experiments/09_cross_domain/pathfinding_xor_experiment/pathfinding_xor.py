#!/usr/bin/env python3
"""
Pathfinding: GF(2) cycle-space XOR vs repeated A*
=================================================

Hipotese: em grafos PLANARES ESPARSOS (mapas em grade com obstaculos), os
ciclos fundamentais sao "corredores" geometricamente locais, e o XOR (soma em
GF(2)) do caminho-base com um ciclo fundamental quase sempre produz um caminho
alternativo s->t VALIDO. Isso contrasta com o TSP (grafo completo denso), onde
so ~3% dos ciclos fundamentais preservavam a hamiltonicidade.

Baseline: A* repetido com penalizacao de arestas (k-shortest aproximado).

Comparacao honesta de:
  - taxa de compatibilidade (XOR simples gera caminho valido?)
  - tempo para achar N=10 caminhos
  - diversidade (distancia de Jaccard media par-a-par)
  - modo de falha dominante (grau vs desconexao)

Uso:
    ../venv/bin/python pathfinding_xor.py
"""

import csv
import json
import os
import time
from collections import defaultdict

import numpy as np
import networkx as nx

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(HERE, "results")
JSON_PATH = os.path.join(RESULTS_DIR, "pathfinding_xor_experiment.json")
CSV_PATH = os.path.join(RESULTS_DIR, "pathfinding_xor_per_instance.csv")

SIZES = [20, 30, 50]
N_SEEDS = 20
OBSTACLE_DENSITY = 0.30
N_PATHS = 10          # quantos caminhos alternativos buscar
MAX_PAIRS = 500       # limite de pares de ciclos para XOR duplo
LOCAL_LEN = 10        # ciclo "local" = ate 10 arestas
TSP_BASELINE = 3.0    # % de compatibilidade observado no experimento TSP


# --------------------------------------------------------------------------
# Arestas como frozenset({u, v}); caminho/ciclo como conjunto de arestas.
# --------------------------------------------------------------------------
def E(u, v):
    return frozenset((u, v))


def seq_to_edges(seq):
    return {E(seq[k], seq[k + 1]) for k in range(len(seq) - 1)}


def manhattan(a, b):
    return abs(a[0] - b[0]) + abs(a[1] - b[1])


# --------------------------------------------------------------------------
# STEP 1 - geracao de mapa e grafo
# --------------------------------------------------------------------------
def gen_graph(n, seed):
    """Gera grade n x n com 30% de celulas bloqueadas, garantindo caminho
    s->t. Se nao houver caminho, tenta seed+1, seed+2, ... Retorna
    (Gc, S, T, seed_usado), onde Gc e a componente conexa que contem S e T."""
    S = (0, 0)
    T = (n - 1, n - 1)
    attempt = 0
    while True:
        cur_seed = seed + attempt
        rng = np.random.default_rng(cur_seed)
        blocked = rng.random((n, n)) < OBSTACLE_DENSITY
        blocked[S] = False
        blocked[T] = False

        G = nx.grid_2d_graph(n, n)
        G.remove_nodes_from(
            [(r, c) for r in range(n) for c in range(n) if blocked[r, c]]
        )
        if T in G and nx.has_path(G, S, T):
            comp = nx.node_connected_component(G, S)
            Gc = G.subgraph(comp).copy()
            return Gc, S, T, cur_seed
        attempt += 1


# --------------------------------------------------------------------------
# STEP 2 - caminho base via A*
# --------------------------------------------------------------------------
def astar_path(G, S, T):
    return nx.astar_path(G, S, T, heuristic=manhattan)


# --------------------------------------------------------------------------
# STEP 3 - ciclos fundamentais via arvore DFS + back edges
# --------------------------------------------------------------------------
def fundamental_cycles(G, S):
    """DFS a partir de S. Cada back edge (u,v) define um ciclo fundamental:
    caminho-na-arvore(u->v) + aresta(u,v). Retorna lista de frozenset de
    arestas (frozenset of frozenset)."""
    Tdfs = nx.dfs_tree(G, source=S)
    parent = {}
    for p, c in Tdfs.edges():        # arestas dirigidas pai->filho
        parent[c] = p

    # profundidade via BFS na arvore (a partir da raiz S)
    depth = {S: 0}
    stack = [S]
    children = defaultdict(list)
    for c, p in parent.items():
        children[p].append(c)
    while stack:
        u = stack.pop()
        for w in children[u]:
            depth[w] = depth[u] + 1
            stack.append(w)

    tree_edges = {E(p, c) for c, p in parent.items()}

    def tree_path_edges(u, v):
        a, b = u, v
        edges = []
        while depth[a] > depth[b]:
            edges.append(E(a, parent[a])); a = parent[a]
        while depth[b] > depth[a]:
            edges.append(E(b, parent[b])); b = parent[b]
        while a != b:
            edges.append(E(a, parent[a])); a = parent[a]
            edges.append(E(b, parent[b])); b = parent[b]
        return edges

    cycles = []
    for u, v in G.edges():
        e = E(u, v)
        if e in tree_edges:
            continue
        cyc = set(tree_path_edges(u, v))
        cyc.add(e)
        cycles.append(frozenset(cyc))
    return cycles


# --------------------------------------------------------------------------
# STEP 4 - validade de um conjunto de arestas como caminho s->t
# --------------------------------------------------------------------------
def check_path(edges, S, T):
    """Retorna (valido, motivo). motivo in {None,'degree','disconnected',
    'both'}. Caminho s->t valido = S,T grau 1; demais grau 2; conexo."""
    deg = defaultdict(int)
    adj = defaultdict(list)
    for ed in edges:
        a, b = tuple(ed)
        deg[a] += 1
        deg[b] += 1
        adj[a].append(b)
        adj[b].append(a)

    degree_ok = True
    if deg.get(S, 0) != 1 or deg.get(T, 0) != 1:
        degree_ok = False
    else:
        for v, d in deg.items():
            want = 1 if (v == S or v == T) else 2
            if d != want:
                degree_ok = False
                break

    conn_ok = True
    if deg:
        start = S if S in deg else next(iter(deg))
        seen = {start}
        st = [start]
        while st:
            u = st.pop()
            for w in adj[u]:
                if w not in seen:
                    seen.add(w)
                    st.append(w)
        conn_ok = (len(seen) == len(deg))
    else:
        conn_ok = False

    valid = degree_ok and conn_ok
    if valid:
        reason = None
    elif (not degree_ok) and (not conn_ok):
        reason = "both"
    elif not degree_ok:
        reason = "degree"
    else:
        reason = "disconnected"
    return valid, reason


def xor(a, b):
    return a ^ b


# --------------------------------------------------------------------------
# STEP 5 - baseline: A* repetido com penalizacao de arestas
# --------------------------------------------------------------------------
def repeated_astar(G, S, T, n_paths=N_PATHS):
    """k-shortest aproximado: acha P1; remove a aresta do meio de cada caminho
    e re-roda A*. Retorna (lista de edge-sets, time_ms, success)."""
    t0 = time.perf_counter()
    Gw = G.copy()
    try:
        p = astar_path(Gw, S, T)
    except nx.NetworkXNoPath:
        return [], (time.perf_counter() - t0) * 1000.0, False

    vseqs = [p]
    edge_sets = [seq_to_edges(p)]
    seen = {frozenset(edge_sets[0])}

    while len(edge_sets) < n_paths:
        cur = vseqs[-1]
        m = len(cur) - 1                     # numero de arestas
        if m <= 0:
            break
        # ordem meio-para-fora dos indices de aresta
        mid = m // 2
        order = sorted(range(m), key=lambda k: abs(k - mid))
        removed = False
        for k in order:
            u, v = cur[k], cur[k + 1]
            if Gw.has_edge(u, v):
                Gw.remove_edge(u, v)
                removed = True
                break
        if not removed:
            break
        try:
            p = astar_path(Gw, S, T)
        except nx.NetworkXNoPath:
            break
        es = seq_to_edges(p)
        fz = frozenset(es)
        if fz in seen:
            vseqs.append(p)      # avanca a frente de remocao mesmo se repetido
            continue
        seen.add(fz)
        edge_sets.append(es)
        vseqs.append(p)

    t_ms = (time.perf_counter() - t0) * 1000.0
    success = len(edge_sets) >= n_paths
    return edge_sets, t_ms, success


# --------------------------------------------------------------------------
# STEP 6 - diversidade (distancia de Jaccard media par-a-par)
# --------------------------------------------------------------------------
def diversity(paths):
    if len(paths) < 2:
        return 0.0
    ds = []
    for i in range(len(paths)):
        for j in range(i + 1, len(paths)):
            inter = len(paths[i] & paths[j])
            union = len(paths[i] | paths[j])
            ds.append(1.0 - inter / union if union else 0.0)
    return float(np.mean(ds))


# --------------------------------------------------------------------------
# Uma instancia completa
# --------------------------------------------------------------------------
def run_instance(n, seed):
    G, S, T, used_seed = gen_graph(n, seed)
    V = G.number_of_nodes()
    Ec = G.number_of_edges()
    n_fund = Ec - V + 1

    # STEP 2 - caminho base
    t0 = time.perf_counter()
    base_seq = astar_path(G, S, T)
    base_astar_ms = (time.perf_counter() - t0) * 1000.0
    base = seq_to_edges(base_seq)
    base_len = len(base)

    # STEP 3 + 4 - ciclos fundamentais + XOR (tudo cronometrado junto)
    t0 = time.perf_counter()
    cycles = fundamental_cycles(G, S)
    cyc_lens = np.array([len(c) for c in cycles]) if cycles else np.array([0])
    n_local = int(np.sum(cyc_lens <= LOCAL_LEN)) if cycles else 0
    frac_local = n_local / len(cycles) if cycles else 0.0

    # XOR simples
    valid_single = []
    single_valid = 0
    fail_counts = {"degree": 0, "disconnected": 0, "both": 0}
    for c in cycles:
        cand = xor(base, set(c))
        ok, reason = check_path(cand, S, T)
        if ok:
            single_valid += 1
            valid_single.append(cand)
        else:
            fail_counts[reason] += 1
    single_total = len(cycles)
    single_rate = single_valid / single_total if single_total else 0.0

    # XOR duplo (amostra ate MAX_PAIRS pares)
    pair_valid = 0
    pair_total = 0
    valid_pairs = []
    k = len(cycles)
    if k >= 2:
        rng = np.random.default_rng(used_seed)
        all_pairs = k * (k - 1) // 2
        if all_pairs <= MAX_PAIRS:
            pair_iter = [(i, j) for i in range(k) for j in range(i + 1, k)]
        else:
            seen_p = set()
            pair_iter = []
            while len(pair_iter) < MAX_PAIRS:
                i = int(rng.integers(0, k))
                j = int(rng.integers(0, k))
                if i == j:
                    continue
                a, b = (i, j) if i < j else (j, i)
                if (a, b) in seen_p:
                    continue
                seen_p.add((a, b))
                pair_iter.append((a, b))
        for (i, j) in pair_iter:
            cand = xor(xor(base, set(cycles[i])), set(cycles[j]))
            ok, _ = check_path(cand, S, T)
            pair_total += 1
            if ok:
                pair_valid += 1
                valid_pairs.append(cand)
    pair_rate = pair_valid / pair_total if pair_total else 0.0
    time_xor_ms = (time.perf_counter() - t0) * 1000.0

    # caminhos validos unicos via XOR (single + pairs)
    unique = {}
    for p in valid_single + valid_pairs:
        unique[frozenset(p)] = p
    valid_xor_total = len(unique)

    # STEP 6 - conjuntos para diversidade (incluem o caminho base como rota 1)
    xor_set = [base]
    for p in unique.values():
        if frozenset(p) != frozenset(base):
            xor_set.append(p)
        if len(xor_set) >= N_PATHS:
            break
    div_xor = diversity(xor_set)

    # STEP 5 - baseline A* repetido
    astar_paths, astar_ms, astar_success = repeated_astar(G, S, T, N_PATHS)
    div_astar = diversity(astar_paths)

    rec = {
        "grid_size": n,
        "seed": seed,
        "used_seed": used_seed,
        "V": V,
        "E": Ec,
        "fundamental_cycles_total": n_fund,
        "cycles_collected": len(cycles),
        "cycle_len_min": int(cyc_lens.min()),
        "cycle_len_max": int(cyc_lens.max()),
        "cycle_len_median": float(np.median(cyc_lens)),
        "cycle_len_mean": float(cyc_lens.mean()),
        "frac_local_cycles": frac_local,
        "base_path_length": base_len,
        "base_astar_ms": base_astar_ms,
        "single_xor_compatibility_rate": single_rate,
        "valid_paths_xor_single": single_valid,
        "pair_xor_compatibility_rate": pair_rate,
        "valid_paths_xor_pairs": pair_valid,
        "pairs_tested": pair_total,
        "valid_paths_xor_total": valid_xor_total,
        "time_xor_ms": time_xor_ms,
        "paths_found_astar_repeated": len(astar_paths),
        "time_astar_repeated_ms": astar_ms,
        "astar_success": astar_success,
        "diversity_xor": div_xor,
        "diversity_astar": div_astar,
        "fail_degree": fail_counts["degree"],
        "fail_disconnected": fail_counts["disconnected"],
        "fail_both": fail_counts["both"],
    }
    return rec


# --------------------------------------------------------------------------
# Persistencia incremental
# --------------------------------------------------------------------------
CSV_FIELDS = [
    "grid_size", "seed", "used_seed", "V", "E",
    "fundamental_cycles_total", "cycle_len_min", "cycle_len_max",
    "cycle_len_median", "cycle_len_mean", "frac_local_cycles",
    "base_path_length", "single_xor_compatibility_rate",
    "valid_paths_xor_single", "pair_xor_compatibility_rate",
    "valid_paths_xor_pairs", "valid_paths_xor_total", "time_xor_ms",
    "paths_found_astar_repeated", "time_astar_repeated_ms", "astar_success",
    "diversity_xor", "diversity_astar",
    "fail_degree", "fail_disconnected", "fail_both",
]


def save_all(records):
    os.makedirs(RESULTS_DIR, exist_ok=True)
    with open(JSON_PATH, "w") as f:
        json.dump(records, f, indent=2)
    with open(CSV_PATH, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=CSV_FIELDS, extrasaction="ignore")
        w.writeheader()
        for r in records:
            w.writerow(r)


# --------------------------------------------------------------------------
# Relatorio agregado
# --------------------------------------------------------------------------
def report(records):
    print("\n=== PATHFINDING: GF(2) XOR vs REPEATED A* ===")
    for n in SIZES:
        rs = [r for r in records if r["grid_size"] == n]
        if not rs:
            continue
        comp = np.array([r["single_xor_compatibility_rate"] for r in rs]) * 100
        vx = np.array([r["valid_paths_xor_total"] for r in rs])
        va = np.array([r["paths_found_astar_repeated"] for r in rs])
        tx = np.array([r["time_xor_ms"] for r in rs])
        ta = np.array([r["time_astar_repeated_ms"] for r in rs])
        dx = np.array([r["diversity_xor"] for r in rs])
        da = np.array([r["diversity_astar"] for r in rs])
        speedup = ta.mean() / tx.mean() if tx.mean() > 0 else float("inf")
        print(f"\nGrid {n}x{n} ({len(rs)} instances):")
        print(f"  Compatibility rate (single XOR): {comp.mean():.1f}% "
              f"+/- {comp.std():.1f}%   [TSP baseline was {TSP_BASELINE:.0f}%]")
        print(f"  Valid paths found XOR: {vx.mean():.1f} avg "
              f"({vx.min()} min, {vx.max()} max)")
        print(f"  Valid paths found A*:  {va.mean():.1f} avg")
        print(f"  Time XOR: {tx.mean():.1f}ms avg | "
              f"Time A* repeated: {ta.mean():.1f}ms avg")
        print(f"  Speedup XOR vs A*: {speedup:.2f}x")
        print(f"  Diversity XOR: {dx.mean():.3f} | Diversity A*: {da.mean():.3f}")

    # agregados globais
    comp_all = np.array([r["single_xor_compatibility_rate"] for r in records]) * 100
    tx_all = np.array([r["time_xor_ms"] for r in records])
    ta_all = np.array([r["time_astar_repeated_ms"] for r in records])
    dx_all = np.array([r["diversity_xor"] for r in records])
    da_all = np.array([r["diversity_astar"] for r in records])

    # correlacao densidade-real de obstaculos vs taxa de compatibilidade.
    # densidade efetiva = 1 - V/(n^2)
    dens = np.array([1.0 - r["V"] / (r["grid_size"] ** 2) for r in records])
    comp_frac = np.array([r["single_xor_compatibility_rate"] for r in records])
    corr_dens = float(np.corrcoef(dens, comp_frac)[0, 1]) if len(records) > 1 else 0.0
    mean_len = np.array([r["cycle_len_mean"] for r in records])
    corr_len = float(np.corrcoef(mean_len, comp_frac)[0, 1]) if len(records) > 1 else 0.0

    faster = ta_all.mean() > tx_all.mean()
    more_div = dx_all.mean() > da_all.mean()
    hyp_conf = comp_all.mean() > TSP_BASELINE

    # modo de falha agregado
    fd = sum(r["fail_degree"] for r in records)
    fdis = sum(r["fail_disconnected"] for r in records)
    fb = sum(r["fail_both"] for r in records)
    ftot = fd + fdis + fb

    print("\nKEY FINDINGS:")
    print(f"  Hypothesis (compat >> 3%): "
          f"{'CONFIRMED' if hyp_conf else 'REFUTED'} "
          f"(mean compat = {comp_all.mean():.1f}%)")
    print(f"  XOR faster than repeated A* for N={N_PATHS} paths: "
          f"{'YES' if faster else 'NO'} "
          f"(XOR {tx_all.mean():.1f}ms vs A* {ta_all.mean():.1f}ms)")
    print(f"  XOR produces more diverse paths than repeated A*: "
          f"{'YES' if more_div else 'NO'} "
          f"(XOR {dx_all.mean():.3f} vs A* {da_all.mean():.3f})")
    print(f"  Correlation obstacle_density vs compatibility_rate: {corr_dens:+.3f}")
    print(f"  Correlation cycle_mean_length vs compatibility_rate: {corr_len:+.3f}")

    print("\nFAILURE MODE (single-XOR invalid candidates):")
    if ftot:
        print(f"  degree only:    {100*fd/ftot:5.1f}%")
        print(f"  disconnected:   {100*fdis/ftot:5.1f}%")
        print(f"  both:           {100*fb/ftot:5.1f}%")
        print(f"  (total invalid candidates: {ftot})")
    else:
        print("  (nenhum candidato invalido)")

    print("\n--- KEY QUESTIONS ---")
    print(f"1. Single-XOR compatibility rate: {comp_all.mean():.1f}% "
          f"(TSP era ~{TSP_BASELINE:.0f}%) -> "
          f"{'>> 3%, viavel' if hyp_conf else 'NAO >> 3%'}")
    print(f"2. XOR mais rapido que A* repetido p/ {N_PATHS} caminhos: "
          f"{'SIM' if faster else 'NAO'}")
    print(f"3. XOR mais diverso que penalizacao A*: {'SIM' if more_div else 'NAO'}")
    print(f"4. Correlacao densidade={corr_dens:+.3f}, comprimento_ciclo={corr_len:+.3f}")
    dominant = max([("degree", fd), ("disconnected", fdis), ("both", fb)],
                   key=lambda x: x[1])[0] if ftot else "n/a"
    print(f"5. Modo de falha dominante: {dominant}")

    print(f"\nResultados salvos em:\n  {JSON_PATH}\n  {CSV_PATH}")


def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)
    records = []
    total = len(SIZES) * N_SEEDS
    done = 0
    for n in SIZES:
        for s in range(N_SEEDS):
            rec = run_instance(n, s)
            records.append(rec)
            done += 1
            print(f"  [{done:2d}/{total}] n={n} seed={s} "
                  f"(used {rec['used_seed']}): V={rec['V']} E={rec['E']} "
                  f"cyc={rec['cycles_collected']} "
                  f"compat={rec['single_xor_compatibility_rate']*100:5.1f}% "
                  f"xorPaths={rec['valid_paths_xor_total']:3d} "
                  f"A*={rec['paths_found_astar_repeated']} "
                  f"tX={rec['time_xor_ms']:.0f}ms tA={rec['time_astar_repeated_ms']:.0f}ms")
        save_all(records)   # salvamento incremental por tamanho de grade

    save_all(records)
    report(records)


if __name__ == "__main__":
    main()
