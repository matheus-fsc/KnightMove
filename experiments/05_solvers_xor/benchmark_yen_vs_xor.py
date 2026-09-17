#!/usr/bin/env python3
"""
Benchmark: Yen's k-shortest paths  vs  GF(2) Cycle-Space XOR Pathfinding
=========================================================================

Cenários de avaliação:
  A) Labirinto Confinado (20×20, árvore geradora → corredores estreitos)
  B) Escalonamento de K  (20×20, 30% obstáculos, k ∈ {5, 20, 50, 100})
  C) Vértices Compulsórios (grafo planar, 2 via-points obrigatórios)

Métricas:
  - Tempo de Execução Total (ms)
  - Taxa de Validade (% caminhos XOR válidos s→t)
  - Diversidade Média (Jaccard distance par-a-par)
  - Custo de Extensão (comprimento médio rotas / comprimento mínimo)

Uso:
    python benchmark_yen_vs_xor.py
    # ou com venv:
    ../venv/bin/python benchmark_yen_vs_xor.py

Autor: Benchmark gerado para pesquisa em Teoria dos Grafos / XOR Pathfinding.
"""

import json
import os
import signal
import sys
import time
from collections import defaultdict
from itertools import combinations, islice

import matplotlib
matplotlib.use("Agg")  # backend sem GUI para salvar PNGs
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import networkx as nx

# ---------------------------------------------------------------------------
# Constantes globais
# ---------------------------------------------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(HERE, "benchmark_yen_vs_xor_results")

GRID_SIZE = 20
OBSTACLE_DENSITY = 0.30
K_VALUES = [5, 20, 50, 100]
N_SEEDS = 5                 # instâncias por cenário
YEN_TIMEOUT_S = 60           # timeout por chamada ao Yen (segundos)
MAX_XOR_PAIRS = 2000         # limite de pares para XOR duplo

# ---------------------------------------------------------------------------
# Utilidades de aresta
# ---------------------------------------------------------------------------
def E(u, v):
    """Aresta como frozenset para grafos não-dirigidos."""
    return frozenset((u, v))


def seq_to_edges(seq):
    """Converte sequência de nós em conjunto de arestas frozenset."""
    return {E(seq[i], seq[i + 1]) for i in range(len(seq) - 1)}


def edges_to_path(edges, s, t):
    """Reconstrói a sequência de nós de um caminho s→t a partir de arestas.
    Retorna None se não for um caminho simples válido."""
    if not edges:
        return None
    adj = defaultdict(list)
    for e in edges:
        a, b = tuple(e)
        adj[a].append(b)
        adj[b].append(a)
    # Verificar graus: s e t devem ter grau 1, todos os outros grau 2
    for v, nbrs in adj.items():
        expected = 1 if v in (s, t) else 2
        if len(nbrs) != expected:
            return None
    # Percorrer o caminho
    path = [s]
    visited = {s}
    current = s
    while current != t:
        found_next = False
        for nbr in adj[current]:
            if nbr not in visited:
                visited.add(nbr)
                path.append(nbr)
                current = nbr
                found_next = True
                break
        if not found_next:
            return None
    if len(path) != len(edges) + 1:
        return None
    return path


def manhattan(a, b):
    """Heurística Manhattan para A*."""
    return abs(a[0] - b[0]) + abs(a[1] - b[1])


# ---------------------------------------------------------------------------
# Timeout via signal (Unix)
# ---------------------------------------------------------------------------
class TimeoutError(Exception):
    pass


def _timeout_handler(signum, frame):
    raise TimeoutError("Yen timeout exceeded")


# ---------------------------------------------------------------------------
# Validação de caminho s→t
# ---------------------------------------------------------------------------
def check_path_valid(edges, s, t):
    """Verifica se um conjunto de arestas forma um caminho simples s→t.
    Retorna (válido: bool, motivo: str|None)."""
    if not edges:
        return False, "empty"
    deg = defaultdict(int)
    adj = defaultdict(set)
    for e in edges:
        a, b = tuple(e)
        deg[a] += 1
        deg[b] += 1
        adj[a].add(b)
        adj[b].add(a)

    # Graus corretos
    degree_ok = True
    if deg.get(s, 0) != 1 or deg.get(t, 0) != 1:
        degree_ok = False
    else:
        for v, d in deg.items():
            want = 1 if v in (s, t) else 2
            if d != want:
                degree_ok = False
                break

    # Conectividade
    if deg:
        start = s if s in deg else next(iter(deg))
        seen = {start}
        stack = [start]
        while stack:
            u = stack.pop()
            for w in adj[u]:
                if w not in seen:
                    seen.add(w)
                    stack.append(w)
        conn_ok = len(seen) == len(deg)
    else:
        conn_ok = False

    if degree_ok and conn_ok:
        return True, None
    elif not degree_ok and not conn_ok:
        return False, "both"
    elif not degree_ok:
        return False, "degree"
    else:
        return False, "disconnected"


def xor_edges(a, b):
    """Diferença simétrica de dois conjuntos de arestas (XOR em GF(2))."""
    return a ^ b


# ---------------------------------------------------------------------------
# Diversidade: distância de Jaccard média par-a-par
# ---------------------------------------------------------------------------
def jaccard_diversity(paths):
    """Calcula a distância de Jaccard média entre todos os pares de caminhos.
    Cada caminho é representado como um set de arestas (frozenset)."""
    if len(paths) < 2:
        return 0.0
    dists = []
    for i in range(len(paths)):
        for j in range(i + 1, len(paths)):
            inter = len(paths[i] & paths[j])
            union = len(paths[i] | paths[j])
            dists.append(1.0 - inter / union if union > 0 else 0.0)
    return float(np.mean(dists))


# ---------------------------------------------------------------------------
# Custo de extensão
# ---------------------------------------------------------------------------
def cost_extension(paths_edges, shortest_len):
    """Razão entre o comprimento médio das k rotas e o caminho mínimo."""
    if not paths_edges or shortest_len == 0:
        return 0.0
    avg_len = np.mean([len(p) for p in paths_edges])
    return float(avg_len / shortest_len)


# ============================================================================
# ALGORITMO 1: Yen (via NetworkX shortest_simple_paths)
# ============================================================================
def run_yen(G, s, t, k, timeout_s=YEN_TIMEOUT_S):
    """Executa Yen's k-shortest simple paths com timeout.
    Retorna (list[set de arestas], tempo_ms, n_encontrados)."""
    paths_edges = []
    paths_seqs = []
    t0 = time.perf_counter()

    # Configurar timeout (só funciona no thread principal em Unix)
    old_handler = None
    try:
        old_handler = signal.signal(signal.SIGALRM, _timeout_handler)
        signal.alarm(timeout_s)
    except (ValueError, OSError, AttributeError):
        pass  # Não está no thread principal ou não é Unix

    try:
        gen = nx.shortest_simple_paths(G, s, t)
        for path_seq in islice(gen, k):
            paths_seqs.append(path_seq)
            paths_edges.append(seq_to_edges(path_seq))
    except TimeoutError:
        pass  # Retorna o que conseguiu
    except nx.NetworkXNoPath:
        pass
    finally:
        try:
            signal.alarm(0)
            if old_handler is not None:
                signal.signal(signal.SIGALRM, old_handler)
        except (ValueError, OSError, AttributeError):
            pass

    t_ms = (time.perf_counter() - t0) * 1000.0
    return paths_edges, t_ms, len(paths_edges)


# ============================================================================
# ALGORITMO 2: XOR Cycle-Space Pathfinding
# ============================================================================
def fundamental_cycles_dfs(G, source):
    """Calcula ciclos fundamentais via DFS tree + back edges.
    Retorna lista de frozenset de arestas (cada ciclo)."""
    Tdfs = nx.dfs_tree(G, source=source)
    parent = {}
    for p, c in Tdfs.edges():
        parent[c] = p

    # Calcular profundidade
    depth = {source: 0}
    children = defaultdict(list)
    for c, p in parent.items():
        children[p].append(c)
    stack = [source]
    while stack:
        u = stack.pop()
        for w in children[u]:
            depth[w] = depth[u] + 1
            stack.append(w)

    tree_edges = {E(p, c) for c, p in parent.items()}

    def tree_path_edges(u, v):
        """Arestas no caminho u→v na DFS tree."""
        a, b = u, v
        edges = []
        while depth.get(a, 0) > depth.get(b, 0):
            edges.append(E(a, parent[a]))
            a = parent[a]
        while depth.get(b, 0) > depth.get(a, 0):
            edges.append(E(b, parent[b]))
            b = parent[b]
        while a != b:
            edges.append(E(a, parent[a]))
            a = parent[a]
            edges.append(E(b, parent[b]))
            b = parent[b]
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


def run_xor_pathfinding(G, s, t, k, cycles=None, required_nodes=None):
    """Gerador XOR: a partir de caminho base + ciclos fundamentais, gera
    alternativas via diferença simétrica.

    Args:
        G: grafo NetworkX
        s, t: origem e destino
        k: número de caminhos desejados
        cycles: ciclos pré-computados (opcional)
        required_nodes: set de nós obrigatórios (via-points) no Cenário C

    Retorna (all_paths_edges, valid_paths_edges, tempo_ms, stats).
    """
    t0 = time.perf_counter()

    # 1) Caminho base via A*
    try:
        base_seq = nx.astar_path(G, s, t, heuristic=manhattan)
    except nx.NetworkXNoPath:
        t_ms = (time.perf_counter() - t0) * 1000.0
        return [], [], t_ms, {"valid": 0, "total_candidates": 0}

    base_edges = seq_to_edges(base_seq)
    base_nodes = set(base_seq)

    # 2) Ciclos fundamentais
    if cycles is None:
        cycles = fundamental_cycles_dfs(G, s)

    # 3) Filtrar ciclos que intersectam o caminho base (overlap >= 1 aresta)
    overlapping = []
    for c in cycles:
        if len(c & base_edges) >= 1:
            overlapping.append(c)

    # No cenário C, filtrar ciclos para forçar inclusão dos via-points
    if required_nodes:
        def cycle_nodes(c):
            nodes = set()
            for e in c:
                nodes.update(e)
            return nodes

        # Priorizar ciclos que passem por via-points
        via_cycles = []
        non_via_cycles = []
        for c in overlapping:
            cn = cycle_nodes(c)
            if any(vp in cn for vp in required_nodes):
                via_cycles.append(c)
            else:
                non_via_cycles.append(c)
        # Colocar ciclos com via-points primeiro
        overlapping = via_cycles + non_via_cycles

    # 4) XOR simples: base ⊕ ciclo
    unique_paths = {}
    unique_paths[frozenset(base_edges)] = base_edges  # caminho base como rota 1

    total_candidates = 0
    valid_count = 0

    for c in overlapping:
        if len(unique_paths) >= k:
            break
        cand = xor_edges(base_edges, set(c))
        total_candidates += 1
        ok, _ = check_path_valid(cand, s, t)
        if ok:
            fz = frozenset(cand)
            if fz not in unique_paths:
                # Se Cenário C, verificar se contém via-points
                if required_nodes:
                    path_seq = edges_to_path(cand, s, t)
                    if path_seq and required_nodes.issubset(set(path_seq)):
                        unique_paths[fz] = cand
                        valid_count += 1
                else:
                    unique_paths[fz] = cand
                    valid_count += 1

    # 5) XOR duplo: base ⊕ c_i ⊕ c_j (se ainda precisar de mais caminhos)
    if len(unique_paths) < k and len(overlapping) >= 2:
        rng = np.random.default_rng(42)
        n_cyc = len(overlapping)
        all_pairs_count = n_cyc * (n_cyc - 1) // 2

        if all_pairs_count <= MAX_XOR_PAIRS:
            pair_iter = list(combinations(range(n_cyc), 2))
        else:
            seen_pairs = set()
            pair_iter = []
            attempts = 0
            while len(pair_iter) < MAX_XOR_PAIRS and attempts < MAX_XOR_PAIRS * 3:
                i = int(rng.integers(0, n_cyc))
                j = int(rng.integers(0, n_cyc))
                if i == j:
                    attempts += 1
                    continue
                pair = (min(i, j), max(i, j))
                if pair not in seen_pairs:
                    seen_pairs.add(pair)
                    pair_iter.append(pair)
                attempts += 1

        for (i, j) in pair_iter:
            if len(unique_paths) >= k:
                break
            cand = xor_edges(xor_edges(base_edges, set(overlapping[i])),
                             set(overlapping[j]))
            total_candidates += 1
            ok, _ = check_path_valid(cand, s, t)
            if ok:
                fz = frozenset(cand)
                if fz not in unique_paths:
                    if required_nodes:
                        path_seq = edges_to_path(cand, s, t)
                        if path_seq and required_nodes.issubset(set(path_seq)):
                            unique_paths[fz] = cand
                            valid_count += 1
                    else:
                        unique_paths[fz] = cand
                        valid_count += 1

    t_ms = (time.perf_counter() - t0) * 1000.0

    all_edges = list(unique_paths.values())
    stats = {
        "valid": valid_count,
        "total_candidates": total_candidates,
        "overlapping_cycles": len(overlapping),
        "total_cycles": len(cycles),
    }
    return all_edges, [p for p in all_edges], t_ms, stats


# ============================================================================
# CENÁRIO A: Labirinto Confinado (Maze via Spanning Tree)
# ============================================================================
def gen_maze_graph(n, seed):
    """Gera um labirinto n×n via árvore geradora aleatória (DFS randomizada).
    O grafo resultante é uma árvore (corredor estreito) com back-edges extras
    adicionadas para criar ciclos locais (~15% das arestas removidas de volta).
    """
    rng = np.random.default_rng(seed)

    # Grade completa
    G_full = nx.grid_2d_graph(n, n)

    # Spanning tree via DFS randomizada (cria corredores)
    visited = set()
    stack = [(0, 0)]
    tree_edges = []

    visited.add((0, 0))
    while stack:
        current = stack[-1]
        neighbors = [nb for nb in G_full.neighbors(current) if nb not in visited]
        if neighbors:
            rng.shuffle(neighbors)
            next_node = neighbors[0]
            visited.add(next_node)
            tree_edges.append((current, next_node))
            stack.append(next_node)
        else:
            stack.pop()

    # Criar grafo-labirinto
    G_maze = nx.Graph()
    G_maze.add_nodes_from(G_full.nodes())
    G_maze.add_edges_from(tree_edges)

    # Adicionar ~15% das arestas não-tree de volta (criar ciclos locais)
    non_tree = [(u, v) for u, v in G_full.edges()
                if not G_maze.has_edge(u, v)]
    n_extra = max(1, int(0.15 * len(non_tree)))
    if non_tree:
        idx = rng.choice(len(non_tree), size=min(n_extra, len(non_tree)),
                         replace=False)
        for i in idx:
            G_maze.add_edge(*non_tree[i])

    s = (0, 0)
    t = (n - 1, n - 1)
    return G_maze, s, t


def run_scenario_a(n=GRID_SIZE, n_seeds=N_SEEDS, k=10):
    """Cenário A: Labirinto Confinado. Retorna lista de resultados."""
    print(f"\n{'='*70}")
    print(f"CENÁRIO A: Labirinto Confinado {n}×{n} (k={k})")
    print(f"{'='*70}")
    results = []

    for seed in range(n_seeds):
        G, s, t = gen_maze_graph(n, seed)
        shortest = nx.shortest_path(G, s, t)
        sp_len = len(shortest) - 1  # número de arestas

        print(f"  [A] seed={seed}: V={G.number_of_nodes()} "
              f"E={G.number_of_edges()} sp_len={sp_len}")

        # --- Yen ---
        yen_paths, yen_ms, yen_found = run_yen(G, s, t, k)
        yen_div = jaccard_diversity(yen_paths) if yen_paths else 0.0
        yen_cost = cost_extension(yen_paths, sp_len) if yen_paths else 0.0

        # --- XOR ---
        xor_all, xor_valid, xor_ms, xor_stats = run_xor_pathfinding(G, s, t, k)
        xor_div = jaccard_diversity(xor_valid) if xor_valid else 0.0
        xor_cost = cost_extension(xor_valid, sp_len) if xor_valid else 0.0
        xor_validity = (xor_stats["valid"] / xor_stats["total_candidates"] * 100
                        if xor_stats["total_candidates"] > 0 else 0.0)

        rec = {
            "scenario": "A",
            "seed": seed,
            "grid_size": n,
            "k": k,
            "V": G.number_of_nodes(),
            "E": G.number_of_edges(),
            "shortest_path_len": sp_len,
            "yen_time_ms": yen_ms,
            "yen_paths_found": yen_found,
            "yen_diversity": yen_div,
            "yen_cost_extension": yen_cost,
            "xor_time_ms": xor_ms,
            "xor_paths_found": len(xor_valid),
            "xor_validity_rate": xor_validity,
            "xor_diversity": xor_div,
            "xor_cost_extension": xor_cost,
            "xor_total_candidates": xor_stats["total_candidates"],
            "xor_valid_candidates": xor_stats["valid"],
            "xor_overlapping_cycles": xor_stats["overlapping_cycles"],
            "xor_total_cycles": xor_stats["total_cycles"],
        }
        results.append(rec)
        print(f"       Yen: {yen_found} caminhos em {yen_ms:.1f}ms | "
              f"XOR: {len(xor_valid)} caminhos em {xor_ms:.1f}ms | "
              f"η={xor_validity:.1f}%")

    return results


# ============================================================================
# CENÁRIO B: Escalonamento de K (grade 20×20, 30% obstáculos)
# ============================================================================
def gen_obstacle_grid(n, density, seed):
    """Gera grade n×n com 'density' de obstáculos aleatórios."""
    s = (0, 0)
    t = (n - 1, n - 1)
    attempt = 0
    while True:
        cur_seed = seed + attempt
        rng = np.random.default_rng(cur_seed)
        blocked = rng.random((n, n)) < density
        blocked[s] = False
        blocked[t] = False
        G = nx.grid_2d_graph(n, n)
        G.remove_nodes_from([(r, c) for r in range(n) for c in range(n)
                             if blocked[r, c]])
        if t in G and nx.has_path(G, s, t):
            comp = nx.node_connected_component(G, s)
            return G.subgraph(comp).copy(), s, t, cur_seed
        attempt += 1
        if attempt > 200:
            raise RuntimeError(f"Não foi possível gerar grade {n}×{n} "
                               f"com densidade {density} e seed base {seed}")


def run_scenario_b(n=GRID_SIZE, density=OBSTACLE_DENSITY,
                   k_values=None, n_seeds=N_SEEDS):
    """Cenário B: Escalonamento de K. Retorna lista de resultados."""
    if k_values is None:
        k_values = K_VALUES

    print(f"\n{'='*70}")
    print(f"CENÁRIO B: Escalonamento de K ({n}×{n}, {density*100:.0f}% obst.)")
    print(f"  k ∈ {k_values}")
    print(f"{'='*70}")
    results = []

    for seed in range(n_seeds):
        G, s, t, used_seed = gen_obstacle_grid(n, density, seed * 1000)
        shortest = nx.shortest_path(G, s, t)
        sp_len = len(shortest) - 1

        # Pré-computar ciclos uma vez por grafo
        cycles = fundamental_cycles_dfs(G, s)

        print(f"  [B] seed={seed} (used={used_seed}): V={G.number_of_nodes()} "
              f"E={G.number_of_edges()} cycles={len(cycles)} sp_len={sp_len}")

        for k in k_values:
            # --- Yen ---
            yen_paths, yen_ms, yen_found = run_yen(G, s, t, k)
            yen_div = jaccard_diversity(yen_paths) if yen_paths else 0.0
            yen_cost = cost_extension(yen_paths, sp_len) if yen_paths else 0.0

            # --- XOR ---
            xor_all, xor_valid, xor_ms, xor_stats = run_xor_pathfinding(
                G, s, t, k, cycles=cycles)
            xor_div = jaccard_diversity(xor_valid) if xor_valid else 0.0
            xor_cost = cost_extension(xor_valid, sp_len) if xor_valid else 0.0
            xor_validity = (xor_stats["valid"] / xor_stats["total_candidates"] * 100
                            if xor_stats["total_candidates"] > 0 else 0.0)

            rec = {
                "scenario": "B",
                "seed": seed,
                "grid_size": n,
                "k": k,
                "V": G.number_of_nodes(),
                "E": G.number_of_edges(),
                "shortest_path_len": sp_len,
                "yen_time_ms": yen_ms,
                "yen_paths_found": yen_found,
                "yen_diversity": yen_div,
                "yen_cost_extension": yen_cost,
                "xor_time_ms": xor_ms,
                "xor_paths_found": len(xor_valid),
                "xor_validity_rate": xor_validity,
                "xor_diversity": xor_div,
                "xor_cost_extension": xor_cost,
                "xor_total_candidates": xor_stats["total_candidates"],
                "xor_valid_candidates": xor_stats["valid"],
                "xor_overlapping_cycles": xor_stats["overlapping_cycles"],
                "xor_total_cycles": xor_stats["total_cycles"],
            }
            results.append(rec)
            print(f"    k={k:3d}: Yen {yen_found:3d} em {yen_ms:8.1f}ms | "
                  f"XOR {len(xor_valid):3d} em {xor_ms:8.1f}ms | "
                  f"η={xor_validity:5.1f}%")

    return results


# ============================================================================
# CENÁRIO C: Vértices Compulsórios (Via-Points)
# ============================================================================
def gen_planar_with_waypoints(n, seed):
    """Gera um grafo planar (grade com ~20% de obstáculos) e seleciona
    2 via-points intermediários no caminho mais curto."""
    s = (0, 0)
    t = (n - 1, n - 1)
    attempt = 0
    while True:
        cur_seed = seed + attempt
        rng = np.random.default_rng(cur_seed)
        blocked = rng.random((n, n)) < 0.20  # menor densidade para garantir flexibilidade
        blocked[s] = False
        blocked[t] = False
        G = nx.grid_2d_graph(n, n)
        G.remove_nodes_from([(r, c) for r in range(n) for c in range(n)
                             if blocked[r, c]])
        if t in G and nx.has_path(G, s, t):
            comp = nx.node_connected_component(G, s)
            Gc = G.subgraph(comp).copy()
            # Selecionar via-points no terço 1/3 e 2/3 do shortest path
            sp = nx.shortest_path(Gc, s, t)
            if len(sp) >= 5:
                wp1 = sp[len(sp) // 3]
                wp2 = sp[2 * len(sp) // 3]
                if wp1 != s and wp1 != t and wp2 != s and wp2 != t and wp1 != wp2:
                    return Gc, s, t, {wp1, wp2}, cur_seed
        attempt += 1
        if attempt > 300:
            raise RuntimeError(f"Não foi possível gerar grafo planar com "
                               f"via-points válidos (seed base {seed})")


def yen_with_waypoints(G, s, t, waypoints, k, timeout_s=YEN_TIMEOUT_S):
    """Yen adaptado para vértices compulsórios: decompõe o problema em
    segmentos s→wp1→wp2→t e combina os k-shortest de cada segmento.
    """
    wp_list = sorted(waypoints, key=lambda wp: nx.shortest_path_length(G, s, wp))
    segments = [s] + wp_list + [t]

    t0 = time.perf_counter()
    segment_paths = []

    for i in range(len(segments) - 1):
        src, dst = segments[i], segments[i + 1]
        paths_seg, _, n_found = run_yen(G, src, dst, k, timeout_s=timeout_s)
        if n_found == 0:
            t_ms = (time.perf_counter() - t0) * 1000.0
            return [], t_ms, 0
        segment_paths.append(paths_seg)

    # Combinar: para cada combinação de segmentos, gerar rota completa
    combined = []
    # Gerar até k combinações via round-robin
    max_per_seg = [len(sp) for sp in segment_paths]
    count = 0
    for combo_idx in range(k * 2):  # tentar mais combinações
        if count >= k:
            break
        indices = []
        temp = combo_idx
        valid_combo = True
        for seg_i in range(len(segment_paths)):
            idx = temp % max_per_seg[seg_i]
            indices.append(idx)
            temp //= max_per_seg[seg_i]

        # Montar caminho completo (edges)
        full_edges = set()
        for seg_i in range(len(segment_paths)):
            full_edges |= segment_paths[seg_i][indices[seg_i]]

        fz = frozenset(full_edges)
        # Verificar se é caminho válido s→t
        ok, _ = check_path_valid(full_edges, s, t)
        if ok:
            # Verificar waypoints
            path_seq = edges_to_path(full_edges, s, t)
            if path_seq and waypoints.issubset(set(path_seq)):
                if fz not in {frozenset(c) for c in combined}:
                    combined.append(full_edges)
                    count += 1

    t_ms = (time.perf_counter() - t0) * 1000.0
    return combined, t_ms, len(combined)


def run_scenario_c(n=GRID_SIZE, n_seeds=N_SEEDS, k=20):
    """Cenário C: Vértices Compulsórios. Retorna lista de resultados."""
    print(f"\n{'='*70}")
    print(f"CENÁRIO C: Vértices Compulsórios {n}×{n} (k={k})")
    print(f"{'='*70}")
    results = []

    for seed in range(n_seeds):
        try:
            G, s, t, waypoints, used_seed = gen_planar_with_waypoints(n, seed * 1000)
        except RuntimeError as e:
            print(f"  [C] seed={seed}: SKIP - {e}")
            continue

        shortest = nx.shortest_path(G, s, t)
        sp_len = len(shortest) - 1

        print(f"  [C] seed={seed} (used={used_seed}): V={G.number_of_nodes()} "
              f"E={G.number_of_edges()} sp_len={sp_len} "
              f"waypoints={waypoints}")

        # --- Yen com waypoints ---
        yen_paths, yen_ms, yen_found = yen_with_waypoints(G, s, t, waypoints, k)
        yen_div = jaccard_diversity(yen_paths) if yen_paths else 0.0
        yen_cost = cost_extension(yen_paths, sp_len) if yen_paths else 0.0

        # --- XOR com via-points obrigatórios ---
        cycles = fundamental_cycles_dfs(G, s)
        xor_all, xor_valid, xor_ms, xor_stats = run_xor_pathfinding(
            G, s, t, k, cycles=cycles, required_nodes=waypoints)
        xor_div = jaccard_diversity(xor_valid) if xor_valid else 0.0
        xor_cost = cost_extension(xor_valid, sp_len) if xor_valid else 0.0
        xor_validity = (xor_stats["valid"] / xor_stats["total_candidates"] * 100
                        if xor_stats["total_candidates"] > 0 else 0.0)

        rec = {
            "scenario": "C",
            "seed": seed,
            "grid_size": n,
            "k": k,
            "V": G.number_of_nodes(),
            "E": G.number_of_edges(),
            "shortest_path_len": sp_len,
            "waypoints": [list(wp) for wp in waypoints],
            "yen_time_ms": yen_ms,
            "yen_paths_found": yen_found,
            "yen_diversity": yen_div,
            "yen_cost_extension": yen_cost,
            "xor_time_ms": xor_ms,
            "xor_paths_found": len(xor_valid),
            "xor_validity_rate": xor_validity,
            "xor_diversity": xor_div,
            "xor_cost_extension": xor_cost,
            "xor_total_candidates": xor_stats["total_candidates"],
            "xor_valid_candidates": xor_stats["valid"],
            "xor_overlapping_cycles": xor_stats["overlapping_cycles"],
            "xor_total_cycles": xor_stats["total_cycles"],
        }
        results.append(rec)
        print(f"       Yen: {yen_found} caminhos em {yen_ms:.1f}ms | "
              f"XOR: {len(xor_valid)} caminhos em {xor_ms:.1f}ms | "
              f"η={xor_validity:.1f}%")

    return results


# ============================================================================
# AGREGAÇÃO E RELATÓRIO
# ============================================================================
def aggregate_results(all_results):
    """Agrega os resultados por cenário para o relatório."""
    summary = {}

    for scenario in ["A", "B", "C"]:
        recs = [r for r in all_results if r["scenario"] == scenario]
        if not recs:
            continue

        if scenario == "B":
            # Agrupar por k
            by_k = defaultdict(list)
            for r in recs:
                by_k[r["k"]].append(r)

            summary[scenario] = {}
            for k_val, rs in sorted(by_k.items()):
                summary[scenario][f"k={k_val}"] = _aggregate_group(rs)
        else:
            summary[scenario] = _aggregate_group(recs)

    return summary


def _aggregate_group(recs):
    """Agrega métricas de um grupo de registros."""
    return {
        "n_instances": len(recs),
        "yen_time_ms_mean": float(np.mean([r["yen_time_ms"] for r in recs])),
        "yen_time_ms_std": float(np.std([r["yen_time_ms"] for r in recs])),
        "yen_paths_mean": float(np.mean([r["yen_paths_found"] for r in recs])),
        "xor_time_ms_mean": float(np.mean([r["xor_time_ms"] for r in recs])),
        "xor_time_ms_std": float(np.std([r["xor_time_ms"] for r in recs])),
        "xor_paths_mean": float(np.mean([r["xor_paths_found"] for r in recs])),
        "xor_validity_rate_mean": float(np.mean([r["xor_validity_rate"]
                                                  for r in recs])),
        "yen_diversity_mean": float(np.mean([r["yen_diversity"] for r in recs])),
        "xor_diversity_mean": float(np.mean([r["xor_diversity"] for r in recs])),
        "yen_cost_ext_mean": float(np.mean([r["yen_cost_extension"]
                                            for r in recs])),
        "xor_cost_ext_mean": float(np.mean([r["xor_cost_extension"]
                                            for r in recs])),
        "speedup_xor_vs_yen": (
            float(np.mean([r["yen_time_ms"] for r in recs]) /
                  np.mean([r["xor_time_ms"] for r in recs]))
            if np.mean([r["xor_time_ms"] for r in recs]) > 0 else float("inf")
        ),
    }


def print_report(summary, all_results):
    """Imprime relatório final no console."""
    print(f"\n{'='*70}")
    print("RELATÓRIO FINAL: Yen vs XOR Cycle-Space")
    print(f"{'='*70}")

    for scenario in ["A", "B", "C"]:
        if scenario not in summary:
            continue
        label = {"A": "Labirinto Confinado",
                 "B": "Escalonamento de K",
                 "C": "Vértices Compulsórios"}[scenario]
        print(f"\n--- Cenário {scenario}: {label} ---")

        if scenario == "B":
            for k_label, agg in summary[scenario].items():
                _print_agg(k_label, agg)
        else:
            _print_agg(f"Cenário {scenario}", summary[scenario])


def _print_agg(label, agg):
    """Imprime uma agregação."""
    print(f"\n  [{label}] ({agg['n_instances']} instâncias)")
    print(f"    Yen:  {agg['yen_time_ms_mean']:8.1f}ms ± "
          f"{agg['yen_time_ms_std']:6.1f}ms | "
          f"{agg['yen_paths_mean']:.1f} caminhos | "
          f"div={agg['yen_diversity_mean']:.3f} | "
          f"custo={agg['yen_cost_ext_mean']:.3f}")
    print(f"    XOR:  {agg['xor_time_ms_mean']:8.1f}ms ± "
          f"{agg['xor_time_ms_std']:6.1f}ms | "
          f"{agg['xor_paths_mean']:.1f} caminhos | "
          f"div={agg['xor_diversity_mean']:.3f} | "
          f"custo={agg['xor_cost_ext_mean']:.3f}")
    print(f"    Speedup XOR/Yen: {agg['speedup_xor_vs_yen']:.2f}x | "
          f"Validade XOR: {agg['xor_validity_rate_mean']:.1f}%")


# ============================================================================
# PLOTTING
# ============================================================================
def setup_plot_style():
    """Configura estilo premium para os gráficos."""
    plt.rcParams.update({
        "figure.facecolor": "#0d1117",
        "axes.facecolor": "#161b22",
        "axes.edgecolor": "#30363d",
        "axes.labelcolor": "#c9d1d9",
        "axes.titlesize": 14,
        "axes.labelsize": 12,
        "xtick.color": "#8b949e",
        "ytick.color": "#8b949e",
        "text.color": "#c9d1d9",
        "grid.color": "#21262d",
        "grid.alpha": 0.7,
        "legend.facecolor": "#161b22",
        "legend.edgecolor": "#30363d",
        "legend.fontsize": 10,
        "font.family": "sans-serif",
        "font.size": 11,
    })


def plot_scenario_b_time(results_b, output_dir):
    """Gráfico de linha: Tempo (ms) vs K para Yen e XOR."""
    setup_plot_style()
    fig, ax = plt.subplots(figsize=(10, 6))

    by_k = defaultdict(lambda: {"yen": [], "xor": []})
    for r in results_b:
        by_k[r["k"]]["yen"].append(r["yen_time_ms"])
        by_k[r["k"]]["xor"].append(r["xor_time_ms"])

    k_vals = sorted(by_k.keys())
    yen_means = [np.mean(by_k[k]["yen"]) for k in k_vals]
    yen_stds = [np.std(by_k[k]["yen"]) for k in k_vals]
    xor_means = [np.mean(by_k[k]["xor"]) for k in k_vals]
    xor_stds = [np.std(by_k[k]["xor"]) for k in k_vals]

    ax.errorbar(k_vals, yen_means, yerr=yen_stds,
                marker="o", linewidth=2.5, markersize=8,
                color="#ff6b6b", markerfacecolor="#ff6b6b",
                capsize=5, capthick=1.5,
                label="Yen (k-shortest paths)", zorder=5)
    ax.errorbar(k_vals, xor_means, yerr=xor_stds,
                marker="s", linewidth=2.5, markersize=8,
                color="#4ecdc4", markerfacecolor="#4ecdc4",
                capsize=5, capthick=1.5,
                label="XOR Cycle-Space", zorder=5)

    ax.set_xlabel("k (número de caminhos)", fontweight="bold")
    ax.set_ylabel("Tempo de Execução (ms)", fontweight="bold")
    ax.set_title("Cenário B: Tempo de Execução vs k\n"
                 "(Grade 20×20, 30% obstáculos)",
                 fontweight="bold", fontsize=14)
    ax.set_yscale("log")
    ax.legend(framealpha=0.9)
    ax.grid(True, alpha=0.3)
    ax.set_xticks(k_vals)

    fig.tight_layout()
    path = os.path.join(output_dir, "cenario_b_tempo_vs_k.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ Salvo: {path}")
    return path


def plot_scenario_b_diversity(results_b, output_dir):
    """Gráfico de linha: Diversidade de Jaccard vs K."""
    setup_plot_style()
    fig, ax = plt.subplots(figsize=(10, 6))

    by_k = defaultdict(lambda: {"yen": [], "xor": []})
    for r in results_b:
        by_k[r["k"]]["yen"].append(r["yen_diversity"])
        by_k[r["k"]]["xor"].append(r["xor_diversity"])

    k_vals = sorted(by_k.keys())
    yen_means = [np.mean(by_k[k]["yen"]) for k in k_vals]
    xor_means = [np.mean(by_k[k]["xor"]) for k in k_vals]

    ax.plot(k_vals, yen_means, marker="o", linewidth=2.5, markersize=8,
            color="#ff6b6b", label="Yen", zorder=5)
    ax.plot(k_vals, xor_means, marker="s", linewidth=2.5, markersize=8,
            color="#4ecdc4", label="XOR", zorder=5)

    ax.set_xlabel("k (número de caminhos)", fontweight="bold")
    ax.set_ylabel("Diversidade Média (Jaccard)", fontweight="bold")
    ax.set_title("Cenário B: Diversidade vs k", fontweight="bold", fontsize=14)
    ax.legend(framealpha=0.9)
    ax.grid(True, alpha=0.3)
    ax.set_xticks(k_vals)
    ax.set_ylim(0, 1.05)

    fig.tight_layout()
    path = os.path.join(output_dir, "cenario_b_diversidade_vs_k.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ Salvo: {path}")
    return path


def plot_structural_bars(all_results, output_dir):
    """Gráfico de barras agrupadas: comparação entre cenários A e C
    (métricas estruturais: tempo, diversidade, custo, validade)."""
    setup_plot_style()
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    scenario_data = {}
    for sc in ["A", "C"]:
        recs = [r for r in all_results if r["scenario"] == sc]
        if recs:
            scenario_data[sc] = {
                "yen_time": np.mean([r["yen_time_ms"] for r in recs]),
                "xor_time": np.mean([r["xor_time_ms"] for r in recs]),
                "yen_div": np.mean([r["yen_diversity"] for r in recs]),
                "xor_div": np.mean([r["xor_diversity"] for r in recs]),
                "yen_cost": np.mean([r["yen_cost_extension"] for r in recs]),
                "xor_cost": np.mean([r["xor_cost_extension"] for r in recs]),
                "xor_validity": np.mean([r["xor_validity_rate"] for r in recs]),
            }

    # Cenário B com k=50 (valor representativo)
    recs_b50 = [r for r in all_results
                if r["scenario"] == "B" and r["k"] == 50]
    if recs_b50:
        scenario_data["B(k=50)"] = {
            "yen_time": np.mean([r["yen_time_ms"] for r in recs_b50]),
            "xor_time": np.mean([r["xor_time_ms"] for r in recs_b50]),
            "yen_div": np.mean([r["yen_diversity"] for r in recs_b50]),
            "xor_div": np.mean([r["xor_diversity"] for r in recs_b50]),
            "yen_cost": np.mean([r["yen_cost_extension"] for r in recs_b50]),
            "xor_cost": np.mean([r["xor_cost_extension"] for r in recs_b50]),
            "xor_validity": np.mean([r["xor_validity_rate"] for r in recs_b50]),
        }

    labels = list(scenario_data.keys())
    x = np.arange(len(labels))
    width = 0.35

    # 1. Tempo de execução
    ax = axes[0, 0]
    yen_vals = [scenario_data[l]["yen_time"] for l in labels]
    xor_vals = [scenario_data[l]["xor_time"] for l in labels]
    bars1 = ax.bar(x - width / 2, yen_vals, width, color="#ff6b6b",
                   label="Yen", edgecolor="#30363d", linewidth=0.5)
    bars2 = ax.bar(x + width / 2, xor_vals, width, color="#4ecdc4",
                   label="XOR", edgecolor="#30363d", linewidth=0.5)
    ax.set_title("Tempo de Execução (ms)", fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.legend()
    ax.grid(axis="y", alpha=0.3)
    ax.set_yscale("log")
    # Adicionar valores nas barras
    for bar in bars1:
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                f"{bar.get_height():.0f}", ha="center", va="bottom",
                fontsize=8, color="#c9d1d9")
    for bar in bars2:
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                f"{bar.get_height():.0f}", ha="center", va="bottom",
                fontsize=8, color="#c9d1d9")

    # 2. Diversidade de Jaccard
    ax = axes[0, 1]
    yen_vals = [scenario_data[l]["yen_div"] for l in labels]
    xor_vals = [scenario_data[l]["xor_div"] for l in labels]
    ax.bar(x - width / 2, yen_vals, width, color="#ff6b6b",
           label="Yen", edgecolor="#30363d", linewidth=0.5)
    ax.bar(x + width / 2, xor_vals, width, color="#4ecdc4",
           label="XOR", edgecolor="#30363d", linewidth=0.5)
    ax.set_title("Diversidade Média (Jaccard)", fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.legend()
    ax.grid(axis="y", alpha=0.3)
    ax.set_ylim(0, 1.05)

    # 3. Custo de extensão
    ax = axes[1, 0]
    yen_vals = [scenario_data[l]["yen_cost"] for l in labels]
    xor_vals = [scenario_data[l]["xor_cost"] for l in labels]
    ax.bar(x - width / 2, yen_vals, width, color="#ff6b6b",
           label="Yen", edgecolor="#30363d", linewidth=0.5)
    ax.bar(x + width / 2, xor_vals, width, color="#4ecdc4",
           label="XOR", edgecolor="#30363d", linewidth=0.5)
    ax.set_title("Custo de Extensão (média / mínimo)", fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.legend()
    ax.grid(axis="y", alpha=0.3)

    # 4. Taxa de Validade XOR + Speedup
    ax = axes[1, 1]
    validity = [scenario_data[l]["xor_validity"] for l in labels]
    speedup = [scenario_data[l]["yen_time"] / scenario_data[l]["xor_time"]
               if scenario_data[l]["xor_time"] > 0 else 0
               for l in labels]
    bars_v = ax.bar(x - width / 2, validity, width, color="#45b7d1",
                    label="Validade XOR (%)", edgecolor="#30363d", linewidth=0.5)
    ax2 = ax.twinx()
    bars_s = ax2.bar(x + width / 2, speedup, width, color="#f7dc6f",
                     label="Speedup (×)", edgecolor="#30363d", linewidth=0.5)
    ax.set_title("Validade XOR & Speedup", fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Validade (%)", color="#45b7d1")
    ax2.set_ylabel("Speedup (×)", color="#f7dc6f")
    ax.grid(axis="y", alpha=0.3)
    # Legenda combinada
    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, loc="upper left")
    ax2.spines["right"].set_color("#f7dc6f")

    fig.suptitle("Benchmark Estrutural: Yen vs XOR Cycle-Space",
                 fontsize=16, fontweight="bold", y=1.02)
    fig.tight_layout()
    path = os.path.join(output_dir, "comparacao_estrutural.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ Salvo: {path}")
    return path


def plot_validity_heatmap(results_b, output_dir):
    """Gráfico adicional: Validade XOR e custo de extensão vs k."""
    setup_plot_style()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    by_k = defaultdict(lambda: {"validity": [], "cost_yen": [], "cost_xor": []})
    for r in results_b:
        by_k[r["k"]]["validity"].append(r["xor_validity_rate"])
        by_k[r["k"]]["cost_yen"].append(r["yen_cost_extension"])
        by_k[r["k"]]["cost_xor"].append(r["xor_cost_extension"])

    k_vals = sorted(by_k.keys())

    # Validade vs K
    means = [np.mean(by_k[k]["validity"]) for k in k_vals]
    stds = [np.std(by_k[k]["validity"]) for k in k_vals]
    ax1.bar(range(len(k_vals)), means, yerr=stds,
            color="#4ecdc4", edgecolor="#30363d", capsize=5,
            alpha=0.85, linewidth=0.5)
    ax1.set_xticks(range(len(k_vals)))
    ax1.set_xticklabels([str(k) for k in k_vals])
    ax1.set_xlabel("k", fontweight="bold")
    ax1.set_ylabel("Taxa de Validade XOR (%)", fontweight="bold")
    ax1.set_title("Validade XOR vs k", fontweight="bold")
    ax1.grid(axis="y", alpha=0.3)

    # Custo de extensão vs K
    yen_costs = [np.mean(by_k[k]["cost_yen"]) for k in k_vals]
    xor_costs = [np.mean(by_k[k]["cost_xor"]) for k in k_vals]
    x = np.arange(len(k_vals))
    width = 0.35
    ax2.bar(x - width / 2, yen_costs, width, color="#ff6b6b",
            label="Yen", edgecolor="#30363d", linewidth=0.5)
    ax2.bar(x + width / 2, xor_costs, width, color="#4ecdc4",
            label="XOR", edgecolor="#30363d", linewidth=0.5)
    ax2.set_xticks(x)
    ax2.set_xticklabels([str(k) for k in k_vals])
    ax2.set_xlabel("k", fontweight="bold")
    ax2.set_ylabel("Custo de Extensão (avg / min)", fontweight="bold")
    ax2.set_title("Custo de Extensão vs k", fontweight="bold")
    ax2.legend()
    ax2.grid(axis="y", alpha=0.3)

    fig.suptitle("Cenário B: Métricas Adicionais",
                 fontsize=14, fontweight="bold")
    fig.tight_layout()
    path = os.path.join(output_dir, "cenario_b_metricas_adicionais.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ Salvo: {path}")
    return path


# ============================================================================
# MAIN
# ============================================================================
def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)
    print("="*70)
    print("BENCHMARK: Yen's k-shortest paths vs GF(2) XOR Cycle-Space")
    print("="*70)
    print(f"Diretório de saída: {RESULTS_DIR}")
    print(f"Grade: {GRID_SIZE}×{GRID_SIZE} | Obstáculos: {OBSTACLE_DENSITY*100:.0f}%")
    print(f"Valores de k: {K_VALUES}")
    print(f"Seeds por cenário: {N_SEEDS}")
    print(f"Timeout Yen: {YEN_TIMEOUT_S}s")

    all_results = []

    # --- Cenário A ---
    try:
        results_a = run_scenario_a()
        all_results.extend(results_a)
    except Exception as e:
        print(f"\n[ERRO] Cenário A falhou: {e}")
        import traceback
        traceback.print_exc()

    # --- Cenário B ---
    try:
        results_b = run_scenario_b()
        all_results.extend(results_b)
    except Exception as e:
        print(f"\n[ERRO] Cenário B falhou: {e}")
        import traceback
        traceback.print_exc()

    # --- Cenário C ---
    try:
        results_c = run_scenario_c()
        all_results.extend(results_c)
    except Exception as e:
        print(f"\n[ERRO] Cenário C falhou: {e}")
        import traceback
        traceback.print_exc()

    # --- Salvar JSON ---
    json_path = os.path.join(RESULTS_DIR, "benchmark_results.json")
    summary = aggregate_results(all_results)
    output = {
        "metadata": {
            "grid_size": GRID_SIZE,
            "obstacle_density": OBSTACLE_DENSITY,
            "k_values": K_VALUES,
            "n_seeds": N_SEEDS,
            "yen_timeout_s": YEN_TIMEOUT_S,
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        },
        "summary": summary,
        "raw_results": all_results,
    }
    with open(json_path, "w") as f:
        json.dump(output, f, indent=2, default=str)
    print(f"\n✓ Resultados JSON salvos em: {json_path}")

    # --- Relatório ---
    print_report(summary, all_results)

    # --- Plots ---
    print(f"\n{'='*70}")
    print("GERANDO GRÁFICOS...")
    print(f"{'='*70}")

    results_b = [r for r in all_results if r["scenario"] == "B"]
    if results_b:
        plot_scenario_b_time(results_b, RESULTS_DIR)
        plot_scenario_b_diversity(results_b, RESULTS_DIR)
        plot_validity_heatmap(results_b, RESULTS_DIR)

    if any(r["scenario"] in ("A", "C") for r in all_results):
        plot_structural_bars(all_results, RESULTS_DIR)

    print(f"\n{'='*70}")
    print("BENCHMARK COMPLETO")
    print(f"Todos os resultados em: {RESULTS_DIR}/")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
