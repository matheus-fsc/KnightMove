"""Q4 — Mapa geométrico: onde D&C ganha vs BT direto.

Medições:
  1. BT direto 12×12: nós, profundidade, backtracks, status terminal por nó
  2. D&C: rejeições por bloco no cycle; ordem alternativa
  3. Tamanho efetivo do espaço de busca: combinações vs nós BT
"""

# --- bootstrap: torna core/ importavel a partir de qualquer subpasta ---
import sys as _sys, pathlib as _pathlib
for _anc in _pathlib.Path(__file__).resolve().parents:
    if (_anc / "core").is_dir() and (_anc / "pyproject.toml").exists():
        _sys.path.insert(0, str(_anc / "core"))
        break
# --- fim do bootstrap ---

import os
import time
import pickle
from collections import defaultdict, Counter

import numpy as np

import knight_tours as kt
import knight_tours_dnc as dnc


# ---------------------------------------------------------------------------
# Instrumentação do BT direto
# ---------------------------------------------------------------------------

class BTMetrics:
    __slots__ = ('nodes', 'leaves_tour', 'leaves_invalid',
                 'contradictions', 'subtours', 'completes',
                 'depth_at_tour', 'depth_at_failure', 'max_depth')

    def __init__(self):
        self.nodes = 0
        self.leaves_tour = 0
        self.leaves_invalid = 0
        self.contradictions = 0
        self.subtours = 0
        self.completes = 0
        self.depth_at_tour = []
        self.depth_at_failure = []
        self.max_depth = 0


def _backtrack_instr(state, ctx, rng, tours, K, metrics, depth=0):
    metrics.nodes += 1
    metrics.max_depth = max(metrics.max_depth, depth)
    if len(tours) >= K:
        return True

    edge, first_val = kt._choose_next_edge(state, ctx, rng)
    if edge == -1:
        if kt._is_complete_tour(state, ctx):
            tours.append(kt._extract_tour(state, ctx))
            metrics.leaves_tour += 1
            metrics.depth_at_tour.append(depth)
        else:
            metrics.leaves_invalid += 1
            metrics.depth_at_failure.append(depth)
        return len(tours) >= K

    for val in (first_val, 1 - first_val):
        snap = state.snapshot()
        status = kt.fix_and_propagate(state, ctx, edge, val)
        if status == kt.OK:
            if _backtrack_instr(state, ctx, rng, tours, K, metrics, depth + 1):
                state.restore(snap)
                return True
        elif status == kt.COMPLETE_TOUR:
            metrics.completes += 1
            if kt._is_complete_tour(state, ctx):
                tours.append(kt._extract_tour(state, ctx))
                metrics.leaves_tour += 1
                metrics.depth_at_tour.append(depth + 1)
            if len(tours) >= K:
                state.restore(snap)
                return True
        elif status == kt.CONTRADICTION:
            metrics.contradictions += 1
            metrics.depth_at_failure.append(depth + 1)
        elif status == kt.SUBTOUR:
            metrics.subtours += 1
            metrics.depth_at_failure.append(depth + 1)
        state.restore(snap)
    return len(tours) >= K


def knight_tours_instrumented(n, K, seed=0):
    if K == 0:
        return [], BTMetrics()
    ctx = kt.build_graph(n)
    rng = np.random.default_rng(seed)
    state = kt.State(ctx['V'], ctx['E'])
    status = kt._propagate_initial(state, ctx)
    if status in (kt.CONTRADICTION, kt.SUBTOUR):
        return [], BTMetrics()
    metrics = BTMetrics()
    tours = []
    if status == kt.COMPLETE_TOUR:
        if kt._is_complete_tour(state, ctx):
            tours.append(kt._extract_tour(state, ctx))
            metrics.leaves_tour += 1
    else:
        _backtrack_instr(state, ctx, rng, tours, K, metrics, 0)
    return tours[:K], metrics


# ---------------------------------------------------------------------------
# D&C com contagem de rejeições por bloco
# ---------------------------------------------------------------------------

def matching_with_block_rejections(catalog, cross, K, rng_seed=0,
                                    order='perimeter', max_attempts=None):
    """Mede em qual bloco do cycle a rejeição mais ocorre.

    order ∈ {'perimeter': BL→BR→TR→TL, 'alt': BL→TL→TR→BR}.
    """
    rng = np.random.default_rng(rng_seed)
    if order == 'perimeter':
        cycle = ['BL', 'BR', 'TR', 'TL']
    elif order == 'alt':
        cycle = ['BL', 'TL', 'TR', 'BR']
    else:
        raise ValueError(order)

    n_blocks = len(cycle)
    fb = {b: set(dnc.get_frontier_vertices(b)) for b in cycle}

    # Cross indices direcionais — populando ambas direções para todos
    # os pares de blocos no ciclo
    cross_idx = {}
    for i in range(n_blocks):
        for j in range(n_blocks):
            if i == j:
                continue
            a = cycle[i]
            b = cycle[j]
            if (a, b) in cross_idx:
                continue
            edges_global = dnc.get_cross_edges(a, b)
            idx = defaultdict(list)
            for (ga, gb) in edges_global:
                la = dnc.global_to_local(ga, a)
                lb = dnc.global_to_local(gb, b)
                idx[la].append(lb)
            cross_idx[(a, b)] = idx

    def frontier_to_block(block_a, block_b):
        return set(cross_idx[(block_a, block_b)].keys())

    pairs_by_block = {}
    by_start_by_block = {}
    for i in range(n_blocks):
        a = cycle[i]
        b_prev = cycle[(i - 1) % n_blocks]
        b_next = cycle[(i + 1) % n_blocks]
        fr_p = frontier_to_block(a, b_prev)
        fr_n = frontier_to_block(a, b_next)
        pairs = [(s, e) for (s, e) in catalog
                 if s in fr_p and e in fr_n and catalog[(s, e)]]
        pairs_by_block[a] = pairs
        bs = defaultdict(list)
        for (s, e) in pairs:
            bs[s].append(e)
        by_start_by_block[a] = bs

    tours = []
    attempts = 0
    rejections = Counter()  # bloco → count
    if max_attempts is None:
        max_attempts = K * 100

    while len(tours) < K and attempts < max_attempts:
        attempts += 1
        # bloco 0
        a0 = cycle[0]
        pairs0 = pairs_by_block[a0]
        if not pairs0:
            rejections[a0] += 1
            continue
        s0, e0 = pairs0[rng.integers(len(pairs0))]
        endpoints = [(s0, e0)]
        feasible = True
        for i in range(1, n_blocks):
            a_prev = cycle[i - 1]
            a_curr = cycle[i]
            e_prev = endpoints[-1][1]
            cand_s = [v for v in cross_idx[(a_prev, a_curr)][e_prev]
                      if v in fb[a_curr]]
            if not cand_s:
                rejections[a_curr] += 1
                feasible = False
                break
            s_curr = cand_s[rng.integers(len(cand_s))]
            if i == n_blocks - 1:
                a_close = cycle[0]
                cand_e = [
                    e for e in by_start_by_block[a_curr].get(s_curr, [])
                    if s0 in cross_idx[(a_curr, a_close)].get(e, [])
                ]
            else:
                cand_e = by_start_by_block[a_curr].get(s_curr, [])
            if not cand_e:
                rejections[a_curr] += 1
                feasible = False
                break
            e_curr = cand_e[rng.integers(len(cand_e))]
            endpoints.append((s_curr, e_curr))
        if not feasible:
            continue

        # Montar tour
        block_paths = []
        for i, b in enumerate(cycle):
            s, e = endpoints[i]
            paths = catalog[(s, e)]
            p = paths[rng.integers(len(paths))]
            block_paths.append(p)

        # Concatenar paths na ordem do cycle, traduzindo para global
        seq = []
        for i, b in enumerate(cycle):
            for v in block_paths[i]:
                seq.append(dnc.local_to_global(v, b))
        tour = np.asarray(seq, dtype=np.int32)
        if kt.verify_tour(tour, dnc.N_GLOBAL):
            tours.append(tour)
        else:
            rejections['verify'] += 1

    return tours, attempts, rejections


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    print("=" * 70)
    print(" Q4 — Mapa geométrico D&C vs BT")
    print("=" * 70)

    # ---- BT direto 12×12 instrumentado
    K_bt = 50
    print(f"\n[1] BT direto 12×12 (K={K_bt}, seed=0)")
    t0 = time.time()
    tours_bt, metrics = knight_tours_instrumented(12, K_bt, seed=0)
    t1 = time.time()
    print(f"  tempo: {t1 - t0:.3f}s")
    print(f"  tours: {len(tours_bt)}, válidos: "
          f"{sum(1 for t in tours_bt if kt.verify_tour(t, 12))}")
    print(f"  nós: {metrics.nodes}")
    print(f"  folhas-tour: {metrics.leaves_tour}, "
          f"folhas-invalidas: {metrics.leaves_invalid}")
    print(f"  contradições: {metrics.contradictions}, "
          f"subtours: {metrics.subtours}, "
          f"complete-via-prop: {metrics.completes}")
    print(f"  max_depth: {metrics.max_depth}")
    if metrics.depth_at_tour:
        d = metrics.depth_at_tour
        print(f"  profundidade média ao encontrar tour: "
              f"{sum(d)/len(d):.2f} (min={min(d)}, max={max(d)})")
    if metrics.depth_at_failure:
        d = metrics.depth_at_failure
        print(f"  profundidade média de falha: "
              f"{sum(d)/len(d):.2f} (min={min(d)}, max={max(d)})")
    # Por bin
    if metrics.depth_at_failure:
        bins = Counter()
        for d in metrics.depth_at_failure:
            bins[d] += 1
        if max(bins.keys()) <= 10:
            for d in sorted(bins.keys()):
                print(f"    depth={d}: {bins[d]} falhas")

    # ---- D&C com order=perimeter e order=alt
    if not os.path.exists(dnc.CATALOG_FILE):
        raise SystemExit("Catálogo não existe.")
    with open(dnc.CATALOG_FILE, 'rb') as f:
        obj = pickle.load(f)
    catalog = obj['catalog']
    cross = dnc.precompute_cross_edge_index()

    K_dnc = 500
    print(f"\n[2] D&C cycle PERIMETER (BL→BR→TR→TL), K={K_dnc}")
    t0 = time.time()
    tours, attempts, rej = matching_with_block_rejections(
        catalog, cross, K_dnc, rng_seed=0, order='perimeter')
    t1 = time.time()
    print(f"  tempo: {t1 - t0:.3f}s | tours: {len(tours)}, "
          f"tentativas: {attempts}")
    print(f"  rejeições por bloco/etapa: {dict(rej)}")

    print(f"\n[3] D&C cycle ALT (BL→TL→TR→BR), K={K_dnc}")
    t0 = time.time()
    tours_alt, attempts_alt, rej_alt = matching_with_block_rejections(
        catalog, cross, K_dnc, rng_seed=0, order='alt')
    t1 = time.time()
    print(f"  tempo: {t1 - t0:.3f}s | tours: {len(tours_alt)}, "
          f"tentativas: {attempts_alt}")
    print(f"  rejeições por bloco/etapa: {dict(rej_alt)}")

    # ---- Comparação Q1 vs Q4: papel do pré-filtro
    print(f"\n[3b] Comparação com filtro 'frouxo' (s_BL não restrito a F_to_TL)")
    import dnc_q1_diversity
    t0 = time.time()
    t_lous, attempts_lous, fails_lous = dnc_q1_diversity.random_matching(
        catalog, cross, K_dnc, rng_seed=0)
    t1 = time.time()
    print(f"  tours: {len(t_lous)}, tentativas: {attempts_lous}, "
          f"falhas: {fails_lous} (rejeição={fails_lous/attempts_lous:.1%})")
    print(f"  ⇒ Pré-filtro direcionado (Q4) elimina ~{fails_lous} rejeições.")

    # ---- Tamanho do espaço de busca
    print(f"\n[4] Espaço de busca")
    print(f"  BT direto: {metrics.nodes} nós para {len(tours_bt)} tours "
          f"= {metrics.nodes/max(len(tours_bt),1):.1f} nós/tour")
    print(f"  Profundidade média de cada nó BT: ~{(sum(metrics.depth_at_tour)+sum(metrics.depth_at_failure))/(len(metrics.depth_at_tour)+max(len(metrics.depth_at_failure),1)):.0f} arestas decididas")
    print(f"  D&C: catalog tem 25600 paths; 4 blocos com ~70 pares cada;")
    print(f"        combinações teóricas brutas ≈ 25600^4 ≈ 4.3e17,")
    print(f"        mas filtro de fronteira direcionado reduz a {attempts}/K={K_dnc} = "
          f"{attempts/K_dnc:.1f} tentativas/tour")
    print(f"  Taxa de sucesso D&C: {len(tours)/attempts:.1%}")
    print()
    print(f"  Razão exploração:")
    print(f"    BT direto:   {metrics.nodes/max(len(tours_bt),1):.2f} nós/tour × ~{metrics.max_depth} níveis")
    print(f"    D&C filtrado: {attempts/K_dnc:.2f} tentativas/tour × 4 blocos × O(1) lookup")
