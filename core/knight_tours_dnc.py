"""knight_tours_dnc.py — Divide-and-Conquer 12×12 com blocos 6×6.

Esquema:
    FASE 0: Setup do tabuleiro 12×12 em 4 blocos 6×6 (BL, BR, TR, TL).
    FASE 1: Pré-computação de paths hamiltonianos abertos no 6×6
            isolado, indexados por par (s, e) de vértices de fronteira.
    FASE 2: Combinação dos paths dos 4 blocos via knight edges
            cruzando fronteira, com critério Union-Find de ciclo único.
    FASE 3: Benchmark vs backtracking direto em 12×12.

Geometria:
    BL: linhas 0-5, colunas 0-5
    BR: linhas 0-5, colunas 6-11
    TL: linhas 6-11, colunas 0-5
    TR: linhas 6-11, colunas 6-11

Convenção de coordenadas:
    Vértice local v ∈ {0,...,35} mapeia para (lr, lc) = divmod(v, 6).
    Vértice global g ∈ {0,...,143} mapeia para (gr, gc) = divmod(g, 12).
"""

# --- bootstrap: torna core/ importavel a partir de qualquer subpasta ---
import sys as _sys, pathlib as _pathlib
for _anc in _pathlib.Path(__file__).resolve().parents:
    if (_anc / "core").is_dir() and (_anc / "pyproject.toml").exists():
        _sys.path.insert(0, str(_anc / "core"))
        break
# --- fim do bootstrap ---

import sys
import time
import os
import pickle
from collections import defaultdict

import numpy as np

import knight_tours as kt

sys.setrecursionlimit(20000)


KNIGHT_MOVES = ((-2, -1), (-2, 1), (-1, -2), (-1, 2),
                (1, -2), (1, 2), (2, -1), (2, 1))

N_BLOCK = 6
N_GLOBAL = 12
V_BLOCK = N_BLOCK * N_BLOCK
V_GLOBAL = N_GLOBAL * N_GLOBAL

BLOCKS = {
    'BL': (0, 0),
    'BR': (0, 6),
    'TL': (6, 0),
    'TR': (6, 6),
}

# Ordem cíclica perimétrica (sentido horário a partir de BL)
PERIMETER_ORDER = ['BL', 'BR', 'TR', 'TL']


# ---------------------------------------------------------------------------
# FASE 0 — Geometria dos blocos
# ---------------------------------------------------------------------------

def local_to_global(local_v, block_name):
    or_, oc_ = BLOCKS[block_name]
    lr, lc = divmod(local_v, N_BLOCK)
    return (or_ + lr) * N_GLOBAL + (oc_ + lc)


def global_to_block(global_v):
    gr, gc = divmod(global_v, N_GLOBAL)
    br = (gr // N_BLOCK) * N_BLOCK
    bc = (gc // N_BLOCK) * N_BLOCK
    for name, (or_, oc_) in BLOCKS.items():
        if or_ == br and oc_ == bc:
            return name
    return None


def global_to_local(global_v, block_name):
    or_, oc_ = BLOCKS[block_name]
    gr, gc = divmod(global_v, N_GLOBAL)
    return (gr - or_) * N_BLOCK + (gc - oc_)


def get_frontier_vertices(block_name):
    """Vértices locais com ≥1 movimento de cavalo saindo do bloco
    e dentro do tabuleiro global."""
    or_, oc_ = BLOCKS[block_name]
    frontier = set()
    for r in range(N_BLOCK):
        for c in range(N_BLOCK):
            for dr, dc in KNIGHT_MOVES:
                gr = or_ + r + dr
                gc = oc_ + c + dc
                if not (0 <= gr < N_GLOBAL and 0 <= gc < N_GLOBAL):
                    continue
                lr = gr - or_
                lc = gc - oc_
                if 0 <= lr < N_BLOCK and 0 <= lc < N_BLOCK:
                    continue
                frontier.add(r * N_BLOCK + c)
                break
    return sorted(frontier)


def get_cross_edges(block_a, block_b):
    """Knight edges globais conectando block_a a block_b.
    Cada aresta é (v_global_em_A, v_global_em_B)."""
    or_a, oc_a = BLOCKS[block_a]
    or_b, oc_b = BLOCKS[block_b]
    edges = []
    for ra in range(N_BLOCK):
        for ca in range(N_BLOCK):
            gra = or_a + ra
            gca = oc_a + ca
            for dr, dc in KNIGHT_MOVES:
                grb = gra + dr
                gcb = gca + dc
                if not (0 <= grb < N_GLOBAL and 0 <= gcb < N_GLOBAL):
                    continue
                if or_b <= grb < or_b + N_BLOCK and oc_b <= gcb < oc_b + N_BLOCK:
                    edges.append((gra * N_GLOBAL + gca,
                                  grb * N_GLOBAL + gcb))
    return edges


# ---------------------------------------------------------------------------
# FASE 1 — Enumerador de paths hamiltonianos no 6×6 (R2 + UF via kt)
# ---------------------------------------------------------------------------

def _color(v):
    r, c = divmod(v, N_BLOCK)
    return (r + c) & 1


def _build_modified_ctx(s, e, n=N_BLOCK):
    """Retorna (ctx, virtual_edge_index). Acrescenta aresta virtual (s,e) se
    não existir; ctx é uma cópia mutável independente."""
    base = kt.build_graph(n)
    V = base['V']
    E = base['E']

    s2, e2 = min(s, e), max(s, e)
    existing = -1
    eps = base['edge_endpoints']
    for ei in range(E):
        if int(eps[ei, 0]) == s2 and int(eps[ei, 1]) == e2:
            existing = ei
            break

    if existing >= 0:
        return base, existing

    # Adiciona aresta virtual
    new_eps = np.vstack(
        [eps, np.array([[s2, e2]], dtype=np.int32)])
    new_adj = list(base['adj_edges'])
    new_adj[s] = np.append(new_adj[s], np.int32(E))
    new_adj[e] = np.append(new_adj[e], np.int32(E))
    new_inc = base['total_incident'].copy()
    new_inc[s] += 1
    new_inc[e] += 1
    ctx = {
        'n': n,
        'V': V,
        'E': E + 1,
        'edge_endpoints': new_eps,
        'adj_edges': new_adj,
        'total_incident': new_inc,
    }
    return ctx, E


def _cycle_to_path(tour, s, e):
    """Dado um ciclo fechado contendo (virtualmente) a aresta (s,e), retorna a
    tupla de vértices do caminho aberto s → ... → e."""
    V = len(tour)
    idx_s = -1
    for i in range(V):
        if int(tour[i]) == s:
            idx_s = i
            break
    nxt = int(tour[(idx_s + 1) % V])
    prv = int(tour[(idx_s - 1) % V])
    if prv == e:
        # ciclo: ... → e → s → next → ... → e (virtual edge é (e,s))
        # removendo: s → next → ... → e
        path = [int(tour[(idx_s + i) % V]) for i in range(V)]
        return tuple(path)
    if nxt == e:
        # ciclo: s → e → ... → s (virtual edge é (s,e))
        # removendo: caminho deveria ser s → prv → ... → e (traverse backward)
        path = [int(tour[(idx_s - i) % V]) for i in range(V)]
        return tuple(path)
    raise RuntimeError(f"e={e} não adjacente a s={s} no tour {list(map(int, tour))}")


def _enum_backtrack(state, ctx, rng, paths, s, e, max_paths, nodes):
    if max_paths is not None and len(paths) >= max_paths:
        return True
    nodes[0] += 1
    edge, first_val = kt._choose_next_edge(state, ctx, rng)
    if edge == -1:
        if kt._is_complete_tour(state, ctx):
            tour = kt._extract_tour(state, ctx)
            paths.append(_cycle_to_path(tour, s, e))
        return max_paths is not None and len(paths) >= max_paths

    for val in (first_val, 1 - first_val):
        snap = state.snapshot()
        status = kt.fix_and_propagate(state, ctx, edge, val)
        if status == kt.OK:
            if _enum_backtrack(state, ctx, rng, paths, s, e, max_paths, nodes):
                state.restore(snap)
                return True
        elif status == kt.COMPLETE_TOUR:
            if kt._is_complete_tour(state, ctx):
                tour = kt._extract_tour(state, ctx)
                paths.append(_cycle_to_path(tour, s, e))
        state.restore(snap)
        if max_paths is not None and len(paths) >= max_paths:
            return True
    return False


def enumerate_block_paths(s, e, max_paths=None, seed=0):
    """Enumera HPs de s a e no 6×6 isolado, via aresta virtual + R2 + UF do kt."""
    if _color(s) == _color(e):
        return []  # impossível em bipartido com V par
    ctx, virtual = _build_modified_ctx(s, e, n=N_BLOCK)
    V = ctx['V']
    state = kt.State(V, ctx['E'])

    # Força aresta virtual ACTIVE e propaga R2 dos cantos
    queue = [(virtual, 1)]
    total_inc = ctx['total_incident']
    for w in range(V):
        if int(total_inc[w]) == 2:
            for ei in ctx['adj_edges'][w]:
                if int(ei) != virtual:
                    queue.append((int(ei), 1))
    status = kt._process_queue(state, ctx, queue)
    if status in (kt.CONTRADICTION, kt.SUBTOUR):
        return []
    paths = []
    if status == kt.COMPLETE_TOUR:
        if kt._is_complete_tour(state, ctx):
            tour = kt._extract_tour(state, ctx)
            paths.append(_cycle_to_path(tour, s, e))
        return paths

    rng = np.random.default_rng(seed)
    nodes = [0]
    _enum_backtrack(state, ctx, rng, paths, s, e, max_paths, nodes)
    return paths


def build_catalog(frontier_union, max_paths_per_pair=None, verbose=True):
    """Constrói catálogo de HPs entre todos os pares (s, e) com s, e ∈ frontier_union
    e cor(s) ≠ cor(e).

    Retorna dict: (s, e) → lista de tuplas de vértices.
    """
    catalog = {}
    pairs_skipped_color = 0
    pairs_empty = 0
    pairs_full = 0
    total_paths = 0

    pairs = []
    for s in frontier_union:
        for e in frontier_union:
            if s == e:
                continue
            if _color(s) == _color(e):
                pairs_skipped_color += 1
                continue
            pairs.append((s, e))

    n_pairs = len(pairs)
    t0 = time.time()
    for idx, (s, e) in enumerate(pairs):
        # Para evitar duplicar (s, e) e (e, s): cada par dirigido tem HPs distintos
        # (reverso do path do outro), mas armazenamos ambos para simplicidade de busca.
        # Otimização: enumerar apenas s < e e armazenar reversos para (e, s).
        if s < e:
            paths = enumerate_block_paths(s, e, max_paths=max_paths_per_pair)
            catalog[(s, e)] = paths
            catalog[(e, s)] = [p[::-1] for p in paths]
            total_paths += 2 * len(paths)
            if len(paths) == 0:
                pairs_empty += 2
            if max_paths_per_pair is not None and len(paths) >= max_paths_per_pair:
                pairs_full += 2
        if verbose and (idx + 1) % 50 == 0:
            elapsed = time.time() - t0
            print(f"  [{idx+1}/{n_pairs}] (s,e)=({s},{e}) paths_acum={total_paths} "
                  f"tempo={elapsed:.1f}s")
    return catalog, {
        'n_pairs_processed': len(catalog),
        'pairs_skipped_color': pairs_skipped_color,
        'pairs_empty': pairs_empty,
        'pairs_full': pairs_full,
        'total_paths': total_paths,
    }


# ---------------------------------------------------------------------------
# FASE 2 — Combinação de blocos
# ---------------------------------------------------------------------------

def precompute_cross_edge_index():
    """Para cada par ordenado (block_a, block_b) ≠, devolve:
       cross[(block_a, block_b)] = lista de (s_local_em_a, e_local_em_b)
       de pares conectados por knight edge cruzando fronteira.
    """
    cross = {}
    for a in BLOCKS:
        for b in BLOCKS:
            if a == b:
                continue
            edges_global = get_cross_edges(a, b)
            edges_local = [(global_to_local(ga, a), global_to_local(gb, b))
                           for (ga, gb) in edges_global]
            cross[(a, b)] = edges_local
    return cross


def combine_perimeter(catalog, cross, max_tours=50, verbose=True):
    """Busca tours fechados no 12×12 combinando 4 paths internos
    (um por bloco) na ordem perimétrica BL → BR → TR → TL → BL.

    Métrica: # de "nós explorados" = # de combinações (e_BL, s_BR, e_BR, s_TR,
    e_TR, s_TL, e_TL, s_BL_close) testadas.
    """
    tours = []
    nodes_explored = 0
    order = PERIMETER_ORDER  # [BL, BR, TR, TL]

    # Pré-computar índices por bloco de pares (s, e) disponíveis no catálogo
    # Como o catálogo é universal (paths no 6×6 isolado), apenas filtramos por
    # frontiers do bloco em questão.
    frontier_per_block = {b: set(get_frontier_vertices(b)) for b in BLOCKS}

    # Pre-index: paths por start vertex em cada bloco
    paths_by_start = defaultdict(list)
    for (s, e), plist in catalog.items():
        if not plist:
            continue
        paths_by_start[s].append((e, plist))

    # Cross edges
    bl_to_br = cross[('BL', 'BR')]
    br_to_tr = cross[('BR', 'TR')]
    tr_to_tl = cross[('TR', 'TL')]
    tl_to_bl = cross[('TL', 'BL')]

    # Construir índices: dado vértice local-saída em A, vértices locais de entrada em B
    def index_cross(edges):
        idx = defaultdict(list)
        for (sa, sb) in edges:
            idx[sa].append(sb)
        return idx

    idx_BL_BR = index_cross(bl_to_br)
    idx_BR_TR = index_cross(br_to_tr)
    idx_TR_TL = index_cross(tr_to_tl)
    idx_TL_BL = index_cross(tl_to_bl)

    fb_BL = frontier_per_block['BL']
    fb_BR = frontier_per_block['BR']
    fb_TR = frontier_per_block['TR']
    fb_TL = frontier_per_block['TL']

    # Busca: iterar s_BL ∈ frontier_BL (a entrada de BL vem do último cross-edge,
    # TL → BL). Para cada s_BL, enumerar paths internos em BL.
    for s_BL in fb_BL:
        if s_BL not in paths_by_start:
            continue
        for (e_BL, paths_BL) in paths_by_start[s_BL]:
            if e_BL not in fb_BL or e_BL not in idx_BL_BR:
                continue
            # paths_BL[i] usa este (s_BL, e_BL); todos têm 36 vértices
            if not paths_BL:
                continue

            # Cross BL → BR
            for s_BR in idx_BL_BR[e_BL]:
                if s_BR not in fb_BR or s_BR not in paths_by_start:
                    continue
                for (e_BR, paths_BR) in paths_by_start[s_BR]:
                    if e_BR not in fb_BR or e_BR not in idx_BR_TR:
                        continue
                    if not paths_BR:
                        continue

                    # Cross BR → TR
                    for s_TR in idx_BR_TR[e_BR]:
                        if s_TR not in fb_TR or s_TR not in paths_by_start:
                            continue
                        for (e_TR, paths_TR) in paths_by_start[s_TR]:
                            if e_TR not in fb_TR or e_TR not in idx_TR_TL:
                                continue
                            if not paths_TR:
                                continue

                            # Cross TR → TL
                            for s_TL in idx_TR_TL[e_TR]:
                                if s_TL not in fb_TL or s_TL not in paths_by_start:
                                    continue
                                for (e_TL, paths_TL) in paths_by_start[s_TL]:
                                    if e_TL not in fb_TL:
                                        continue
                                    # Fechar: TL → BL deve usar (e_TL, s_BL)
                                    if e_TL not in idx_TL_BL:
                                        continue
                                    if s_BL not in idx_TL_BL[e_TL]:
                                        continue
                                    if not paths_TL:
                                        continue

                                    nodes_explored += 1

                                    # Encontrado um esqueleto válido!
                                    # Pegar UM path por bloco e montar tour.
                                    p_BL = paths_BL[0]
                                    p_BR = paths_BR[0]
                                    p_TR = paths_TR[0]
                                    p_TL = paths_TL[0]

                                    tour = _assemble_tour(
                                        p_BL, p_BR, p_TR, p_TL)
                                    if tour is not None:
                                        tours.append(tour)
                                        if verbose and len(tours) <= 10:
                                            print(f"  Tour {len(tours)}: "
                                                  f"BL({s_BL}→{e_BL}), "
                                                  f"BR({s_BR}→{e_BR}), "
                                                  f"TR({s_TR}→{e_TR}), "
                                                  f"TL({s_TL}→{e_TL})")
                                        if len(tours) >= max_tours:
                                            return tours, nodes_explored
    return tours, nodes_explored


def _assemble_tour(p_BL, p_BR, p_TR, p_TL):
    """Concatena paths locais traduzidos para global na ordem
    BL → BR → TR → TL. Retorna array global ou None se algo falhou."""
    seq = []
    for v in p_BL:
        seq.append(local_to_global(v, 'BL'))
    for v in p_BR:
        seq.append(local_to_global(v, 'BR'))
    for v in p_TR:
        seq.append(local_to_global(v, 'TR'))
    for v in p_TL:
        seq.append(local_to_global(v, 'TL'))
    return np.asarray(seq, dtype=np.int32)


# ---------------------------------------------------------------------------
# Cache em disco
# ---------------------------------------------------------------------------

CATALOG_FILE = 'dnc_catalog_6x6.pkl'


def load_or_build_catalog(frontier_union, max_paths_per_pair=None,
                           force=False, verbose=True):
    if not force and os.path.exists(CATALOG_FILE):
        with open(CATALOG_FILE, 'rb') as f:
            obj = pickle.load(f)
        if (obj.get('max_paths_per_pair') == max_paths_per_pair and
                obj.get('frontier_union') == frontier_union):
            if verbose:
                print(f"  catálogo carregado de {CATALOG_FILE}")
            return obj['catalog'], obj['stats']
    if verbose:
        print(f"  construindo catálogo (max_paths_per_pair={max_paths_per_pair}) ...")
    t0 = time.time()
    catalog, stats = build_catalog(frontier_union, max_paths_per_pair, verbose=verbose)
    t1 = time.time()
    stats['build_time_s'] = t1 - t0
    with open(CATALOG_FILE, 'wb') as f:
        pickle.dump({
            'frontier_union': frontier_union,
            'max_paths_per_pair': max_paths_per_pair,
            'catalog': catalog,
            'stats': stats,
        }, f)
    if verbose:
        print(f"  catálogo construído em {t1 - t0:.1f}s, salvo em {CATALOG_FILE}")
    return catalog, stats


# ---------------------------------------------------------------------------
# FASE 3 — Benchmark
# ---------------------------------------------------------------------------

def verify_assembled_tour(tour, n=N_GLOBAL):
    """Wrapper sobre kt.verify_tour."""
    return kt.verify_tour(tour, n)


def benchmark_dnc(K=50, max_paths_per_pair=200, verbose=True):
    """Pipeline completo: catalog + matching + extract K tours."""
    print("=" * 70)
    print(" FASE 0 — Geometria")
    print("=" * 70)
    fronts = {b: get_frontier_vertices(b) for b in BLOCKS}
    for b, f in fronts.items():
        print(f"  {b}: |F|={len(f)} frontiers={f}")
    cross = precompute_cross_edge_index()
    for (a, b), eds in sorted(cross.items()):
        if a in PERIMETER_ORDER and b in PERIMETER_ORDER:
            i = PERIMETER_ORDER.index(a)
            j = PERIMETER_ORDER.index(b)
            if (j - i) % 4 == 1:
                print(f"  cross[{a}→{b}]: {len(eds)} arestas")
    print()

    # União das fronteiras
    frontier_union = sorted(set().union(*fronts.values()))
    print(f"  União de fronteiras: {len(frontier_union)} vértices")
    print(f"  Pares ordenados candidatos (cores opostas): "
          f"{sum(1 for s in frontier_union for e in frontier_union if s != e and _color(s) != _color(e))}")
    print()

    print("=" * 70)
    print(" FASE 1 — Catálogo de paths no 6×6")
    print("=" * 70)
    t_phase1 = time.time()
    catalog, stats = load_or_build_catalog(
        frontier_union, max_paths_per_pair=max_paths_per_pair, verbose=verbose)
    t_phase1 = time.time() - t_phase1
    print(f"  Vértices de fronteira do 6×6 (união 4 blocos): "
          f"{len(frontier_union)} de {V_BLOCK}")
    print(f"  Pares (s,e) com cores opostas catalogados: {stats['n_pairs_processed']}")
    print(f"  Pares vazios (sem HP): {stats['pairs_empty']}")
    print(f"  Pares saturados (atingiram cap={max_paths_per_pair}): {stats['pairs_full']}")
    print(f"  Total de paths catalogados: {stats['total_paths']}")
    print(f"  Tempo Fase 1: {t_phase1:.2f}s")
    print()

    print("=" * 70)
    print(f" FASE 2 — Matching perimetral (alvo K={K})")
    print("=" * 70)
    t_phase2 = time.time()
    tours, nodes = combine_perimeter(catalog, cross, max_tours=K, verbose=verbose)
    t_phase2 = time.time() - t_phase2
    print(f"  Tours encontrados: {len(tours)}")
    print(f"  Combinações testadas (esqueletos): {nodes}")
    print(f"  Tempo Fase 2: {t_phase2:.2f}s")

    # Validar
    valid = sum(1 for t in tours if verify_assembled_tour(t))
    print(f"  Tours válidos: {valid}/{len(tours)}")
    print()

    return {
        'tours': tours,
        'valid': valid,
        'nodes_dnc': nodes,
        't_phase1': t_phase1,
        't_phase2': t_phase2,
        'stats_catalog': stats,
    }


def benchmark_backtracking(K=50, seed=0):
    print("=" * 70)
    print(f" Backtracking direto 12×12 (kt.knight_tours, K={K})")
    print("=" * 70)
    t0 = time.time()
    tours = kt.knight_tours(N_GLOBAL, K, seed=seed)
    t1 = time.time()
    valid = sum(1 for t in tours if kt.verify_tour(t, N_GLOBAL))
    print(f"  Tours: {len(tours)} (válidos: {valid})")
    print(f"  Tempo: {t1 - t0:.2f}s")
    return {'tours': tours, 'time_s': t1 - t0}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--K', type=int, default=50)
    parser.add_argument('--max-paths-per-pair', type=int, default=200)
    parser.add_argument('--force-rebuild', action='store_true')
    parser.add_argument('--skip-bt', action='store_true',
                        help='Pular benchmark de backtracking')
    parser.add_argument('--quiet', action='store_true')
    args = parser.parse_args()

    if args.force_rebuild and os.path.exists(CATALOG_FILE):
        os.remove(CATALOG_FILE)

    verbose = not args.quiet

    res_dnc = benchmark_dnc(K=args.K, max_paths_per_pair=args.max_paths_per_pair,
                            verbose=verbose)

    if not args.skip_bt:
        res_bt = benchmark_backtracking(K=args.K, seed=0)
        print()
        print("=" * 70)
        print(" FASE 3 — Comparação")
        print("=" * 70)
        time_dnc_total = res_dnc['t_phase1'] + res_dnc['t_phase2']
        time_dnc_match = res_dnc['t_phase2']
        time_bt = res_bt['time_s']

        print(f"  | Método              | tempo total (s) | tours | nós explorados |")
        print(f"  |---------------------|-----------------|-------|----------------|")
        print(f"  | BT direto 12×12     | {time_bt:>15.2f} | {len(res_bt['tours']):>5} | {'(~5·K)':>14} |")
        print(f"  | D&C (precomp+match) | {time_dnc_total:>15.2f} | {len(res_dnc['tours']):>5} | "
              f"{res_dnc['nodes_dnc']:>14} |")
        print(f"  | D&C (só match)      | {time_dnc_match:>15.2f} | {len(res_dnc['tours']):>5} | "
              f"{res_dnc['nodes_dnc']:>14} |")
        print()
        if time_dnc_match > 0:
            speedup_amort = time_bt / time_dnc_match
            print(f"  Speedup match-only (catálogo amortizado): {speedup_amort:.2f}×")
        if time_dnc_total > 0:
            speedup_full = time_bt / time_dnc_total
            print(f"  Speedup full (cold-start): {speedup_full:.2f}×")
