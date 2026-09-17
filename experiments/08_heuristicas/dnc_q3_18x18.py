"""Q3 — Estender D&C para 18×18 com grade 3×3 de blocos 6×6.

Esquema:
  - 9 blocos, etiquetados como (br, bc) ∈ {0,1,2}×{0,1,2}
  - Cycle no meta-grafo de blocos: 9 blocos visitados em ordem cíclica
  - Reusa catálogo dnc_catalog_6x6.pkl (paths internos no 6×6 isolado)

Cycle escolhido (Hamiltoniano em meta-grafo king 3×3):
  (0,0) → (0,1) → (0,2) → (1,2) → (2,2) → (2,1) → (2,0) → (1,0) → (1,1) → (0,0)

Cada bloco tem (predecessor, sucessor) no ciclo; entry/exit vertices
devem estar nas fronteiras toward esses dois blocos.
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
from collections import defaultdict

import numpy as np

import knight_tours as kt
import knight_tours_dnc as dnc


N_BLOCK = 6
N_GLOBAL_18 = 18

# Blocos no 18×18, indexados por (br, bc)
BLOCKS_18 = {(br, bc): (br * N_BLOCK, bc * N_BLOCK)
             for br in range(3) for bc in range(3)}

KNIGHT_MOVES = dnc.KNIGHT_MOVES

# Cycle de 9 blocos (visitando todos)
CYCLE_18 = [(0, 0), (0, 1), (0, 2), (1, 2), (2, 2),
            (2, 1), (2, 0), (1, 0), (1, 1)]


def local_to_global_18(local_v, block):
    or_, oc_ = BLOCKS_18[block]
    lr, lc = divmod(local_v, N_BLOCK)
    return (or_ + lr) * N_GLOBAL_18 + (oc_ + lc)


def get_frontier_vertices_18(block):
    or_, oc_ = BLOCKS_18[block]
    fr = set()
    for r in range(N_BLOCK):
        for c in range(N_BLOCK):
            for dr, dc in KNIGHT_MOVES:
                gr = or_ + r + dr
                gc = oc_ + c + dc
                if not (0 <= gr < N_GLOBAL_18 and 0 <= gc < N_GLOBAL_18):
                    continue
                lr = gr - or_
                lc = gc - oc_
                if 0 <= lr < N_BLOCK and 0 <= lc < N_BLOCK:
                    continue
                fr.add(r * N_BLOCK + c)
                break
    return fr


def frontier_to(block_a, block_b):
    """Vértices locais de block_a com knight move para block_b."""
    or_a, oc_a = BLOCKS_18[block_a]
    or_b, oc_b = BLOCKS_18[block_b]
    fr = set()
    for r in range(N_BLOCK):
        for c in range(N_BLOCK):
            gra = or_a + r
            gca = oc_a + c
            for dr, dc in KNIGHT_MOVES:
                grb = gra + dr
                gcb = gca + dc
                if not (0 <= grb < N_GLOBAL_18 and 0 <= gcb < N_GLOBAL_18):
                    continue
                if or_b <= grb < or_b + N_BLOCK and oc_b <= gcb < oc_b + N_BLOCK:
                    fr.add(r * N_BLOCK + c)
                    break
    return fr


def cross_edges_18(block_a, block_b):
    """Knight edges globais entre block_a e block_b. Cada par (g_a, g_b)."""
    or_a, oc_a = BLOCKS_18[block_a]
    or_b, oc_b = BLOCKS_18[block_b]
    edges = []
    for r in range(N_BLOCK):
        for c in range(N_BLOCK):
            gra = or_a + r
            gca = oc_a + c
            for dr, dc in KNIGHT_MOVES:
                grb = gra + dr
                gcb = gca + dc
                if not (0 <= grb < N_GLOBAL_18 and 0 <= gcb < N_GLOBAL_18):
                    continue
                if or_b <= grb < or_b + N_BLOCK and oc_b <= gcb < oc_b + N_BLOCK:
                    edges.append((gra * N_GLOBAL_18 + gca, grb * N_GLOBAL_18 + gcb))
    return edges


def global_to_local_18(global_v, block):
    or_, oc_ = BLOCKS_18[block]
    gr, gc = divmod(global_v, N_GLOBAL_18)
    return (gr - or_) * N_BLOCK + (gc - oc_)


# ---------------------------------------------------------------------------
# Matching no 18×18
# ---------------------------------------------------------------------------

def random_matching_18(catalog, K, rng_seed=0, max_attempts=None,
                       verbose=False):
    """Tenta K tours no 18×18 via matching aleatório seguindo CYCLE_18.

    Para cada bloco A_i no ciclo:
      s_i ∈ frontier_to(A_i, A_{i-1})
      e_i ∈ frontier_to(A_i, A_{i+1})

    Conexão: e_i → s_{i+1} via knight edge na cross[A_i, A_{i+1}].
    Fecha: e_9 → s_0.
    """
    rng = np.random.default_rng(rng_seed)
    n_blocks = len(CYCLE_18)

    # Pré-computar fronteiras direcionais e cross-indices
    frontier_prev = {}  # block -> frontier vertices toward predecessor
    frontier_next = {}  # block -> frontier vertices toward successor
    cross_idx_next = {}  # block_i -> dict[v_local_e_i] -> [v_local_s_{i+1}]
    for i, b in enumerate(CYCLE_18):
        b_prev = CYCLE_18[(i - 1) % n_blocks]
        b_next = CYCLE_18[(i + 1) % n_blocks]
        frontier_prev[b] = frontier_to(b, b_prev)
        frontier_next[b] = frontier_to(b, b_next)
        # cross[A_i, A_{i+1}]: indexar por v_local em A_i
        idx = defaultdict(list)
        for (ga, gb) in cross_edges_18(b, b_next):
            la = global_to_local_18(ga, b)
            lb = global_to_local_18(gb, b_next)
            idx[la].append(lb)
        cross_idx_next[b] = idx

    if verbose:
        print(f"  Fronteiras direcionais (sample):")
        for i, b in enumerate(CYCLE_18):
            b_prev = CYCLE_18[(i - 1) % n_blocks]
            b_next = CYCLE_18[(i + 1) % n_blocks]
            print(f"    {b}: prev={b_prev}, next={b_next}, "
                  f"|F_prev|={len(frontier_prev[b])}, |F_next|={len(frontier_next[b])}")

    # Pré-computar pares (s, e) por bloco filtrados pelas fronteiras direcionais
    pairs_by_block = {}
    by_start_by_block = {}
    for b in CYCLE_18:
        fr_p = frontier_prev[b]
        fr_n = frontier_next[b]
        pairs = [(s, e) for (s, e) in catalog
                 if s in fr_p and e in fr_n and catalog[(s, e)]]
        # Filtrar: e deve ter pelo menos um cross para o próximo bloco
        idx_next = cross_idx_next[b]
        pairs = [(s, e) for (s, e) in pairs if e in idx_next]
        pairs_by_block[b] = pairs
        bs = defaultdict(list)
        for (s, e) in pairs:
            bs[s].append(e)
        by_start_by_block[b] = bs

    if verbose:
        for b in CYCLE_18:
            print(f"  {b}: {len(pairs_by_block[b])} pares (s,e) válidos")

    tours = []
    attempts = 0
    construct_fail = 0
    verify_fail = 0
    if max_attempts is None:
        max_attempts = K * 1000

    while len(tours) < K and attempts < max_attempts:
        attempts += 1
        # Amostrar primeiro bloco
        b0 = CYCLE_18[0]
        if not pairs_by_block[b0]:
            break
        s0, e0 = pairs_by_block[b0][rng.integers(len(pairs_by_block[b0]))]

        # Propagar: blocks_paths[(b, s, e)]; armazenar sequência
        endpoints = [(s0, e0)]  # endpoint do bloco 0
        feasible = True
        for i in range(1, n_blocks):
            b_prev = CYCLE_18[i - 1]
            b_curr = CYCLE_18[i]
            e_prev = endpoints[-1][1]
            # s_curr ∈ knight_dest(e_prev) ∩ frontier_prev(b_curr)
            cand_s = [v for v in cross_idx_next[b_prev][e_prev]
                      if v in frontier_prev[b_curr]]
            if not cand_s:
                construct_fail += 1
                feasible = False
                break
            s_curr = cand_s[rng.integers(len(cand_s))]
            # e_curr ∈ by_start_by_block[b_curr][s_curr]
            if i == n_blocks - 1:
                # Último bloco: e_curr deve fechar para s0 via cross
                cand_e = [
                    e for e in by_start_by_block[b_curr].get(s_curr, [])
                    if s0 in cross_idx_next[b_curr][e]
                ]
            else:
                cand_e = by_start_by_block[b_curr].get(s_curr, [])
            if not cand_e:
                construct_fail += 1
                feasible = False
                break
            e_curr = cand_e[rng.integers(len(cand_e))]
            endpoints.append((s_curr, e_curr))

        if not feasible:
            continue

        # Montar tour
        block_paths = []
        for i, b in enumerate(CYCLE_18):
            s, e = endpoints[i]
            paths = catalog[(s, e)]
            p = paths[rng.integers(len(paths))]
            block_paths.append((b, p))

        seq = []
        for (b, p) in block_paths:
            for v in p:
                seq.append(local_to_global_18(v, b))
        tour = np.asarray(seq, dtype=np.int32)

        if kt.verify_tour(tour, N_GLOBAL_18):
            tours.append(tour)
        else:
            verify_fail += 1

    return {
        'tours': tours,
        'attempts': attempts,
        'construct_fail': construct_fail,
        'verify_fail': verify_fail,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    print("=" * 70)
    print(" Q3 — D&C 18×18 (grade 3×3 de blocos 6×6)")
    print("=" * 70)

    # Carrega catálogo do 12×12
    if not os.path.exists(dnc.CATALOG_FILE):
        raise SystemExit("Catálogo não existe. Rode knight_tours_dnc.py primeiro.")
    with open(dnc.CATALOG_FILE, 'rb') as f:
        obj = pickle.load(f)
    catalog = obj['catalog']
    print(f"  Catálogo carregado: {len(catalog)} pares, "
          f"{sum(len(v) for v in catalog.values())} paths")

    # Diagnóstico das fronteiras direcionais
    print("\n  Cycle:")
    for b in CYCLE_18:
        print(f"    {b}")

    # Diagnóstico do alcance do catálogo:
    # quais fronteiras direcionais existem
    print("\n  Fronteiras direcionais por bloco (no cycle):")
    for i, b in enumerate(CYCLE_18):
        b_prev = CYCLE_18[(i - 1) % len(CYCLE_18)]
        b_next = CYCLE_18[(i + 1) % len(CYCLE_18)]
        fp = frontier_to(b, b_prev)
        fn = frontier_to(b, b_next)
        print(f"    {b} ← {b_prev}: |F|={len(fp)}; → {b_next}: |F|={len(fn)}")

    # Diagnóstico de catálogo: alguma fronteira não está nele?
    cat_vertices = set()
    for (s, e) in catalog.keys():
        cat_vertices.add(s)
        cat_vertices.add(e)
    print(f"\n  Vértices presentes no catálogo: {len(cat_vertices)}")
    needed_vertices = set()
    for i, b in enumerate(CYCLE_18):
        b_prev = CYCLE_18[(i - 1) % len(CYCLE_18)]
        b_next = CYCLE_18[(i + 1) % len(CYCLE_18)]
        needed_vertices.update(frontier_to(b, b_prev))
        needed_vertices.update(frontier_to(b, b_next))
    print(f"  Vértices necessários: {len(needed_vertices)}")
    missing = needed_vertices - cat_vertices
    print(f"  Faltando no catálogo: {len(missing)} {sorted(missing) if missing else ''}")

    # Tentar matching
    K = 20
    print(f"\n  Tentando K={K} tours no 18×18 (D&C, timeout 60s) ...")
    t0 = time.time()
    res = random_matching_18(catalog, K=K, rng_seed=0,
                             max_attempts=200000, verbose=True)
    t1 = time.time()
    print(f"\n  D&C 18×18: tours={len(res['tours'])}/{K}, "
          f"tentativas={res['attempts']}, "
          f"construct_fail={res['construct_fail']}, "
          f"verify_fail={res['verify_fail']} | {t1 - t0:.2f}s")
    if res['tours']:
        valid = sum(1 for t in res['tours'] if kt.verify_tour(t, N_GLOBAL_18))
        print(f"  Tours válidos: {valid}/{len(res['tours'])}")

    # Comparar com BT direto 18×18
    print(f"\n  Backtracking direto 18×18 (kt.knight_tours, K={K}, timeout 60s)...")
    import signal

    class TimeoutError(Exception):
        pass

    def _to_handler(signum, frame):
        raise TimeoutError()

    signal.signal(signal.SIGALRM, _to_handler)
    signal.alarm(60)
    t0 = time.time()
    try:
        tours_bt = kt.knight_tours(N_GLOBAL_18, K, seed=0)
        t1 = time.time()
        signal.alarm(0)
        print(f"  BT 18×18: tours={len(tours_bt)} em {t1 - t0:.2f}s")
        valid = sum(1 for t in tours_bt if kt.verify_tour(t, N_GLOBAL_18))
        print(f"  Tours válidos: {valid}/{len(tours_bt)}")
    except TimeoutError:
        t1 = time.time()
        signal.alarm(0)
        print(f"  BT 18×18: TIMEOUT após {t1 - t0:.2f}s")

    print()
    print("=" * 70)
    print(" Reusabilidade do catálogo")
    print("=" * 70)
    print("  O catálogo é UNIVERSAL: paths são no 6×6 isolado, idênticos")
    print("  independente da posição do bloco no tabuleiro maior.")
    print("  O que muda por posição é o FILTRO de fronteira (quais (s,e)")
    print("  são utilizáveis para esse bloco-no-cycle).")
