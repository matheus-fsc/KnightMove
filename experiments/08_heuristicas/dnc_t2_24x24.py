"""T2 — D&C para 24×24 (grade 4×4 de blocos 6×6).

Hamilton cycle no meta-grafo king 4×4 (16 nós, 16 arestas).

Cycle escolhido:
  (0,0)→(1,0)→(2,0)→(3,0)→(3,1)→(3,2)→(3,3)→(2,3)
       →(2,2)→(1,3)→(0,3)→(0,2)→(1,2)→(2,1)→(1,1)→(0,1)→(0,0)

Verificação de adjacências (king: |dr|,|dc| ≤ 1, não ambos zero):
  (0,0)-(1,0) adj
  (1,0)-(2,0) adj
  (2,0)-(3,0) adj
  (3,0)-(3,1) adj
  (3,1)-(3,2) adj
  (3,2)-(3,3) adj
  (3,3)-(2,3) adj
  (2,3)-(2,2) adj
  (2,2)-(1,3) diag
  (1,3)-(0,3) adj
  (0,3)-(0,2) adj
  (0,2)-(1,2) adj
  (1,2)-(2,1) diag
  (2,1)-(1,1) adj
  (1,1)-(0,1) adj
  (0,1)-(0,0) adj
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
import signal
from collections import defaultdict

import numpy as np

import knight_tours as kt
import knight_tours_dnc as dnc


N_BLOCK = 6
N_GLOBAL_24 = 24

BLOCKS_24 = {(br, bc): (br * N_BLOCK, bc * N_BLOCK)
             for br in range(4) for bc in range(4)}

CYCLE_24 = [(0, 0), (1, 0), (2, 0), (3, 0),
            (3, 1), (3, 2), (3, 3), (2, 3),
            (2, 2), (1, 3), (0, 3), (0, 2),
            (1, 2), (2, 1), (1, 1), (0, 1)]


def local_to_global_24(local_v, block):
    or_, oc_ = BLOCKS_24[block]
    lr, lc = divmod(local_v, N_BLOCK)
    return (or_ + lr) * N_GLOBAL_24 + (oc_ + lc)


def global_to_local_24(global_v, block):
    or_, oc_ = BLOCKS_24[block]
    gr, gc = divmod(global_v, N_GLOBAL_24)
    return (gr - or_) * N_BLOCK + (gc - oc_)


def frontier_to_24(block_a, block_b):
    or_a, oc_a = BLOCKS_24[block_a]
    or_b, oc_b = BLOCKS_24[block_b]
    fr = set()
    for r in range(N_BLOCK):
        for c in range(N_BLOCK):
            gra = or_a + r
            gca = oc_a + c
            for dr, dc in dnc.KNIGHT_MOVES:
                grb = gra + dr
                gcb = gca + dc
                if not (0 <= grb < N_GLOBAL_24 and 0 <= gcb < N_GLOBAL_24):
                    continue
                if or_b <= grb < or_b + N_BLOCK and oc_b <= gcb < oc_b + N_BLOCK:
                    fr.add(r * N_BLOCK + c)
                    break
    return fr


def cross_edges_24(block_a, block_b):
    or_a, oc_a = BLOCKS_24[block_a]
    or_b, oc_b = BLOCKS_24[block_b]
    edges = []
    for r in range(N_BLOCK):
        for c in range(N_BLOCK):
            gra = or_a + r
            gca = oc_a + c
            for dr, dc in dnc.KNIGHT_MOVES:
                grb = gra + dr
                gcb = gca + dc
                if not (0 <= grb < N_GLOBAL_24 and 0 <= gcb < N_GLOBAL_24):
                    continue
                if or_b <= grb < or_b + N_BLOCK and oc_b <= gcb < oc_b + N_BLOCK:
                    edges.append((gra * N_GLOBAL_24 + gca,
                                  grb * N_GLOBAL_24 + gcb))
    return edges


def random_matching_24(catalog, K, rng_seed=0, max_attempts=None,
                       verbose=False):
    rng = np.random.default_rng(rng_seed)
    n_blocks = len(CYCLE_24)

    frontier_prev = {}
    frontier_next = {}
    cross_idx_next = {}
    for i, b in enumerate(CYCLE_24):
        b_prev = CYCLE_24[(i - 1) % n_blocks]
        b_next = CYCLE_24[(i + 1) % n_blocks]
        frontier_prev[b] = frontier_to_24(b, b_prev)
        frontier_next[b] = frontier_to_24(b, b_next)
        idx = defaultdict(list)
        for (ga, gb) in cross_edges_24(b, b_next):
            la = global_to_local_24(ga, b)
            lb = global_to_local_24(gb, b_next)
            idx[la].append(lb)
        cross_idx_next[b] = idx

    pairs_by_block = {}
    by_start_by_block = {}
    for b in CYCLE_24:
        fr_p = frontier_prev[b]
        fr_n = frontier_next[b]
        idx_next = cross_idx_next[b]
        pairs = [(s, e) for (s, e) in catalog
                 if s in fr_p and e in fr_n and catalog[(s, e)]
                 and e in idx_next]
        pairs_by_block[b] = pairs
        bs = defaultdict(list)
        for (s, e) in pairs:
            bs[s].append(e)
        by_start_by_block[b] = bs

    if verbose:
        print(f"  Pares válidos por bloco no cycle 24×24:")
        for i, b in enumerate(CYCLE_24):
            b_prev = CYCLE_24[(i - 1) % n_blocks]
            b_next = CYCLE_24[(i + 1) % n_blocks]
            print(f"    {b} (prev={b_prev}, next={b_next}): "
                  f"|F_prev|={len(frontier_prev[b])}, "
                  f"|F_next|={len(frontier_next[b])}, "
                  f"|pares|={len(pairs_by_block[b])}")

    tours = []
    attempts = 0
    construct_fail = 0
    verify_fail = 0
    if max_attempts is None:
        max_attempts = K * 2000

    while len(tours) < K and attempts < max_attempts:
        attempts += 1
        b0 = CYCLE_24[0]
        if not pairs_by_block[b0]:
            construct_fail += 1
            break
        s0, e0 = pairs_by_block[b0][rng.integers(len(pairs_by_block[b0]))]
        endpoints = [(s0, e0)]
        feasible = True
        for i in range(1, n_blocks):
            b_prev = CYCLE_24[i - 1]
            b_curr = CYCLE_24[i]
            e_prev = endpoints[-1][1]
            cand_s = [v for v in cross_idx_next[b_prev][e_prev]
                      if v in frontier_prev[b_curr]]
            if not cand_s:
                construct_fail += 1
                feasible = False
                break
            s_curr = cand_s[rng.integers(len(cand_s))]
            if i == n_blocks - 1:
                cand_e = [e for e in by_start_by_block[b_curr].get(s_curr, [])
                          if s0 in cross_idx_next[b_curr].get(e, [])]
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

        seq = []
        for i, b in enumerate(CYCLE_24):
            s, e = endpoints[i]
            p = catalog[(s, e)][rng.integers(len(catalog[(s, e)]))]
            for v in p:
                seq.append(local_to_global_24(v, b))
        tour = np.asarray(seq, dtype=np.int32)
        if kt.verify_tour(tour, N_GLOBAL_24):
            tours.append(tour)
        else:
            verify_fail += 1

    return {
        'tours': tours,
        'attempts': attempts,
        'construct_fail': construct_fail,
        'verify_fail': verify_fail,
    }


def edge_set(tour):
    V = len(tour)
    return frozenset(frozenset((int(tour[i]), int(tour[(i + 1) % V])))
                     for i in range(V))


if __name__ == '__main__':
    print("=" * 70)
    print(" T2 — D&C 24×24 (grade 4×4 = 16 blocos 6×6)")
    print("=" * 70)

    with open(dnc.CATALOG_FILE, 'rb') as f:
        obj = pickle.load(f)
    catalog = obj['catalog']
    print(f"  Catálogo reusado: {sum(len(v) for v in catalog.values())} paths")

    # Validar Hamilton cycle
    print(f"\n  Hamilton cycle no king 4×4 (16 blocos):")
    for i, b in enumerate(CYCLE_24):
        b_next = CYCLE_24[(i + 1) % len(CYCLE_24)]
        dr = abs(b[0] - b_next[0])
        dc = abs(b[1] - b_next[1])
        adj = (dr <= 1 and dc <= 1 and (dr + dc) > 0)
        print(f"    {b} → {b_next}  king_adj={adj}")

    K = 20
    print(f"\n  Tentando K={K} tours via D&C 24×24...")
    t0 = time.time()
    res = random_matching_24(catalog, K=K, rng_seed=0, verbose=True)
    t1 = time.time()
    print(f"\n  tours={len(res['tours'])}/{K}, "
          f"attempts={res['attempts']}, "
          f"construct_fail={res['construct_fail']}, "
          f"verify_fail={res['verify_fail']} | {t1 - t0:.3f}s")
    if res['tours']:
        valid = sum(1 for t in res['tours'] if kt.verify_tour(t, N_GLOBAL_24))
        hashes = set(hash(edge_set(t)) for t in res['tours'])
        print(f"  Tours válidos: {valid}/{len(res['tours'])}")
        print(f"  Tours únicos: {len(hashes)}/{len(res['tours'])}")

    # BT direto 24×24 com timeout 120s
    print(f"\n  Backtracking direto 24×24 (kt.knight_tours, K={K}, timeout 120s)...")

    class TimeoutErr(Exception):
        pass

    def _handler(signum, frame):
        raise TimeoutErr()

    signal.signal(signal.SIGALRM, _handler)
    signal.alarm(120)
    t0 = time.time()
    try:
        tours_bt = kt.knight_tours(N_GLOBAL_24, K, seed=0)
        t1 = time.time()
        signal.alarm(0)
        valid_bt = sum(1 for t in tours_bt if kt.verify_tour(t, N_GLOBAL_24))
        print(f"  BT 24×24: tours={len(tours_bt)} válidos={valid_bt} em {t1-t0:.2f}s")
    except TimeoutErr:
        t1 = time.time()
        signal.alarm(0)
        print(f"  BT 24×24: TIMEOUT após {t1-t0:.2f}s (0 tours)")
