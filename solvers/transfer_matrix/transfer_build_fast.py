"""
Versao otimizada do build de T_bulk para n=8 e n=10.

Otimizacoes:
  1. Estado codificado como int (base-3 para deg_c, base-2 (wait base-3) para deg_cp1).
  2. Pre-computa, para cada deg_c, a lista de tuplas (cp1_inc, cp2_inc)
     possiveis em E_c -- depois para cada estado iteramos apenas filtrando
     deg_cp1 + cp1_inc <= 2.
  3. Usa listas de inteiros em vez de dicts onde possivel.

State id = deg_c_int * BASE_CP1 + deg_cp1_int
  deg_c_int  in [0, 3^n)
  deg_cp1_int in [0, 3^n)
  total raw state id in [0, 9^n)
"""

import json
import sys
import time
from itertools import combinations, product
from pathlib import Path

import numpy as np
import scipy.sparse as sp

THIS = Path(__file__).parent
DATA = THIS / "data"
DATA.mkdir(exist_ok=True)


def forward_options(n: int, r: int):
    opts = []
    for dr in (-2, 2):
        if 0 <= r + dr < n:
            opts.append((1, r + dr))
    for dr in (-1, 1):
        if 0 <= r + dr < n:
            opts.append((2, r + dr))
    return opts


# ---------------------------------------------------------------------------
# Encoding / decoding state
# ---------------------------------------------------------------------------

def to_int_base3(arr, n):
    v = 0
    for x in arr:
        v = v * 3 + int(x)
    return v


def from_int_base3(v, n):
    out = [0] * n
    for i in range(n - 1, -1, -1):
        out[i] = v % 3
        v //= 3
    return tuple(out)


def state_to_int(state, n, BASE):
    deg_c, deg_cp1 = state
    return to_int_base3(deg_c, n) * BASE + to_int_base3(deg_cp1, n)


def int_to_state(idx, n, BASE):
    deg_c_int = idx // BASE
    deg_cp1_int = idx % BASE
    return (from_int_base3(deg_c_int, n), from_int_base3(deg_cp1_int, n))


# ---------------------------------------------------------------------------
# Precomputa transicoes por deg_c
# ---------------------------------------------------------------------------

def precompute_transitions_by_deg_c(n, forbid_cp2=False):
    """Para cada deg_c, gera lista de (cp1_inc_tuple, cp2_inc_tuple) usando
    BACKTRACKING (com pruning), evitando o blowup do produto cartesiano.

    Pruning:
      - durante a recursao por linha, cp1_inc[i] e cp2_inc[i] nao podem exceder 2
        (deg_cp1[i] ainda nao verificado -- so' verificamos cp1_inc, mais permissivo)
    """
    transitions = {}

    for deg_c_int in range(3 ** n):
        deg_c = from_int_base3(deg_c_int, n)
        # checar viabilidade
        row_opts = []
        valid_state = True
        for r in range(n):
            need = 2 - deg_c[r]
            if need < 0:
                valid_state = False
                break
            opts = forward_options(n, r)
            if forbid_cp2:
                opts = [o for o in opts if o[0] == 1]
            if need > len(opts):
                valid_state = False
                break
            row_opts.append((need, opts))
        if not valid_state:
            continue

        outs = []
        cp1_inc = [0] * n
        cp2_inc = [0] * n

        def rec(r):
            if r == n:
                outs.append((tuple(cp1_inc), tuple(cp2_inc)))
                return
            need, opts = row_opts[r]
            if need == 0:
                rec(r + 1)
                return
            for combo in combinations(range(len(opts)), need):
                back1 = []
                back2 = []
                ok = True
                for idx in combo:
                    co, dst = opts[idx]
                    if co == 1:
                        cp1_inc[dst] += 1
                        if cp1_inc[dst] > 2:
                            ok = False
                            cp1_inc[dst] -= 1
                            break
                        back1.append(dst)
                    else:
                        cp2_inc[dst] += 1
                        if cp2_inc[dst] > 2:
                            ok = False
                            cp2_inc[dst] -= 1
                            break
                        back2.append(dst)
                if ok:
                    rec(r + 1)
                for d in back1:
                    cp1_inc[d] -= 1
                for d in back2:
                    cp2_inc[d] -= 1

        rec(0)
        if outs:
            transitions[deg_c_int] = outs
    return transitions


# ---------------------------------------------------------------------------
# BFS de estados alcancaveis usando precomputado
# ---------------------------------------------------------------------------

def build_T_fast(n, forbid_cp2=False, verbose=True):
    BASE = 3 ** n

    if verbose:
        print(f"  precomputando transicoes por deg_c (3^{n}={3**n} entries)...")
    t0 = time.time()
    trans_by_deg_c = precompute_transitions_by_deg_c(n, forbid_cp2=forbid_cp2)
    t1 = time.time()
    if verbose:
        print(f"    {len(trans_by_deg_c)} deg_c validos, "
              f"avg {sum(len(v) for v in trans_by_deg_c.values()) / max(1, len(trans_by_deg_c)):.1f} trans/deg_c, "
              f"tempo {t1 - t0:.2f}s")

    if verbose:
        print(f"  BFS dos estados alcancaveis...")

    # Estado inicial: deg_c = deg_cp1 = 0
    s0 = (0, 0)  # (deg_c_int=0, deg_cp1_int=0)
    s0_id = 0 * BASE + 0
    idx = {s0_id: 0}
    queue = [s0_id]
    rows, cols, data = [], [], []

    head = 0
    last_report = time.time()
    while head < len(queue):
        sid = queue[head]
        head += 1
        i = idx[sid]
        deg_c_int = sid // BASE
        deg_cp1_int = sid % BASE
        deg_cp1 = from_int_base3(deg_cp1_int, n)

        trans = trans_by_deg_c.get(deg_c_int)
        if not trans:
            continue

        for cp1_inc, cp2_inc in trans:
            # check deg_cp1[r] + cp1_inc[r] <= 2
            new_cp1 = [0] * n
            valid = True
            for r in range(n):
                v = deg_cp1[r] + cp1_inc[r]
                if v > 2:
                    valid = False
                    break
                new_cp1[r] = v
            if not valid:
                continue
            new_cp1_int = to_int_base3(new_cp1, n)
            new_cp2_int = to_int_base3(cp2_inc, n)
            new_sid = new_cp1_int * BASE + new_cp2_int

            j = idx.get(new_sid)
            if j is None:
                j = len(idx)
                idx[new_sid] = j
                queue.append(new_sid)
            rows.append(i)
            cols.append(j)
            data.append(1)

        if verbose and time.time() - last_report > 30:
            print(f"    BFS: |idx|={len(idx)}, head={head}, queue_len={len(queue)}, "
                  f"transicoes acum={len(rows)}")
            last_report = time.time()

    K = len(idx)
    if verbose:
        print(f"  K = {K} estados alcancaveis")
        print(f"  nnz total = {len(rows)}")

    T = sp.csr_matrix(
        (data, (rows, cols)), shape=(K, K), dtype=np.int64
    )
    return T, idx, trans_by_deg_c


# ---------------------------------------------------------------------------
# Contagem via DP de produto matricial
# ---------------------------------------------------------------------------

def count_2factors_via_T(n, T_bulk, T_nocp2, idx_ext, BASE):
    """count = <s_0 | T_bulk^{n-2} . T_nocp2 | s_target>
    onde s_target = (deg_c=(2,...,2), deg_cp1=(0,...,0)).
    """
    K = T_bulk.shape[0]
    s0_id = 0
    s_target_int = to_int_base3([2] * n, n) * BASE + 0
    if s_target_int not in idx_ext:
        return None
    target = idx_ext[s_target_int]

    v = np.zeros(K, dtype=np.float64)
    v[idx_ext[s0_id]] = 1.0
    for _ in range(n - 2):
        v = T_bulk.T @ v
    v = T_nocp2.T @ v
    return int(round(v[target]))


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def run_for(n: int, save: bool = True):
    print(f"\n=== n = {n} ===")

    BASE = 3 ** n

    t0 = time.time()
    T_bulk, idx_bulk, _ = build_T_fast(n, forbid_cp2=False)
    t1 = time.time()
    print(f"  build T_bulk: {t1 - t0:.2f}s")

    # T_nocp2 sobre o mesmo idx + extensoes
    t2 = time.time()
    print(f"  precomputando trans nocp2...")
    trans_nocp2_by_deg_c = precompute_transitions_by_deg_c(n, forbid_cp2=True)
    idx_ext = dict(idx_bulk)
    rows, cols, data = [], [], []
    for sid, i in list(idx_bulk.items()):
        deg_c_int = sid // BASE
        deg_cp1_int = sid % BASE
        deg_cp1 = from_int_base3(deg_cp1_int, n)
        trans = trans_nocp2_by_deg_c.get(deg_c_int)
        if not trans:
            continue
        for cp1_inc, cp2_inc in trans:
            new_cp1 = [0] * n
            valid = True
            for r in range(n):
                v = deg_cp1[r] + cp1_inc[r]
                if v > 2:
                    valid = False
                    break
                new_cp1[r] = v
            if not valid:
                continue
            new_cp1_int = to_int_base3(new_cp1, n)
            new_cp2_int = to_int_base3(cp2_inc, n)
            new_sid = new_cp1_int * BASE + new_cp2_int
            j = idx_ext.get(new_sid)
            if j is None:
                j = len(idx_ext)
                idx_ext[new_sid] = j
            rows.append(i)
            cols.append(j)
            data.append(1)
    K_ext = len(idx_ext)
    T_nocp2 = sp.csr_matrix(
        (data, (rows, cols)), shape=(K_ext, K_ext), dtype=np.int64
    )
    t3 = time.time()
    print(f"  build T_nocp2: {t3 - t2:.2f}s")
    print(f"  K_ext={K_ext}, nnz_nocp2={T_nocp2.nnz}")

    # Re-shape T_bulk para idx_ext (com zeros nas linhas/colunas extras)
    if K_ext > T_bulk.shape[0]:
        # extend by zero rows/cols
        T_bulk_full = sp.lil_matrix((K_ext, K_ext), dtype=np.int64)
        T_coo = T_bulk.tocoo()
        for r, c, d in zip(T_coo.row, T_coo.col, T_coo.data):
            T_bulk_full[r, c] = d
        T_bulk_full = T_bulk_full.tocsr()
    else:
        T_bulk_full = T_bulk

    print(f"  contagem via produto matricial...")
    t4 = time.time()
    cnt = count_2factors_via_T(n, T_bulk_full, T_nocp2, idx_ext, BASE)
    t5 = time.time()
    print(f"  #2-fatores(n={n}) = {cnt}  (tempo {t5 - t4:.2f}s)")

    if save:
        sp.save_npz(DATA / f"transfer_n{n}_bulk.npz", T_bulk)
        sp.save_npz(DATA / f"transfer_n{n}_nocp2.npz", T_nocp2)
        # Salvar idx mapping (compatibilidade com analise legada: precisamos
        # dos estados como tuples para load_idx em spectral_analysis e
        # compare_finf)
        idx_tuple_repr = {}
        for sid, i in idx_ext.items():
            deg_c_int = sid // BASE
            deg_cp1_int = sid % BASE
            deg_c = from_int_base3(deg_c_int, n)
            deg_cp1 = from_int_base3(deg_cp1_int, n)
            key = f"{deg_c}|{deg_cp1}"
            idx_tuple_repr[key] = i
        with open(DATA / f"transfer_n{n}_idx.json", "w") as f:
            json.dump(idx_tuple_repr, f)
        print(f"  salvo em data/transfer_n{n}_*.{{npz,json}}")

    info = {
        "n": n,
        "K": len(idx_bulk),
        "K_ext": K_ext,
        "nnz_bulk": int(T_bulk.nnz),
        "nnz_nocp2": int(T_nocp2.nnz),
        "count_dp": None,
        "count_matrix": cnt,
        "time_build_s": t1 - t0,
        "time_nocp2_s": t3 - t2,
    }
    sc_path = DATA / "state_count.json"
    existing = []
    if sc_path.exists():
        with open(sc_path) as f:
            existing = json.load(f)
    by_n = {r["n"]: r for r in existing}
    by_n[n] = info
    with open(sc_path, "w") as f:
        json.dump(sorted(by_n.values(), key=lambda r: r["n"]), f, indent=2)

    return info


if __name__ == "__main__":
    ns = [int(x) for x in sys.argv[1:]] if len(sys.argv) > 1 else [6]
    for n in ns:
        run_for(n)
