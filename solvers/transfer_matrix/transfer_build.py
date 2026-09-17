"""
Constroi a transfer matrix (column-sweep) do passeio do cavalo n×n
para contagem de 2-fatores.

Estado: s = (deg_c, deg_{c+1}) com deg_c[r], deg_{c+1}[r] em {0,1,2}.
Passo c: decide arestas com extremo esquerdo em col c:
  - (c, c+1) via knight (col_diff=1, row_diff=±2)
  - (c, c+2) via knight (col_diff=2, row_diff=±1)

Saidas:
  - data/transfer_n{n}.npz  (T_bulk esparsa + state index)
  - data/state_count.json  (K(n) e nnz por n)
  - log no stdout
"""

import json
import time
from itertools import combinations
from pathlib import Path

import numpy as np
import scipy.sparse as sp

DATA = Path(__file__).parent / "data"
DATA.mkdir(exist_ok=True)


# ---------------------------------------------------------------------------
# Geometria das arestas forward
# ---------------------------------------------------------------------------

def forward_options(n: int, r: int):
    """Arestas forward saindo de (r, c). Retorna lista de (col_off, row_dst)."""
    opts = []
    for dr in (-2, 2):
        if 0 <= r + dr < n:
            opts.append((1, r + dr))  # para col c+1, row r+dr
    for dr in (-1, 1):
        if 0 <= r + dr < n:
            opts.append((2, r + dr))  # para col c+2, row r+dr
    return opts


# ---------------------------------------------------------------------------
# Enumeracao de transicoes a partir de um estado
# ---------------------------------------------------------------------------

def transitions(state, n: int, forbid_cp2: bool = False):
    """Itera (s_next,) para cada escolha valida de E_c.

    Multiplicidade e sempre 1: cada conjunto de arestas E_c distinto produz
    exatamente uma transicao.
    """
    deg_c, deg_cp1 = state

    # Por linha, precisamos escolher exatamente 'need' arestas forward
    row_choices = []
    for r in range(n):
        need = 2 - deg_c[r]
        if need < 0:
            return
        opts = forward_options(n, r)
        if forbid_cp2:
            opts = [o for o in opts if o[0] == 1]
        if need > len(opts):
            return  # impossivel fechar col c
        row_choices.append((need, opts))

    # Backtracking por linha
    def rec(r, cp1_inc, cp2_inc):
        if r == n:
            new_cp1 = tuple(deg_cp1[i] + cp1_inc[i] for i in range(n))
            new_cp2 = tuple(cp2_inc[i] for i in range(n))
            yield (new_cp1, new_cp2)
            return

        need, opts = row_choices[r]
        if need == 0:
            yield from rec(r + 1, cp1_inc, cp2_inc)
            return

        for combo in combinations(range(len(opts)), need):
            # tentar adicionar combo
            cp1_back = []
            cp2_back = []
            valid = True
            for idx in combo:
                col_off, dst_r = opts[idx]
                if col_off == 1:
                    if deg_cp1[dst_r] + cp1_inc[dst_r] + 1 > 2:
                        valid = False
                        break
                    cp1_inc[dst_r] += 1
                    cp1_back.append(dst_r)
                else:
                    if cp2_inc[dst_r] + 1 > 2:
                        valid = False
                        break
                    cp2_inc[dst_r] += 1
                    cp2_back.append(dst_r)
            if valid:
                yield from rec(r + 1, cp1_inc, cp2_inc)
            # desfazer
            for dst_r in cp1_back:
                cp1_inc[dst_r] -= 1
            for dst_r in cp2_back:
                cp2_inc[dst_r] -= 1

    yield from rec(0, [0] * n, [0] * n)


# ---------------------------------------------------------------------------
# Construcao da matriz e BFS de estados alcancaveis
# ---------------------------------------------------------------------------

def build_T(n: int):
    """BFS a partir de s_0 = (0,0) usando transicoes BULK (forbid_cp2=False).

    Retorna (T_csr, idx) onde idx[state]=i e T_csr eh K x K.
    """
    s0 = ((0,) * n, (0,) * n)
    idx = {s0: 0}
    queue = [s0]
    rows, cols, data = [], [], []

    head = 0
    while head < len(queue):
        s = queue[head]
        head += 1
        i = idx[s]
        for s_next in transitions(s, n, forbid_cp2=False):
            j = idx.get(s_next)
            if j is None:
                j = len(idx)
                idx[s_next] = j
                queue.append(s_next)
            rows.append(i)
            cols.append(j)
            data.append(1)

    K = len(idx)
    T = sp.csr_matrix(
        (data, (rows, cols)), shape=(K, K), dtype=np.int64
    )
    return T, idx


def build_T_no_cp2(n: int, idx_known):
    """Constroi T_no_cp2 usando o mesmo indexamento (idx_known).

    Estados sucessores podem ser novos -- vamos estende-los em idx_known.
    Retorna (T_csr, idx_extended).
    """
    idx = dict(idx_known)
    rows, cols, data = [], [], []
    # Iterar sobre todos os estados ja conhecidos
    for s, i in list(idx.items()):
        for s_next in transitions(s, n, forbid_cp2=True):
            j = idx.get(s_next)
            if j is None:
                j = len(idx)
                idx[s_next] = j
            rows.append(i)
            cols.append(j)
            data.append(1)

    K = len(idx)
    T = sp.csr_matrix(
        (data, (rows, cols)), shape=(K, K), dtype=np.int64
    )
    return T, idx


# ---------------------------------------------------------------------------
# Contagem de 2-fatores via DP
# ---------------------------------------------------------------------------

def count_2factors_dp(n: int):
    """Conta 2-fatores via DP step-by-step sem materializar T (sanity)."""
    s0 = ((0,) * n, (0,) * n)
    v = {s0: 1}

    for c in range(n):
        last = (c == n - 1)
        second_last = (c == n - 2)
        forbid_cp2 = second_last or last

        new_v = {}
        if last:
            # Apenas E_c = vazio; requer deg_c = (2,...,2).
            for s, w in v.items():
                deg_c, deg_cp1 = s
                if all(d == 2 for d in deg_c):
                    s_next = (deg_cp1, (0,) * n)
                    new_v[s_next] = new_v.get(s_next, 0) + w
        else:
            for s, w in v.items():
                for s_next in transitions(s, n, forbid_cp2=forbid_cp2):
                    new_v[s_next] = new_v.get(s_next, 0) + w
        v = new_v

    return v.get(((0,) * n, (0,) * n), 0)


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def run_for(n: int, save: bool = True):
    print(f"\n=== n = {n} ===")

    t0 = time.time()
    T, idx = build_T(n)
    t1 = time.time()
    K = len(idx)
    nnz = T.nnz
    print(f"  K (estados alcancaveis a partir de s_0) = {K}")
    print(f"  nnz(T_bulk) = {nnz}")
    print(f"  tempo build T_bulk: {t1 - t0:.2f}s")

    t2 = time.time()
    T_nc2, idx_ext = build_T_no_cp2(n, idx)
    t3 = time.time()
    K_ext = len(idx_ext)
    print(f"  K_ext (incluindo successores no_cp2) = {K_ext}")
    print(f"  nnz(T_no_cp2) = {T_nc2.nnz}")
    print(f"  tempo build T_no_cp2: {t3 - t2:.2f}s")

    # Verificacao via DP (sanity)
    print(f"  computando contagem por DP...")
    t4 = time.time()
    if n <= 8:
        cnt = count_2factors_dp(n)
    else:
        cnt = None  # caro, pulamos
    t5 = time.time()
    if cnt is not None:
        print(f"  #2-fatores(n={n}) = {cnt}")
        print(f"  tempo DP: {t5 - t4:.2f}s")

    # Verificacao via T (multiplicacao vetor-matriz)
    # v_{n-1} = e_{s_0} * T_bulk^{n-2} * T_no_cp2
    # Alvo: estado ((2,...,2), (0,...,0))
    print(f"  verificando via produto matricial...")
    # Re-indexar T_bulk no espaco idx_ext (estendido)
    if K_ext > K:
        # remap T_bulk para idx_ext
        T_bulk_ext = sp.lil_matrix((K_ext, K_ext), dtype=np.int64)
        T_coo = T.tocoo()
        for r, c_, v_ in zip(T_coo.row, T_coo.col, T_coo.data):
            T_bulk_ext[r, c_] = v_
        T_bulk_ext = T_bulk_ext.tocsr()
    else:
        T_bulk_ext = T

    s0 = ((0,) * n, (0,) * n)
    s_target = ((2,) * n, (0,) * n)
    if s_target not in idx_ext:
        print(f"  [WARN] estado alvo nao alcancavel em idx_ext.")
        cnt_mat = None
    else:
        # Usa float64 para evitar overflow; converte para int no final.
        v = np.zeros(K_ext, dtype=np.float64)
        v[idx_ext[s0]] = 1.0
        # n-2 bulk steps. v_{t+1}[s'] = sum_s v_t[s] * T[s, s'] = (T.T @ v)
        for _ in range(n - 2):
            v = T_bulk_ext.T @ v
        # 1 no_cp2 step
        v = T_nc2.T @ v
        cnt_mat = int(round(v[idx_ext[s_target]]))
        print(f"  contagem via T = {cnt_mat}")
        if cnt is not None and cnt_mat != cnt:
            print(f"  [ERRO] DP={cnt}, matriz={cnt_mat} -- INCONSISTENCIA")

    if save:
        # Salvar T_bulk (sobre idx original) e idx mapping
        sp.save_npz(DATA / f"transfer_n{n}_bulk.npz", T)
        sp.save_npz(DATA / f"transfer_n{n}_nocp2.npz", T_nc2)
        # idx: salvar como JSON com keys como strings
        idx_serialized = {
            f"{s[0]}|{s[1]}": i for s, i in idx_ext.items()
        }
        with open(DATA / f"transfer_n{n}_idx.json", "w") as f:
            json.dump(idx_serialized, f)
        print(f"  salvo em data/transfer_n{n}_*.npz")

    return {
        "n": n,
        "K": K,
        "K_ext": K_ext,
        "nnz_bulk": int(nnz),
        "nnz_nocp2": int(T_nc2.nnz),
        "count_dp": cnt,
        "count_matrix": cnt_mat,
        "time_build_s": t1 - t0,
        "time_nocp2_s": t3 - t2,
        "time_dp_s": (t5 - t4) if cnt is not None else None,
    }


if __name__ == "__main__":
    import sys

    ns = [int(x) for x in sys.argv[1:]] if len(sys.argv) > 1 else [6]
    summary = []
    for n in ns:
        info = run_for(n, save=True)
        summary.append(info)

    # Append/merge em state_count.json
    sc_path = DATA / "state_count.json"
    existing = []
    if sc_path.exists():
        with open(sc_path) as f:
            existing = json.load(f)
    existing_by_n = {r["n"]: r for r in existing}
    for r in summary:
        existing_by_n[r["n"]] = r
    merged = sorted(existing_by_n.values(), key=lambda r: r["n"])
    with open(sc_path, "w") as f:
        json.dump(merged, f, indent=2)
    print("\nResumo salvo em data/state_count.json")
