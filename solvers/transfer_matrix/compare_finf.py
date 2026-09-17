"""
Compara o autovetor dominante v_1 de T_bulk com a frequencia empirica
f_inf(L) medida em incremental_subtour/data/edge_freq_scaling.json.

Marginal de uma aresta no bulk:
  P(e in E_c) = (sum_{transicoes contendo e} v_left[s] * v_right[s'])
              / (lambda_1 * v_left . v_right)

Cada aresta tem "shape" (col_off, src_row, dst_row). Para um tabuleiro
nxn, src_row varia em 0..n-1 e dst_row eh determinado por col_off + delta_r.

Saidas:
  - data/eigenvector_marginals.json
  - data/plots/eigenvector_vs_finf.png
"""

import json
import sys
import time
from collections import defaultdict
from itertools import combinations
from pathlib import Path

import numpy as np
import scipy.sparse as sp

THIS = Path(__file__).parent
DATA = THIS / "data"
PLOTS = DATA / "plots"
PLOTS.mkdir(exist_ok=True)

# Reusar a funcao forward_options e transitions de transfer_build.py
sys.path.insert(0, str(THIS))
from transfer_build import forward_options


# ---------------------------------------------------------------------------
# transicoes que tambem retornam o E_c (lista de edges) usado
# ---------------------------------------------------------------------------

def transitions_with_edges(state, n: int, forbid_cp2: bool = False):
    """Yield (s_next, E_c) onde E_c eh tupla de ((src_row, 'cp1'|'cp2', dst_row), ...)."""
    deg_c, deg_cp1 = state

    row_choices = []
    for r in range(n):
        need = 2 - deg_c[r]
        if need < 0:
            return
        opts = forward_options(n, r)
        if forbid_cp2:
            opts = [o for o in opts if o[0] == 1]
        if need > len(opts):
            return
        row_choices.append((need, opts))

    def rec(r, cp1_inc, cp2_inc, partial_edges):
        if r == n:
            new_cp1 = tuple(deg_cp1[i] + cp1_inc[i] for i in range(n))
            new_cp2 = tuple(cp2_inc[i] for i in range(n))
            yield (new_cp1, new_cp2), tuple(partial_edges)
            return

        need, opts = row_choices[r]
        if need == 0:
            yield from rec(r + 1, cp1_inc, cp2_inc, partial_edges)
            return

        for combo in combinations(range(len(opts)), need):
            new_edges = []
            valid = True
            cp1_back, cp2_back = [], []
            for idx in combo:
                col_off, dst_r = opts[idx]
                if col_off == 1:
                    if deg_cp1[dst_r] + cp1_inc[dst_r] + 1 > 2:
                        valid = False
                        break
                    cp1_inc[dst_r] += 1
                    cp1_back.append(dst_r)
                    new_edges.append((r, col_off, dst_r))
                else:
                    if cp2_inc[dst_r] + 1 > 2:
                        valid = False
                        break
                    cp2_inc[dst_r] += 1
                    cp2_back.append(dst_r)
                    new_edges.append((r, col_off, dst_r))
            if valid:
                yield from rec(r + 1, cp1_inc, cp2_inc, partial_edges + new_edges)
            for dst_r in cp1_back:
                cp1_inc[dst_r] -= 1
            for dst_r in cp2_back:
                cp2_inc[dst_r] -= 1

    yield from rec(0, [0] * n, [0] * n, [])


# ---------------------------------------------------------------------------
# Carregamento
# ---------------------------------------------------------------------------

def load_idx(n):
    with open(DATA / f"transfer_n{n}_idx.json") as f:
        idx_raw = json.load(f)
    idx = {}
    for k, v in idx_raw.items():
        a, b = k.split("|")
        a_t = tuple(int(x) for x in a.strip("()").split(",") if x.strip())
        b_t = tuple(int(x) for x in b.strip("()").split(",") if x.strip())
        idx[(a_t, b_t)] = v
    return idx


def load_eigenvectors(n):
    z = np.load(DATA / f"eigenvector_dominant_n{n}.npz")
    v_right = z["v_real"] + 1j * z["v_imag"]
    v_left = z["w_real"] + 1j * z["w_imag"]
    lam1 = complex(z["lambda1_real"], z["lambda1_imag"])
    return v_right, v_left, lam1


# ---------------------------------------------------------------------------
# Computa marginais de aresta a partir do autovetor
# ---------------------------------------------------------------------------

def edge_marginals(n, idx, v_right, v_left, lam1):
    """Para cada edge shape (src_row, col_off, dst_row), computa P(e in E_c)
    no bulk usando v_left[s] * v_right[s'] como pesos.

    Retorna dict shape -> float (real).
    """
    # Normalizacao: lam1 * (v_left . v_right)
    norm_denom = lam1 * (v_left @ v_right)
    # numerador acumulado
    edge_w = defaultdict(complex)
    total_w = 0.0 + 0j

    K = len(idx)
    print(f"  iterando {K} estados...")
    t0 = time.time()
    for s, i in idx.items():
        wl = v_left[i]
        if wl == 0:
            continue
        for s_next, E_c in transitions_with_edges(s, n, forbid_cp2=False):
            j = idx.get(s_next)
            if j is None:
                continue
            w = wl * v_right[j]
            total_w += w
            for e in E_c:
                edge_w[e] += w
    t1 = time.time()
    print(f"  {t1 - t0:.2f}s")

    # P(e) = edge_w[e] / (lam1 * (vl.vr))
    # Verificar: total_w deve = lam1 * (vl.vr)
    rel_err = abs(total_w - norm_denom) / abs(norm_denom)
    print(f"  total_w = {total_w.real:+.4f}{total_w.imag:+.4f}j")
    print(f"  norm_denom = lam1 * (vl.vr) = {norm_denom.real:+.4f}{norm_denom.imag:+.4f}j")
    print(f"  rel_err = {rel_err:.3e}")

    marginals = {}
    for e, m in edge_w.items():
        # tomar real part; se imag for grande, alertar
        if abs(m.imag) > 1e-6 * abs(m):
            pass  # ignora silenciosamente; warm warning so soft.
        marginals[e] = (m / norm_denom).real

    return marginals


# ---------------------------------------------------------------------------
# Carregar f_inf empirico
# ---------------------------------------------------------------------------

def load_finf():
    path = THIS.parent / "incremental_subtour" / "data" / "edge_freq_scaling.json"
    with open(path) as f:
        d = json.load(f)
    # Coletar por nivel L, agregando across n
    finf_by_n = {}
    for r in d["resultados"]:
        n = r["n"]
        finf_by_n[n] = {int(L): st["freq_mean"] for L, st in r["por_nivel"].items()}
    return finf_by_n


# ---------------------------------------------------------------------------
# Mapeamento aresta -> nivel L
# ---------------------------------------------------------------------------

def edge_level(n, r1, c1, r2, c2):
    """Nivel L de uma aresta = min(min(r1, n-1-r1, c1, n-1-c1), idem r2,c2).

    Conforme convencao em incremental_subtour: L(v) = dist a' borda (anel).
    L(e) = min(L(v1), L(v2)).
    """
    L1 = min(r1, n - 1 - r1, c1, n - 1 - c1)
    L2 = min(r2, n - 1 - r2, c2, n - 1 - c2)
    return min(L1, L2)


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def run(n):
    print(f"\n=== n = {n} comparison ===")
    idx = load_idx(n)
    v_right, v_left, lam1 = load_eigenvectors(n)
    print(f"  K = {len(idx)}")
    print(f"  lambda_1 = {lam1.real:+.6f}{lam1.imag:+.6f}j")

    shape_marginals = edge_marginals(n, idx, v_right, v_left, lam1)
    print(f"  {len(shape_marginals)} edge shapes encontrados")

    # Converter shape -> nivel L para arestas no INTERIOR (col c em [2, n-3])
    # mas como o "shape" eh (src_row, col_off, dst_row) sem coluna explicita,
    # interpretamos como aresta no bulk: src_col = (n-1)/2 (centro) e
    # avaliamos L com isso.
    #
    # Melhor: enumerar todas as arestas reais do knight graph e atribuir
    # a marginal da shape correspondente.

    # Enumerar arestas reais
    edges_real = []  # list of dict
    for c1 in range(n):
        for r1 in range(n):
            for col_off, dst_r in forward_options(n, r1):
                c2 = c1 + col_off
                r2 = dst_r
                if c2 >= n:
                    continue  # fora do tabuleiro
                shape = (r1, col_off, r2)
                marg_eigen = shape_marginals.get(shape, 0.0)
                L = edge_level(n, r1, c1, r2, c2)
                # is it a "boundary column" edge? c1 < 2 ou c2 > n-3
                is_bulk_col = (c1 >= 2 and c2 <= n - 3)
                edges_real.append({
                    "src": (r1, c1),
                    "dst": (r2, c2),
                    "shape": shape,
                    "level_L": L,
                    "marginal_eigenvector": float(marg_eigen),
                    "bulk_column": bool(is_bulk_col),
                })

    # Agregar por nivel
    by_L = defaultdict(list)
    by_L_bulk = defaultdict(list)
    for e in edges_real:
        by_L[e["level_L"]].append(e["marginal_eigenvector"])
        if e["bulk_column"]:
            by_L_bulk[e["level_L"]].append(e["marginal_eigenvector"])

    finf_by_n = load_finf()
    finf_target = finf_by_n.get(n, {})

    print(f"\n  Marginal por nivel L (todas as arestas):")
    print(f"  {'L':>3} {'n_edges':>8} {'eig_mean':>12} {'eig_std':>10} {'f_inf':>10} {'rel_err':>10}")
    table = []
    for L in sorted(by_L):
        arr = np.array(by_L[L])
        mean = arr.mean()
        std = arr.std()
        finf = finf_target.get(L, None)
        if finf is None or finf == 0:
            rel = float("nan")
        else:
            rel = abs(mean - finf) / finf
        print(f"  {L:>3} {len(arr):>8} {mean:>12.4f} {std:>10.4f} "
              f"{finf if finf else 'NA':>10} {rel:>10.4f}")
        table.append({"L": L, "n_edges": len(arr), "eig_mean": float(mean),
                      "eig_std": float(std), "finf_empirical": finf,
                      "rel_err": float(rel) if not np.isnan(rel) else None})

    print(f"\n  Marginal por nivel L (apenas colunas bulk c1>=2, c2<=n-3):")
    print(f"  {'L':>3} {'n_edges':>8} {'eig_mean':>12} {'eig_std':>10} {'f_inf':>10} {'rel_err':>10}")
    table_bulk = []
    for L in sorted(by_L_bulk):
        arr = np.array(by_L_bulk[L])
        mean = arr.mean()
        std = arr.std()
        finf = finf_target.get(L, None)
        if finf is None or finf == 0:
            rel = float("nan")
        else:
            rel = abs(mean - finf) / finf
        print(f"  {L:>3} {len(arr):>8} {mean:>12.4f} {std:>10.4f} "
              f"{finf if finf else 'NA':>10} {rel:>10.4f}")
        table_bulk.append({"L": L, "n_edges": len(arr), "eig_mean": float(mean),
                           "eig_std": float(std), "finf_empirical": finf,
                           "rel_err": float(rel) if not np.isnan(rel) else None})

    # Salvar
    out = {
        "n": n,
        "lambda_1": [lam1.real, lam1.imag],
        "edges": [
            {
                "src": list(e["src"]),
                "dst": list(e["dst"]),
                "level_L": e["level_L"],
                "marginal_eigenvector": e["marginal_eigenvector"],
                "bulk_column": e["bulk_column"],
            } for e in edges_real
        ],
        "by_level_all": table,
        "by_level_bulk": table_bulk,
    }
    with open(DATA / f"eigenvector_marginals_n{n}.json", "w") as f:
        json.dump(out, f, indent=2)
    print(f"\n  salvo em data/eigenvector_marginals_n{n}.json")

    return out


if __name__ == "__main__":
    ns = [int(x) for x in sys.argv[1:]] if len(sys.argv) > 1 else [6]
    for n in ns:
        run(n)
