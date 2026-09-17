"""
TAREFA 3 - testar 4 hipoteses para a obstrucao Ham != H_1.

H1  paridade bipartida
H2  obstrucao por corte (existe corte S t.q. todo tour tem |cut(S) ∩ tour| par
    mas f_i da impar)
H3  obstrucao por sub-tour (S_i nao cobre os 36 vertices)
H4  obstrucao por projecao: o subespaco gerado por TODOS os 2-fatores
    (uniao disjunta de ciclos cobrindo V, sem necessidade de conectividade)
    tem rank 45?
"""
from __future__ import annotations
import json
import os
from itertools import combinations, product
import numpy as np

THIS = os.path.dirname(__file__)
RES = os.path.join(THIS, "data", "results")
HO = os.path.join(THIS, "..", "6x6_higher_order", "data")


def gf2_rank(M: np.ndarray) -> int:
    A = (M.astype(np.uint8) & 1).copy()
    rows, cols = A.shape
    r = 0
    for c in range(cols):
        piv = None
        for i in range(r, rows):
            if A[i, c]:
                piv = i
                break
        if piv is None:
            continue
        if piv != r:
            A[[r, piv]] = A[[piv, r]]
        mask = A[:, c].astype(bool)
        mask[r] = False
        if mask.any():
            A[mask] ^= A[r]
        r += 1
        if r == rows:
            break
    return r


def load_all():
    F = np.load(os.path.join(RES, "forbidden_cycles.npy"))
    T = np.load(os.path.join(HO, "incidence_matrix_6x6.npy"))
    inc = np.load(os.path.join(RES, "incidence_vertex_edge.npy"))
    with open(os.path.join(HO, "edges_6x6.json")) as fh:
        meta = json.load(fh)
    return F, T, inc, [tuple(e) for e in meta["edges_uv"]], meta["edge_labels"]


# ----------------------------------------------------------------------
# H1  Paridade bipartida
# ----------------------------------------------------------------------
def hypothesis_paridade(F: np.ndarray, edges_uv) -> dict:
    """
    O grafo do cavalo eh bipartido: cor(v) = (row + col) mod 2.
    Como tem-se a bipartite, TODA aresta (u,v) liga cores opostas - nao ha
    arestas "mesmo lado". Logo o teste da forma literal e degenerado.

    Reformulamos: para QUALQUER 2-coloracao dos vertices c: V -> F_2,
    associamos o numero de arestas "mesma cor" -- isto e linear em GF(2).
    Procurar c tal que f_i da somatorio impar mas TODA linha de T da par
    fornece um corte invalidante.

    Vamos exaustivamente testar c via vetor x in F_2^36 (excluindo 0 e all-1).
    Para cada x, definir edge_val(e) = x_u XOR x_v (= 0 sse mesma cor).
    "edges_same_color"(c) = #{e: x_u==x_v}.
    """
    # Para cada x in F_2^36, gerar vetor coloracao e testar paridade
    # Mas 2^36 e demais. Vamos restringir a bipartites estruturais relevantes:
    # cor xadrez, paridade de linha, paridade de coluna, etc.
    V = 36
    colorings: list[tuple[str, np.ndarray]] = []
    # cor xadrez
    chess = np.array([(v // 6 + v % 6) % 2 for v in range(V)], dtype=np.uint8)
    colorings.append(("chess", chess))
    colorings.append(("row_parity", np.array([(v // 6) % 2 for v in range(V)], dtype=np.uint8)))
    colorings.append(("col_parity", np.array([(v % 6) % 2 for v in range(V)], dtype=np.uint8)))
    # quadrantes
    colorings.append(("top_half", np.array([1 if (v // 6) < 3 else 0 for v in range(V)], dtype=np.uint8)))
    colorings.append(("left_half", np.array([1 if (v % 6) < 3 else 0 for v in range(V)], dtype=np.uint8)))

    results = []
    for name, col in colorings:
        # vetor por aresta: same_color(e) = 1 - (col[u] XOR col[v])
        same = np.array([1 - int(col[u] ^ col[v]) for (u, v) in edges_uv], dtype=np.uint8)
        # paridade do F_i * same eh um ESCALAR em F_2
        Fparity = (F @ same) % 2
        results.append({"coloring": name, "F_parity": [int(x) for x in Fparity]})

    return {
        "note": (
            "Cavalo eh bipartido pela cor xadrez -> arestas 'mesma cor' = 0 "
            "(coloracao chess). Teste so e nao-trivial para 2-coloracoes "
            "diferentes da bipartite natural."
        ),
        "tests": results,
    }


# ----------------------------------------------------------------------
# H2  Obstrucao por corte
# ----------------------------------------------------------------------
def cut_vector(S: set, edges_uv, n_edges: int) -> np.ndarray:
    """Indicador das arestas com exatamente um endpoint em S."""
    v = np.zeros(n_edges, dtype=np.uint8)
    for ei, (u, w) in enumerate(edges_uv):
        a = int(u in S); b = int(w in S)
        if a != b:
            v[ei] = 1
    return v


def hypothesis_corte(F: np.ndarray, T: np.ndarray, edges_uv) -> dict:
    """
    Para que f_i seja 'obstruido por corte', precisamos de S subset V t.q.
       (T @ cut(S)) % 2  = 0 (vetor)        --> todo tour cruza S par num
       (F[i] @ cut(S)) % 2 = 1               --> f_i cruza S impar

    Como dim(Ham) = 42 em F_2^E, o anulador de Ham em (F_2^E)^* (= row-space
    de T^perp) tem dim E - 42 = 38. Mas nem todo vetor de F_2^E e um cut(S).
    Vamos enumerar S subset V via 2^36 (impossivel). Em vez disso,
    enumeramos subconjuntos pequenos S (|S| <= k) e tambem cuts induzidos
    por classes simples (linhas, colunas, quadrantes, etc.).
    """
    V = 36
    E = len(edges_uv)
    bad_cuts_per_fi = [[] for _ in range(F.shape[0])]

    # exaustivo |S|=1..3 (manageable)
    pools: list[tuple[str, list[set]]] = []
    pools.append(("|S|=1", [{v} for v in range(V)]))
    pools.append(("|S|=2", [set(c) for c in combinations(range(V), 2)]))
    pools.append(("|S|=3", [set(c) for c in combinations(range(V), 3)]))

    # cortes "estruturais": linhas, colunas, quadrantes, casas escuras
    structural: list[tuple[str, set]] = []
    for r in range(6):
        structural.append((f"row_eq_{6-r}", {v for v in range(V) if v // 6 == r}))
    for c in range(6):
        structural.append((f"col_eq_{chr(65+c)}", {v for v in range(V) if v % 6 == c}))
    structural.append(("top_half", {v for v in range(V) if v // 6 < 3}))
    structural.append(("bottom_half", {v for v in range(V) if v // 6 >= 3}))
    structural.append(("left_half", {v for v in range(V) if v % 6 < 3}))
    structural.append(("right_half", {v for v in range(V) if v % 6 >= 3}))
    structural.append(("chess_white", {v for v in range(V) if (v // 6 + v % 6) % 2 == 0}))
    structural.append(("dark_4_center", {14, 15, 20, 21}))  # C4 D4 C3 D3
    structural.append(("corners", {0, 5, 30, 35}))

    summary = {"cuts_blocked_by_Ham": 0, "cuts_blocking_some_fi": 0, "examples": [], "tests_run": 0}

    def test_cut(name: str, S: set):
        nonlocal summary
        cv = cut_vector(S, edges_uv, E)
        if not cv.any():
            return
        # todo tour par?
        if ((T @ cv) % 2).any():
            return
        summary["cuts_blocked_by_Ham"] += 1
        # f_i quebra?
        Fp = (F @ cv) % 2
        if Fp.any():
            summary["cuts_blocking_some_fi"] += 1
            for i in range(F.shape[0]):
                if Fp[i]:
                    bad_cuts_per_fi[i].append(name)
            if len(summary["examples"]) < 12:
                summary["examples"].append({
                    "S_name": name, "size_S": len(S),
                    "S": sorted(S),
                    "F_parity": [int(x) for x in Fp],
                })

    for name, S in structural:
        summary["tests_run"] += 1
        test_cut(name, S)

    for label, lst in pools:
        for S in lst:
            summary["tests_run"] += 1
            test_cut(f"{label}:{sorted(S)}", S)

    summary["bad_cuts_per_fi"] = [len(x) for x in bad_cuts_per_fi]
    summary["sample_bad_cuts_f0"] = bad_cuts_per_fi[0][:6]
    summary["sample_bad_cuts_f1"] = bad_cuts_per_fi[1][:6]
    summary["sample_bad_cuts_f2"] = bad_cuts_per_fi[2][:6]
    return summary


# ----------------------------------------------------------------------
# H3  Obstrucao por sub-tour
# ----------------------------------------------------------------------
def hypothesis_subtour(F: np.ndarray, edges_uv) -> dict:
    """ S_i cobre todos os 36 vertices? se NAO, eh impossivel como
    subgrafo de tour hamiltoniano. """
    V = 36
    out = []
    for i, f in enumerate(F):
        verts = set()
        for ei, (u, v) in enumerate(edges_uv):
            if f[ei]:
                verts.add(u); verts.add(v)
        out.append({
            "f_i": i, "n_vertices_covered": len(verts),
            "covers_all_36": len(verts) == V,
            "missing_vertices": sorted(set(range(V)) - verts) if len(verts) < V else [],
        })
    return {"per_fi": out}


# ----------------------------------------------------------------------
# H4  Obstrucao por projecao (2-fatores)
# ----------------------------------------------------------------------
def hypothesis_projecao(F: np.ndarray, T: np.ndarray, inc: np.ndarray,
                       edges_uv, sample_cap: int = 200_000) -> dict:
    """
    2-fator do grafo do cavalo 6x6 = subconjunto de arestas onde TODO vertice
    tem grau exatamente 2. Inclui tours hamiltonianos + unioes disjuntas de
    ciclos. Enumera-los exaustivamente eh tipicamente vavavivel para 6x6
    (centenas a milhares); mas como tambem queremos so estimativa do rank,
    podemos:
      a) usar enumeracao backtracking com poda de grau
      b) amostrar via SAT / random
      c) usar o resultado conhecido: para 6x6, # 2-fatores eh moderado
    Implementamos (a) limitado por sample_cap.
    """
    V = 36
    E = len(edges_uv)
    # adjacencia: para cada vertice, lista de (vizinho, edge_idx)
    adj = [[] for _ in range(V)]
    for ei, (u, v) in enumerate(edges_uv):
        adj[u].append((v, ei))
        adj[v].append((u, ei))

    # backtracking: para vertices em ordem, escolher 2 das suas arestas restantes.
    # Usamos abordagem por "match perfeito em 2-regular": dificil de baixo nivel.
    # Reduzimos a problema: para cada vertice v, exatamente 2 vizinhos sao escolhidos.
    # Implementamos via backtracking decidindo arestas em ordem; corte por grau.

    deg = [0] * V
    used = np.zeros(E, dtype=np.uint8)
    factors = []   # lista de vetores indicadores
    found = [0]

    # ordenar vertices pelo grau crescente
    sorted_v = sorted(range(V), key=lambda v: len(adj[v]))
    v_order = {v: i for i, v in enumerate(sorted_v)}
    # arestas: para cada par (u,v), ei. Vamos decidir em ordem fixa de arestas.

    def can_extend():
        # corte: para todo vertice, deg restante (incluindo arestas livres) >= 2
        # pega arestas restantes (nao decididas) - mas backtracking por aresta:
        return True

    def feasible(deg, edge_decided, edges_left_per_v):
        for v in range(V):
            if deg[v] > 2:
                return False
            if deg[v] + edges_left_per_v[v] < 2:
                return False
        return True

    edges_left_per_v = [len(adj[v]) for v in range(V)]

    def backtrack(idx):
        if found[0] >= sample_cap:
            return
        if idx == E:
            if all(d == 2 for d in deg):
                vec = used.copy()
                factors.append(vec)
                found[0] += 1
            return
        u, v = edges_uv[idx]
        # opcao 1: nao usar
        edges_left_per_v[u] -= 1
        edges_left_per_v[v] -= 1
        if feasible(deg, idx + 1, edges_left_per_v):
            backtrack(idx + 1)
        edges_left_per_v[u] += 1
        edges_left_per_v[v] += 1
        if found[0] >= sample_cap:
            return
        # opcao 2: usar
        if deg[u] < 2 and deg[v] < 2:
            used[idx] = 1
            deg[u] += 1; deg[v] += 1
            edges_left_per_v[u] -= 1
            edges_left_per_v[v] -= 1
            if feasible(deg, idx + 1, edges_left_per_v):
                backtrack(idx + 1)
            edges_left_per_v[u] += 1
            edges_left_per_v[v] += 1
            deg[u] -= 1; deg[v] -= 1
            used[idx] = 0

    backtrack(0)

    n_factors = len(factors)
    if n_factors == 0:
        return {"error": "no 2-factor found", "n_factors": 0}

    M = np.array(factors, dtype=np.uint8)
    # confirma: rows sao ciclos? (boundary = 0)
    bd = (M @ inc.T) % 2
    assert not bd.any()
    rank2f = gf2_rank(M)
    # rank(M ∪ T) = ?
    combined = np.vstack([T, M])
    rank_comb = gf2_rank(combined)
    # M contem Forbidden?
    F_inside = []
    # para cada f_i, testar se esta em rowspan(M)
    # mas matriz F + M tem rank == rank(M) sse F subset rowspan(M)
    F_count = 0
    F_load = np.load(os.path.join(RES, "forbidden_cycles.npy"))
    for i in range(F_load.shape[0]):
        rk = gf2_rank(np.vstack([M, F_load[i:i+1]]))
        if rk == rank2f:
            F_count += 1
    return {
        "n_2_factors_enumerated": n_factors,
        "rank_2_factors": int(rank2f),
        "rank_T_union_2f": int(rank_comb),
        "rank_T": int(gf2_rank(T)),
        "forbidden_in_2f_span": int(F_count),
        "interpretation": "se rank_2_factors == 45 e forbidden_in_2f_span == 3 -> H4 CONFIRMADA",
    }


# ----------------------------------------------------------------------
# main
# ----------------------------------------------------------------------
def main():
    F, T, inc, edges_uv, labels = load_all()
    print(f"F shape {F.shape}, T shape {T.shape}")

    print()
    print("=" * 64)
    print("H1  paridade bipartida")
    print("=" * 64)
    h1 = hypothesis_paridade(F, edges_uv)
    print(h1["note"])
    for t in h1["tests"]:
        print(f"  coloring={t['coloring']:>12s}  F_parity={t['F_parity']}")
    # se TODO tour da paridade par (= invariante), mas algum f_i da impar -> obstrucao
    # como T eh GF(2)-linear, T@same eh um vetor; "todo tour da par" sse esse vetor for 0.
    # vamos calcular e reportar:
    print()
    h1_verdict_lines = []
    h1_genuine = []
    for t in h1["tests"]:
        name = t['coloring']
        col = None
        for n2, c2 in [("chess", "chess"),
                       ("row_parity", "row_parity"),
                       ("col_parity", "col_parity"),
                       ("top_half", "top_half"),
                       ("left_half", "left_half")]:
            if n2 == name:
                col = c2
        # recomputar coluna
        V = 36
        if name == "chess":
            cv = np.array([(v // 6 + v % 6) % 2 for v in range(V)], dtype=np.uint8)
        elif name == "row_parity":
            cv = np.array([(v // 6) % 2 for v in range(V)], dtype=np.uint8)
        elif name == "col_parity":
            cv = np.array([(v % 6) % 2 for v in range(V)], dtype=np.uint8)
        elif name == "top_half":
            cv = np.array([1 if v // 6 < 3 else 0 for v in range(V)], dtype=np.uint8)
        elif name == "left_half":
            cv = np.array([1 if v % 6 < 3 else 0 for v in range(V)], dtype=np.uint8)
        same = np.array([1 - int(cv[u] ^ cv[v]) for (u, v) in edges_uv], dtype=np.uint8)
        T_par = (T @ same) % 2
        all_par_for_tours = not T_par.any()
        F_par = t["F_parity"]
        h1_verdict_lines.append(
            f"  {name:>12s}: todo tour par? {all_par_for_tours}  "
            f"f_parity={F_par}"
        )
        if all_par_for_tours and any(F_par):
            h1_genuine.append(name)
    print("\n".join(h1_verdict_lines))
    h1_status = ("CONFIRMADA (parcial)" if h1_genuine else "REFUTADA")
    print(f"--> H1 status: {h1_status}  obstrucoes: {h1_genuine}")

    print()
    print("=" * 64)
    print("H2  obstrucao por corte")
    print("=" * 64)
    h2 = hypothesis_corte(F, T, edges_uv)
    print(f"cortes testados : {h2['tests_run']}")
    print(f"cortes 'pares p/ tours': {h2['cuts_blocked_by_Ham']}")
    print(f"cortes que ainda separam algum f_i: {h2['cuts_blocking_some_fi']}")
    print(f"distribuicao por f_i: {h2['bad_cuts_per_fi']}")
    if h2["examples"]:
        print("exemplos:")
        for ex in h2["examples"]:
            print(f"  {ex['S_name']:>40s}  S={ex['S']}  Fpar={ex['F_parity']}")
    if h2["cuts_blocking_some_fi"] > 0:
        h2_status = "CONFIRMADA"
    else:
        h2_status = "INCONCLUSIVA (testes limitados)"
    print(f"--> H2 status: {h2_status}")

    print()
    print("=" * 64)
    print("H3  obstrucao por sub-tour (cobertura de vertices)")
    print("=" * 64)
    h3 = hypothesis_subtour(F, edges_uv)
    for r in h3["per_fi"]:
        print(f"  f_{r['f_i']}  cobre {r['n_vertices_covered']}/36 vertices  "
              f"all? {r['covers_all_36']}")
        if r["missing_vertices"]:
            mv = r["missing_vertices"]
            # converter para labels
            mv_lab = [f"{chr(65 + v%6)}{6 - v//6}" for v in mv]
            print(f"    faltantes ({len(mv)}): {mv_lab}")
    h3_status = "CONFIRMADA (todos f_i deixam vertices nao-cobertos)" \
        if all(not r["covers_all_36"] for r in h3["per_fi"]) else "REFUTADA"
    print(f"--> H3 status: {h3_status}")

    print()
    print("=" * 64)
    print("H4  obstrucao por projecao (2-fatores)")
    print("=" * 64)
    print("enumerando 2-fatores (uniao disjunta de ciclos cobrindo V)...")
    h4 = hypothesis_projecao(F, T, inc, edges_uv, sample_cap=200_000)
    for k, v in h4.items():
        print(f"  {k}: {v}")
    if "error" in h4:
        h4_status = "ERRO"
    else:
        h4_status = ("CONFIRMADA" if h4["rank_2_factors"] == 45
                     and h4["forbidden_in_2f_span"] == 3
                     else "REFUTADA")
    print(f"--> H4 status: {h4_status}")

    print()
    print("=" * 64)
    print("BONUS  funcionais duais (3 invariantes lineares que detectam f_i)")
    print("=" * 64)
    # Procuramos 3 vetores w_1, w_2, w_3 in F_2^80 tais que:
    #   <t, w_i> = 0 para todo tour t
    #   <f_j, w_i> = delta_{ij}  (ou equivalente: posto cheio em F)
    # i.e. base do anulador(Ham) modulo anulador(H_1).
    # Algoritmo:
    #   Ann(Ham) = ker(T) em F_2^80 (vetores w t.q. T w = 0).
    #   dentro de Ann(Ham), achar 3 vetores que SAO LINEARMENTE INDEPENDENTES
    #   contra a base Forbidden (i.e. matriz F @ W e invertivel 3x3).
    # Implementacao: encontramos ker(T) por row reduction da transposta.
    E = T.shape[1]
    # ker(T): solucionar T w = 0. equivalente a achar base do nullspace.
    # row-reduce T -> identificar colunas livres
    A = (T.copy().astype(np.uint8)) & 1
    rows, cols = A.shape
    r = 0
    pivot_col_of_row = []
    pivot_cols = []
    free_cols = []
    col_to_pivot_row = [-1] * cols
    Aw = A.copy()
    for c in range(cols):
        piv = None
        for i in range(r, rows):
            if Aw[i, c]:
                piv = i; break
        if piv is None:
            free_cols.append(c)
            continue
        if piv != r:
            Aw[[r, piv]] = Aw[[piv, r]]
        mask = Aw[:, c].astype(bool); mask[r] = False
        if mask.any():
            Aw[mask] ^= Aw[r]
        pivot_cols.append(c)
        col_to_pivot_row[c] = r
        r += 1
        if r == rows:
            # restantes sao livres
            for cc in range(c + 1, cols):
                free_cols.append(cc)
            break
    print(f"rank T = {r}  pivots={len(pivot_cols)}  free={len(free_cols)}  E={E}")
    # kernel basis: 1 vetor por coluna livre
    kernel = []
    for fc in free_cols:
        w = np.zeros(E, dtype=np.uint8)
        w[fc] = 1
        for pc in pivot_cols:
            row_idx = col_to_pivot_row[pc]
            if Aw[row_idx, fc]:
                w[pc] = 1
        kernel.append(w)
    Kmat = np.array(kernel, dtype=np.uint8)  # (E - rank) x E
    # sanity: T @ K^T = 0 mod 2
    check = (T @ Kmat.T) % 2
    assert not check.any(), "kernel construido errado"
    print(f"kernel(T) basis: shape {Kmat.shape}")

    # projetar nos 3 funcionais detectores de f_i: queremos W (3xE) com F @ W^T = I_3
    # F @ K^T tem shape (3, dim_ker) - precisamos selecionar/transformar para obter um 3x3 invertivel
    FK = (F @ Kmat.T) % 2  # 3 x dim_ker
    # gaussian elim sobre colunas para extrair 3x3 invertivel
    FKc = FK.copy()
    chosen_cols = []
    used = [False] * FKc.shape[1]
    for row in range(3):
        # achar coluna livre com 1 nessa linha
        sel = None
        for c in range(FKc.shape[1]):
            if used[c]: continue
            if FKc[row, c] == 1:
                sel = c; break
        if sel is None:
            print(f"!! falha: nao ha funcional em kernel(T) que detecte f_{row}")
            chosen_cols = None; break
        used[sel] = True
        chosen_cols.append(sel)
        # eliminar abaixo na mesma linha por XOR de colunas
        for c in range(FKc.shape[1]):
            if c == sel or used[c]:
                continue
            if FKc[row, c] == 1:
                FKc[:, c] ^= FKc[:, sel]
    detectors = None
    if chosen_cols is not None:
        # W tem 3 linhas, cada uma = combinacao das colunas escolhidas de K
        Wbase = Kmat[chosen_cols]   # 3 x E
        # Aplica transformacao para que F @ W^T = I_3. Note que apos elim, FK[:, chosen_cols]
        # eh upper triangular com 1's na diagonal; faltam zerar acima da diagonal.
        # como ja estamos em diagonal-1's-only (back-substitution), agora torna identidade:
        # vou recombinar manualmente:
        # construir 3x3 = F @ Wbase^T
        M3 = (F @ Wbase.T) % 2
        # inverter M3 em F_2 (M3 e quadrada 3x3)
        def inv_gf2(M):
            n = M.shape[0]
            A = np.hstack([M.copy() & 1, np.eye(n, dtype=np.uint8)]).astype(np.uint8)
            for c in range(n):
                piv = None
                for i in range(c, n):
                    if A[i, c]:
                        piv = i; break
                if piv is None:
                    return None
                if piv != c:
                    A[[c, piv]] = A[[piv, c]]
                for i in range(n):
                    if i != c and A[i, c]:
                        A[i] ^= A[c]
            return A[:, n:]
        inv = inv_gf2(M3)
        if inv is not None:
            # queremos D 3xE com F @ D^T = I.  D = X @ Wbase => F W^T X^T = I
            # => X^T = M3^{-1} => X = (M3^{-1})^T  => D = (M3^{-1})^T @ Wbase
            detectors = (inv.T @ Wbase) % 2  # 3 x E
            verify = (F @ detectors.T) % 2
            print(f"  F @ detectors^T = (deve ser I_3):\n{verify}")
            tv = (T @ detectors.T) % 2
            print(f"  todos tours -> detectors aplicado:  zeros? {not tv.any()}")
            for i in range(3):
                supp = [int(j) for j in np.where(detectors[i] == 1)[0]]
                supp_lab = [labels[j] for j in supp]
                print(f"  detector phi_{i}: |supp|={len(supp)}  arestas={supp_lab}")
    if detectors is not None:
        np.save(os.path.join(RES, "dual_detectors.npy"), detectors)
        det_json = []
        for i in range(3):
            supp = [int(j) for j in np.where(detectors[i] == 1)[0]]
            det_json.append({
                "phi_i": i, "support_size": len(supp),
                "edge_indices": supp,
                "edge_labels": [labels[j] for j in supp],
            })
        with open(os.path.join(RES, "dual_detectors.json"), "w") as fh:
            json.dump(det_json, fh, indent=2)
        print(f"  -> dual_detectors.{{npy,json}}")

    # ---- salvar
    out = {
        "H1": {"status": h1_status, "details": h1, "genuine_colorings": h1_genuine},
        "H2": {"status": h2_status, "details": h2},
        "H3": {"status": h3_status, "details": h3,
               "note": "tecnicamente trivial: qualquer vetor F_2 nao-tour deixa "
                       "vertices descobertos. nao distingue proibidos de outros."},
        "H4": {"status": h4_status, "details": h4},
        "dual_detectors_summary": ("3 funcionais lineares phi_i in F_2^E "
                                   "salvos em dual_detectors.npy "
                                   "(phi_i(t)=0 forall tour, phi_i(f_j)=delta_ij)"),
    }
    out_path = os.path.join(RES, "obstruction_hypotheses.json")
    with open(out_path, "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"\n--> {out_path}")


if __name__ == "__main__":
    main()
