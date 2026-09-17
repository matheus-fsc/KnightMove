"""
TAREFAS 0 e 1 — identificar as 3 dimensões proibidas em H_1(G, F_2)
================================================================

T  : 9862 x 80  (incidência tour-aresta, GF(2))
G  : grafo do cavalo 6x6, V=36, E=80, beta_1 = 80 - 36 + 1 = 45
Ham: span GF(2) das linhas de T;  rank(T) = 42

Forbidden = base de um complemento de Ham em H_1, dim = 3.
"""
from __future__ import annotations
import json
import os
import sys
import numpy as np

DATA_IN = os.path.join(os.path.dirname(__file__), "..", "6x6_higher_order", "data")
DATA_OUT = os.path.join(os.path.dirname(__file__), "data", "results")
os.makedirs(DATA_OUT, exist_ok=True)


# ----------------------------------------------------------------------
# GF(2) helpers (matrizes pequenas: usamos bools + xor)
# ----------------------------------------------------------------------

def gf2_rank(M: np.ndarray) -> int:
    """rank de M sobre GF(2) por eliminacao gaussiana in-place."""
    A = M.astype(np.uint8).copy() & 1
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


def gf2_row_reduce(M: np.ndarray):
    """Retorna (rref, pivots, rank). M nao e modificado."""
    A = (M.astype(np.uint8) & 1).copy()
    rows, cols = A.shape
    pivots: list[int] = []
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
        pivots.append(c)
        r += 1
        if r == rows:
            break
    return A[:r], pivots, r


def gf2_in_span(basis_rref: np.ndarray, pivots: list[int], v: np.ndarray) -> bool:
    """v esta em rowspan(basis_rref) se reduzir-o por basis_rref der 0."""
    w = (v.astype(np.uint8) & 1).copy()
    for r, pc in enumerate(pivots):
        if w[pc]:
            w ^= basis_rref[r]
    return not w.any()


# ----------------------------------------------------------------------
# Grafo do cavalo 6x6
# ----------------------------------------------------------------------

def load_graph():
    with open(os.path.join(DATA_IN, "edges_6x6.json")) as f:
        meta = json.load(f)
    edges_uv = [tuple(e) for e in meta["edges_uv"]]
    labels = meta["edge_labels"]
    n_edges = meta["n_edges"]
    V = 36
    assert n_edges == len(edges_uv) == len(labels) == 80
    # mapa aresta -> indice
    edge_index = {e: i for i, e in enumerate(edges_uv)}
    # lista de adjacencias com indices das arestas
    adj: list[list[tuple[int, int]]] = [[] for _ in range(V)]
    for i, (u, v) in enumerate(edges_uv):
        adj[u].append((v, i))
        adj[v].append((u, i))
    return V, edges_uv, labels, edge_index, adj


def spanning_tree_bfs(V: int, adj):
    """BFS de 0. Retorna (parent, parent_edge, tree_edges_idx, nontree_edges_idx)."""
    parent = [-1] * V
    parent_edge = [-1] * V
    visited = [False] * V
    visited[0] = True
    queue = [0]
    tree_edges = set()
    while queue:
        nq = []
        for u in queue:
            for (v, ei) in adj[u]:
                if not visited[v]:
                    visited[v] = True
                    parent[v] = u
                    parent_edge[v] = ei
                    tree_edges.add(ei)
                    nq.append(v)
        queue = nq
    assert all(visited), "grafo desconexo? impossivel para cavalo 6x6"
    nontree = [i for i in range(80) if i not in tree_edges]
    return parent, parent_edge, sorted(tree_edges), nontree


def path_to_root(u: int, parent, parent_edge):
    """retorna lista de arestas no caminho u -> 0 na spanning tree."""
    out = []
    while parent[u] != -1:
        out.append(parent_edge[u])
        u = parent[u]
    return out


def fundamental_cycle(e_idx: int, edges_uv, parent, parent_edge, n_edges=80) -> np.ndarray:
    """Ciclo fundamental: aresta e_idx + caminho na arvore entre seus extremos."""
    u, v = edges_uv[e_idx]
    pu = path_to_root(u, parent, parent_edge)
    pv = path_to_root(v, parent, parent_edge)
    # XOR dos caminhos = caminho de u a v na arvore
    su = set(pu); sv = set(pv)
    sym = (su ^ sv)  # XOR simetrica
    cyc = np.zeros(n_edges, dtype=np.uint8)
    cyc[e_idx] = 1
    for ei in sym:
        cyc[ei] = 1
    return cyc


# ----------------------------------------------------------------------
# Pipeline principal
# ----------------------------------------------------------------------

def main():
    # TAREFA 0 ----------------------------------------------------------
    print("=" * 64)
    print("TAREFA 0  Carregamento + sanity check")
    print("=" * 64)

    T = np.load(os.path.join(DATA_IN, "incidence_matrix_6x6.npy"))
    print(f"T shape       : {T.shape}  dtype={T.dtype}")

    V, edges_uv, labels, edge_index, adj = load_graph()
    E = len(edges_uv)
    beta1 = E - V + 1
    print(f"V             : {V}")
    print(f"E             : {E}")
    print(f"beta_1 = E-V+1: {beta1}")

    rank_T = gf2_rank(T)
    print(f"rank_GF2(T)   : {rank_T}")
    if rank_T != 42:
        print(f"!! rank inesperado: esperado 42, obtido {rank_T}")
        sys.exit(1)
    if beta1 != 45:
        print(f"!! beta_1 inesperado: esperado 45, obtido {beta1}")
        sys.exit(1)
    print("OK: rank(T)=42  beta_1=45  deficit = 3")

    # ---- verifica ainda: cada linha de T tem grau par em cada vertice
    inc = np.zeros((V, E), dtype=np.uint8)
    for i, (u, v) in enumerate(edges_uv):
        inc[u, i] = 1
        inc[v, i] = 1
    deg = (T @ inc.T) % 2  # 9862 x 36
    if deg.any():
        bad = int(deg.any(axis=1).sum())
        print(f"!! {bad} linhas de T nao sao ciclos (boundary != 0)")
        sys.exit(1)
    print("OK: todas as linhas de T sao ciclos (boundary = 0)")

    # TAREFA 1.1 --------------------------------------------------------
    print()
    print("=" * 64)
    print("TAREFA 1.1  Base de H_1(G, F_2) via spanning tree")
    print("=" * 64)
    parent, parent_edge, tree_edges, nontree = spanning_tree_bfs(V, adj)
    print(f"spanning tree : {len(tree_edges)} arestas  (esperado V-1 = 35)")
    print(f"co-arvore     : {len(nontree)} arestas    (esperado beta_1 = 45)")
    assert len(tree_edges) == V - 1
    assert len(nontree) == beta1

    Base_H1 = np.zeros((beta1, E), dtype=np.uint8)
    for k, e in enumerate(nontree):
        Base_H1[k] = fundamental_cycle(e, edges_uv, parent, parent_edge, n_edges=E)
    # cada ciclo fundamental e mesmo um ciclo (boundary=0)?
    bH1 = (Base_H1 @ inc.T) % 2
    assert not bH1.any(), "ciclo fundamental com boundary != 0"
    rank_H1 = gf2_rank(Base_H1)
    print(f"rank(Base_H1) : {rank_H1}")
    if rank_H1 != beta1:
        print("!! Base_H1 nao tem rank beta_1 -- bug na spanning tree")
        sys.exit(1)
    print("OK")

    # TAREFA 1.2 --------------------------------------------------------
    print()
    print("=" * 64)
    print("TAREFA 1.2  Base de Ham (row reduction de T)")
    print("=" * 64)
    Base_Ham_rref, pivots_ham, rank_ham = gf2_row_reduce(T)
    print(f"Base_Ham shape: {Base_Ham_rref.shape}   rank={rank_ham}")
    assert rank_ham == 42
    # tambem deve ser subespaco de H1 (cada linha tem boundary=0)
    bH = (Base_Ham_rref @ inc.T) % 2
    assert not bH.any()

    # TAREFA 1.3 --------------------------------------------------------
    print()
    print("=" * 64)
    print("TAREFA 1.3  Complemento Forbidden = base de H_1 / Ham")
    print("=" * 64)

    # Vamos manter dois "estados": rref incremental do span Ham+coisas-adicionadas,
    # comecando de Base_Ham_rref, e a lista das adicoes em forma 'original'.
    current = Base_Ham_rref.copy()
    current_pivots = list(pivots_ham)
    forbidden_rows = []          # vetores f_i originais (antes de reducao)
    forbidden_source_edge = []   # qual aresta nao-arvore gerou cada f_i

    for k, e in enumerate(nontree):
        v = Base_H1[k].copy()
        if gf2_in_span(current, current_pivots, v):
            continue
        # ainda nao esta no span -> adicionar
        forbidden_rows.append(v.copy())
        forbidden_source_edge.append(e)
        # incorporar em current rref
        w = v.copy()
        for r, pc in enumerate(current_pivots):
            if w[pc]:
                w ^= current[r]
        # w tem pelo menos um 1 (nao estava no span). Achar primeiro 1:
        new_pivot = int(np.argmax(w))
        # inserir em current na posicao certa (para manter ordenado por pivot)
        insert_at = 0
        for i, pc in enumerate(current_pivots):
            if pc > new_pivot:
                insert_at = i
                break
            insert_at = i + 1
        # eliminar new_pivot nas linhas existentes acima
        for r in range(current.shape[0]):
            if r != insert_at and current[r, new_pivot]:
                current[r] ^= w
        current = np.insert(current, insert_at, w, axis=0)
        current_pivots.insert(insert_at, new_pivot)
        if len(forbidden_rows) == 3:
            # confere se ja completamos a base de H_1
            if current.shape[0] == 45:
                # ainda assim, continuar varrendo seria redundante
                break

    n_added = len(forbidden_rows)
    print(f"forbidden adicionados: {n_added} (esperado 3)")
    if n_added != 3:
        print(f"!! adicionados {n_added} vetores, esperado 3")
        sys.exit(1)

    Forbidden = np.array(forbidden_rows, dtype=np.uint8)
    # Verificacoes finais
    combined = np.vstack([Base_Ham_rref, Forbidden])
    rank_comb = gf2_rank(combined)
    print(f"rank([Base_Ham | Forbidden])  : {rank_comb}  (esperado {beta1})")
    assert rank_comb == beta1

    # Nenhum f_i no span de Base_Ham
    for i, f in enumerate(Forbidden):
        assert not gf2_in_span(Base_Ham_rref, pivots_ham, f), f"f{i} cai no span de Ham"
    print("OK: nenhum f_i no span(Ham)")

    # Span(Ham) interseccao Span(Forbidden) = {0}
    # Equivalente a: rank(Ham) + rank(Forbidden) == rank(Ham ∪ Forbidden) == 45
    rank_F = gf2_rank(Forbidden)
    print(f"rank(Forbidden)               : {rank_F}")
    assert rank_F == 3
    assert rank_ham + rank_F == rank_comb, "intersecao nao trivial!"
    print("OK: Ham ∩ Forbidden = {0}")

    # ---- salvar resultados
    np.save(os.path.join(DATA_OUT, "forbidden_cycles.npy"), Forbidden)
    np.save(os.path.join(DATA_OUT, "base_h1.npy"), Base_H1)
    np.save(os.path.join(DATA_OUT, "base_ham_rref.npy"), Base_Ham_rref)
    np.save(os.path.join(DATA_OUT, "incidence_vertex_edge.npy"), inc)

    json_out = {
        "V": V, "E": E, "beta_1": beta1,
        "rank_T": int(rank_T),
        "deficit": int(beta1 - rank_T),
        "forbidden": [],
    }
    for i, (f, src) in enumerate(zip(Forbidden, forbidden_source_edge)):
        active = [int(j) for j in np.where(f == 1)[0]]
        json_out["forbidden"].append({
            "index": i,
            "source_nontree_edge": int(src),
            "source_edge_label": labels[src],
            "support_size": len(active),
            "edge_indices": active,
            "edge_labels": [labels[j] for j in active],
        })
    with open(os.path.join(DATA_OUT, "forbidden_cycles.json"), "w") as f:
        json.dump(json_out, f, indent=2)

    print()
    print("=" * 64)
    print("RESUMO")
    print("=" * 64)
    for entry in json_out["forbidden"]:
        print(f"f_{entry['index']}: |S|={entry['support_size']:2d}  "
              f"src={entry['source_edge_label']}")
    print()
    print(f"escritos:")
    print(f"  {DATA_OUT}/forbidden_cycles.npy")
    print(f"  {DATA_OUT}/forbidden_cycles.json")
    print(f"  {DATA_OUT}/base_h1.npy")
    print(f"  {DATA_OUT}/base_ham_rref.npy")
    print(f"  {DATA_OUT}/incidence_vertex_edge.npy")


if __name__ == "__main__":
    main()
