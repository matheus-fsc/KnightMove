"""
TAREFA 2 - Geometria dos ciclos proibidos.

Para cada f_i:
  2.1  plot 6x6 destacando S_i = {e : f_i[e]=1}
  2.2  topologia: |S_i|, graus, componentes (= ciclos simples que compoem f_i)
  2.3  posicao: linhas/colunas tocadas, simetria D_4, winding em torno do centro
  2.4  relacao com obrigatorias (8 arestas) e exclusoes minimais (88 / 1776 / 17004)
"""
from __future__ import annotations
import json
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

THIS = os.path.dirname(__file__)
RES = os.path.join(THIS, "data", "results")
PLOTS = os.path.join(THIS, "data", "plots")
HO = os.path.join(THIS, "..", "6x6_higher_order", "data")
os.makedirs(PLOTS, exist_ok=True)


# ----------------------------------------------------------------------
# Coordenadas do tabuleiro 6x6
#   vertice v = row*6 + col, com row=0 = linha 6 (topo).
#   A6=0, F1=35
# ----------------------------------------------------------------------

def v_to_xy(v: int) -> tuple[int, int]:
    """retorna (x=col, y=rank) onde rank in 1..6 (1 = base)."""
    row = v // 6
    col = v % 6
    rank = 6 - row
    return col, rank


def v_label(v: int) -> str:
    row = v // 6
    col = v % 6
    return f"{chr(ord('A') + col)}{6 - row}"


def cell_color(v: int) -> int:
    """retorna 0 (escura) ou 1 (clara). A1 e tipicamente escura no xadrez."""
    row = v // 6
    col = v % 6
    return (row + col) % 2


# ----------------------------------------------------------------------
# Carregar tudo
# ----------------------------------------------------------------------

def load_all():
    Forbidden = np.load(os.path.join(RES, "forbidden_cycles.npy"))
    Base_Ham = np.load(os.path.join(RES, "base_ham_rref.npy"))
    inc = np.load(os.path.join(RES, "incidence_vertex_edge.npy"))
    with open(os.path.join(HO, "edges_6x6.json")) as f:
        meta = json.load(f)
    edges_uv = [tuple(e) for e in meta["edges_uv"]]
    labels = meta["edge_labels"]
    with open(os.path.join(RES, "forbidden_cycles.json")) as f:
        fjson = json.load(f)
    return Forbidden, Base_Ham, inc, edges_uv, labels, fjson


# ----------------------------------------------------------------------
# 2.1  plotagem
# ----------------------------------------------------------------------

def plot_cycle(f: np.ndarray, edges_uv, idx: int, info_lines: list[str], out_path: str):
    fig, ax = plt.subplots(figsize=(6.5, 6.8))
    # tabuleiro (xadrez)
    for v in range(36):
        x, y = v_to_xy(v)
        c = cell_color(v)
        face = "#e9e1c8" if c else "#7d8a64"
        ax.add_patch(plt.Rectangle((x - 0.5, y - 0.5), 1, 1,
                                   facecolor=face, edgecolor="none", zorder=0))
    # arestas cinza claro (todas) - desenhar arestas fora de S como linhas finas
    gray_lines = []
    red_lines = []
    for ei, (u, v) in enumerate(edges_uv):
        xu, yu = v_to_xy(u)
        xv, yv = v_to_xy(v)
        seg = [(xu, yu), (xv, yv)]
        if f[ei]:
            red_lines.append(seg)
        else:
            gray_lines.append(seg)
    ax.add_collection(LineCollection(gray_lines, colors="#bbbbbb",
                                     linewidths=0.6, alpha=0.55, zorder=1))
    ax.add_collection(LineCollection(red_lines, colors="#c0212a",
                                     linewidths=2.4, alpha=0.95, zorder=3))
    # vertices
    for v in range(36):
        x, y = v_to_xy(v)
        # marcar diferente se vertice participa do ciclo
        deg_S = 0
        for ei, (a, b) in enumerate(edges_uv):
            if f[ei] and (a == v or b == v):
                deg_S += 1
        if deg_S > 0:
            ax.plot(x, y, marker="o", color="#222", markersize=7, zorder=4)
        ax.text(x, y - 0.36, v_label(v), ha="center", va="center",
                fontsize=7, color="#333", zorder=5)
    ax.set_xlim(-0.6, 5.6)
    ax.set_ylim(0.4, 6.6)
    ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(f"Forbidden cycle f_{idx}    |S|={int(f.sum())}", fontsize=11)
    txt = "\n".join(info_lines)
    fig.text(0.02, 0.02, txt, fontsize=8, family="monospace", va="bottom")
    fig.tight_layout()
    fig.savefig(out_path, dpi=160)
    plt.close(fig)


# ----------------------------------------------------------------------
# 2.2  topologia: componentes do subgrafo S
# ----------------------------------------------------------------------

def cycle_components(f: np.ndarray, edges_uv):
    """retorna lista de componentes conexas (cada uma = lista de vertices ordenados)."""
    V = 36
    adj = {v: [] for v in range(V)}
    for ei, (u, v) in enumerate(edges_uv):
        if f[ei]:
            adj[u].append((v, ei))
            adj[v].append((u, ei))
    seen = set()
    comps = []
    for start in range(V):
        if start in seen or not adj[start]:
            continue
        # passeio: como grau eh 2, eh um ciclo simples
        comp = []
        cur = start
        prev = -1
        first_edge = adj[cur][0][1]
        steps = 0
        while True:
            comp.append(cur)
            seen.add(cur)
            # proximo vizinho diferente de 'prev'
            nxts = [nb for nb in adj[cur] if nb[0] != prev or len(adj[cur]) == 1]
            # escolhe consistentemente
            nxt = None
            for (w, ei) in adj[cur]:
                if w != prev:
                    nxt = (w, ei)
                    break
            if nxt is None:
                break
            prev = cur
            cur = nxt[0]
            steps += 1
            if cur == start:
                break
            if steps > V * 2:
                break
        comps.append(comp)
    return comps


def vertex_degrees(f: np.ndarray, edges_uv):
    deg = np.zeros(36, dtype=int)
    for ei, (u, v) in enumerate(edges_uv):
        if f[ei]:
            deg[u] += 1; deg[v] += 1
    return deg


# ----------------------------------------------------------------------
# 2.3  simetria D_4 e winding
# ----------------------------------------------------------------------

def d4_action(perm_name: str):
    """Retorna funcao v -> v' que aplica a transformacao do D_4 ao tabuleiro 6x6.
    perm_name in {r90, r180, r270, fx, fy, fd, fa}."""
    def to_rc(v):
        return v // 6, v % 6
    def from_rc(r, c):
        return r * 6 + c

    def f(v):
        r, c = to_rc(v)
        if perm_name == "r90":   # 90 deg anti-hor: (r,c) -> (5-c, r)
            return from_rc(5 - c, r)
        if perm_name == "r180":
            return from_rc(5 - r, 5 - c)
        if perm_name == "r270":
            return from_rc(c, 5 - r)
        if perm_name == "fx":    # flip vertical (sobre eixo horizontal central)
            return from_rc(5 - r, c)
        if perm_name == "fy":    # flip horizontal (sobre eixo vertical central)
            return from_rc(r, 5 - c)
        if perm_name == "fd":    # diagonal A1-F6 (r,c) -> (c,r) (apos ajustes)
            return from_rc(c, r)
        if perm_name == "fa":    # antidiagonal
            return from_rc(5 - c, 5 - r)
        raise ValueError(perm_name)
    return f


def apply_perm_to_vec(f: np.ndarray, edges_uv, perm) -> np.ndarray:
    """Retorna vetor indicador transformado."""
    edge_idx = {}
    for i, (u, v) in enumerate(edges_uv):
        a, b = (u, v) if u < v else (v, u)
        edge_idx[(a, b)] = i
    out = np.zeros_like(f)
    for ei, (u, v) in enumerate(edges_uv):
        if f[ei]:
            pu, pv = perm(u), perm(v)
            key = (pu, pv) if pu < pv else (pv, pu)
            out[edge_idx[key]] = 1
    return out


def in_coset(vec: np.ndarray, base_ham_rref: np.ndarray, pivots, forbidden: np.ndarray) -> int:
    """Retorna i in {0,1,2,3..} indicando a classe (com 0 = Ham, 1..7 = cosets nao-trivais).
    Como Forbidden tem dim 3, ha 8 cosets em H_1/Ham. Vamos retornar como uma tripla.
    """
    # reduce vec by base_ham
    w = (vec.astype(np.uint8) & 1).copy()
    for r, pc in enumerate(pivots):
        if w[pc]:
            w ^= base_ham_rref[r]
    if not w.any():
        return (0, 0, 0)
    # agora w esta em forbidden_span. expressar w como combinacao de f_0,f_1,f_2
    # construir matriz: [forbidden_rref] e reduzir w
    # mais simples: tentar as 8 combinacoes
    F = forbidden
    for a in range(2):
        for b in range(2):
            for c in range(2):
                if (a, b, c) == (0, 0, 0):
                    continue
                combo = (a * F[0] ^ b * F[1] ^ c * F[2]) & 1
                test = (w ^ combo) & 1
                # reduzir test por base_ham
                w2 = test.copy()
                for r, pc in enumerate(pivots):
                    if w2[pc]:
                        w2 ^= base_ham_rref[r]
                if not w2.any():
                    return (a, b, c)
    return None  # nao devia chegar aqui


def winding_around_center(f: np.ndarray, edges_uv) -> int:
    """Winding number discreto do subgrafo S em torno do centro do tabuleiro.

    O centro 6x6 fica entre as casas C3, C4, D3, D4 (entre row 2/3 e col 2/3).
    Como o subgrafo eh um ciclo (ou uniao de ciclos) em GF(2), podemos contar
    quantas arestas de S cruzam o eixo vertical x=2.5 indo (cima->baixo) menos
    (baixo->cima). Para uniao de ciclos isso da winding total mod 2 == 0; o
    valor inteiro depende da orientacao escolhida. Como nao temos orientacao
    canonica, retornamos:
      - n_cross_horizontal_axis: arestas cruzando y=3.5
      - n_cross_vertical_axis:   arestas cruzando x=2.5
    Sao invariantes do conjunto de arestas (independem de orientacao).
    """
    n_h = 0; n_v = 0
    for ei, (u, v) in enumerate(edges_uv):
        if not f[ei]:
            continue
        xu, yu = v_to_xy(u)
        xv, yv = v_to_xy(v)
        # cruza y=3.5?
        if (yu - 3.5) * (yv - 3.5) < 0:
            n_h += 1
        if (xu - 2.5) * (xv - 2.5) < 0:
            n_v += 1
    return n_h, n_v


# ----------------------------------------------------------------------
# 2.4  relacao com obrigatorias e exclusoes
# ----------------------------------------------------------------------

def load_exclusions():
    out = {}
    for k, fname in [(2, "exclusions_pairs.json"),
                     (3, "exclusions_triples.json"),
                     (4, "exclusions_quads.json")]:
        path = os.path.join(HO, fname)
        if os.path.exists(path):
            with open(path) as fh:
                data = json.load(fh)
            out[k] = data
    return out


def mandatory_edges(T: np.ndarray):
    """arestas em todos os tours."""
    freq = T.mean(axis=0)
    return [int(i) for i in np.where(freq == 1.0)[0]]


def _extract_tuple_list(payload):
    """Aceita formatos:
       - lista direta: [[e1,e2,...], ...]
       - dict com chave 'minimal' / 'minimals' / 'tuples'
       Cada tupla pode ser lista de int (idx) ou de str (label).
    """
    if isinstance(payload, list):
        return payload
    if isinstance(payload, dict):
        for k in ("minimal", "minimals", "tuples", "exclusions"):
            if k in payload:
                return payload[k]
        # fallback: primeiro valor que e lista
        for v in payload.values():
            if isinstance(v, list):
                return v
    return []


def exclusion_contained(S_set, excl_payload, label_to_idx):
    """Retorna nº de tuplas excluidas (k) totalmente contidas em S_set (como conjunto de aresta_idx)."""
    tuples = _extract_tuple_list(excl_payload)
    count = 0
    total = 0
    for t in tuples:
        if isinstance(t, dict):
            t = t.get("edges") or t.get("tuple") or t.get("indices") or []
        if not t:
            continue
        total += 1
        # converte para idx
        if all(isinstance(x, int) for x in t):
            idxs = t
        else:
            try:
                idxs = [label_to_idx[x] for x in t]
            except KeyError:
                continue
        if set(idxs).issubset(S_set):
            count += 1
    return count, total


# ----------------------------------------------------------------------
# Pipeline
# ----------------------------------------------------------------------

def main():
    Forbidden, Base_Ham, inc, edges_uv, labels, fjson = load_all()
    T = np.load(os.path.join(HO, "incidence_matrix_6x6.npy"))
    excl = load_exclusions()
    label_to_idx = {lab: i for i, lab in enumerate(labels)}
    mand_idx = set(mandatory_edges(T))
    mand_labels = [labels[i] for i in sorted(mand_idx)]
    print(f"arestas obrigatorias (P=1) : {len(mand_idx)}   -> {mand_labels}")

    # pivots de Base_Ham para uso em testes de simetria
    pivots_ham = []
    for r in range(Base_Ham.shape[0]):
        idx = int(np.argmax(Base_Ham[r]))
        pivots_ham.append(idx)

    report = {"forbidden": []}

    for i, f in enumerate(Forbidden):
        print()
        print("=" * 64)
        print(f"Ciclo proibido f_{i}")
        print("=" * 64)
        active_idx = [int(j) for j in np.where(f == 1)[0]]
        S_set = set(active_idx)
        active_labels = [labels[j] for j in active_idx]
        print(f"|S| = {len(active_idx)}")
        print(f"arestas: {active_labels}")

        # graus
        deg = vertex_degrees(f, edges_uv)
        odd = np.where(deg % 2 == 1)[0]
        assert odd.size == 0, f"f_{i} tem vertice de grau impar! nao e um ciclo"
        nz = [int(v) for v in np.where(deg > 0)[0]]
        deg_counts = {int(d): int((deg == d).sum()) for d in sorted(set(deg.tolist())) if d > 0}
        print(f"#vertices em S: {len(nz)}  | distribuicao de grau: {deg_counts}")

        # componentes (ciclos simples que compoem f)
        comps = cycle_components(f, edges_uv)
        print(f"componentes (ciclos simples): {len(comps)}")
        for k, comp in enumerate(comps):
            print(f"  c_{k} ({len(comp)} vertices): {[v_label(v) for v in comp]}")

        # cobertura
        rows = sorted(set(v // 6 for v in nz))
        cols = sorted(set(v % 6 for v in nz))
        ranks = [6 - r for r in rows]
        col_letters = [chr(ord('A') + c) for c in cols]
        print(f"linhas (rank) tocadas : {ranks}")
        print(f"colunas tocadas       : {col_letters}")

        # simetria D_4 sobre o subespaco (classe!)
        # ou seja: aplicar perm em f e ver em que coset cai
        sym_table = {}
        for pname in ["r90", "r180", "r270", "fx", "fy", "fd", "fa"]:
            perm = d4_action(pname)
            fp = apply_perm_to_vec(f, edges_uv, perm)
            coset = in_coset(fp, Base_Ham, pivots_ham, Forbidden)
            sym_table[pname] = coset
        print(f"D_4 acao (classe em H_1/Ham, tripla = (a,b,c) em f_0+f_1+f_2):")
        for k, vcoset in sym_table.items():
            print(f"  {k:>4}: {vcoset}")

        # winding aproximado
        nh, nv = winding_around_center(f, edges_uv)
        print(f"cruzamentos do eixo horizontal central (y=3.5): {nh}")
        print(f"cruzamentos do eixo vertical central   (x=2.5): {nv}")

        # interacao com obrigatorias
        S_mand = sorted(S_set & mand_idx)
        print(f"arestas obrigatorias em S: {len(S_mand)}  -> {[labels[j] for j in S_mand]}")

        # interacao com exclusoes minimais
        excl_in_S = {}
        for k, payload in excl.items():
            c, total = exclusion_contained(S_set, payload, label_to_idx)
            excl_in_S[k] = (c, total)
            print(f"exclusoes minimais ordem {k} contidas em S: {c} / {total}")

        # plot
        info_lines = [
            f"|S| = {len(active_idx)}",
            f"componentes: {len(comps)}",
            f"  " + " | ".join(f"c_{k}={len(comp)}" for k, comp in enumerate(comps)),
            f"linhas: {ranks}",
            f"colunas: {col_letters}",
            f"cruz eixo horizontal y=3.5: {nh}",
            f"cruz eixo vertical   x=2.5: {nv}",
            f"obrigatorias em S: {len(S_mand)}",
            f"exclusoes contidas (pares/triplas/quadras): "
            f"{excl_in_S.get(2, (0,0))[0]} / {excl_in_S.get(3, (0,0))[0]} / {excl_in_S.get(4, (0,0))[0]}",
        ]
        out_png = os.path.join(PLOTS, f"forbidden_cycle_{i}.png")
        plot_cycle(f, edges_uv, i, info_lines, out_png)
        print(f"plot -> {out_png}")

        report["forbidden"].append({
            "index": i,
            "support_size": len(active_idx),
            "vertices_in_S": len(nz),
            "degree_distribution": deg_counts,
            "n_components": len(comps),
            "component_sizes": [len(c) for c in comps],
            "component_vertices": [[v_label(v) for v in comp] for comp in comps],
            "rows_touched": ranks,
            "cols_touched": col_letters,
            "horizontal_axis_crossings": nh,
            "vertical_axis_crossings": nv,
            "mandatory_in_S": [labels[j] for j in S_mand],
            "minimal_pairs_in_S": excl_in_S.get(2, (0, 0))[0],
            "minimal_triples_in_S": excl_in_S.get(3, (0, 0))[0],
            "minimal_quads_in_S": excl_in_S.get(4, (0, 0))[0],
            "d4_coset_action": {k: list(map(int, val)) for k, val in sym_table.items()},
        })

    # ---- analise conjunta D_4 sobre o subespaco Forbidden
    print()
    print("=" * 64)
    print("ACAO D_4 no QUOCIENTE H_1/Ham (subespaco F_2^3)")
    print("=" * 64)
    # Para cada perm, qual sao as imagens dos 3 geradores?
    # se a tripla (a,b,c) -> (a',b',c') for linear, ganhamos um endomorfismo.
    for pname in ["r90", "r180", "r270", "fx", "fy", "fd", "fa"]:
        perm = d4_action(pname)
        cols = []
        for fi in Forbidden:
            fp = apply_perm_to_vec(fi, edges_uv, perm)
            c = in_coset(fp, Base_Ham, pivots_ham, Forbidden)
            cols.append(c)
        mat = np.array(cols).T  # 3x3 (col j = imagem do f_j)
        # determinante mod 2 -> invertivel se 1
        det = int(round(np.linalg.det(mat))) % 2
        print(f"{pname:>4} -> matriz acao em F_2^3:")
        print("       " + str(mat.tolist()))
        print(f"       det mod 2 = {det}")
    out_json = os.path.join(RES, "geometry_summary.json")
    with open(out_json, "w") as fh:
        json.dump(report, fh, indent=2)
    print(f"\nresumo -> {out_json}")


if __name__ == "__main__":
    main()
