"""
d4_orbits.py
============
Grupo D₄ (8 simetrias do quadrado) sobre o tabuleiro n×n do cavalo.

Aplicado a:
  - vértices  (r, c)
  - arestas   ((r1,c1), (r2,c2))   sempre tuple-sorted
  - pares de arestas (edge_a, edge_b)  não-ordenado, canonicalizado lex-min

Fronteira de I/O: labels em notação chess "A8-C7" (ou "A6-..." em 6×6) usados pelo engine.

Tamanho do tabuleiro: configurável via `set_board(n)`. Default 8 (compat 8×8 legado).
Todas as funções consultam BOARD/N dinamicamente, então basta um set_board(6) no início
do script para alternar para 6×6.
"""

BOARD = 8
N = BOARD - 1  # índice máximo


def set_board(n):
    """Reconfigura o tamanho do tabuleiro (afeta encoding de labels e D₄)."""
    global BOARD, N
    BOARD = int(n)
    N = BOARD - 1


# ── conversão label ↔ tupla ──────────────────────────────────────────

def label_to_node(s):
    """'A8' -> (0, 0).  Convenção: rank 8 = row 0 (topo)."""
    c = ord(s[0]) - ord('A')
    r = BOARD - int(s[1:])
    return (r, c)

def node_to_label(rc):
    r, c = rc
    return chr(ord('A') + c) + str(BOARD - r)

def label_to_edge(s):
    """'A8-C7' -> ((0,0), (1,2))  já tuple-sorted."""
    a, b = s.split("-")
    u, v = label_to_node(a), label_to_node(b)
    return (u, v) if u <= v else (v, u)

def edge_to_label(e):
    """((0,0), (1,2)) -> 'A8-C7'  (canonical, tuple-sorted)."""
    u, v = e
    if u > v:
        u, v = v, u
    return f"{node_to_label(u)}-{node_to_label(v)}"


# ── ações de D₄ ──────────────────────────────────────────────────────

def d4_apply_node(rc, t):
    """Uma das 8 transformações de D₄ aplicada a (r, c)."""
    r, c = rc
    if t == 0: return (r,     c    )   # identidade
    if t == 1: return (c,     N - r)   # rot 90° CW
    if t == 2: return (N - r, N - c)   # rot 180°
    if t == 3: return (N - c, r    )   # rot 270° CW
    if t == 4: return (r,     N - c)   # reflexão horizontal
    if t == 5: return (N - r, c    )   # reflexão vertical
    if t == 6: return (c,     r    )   # transposição (diag principal)
    if t == 7: return (N - c, N - r)   # anti-transposição (diag secundária)
    raise ValueError(f"t fora de [0,7]: {t}")

def d4_apply_edge(e, t):
    """Aplica D₄ a uma aresta; retorna tuple-sorted."""
    u, v = e
    u2 = d4_apply_node(u, t)
    v2 = d4_apply_node(v, t)
    return (u2, v2) if u2 <= v2 else (v2, u2)


# ── órbitas ───────────────────────────────────────────────────────────

def edge_orbit(e):
    """Conjunto das 8 imagens D₄ da aresta e (size pode ser <8 se autoinvariante)."""
    return {d4_apply_edge(e, t) for t in range(8)}

def edge_orbits(edges_iterable):
    """Particiona arestas em órbitas D₄. Retorna lista ordenada de listas ordenadas."""
    edge_set = set(edges_iterable)
    seen = set()
    orbits = []
    for e in sorted(edge_set):
        if e in seen:
            continue
        orb = edge_orbit(e) & edge_set
        orbits.append(sorted(orb))
        seen |= orb
    return orbits

def edge_to_orbit_id(orbits):
    m = {}
    for i, orb in enumerate(orbits):
        for e in orb:
            m[e] = i
    return m

def pair_d4_canonical(edge_a, edge_b):
    """
    Representante canônico (lex-min) da órbita D₄ do par não-ordenado
    {edge_a, edge_b}. Usado para colapsar correlações em invariantes
    topológicos distintos.
    """
    images = set()
    for t in range(8):
        a = d4_apply_edge(edge_a, t)
        b = d4_apply_edge(edge_b, t)
        images.add(tuple(sorted([a, b])))
    return min(images)


# ── smoke test ────────────────────────────────────────────────────────

if __name__ == "__main__":
    # Reconstrói as 168 arestas do grafo 8×8 só para conferir contagens.
    deltas = [(1,2),(2,1),(-1,2),(-2,1),(1,-2),(2,-1),(-1,-2),(-2,-1)]
    edges = set()
    for r in range(BOARD):
        for c in range(BOARD):
            for dr, dc in deltas:
                nr, nc = r+dr, c+dc
                if 0 <= nr < BOARD and 0 <= nc < BOARD:
                    a, b = (r, c), (nr, nc)
                    edges.add((a, b) if a <= b else (b, a))

    print(f"Arestas no grafo do cavalo 8×8: {len(edges)}")
    orbits = edge_orbits(edges)
    print(f"Órbitas D₄ de arestas:          {len(orbits)}")
    sizes = {}
    for o in orbits:
        sizes[len(o)] = sizes.get(len(o), 0) + 1
    for sz, cnt in sorted(sizes.items()):
        print(f"  {cnt} órbita(s) de tamanho {sz}  →  {sz*cnt} arestas")
    print()
    print("Representantes (lex-min) das primeiras órbitas:")
    for i, o in enumerate(orbits[:5]):
        print(f"  O{i:02d} (|{len(o)}|): {edge_to_label(o[0])}")

    # Conferência de pares: a aresta-do-canto e sua imagem D₄ formam órbita-de-par
    e_corner = ((0, 0), (1, 2))
    e_other  = ((0, 0), (2, 1))   # outra opção do canto
    canon = pair_d4_canonical(e_corner, e_other)
    print()
    print(f"Par canônico de {edge_to_label(e_corner)} ↔ {edge_to_label(e_other)}:")
    print(f"  {edge_to_label(canon[0])} ↔ {edge_to_label(canon[1])}")
