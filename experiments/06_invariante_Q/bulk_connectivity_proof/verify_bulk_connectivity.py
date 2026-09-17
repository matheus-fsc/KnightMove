#!/usr/bin/env python3
r"""
Verificação computacional do Lema de Conexidade do Bulk e das sub-afirmações
da prova INDUTIVA (passo n -> n+2).

Statement alvo:
    Para todo n >= 6, Bulk(n) := G_n[V_n \ Corners(n)] é conexo,
    onde G_n é o grafo do cavalo no tabuleiro n x n e Corners(n) são as 4
    casas de canto (cada uma de grau 2).

Esta script:
  (1) Verifica diretamente a conexidade de Bulk(n) por BFS, para n grande.
  (2) Verifica EXPLICITAMENTE cada sub-afirmação da prova indutiva:
        - Claim I  : a imersão de B = {1<=i,j<=n} no tabuleiro (n+2) induz
                     um subgrafo cujo conjunto de arestas, retirados os 4
                     cantos imersos, é EXATAMENTE Bulk(n) (invariância por
                     translação).
        - Claim II : cada canto imerso (1,1),(1,n),(n,1),(n,n) tem >=1 vizinho
                     de cavalo em B que NÃO é canto imerso (anexa-se ao bulk).
        - Claim IV : toda casa da moldura (anel externo de largura 1) não-canto
                     do tabuleiro (n+2) conecta-se a B pelo MOVIMENTO ESPECÍFICO
                     usado na prova; e esse movimento é um movimento de cavalo
                     válido que atersna em B.
"""

import time
from collections import deque

KNIGHT = [(1, 2), (2, 1), (2, -1), (1, -2),
          (-1, -2), (-2, -1), (-2, 1), (-1, 2)]


def knight_neighbors(n, i, j):
    out = []
    for di, dj in KNIGHT:
        a, b = i + di, j + dj
        if 0 <= a < n and 0 <= b < n:
            out.append((a, b))
    return out


def corners(n):
    return {(0, 0), (0, n - 1), (n - 1, 0), (n - 1, n - 1)}


def bulk_bfs(n):
    """BFS em Bulk(n) = G_n menos os 4 cantos. Retorna (alcançados, esperado,
    corners_isolados_bool)."""
    cset = corners(n)
    # vértice inicial não-canto garantido (interior)
    start = (2, 2)
    seen = {start}
    dq = deque([start])
    while dq:
        v = dq.popleft()
        for w in knight_neighbors(n, *v):
            if w in cset or w in seen:
                continue
            seen.add(w)
            dq.append(w)
    reached = len(seen)
    expected = n * n - 4
    # cantos isolados: cada canto tem todos os vizinhos... no Bulk eles foram
    # removidos do conjunto de vértices, então "isolado" significa: cada canto
    # tem grau 2 em G_n (logo após remover Mand fica isolado). Verificamos grau.
    corners_isolated = all(len(knight_neighbors(n, *c)) == 2 for c in cset)
    return reached, expected, corners_isolated


# ---------- Sub-afirmações da prova indutiva (passo n -> n+2) ----------

def embedded_corners(n):
    """Cantos do n-tabuleiro imersos no (n+2)-tabuleiro: deslocamento +1."""
    return {(1, 1), (1, n), (n, 1), (n, n)}


def interior_block(n):
    """B = {(i,j) : 1<=i,j<=n} dentro do tabuleiro (n+2)."""
    return {(i, j) for i in range(1, n + 1) for j in range(1, n + 1)}


def claim_I_embedding_edges(n):
    """A imersão preserva exatamente as arestas: arestas induzidas em B
    (no grafo do cavalo (n+2)) menos os 4 cantos imersos == arestas de Bulk(n).
    Comparamos os conjuntos de arestas (em coordenadas do n-tabuleiro)."""
    m = n + 2
    B = interior_block(n)
    ec = embedded_corners(n)
    # arestas induzidas em B\ec dentro de G_m, traduzidas para coords do n-board
    edges_embedded = set()
    for (i, j) in B:
        if (i, j) in ec:
            continue
        for (a, b) in knight_neighbors(m, i, j):
            if (a, b) in B and (a, b) not in ec:
                e = tuple(sorted([(i - 1, j - 1), (a - 1, b - 1)]))
                edges_embedded.add(e)
    # arestas de Bulk(n) diretamente
    cset = corners(n)
    edges_bulk = set()
    for i in range(n):
        for j in range(n):
            if (i, j) in cset:
                continue
            for (a, b) in knight_neighbors(n, i, j):
                if (a, b) in cset:
                    continue
                e = tuple(sorted([(i, j), (a, b)]))
                edges_bulk.add(e)
    return edges_embedded == edges_bulk


def claim_II_embedded_corner_attaches(n):
    """Cada canto imerso tem >=1 vizinho de cavalo em B que não é canto imerso."""
    m = n + 2
    B = interior_block(n)
    ec = embedded_corners(n)
    for c in ec:
        nbrs = [w for w in knight_neighbors(m, *c) if w in B and w not in ec]
        if len(nbrs) == 0:
            return False, c
    return True, None


def frame_target(n, cell):
    """Movimento ESPECÍFICO da prova: para casa de moldura do (n+2)-tabuleiro,
    devolve a casa-alvo em B (move 'duas para dentro, uma ao longo').
    Coordenadas no (n+2)-tabuleiro (0..n+1). B = 1..n."""
    i, j = cell
    m = n + 2
    # topo: i == 0
    if i == 0:
        k = j  # 1 <= k <= n
        return (2, k + 1) if k <= n - 1 else (2, k - 1)
    # base: i == n+1
    if i == m - 1:
        k = j
        return (n - 1, k + 1) if k <= n - 1 else (n - 1, k - 1)
    # esquerda: j == 0
    if j == 0:
        k = i
        return (k + 1, 2) if k <= n - 1 else (k - 1, 2)
    # direita: j == n+1
    if j == m - 1:
        k = i
        return (k + 1, n - 1) if k <= n - 1 else (k - 1, n - 1)
    raise ValueError("não é casa de moldura")


def is_knight_move(u, v):
    di, dj = abs(u[0] - v[0]), abs(u[1] - v[1])
    return (di, dj) in [(1, 2), (2, 1)]


def claim_IV_frame_connects(n):
    """Toda casa de moldura não-canto do (n+2)-tabuleiro conecta-se a B pelo
    movimento específico, que é (a) movimento de cavalo válido e (b) aterrissa
    em B. Retorna (ok, contraexemplo)."""
    m = n + 2
    B = interior_block(n)
    real_corners = {(0, 0), (0, m - 1), (m - 1, 0), (m - 1, m - 1)}
    frame = []
    for i in range(m):
        for j in range(m):
            if (i in (0, m - 1) or j in (0, m - 1)) and (i, j) not in real_corners:
                frame.append((i, j))
    for cell in frame:
        tgt = frame_target(n, cell)
        if not is_knight_move(cell, tgt):
            return False, ("não-cavalo", cell, tgt)
        if tgt not in B:
            return False, ("fora-de-B", cell, tgt)
    return True, None


def main():
    print("=== BULK CONNECTIVITY VERIFICATION ===")
    print(f"{'n':>3} | {'reached':>8} | {'expected':>8} | {'connected':>9} | "
          f"{'corners iso':>11} | {'time':>8}")
    ns = [6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 25, 30]
    all_conn = True
    for n in ns:
        t0 = time.perf_counter()
        reached, expected, ciso = bulk_bfs(n)
        dt = (time.perf_counter() - t0) * 1000
        conn = (reached == expected)
        all_conn = all_conn and conn and ciso
        print(f"{n:>3} | {reached:>8} | {expected:>8} | "
              f"{'YES' if conn else 'NO':>9} | {'YES' if ciso else 'NO':>11} | "
              f"{dt:>6.1f}ms")

    print()
    print("=== SUB-CLAIM VERIFICATION (inductive step n -> n+2) ===")
    # passo aplicado para cada n cujo n+2 queremos derivar; cobre as duas
    # paridades a partir das bases {6,7}.
    step_ns = [6, 7, 8, 9, 10, 12, 14, 16, 18, 23, 28]
    cI = cII = cIV = True
    for n in step_ns:
        okI = claim_I_embedding_edges(n)
        okII, cexII = claim_II_embedded_corner_attaches(n)
        okIV, cexIV = claim_IV_frame_connects(n)
        cI = cI and okI
        cII = cII and okII
        cIV = cIV and okIV
        flag = "OK" if (okI and okII and okIV) else "FAIL"
        extra = ""
        if not okII:
            extra += f"  II-cex={cexII}"
        if not okIV:
            extra += f"  IV-cex={cexIV}"
        print(f"  n={n:>2}->{n+2:<2}: Claim I (edges)={ok(okI)}  "
              f"Claim II (corner attach)={ok(okII)}  "
              f"Claim IV (frame connects)={ok(okIV)}  [{flag}]{extra}")

    print()
    print("Resumo sub-afirmações:")
    print(f"  Claim I  (imersão preserva arestas == Bulk(n)) : {PASS(cI)}")
    print(f"  Claim II (cantos imersos anexam ao bulk)        : {PASS(cII)}")
    print(f"  Claim IV (moldura conecta via movimento dado)   : {PASS(cIV)}")
    print()
    print(f"Conexidade direta (BFS) até n=30: {PASS(all_conn)}")
    print()
    overall = all_conn and cI and cII and cIV
    print(f"VEREDICTO GERAL: {PASS(overall)}")
    return 0 if overall else 1


def ok(b):
    return "PASS" if b else "FAIL"


def PASS(b):
    return "PASS" if b else "FAIL"


if __name__ == "__main__":
    raise SystemExit(main())
