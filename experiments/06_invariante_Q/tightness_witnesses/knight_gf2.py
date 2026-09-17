#!/usr/bin/env python3
"""
knight_gf2.py
=============
Infraestrutura compartilhada da Fase 2 (mapa dos witnesses duais).

Fornece:
  - grafo do cavalo n x n, arestas, rotulos de xadrez, nivel de anel L
  - enumeracao exaustiva (n=6) / amostragem (n>=8) de tours fechados
  - algebra GF(2): rref, rank, reducao modulo subespaco
  - sindromes de aresta para teste de pertinencia em C_n^perp

Convencoes (iguais a deficit_theorem/verify_small_cases.py):
  vid(r, c) = r * n + c ; label = coluna 'A'+c, linha n-r
"""
from __future__ import annotations

import random
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np

MOVES = [(2, 1), (2, -1), (-2, 1), (-2, -1),
         (1, 2), (1, -2), (-1, 2), (-1, -2)]


# ── grafo ───────────────────────────────────────────────────────────

def vid(r: int, c: int, n: int) -> int:
    return r * n + c


def label(v: int, n: int) -> str:
    r, c = divmod(v, n)
    return f"{chr(ord('A') + c)}{n - r}"


def ring_level(v: int, n: int) -> int:
    """Nivel de anel L: 0 = borda externa, crescendo para o bulk."""
    r, c = divmod(v, n)
    return min(r, c, n - 1 - r, n - 1 - c)


def build_graph(n: int) -> Tuple[List[List[int]], List[Tuple[int, int]]]:
    """Retorna (adj, edges) com edges ordenado, cada aresta (u, v), u < v."""
    total = n * n
    adj: List[List[int]] = [[] for _ in range(total)]
    edge_set = set()
    for r in range(n):
        for c in range(n):
            v = vid(r, c, n)
            for dr, dc in MOVES:
                nr, nc = r + dr, c + dc
                if 0 <= nr < n and 0 <= nc < n:
                    u = vid(nr, nc, n)
                    adj[v].append(u)
                    if u > v:
                        edge_set.add((v, u))
    return adj, sorted(edge_set)


def edge_index(edges: Sequence[Tuple[int, int]]) -> Dict[Tuple[int, int], int]:
    return {e: i for i, e in enumerate(edges)}


def edge_label(e: Tuple[int, int], n: int) -> str:
    return f"{label(e[0], n)}-{label(e[1], n)}"


def edge_level(e: Tuple[int, int], n: int) -> int:
    """Nivel do anel de uma aresta = min dos niveis dos extremos."""
    return min(ring_level(e[0], n), ring_level(e[1], n))


def corners(n: int) -> List[int]:
    return [vid(0, 0, n), vid(0, n - 1, n), vid(n - 1, 0, n), vid(n - 1, n - 1, n)]


def mandatory_edges(n: int, adj, eidx) -> Dict[int, List[int]]:
    """Para cada canto (grau 2), os indices das suas 2 arestas obrigatorias."""
    out: Dict[int, List[int]] = {}
    for c in corners(n):
        nbrs = sorted(set(adj[c]))
        assert len(nbrs) == 2, f"canto {label(c, n)} tem grau {len(nbrs)}"
        out[c] = [eidx[(min(c, u), max(c, u))] for u in nbrs]
    return out


def boundary_matrix(edges, total) -> np.ndarray:
    B = np.zeros((total, len(edges)), dtype=np.uint8)
    for i, (u, v) in enumerate(edges):
        B[u, i] = 1
        B[v, i] = 1
    return B


# ── enumeracao de tours ─────────────────────────────────────────────

def enumerate_tours(n: int, limit: int | None = None,
                    seed: int | None = None, jitter: int = 0,
                    shallow_random: int = 0) -> List[List[int]]:
    """
    Tours hamiltonianos FECHADOS, como listas de indices de aresta.

    Todo tour passa pelo vertice 0 (canto), entao basta enumerar ciclos
    hamiltonianos a partir dele. O vertice 0 tem grau 2, logo AS DUAS arestas
    incidentes estao em todo tour: fixando path[1] = min(nbrs[0]) o outro
    vizinho e' forcosamente o extremo final, e cada tour nao-orientado sai
    exatamente uma vez. (Desempatar a posteriori com path[1] < cur seria
    errado: condena metade da arvore a nao produzir tour nenhum.)

    limit=None  -> exaustivo (viavel para n=6: 9862 tours)
    limit=K     -> para em K tours, com ordem de vizinhos aleatoria (n>=8)
    """
    adj, edges = build_graph(n)
    eidx = edge_index(edges)
    total = n * n
    rng = random.Random(seed)

    nbrs = [sorted(set(a)) for a in adj]
    if limit is not None:
        for a in nbrs:
            rng.shuffle(a)

    deg = [len(a) for a in nbrs]
    visited = bytearray(total)
    path: List[int] = [0]
    visited[0] = 1
    out: List[List[int]] = []

    # grau residual: usado na poda "vertice livre com <2 vizinhos livres"
    free_deg = deg[:]
    for u in nbrs[0]:
        free_deg[u] -= 1

    def prune(cur: int, depth: int) -> bool:
        """
        Poda em dois niveis. O grafo residual e' o induzido pelos vertices
        livres mais {cur, 0} (0 e' o ponto de fechamento do ciclo).

        (a) grau: todo vertice livre precisa de >= 2 saidas no residual;
        (b) conexidade: o residual precisa ser conexo -- sem isto o DFS
            passeia por semanas em n=8 (era o gargalo da versao anterior).
        """
        for v in range(total):
            if visited[v]:
                continue
            avail = 0
            for u in nbrs[v]:
                if not visited[u] or u == cur or u == 0:
                    avail += 1
                    if avail >= 2:
                        break
            if avail < 2:
                return True

        # (b) BFS a partir de cur sobre o residual
        stack = [cur]
        seen = bytearray(total)
        seen[cur] = 1
        reach = 0
        while stack:
            x = stack.pop()
            for u in nbrs[x]:
                if seen[u]:
                    continue
                if not visited[u]:
                    seen[u] = 1
                    reach += 1
                    stack.append(u)
                elif u == 0:
                    seen[u] = 1
        if not seen[0]:
            return True
        return reach != total - depth

    def dfs(cur: int, depth: int) -> bool:
        """Retorna True se atingiu o limite de tours."""
        if depth == total:
            if 0 in nbrs[cur]:
                out.append([eidx[(min(a, b), max(a, b))]
                            for a, b in zip(path, path[1:] + [0])])
                if limit is not None and len(out) >= limit:
                    return True
            return False

        if prune(cur, depth):
            return False

        # ordena por grau residual crescente (Warnsdorff)
        if depth == 1:
            cands = [min(nbrs[0])]          # quebra ida/volta (ver docstring)
        else:
            cands = [u for u in nbrs[cur] if not visited[u]]
        if limit is None:
            cands.sort(key=lambda u: free_deg[u])
        elif depth <= shallow_random:
            # Nos primeiros niveis, ordem UNIFORME. O Warnsdorff quase sempre
            # fecha um tour na primeira tentativa, entao as opcoes seguintes
            # nunca sao exploradas: no 8x8 isso deixava as duas arestas de
            # saida de C7 (= path[1] forcado) com frequencia ZERO, travando o
            # rank de Ham em 100 em vez de 102.
            rng.shuffle(cands)
        else:
            # Warnsdorff com ruido: sem o jitter os tours saem correlacionados.
            cands.sort(key=lambda u: (free_deg[u] + rng.randrange(jitter + 1),
                                      rng.random()))

        for u in cands:
            visited[u] = 1
            for w in nbrs[u]:
                free_deg[w] -= 1
            path.append(u)
            if dfs(u, depth + 1):
                return True
            path.pop()
            for w in nbrs[u]:
                free_deg[w] += 1
            visited[u] = 0
        return False

    dfs(0, 1)
    return out


def sample_tours(n: int, k: int, seed: int = 0,
                 restarts: int = 4000, jitter: int = 2,
                 shallow_random: int = 4) -> List[List[int]]:
    """
    Amostra ~k tours distintos com reinicios aleatorios.

    Diversidade e' o ponto critico: se o rank de Ham nao atingir beta_1 - 3,
    C_n^perp sai grande demais e o mapa de witnesses enche de artefatos.
    Por isso: muitos reinicios curtos + jitter no Warnsdorff.
    """
    seen = set()
    out: List[List[int]] = []
    per = max(1, k // restarts)
    for s in range(restarts):
        for t in enumerate_tours(n, limit=per, seed=seed * 100003 + s,
                                 jitter=jitter,
                                 shallow_random=shallow_random):
            key = frozenset(t)
            if key not in seen:
                seen.add(key)
                out.append(t)
        if len(out) >= k:
            break
    return out[:k]


# ── algebra GF(2) ───────────────────────────────────────────────────

def gf2_rref(M: np.ndarray) -> Tuple[np.ndarray, List[int], int]:
    A = M.copy().astype(np.uint8)
    rows, cols = A.shape
    pivots: List[int] = []
    r = 0
    for c in range(cols):
        if r >= rows:
            break
        pivot = None
        for rr in range(r, rows):
            if A[rr, c]:
                pivot = rr
                break
        if pivot is None:
            continue
        if pivot != r:
            A[[r, pivot]] = A[[pivot, r]]
        for rr in range(rows):
            if rr != r and A[rr, c]:
                A[rr] ^= A[r]
        pivots.append(c)
        r += 1
    return A[:r], pivots, r


def gf2_rank(M: np.ndarray) -> int:
    return gf2_rref(M)[2]


def rows_to_ints(rows: np.ndarray) -> List[int]:
    """Cada linha vira um int (bit i = coluna i)."""
    out = []
    for row in rows:
        x = 0
        for i in np.flatnonzero(row):
            x |= 1 << int(i)
        out.append(x)
    return out


def vec_to_int(support: Iterable[int]) -> int:
    x = 0
    for i in support:
        x ^= 1 << int(i)
    return x


def int_to_support(x: int) -> List[int]:
    out = []
    i = 0
    while x:
        if x & 1:
            out.append(i)
        x >>= 1
        i += 1
    return out


class Reducer:
    """Reducao de vetores (como int) modulo um subespaco em rref."""

    def __init__(self, basis_rows: np.ndarray):
        red, pivots, rank = gf2_rref(basis_rows)
        self.rows = rows_to_ints(red)
        self.pivots = pivots
        self.rank = rank

    def reduce(self, x: int) -> int:
        for p, row in zip(self.pivots, self.rows):
            if (x >> p) & 1:
                x ^= row
        return x

    def contains(self, x: int) -> bool:
        return self.reduce(x) == 0
