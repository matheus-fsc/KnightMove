"""d4_symmetry.py — Simetria diedral D4 pura (sem dependência do motor de busca).

Representa o grupo diedral D4 (8 simetrias do quadrado) agindo sobre um
tabuleiro n×n. Vértice v = row*n + col (mesma convenção de knight_tours.py).

As 8 transformações de coordenada (r, c) -> (r', c'):
  id        (r, c)
  rot90     (c, n-1-r)
  rot180    (n-1-r, n-1-c)
  rot270    (n-1-c, r)
  flip_h    (r, n-1-c)        reflexão no eixo vertical (espelha colunas)
  flip_v    (n-1-r, c)        reflexão no eixo horizontal (espelha linhas)
  diag      (c, r)            reflexão na diagonal principal (transposta)
  anti      (n-1-c, n-1-r)    reflexão na diagonal secundária

Interface pública:
  D4Element / d4_elements()                    — as 8 transformações (parametrizadas por n)
  d4_group(n) -> list[BoundD4]                  — grupo "ligado" a n (callables sobre vértices/arestas)
  orbit_of_vertex(v, n) -> frozenset[int]
  canonical_vertex_representatives(n) -> list[int]
  stabilizer(fixed_edges_1, fixed_edges_0, n) -> list[BoundD4]
  filter_candidates_by_orbit(candidates, stabilizer_group) -> list[edge]
"""
from __future__ import annotations


# ---------------------------------------------------------------------------
# As 8 fórmulas de coordenada de D4
# ---------------------------------------------------------------------------

def _f_id(r, c, n):     return (r, c)
def _f_rot90(r, c, n):  return (c, n - 1 - r)
def _f_rot180(r, c, n): return (n - 1 - r, n - 1 - c)
def _f_rot270(r, c, n): return (n - 1 - c, r)
def _f_flip_h(r, c, n): return (r, n - 1 - c)
def _f_flip_v(r, c, n): return (n - 1 - r, c)
def _f_diag(r, c, n):   return (c, r)
def _f_anti(r, c, n):   return (n - 1 - c, n - 1 - r)


_D4_FORMULAS = [
    ("id", _f_id),
    ("rot90", _f_rot90),
    ("rot180", _f_rot180),
    ("rot270", _f_rot270),
    ("flip_h", _f_flip_h),
    ("flip_v", _f_flip_v),
    ("diag", _f_diag),
    ("anti", _f_anti),
]


class D4Element:
    """Uma das 8 transformações de D4, parametrizada por n em cada chamada."""

    __slots__ = ("name", "_fn")

    def __init__(self, name, fn):
        self.name = name
        self._fn = fn

    def cell(self, r, c, n):
        return self._fn(r, c, n)

    def vertex(self, v, n):
        r, c = divmod(int(v), n)
        r2, c2 = self._fn(r, c, n)
        return r2 * n + c2

    def edge(self, edge, n):
        u, w = edge
        return (self.vertex(u, n), self.vertex(w, n))

    def __repr__(self):
        return f"D4Element({self.name})"


def d4_elements():
    """Retorna as 8 transformações de D4 (não ligadas a um n específico)."""
    return [D4Element(name, fn) for name, fn in _D4_FORMULAS]


class BoundD4:
    """Transformação de D4 já ligada a um tabuleiro n×n.

    Callable sobre vértices; também mapeia arestas. Não precisa mais de n.
    """

    __slots__ = ("name", "n", "_elem")

    def __init__(self, elem: D4Element, n: int):
        self.name = elem.name
        self.n = n
        self._elem = elem

    def __call__(self, v):
        return self._elem.vertex(v, self.n)

    def vertex(self, v):
        return self._elem.vertex(v, self.n)

    def edge(self, edge):
        return self._elem.edge(edge, self.n)

    def __repr__(self):
        return f"BoundD4({self.name}, n={self.n})"


def d4_group(n):
    """Grupo D4 ligado ao tabuleiro n×n (8 BoundD4)."""
    return [BoundD4(e, n) for e in d4_elements()]


# ---------------------------------------------------------------------------
# Órbitas de vértices
# ---------------------------------------------------------------------------

def orbit_of_vertex(v, n):
    """Órbita do vértice v sob D4: frozenset dos vértices-imagem."""
    return frozenset(e.vertex(v, n) for e in d4_elements())


def canonical_vertex_representatives(n):
    """Um representante por órbita de vértice sob D4 (o de menor índice).

    Retorna a lista ordenada de representantes.
    Para n=8 deve devolver 10 órbitas (4 de tamanho 4 nas diagonais +
    6 de tamanho 8 fora delas = 64 vértices).
    """
    reps = []
    for v in range(n * n):
        if v == min(orbit_of_vertex(v, n)):
            reps.append(v)
    return reps


# ---------------------------------------------------------------------------
# Estabilizador residual e filtro de candidatos
# ---------------------------------------------------------------------------

def _norm_edge(edge):
    """Forma canônica de uma aresta não-orientada: (min, max)."""
    u, w = edge
    return (u, w) if u <= w else (w, u)


def _norm_edge_set(edges):
    return frozenset(_norm_edge(e) for e in edges)


def stabilizer(fixed_edges_1, fixed_edges_0, n):
    """Subgrupo de D4 que preserva SIMULTANEAMENTE E1 e E0 (como conjuntos).

    fixed_edges_1: iterável de arestas (pares de vértices) fixadas em 1.
    fixed_edges_0: iterável de arestas fixadas em 0.
    Retorna a lista de BoundD4 g tais que g(E1)=E1 e g(E0)=E0.
    (Sempre contém a identidade.)
    """
    e1 = _norm_edge_set(fixed_edges_1)
    e0 = _norm_edge_set(fixed_edges_0)
    out = []
    for g in d4_group(n):
        g_e1 = frozenset(_norm_edge(g.edge(e)) for e in e1)
        if g_e1 != e1:
            continue
        g_e0 = frozenset(_norm_edge(g.edge(e)) for e in e0)
        if g_e0 != e0:
            continue
        out.append(g)
    return out


def filter_candidates_by_orbit(candidates, stabilizer_group):
    """Reduz `candidates` a um representante por órbita sob o subgrupo dado.

    candidates: lista de arestas (pares de vértices).
    stabilizer_group: lista de BoundD4 (o estabilizador residual — NÃO o D4
        inteiro).
    Retorna, na ordem de entrada, o primeiro candidato de cada órbita.
    """
    covered = set()
    reps = []
    for c in candidates:
        key = _norm_edge(c)
        if key in covered:
            continue
        orbit = {_norm_edge(g.edge(c)) for g in stabilizer_group}
        covered |= orbit
        reps.append(c)
    return reps
