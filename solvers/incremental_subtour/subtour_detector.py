"""
Detector incremental de sub-ciclos usando Union-Find.

Mantido em sincronia com o estado de arestas fixadas em 1 durante o
backtracking. Detecta SUB-CICLOS ASSIM QUE FORMADOS — antes da folha.

Status retornados em fix(e, 1):
  OK              — aresta absorvida, sem ciclo fechado
  SUBTOUR         — ciclo fechado de comprimento < n² → poda
  COMPLETE_TOUR   — ciclo fechado de comprimento = n² → candidato a tour
  CONTRADICTION   — algum vértice atingiria grau > 2

Quando uma aresta é fixada em 0, o detector ignora (só monitora =1).

Duas variantes de undo:
  - SubtourDetectorRollback: pilha de operações revertidas em rollback()
  - SubtourDetectorCopy: copia arrays a cada nó (mais simples, mais lento)
Ambas implementam a mesma API; o backtracking escolhe via flag.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np


OK = "OK"
SUBTOUR = "SUBTOUR"
COMPLETE_TOUR = "COMPLETE_TOUR"
CONTRADICTION = "CONTRADICTION"


@dataclass
class _OpDeg:
    """Registro: degree[v] foi incrementado em delta (sempre +1 em fix)."""
    v: int
    delta: int


@dataclass
class _OpUnion:
    """Registro: realizamos union (child, parent), antigo parent[child]=old_p,
    antigo rank[parent]=old_r, antigo size[parent]=old_s,
    antigo closed_count delta=closed_delta."""
    child: int
    parent: int
    old_parent_of_child: int
    old_rank_of_parent: int
    old_size_of_parent: int
    closed_delta_revert: int  # quanto subtrair de n_closed para reverter


@dataclass
class _OpCloseSelf:
    """Registro: aresta fechou ciclo no mesmo componente (não há union).
    Não muda parent/rank/size — só degrees e closed_count."""
    closed_delta_revert: int


class SubtourDetectorRollback:
    """
    Union-Find com rollback explícito. Não usa path compression
    (incompatível com rollback simples).

    edges_uv: lista de (u, v) por idx de aresta.
    n_vertices: número total de vértices.
    """

    def __init__(self, edges_uv, n_vertices: int):
        self.edges = edges_uv
        self.V = n_vertices
        self.parent = np.arange(n_vertices, dtype=np.int32)
        self.rank = np.zeros(n_vertices, dtype=np.int8)
        self.size = np.ones(n_vertices, dtype=np.int32)
        self.degree = np.zeros(n_vertices, dtype=np.int8)
        # closed_count_of_root: para cada raiz r, quantos vértices do
        # componente têm degree 2 (degree[v] == 2). Mantido aproximado:
        # rastreamos n_closed_components (componentes inteiramente fechados).
        # Para detectar fechamento incremental, mantemos por componente
        # n_deg2 = nº de vértices com degree 2 no componente.
        self.n_deg2_of_root = np.zeros(n_vertices, dtype=np.int32)
        self.ops: list = []

    # --- find sem path compression (necessário para rollback) ---
    def find(self, v: int) -> int:
        while self.parent[v] != v:
            v = self.parent[v]
        return int(v)

    def checkpoint(self) -> int:
        return len(self.ops)

    def rollback(self, cp: int) -> None:
        while len(self.ops) > cp:
            op = self.ops.pop()
            if isinstance(op, _OpDeg):
                v = op.v
                # antes de reverter o degree, reverter o efeito em n_deg2
                # do componente. Como rollback ocorre LIFO em relação aos
                # demais ops, o estado parent/size já foi (ou será)
                # consistentemente revertido pelos outros ops.
                # n_deg2 do root atual:
                r = self.find(v)
                # se degree[v] era 2 e agora vai voltar a (2 - delta)
                novo = int(self.degree[v]) - op.delta
                if self.degree[v] == 2 and novo != 2:
                    self.n_deg2_of_root[r] -= 1
                elif self.degree[v] != 2 and novo == 2:
                    # caso degenerado (não ocorre em fix(+1))
                    self.n_deg2_of_root[r] += 1
                self.degree[v] = np.int8(novo)
            elif isinstance(op, _OpUnion):
                # desfazer union: parent[child] volta ao próprio child;
                # parent (raiz antiga) tem size/rank/n_deg2 restaurados
                self.parent[op.child] = op.old_parent_of_child
                self.rank[op.parent] = op.old_rank_of_parent
                self.size[op.parent] = op.old_size_of_parent
                # n_deg2_of_root[child] já estava em n_deg2_of_root[child]
                # antes do union; foi adicionado ao parent. Recomputar a
                # partir de closed_delta_revert:
                self.n_deg2_of_root[op.parent] -= op.closed_delta_revert
                self.n_deg2_of_root[op.child] = op.closed_delta_revert
            elif isinstance(op, _OpCloseSelf):
                # nada de parent/size/rank mudou
                pass
            else:  # pragma: no cover
                raise RuntimeError(f"op desconhecido: {op}")

    def fix(self, e_idx: int) -> str:
        """Absorve aresta fixada em 1. Retorna OK / SUBTOUR / COMPLETE / CONTRA."""
        u, v = self.edges[e_idx]
        # 1) checar grau prospectivo
        if self.degree[u] >= 2 or self.degree[v] >= 2:
            return CONTRADICTION

        # 2) incrementar degrees (com log)
        # antes de mudar, lembrar estado para incrementar n_deg2:
        ru = self.find(u)
        rv = self.find(v)

        # incrementa u
        old_du = int(self.degree[u])
        self.degree[u] = np.int8(old_du + 1)
        if old_du + 1 == 2:
            self.n_deg2_of_root[ru] += 1
        self.ops.append(_OpDeg(v=u, delta=1))

        # incrementa v
        old_dv = int(self.degree[v])
        self.degree[v] = np.int8(old_dv + 1)
        # cuidado: se ru == rv, n_deg2_of_root[ru] cobre v também
        if old_dv + 1 == 2:
            self.n_deg2_of_root[rv] += 1
        self.ops.append(_OpDeg(v=v, delta=1))

        # 3) union (se necessário) e detecção
        if ru == rv:
            # adicionar aresta dentro do mesmo componente FECHA um ciclo
            # cobrindo todos os 'size[ru]' vértices se e somente se o
            # subgrafo for um ciclo simples sobre eles (todos com grau 2).
            self.ops.append(_OpCloseSelf(closed_delta_revert=0))
            size_c = int(self.size[ru])
            n_deg2 = int(self.n_deg2_of_root[ru])
            if n_deg2 == size_c:
                if size_c == self.V:
                    return COMPLETE_TOUR
                return SUBTOUR
            # caso degenerado: aresta entre 2 nós já no mesmo componente
            # mas o componente não é um ciclo simples — significa que
            # algum vértice teria grau >2 (já bloqueado acima) ou o
            # componente tem desvios. Mantemos OK; folha cuidará.
            return OK

        # union by rank
        if self.rank[ru] < self.rank[rv]:
            ru, rv = rv, ru
        child, parent = rv, ru
        # Antes da fusão, lembrar antigo n_deg2_of_root[child] para
        # rollback. Note: durante a fusão, contagens de degree-2 do child
        # serão somadas ao parent.
        n_deg2_child = int(self.n_deg2_of_root[child])
        old_parent_of_child = int(self.parent[child])
        old_rank_parent = int(self.rank[parent])
        old_size_parent = int(self.size[parent])

        self.parent[child] = parent
        self.size[parent] += self.size[child]
        if self.rank[ru] == self.rank[rv]:
            self.rank[parent] += 1
        self.n_deg2_of_root[parent] += n_deg2_child
        self.n_deg2_of_root[child] = 0  # zera (raiz inativa)

        self.ops.append(_OpUnion(
            child=child, parent=parent,
            old_parent_of_child=old_parent_of_child,
            old_rank_of_parent=old_rank_parent,
            old_size_of_parent=old_size_parent,
            closed_delta_revert=n_deg2_child,
        ))

        # Após união, verificar se o novo componente está FECHADO:
        size_c = int(self.size[parent])
        n_deg2 = int(self.n_deg2_of_root[parent])
        if n_deg2 == size_c:
            if size_c == self.V:
                return COMPLETE_TOUR
            return SUBTOUR
        return OK


class SubtourDetectorCopy:
    """Variante 'copy-on-branch': sem rollback, faz cópia do estado.

    Mesma API que SubtourDetectorRollback. checkpoint() retorna um snapshot
    (tupla de arrays). rollback(cp) restaura tudo.
    Mais simples; arrays pequenos (V=100 → ~400 B por snapshot).
    """

    def __init__(self, edges_uv, n_vertices: int):
        self.edges = edges_uv
        self.V = n_vertices
        self.parent = np.arange(n_vertices, dtype=np.int32)
        self.rank = np.zeros(n_vertices, dtype=np.int8)
        self.size = np.ones(n_vertices, dtype=np.int32)
        self.degree = np.zeros(n_vertices, dtype=np.int8)
        self.n_deg2_of_root = np.zeros(n_vertices, dtype=np.int32)

    def find(self, v: int) -> int:
        # path compression OK aqui — rollback é por cópia completa
        root = v
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[v] != root:
            self.parent[v], v = root, self.parent[v]
        return int(root)

    def checkpoint(self):
        return (
            self.parent.copy(),
            self.rank.copy(),
            self.size.copy(),
            self.degree.copy(),
            self.n_deg2_of_root.copy(),
        )

    def rollback(self, cp) -> None:
        self.parent, self.rank, self.size, self.degree, self.n_deg2_of_root = (
            cp[0].copy(), cp[1].copy(), cp[2].copy(), cp[3].copy(), cp[4].copy()
        )

    def fix(self, e_idx: int) -> str:
        u, v = self.edges[e_idx]
        if self.degree[u] >= 2 or self.degree[v] >= 2:
            return CONTRADICTION

        ru = self.find(u)
        rv = self.find(v)

        old_du = int(self.degree[u])
        self.degree[u] = np.int8(old_du + 1)
        if old_du + 1 == 2:
            self.n_deg2_of_root[ru] += 1

        old_dv = int(self.degree[v])
        self.degree[v] = np.int8(old_dv + 1)
        # Se v compartilha componente com u (ru==rv) e ru já incrementado:
        # mas n_deg2_of_root[rv] == n_deg2_of_root[ru] (mesmo array slot)
        if old_dv + 1 == 2:
            self.n_deg2_of_root[rv] += 1

        if ru == rv:
            size_c = int(self.size[ru])
            n_deg2 = int(self.n_deg2_of_root[ru])
            if n_deg2 == size_c:
                if size_c == self.V:
                    return COMPLETE_TOUR
                return SUBTOUR
            return OK

        # union by rank
        if self.rank[ru] < self.rank[rv]:
            ru, rv = rv, ru
        child, parent = rv, ru
        n_deg2_child = int(self.n_deg2_of_root[child])
        self.parent[child] = parent
        self.size[parent] += self.size[child]
        if self.rank[ru] == self.rank[rv]:
            self.rank[parent] += 1
        self.n_deg2_of_root[parent] += n_deg2_child
        self.n_deg2_of_root[child] = 0

        size_c = int(self.size[parent])
        n_deg2 = int(self.n_deg2_of_root[parent])
        if n_deg2 == size_c:
            if size_c == self.V:
                return COMPLETE_TOUR
            return SUBTOUR
        return OK
