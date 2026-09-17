#!/usr/bin/env python3
"""
flip_paths.py
=============
Decidir flip-realizabilidade de um hexagono por DECOMPOSICAO EM CAMINHOS, em
vez de busca hamiltoniana com arestas forcadas.

## A reformulacao

Se `tau` e `tau XOR C` sao ambos tours, escreva `F = tau \\ C`. Entao:

  - todo vertice de C tem grau 1 em F (uma aresta de C, uma fora);
  - todo outro vertice tem grau 2 em F;

logo **F e' uma uniao de 3 caminhos disjuntos que cobrem o tabuleiro, com
extremos exatamente nos 6 vertices de C**, e nenhum caminho PASSA por um
vertice de C.

Seja `sigma` o emparelhamento dos 6 vertices induzido pelos caminhos. Entao

    F ∪ m  e' um unico ciclo  <=>  sigma ∪ m  e' um 6-ciclo,

e o mesmo para `m'`. Portanto C e' flip-realizavel se e so se existe um tal
sistema de caminhos cujo `sigma` satisfaz **ambas** as condicoes.

Dos 15 emparelhamentos perfeitos de 6 pontos, apenas **4** satisfazem as duas
(8 satisfazem so' `m`). Isso da' tres ganhos sobre a busca anterior:

  1. as 6 arestas de C ficam proibidas, nao 3;
  2. nenhum caminho pode ATRAVESSAR um vertice de C -- entrar num deles fecha
     o caminho na hora, que e' uma poda forte e imediata;
  3. o teste do flip sai da folha e vira escolha de `sigma` na raiz: nao se
     gasta busca em ramos que so' serviriam a `m`.

A busca e' canonica (caminhos ordenados pelo extremo menor, percorridos do
menor para o maior extremo), logo sem redundancia de ordem/orientacao.
"""
from __future__ import annotations

from typing import Dict, List, Sequence, Tuple

from flip_graph import cycle_mask


class Budget(Exception):
    pass


def _matchings(vs: List[int]):
    if not vs:
        yield frozenset()
        return
    a = vs[0]
    for i in range(1, len(vs)):
        b = vs[i]
        for r in _matchings(vs[1:i] + vs[i + 1:]):
            yield r | {(a, b)}


def _is_6cycle(es, nodes) -> bool:
    nb: Dict[int, List[int]] = {v: [] for v in nodes}
    for a, b in es:
        nb[a].append(b)
        nb[b].append(a)
    if any(len(x) != 2 for x in nb.values()):
        return False
    s = nodes[0]
    prev, cur, cnt = -1, s, 0
    while True:
        x, y = nb[cur]
        nxt = x if x != prev else y
        prev, cur = cur, nxt
        cnt += 1
        if cur == s:
            break
        if cnt > len(nodes):
            return False
    return cnt == len(nodes)


def valid_sigmas(cyc: Sequence[int]) -> List[List[Tuple[int, int]]]:
    """Os emparelhamentos sigma dos 6 vertices de `cyc` (ordem ciclica) tais
    que sigma∪m e sigma∪m' sao ambos 6-ciclos. Sempre 4."""
    v = list(cyc)
    nodes = sorted(v)
    m = [(v[0], v[1]), (v[2], v[3]), (v[4], v[5])]
    mp = [(v[1], v[2]), (v[3], v[4]), (v[5], v[0])]
    out = []
    for s in _matchings(v):
        se = [tuple(e) for e in s]
        if _is_6cycle(se + m, nodes) and _is_6cycle(se + mp, nodes):
            out.append([(min(a, b), max(a, b)) for a, b in se])
    return out


def flip_realizable(V: int, adj, edges, eidx, cyc, budget: int) -> Tuple[bool, bool]:
    """(achou, exhausted). `exhausted=False` so' quando o orcamento estourou."""
    cm = cycle_mask(cyc, eidx)
    cverts = set(cyc)
    cedges = {(min(a, b), max(a, b))
              for a, b in zip(cyc, tuple(cyc[1:]) + (cyc[0],))}
    allowed = [[u for u in adj[v] if (min(v, u), max(v, u)) not in cedges]
               for v in range(V)]

    nodes_used = 0

    def try_sigma(sigma) -> bool:
        nonlocal nodes_used
        partner: Dict[int, int] = {}
        for a, b in sigma:
            partner[a] = b
            partner[b] = a
        used = [False] * V
        # caminhos em ordem canonica: pelo extremo MENOR, crescente
        starts = sorted({min(a, b) for a, b in sigma})
        found = False

        def viable(cur: int, k: int) -> bool:
            """Poda CONSERVADORA: so' rejeita o que nao pode dar solucao.

            (a) alcancabilidade -- todo vertice livre tem de ser alcancavel a
                partir da cabeca do caminho atual OU de algum extremo de um
                caminho ainda por percorrer (os dois extremos contam), andando
                so' por vertices livres e sem ATRAVESSAR vertice de C;
            (b) grau -- todo vertice livre fora de C precisa de >= 2 vizinhos
                disponiveis (livres, ou a propria cabeca `cur`); todo extremo
                ainda livre precisa de >= 1.
            """
            free = [v for v in range(V) if not used[v]]
            if not free:
                return True
            srcs = [cur]
            for j in range(k + 1, len(starts)):
                s2 = starts[j]
                srcs.append(s2)
                srcs.append(partner[s2])
            seen = set()
            stack = []
            for s2 in srcs:
                if s2 not in seen:
                    seen.add(s2)
                    stack.append(s2)
            while stack:
                x = stack.pop()
                for u in allowed[x]:
                    if used[u] or u in seen:
                        continue
                    seen.add(u)
                    if u not in cverts:
                        stack.append(u)
            for v in free:
                if v not in seen:
                    return False
            for v in free:
                need = 1 if v in cverts else 2
                c = 0
                for u in allowed[v]:
                    if (not used[u]) or u == cur:
                        c += 1
                        if c >= need:
                            break
                if c < need:
                    return False
            return True

        def walk(cur: int, target: int, k: int) -> None:
            """Estende o caminho k (que deve terminar em `target`)."""
            nonlocal found, nodes_used
            if found:
                return
            nodes_used += 1
            if nodes_used > budget:
                raise Budget
            if cur == target:
                if k + 1 == len(starts):
                    if all(used):
                        found = True
                    return
                nxt_start = starts[k + 1]
                used[nxt_start] = True
                walk(nxt_start, partner[nxt_start], k + 1)
                used[nxt_start] = False
                return
            for u in allowed[cur]:
                if used[u]:
                    continue
                if u in cverts and u != target:
                    continue          # caminho nao atravessa vertice de C
                used[u] = True
                if u == target or viable(u, k):
                    walk(u, target, k)
                used[u] = False
                if found:
                    return

        s0 = starts[0]
        used[s0] = True
        walk(s0, partner[s0], 0)
        return found

    try:
        for sigma in valid_sigmas(cyc):
            if try_sigma(sigma):
                return True, True
        return False, True
    except Budget:
        return False, False
