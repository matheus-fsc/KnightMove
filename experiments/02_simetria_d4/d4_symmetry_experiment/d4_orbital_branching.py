"""d4_orbital_branching.py — Enumeração de tours com redução por simetria D4.

Reusa o motor de knight_tours.py (SEM modificá-lo): build_graph, State,
fix_and_propagate, _choose_next_edge, _propagate_initial, _is_complete_tour
e as constantes de status. Importa d4_symmetry.py para a álgebra de D4.

Duas enumerações COMPLETAS (K=infinito), compartilhando a mesma propagação:

  enumerate_tours(n, use_symmetry=False)  -> conta cada tour fechado uma vez.
  enumerate_tours(n, use_symmetry=True)   -> enumera UM representante por
      D4-órbita (poda lex-leader canônica) e reconstrói a contagem total via
      soma dos tamanhos de órbita  Σ 8/|Stab(S)|  (órbita-estabilizador).

Correção da contagem sob simetria (ponto crítico):
  Cada tour = conjunto S de arestas ativas. D4 age sobre S (o grafo do
  cavalo é D4-invariante). S é *canônico* sse o bitstring de arestas (na
  ordem de índice de aresta) é o lex-mínimo da sua órbita. A poda mantém,
  a cada nó, o conjunto de simetrias "ainda vivas" (o estabilizador residual
  do prefixo decidido) e poda o nó assim que alguma simetria PROVA que ele
  não pode completar num canônico. As folhas sobreviventes são exatamente
  uma por órbita; o total é Σ_{reps} 8/|Stab(S)|.
  A poda é *sound* (nunca elimina o canônico de uma órbita); a corretude da
  contagem NÃO depende da força da poda, apenas do teste canônico na folha.
"""
from __future__ import annotations

# --- bootstrap: torna core/ importavel a partir de qualquer subpasta ---
import sys as _sys, pathlib as _pathlib
for _anc in _pathlib.Path(__file__).resolve().parents:
    if (_anc / "core").is_dir() and (_anc / "pyproject.toml").exists():
        _sys.path.insert(0, str(_anc / "core"))
        break
# --- fim do bootstrap ---

import os
import sys
import time

import numpy as np

# --- importa o motor existente sem copiá-lo ---
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

import knight_tours as kt  # noqa: E402
from knight_tours import (  # noqa: E402
    build_graph, State, fix_and_propagate, _choose_next_edge,
    _propagate_initial, _is_complete_tour,
    FREE, ACTIVE, INACTIVE, OK, SUBTOUR, CONTRADICTION, COMPLETE_TOUR,
)
import d4_symmetry as d4  # noqa: E402


# ---------------------------------------------------------------------------
# Permutações de arestas induzidas por D4
# ---------------------------------------------------------------------------

def build_edge_perms(ctx):
    """Constrói, para cada g em D4, a permutação sobre ÍNDICES de aresta.

    Retorna dict com:
      gperm    : int32 (8, E)  gperm[k][e] = índice da aresta g_k(aresta e)
      gperm_inv: int32 (8, E)  inversa de gperm[k]
      names    : lista de nomes dos 8 elementos (identidade em k=0)
    """
    n = ctx['n']
    E = ctx['E']
    ep = ctx['edge_endpoints']
    edge_index = {}
    for e in range(E):
        u = int(ep[e, 0]); w = int(ep[e, 1])
        edge_index[(u, w) if u <= w else (w, u)] = e

    elems = d4.d4_elements()
    gperm = np.zeros((8, E), dtype=np.int32)
    for k, g in enumerate(elems):
        for e in range(E):
            u = int(ep[e, 0]); w = int(ep[e, 1])
            u2 = g.vertex(u, n); w2 = g.vertex(w, n)
            key = (u2, w2) if u2 <= w2 else (w2, u2)
            gperm[k, e] = edge_index[key]

    gperm_inv = np.zeros((8, E), dtype=np.int32)
    ar = np.arange(E, dtype=np.int32)
    for k in range(8):
        gperm_inv[k, gperm[k]] = ar

    return {'gperm': gperm, 'gperm_inv': gperm_inv,
            'names': [g.name for g in elems]}


def _active_set(state, E):
    return frozenset(int(e) for e in range(E) if state.fixed[e] == ACTIVE)


def _key(S, E):
    """Bitstring do conjunto de arestas ativas (ordem de índice)."""
    b = bytearray(E)
    for e in S:
        b[e] = 1
    return bytes(b)


def _stab_size_and_canonical(S, perms, E):
    """Retorna (|Stab_D4(S)|, is_canonical)."""
    gperm = perms['gperm']
    base = _key(S, E)
    stab = 0
    canonical = True
    for k in range(8):
        gk = gperm[k]
        gS = frozenset(int(gk[e]) for e in S)
        kk = _key(gS, E)
        if kk == base:
            stab += 1
        elif kk < base:
            canonical = False
    return stab, canonical


# Classificação de uma simetria g no nó atual (poda lex-leader sound)
_PRUNE = 0   # provado g·b < b  -> nó não pode ser canônico -> podar
_DEAD = 1    # provado b < g·b   -> g nunca mais ameaça -> descartar
_ALIVE = 2   # empate ou lacuna  -> manter para descendentes


def _classify(fixed, src, E):
    """Classifica uma simetria (via seu array inverso `src`) no estado atual.

    src[p] = índice da aresta g^{-1}(p), de modo que (g·b)[p] = b[src[p]].
    Percorre posições em ordem de índice; para na 1ª diferença decidida ou
    numa lacuna (posição indecidida em qualquer dos lados).
    """
    for p in range(E):
        xp = fixed[p]
        yp = fixed[src[p]]
        if xp == FREE or yp == FREE:
            return _ALIVE            # lacuna: não dá pra provar nada ainda
        xv = 1 if xp == ACTIVE else 0
        yv = 1 if yp == ACTIVE else 0
        if xv == yv:
            continue
        return _PRUNE if yv < xv else _DEAD
    return _ALIVE                    # tudo decidido e igual (estabilizador)


# ---------------------------------------------------------------------------
# Enumeração
# ---------------------------------------------------------------------------

def enumerate_tours(n, use_symmetry=False, seed=0, collect_sets=False):
    """Enumera COMPLETAMENTE os tours fechados do cavalo em n×n.

    use_symmetry=False: conta cada tour uma vez (contagem bruta).
    use_symmetry=True : enumera um representante por D4-órbita e reconstrói
                        a contagem total.

    Retorna dict:
      count           — contagem total de tours (bruta ou reconstruída)
      nodes           — nós de busca explorados
      seconds         — tempo de parede
      raw_leaves      — nº de folhas-tour visitadas (sem simetria == count;
                        com simetria == nº de órbitas/representantes)
      orbit_sizes     — (só com simetria) lista dos tamanhos de órbita
      sets            — (se collect_sets) lista de frozensets de arestas ativas
    """
    ctx = build_graph(n)
    E = ctx['E']
    rng = np.random.default_rng(seed)
    state = State(ctx['V'], ctx['E'])

    perms = build_edge_perms(ctx) if use_symmetry else None

    stats = {'nodes': 0, 'raw_leaves': 0, 'count': 0}
    orbit_sizes = []
    sets = [] if collect_sets else None

    def record_leaf():
        stats['raw_leaves'] += 1
        if collect_sets:
            sets.append(_active_set(state, E))
        if not use_symmetry:
            stats['count'] += 1

    def record_leaf_sym(S):
        stab, canonical = _stab_size_and_canonical(S, perms, E)
        if not canonical:
            return
        stats['raw_leaves'] += 1
        osize = 8 // stab
        orbit_sizes.append(osize)
        stats['count'] += osize
        if collect_sets:
            sets.append(S)

    t0 = time.time()

    status = _propagate_initial(state, ctx)
    if status in (CONTRADICTION, SUBTOUR):
        out = {'count': 0, 'nodes': 0, 'seconds': time.time() - t0,
               'raw_leaves': 0, 'orbit_sizes': [], 'sets': sets}
        return out
    if status == COMPLETE_TOUR:
        if _is_complete_tour(state, ctx):
            if use_symmetry:
                record_leaf_sym(_active_set(state, E))
            else:
                record_leaf()
        return {'count': stats['count'], 'nodes': 1,
                'seconds': time.time() - t0, 'raw_leaves': stats['raw_leaves'],
                'orbit_sizes': orbit_sizes, 'sets': sets}

    # --- DFS sem simetria ---
    def dfs_plain():
        stats['nodes'] += 1
        e, first = _choose_next_edge(state, ctx, rng)
        if e == -1:
            if _is_complete_tour(state, ctx):
                record_leaf()
            return
        for val in (first, 1 - first):
            snap = state.snapshot()
            st = fix_and_propagate(state, ctx, e, val)
            if st == OK:
                dfs_plain()
            elif st == COMPLETE_TOUR:
                if _is_complete_tour(state, ctx):
                    record_leaf()
            state.restore(snap)

    # --- DFS com poda lex-leader D4 ---
    gperm_inv = perms['gperm_inv'] if use_symmetry else None

    def dfs_sym(alive):
        stats['nodes'] += 1
        fixed = state.fixed
        next_alive = []
        for k in alive:
            res = _classify(fixed, gperm_inv[k], E)
            if res == _PRUNE:
                return                 # nó não-canônico: poda a subárvore
            if res == _ALIVE:
                next_alive.append(k)
            # _DEAD: descartado
        e, first = _choose_next_edge(state, ctx, rng)
        if e == -1:
            if _is_complete_tour(state, ctx):
                record_leaf_sym(_active_set(state, E))
            return
        for val in (first, 1 - first):
            snap = state.snapshot()
            st = fix_and_propagate(state, ctx, e, val)
            if st == OK:
                dfs_sym(next_alive)
            elif st == COMPLETE_TOUR:
                if _is_complete_tour(state, ctx):
                    record_leaf_sym(_active_set(state, E))
            state.restore(snap)

    if use_symmetry:
        dfs_sym(list(range(1, 8)))     # identidade (k=0) nunca poda
    else:
        dfs_plain()

    return {'count': stats['count'], 'nodes': stats['nodes'],
            'seconds': time.time() - t0, 'raw_leaves': stats['raw_leaves'],
            'orbit_sizes': orbit_sizes, 'sets': sets}


if __name__ == '__main__':
    for sym in (False, True):
        r = enumerate_tours(6, use_symmetry=sym, seed=0)
        tag = 'D4 ' if sym else 'raw'
        print(f"[{tag}] count={r['count']} nodes={r['nodes']} "
              f"leaves={r['raw_leaves']} t={r['seconds']:.2f}s")
