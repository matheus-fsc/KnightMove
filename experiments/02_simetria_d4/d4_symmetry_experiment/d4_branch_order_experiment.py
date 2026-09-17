"""d4_branch_order_experiment.py — A ordem de ramificação limita a poda D4?

Hipótese falsificável: os NÓS D4 caem só 1.85x (vs 7.92x nas folhas) porque a
ordem de ramificação DINÂMICA (pressão de vértice, _choose_next_edge) decide
arestas fora da ordem de índice, criando "lacunas" que impedem o _classify
(poda lex-leader, que compara bitstrings em ordem de índice) de podar cedo.

Este arquivo reimplementa a enumeração D4 (copiando a lógica de
d4_orbital_branching.py para poder parametrizar a ordem) com DOIS modos de
seleção de aresta de ramificação, mantendo TODO o resto idêntico:

  order="dynamic" : _choose_next_edge do motor (pressão de vértice)
  order="index"   : sempre a aresta LIVRE de menor índice (alinhada ao _classify)

Corretude: qualquer ordem tem de dar 9862 (n=6) / 1245 órbitas. Se não der,
é bug da variante de ordem — reportar, não mascarar.

Nada dos módulos existentes é modificado; tudo é importado.
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

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

import knight_tours as kt
from knight_tours import (
    build_graph, State, fix_and_propagate, _choose_next_edge,
    _propagate_initial, _is_complete_tour, edge_level, edge_priority,
    FREE, ACTIVE, INACTIVE, OK, SUBTOUR, CONTRADICTION, COMPLETE_TOUR,
)
from d4_orbital_branching import (
    build_edge_perms, _active_set, _stab_size_and_canonical, _classify,
    _PRUNE, _ALIVE,
)


class _Timeout(Exception):
    pass


def _make_chooser(order, ctx):
    """Retorna choose(state, rng) -> (edge_idx, first_val) para o modo dado."""
    if order == "dynamic":
        def choose(state, rng):
            return _choose_next_edge(state, ctx, rng)
        return choose
    if order == "index":
        E = ctx['E']
        n = ctx['n']
        ep = ctx['edge_endpoints']
        # valor inicial por aresta (mesma heurística edge_priority do motor),
        # pré-computado — só a SELEÇÃO da aresta muda entre os modos.
        first_of = np.zeros(E, dtype=np.int8)
        for e in range(E):
            a = int(ep[e, 0]); b = int(ep[e, 1])
            _, fv = edge_priority(edge_level(a, b, n))
            first_of[e] = fv

        def choose(state, rng):
            fixed = state.fixed
            for e in range(E):
                if fixed[e] == FREE:
                    return e, int(first_of[e])
            return -1, 0
        return choose
    raise ValueError(f"ordem desconhecida: {order}")


def enumerate_variant(n, use_symmetry, order, seed=0, time_budget=None):
    """Enumeração COMPLETA parametrizada pela ordem de ramificação.

    Retorna dict com: count, orbits (=raw_leaves), nodes, seconds, completed.
    completed=False sinaliza timeout (números são parciais e NÃO devem ser
    interpretados como contagem final).
    """
    ctx = build_graph(n)
    E = ctx['E']
    rng = np.random.default_rng(seed)
    state = State(ctx['V'], ctx['E'])
    choose = _make_chooser(order, ctx)

    perms = build_edge_perms(ctx) if use_symmetry else None
    gperm_inv = perms['gperm_inv'] if use_symmetry else None

    stats = {'nodes': 0, 'raw_leaves': 0, 'count': 0}
    orbit_sizes = []
    deadline = (time.time() + time_budget) if time_budget else None
    t0 = time.time()

    def check_deadline():
        if deadline and (stats['nodes'] & 0x3FFF) == 0 and time.time() > deadline:
            raise _Timeout()

    def record_plain():
        stats['raw_leaves'] += 1
        stats['count'] += 1

    def record_sym(S):
        stab, canonical = _stab_size_and_canonical(S, perms, E)
        if not canonical:
            return
        stats['raw_leaves'] += 1
        osize = 8 // stab
        orbit_sizes.append(osize)
        stats['count'] += osize

    def dfs_plain():
        stats['nodes'] += 1
        check_deadline()
        e, first = choose(state, rng)
        if e == -1:
            if _is_complete_tour(state, ctx):
                record_plain()
            return
        for val in (first, 1 - first):
            snap = state.snapshot()
            st = fix_and_propagate(state, ctx, e, val)
            if st == OK:
                dfs_plain()
            elif st == COMPLETE_TOUR:
                if _is_complete_tour(state, ctx):
                    record_plain()
            state.restore(snap)

    def dfs_sym(alive):
        stats['nodes'] += 1
        check_deadline()
        fixed = state.fixed
        next_alive = []
        for k in alive:
            res = _classify(fixed, gperm_inv[k], E)
            if res == _PRUNE:
                return
            if res == _ALIVE:
                next_alive.append(k)
        e, first = choose(state, rng)
        if e == -1:
            if _is_complete_tour(state, ctx):
                record_sym(_active_set(state, E))
            return
        for val in (first, 1 - first):
            snap = state.snapshot()
            st = fix_and_propagate(state, ctx, e, val)
            if st == OK:
                dfs_sym(next_alive)
            elif st == COMPLETE_TOUR:
                if _is_complete_tour(state, ctx):
                    record_sym(_active_set(state, E))
            state.restore(snap)

    completed = True
    try:
        status = _propagate_initial(state, ctx)
        if status in (CONTRADICTION, SUBTOUR):
            pass
        elif status == COMPLETE_TOUR:
            if _is_complete_tour(state, ctx):
                if use_symmetry:
                    record_sym(_active_set(state, E))
                else:
                    record_plain()
            stats['nodes'] = 1
        else:
            if use_symmetry:
                dfs_sym(list(range(1, 8)))
            else:
                dfs_plain()
    except _Timeout:
        completed = False

    return {'count': stats['count'], 'orbits': stats['raw_leaves'],
            'nodes': stats['nodes'], 'seconds': time.time() - t0,
            'completed': completed}


# ---------------------------------------------------------------------------
# Bancada de medição
# ---------------------------------------------------------------------------

def run_grid(n, seeds, time_budget=None, expect_count=None, expect_orbits=None):
    print(f"\n########## n={n} ##########")
    header = (f"{'modo':>8} {'sym':>4} {'seed':>4} | {'count':>12} "
              f"{'orbits':>8} {'nodes':>10} {'leaves':>8} {'t(s)':>7} {'ok':>4}")
    print(header)
    print("-" * len(header))
    results = {}
    for order in ("dynamic", "index"):
        for sym in (False, True):
            key = (order, sym)
            results[key] = []
            for s in seeds:
                r = enumerate_variant(n, use_symmetry=sym, order=order,
                                      seed=s, time_budget=time_budget)
                results[key].append(r)
                leaves = r['orbits']
                okstr = "OK" if r['completed'] else "T/O"
                # checagens de corretude (só quando completou)
                if r['completed'] and expect_count is not None:
                    if r['count'] != expect_count:
                        okstr = "BUG!"
                    if sym and expect_orbits is not None and r['orbits'] != expect_orbits:
                        okstr = "BUG!"
                print(f"{order:>8} {str(sym):>4} {s:>4} | {r['count']:>12} "
                      f"{leaves:>8} {r['nodes']:>10} {leaves:>8} "
                      f"{r['seconds']:>7.2f} {okstr:>4}")
    return results


def summarize(n, results, expect_count, expect_orbits):
    def avg(key, field):
        rs = [r for r in results[key] if r['completed']]
        return np.mean([r[field] for r in rs]) if rs else float('nan')

    print(f"\n--- Análise n={n} (médias sobre seeds completos) ---")
    for order in ("dynamic", "index"):
        raw_nodes = avg((order, False), 'nodes')
        d4_nodes = avg((order, True), 'nodes')
        raw_t = avg((order, False), 'seconds')
        d4_t = avg((order, True), 'seconds')
        raw_leaves = avg((order, False), 'orbits')
        d4_leaves = avg((order, True), 'orbits')
        xn = raw_nodes / d4_nodes if d4_nodes else float('nan')
        xt = raw_t / d4_t if d4_t else float('nan')
        xl = raw_leaves / d4_leaves if d4_leaves else float('nan')
        print(f"  [{order:>7}] bruto: nodes={raw_nodes:.0f} t={raw_t:.2f}s | "
              f"D4: nodes={d4_nodes:.0f} t={d4_t:.2f}s | "
              f"poda D4 nós={xn:.2f}x folhas={xl:.2f}x tempo={xt:.2f}x")

    # efeito (a): ordem sozinha na busca BRUTA
    rd = avg(("dynamic", False), 'nodes')
    ri = avg(("index", False), 'nodes')
    td = avg(("dynamic", False), 'seconds')
    ti = avg(("index", False), 'seconds')
    print(f"  (a) ordem na busca BRUTA: nodes dynamic={rd:.0f} vs index={ri:.0f} "
          f"(index/dynamic={ri/rd:.2f}x nós, {ti/td:.2f}x tempo)")
    # efeito (b): poda D4 em cada ordem — já impresso acima
    d4d = avg(("dynamic", True), 'nodes')
    d4i = avg(("index", True), 'nodes')
    print(f"  (b) nós D4: dynamic={d4d:.0f} vs index={d4i:.0f} "
          f"(index/dynamic={d4i/d4d:.2f}x)")


if __name__ == '__main__':
    seeds = [0, 1, 2]

    r6 = run_grid(6, seeds, time_budget=None,
                  expect_count=9862, expect_orbits=1245)
    summarize(6, r6, 9862, 1245)

    # n=8: enumeração completa é astronômica (>10^13 tours). Orçamento de tempo
    # honesto; espera-se TIMEOUT. Só 1 seed (index é determinístico de qq modo).
    print("\n(n=8: enumeração completa é inviável — ~1.3e13 tours fechados. "
          "Rodando com orçamento de 120s por config só para registrar o T/O.)")
    r8 = run_grid(8, [0], time_budget=120.0,
                  expect_count=None, expect_orbits=None)
