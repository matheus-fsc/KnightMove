"""Investiga por que o patch nunca dispara: quais filtros falham?"""

# --- bootstrap: torna core/ importavel a partir de qualquer subpasta ---
import sys as _sys, pathlib as _pathlib
for _anc in _pathlib.Path(__file__).resolve().parents:
    if (_anc / "core").is_dir() and (_anc / "pyproject.toml").exists():
        _sys.path.insert(0, str(_anc / "core"))
        break
# --- fim do bootstrap ---
import numpy as np
from collections import Counter
import knight_tours_patch as kt
from knight_tours_patch import (
    build_graph, State, _propagate_initial, _process_queue,
    fix_and_propagate, _choose_next_edge, _is_complete_tour, _extract_tour,
    uf_find, FREE, ACTIVE, INACTIVE, OK, CONTRADICTION, SUBTOUR, COMPLETE_TOUR,
)

# Versão instrumentada: ao invés de aplicar patch, conta motivos de falha
FAILURE_REASONS = Counter()
NO_LUT_ENTRY = Counter()  # quantas vezes LUT[e1] estava vazia


def _process_queue_instrumented(state, ctx, queue):
    n2 = ctx['n2']
    total_inc = ctx['total_incident']
    adj_edges = ctx['adj_edges']
    edge_endpoints = ctx['edge_endpoints']
    patch_LUT = ctx['patch_LUT']

    while queue:
        e, val = queue.pop()
        cur = state.fixed[e]
        if cur != FREE:
            if (cur == ACTIVE and val == 1) or (cur == INACTIVE and val == 0):
                continue
            return CONTRADICTION

        u = int(edge_endpoints[e, 0])
        v = int(edge_endpoints[e, 1])

        if val == 1:
            state.fixed[e] = ACTIVE
            state.n_free -= 1
            state.degree[u] += 1
            state.degree[v] += 1
            if state.degree[u] > 2 or state.degree[v] > 2:
                return CONTRADICTION

            ru = uf_find(state, u)
            rv = uf_find(state, v)
            new_deg2 = (1 if state.degree[u] == 2 else 0) + \
                       (1 if state.degree[v] == 2 else 0)

            if ru == rv:
                state.uf_deg2[ru] += new_deg2
                if state.uf_size[ru] == n2:
                    return COMPLETE_TOUR

                # Contar motivo da falha de cada candidato no LUT
                if len(patch_LUT[e]) == 0:
                    NO_LUT_ENTRY['empty'] += 1

                worst_reason = None
                for (e2, e3, e4) in patch_LUT[e]:
                    if state.fixed[e2] != ACTIVE:
                        FAILURE_REASONS['e2_not_active'] += 1
                        continue
                    if state.fixed[e3] != FREE:
                        FAILURE_REASONS['e3_not_free'] += 1
                        continue
                    if state.fixed[e4] != FREE:
                        FAILURE_REASONS['e4_not_free'] += 1
                        continue
                    u2 = int(edge_endpoints[e2, 0])
                    if uf_find(state, u2) == ru:
                        FAILURE_REASONS['e2_same_component'] += 1
                        continue
                    FAILURE_REASONS['would_apply'] += 1  # nunca chega aqui no original
                    break

                return SUBTOUR

            if state.uf_size[ru] < state.uf_size[rv]:
                ru, rv = rv, ru
            state.uf_parent[rv] = ru
            state.uf_size[ru] += state.uf_size[rv]
            state.uf_deg2[ru] += state.uf_deg2[rv] + new_deg2
        else:
            state.fixed[e] = INACTIVE
            state.n_free -= 1
            state.n_inactive[u] += 1
            state.n_inactive[v] += 1
            if total_inc[u] - state.n_inactive[u] < 2:
                return CONTRADICTION
            if total_inc[v] - state.n_inactive[v] < 2:
                return CONTRADICTION

        for w in (u, v):
            dw = int(state.degree[w])
            avail = int(total_inc[w] - state.n_inactive[w])
            if dw == 2:
                for ei in adj_edges[w]:
                    if state.fixed[ei] == FREE:
                        queue.append((int(ei), 0))
            elif avail == 2 and dw < 2:
                for ei in adj_edges[w]:
                    if state.fixed[ei] == FREE:
                        queue.append((int(ei), 1))

    return OK


# Monkey-patch
kt._process_queue = _process_queue_instrumented


def run(n, K, seed=0):
    FAILURE_REASONS.clear()
    NO_LUT_ENTRY.clear()
    tours = kt.knight_tours_patch(n, K, seed=seed)
    return tours, dict(FAILURE_REASONS), dict(NO_LUT_ENTRY)


if __name__ == '__main__':
    print("=" * 60)
    print("Instrumentação: por que o patch nunca dispara?")
    print("=" * 60)

    for n, K in [(6, 10**6), (10, 200)]:
        tours, fail, empty = run(n, K, seed=0)
        total_fail = sum(fail.values())
        total_subtour = total_fail + empty.get('empty', 0)  # aproximado
        print(f"\nn={n}, K={K}: {len(tours)} tours")
        print(f"  Total checagens de candidato: {total_fail}")
        print(f"  Sub-tours com LUT[e1] vazia : {empty.get('empty', 0)}")
        print(f"  Motivos da falha por candidato:")
        for k, v in sorted(fail.items(), key=lambda kv: -kv[1]):
            pct = 100 * v / total_fail if total_fail else 0
            print(f"    {k:>22s} : {v:>7d}  ({pct:5.1f}%)")
