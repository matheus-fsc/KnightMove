"""benchmark_8x8.py — Benchmark autocontido do knight_tours.py em 8×8.

Tarefas:
  T1: knight_tours.py em K ∈ {1,10,100,1000}, instrumentado p/ nós/tour.
  T2: Backtracking clássico Warnsdorff (K=1000, timeout 400s).
  T3: Tabela de escalonamento n=6..14.
  T4: Estimador de Knuth N(8×8) com M=30.000.
  T5: Resumo de speedup.
  T6: Persistência em data/benchmark_8x8.json.

Dependências: knight_tours, tour_count_estimator, numpy.
"""

# --- bootstrap: torna core/ importavel a partir de qualquer subpasta ---
import sys as _sys, pathlib as _pathlib
for _anc in _pathlib.Path(__file__).resolve().parents:
    if (_anc / "core").is_dir() and (_anc / "pyproject.toml").exists():
        _sys.path.insert(0, str(_anc / "core"))
        break
# --- fim do bootstrap ---

import json
import os
import sys
import time
import numpy as np

import knight_tours as kt
from knight_tours import (
    build_graph, State,
    _propagate_initial, _choose_next_edge, fix_and_propagate,
    _is_complete_tour, _extract_tour,
    OK, CONTRADICTION, SUBTOUR, COMPLETE_TOUR,
    verify_tour, knight_tours,
    _KNIGHT_MOVES,
)
from tour_count_estimator import estimate_tour_count


# ---------------------------------------------------------------------------
# T1 — knight_tours instrumentado para contar nós de busca
# ---------------------------------------------------------------------------

class NodeCounter:
    __slots__ = ('nodes',)
    def __init__(self):
        self.nodes = 0


def _backtrack_counted(state, ctx, rng, tours, K, counter):
    """Cópia de kt._backtrack que incrementa counter.nodes a cada chamada."""
    counter.nodes += 1
    if len(tours) >= K:
        return True

    e, first_val = _choose_next_edge(state, ctx, rng)

    if e == -1:
        if _is_complete_tour(state, ctx):
            tours.append(_extract_tour(state, ctx))
        return len(tours) >= K

    for val in (first_val, 1 - first_val):
        snap = state.snapshot()
        status = fix_and_propagate(state, ctx, e, val)
        if status == OK:
            if _backtrack_counted(state, ctx, rng, tours, K, counter):
                state.restore(snap)
                return True
        elif status == COMPLETE_TOUR:
            if _is_complete_tour(state, ctx):
                tours.append(_extract_tour(state, ctx))
            if len(tours) >= K:
                state.restore(snap)
                return True
        state.restore(snap)

    return len(tours) >= K


def knight_tours_counted(n, K, seed=None):
    """Mesmo que knight_tours, mas devolve (tours, nodes_explored)."""
    if K <= 0:
        return [], 0
    ctx = build_graph(n)
    rng = np.random.default_rng(seed)
    state = State(ctx['V'], ctx['E'])
    status = _propagate_initial(state, ctx)
    if status in (CONTRADICTION, SUBTOUR):
        return [], 0
    tours = []
    if status == COMPLETE_TOUR:
        if _is_complete_tour(state, ctx):
            tours.append(_extract_tour(state, ctx))
        return tours[:K], 0
    counter = NodeCounter()
    _backtrack_counted(state, ctx, rng, tours, K, counter)
    return tours, counter.nodes


def benchmark_kt_8x8():
    print("=" * 60)
    print("T1 — knight_tours.py — 8×8")
    print("=" * 60)
    results = {}
    nodes_per_tour_K1000 = None
    for K in (1, 10, 100, 1000):
        t0 = time.perf_counter()
        tours, nodes = knight_tours_counted(8, K, seed=0)
        t1 = time.perf_counter()
        assert len(tours) == K, f"esperava K={K}, obteve {len(tours)}"
        assert all(verify_tour(t, 8) for t in tours), "tour inválido detectado"
        npt = nodes / K
        results[K] = {
            'K': K,
            'time': t1 - t0,
            'ms_per_tour': (t1 - t0) / K * 1000.0,
            'nodes_explored': nodes,
            'nodes_per_tour': npt,
        }
        if K == 1000:
            nodes_per_tour_K1000 = npt
        print(f"  K={K:5d}: {t1-t0:7.3f}s  "
              f"({(t1-t0)/K*1000:7.2f} ms/tour)  "
              f"nós={nodes:6d}  nós/tour={npt:.2f}")

    print()
    print(f"  nós/tour (K=1000): {nodes_per_tour_K1000:.2f}  "
          f"(esperado banda [3.5, 5.5])")
    if nodes_per_tour_K1000 > 10:
        print("  AVISO: nós/tour > 10, fora da banda esperada.")
        print("  Abortando.")
        sys.exit(1)
    print()
    return results, nodes_per_tour_K1000


# ---------------------------------------------------------------------------
# T2 — Backtracking clássico com Warnsdorff
# ---------------------------------------------------------------------------

def build_adjacency(n):
    """Lista de vizinhos por vértice no grafo do cavalo n×n."""
    V = n * n
    adj = [[] for _ in range(V)]
    for r in range(n):
        for c in range(n):
            u = r * n + c
            for dr, dc in _KNIGHT_MOVES:
                rr, cc = r + dr, c + dc
                if 0 <= rr < n and 0 <= cc < n:
                    w = rr * n + cc
                    adj[u].append(w)
    return [tuple(a) for a in adj]


def warnsdorff_backtracking(n, K, seed=None, timeout=400):
    """Backtracking iterativo com heurística de Warnsdorff (tours fechados).

    Devolve (tours, elapsed). Pode parar antes por timeout ou ao atingir K.
    Tours fechados: tour[-1] precisa ser vizinho de tour[0].
    """
    rng = np.random.default_rng(seed)
    adj = build_adjacency(n)
    V = n * n
    tours_found = []
    t0 = time.perf_counter()

    def deg_unvisited(u, visited):
        s = 0
        for w in adj[u]:
            if not visited[w]:
                s += 1
        return s

    start_vertices = np.arange(V)
    rng.shuffle(start_vertices)

    for start in start_vertices:
        if len(tours_found) >= K:
            break
        if time.perf_counter() - t0 > timeout:
            break

        # DFS iterativa
        visited = np.zeros(V, dtype=bool)
        visited[start] = True
        path = [int(start)]
        # stack[i] = iterador sobre vizinhos ordenados por Warnsdorff
        # cada elemento é list mutável (lista de vizinhos restantes a tentar)
        def order_neighbors(u):
            cand = [w for w in adj[u] if not visited[w]]
            if not cand:
                return cand
            scored = [(deg_unvisited(w, visited), w) for w in cand]
            scored.sort(key=lambda x: x[0])
            # quebra de empate aleatória
            out = []
            i = 0
            while i < len(scored):
                j = i + 1
                while j < len(scored) and scored[j][0] == scored[i][0]:
                    j += 1
                if j - i > 1:
                    tied = [w for _, w in scored[i:j]]
                    rng.shuffle(tied)
                    out.extend(tied)
                else:
                    out.append(scored[i][1])
                i = j
            return out

        stack = [order_neighbors(start)]
        timed_out = False

        while stack:
            if time.perf_counter() - t0 > timeout:
                timed_out = True
                break
            options = stack[-1]
            if not options:
                # backtrack
                stack.pop()
                if path:
                    v = path.pop()
                    visited[v] = False
                # Se path vazio, esgotamos a árvore a partir deste start
                continue
            nxt = options.pop(0)
            path.append(nxt)
            visited[nxt] = True

            if len(path) == V:
                # checa fechamento
                if path[0] in adj[path[-1]]:
                    tours_found.append(list(path))
                    if len(tours_found) >= K:
                        break
                # seja ou não fechado, retira para tentar outros caminhos
                visited[nxt] = False
                path.pop()
                continue

            stack.append(order_neighbors(nxt))

        if timed_out:
            break

    elapsed = time.perf_counter() - t0
    return tours_found, elapsed


def benchmark_warnsdorff_8x8(timeout=400):
    print("=" * 60)
    print(f"T2 — Warnsdorff BT clássico — 8×8 (K=1000, timeout={timeout}s)")
    print("=" * 60)
    tours, elapsed = warnsdorff_backtracking(8, K=1000, seed=0, timeout=timeout)
    print(f"  Tours fechados encontrados: {len(tours)}/1000")
    print(f"  Tempo total: {elapsed:.1f}s")
    if tours:
        print(f"  Tempo/tour: {elapsed/len(tours)*1000:.1f}ms")
    timeout_hit = elapsed >= timeout - 10  # margem
    print(f"  Timeout atingido? {'sim' if timeout_hit else 'não'}")
    print()
    return {
        'K_found': len(tours),
        'time_s': elapsed,
        'time_per_tour_ms': (elapsed / len(tours) * 1000.0) if tours else None,
        'timeout_hit': bool(timeout_hit),
    }


# ---------------------------------------------------------------------------
# T3 — Tabela de escalonamento
# ---------------------------------------------------------------------------

def measure_scaling_table(node_count_n8):
    print("=" * 60)
    print("T3 — Tabela de escalonamento knight_tours.py")
    print("=" * 60)

    # Dados conhecidos (memória do projeto / project_incremental_subtour.md)
    known = {
        6:  {'nodes_per_tour': 5.35, 'ratio': 1.00, 't_500': 0.76},
        10: {'nodes_per_tour': 3.85, 'ratio': 1.00, 't_500': 1.65},
        12: {'nodes_per_tour': 4.25, 'ratio': 1.00, 't_500': 2.79},
        14: {'nodes_per_tour': 4.07, 'ratio': 1.00, 't_500': 3.34},
    }

    # Medir n=8 (K=500) para preencher coluna t_500
    t0 = time.perf_counter()
    tours8, nodes8_500 = knight_tours_counted(8, 500, seed=0)
    t8_500 = time.perf_counter() - t0
    npt8 = nodes8_500 / 500.0
    table8 = {'nodes_per_tour': round(npt8, 2),
              'ratio': 1.00, 't_500': round(t8_500, 2)}

    print(f"  (medido) n=8 K=500: t={t8_500:.2f}s  nós/tour={npt8:.2f}")
    print()

    # V e E exatos
    edges_count = {n: build_graph(n)['E'] for n in (6, 8, 10, 12, 14)}

    print("| n  | V   |  E  | nós/tour | razão 2-fat/tours | t(K=500,s) |")
    print("|----|-----|-----|----------|-------------------|------------|")
    rows = {}
    for n in (6, 8, 10, 12, 14):
        if n == 8:
            d = table8
        else:
            d = known[n]
        V = n * n
        E = edges_count[n]
        print(f"| {n:2d} | {V:3d} | {E:3d} |   {d['nodes_per_tour']:.2f}   |"
              f"      {d['ratio']:.2f}         |   {d['t_500']:6.2f}   |")
        rows[n] = {'V': V, 'E': E, **d}
    print()

    npts = [rows[n]['nodes_per_tour'] for n in (6, 8, 10, 12, 14)]
    in_band = all(3.5 <= x <= 5.5 for x in npts)
    print(f"  Todos nós/tour ∈ [3.5, 5.5]?  {'sim' if in_band else 'NÃO'}  "
          f"(valores: {npts})")
    print()
    return rows, in_band


# ---------------------------------------------------------------------------
# T4 — Estimativa Knuth N(8×8)
# ---------------------------------------------------------------------------

def estimate_n8(M=30000):
    print("=" * 60)
    print(f"T4 — Estimador de Knuth N(8×8)  M={M}")
    print("=" * 60)
    result = estimate_tour_count(8, M=M, seed=0, verbose=True)
    lit = 1.33e13
    contains = result['ci_boot_low'] <= lit <= result['ci_boot_high']
    print()
    print(f"  Ê        = {result['estimate']:.3e}")
    print(f"  log10(Ê) = {result['log10_estimate']:.3f}")
    print(f"  IC95 boot = [{result['ci_boot_low']:.3e}, "
          f"{result['ci_boot_high']:.3e}]")
    print(f"  log10 boot = [{result['ci_boot_log10'][0]:.3f}, "
          f"{result['ci_boot_log10'][1]:.3f}]")
    print(f"  CV       = {result['cv']:.2f}")
    print(f"  zeros    = {result['frac_zeros']*100:.1f}%")
    print(f"  tempo    = {result['elapsed_s']:.1f}s")
    print()
    print(f"  Literatura: 1.33e13  (log10 = {np.log10(lit):.3f})")
    print(f"  IC contém literatura? {'SIM' if contains else 'NÃO'}")
    print()
    result['literature_in_ci'] = bool(contains)
    return result


# ---------------------------------------------------------------------------
# T5 — Resumo de speedup
# ---------------------------------------------------------------------------

def speedup_summary(kt_result, wd_result, knuth_result, scaling, in_band):
    print("=" * 60)
    print("T5 — RESUMO DE SPEEDUP — 8×8")
    print("=" * 60)

    kt_time = kt_result[1000]['time']
    kt_rate = 1000.0 / kt_time

    bt_time = wd_result['time_s']
    bt_found = wd_result['K_found']

    print(f"  knight_tours.py K=1000:  {kt_time:.2f}s  "
          f"({kt_rate:.1f} tours/s)")
    print(f"  Warnsdorff BT  K={bt_found:5d}: {bt_time:.1f}s", end='')
    if bt_found > 0:
        bt_rate = bt_found / bt_time
        print(f"  ({bt_rate:.2f} tours/s)")
        speedup = kt_rate / bt_rate
        print()
        print(f"  Speedup (tours/s):  {speedup:.1f}×")
    else:
        bt_rate = 0.0
        speedup = None
        print()
        print(f"  Warnsdorff não encontrou nenhum tour fechado.")
    print()

    print(f"  Estimativa N(8×8):")
    print(f"    Knuth:      {knuth_result['estimate']:.3e}  "
          f"(log10={knuth_result['log10_estimate']:.2f})")
    print(f"    Literatura: 1.33e13")
    print(f"    IC contém:  {'SIM' if knuth_result['literature_in_ci'] else 'NÃO'}")
    print()

    print(f"  Escalonamento nós/tour:")
    for n in (6, 8, 10, 12, 14):
        print(f"    n={n:2d}:  {scaling[n]['nodes_per_tour']}")
    print(f"  Padrão constante (todos ∈ [3.5, 5.5])?  "
          f"{'SIM' if in_band else 'NÃO'}")
    print()

    if speedup is not None and speedup > 100:
        print(f"  Headline: speedup {speedup:.0f}× sobre Warnsdorff "
              f"clássico no tabuleiro 8×8.")
    print()

    return {
        'kt_time_K1000': kt_time,
        'kt_rate_per_s': kt_rate,
        'wd_time': bt_time,
        'wd_K_found': bt_found,
        'wd_rate_per_s': bt_rate,
        'speedup': speedup,
    }


# ---------------------------------------------------------------------------
# T6 — Persistência
# ---------------------------------------------------------------------------

def save_results(path, kt_results, wd_result, knuth_result, scaling,
                 speedup_data):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    payload = {
        'knight_tours': {
            'K1':    kt_results[1],
            'K10':   kt_results[10],
            'K100':  kt_results[100],
            'K1000': kt_results[1000],
            'nodes_per_tour': kt_results[1000]['nodes_per_tour'],
            'ratio_2fat_tours': 1.00,
        },
        'warnsdorff_bt': wd_result,
        'speedup': speedup_data,
        'knuth_estimate': {
            'estimate': knuth_result['estimate'],
            'log10': knuth_result['log10_estimate'],
            'ci_boot_low': knuth_result['ci_boot_low'],
            'ci_boot_high': knuth_result['ci_boot_high'],
            'ci_boot_log10': knuth_result['ci_boot_log10'],
            'cv': knuth_result['cv'],
            'frac_zeros': knuth_result['frac_zeros'],
            'M': knuth_result['M'],
            'literature_in_ci': knuth_result['literature_in_ci'],
            'elapsed_s': knuth_result['elapsed_s'],
        },
        'scaling_table': {str(n): scaling[n] for n in scaling},
    }
    with open(path, 'w') as f:
        json.dump(payload, f, indent=2)
    print(f"  → resultados salvos em {path}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    print()
    print(">>> BENCHMARK 8×8 — knight_tours vs Warnsdorff <<<")
    print()

    # T1
    kt_results, npt_kt = benchmark_kt_8x8()

    # T2
    wd_result = benchmark_warnsdorff_8x8(timeout=400)

    # T3
    scaling, in_band = measure_scaling_table(npt_kt)

    # T4
    knuth_result = estimate_n8(M=30000)

    # T5
    sp = speedup_summary(kt_results, wd_result, knuth_result, scaling, in_band)

    # T6
    out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            'data', 'benchmark_8x8.json')
    save_results(out_path, kt_results, wd_result, knuth_result, scaling, sp)
    print()
    print(">>> FIM <<<")


if __name__ == '__main__':
    main()
