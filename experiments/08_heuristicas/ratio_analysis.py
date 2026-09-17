"""ratio_analysis.py — Análise empírica da razão tours/2-fatores.

================================================================
DEFINIÇÕES FORMAIS
================================================================

Grafo do cavalo G_n:
  V = {0, ..., n² - 1}     (uma casa por vértice)
  E = pares de casas conectadas por um movimento válido de cavalo

2-fator (de G_n):
  Subconjunto F ⊆ E tal que todo vértice tem grau exatamente 2 em F.
  Equivalente a uma cobertura por ciclos simples disjuntos cobrindo V.

Tour (fechado, não-direcionado):
  2-fator F que é CONEXO (um único ciclo hamiltoniano).

Razão investigada:
  r(n) = N(tours) / N(2-fatores)

Hipótese central:
  Se r(n) segue um padrão (constante, polinomial, ou ditado por
  deficit=3), então N(tours) = N(2-fatores) × r(n) — e como
  N(2-fatores) é calculável via transfer matrix (estrutura linear),
  a decomposição resolve o problema de contagem.

Método (B):
  Backtracking de knight_tours.py SEM o detector de sub-ciclos —
  a propagação R2 é mantida, mas fechamentos prematuros de ciclo
  NÃO podam a árvore. Folhas correspondem a 2-fatores; conta-se
  quantas têm 1 componente (tours).

================================================================
"""

# --- bootstrap: torna core/ importavel a partir de qualquer subpasta ---
import sys as _sys, pathlib as _pathlib
for _anc in _pathlib.Path(__file__).resolve().parents:
    if (_anc / "core").is_dir() and (_anc / "pyproject.toml").exists():
        _sys.path.insert(0, str(_anc / "core"))
        break
# --- fim do bootstrap ---

import sys
import time
import json
import os
import numpy as np

sys.setrecursionlimit(200000)

from knight_tours import (
    build_graph, vertex_level, edge_level, edge_priority,
    FREE, ACTIVE, INACTIVE, OK, CONTRADICTION,
    State, uf_find,
)


# ---------------------------------------------------------------------------
# Backtracking modificado — sem prune por sub-ciclo
# ---------------------------------------------------------------------------

def _process_queue_2fac(state, ctx, queue):
    """Versão de _process_queue do knight_tours.py SEM podar sub-ciclos.

    Quando uma aresta ACTIVE fecha um ciclo dentro do mesmo componente,
    nós apenas registramos (não retornamos SUBTOUR). Isto permite a
    construção de 2-fatores multi-componente (válidos como 2-fatores,
    mas não como tours).

    Retorna OK ou CONTRADICTION.
    """
    total_inc = ctx['total_incident']
    adj_edges = ctx['adj_edges']
    edge_endpoints = ctx['edge_endpoints']

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
                # Fechamento de sub-ciclo: PERMITIDO. Só atualiza contador.
                state.uf_deg2[ru] += new_deg2
            else:
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

        # R2 nos dois endpoints
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


def _fix_and_propagate_2fac(state, ctx, e, val):
    return _process_queue_2fac(state, ctx, [(e, val)])


def _propagate_initial_2fac(state, ctx):
    V = ctx['V']
    total_inc = ctx['total_incident']
    adj_edges = ctx['adj_edges']
    queue = []
    for v in range(V):
        if total_inc[v] - state.n_inactive[v] == 2 and state.degree[v] < 2:
            for ei in adj_edges[v]:
                if state.fixed[ei] == FREE:
                    queue.append((int(ei), 1))
    return _process_queue_2fac(state, ctx, queue)


# ---------------------------------------------------------------------------
# Seleção de variável — pressão de vértice + randomização opcional
# ---------------------------------------------------------------------------

def _choose_next_edge_2fac(state, ctx, rng, use_heuristic=True,
                            randomize_branch=True):
    """Escolhe próxima aresta a fixar.

    use_heuristic: usa f∞(L) para escolher o valor inicial.
    randomize_branch: 50% de chance de inverter a ordem dos valores
      (útil para amostragem mais diversa em K << total).
    """
    V = ctx['V']
    n = ctx['n']
    total_inc = ctx['total_incident']
    adj_edges = ctx['adj_edges']
    edge_endpoints = ctx['edge_endpoints']

    best_score = -1
    candidates = []
    for v in range(V):
        free_v = int(total_inc[v]) - int(state.degree[v]) - int(state.n_inactive[v])
        if free_v == 0:
            continue
        score = int(state.degree[v]) * 100 + \
                int(total_inc[v] - state.n_inactive[v] - 2)
        if score > best_score:
            best_score = score
            candidates = [v]
        elif score == best_score:
            candidates.append(v)

    if not candidates:
        return -1, 0

    if len(candidates) == 1:
        v_star = candidates[0]
    else:
        v_star = int(candidates[int(rng.integers(len(candidates)))])

    if use_heuristic:
        best_pri = -1.0
        best_edge = -1
        best_first = 0
        for ei in adj_edges[v_star]:
            ei_i = int(ei)
            if state.fixed[ei_i] != FREE:
                continue
            a = int(edge_endpoints[ei_i, 0])
            b = int(edge_endpoints[ei_i, 1])
            L = edge_level(a, b, n)
            pri, fv = edge_priority(L)
            if pri > best_pri:
                best_pri = pri
                best_edge = ei_i
                best_first = fv
    else:
        # Sem heurística — primeira aresta livre, valor inicial aleatório
        best_edge = -1
        for ei in adj_edges[v_star]:
            if state.fixed[ei] == FREE:
                best_edge = int(ei)
                break
        best_first = int(rng.integers(2))

    if randomize_branch and rng.random() < 0.5:
        best_first = 1 - best_first

    return best_edge, best_first


# ---------------------------------------------------------------------------
# Contagem de componentes
# ---------------------------------------------------------------------------

def _count_components(state, ctx):
    """Conta componentes conexos via raízes UF.

    Pressuposto: grau ativo = 2 em todo vértice (verdadeiro em folhas
    válidas após R2). Cada componente é um ciclo simples.
    """
    V = ctx['V']
    roots = set()
    for v in range(V):
        roots.add(uf_find(state, v))
    return len(roots)


# ---------------------------------------------------------------------------
# Backtracking principal
# ---------------------------------------------------------------------------

def _backtrack_2fac(state, ctx, rng, stats, K_target, t_start, time_limit,
                     use_heuristic, randomize_branch):
    if stats['n_2fac'] >= K_target:
        return True
    if time_limit is not None and (time.time() - t_start) > time_limit:
        stats['stopped_by_time'] = True
        return True

    e, first_val = _choose_next_edge_2fac(
        state, ctx, rng, use_heuristic=use_heuristic,
        randomize_branch=randomize_branch)

    if e == -1:
        # Folha — todos os vértices têm grau 2 (R2 garante)
        stats['n_2fac'] += 1
        n_comp = _count_components(state, ctx)
        stats['comp_dist'][n_comp] = stats['comp_dist'].get(n_comp, 0) + 1
        if n_comp == 1:
            stats['n_tours'] += 1
        return stats['n_2fac'] >= K_target

    for val in (first_val, 1 - first_val):
        snap = state.snapshot()
        status = _fix_and_propagate_2fac(state, ctx, e, val)
        if status == OK:
            if _backtrack_2fac(state, ctx, rng, stats, K_target, t_start,
                                time_limit, use_heuristic, randomize_branch):
                state.restore(snap)
                return True
        state.restore(snap)

    return stats['n_2fac'] >= K_target


# ---------------------------------------------------------------------------
# Interface pública
# ---------------------------------------------------------------------------

def sample_2factor_restart(n, K_target=10000, seed=42, time_limit=None,
                             use_heuristic=False, max_restarts_factor=200):
    """Amostragem por RESTART: cada amostra é uma descida aleatória independente.

    Cada amostra:
      1. Restaura o estado inicial (após propagação inicial)
      2. Faz escolhas aleatórias até folha (sucesso) OU contradição (descarte)
      3. Conta componentes na folha

    Diferente do DFS, isto NÃO compartilha prefixo entre amostras, evitando
    o viés do "primeiro K leaves" do DFS.

    Para K_target amostras VÁLIDAS, faz até K_target × max_restarts_factor
    restarts (descartando contradições). Reporta também n_restarts e
    success_rate.
    """
    ctx = build_graph(n)
    rng = np.random.default_rng(seed)

    init_state = State(ctx['V'], ctx['E'])
    status = _propagate_initial_2fac(init_state, ctx)
    if status == CONTRADICTION:
        return {
            'n': n, 'n_2fatores': 0, 'n_tours': 0, 'ratio': 0.0,
            'component_distribution': {}, 'K_target': K_target,
            'time_limit': time_limit, 'use_heuristic': use_heuristic,
            'mode': 'restart', 'n_restarts': 0, 'success_rate': 0.0,
            'elapsed_s': 0.0, 'stopped_by_time': False,
        }
    init_snap = init_state.snapshot()

    state = State(ctx['V'], ctx['E'])
    n_2fac = 0
    n_tours = 0
    comp_dist = {}
    n_restarts = 0
    stopped_by_time = False
    t_start = time.time()
    max_restarts = K_target * max_restarts_factor

    while n_2fac < K_target and n_restarts < max_restarts:
        if time_limit is not None and (time.time() - t_start) > time_limit:
            stopped_by_time = True
            break
        state.restore(init_snap)
        ok = _descent(state, ctx, rng, use_heuristic)
        n_restarts += 1
        if ok:
            n_2fac += 1
            n_comp = _count_components(state, ctx)
            comp_dist[n_comp] = comp_dist.get(n_comp, 0) + 1
            if n_comp == 1:
                n_tours += 1

    elapsed = time.time() - t_start
    ratio = n_tours / n_2fac if n_2fac > 0 else 0.0
    return {
        'n': n, 'n_2fatores': n_2fac, 'n_tours': n_tours, 'ratio': ratio,
        'component_distribution': comp_dist, 'K_target': K_target,
        'time_limit': time_limit, 'use_heuristic': use_heuristic,
        'mode': 'restart',
        'n_restarts': n_restarts,
        'success_rate': n_2fac / n_restarts if n_restarts > 0 else 0.0,
        'elapsed_s': elapsed,
        'stopped_by_time': stopped_by_time,
    }


def _descent(state, ctx, rng, use_heuristic):
    """Uma descida aleatória INDEPENDENTE até folha ou contradição.

    Em cada passo: escolhe a aresta via pressão de vértice (MRV) e fixa
    um valor binário aleatório uniforme. Sem retry — primeira contradição
    descarta a descida. Cada chamada bem-sucedida é uma amostra
    estatisticamente independente das outras.
    """
    while True:
        e, _ = _choose_next_edge_2fac(
            state, ctx, rng, use_heuristic=use_heuristic,
            randomize_branch=False)
        if e == -1:
            return True
        val = int(rng.integers(2))
        status = _fix_and_propagate_2fac(state, ctx, e, val)
        if status == CONTRADICTION:
            return False


def sample_2factors_and_tours(n, K_target=10000, seed=42, time_limit=None,
                               use_heuristic=True, randomize_branch=True):
    """Conta 2-fatores e tours via backtracking sem prune por subtour.

    Parâmetros:
      n             : lado do tabuleiro
      K_target      : alvo de 2-fatores; busca para ao atingir K_target
      seed          : semente do RNG
      time_limit    : limite em segundos (None = sem limite)
      use_heuristic : usa f∞(L) para ordenar valor inicial
      randomize_branch: 50% de chance de inverter ordem (diversifica K << total)

    Retorna dict com chaves:
      n, n_2fatores, n_tours, ratio, component_distribution,
      K_target, time_limit, use_heuristic, randomize_branch,
      elapsed_s, stopped_by_time
    """
    ctx = build_graph(n)
    rng = np.random.default_rng(seed)
    state = State(ctx['V'], ctx['E'])

    stats = {
        'n_2fac': 0,
        'n_tours': 0,
        'comp_dist': {},
        'stopped_by_time': False,
    }

    t_start = time.time()
    status = _propagate_initial_2fac(state, ctx)

    if status == CONTRADICTION:
        elapsed = time.time() - t_start
        return {
            'n': n, 'n_2fatores': 0, 'n_tours': 0, 'ratio': 0.0,
            'component_distribution': {}, 'K_target': K_target,
            'time_limit': time_limit, 'use_heuristic': use_heuristic,
            'randomize_branch': randomize_branch,
            'elapsed_s': elapsed, 'stopped_by_time': False,
        }

    # Verifica se já é folha após propagação inicial
    e, _ = _choose_next_edge_2fac(state, ctx, rng,
                                    use_heuristic=use_heuristic,
                                    randomize_branch=False)
    if e == -1:
        stats['n_2fac'] += 1
        n_comp = _count_components(state, ctx)
        stats['comp_dist'][n_comp] = stats['comp_dist'].get(n_comp, 0) + 1
        if n_comp == 1:
            stats['n_tours'] += 1
    else:
        _backtrack_2fac(state, ctx, rng, stats, K_target, t_start,
                         time_limit, use_heuristic, randomize_branch)

    elapsed = time.time() - t_start
    ratio = stats['n_tours'] / stats['n_2fac'] if stats['n_2fac'] > 0 else 0.0

    return {
        'n': n,
        'n_2fatores': stats['n_2fac'],
        'n_tours': stats['n_tours'],
        'ratio': ratio,
        'component_distribution': dict(stats['comp_dist']),
        'K_target': K_target,
        'time_limit': time_limit,
        'use_heuristic': use_heuristic,
        'randomize_branch': randomize_branch,
        'elapsed_s': elapsed,
        'stopped_by_time': stats['stopped_by_time'],
    }


# ---------------------------------------------------------------------------
# Análises (T3..T6)
# ---------------------------------------------------------------------------

def analyze_components(results):
    print("\n=== Distribuição de componentes por n ===\n")
    print(f"{'n':>4} | {'comp=1':>8} | {'comp=2':>8} | {'comp=3':>8} | "
          f"{'comp=4+':>8} | {'razão':>8} | {'n_2fac':>10}")
    print("-" * 75)

    for n in sorted(results.keys()):
        r = results[n]
        dist = r['component_distribution']
        total = r['n_2fatores']
        if total == 0:
            continue
        c1 = dist.get(1, 0) / total
        c2 = dist.get(2, 0) / total
        c3 = dist.get(3, 0) / total
        c4p = sum(v for k, v in dist.items() if k >= 4) / total
        print(f"{n:>4} | {c1:>8.4f} | {c2:>8.4f} | {c3:>8.4f} | "
              f"{c4p:>8.4f} | {r['ratio']:>8.4f} | {total:>10d}")
    print()
    print("Hipótese: se comp=2 domina os não-conexos, falhas são fusão")
    print("de exatamente 2 ciclos → estrutura binária, tratável.")


def analyze_ratio_pattern(results):
    ns = sorted(results.keys())
    ratios = [results[n]['ratio'] for n in ns]
    print("\n=== Análise do padrão de razão(n) ===\n")
    for n, r in zip(ns, ratios):
        if r > 0:
            print(f"  n={n:2d}: razão={r:.4f}  log2(1/razão)={np.log2(1/r):.3f}")
        else:
            print(f"  n={n:2d}: razão=0 (sem dados)")

    valid = [(n, r) for n, r in zip(ns, ratios) if r > 0]
    if len(valid) < 2:
        print("\n  (poucos dados para fit)")
        return {}

    ns_v = np.array([n for n, _ in valid])
    rs_v = np.array([r for _, r in valid])
    log_r = np.log(rs_v)

    slope, intercept = np.polyfit(ns_v, log_r, 1)
    print(f"\nFit exponencial em n: razão(n) ≈ {np.exp(intercept):.4f} "
          f"× exp({slope:.4f}·n)")

    slope2, intercept2 = np.polyfit(ns_v ** 2, log_r, 1)
    print(f"Fit exponencial em n²: razão(n) ≈ {np.exp(intercept2):.4f} "
          f"× exp({slope2:.6f}·n²)")

    print("\nPredições (fit em n):")
    for n_pred in [6, 8, 10, 12, 14]:
        pred = np.exp(intercept + slope * n_pred)
        marker = " ✓" if n_pred in results else ""
        print(f"  n={n_pred:2d}: {pred:.5f}{marker}")

    print("\nClassificação:")
    if abs(slope) < 0.05:
        print("  → CASO 1: razão aproximadamente CONSTANTE em n")
        print("  → decomposição realizável: N(tours) ≈ N(2-fat) × const")
    elif slope < -0.1:
        print("  → CASO 3: razão DECAI exponencialmente em n")
        print("  → conectividade é barreira crescente")
    else:
        print(f"  → CASO 2: razão varia moderadamente (slope={slope:.3f})")

    return {
        'slope_n': float(slope),
        'intercept_n': float(intercept),
        'slope_n2': float(slope2),
        'intercept_n2': float(intercept2),
    }


def analyze_deficit_connection(results):
    print("\n=== Conexão com deficit=3 ===\n")
    print(f"2^(-deficit) = 2^(-3) = {2**(-3):.4f} = 12.5%")
    print()
    print(f"{'n':>4} | {'razão':>8} | {'log2(1/r)':>10} | {'extra':>8}")
    print("-" * 45)
    for n in sorted(results.keys()):
        r = results[n]['ratio']
        if r <= 0:
            continue
        log2_inv = np.log2(1 / r)
        print(f"{n:>4} | {r:>8.4f} | {log2_inv:>10.3f} | "
              f"{log2_inv - 3:>+8.3f}")
    print()
    print("Interpretação:")
    print("  log2(1/r) ≈ 3 e CONSTANTE → deficit=3 explica r(n)")
    print("  log2(1/r) cresce com n   → conectividade > deficit")


def estimate_N_via_ratio(results):
    print("\n=== Estimativa N(tours) via decomposição ===\n")
    n2f_6 = 13422
    tours_6 = 9862
    r6 = tours_6 / n2f_6
    print(f"n=6 (exato):    N(2-fat)={n2f_6:>10d}  ratio={r6:.4f}  "
          f"N(tours)={tours_6:>10d}")
    lit = {8: 1.33e13}
    for n in sorted(results.keys()):
        if n == 6:
            continue
        r = results[n]['ratio']
        if r <= 0:
            continue
        if n in lit:
            n_imp = lit[n] / r
            print(f"n={n:2d}: ratio={r:.4f}  N(tours)_lit={lit[n]:.3e}  "
                  f"→ N(2-fat) implícito = {n_imp:.3e}  "
                  f"log10={np.log10(n_imp):.2f}")
        else:
            print(f"n={n:2d}: ratio={r:.4f}  (sem N(tours) de literatura)")


def generate_plots(results, out_path='data/plots/ratio_analysis.png'):
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib indisponível — pulando plots")
        return

    ns = sorted([n for n, r in results.items() if r['ratio'] > 0])
    ratios = [results[n]['ratio'] for n in ns]
    if not ns:
        print("Sem dados para plotar")
        return

    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    axes[0].plot(ns, ratios, 'o-', color='steelblue', linewidth=2,
                  markersize=10)
    axes[0].axhline(y=0.125, color='red', linestyle='--', alpha=0.6,
                     label='1/8 = 12.5% (deficit=3)')
    axes[0].axhline(y=0.735, color='green', linestyle='--', alpha=0.6,
                     label='73.5% (n=6 exato)')
    axes[0].set_xlabel('n')
    axes[0].set_ylabel('razão = tours / 2-fatores')
    axes[0].set_title('razão(n)')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    log2_inv = [np.log2(1 / r) for r in ratios]
    axes[1].plot(ns, log2_inv, 'o-', color='darkorange', linewidth=2,
                  markersize=10)
    axes[1].axhline(y=3, color='black', linestyle='--', alpha=0.6,
                     label='deficit = 3')
    axes[1].set_xlabel('n')
    axes[1].set_ylabel('log2(1 / razão)')
    axes[1].set_title('dimensões "proibidas" pela conectividade')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)

    largest = ns[-1]
    dist = results[largest]['component_distribution']
    total = results[largest]['n_2fatores']
    if dist and total > 0:
        ks = sorted(dist.keys())
        vals = [dist[k] / total for k in ks]
        axes[2].bar([str(k) for k in ks], vals, color='steelblue')
        axes[2].set_xlabel('número de componentes')
        axes[2].set_ylabel('fração de 2-fatores')
        axes[2].set_title(f'distribuição de componentes (n={largest})')
        axes[2].grid(True, axis='y', alpha=0.3)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    print(f"Plot salvo: {out_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def measure_all_ratios(plan):
    """plan: lista de (n, K_target, time_limit)."""
    results = {}
    for n, K, tl in plan:
        print(f"\n=== n={n}  K_target={K}  time_limit={tl}s ===")
        t0 = time.time()
        r = sample_2factors_and_tours(n, K_target=K, seed=42, time_limit=tl)
        dt = time.time() - t0
        results[n] = r
        print(f"  2-fatores:      {r['n_2fatores']}")
        print(f"  tours:          {r['n_tours']}")
        if r['n_2fatores'] > 0:
            print(f"  razão:          {r['ratio']:.4f}  "
                  f"({r['ratio']*100:.2f}%)")
        print(f"  dist. comp.:    {r['component_distribution']}")
        print(f"  tempo:          {dt:.2f}s  "
              f"(stopped_by_time={r['stopped_by_time']})")
    return results


def main():
    # T1 — validação n=6 (exaustiva)
    print("=" * 60)
    print("T1 — VALIDAÇÃO n=6 (enumeração exaustiva esperada)")
    print("=" * 60)

    r6 = sample_2factors_and_tours(6, K_target=10**9, seed=42,
                                     time_limit=120)
    print(f"  n_2fatores = {r6['n_2fatores']}  (esperado: 13422)")
    print(f"  n_tours    = {r6['n_tours']}  (esperado: 9862)")
    print(f"  ratio      = {r6['ratio']:.4f}  (esperado: ~0.735)")
    print(f"  comp dist  = {r6['component_distribution']}")
    print(f"  elapsed    = {r6['elapsed_s']:.2f}s")

    ratio_err = abs(r6['ratio'] - 9862 / 13422)
    if ratio_err > 0.01:
        print(f"\n!! razão fora de ±1% (erro={ratio_err:.4f}) — BUG provável")
        print("   Abortando T2..T7.")
        return
    if r6['n_2fatores'] != 13422:
        print(f"\n!! n_2fatores={r6['n_2fatores']} ≠ 13422 — enumeração "
              "incompleta ou bug")
        print("   Continuando, mas com cautela.")

    print("\n→ T1 OK, prosseguindo para T2.")

    # T2 — n=8, 10, 12 (amostragem)
    print("\n" + "=" * 60)
    print("T2 — AMOSTRAGEM n ∈ {8, 10, 12}")
    print("=" * 60)

    plan = [
        (8,  10000, 1800),
        (10,  5000, 1800),
        (12,  5000, 1800),
    ]
    results = {6: r6}
    results.update(measure_all_ratios(plan))

    # T3..T7
    analyze_components(results)
    fit = analyze_ratio_pattern(results)
    analyze_deficit_connection(results)
    estimate_N_via_ratio(results)
    generate_plots(results)

    # Salvar JSON
    out = {
        'results': {str(k): v for k, v in results.items()},
        'fit': fit,
        'ground_truth_n6': {'n_2fatores': 13422, 'n_tours': 9862,
                             'ratio': 9862 / 13422},
    }
    os.makedirs('data', exist_ok=True)
    with open('data/ratio_analysis.json', 'w') as f:
        json.dump(out, f, indent=2)
    print("\nResultado salvo: data/ratio_analysis.json")


if __name__ == '__main__':
    main()
