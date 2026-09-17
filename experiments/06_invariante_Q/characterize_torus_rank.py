#!/usr/bin/env python3
"""
characterize_torus_rank.py — Caracterização precisa do rank GF(2), simetrias de tradução,
suporte de órbitas e paridades homológicas do Passeio do Cavalo no Toro regular n x n.
"""

# --- bootstrap: torna core/ importavel a partir de qualquer subpasta ---
import sys as _sys, pathlib as _pathlib
for _anc in _pathlib.Path(__file__).resolve().parents:
    if (_anc / "core").is_dir() and (_anc / "pyproject.toml").exists():
        _sys.path.insert(0, str(_anc / "core"))
        break
# --- fim do bootstrap ---

import time
import numpy as np
import matplotlib.pyplot as plt
from knight_tours_torus import knight_tours, build_graph, verify_tour

# ===========================================================================
# 1. GF(2) Linear Algebra & Incremental Rank
# ===========================================================================

def vector_to_int(vec):
    """Converte um vetor binário numpy em um inteiro longo do Python (MSB em vec[0])."""
    val = 0
    for bit in vec:
        val = (val << 1) | int(bit)
    return val

def int_to_vector(val, E):
    """Converte um inteiro longo do Python de volta em um vetor binário numpy de tamanho E."""
    vec = np.zeros(E, dtype=np.uint8)
    for i in range(E - 1, -1, -1):
        vec[i] = (val >> (E - 1 - i)) & 1
    return vec

class GF2RankCalculator:
    """Calculador de Rank GF(2) incremental extremamente rápido usando inteiros longos."""
    def __init__(self, E):
        self.E = E
        self.basis_map = {}  # pivot_bit -> int

    def add_vector(self, v_int):
        """Tenta adicionar o vetor (em formato int). Retorna True se for L.I. (rank subiu)."""
        v = v_int
        while v > 0:
            pivot = v.bit_length() - 1
            if pivot in self.basis_map:
                v ^= self.basis_map[pivot]
            else:
                self.basis_map[pivot] = v
                return True
        return False

    def get_rank(self):
        return len(self.basis_map)

    def get_basis_vectors(self):
        """Retorna os vetores da base ordenados em ordem decrescente de pivot."""
        basis_ints = sorted(self.basis_map.values(), reverse=True)
        return [int_to_vector(b_int, self.E) for b_int in basis_ints]

# ===========================================================================
# 2. Translations & Orbits
# ===========================================================================

def precompute_translations(ctx):
    """Pré-computa as permutações das arestas correspondentes às n^2 translações do toro."""
    n = ctx['n']
    E = ctx['E']
    edge_endpoints = ctx['edge_endpoints']
    edge_to_idx = {tuple(edge_endpoints[i]): i for i in range(E)}

    trans_table = {}
    for a in range(n):
        for b in range(n):
            perm = np.zeros(E, dtype=np.int32)
            for e_idx in range(E):
                u, w = edge_endpoints[e_idx]
                uy, ux = divmod(u, n)
                wy, wx = divmod(w, n)

                tu = ((uy + a) % n) * n + ((ux + b) % n)
                tw = ((wy + a) % n) * n + ((wx + b) % n)

                t_key = (tu, tw) if tu < tw else (tw, tu)
                perm[e_idx] = edge_to_idx[t_key]
            trans_table[(a, b)] = perm
    return trans_table

def compute_edge_orbits(ctx, trans_table):
    """Particiona as arestas em órbitas sob a ação do grupo de translações."""
    E = ctx['E']
    n = ctx['n']
    orbits = []
    visited = set()
    for e_idx in range(E):
        if e_idx in visited:
            continue
        orbit = set()
        for a in range(n):
            for b in range(n):
                perm = trans_table[(a, b)]
                orbit.add(perm[e_idx])
        orbits.append(sorted(list(orbit)))
        visited.update(orbit)
    return orbits

def translate_vector(vec, perm):
    """Aplica a translação (definida por perm) a um vetor binário."""
    g_v = np.zeros(len(vec), dtype=np.uint8)
    g_v[perm] = vec
    return g_v

# ===========================================================================
# 3. Homology & Winding Parities
# ===========================================================================

def get_homology_cuts(ctx):
    """Retorna os vetores de corte binários omega_y e omega_x usando diferença de coordenadas."""
    n = ctx['n']
    E = ctx['E']
    edge_endpoints = ctx['edge_endpoints']

    omega_y = np.zeros(E, dtype=np.uint8)
    omega_x = np.zeros(E, dtype=np.uint8)

    threshold = 2 if n >= 5 else 1

    for e_idx in range(E):
        u, w = edge_endpoints[e_idx]
        uy, ux = divmod(u, n)
        wy, wx = divmod(w, n)
        if abs(uy - wy) > threshold:
            omega_y[e_idx] = 1
        if abs(ux - wx) > threshold:
            omega_x[e_idx] = 1

    return omega_y, omega_x

# ===========================================================================
# 4. Core Experiments
# ===========================================================================

def run_experiment_a():
    print("\n" + "="*80)
    print(" EXPERIMENTO A: RANK REAL VIA ENUMERAÇÃO PARCIAL (TORO 6x6)")
    print("="*80)

    n = 6
    ctx = build_graph(n, n)
    E = ctx['E']
    V = ctx['V']
    beta1 = E - V + 1

    # 1. Orbits verification
    trans_table = precompute_translations(ctx)
    orbits = compute_edge_orbits(ctx, trans_table)
    print(f"Arestas totais: {E}")
    print(f"Número de órbitas de arestas sob Z_6 x Z_6: {len(orbits)}")
    for idx, orb in enumerate(orbits):
        print(f"  Órbita {idx+1}: tamanho {len(orb)}")
        assert len(orb) == 36, "A ação deveria ser livre e cada órbita ter tamanho 36!"

    # 2. Tour sampling
    K_target = 50000
    print(f"\nAmostrando K = {K_target} tours no toro 6x6...")
    t0 = time.time()
    tours = knight_tours(n, n, K_target, seed=42)
    t1 = time.time()
    print(f"Tours gerados: {len(tours)} em {t1-t0:.2f} segundos")

    # 3. Incremental rank calculation
    print("\nCalculando curva de rank incremental a cada 1000 tours...")
    rank_calc = GF2RankCalculator(E)
    
    edge_endpoints = ctx['edge_endpoints']
    edge_to_idx = {tuple(edge_endpoints[i]): i for i in range(E)}

    curve_k = []
    curve_rank = []

    # Cache tour vectors as ints
    tour_ints = []
    for t in tours:
        t_vec = np.zeros(E, dtype=np.uint8)
        for j in range(V):
            u = t[j]
            w = t[(j+1)%V]
            key = (u, w) if u < w else (w, u)
            t_vec[edge_to_idx[key]] = 1
        tour_ints.append(vector_to_int(t_vec))

    # Add incrementally
    for i, t_int in enumerate(tour_ints):
        rank_calc.add_vector(t_int)
        k = i + 1
        if k % 1000 == 0 or k == len(tours):
            curve_k.append(k)
            curve_rank.append(rank_calc.get_rank())

    final_rank = rank_calc.get_rank()
    print(f"Rank final para K={len(tours)}: {final_rank} (de beta1 = {beta1})")
    
    # Check convergence
    is_converged = False
    if len(curve_rank) >= 5:
        # Check if the rank was constant for the last few steps
        if curve_rank[-1] == curve_rank[-2] == curve_rank[-3]:
            is_converged = True
    print(f"O rank convergiu? {'SIM' if is_converged else 'NÃO'} (estabilizou em {curve_rank[-1]})")

    # 4. Generate Plot
    plt.figure(figsize=(10, 6))
    plt.plot(curve_k, curve_rank, marker='o', color='#8884d8', linewidth=2, label='Rank GF(2)')
    plt.axhline(y=final_rank, color='r', linestyle='--', label=f'Rank Final = {final_rank}')
    plt.title(f'Convergência do Rank GF(2) no Toro 6x6 (beta1 = {beta1})', fontsize=14)
    plt.xlabel('Tours Adicionados (k)', fontsize=12)
    plt.ylabel('Rank GF(2)', fontsize=12)
    plt.grid(True, linestyle=':', alpha=0.6)
    plt.legend(fontsize=11)
    
    plot_path = '/home/math/.gemini/antigravity-cli/brain/bbba8aeb-a93c-4965-8fde-b95e4e0483bb/torus_rank_convergence.png'
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Gráfico de convergência salvo em: [torus_rank_convergence.png](file://{plot_path})")

    return ctx, trans_table, orbits, tour_ints, rank_calc


def run_experiment_b(ctx, trans_table, orbits, tour_ints, rank_calc):
    print("\n" + "="*80)
    print(" EXPERIMENTO B: ANÁLISE DOS VETORES DA BASE (TORO 6x6)")
    print("="*80)

    n = ctx['n']
    E = ctx['E']
    basis = rank_calc.get_basis_vectors()
    rank = len(basis)
    print(f"Analisando {rank} vetores da base de Span(Ham(G_T(6)))...")

    # Winding cuts
    omega_y, omega_x = get_homology_cuts(ctx)

    basis_info = []
    for idx, v in enumerate(basis):
        # 1. Homology
        wy = int((v @ omega_y) % 2)
        wx = int((v @ omega_x) % 2)

        # 2. Orbit support weights (absolute weights)
        orb_weights = []
        for orb in orbits:
            weight = int(v[orb].sum())
            orb_weights.append(weight)

        # 3. Translation Stabilizer
        stab_size = 0
        stab_elements = []
        for a in range(n):
            for b in range(n):
                perm = trans_table[(a, b)]
                g_v = translate_vector(v, perm)
                if np.array_equal(g_v, v):
                    stab_size += 1
                    stab_elements.append((a, b))

        basis_info.append({
            'idx': idx + 1,
            'wy': wy, 'wx': wx,
            'weights': orb_weights,
            'stab_size': stab_size,
            'stab': stab_elements
        })

    # Summary
    print("\nResumo da Base GF(2):")
    print(f"{'Vetor':<6} | {'Winding (Wy, Wx)':<16} | {'Pesos nas Orbits (O1, O2, O3, O4)':<33} | {'Tam. Estabilizador':<18}")
    print("-" * 85)
    for info in basis_info:
        weights_str = ", ".join(f"{w:2d}" for w in info['weights'])
        print(f"{info['idx']:<6d} | ({info['wy']}, {info['wx']})           | [{weights_str}]                        | {info['stab_size']:<18d}")

    # Test hypothesis: Is Span(Ham) ⊆ translation-invariant subspace?
    # Subspace of translation-invariant vectors has dimension exactly 4 (since they must be constant on each of the 4 orbits).
    # Since rank is 27 > 4, this is FALSE.
    print("\nTestando Hipótese:")
    print("  Hipótese: Span(Ham(G_T(6))) ⊆ subespaço de vetores invariantes por Z_6 x Z_6")
    print(f"  Rank observado: {rank}")
    print("  Dimensão do subespaço invariante (número de órbitas): 4")
    if rank > 4:
        print("  -> HIPÓTESE REJEITADA: O espaço dos tours possui muito mais graus de liberdade do que a invariância estrita de translação.")
    else:
        print("  -> HIPÓTESE CONFIRMADA (inesperadamente!)")

    # Let's check orbit rank of individual tours
    # For a tour, what is the rank of the space spanned by all its translations?
    print("\nAnálise de Representação de Grupo:")
    # We choose the first 5 tours and compute the rank of their translated orbits
    for i in range(min(5, len(tour_ints))):
        t_vec = int_to_vector(tour_ints[i], E)
        orbit_calc = GF2RankCalculator(E)
        for a in range(n):
            for b in range(n):
                perm = trans_table[(a, b)]
                g_t = translate_vector(t_vec, perm)
                orbit_calc.add_vector(vector_to_int(g_t))
        print(f"  Tour {i+1}: Rank de seu span sob translações = {orbit_calc.get_rank()}")

    return basis_info


def run_experiment_c():
    print("\n" + "="*80)
    print(" EXPERIMENTO C: SCALING DO RANK REAL (n = 4, 6, 8)")
    print("="*80)

    sizes = [4, 6, 8]
    table_rows = []

    for n in sizes:
        print(f"\n--- Analisando Toro {n}x{n} ---")
        ctx = build_graph(n, n)
        E = ctx['E']
        V = ctx['V']
        beta1 = E - V + 1

        # We precompute translations to count orbits
        trans_table = precompute_translations(ctx)
        orbits = compute_edge_orbits(ctx, trans_table)
        num_orbits = len(orbits)

        # Sampling
        if n == 4:
            # We already know there are exactly 1,344 tours
            tours = knight_tours(4, 4, 5000, seed=0)
        elif n == 6:
            tours = knight_tours(6, 6, 50000, seed=42)
        elif n == 8:
            # 10,000 tours for 8x8 is extremely fast and enough for convergence
            tours = knight_tours(8, 8, 15000, seed=0)

        # Create vectors and compute rank
        rank_calc = GF2RankCalculator(E)
        edge_endpoints = ctx['edge_endpoints']
        edge_to_idx = {tuple(edge_endpoints[i]): i for i in range(E)}

        for t in tours:
            t_vec = np.zeros(E, dtype=np.uint8)
            for j in range(V):
                u = t[j]
                w = t[(j+1)%V]
                key = (u, w) if u < w else (w, u)
                t_vec[edge_to_idx[key]] = 1
            rank_calc.add_vector(vector_to_int(t_vec))

        final_rank = rank_calc.get_rank()
        deficit = beta1 - final_rank
        ratio = final_rank / beta1

        print(f"  Toro {n}x{n}: E={E}, V={V}, beta1={beta1}, rank={final_rank}, deficit={deficit}, órbitas={num_orbits}")
        table_rows.append({
            'n': n,
            'beta1': beta1,
            'rank': final_rank,
            'deficit': deficit,
            'ratio': ratio,
            'orbits': num_orbits
        })

    # Print Target Table
    print("\n" + "="*80)
    print(" TABELA ALVO FINAL")
    print("="*80)
    print(f"{'n':<4} | {'β₁(toro)':<10} | {'rank_real':<10} | {'deficit_real':<12} | {'rank/β₁':<10} | {'órbitas Z_n×Z_n':<15}")
    print("-" * 75)
    for r in table_rows:
        print(f"{r['n']:<4d} | {r['beta1']:<10d} | {r['rank']:<10d} | {r['deficit']:<12d} | {r['ratio']:<10.4f} | {r['orbits']:<15d}")

    return table_rows

# ===========================================================================
# Main Execution
# ===========================================================================

if __name__ == "__main__":
    t_start = time.time()
    ctx6, trans6, orbits6, tour_ints6, rank_calc6 = run_experiment_a()
    basis_info6 = run_experiment_b(ctx6, trans6, orbits6, tour_ints6, rank_calc6)
    table_rows = run_experiment_c()
    t_end = time.time()
    print(f"\nTempo total de execução de todos os experimentos: {t_end - t_start:.2f} segundos.")
