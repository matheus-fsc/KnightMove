"""test_d4_orbital_branching.py — Teste diferencial da redução D4.

Critérios de aceite (todos obrigatórios):
  1. canonical_vertex_representatives(8) == exatamente 10 representantes.
  2. Enumeração BRUTA (sem simetria) de n=6 == 9862 tours (valor do paper).
  3. Enumeração COM redução D4 reconstrói exatamente 9862 (Σ tamanhos de órbita).
  4. Validação INDEPENDENTE por Burnside e por partição direta de órbitas
     (a partir dos conjuntos brutos): nº de órbitas e Σ tamanhos batem com
     a versão D4.
  5. Speedup real: nós explorados e tempo de parede, em >=3 seeds.

Se qualquer verificação falhar, o script PARA (AssertionError) e reporta o
mismatch específico — a implementação é ajustada, nunca o teste.
"""

# --- bootstrap: torna core/ importavel a partir de qualquer subpasta ---
import sys as _sys, pathlib as _pathlib
for _anc in _pathlib.Path(__file__).resolve().parents:
    if (_anc / "core").is_dir() and (_anc / "pyproject.toml").exists():
        _sys.path.insert(0, str(_anc / "core"))
        break
# --- fim do bootstrap ---
import sys

import numpy as np

import d4_symmetry as d4
from d4_orbital_branching import enumerate_tours, build_edge_perms
from knight_tours import build_graph

EXPECTED_N6 = 9862


def test_vertex_orbits():
    print("== Teste 1: órbitas de vértice (D4 puro) ==")
    reps = d4.canonical_vertex_representatives(8)
    print(f"  canonical_vertex_representatives(8) = {reps} (len={len(reps)})")
    assert len(reps) == 10, f"esperava 10 órbitas, obteve {len(reps)}"
    sizes = sorted(len(d4.orbit_of_vertex(r, 8)) for r in reps)
    print(f"  tamanhos de órbita: {sizes} (soma={sum(sizes)})")
    assert sizes.count(4) == 4 and sizes.count(8) == 6, \
        f"esperava 4×tam4 + 6×tam8, obteve {sizes}"
    assert sum(sizes) == 64
    # n=6 também precisa ser consistente (usado na enumeração)
    reps6 = d4.canonical_vertex_representatives(6)
    assert sum(len(d4.orbit_of_vertex(r, 6)) for r in reps6) == 36
    print("  OK\n")


def _orbit_partition(sets, perms):
    """Particiona uma lista de conjuntos de arestas em D4-órbitas.

    Retorna (n_orbitas, lista_de_tamanhos, ok_disjunto_cobre).
    """
    gperm = perms['gperm']
    remaining = set(sets)
    total = len(sets)
    sizes = []
    covered = 0
    while remaining:
        S = next(iter(remaining))
        orbit = set()
        for k in range(8):
            gk = gperm[k]
            orbit.add(frozenset(int(gk[e]) for e in S))
        # toda a órbita precisa estar presente no conjunto bruto
        if not orbit <= remaining:
            return len(sizes), sizes, False
        remaining -= orbit
        sizes.append(len(orbit))
        covered += len(orbit)
    return len(sizes), sizes, (covered == total)


def _burnside(sets, perms):
    """Nº de órbitas via lema de Burnside: (1/8) Σ_g |Fix(g)|."""
    gperm = perms['gperm']
    sset = set(sets)
    total_fixed = 0
    per_g = []
    for k in range(8):
        gk = gperm[k]
        fix = 0
        for S in sset:
            if frozenset(int(gk[e]) for e in S) == S:
                fix += 1
        per_g.append(fix)
        total_fixed += fix
    assert total_fixed % 8 == 0, f"Burnside não inteiro: {total_fixed}/8"
    return total_fixed // 8, per_g


def test_counts_and_orbits():
    print("== Teste 2/3/4: contagem bruta, D4 e validação independente (n=6) ==")
    # Enumeração bruta com conjuntos coletados (para validar órbitas)
    raw = enumerate_tours(6, use_symmetry=False, seed=0, collect_sets=True)
    print(f"  bruto:  count={raw['count']} leaves={raw['raw_leaves']} "
          f"nodes={raw['nodes']} t={raw['seconds']:.2f}s")
    assert raw['count'] == EXPECTED_N6, \
        f"BRUTO diverge: {raw['count']} != {EXPECTED_N6}"
    assert len(set(raw['sets'])) == EXPECTED_N6, \
        "conjuntos de arestas não são únicos (tour contado 2x?)"

    # Enumeração com D4
    sym = enumerate_tours(6, use_symmetry=True, seed=0, collect_sets=True)
    print(f"  D4:     count={sym['count']} orbitas(reps)={sym['raw_leaves']} "
          f"nodes={sym['nodes']} t={sym['seconds']:.2f}s")
    assert sym['count'] == EXPECTED_N6, \
        f"D4 reconstruído diverge: {sym['count']} != {EXPECTED_N6}"

    # Validação independente a partir dos conjuntos brutos
    perms = build_edge_perms(build_graph(6))
    n_orb, sizes, ok = _orbit_partition(raw['sets'], perms)
    print(f"  partição direta: {n_orb} órbitas, Σtam={sum(sizes)}, "
          f"disjunta&cobre={ok}")
    assert ok, "órbitas não disjuntas/não cobrem — ação de D4 inconsistente"
    assert sum(sizes) == EXPECTED_N6
    assert n_orb == sym['raw_leaves'], \
        f"nº de órbitas diverge: partição={n_orb} vs D4={sym['raw_leaves']}"
    assert sorted(sizes) == sorted(sym['orbit_sizes']), \
        "tamanhos de órbita divergem entre partição direta e D4"

    n_burn, per_g = _burnside(raw['sets'], perms)
    print(f"  Burnside: {n_burn} órbitas  (Fix por g = {per_g})")
    assert n_burn == n_orb, \
        f"Burnside diverge: {n_burn} vs partição {n_orb}"

    # os representantes da busca D4 devem ser um por órbita (canônicos, únicos)
    assert len(set(sym['sets'])) == sym['raw_leaves'], \
        "representantes D4 repetidos"
    print("  OK — 9862 confirmado por 3 vias independentes "
          "(bruto, D4, Burnside/partição)\n")
    return raw, sym


def test_speedup():
    print("== Teste 5: speedup real (3 seeds) ==")
    seeds = [0, 1, 2]
    print(f"  {'seed':>4} | {'nodes_raw':>10} {'nodes_D4':>10} {'x_nodes':>8} "
          f"| {'t_raw':>7} {'t_D4':>7} {'x_tempo':>8} | {'orbitas':>8}")
    agg = []
    for s in seeds:
        raw = enumerate_tours(6, use_symmetry=False, seed=s)
        sym = enumerate_tours(6, use_symmetry=True, seed=s)
        assert raw['count'] == EXPECTED_N6, f"seed {s} bruto={raw['count']}"
        assert sym['count'] == EXPECTED_N6, f"seed {s} D4={sym['count']}"
        xn = raw['nodes'] / sym['nodes']
        xt = raw['seconds'] / sym['seconds'] if sym['seconds'] > 0 else float('nan')
        agg.append((xn, xt, sym['raw_leaves']))
        print(f"  {s:>4} | {raw['nodes']:>10} {sym['nodes']:>10} {xn:>7.2f}x "
              f"| {raw['seconds']:>6.2f}s {sym['seconds']:>6.2f}s {xt:>7.2f}x "
              f"| {sym['raw_leaves']:>8}")
    mean_xn = np.mean([a[0] for a in agg])
    mean_xt = np.mean([a[1] for a in agg])
    print(f"  média: speedup nós = {mean_xn:.2f}x | speedup tempo = {mean_xt:.2f}x")
    assert all(a[2] == agg[0][2] for a in agg), \
        "nº de órbitas variou entre seeds (deveria ser invariante)"
    print(f"  nº de órbitas invariante entre seeds: {agg[0][2]}\n")
    return mean_xn, mean_xt


if __name__ == '__main__':
    test_vertex_orbits()
    raw, sym = test_counts_and_orbits()
    mean_xn, mean_xt = test_speedup()

    print("=" * 64)
    print("RESUMO")
    print(f"  tours (n=6): esperado {EXPECTED_N6}, obtido {sym['count']} "
          f"-> {'BATE' if sym['count'] == EXPECTED_N6 else 'DIVERGE'}")
    print(f"  órbitas D4: {sym['raw_leaves']}  "
          f"(fator de compactação {EXPECTED_N6 / sym['raw_leaves']:.2f}x)")
    print(f"  speedup médio: {mean_xn:.2f}x em nós, {mean_xt:.2f}x em tempo")
    print("=" * 64)
    print("TODOS OS TESTES PASSARAM")
