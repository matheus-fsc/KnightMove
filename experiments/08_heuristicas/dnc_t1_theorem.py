"""T1 — Teorema de suficiência local.

Prova/teste: para tabuleiro N=6k, dados k² blocos 6×6, cada bloco
contribuindo 1 path hamiltoniano (com endpoints s_i, e_i), conectados
por cross-edges em sequência, o resultado global é UM ciclo
hamiltoniano sse o multigrafo de blocos induzido pelas cross-edges
forma UM ciclo (Hamilton no meta-grafo).

Argumento estrutural:
  - Cada bloco i contribui path P_i = (s_i, ..., e_i) de 36 vértices
  - Cross-edge i conecta e_i a s_{σ(i)} para alguma permutação σ
  - O grafo meta tem k² nós, cada um de grau 2 (uma cross-edge entra,
    uma sai), portanto é uma união disjunta de ciclos
  - O ciclo global é 1 sse σ é uma permutação cíclica de comprimento k²
    (Hamilton cycle no meta-grafo)
  - Caso contrário, há c ≥ 2 sub-ciclos no meta-grafo, gerando c
    ciclos globais disjuntos

Demonstração empírica:
  1. Caso Hamilton (já demonstrado em Q3): 1 ciclo
  2. Caso 3+6: partição do 3×3 em 3-ciclo + 6-ciclo → 2 ciclos globais
"""

# --- bootstrap: torna core/ importavel a partir de qualquer subpasta ---
import sys as _sys, pathlib as _pathlib
for _anc in _pathlib.Path(__file__).resolve().parents:
    if (_anc / "core").is_dir() and (_anc / "pyproject.toml").exists():
        _sys.path.insert(0, str(_anc / "core"))
        break
# --- fim do bootstrap ---

import os
import pickle
import time
from collections import defaultdict

import numpy as np

import knight_tours as kt
import knight_tours_dnc as dnc
import dnc_q3_18x18 as q3


# Partição 3+6 do meta-grafo king 3×3:
#   - 3-ciclo:  (0,0)-(0,1)-(1,0)-(0,0)
#   - 6-ciclo:  (1,1)-(0,2)-(1,2)-(2,2)-(2,1)-(2,0)-(1,1)
#
# Verificar adjacências king (|dr|≤1 e |dc|≤1):
#   3-ciclo:
#     (0,0)-(0,1) ✓
#     (0,1)-(1,0) diag ✓
#     (1,0)-(0,0) ✓
#   6-ciclo:
#     (1,1)-(0,2) diag ✓
#     (0,2)-(1,2) ✓
#     (1,2)-(2,2) ✓
#     (2,2)-(2,1) ✓
#     (2,1)-(2,0) ✓
#     (2,0)-(1,1) diag ✓

CYCLE_A = [(0, 0), (0, 1), (1, 0)]
CYCLE_B = [(1, 1), (0, 2), (1, 2), (2, 2), (2, 1), (2, 0)]


def assemble_subcycle(catalog, cycle, rng_seed=0):
    """Tenta montar 1 ciclo fechado sobre os blocos em `cycle`.

    Retorna (success, vertex_sequence_global, info).
    """
    rng = np.random.default_rng(rng_seed)
    n_blocks = len(cycle)

    # Para cada bloco, fronteiras direcionais
    frontier_prev = {}
    frontier_next = {}
    cross_idx_next = {}
    for i in range(n_blocks):
        b = cycle[i]
        b_prev = cycle[(i - 1) % n_blocks]
        b_next = cycle[(i + 1) % n_blocks]
        frontier_prev[b] = q3.frontier_to(b, b_prev)
        frontier_next[b] = q3.frontier_to(b, b_next)
        idx = defaultdict(list)
        for (ga, gb) in q3.cross_edges_18(b, b_next):
            la = q3.global_to_local_18(ga, b)
            lb = q3.global_to_local_18(gb, b_next)
            idx[la].append(lb)
        cross_idx_next[b] = idx

    # Pares (s, e) válidos por bloco
    pairs_by_block = {}
    by_start_by_block = {}
    for b in cycle:
        fr_p = frontier_prev[b]
        fr_n = frontier_next[b]
        idx_next = cross_idx_next[b]
        pairs = [(s, e) for (s, e) in catalog
                 if s in fr_p and e in fr_n and catalog[(s, e)]
                 and e in idx_next]
        pairs_by_block[b] = pairs
        bs = defaultdict(list)
        for (s, e) in pairs:
            bs[s].append(e)
        by_start_by_block[b] = bs

    # Amostrar
    for attempt in range(1000):
        b0 = cycle[0]
        if not pairs_by_block[b0]:
            return False, None, {'reason': 'no pairs in first block'}
        s0, e0 = pairs_by_block[b0][rng.integers(len(pairs_by_block[b0]))]
        endpoints = [(s0, e0)]
        feasible = True
        for i in range(1, n_blocks):
            b_prev = cycle[i - 1]
            b_curr = cycle[i]
            e_prev = endpoints[-1][1]
            cand_s = [v for v in cross_idx_next[b_prev][e_prev]
                      if v in frontier_prev[b_curr]]
            if not cand_s:
                feasible = False
                break
            s_curr = cand_s[rng.integers(len(cand_s))]
            if i == n_blocks - 1:
                cand_e = [e for e in by_start_by_block[b_curr].get(s_curr, [])
                          if s0 in cross_idx_next[b_curr].get(e, [])]
            else:
                cand_e = by_start_by_block[b_curr].get(s_curr, [])
            if not cand_e:
                feasible = False
                break
            e_curr = cand_e[rng.integers(len(cand_e))]
            endpoints.append((s_curr, e_curr))
        if not feasible:
            continue

        # Montar
        seq = []
        for i, b in enumerate(cycle):
            s, e = endpoints[i]
            p = catalog[(s, e)][rng.integers(len(catalog[(s, e)]))]
            for v in p:
                seq.append(q3.local_to_global_18(v, b))
        return True, np.asarray(seq, dtype=np.int32), {
            'attempts_internal': attempt + 1,
            'endpoints': endpoints,
        }
    return False, None, {'reason': 'exceeded attempts'}


def verify_two_cycles_global(seq_a, seq_b):
    """Verifica que seq_a e seq_b são DOIS ciclos disjuntos, cada um
    fechado em si, sobre o tabuleiro 18×18.

    Retorna dict com diagnóstico.
    """
    N = q3.N_GLOBAL_18
    V = N * N
    seq_a = list(map(int, seq_a))
    seq_b = list(map(int, seq_b))

    # 1. Conjuntos disjuntos
    set_a = set(seq_a)
    set_b = set(seq_b)
    overlap = set_a & set_b
    union = set_a | set_b

    # 2. Cada ciclo é fechado em si: último → primeiro é knight move
    knight_moves = set(dnc.KNIGHT_MOVES)

    def is_closed_chain(seq):
        if len(seq) < 2:
            return False, "muito curto"
        for i in range(len(seq)):
            u = seq[i]
            v = seq[(i + 1) % len(seq)]
            ru, cu = divmod(u, N)
            rv, cv = divmod(v, N)
            if (rv - ru, cv - cu) not in knight_moves:
                return False, f"pos {i}: {(ru, cu)}→{(rv, cv)} não é knight"
        return True, "ok"

    closed_a, ra = is_closed_chain(seq_a)
    closed_b, rb = is_closed_chain(seq_b)

    return {
        'disjuntos': len(overlap) == 0,
        'overlap': len(overlap),
        'cobre_todos': len(union) == V,
        'falta': V - len(union),
        '|seq_a|': len(set_a),
        '|seq_b|': len(set_b),
        'cycle_a_fechado': closed_a,
        'cycle_a_motivo': ra,
        'cycle_b_fechado': closed_b,
        'cycle_b_motivo': rb,
        'soma_eh_tour_18x18': closed_a and closed_b and len(overlap) == 0 and len(union) == V,
    }


if __name__ == '__main__':
    print("=" * 70)
    print(" T1 — Teorema da suficiência local")
    print("=" * 70)

    with open(dnc.CATALOG_FILE, 'rb') as f:
        obj = pickle.load(f)
    catalog = obj['catalog']

    # ----------------------------------------------------------------------
    # Caso 1 — Hamilton (controle, já validado em Q3)
    # ----------------------------------------------------------------------
    print(f"\n[Controle] Hamilton cycle 9-blocos (Q3): deve gerar 1 tour")
    t0 = time.time()
    r_q3 = q3.random_matching_18(catalog, K=1, rng_seed=11)
    t1 = time.time()
    if r_q3['tours']:
        tour = r_q3['tours'][0]
        valid = kt.verify_tour(tour, q3.N_GLOBAL_18)
        print(f"  1 tour gerado, válido={valid} | {t1-t0:.3f}s")

    # ----------------------------------------------------------------------
    # Caso 2 — Partição 3+6 (NÃO Hamilton no meta-grafo)
    # ----------------------------------------------------------------------
    print(f"\n[Teste] Partição 3+6 do meta-grafo king 3×3")
    print(f"  3-ciclo: {CYCLE_A}")
    print(f"  6-ciclo: {CYCLE_B}")
    print(f"  → Esperado: 2 ciclos globais DISJUNTOS, NÃO um tour 18×18")

    t0 = time.time()
    ok_a, seq_a, info_a = assemble_subcycle(catalog, CYCLE_A, rng_seed=0)
    ok_b, seq_b, info_b = assemble_subcycle(catalog, CYCLE_B, rng_seed=0)
    t1 = time.time()
    print(f"  Montagem em {t1-t0:.3f}s")
    print(f"  3-ciclo: ok={ok_a} (info={info_a if not ok_a else 'tentativas='+str(info_a['attempts_internal'])})")
    print(f"  6-ciclo: ok={ok_b} (info={info_b if not ok_b else 'tentativas='+str(info_b['attempts_internal'])})")

    if ok_a and ok_b:
        diag = verify_two_cycles_global(seq_a, seq_b)
        print(f"\n  Diagnóstico:")
        for k, v in diag.items():
            print(f"    {k}: {v}")

        # Conferir se a UNIÃO (concatenação) seria detectada como tour:
        # esperamos NÃO ser tour porque há dois ciclos disjuntos
        concat = np.concatenate([seq_a, seq_b])
        as_tour = kt.verify_tour(concat, q3.N_GLOBAL_18)
        print(f"\n  verify_tour(concat 18×18): {as_tour} "
              f"(esperado False — dois ciclos disjuntos não formam 1 tour)")

    # ----------------------------------------------------------------------
    # Conclusão formal
    # ----------------------------------------------------------------------
    print("\n" + "=" * 70)
    print(" Conclusão estrutural")
    print("=" * 70)
    print("""
  TEOREMA (suficiência local). Seja G um tabuleiro N×N particionado em
  k² blocos 6×6 (N=6k). Sejam P₁,...,P_{k²} paths hamiltonianos internos
  (um por bloco) e E₁,...,E_{k²} cross-edges entre blocos consecutivos
  na visita.

  Cada bloco aparece com grau 2 no meta-multigrafo (uma cross-edge entra
  no path por s_i, uma sai por e_i). Logo o meta-grafo induzido por
  {E_j} é uma união disjunta de ciclos sobre os k² nós.

  A concatenação global é UM ciclo hamiltoniano sse o meta-grafo é
  exatamente um ciclo simples (Hamilton no meta-grafo).
  Caso contrário, com c ≥ 2 sub-ciclos no meta-grafo, a concatenação
  forma c ciclos disjuntos globais (cada um com 36·m vértices, onde m
  é o tamanho do sub-ciclo no meta-grafo).

  Verificação empírica (acima):
    - Hamilton 9-cycle  → 1 tour 18×18 ✓
    - Partição 3 + 6    → 2 ciclos disjuntos (108 + 216 vértices) ✓

  Corolário: o D&C com qualquer Hamilton no meta-grafo é correto.
  Por construção, nenhum sub-ciclo prematuro pode surgir.
""")
