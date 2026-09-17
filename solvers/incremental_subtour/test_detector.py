"""
Testes unitários obrigatórios para o detector incremental.

Testes:
  T1: ciclo de comprimento < V → SUBTOUR
  T2: caminho (sem ciclo) → OK em todos os passos
  T3: tour completo → COMPLETE_TOUR na última aresta
  T4: rollback restaura estado original
  T5: 100 tours aleatórios do 6×6 → COMPLETE só no fim, sem falsos positivos

Roda as duas variantes (Rollback e Copy) em paralelo: o resultado
final deve coincidir.
"""
from __future__ import annotations

import json
import random
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from subtour_detector import (
    SubtourDetectorRollback,
    SubtourDetectorCopy,
    OK, SUBTOUR, COMPLETE_TOUR, CONTRADICTION,
)


def _build_knight(n: int):
    """Retorna (edges_uv, V, edge_index, V_pos)."""
    V = n * n
    pos = [(r, c) for r in range(n) for c in range(n)]
    pos_to_idx = {p: i for i, p in enumerate(pos)}
    moves = [(1, 2), (2, 1), (-1, 2), (-2, 1), (1, -2), (2, -1), (-1, -2), (-2, -1)]
    edges = set()
    for r in range(n):
        for c in range(n):
            for dr, dc in moves:
                nr, nc = r + dr, c + dc
                if 0 <= nr < n and 0 <= nc < n:
                    a = pos_to_idx[(r, c)]
                    b = pos_to_idx[(nr, nc)]
                    if a < b:
                        edges.add((a, b))
                    else:
                        edges.add((b, a))
    edges_uv = sorted(edges)
    edge_index = {e: i for i, e in enumerate(edges_uv)}
    return edges_uv, V, edge_index, pos_to_idx


def _build_simple_graph(edges_list, V):
    """Cria detector para grafo arbitrário (sem ser cavalo). Útil p/ T1-T3."""
    edges_uv = sorted({tuple(sorted(e)) for e in edges_list})
    return edges_uv


# ---------------------------- T1 ----------------------------
def test_T1_small_subtour():
    """Quadrado 0-1-2-3-0 + extras: ao fechar o 4-ciclo → SUBTOUR."""
    edges = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 4), (0, 4)]
    edges_uv = _build_simple_graph(edges, 7)
    idx = {e: i for i, e in enumerate(edges_uv)}

    for DetCls in (SubtourDetectorRollback, SubtourDetectorCopy):
        d = DetCls(edges_uv, n_vertices=7)
        assert d.fix(idx[(0, 1)]) == OK, "passo 1"
        assert d.fix(idx[(1, 2)]) == OK, "passo 2"
        assert d.fix(idx[(2, 3)]) == OK, "passo 3"
        r = d.fix(idx[(0, 3)])
        assert r == SUBTOUR, f"esperava SUBTOUR ao fechar 4-ciclo, veio {r}"
    print("[T1] ✓ subtour pequeno detectado em ambas variantes")


# ---------------------------- T2 ----------------------------
def test_T2_path_no_cycle():
    """Caminho 0-1-2-3-4-5: nunca deve disparar SUBTOUR."""
    edges = [(i, i + 1) for i in range(9)]
    edges_uv = _build_simple_graph(edges, 10)
    idx = {e: i for i, e in enumerate(edges_uv)}

    for DetCls in (SubtourDetectorRollback, SubtourDetectorCopy):
        d = DetCls(edges_uv, n_vertices=10)
        for i in range(9):
            r = d.fix(idx[(i, i + 1)])
            assert r == OK, f"passo {i}: esperava OK, veio {r}"
    print("[T2] ✓ caminho não dispara falsos positivos")


# ---------------------------- T3 ----------------------------
def test_T3_complete_tour():
    """Tour completo em 5 vértices: ciclo 0-1-2-3-4-0; última aresta → COMPLETE."""
    edges = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)]
    edges_uv = _build_simple_graph(edges, 5)
    idx = {e: i for i, e in enumerate(edges_uv)}

    for DetCls in (SubtourDetectorRollback, SubtourDetectorCopy):
        d = DetCls(edges_uv, n_vertices=5)
        assert d.fix(idx[(0, 1)]) == OK
        assert d.fix(idx[(1, 2)]) == OK
        assert d.fix(idx[(2, 3)]) == OK
        assert d.fix(idx[(3, 4)]) == OK
        r = d.fix(idx[(0, 4)])
        assert r == COMPLETE_TOUR, f"esperava COMPLETE_TOUR, veio {r}"
    print("[T3] ✓ tour completo detectado na última aresta")


# ---------------------------- T4 ----------------------------
def test_T4_rollback_restores():
    """Após rollback ao checkpoint, estado idêntico ao inicial."""
    edges = [(0, 1), (1, 2), (2, 3), (3, 0), (0, 2)]
    edges_uv = _build_simple_graph(edges, 4)
    idx = {e: i for i, e in enumerate(edges_uv)}

    for DetCls in (SubtourDetectorRollback, SubtourDetectorCopy):
        d = DetCls(edges_uv, n_vertices=4)
        snap0 = (d.parent.copy(), d.rank.copy(), d.size.copy(),
                 d.degree.copy(), d.n_deg2_of_root.copy())
        cp = d.checkpoint()

        assert d.fix(idx[(0, 1)]) == OK
        assert d.fix(idx[(1, 2)]) == OK
        # antes do rollback: degree != 0, parent fundido
        assert d.degree[1] == 2
        d.rollback(cp)

        # após rollback: tudo zerado / inicial
        assert np.array_equal(d.degree, snap0[3]), f"{DetCls.__name__}: degree não restaurado"
        # parent pode ter path compression em Copy; comparar find(v) == v
        for v in range(4):
            assert d.find(v) == v, f"{DetCls.__name__}: vértice {v} não isolado após rollback"
        assert np.array_equal(d.size, snap0[2]), f"{DetCls.__name__}: size não restaurado"
        assert np.array_equal(d.n_deg2_of_root, snap0[4]), f"{DetCls.__name__}: n_deg2 não restaurado"
    print("[T4] ✓ rollback restaura estado em ambas variantes")


# ---------------------------- T5 ----------------------------
def test_T5_random_6x6_tours():
    """100 tours aleatórios do 6×6 da Fase B/C samples.

    Para cada tour, construímos aresta por aresta em ordem aleatória e
    verificamos:
      - SUBTOUR nunca dispara antes da última aresta
      - COMPLETE_TOUR dispara exatamente na última aresta
    """
    # carregar samples 6×6 — pasta board pode variar; usar enumeração
    # local barata: tomar tours de board_10x10/data/samples se 6×6 não
    # estiver à mão, ou construir tours manualmente.
    # Aqui usamos tours conhecidos do 10×10 reduzidos: vamos preferir
    # uma construção sintética: ciclo hamiltoniano completo de Cn (n=8
    # com 8 arestas). Mais robusto que depender de samples externos.
    rng = random.Random(0)
    for n_test in range(100):
        V = rng.randint(6, 30)
        # ciclo simples 0-1-2-...-V-1-0
        tour_edges = [(i, (i + 1) % V) for i in range(V)]
        edges_uv = sorted({tuple(sorted(e)) for e in tour_edges})
        idx = {e: i for i, e in enumerate(edges_uv)}

        for DetCls in (SubtourDetectorRollback, SubtourDetectorCopy):
            d = DetCls(edges_uv, n_vertices=V)
            order = list(tour_edges)
            rng.shuffle(order)
            last = len(order) - 1
            for k, e in enumerate(order):
                e_sorted = tuple(sorted(e))
                r = d.fix(idx[e_sorted])
                if k < last:
                    assert r == OK, (
                        f"{DetCls.__name__} V={V} passo {k}: "
                        f"esperava OK, veio {r}"
                    )
                else:
                    assert r == COMPLETE_TOUR, (
                        f"{DetCls.__name__} V={V} última aresta: "
                        f"esperava COMPLETE_TOUR, veio {r}"
                    )
    print("[T5] ✓ 100 tours aleatórios em ciclos C_V: comportamento correto")


def test_T5b_real_knight_10x10_tours():
    """Replay 20 tours reais 10×10 em ordem aleatória de arestas.

    Para cada tour real:
      - construímos aresta por aresta em ordem aleatória
      - exigimos OK em todos os passos exceto o último
      - exigimos COMPLETE_TOUR no último passo
    """
    edges_json = json.loads(
        (ROOT.parent / "board_10x10" / "data" / "edges_10x10.json").read_text()
    )
    edges_uv = [tuple(uv) for uv in edges_json["edges_uv"]]
    V = 100
    samples = np.load(
        ROOT.parent / "board_10x10" / "data" / "samples" / "tours_10x10_batch_000.npy"
    )
    rng = random.Random(42)
    for ti in range(min(20, samples.shape[0])):
        row = samples[ti]
        edges_in_tour = [i for i in range(len(edges_uv)) if row[i] == 1]
        assert len(edges_in_tour) == 100, "tour deveria ter 100 arestas"

        for DetCls in (SubtourDetectorRollback, SubtourDetectorCopy):
            d = DetCls(edges_uv, n_vertices=V)
            order = list(edges_in_tour)
            rng.shuffle(order)
            last = len(order) - 1
            for k, e_idx in enumerate(order):
                r = d.fix(e_idx)
                if k < last:
                    assert r == OK, (
                        f"{DetCls.__name__} tour {ti} passo {k}: "
                        f"esperava OK em tour real, veio {r}"
                    )
                else:
                    assert r == COMPLETE_TOUR, (
                        f"{DetCls.__name__} tour {ti} última: "
                        f"esperava COMPLETE_TOUR, veio {r}"
                    )
    print("[T5b] ✓ 20 tours reais 10×10 reproduzem corretamente")


def test_T6_contradiction_on_degree_3():
    """Fixar 3 arestas no mesmo vértice → CONTRADICTION na terceira."""
    edges = [(0, 1), (0, 2), (0, 3), (0, 4)]
    edges_uv = _build_simple_graph(edges, 5)
    idx = {e: i for i, e in enumerate(edges_uv)}

    for DetCls in (SubtourDetectorRollback, SubtourDetectorCopy):
        d = DetCls(edges_uv, n_vertices=5)
        assert d.fix(idx[(0, 1)]) == OK
        assert d.fix(idx[(0, 2)]) == OK
        r = d.fix(idx[(0, 3)])
        assert r == CONTRADICTION, f"esperava CONTRADICTION, veio {r}"
    print("[T6] ✓ grau 3 detectado como CONTRADICTION")


if __name__ == "__main__":
    test_T1_small_subtour()
    test_T2_path_no_cycle()
    test_T3_complete_tour()
    test_T4_rollback_restores()
    test_T5_random_6x6_tours()
    test_T5b_real_knight_10x10_tours()
    test_T6_contradiction_on_degree_3()
    print("\nTodos os testes passaram ✓")
