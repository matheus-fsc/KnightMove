"""
Benchmark: Z3 puro vs Z3+mandatory vs backtracking+propagação no 10×10.

Compara três métodos para encontrar K=50 tours fechados:
  A) Z3 puro com grau-2 + sub-tour elimination via BFS
  B) Z3 + 8 obrigatórias pré-fixadas (R1)
  C) Backtracking + R1..R6 (de backtracking_10x10.py)

Reaproveita o resultado do C salvo em data/backtracking_results.json.
"""

from __future__ import annotations

import argparse
import json
import time
from collections import deque
from pathlib import Path

import numpy as np

from propagation_engine_10x10 import (
    PropagationEngine10,
    ROOT,
    carregar_dados,
)

try:
    from z3 import (
        And,
        Bool,
        Implies,
        Not,
        Or,
        PbEq,
        Solver,
        Sum,
        sat,
        unsat,
    )
    HAS_Z3 = True
except Exception:
    HAS_Z3 = False


def bfs_tour(arestas_em_um, edges_uv, n_vertices):
    adj = [[] for _ in range(n_vertices)]
    for e in arestas_em_um:
        u, v = edges_uv[e]
        adj[u].append(v)
        adj[v].append(u)
    bfs = deque([0])
    seen = {0}
    while bfs:
        x = bfs.popleft()
        for y in adj[x]:
            if y not in seen:
                seen.add(y)
                bfs.append(y)
    return len(seen) == n_vertices, seen


def z3_busca(dados, alvo: int, com_mandatory: bool, timeout_s: float) -> dict:
    if not HAS_Z3:
        return {"erro": "z3 não disponível"}

    edges_uv = dados["edges_uv"]
    n_arestas = len(edges_uv)
    n_vertices = max(max(u, v) for u, v in edges_uv) + 1

    # v2e
    v2e: dict[int, list[int]] = {}
    for ei, (u, v) in enumerate(edges_uv):
        v2e.setdefault(u, []).append(ei)
        v2e.setdefault(v, []).append(ei)

    x = [Bool(f"x_{i}") for i in range(n_arestas)]
    solver = Solver()
    # grau-2 em cada vértice (PbEq = pseudo-boolean equality)
    for v, es in v2e.items():
        solver.add(PbEq([(x[e], 1) for e in es], 2))

    if com_mandatory:
        for e in dados["mandatory_idx"]:
            solver.add(x[e])

    t0 = time.perf_counter()
    tours: list[list[int]] = []
    t_first = None
    tentativas = 0
    tempo_z3_total = 0.0

    while len(tours) < alvo:
        if time.perf_counter() - t0 > timeout_s:
            break
        tentativas += 1
        ts = time.perf_counter()
        if solver.check() != sat:
            tempo_z3_total += time.perf_counter() - ts
            break
        tempo_z3_total += time.perf_counter() - ts
        modelo = solver.model()
        ativos = []
        for e in range(n_arestas):
            v = modelo[x[e]]
            if v is not None and bool(v):
                ativos.append(e)
        is_conexo, seen = bfs_tour(ativos, edges_uv, n_vertices)
        if is_conexo:
            tours.append(ativos)
            if t_first is None:
                t_first = time.perf_counter() - t0
            # bloquear este tour
            solver.add(Sum([x[e] if e in ativos else Not(x[e]) for e in range(n_arestas)])
                       < n_arestas)
            # forma equivalente mais eficiente: bloquear pelo menos uma diferença
            # (usando que ativos completamente determina)
            # já adicionada acima
        else:
            # bloquear sub-tour: pelo menos uma aresta dentro do componente
            # do BFS deve mudar (impedindo recurrência exata do sub-tour)
            arestas_componente = [
                e for e in range(n_arestas)
                if edges_uv[e][0] in seen and edges_uv[e][1] in seen
                and modelo[x[e]] is not None and bool(modelo[x[e]])
            ]
            if arestas_componente:
                solver.add(Or([Not(x[e]) for e in arestas_componente]))
            else:
                # fallback: bloquear esta solução exata
                solver.add(Not(And([x[e] for e in ativos])))

    t_total = time.perf_counter() - t0
    return {
        "alvo": alvo,
        "tours_encontrados": len(tours),
        "t_first": t_first,
        "t_total": t_total,
        "tempo_z3_check_total_s": tempo_z3_total,
        "tentativas": tentativas,
        "com_mandatory": com_mandatory,
    }


def carregar_resultado_backtracking() -> dict:
    p = ROOT / "data" / "backtracking_results.json"
    if not p.exists():
        return {"erro": "backtracking_results.json ausente — rode backtracking_10x10.py primeiro"}
    r = json.loads(p.read_text())
    # encontrar t_first do log
    t_first = None
    for dt, nos, tours, _ in r.get("log_progresso", []):
        if tours >= 1:
            t_first = dt
            break
    return {
        "alvo": r["alvo_tours"],
        "tours_encontrados": r["n_tours"],
        "t_first": t_first,
        "t_total": r["tempo_s"],
        "nos_explorados": r["n_nos"],
        "nos_por_tour": r["nos_por_tour"],
        "podas": r["podas"],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--alvo", type=int, default=50)
    parser.add_argument("--timeout", type=float, default=1800.0)
    parser.add_argument("--pular_z3", action="store_true")
    args = parser.parse_args()

    dados = carregar_dados(strict_pair_mode="todos")
    saida = {"alvo": args.alvo, "timeout": args.timeout}

    if not args.pular_z3:
        print(f"\n=== A) Z3 puro (alvo={args.alvo}) ===")
        a = z3_busca(dados, alvo=args.alvo, com_mandatory=False, timeout_s=args.timeout)
        print(json.dumps(a, indent=2))
        saida["A_z3_puro"] = a

        print(f"\n=== B) Z3 + 8 obrigatórias (alvo={args.alvo}) ===")
        b = z3_busca(dados, alvo=args.alvo, com_mandatory=True, timeout_s=args.timeout)
        print(json.dumps(b, indent=2))
        saida["B_z3_mandatory"] = b

    print("\n=== C) Backtracking + R1..R6 ===")
    c = carregar_resultado_backtracking()
    print(json.dumps(c, indent=2))
    saida["C_backtracking"] = c

    out = ROOT / "data" / "benchmark_comparison.json"
    out.write_text(json.dumps(saida, indent=2))
    print(f"\nSalvo em {out}")


if __name__ == "__main__":
    main()
