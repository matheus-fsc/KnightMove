#!/usr/bin/env python3
"""
sample_tours.py
===============
Amostragem SAT não-viesada de ciclos hamiltonianos do cavalo 10×10.

Formulação:
  - Variáveis booleanas x_e para cada aresta
  - Restrição grau-2 em cada vértice: Σ_{e ∋ v} x_e == 2
  - Quebra de simetria aleatória por sample (forçar uma aresta a 1 ou 0)
  - Pós-filtro de conectividade via BFS — se a solução grau-2 é união de
    sub-tours, rejeita e adiciona corte (no-good) para essa configuração
  - break_symmetry=False (sem viés)

Cada amostra é serializada como assinatura binária sobre as 288 arestas
(uint8, shape (E,)). Um batch agrega K amostras válidas em (K, E).

Saídas:
  data/samples/tours_10x10_batch_{i:03d}.npy   shape (K_validas, E)
  data/samples/run_log.json                    estatísticas por batch
"""

import argparse
import json
import random
import sys
import time
from pathlib import Path

import numpy as np
from z3 import Bool, Or, PbEq, Solver, sat, is_true

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
SAMPLES = DATA / "samples"

# importa o grafo do módulo irmão
sys.path.insert(0, str(ROOT))
from graph_10x10 import (build_graph, label, vid, TOTAL, BOARD)  # noqa: E402


# ── conectividade via BFS sobre arestas ativas ──────────────────────────

def is_single_tour(active_edges, n_vertices=TOTAL):
    """Verifica se o conjunto de arestas ativas (grau-2 por vértice) forma
    um único ciclo hamiltoniano em vez de união de sub-tours."""
    adj = [[] for _ in range(n_vertices)]
    for u, v in active_edges:
        adj[u].append(v)
        adj[v].append(u)
    if any(len(a) != 2 for a in adj):
        return False, []
    visited = [False] * n_vertices
    tour = []
    cur = 0
    prev = -1
    visited[0] = True
    tour.append(0)
    for _ in range(n_vertices - 1):
        nxt = adj[cur][0] if adj[cur][0] != prev else adj[cur][1]
        if visited[nxt]:
            return False, []
        visited[nxt] = True
        tour.append(nxt)
        prev, cur = cur, nxt
    # fechamento
    return (adj[cur][0] == 0 or adj[cur][1] == 0), tour


# ── construção do solver base (reutilizável entre amostras) ─────────────

def build_base_solver(edges):
    """
    Solver com restrição de grau-2 em cada vértice. Não inclui ainda
    a quebra de simetria aleatória (essa muda por sample) nem cortes.
    Retorna (solver, xvars, vertex_inc).

    xvars: lista de Bool, um por aresta
    vertex_inc: dict v -> lista de índices de arestas incidentes em v
    """
    E = len(edges)
    xvars = [Bool(f"x_{i}") for i in range(E)]

    vertex_inc = {v: [] for v in range(TOTAL)}
    for i, (u, v) in enumerate(edges):
        vertex_inc[u].append(i)
        vertex_inc[v].append(i)

    s = Solver()
    for v, inc in vertex_inc.items():
        s.add(PbEq([(xvars[i], 1) for i in inc], 2))

    return s, xvars, vertex_inc


# ── amostragem ─────────────────────────────────────────────────────────

def sample_batch(edges, batch_size, *, batch_idx, seed,
                 timeout_per_sample_s=30, max_sym_attempts_per_sample=4,
                 verbose=False):
    """Retorna (lista de assinaturas válidas, stats)."""
    rng = random.Random(seed)
    E = len(edges)

    solver, xvars, vertex_inc = build_base_solver(edges)

    sigs = []
    stats = {
        "batch_idx": batch_idx,
        "batch_size_target": batch_size,
        "n_valid": 0,
        "n_rejected_subtours": 0,
        "n_unsat_local": 0,
        "n_attempts": 0,
        "elapsed_s": 0.0,
        "time_per_valid": [],
    }

    t0 = time.perf_counter()

    while len(sigs) < batch_size:
        t_attempt = time.perf_counter()
        stats["n_attempts"] += 1

        # quebra de simetria aleatória: força uma aresta a 1 OU 0 com 50%/50%
        # usamos um push/pop para que a cláusula só valha nesta tentativa
        solver.push()

        # tentativas de simetria — se a primeira der UNSAT, tenta outras
        local_unsat = True
        for sym_try in range(max_sym_attempts_per_sample):
            edge_pick = rng.randrange(E)
            val_pick = rng.randrange(2)
            solver.push()
            solver.add(xvars[edge_pick] == (val_pick == 1))
            solver.set("timeout", int(timeout_per_sample_s * 1000))
            res = solver.check()
            if res == sat:
                local_unsat = False
                break
            solver.pop()  # tira a quebra de simetria que deu UNSAT

        if local_unsat:
            stats["n_unsat_local"] += 1
            solver.pop()  # tira o push externo
            # remove a cláusula que travou — adiciona corte permanente? Não:
            # apenas continua, próxima quebra aleatória pode ser feliz
            continue

        m = solver.model()
        active = []
        sig = np.zeros(E, dtype=np.uint8)
        for i, (u, v) in enumerate(edges):
            if is_true(m.evaluate(xvars[i])):
                active.append((u, v))
                sig[i] = 1

        # pós-filtro de conectividade
        ok, tour = is_single_tour(active)
        # pop as duas camadas (simetria + outer)
        solver.pop()
        solver.pop()

        if ok:
            sigs.append(sig)
            t_now = time.perf_counter() - t_attempt
            stats["time_per_valid"].append(round(t_now, 4))
            stats["n_valid"] += 1
            if verbose and len(sigs) % 25 == 0:
                elapsed = time.perf_counter() - t0
                rej_rate = stats["n_rejected_subtours"] / max(stats["n_attempts"], 1)
                print(f"    [batch {batch_idx}] {len(sigs)}/{batch_size}  "
                      f"attempts={stats['n_attempts']}  "
                      f"reject={rej_rate:.1%}  "
                      f"elapsed={elapsed:.1f}s",
                      flush=True)
            # adiciona corte global de exclusão dessa configuração
            solver.add(Or([xvars[i] != bool(sig[i]) for i in range(E)]))
        else:
            stats["n_rejected_subtours"] += 1
            # corte permanente: essa configuração de sub-tour não pode reaparecer
            solver.add(Or([xvars[i] != bool(sig[i]) for i in range(E)]))

    stats["elapsed_s"] = round(time.perf_counter() - t0, 3)
    if stats["time_per_valid"]:
        stats["mean_time_per_valid_s"] = round(
            sum(stats["time_per_valid"]) / len(stats["time_per_valid"]), 4)
    else:
        stats["mean_time_per_valid_s"] = None
    stats["rejection_rate"] = round(
        stats["n_rejected_subtours"] / max(stats["n_attempts"], 1), 4)

    return sigs, stats


# ── runner ──────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--n-batches", type=int, default=20)
    p.add_argument("--batch-size", type=int, default=500)
    p.add_argument("--timeout-per-sample-s", type=int, default=30)
    p.add_argument("--seed", type=int, default=2026)
    p.add_argument("--test", action="store_true",
                   help="Modo teste: 2 batches de 500 (1000 amostras)")
    p.add_argument("--max-sym-attempts", type=int, default=4)
    args = p.parse_args()

    if args.test:
        args.n_batches = 2

    SAMPLES.mkdir(parents=True, exist_ok=True)

    adj, edges = build_graph()
    E = len(edges)
    print(f"Grafo 10×10  V={TOTAL}  E={E}")
    print(f"Plano: {args.n_batches} batches × {args.batch_size} amostras "
          f"(timeout {args.timeout_per_sample_s}s/sample)")
    print()

    all_stats = []
    total_t = time.perf_counter()
    for i in range(args.n_batches):
        print(f"── batch {i:03d} ─────────────────────────")
        t0 = time.perf_counter()
        sigs, stats = sample_batch(
            edges, args.batch_size, batch_idx=i,
            seed=args.seed + i * 1000,
            timeout_per_sample_s=args.timeout_per_sample_s,
            max_sym_attempts_per_sample=args.max_sym_attempts,
            verbose=True,
        )
        elapsed = time.perf_counter() - t0
        M = np.stack(sigs, axis=0) if sigs else np.empty((0, E), dtype=np.uint8)
        out = SAMPLES / f"tours_10x10_batch_{i:03d}.npy"
        np.save(out, M)

        # avaliação rápida de unicidade local
        n_unique = len({tuple(s) for s in sigs})
        stats["n_unique_local"] = n_unique
        stats["mean_time_per_valid_s"] = stats.get("mean_time_per_valid_s")
        all_stats.append(stats)

        print(f"  → {len(sigs)} válidas  únicas={n_unique}  "
              f"rej={stats['rejection_rate']:.1%}  "
              f"unsat_local={stats['n_unsat_local']}  "
              f"tempo={elapsed:.1f}s  "
              f"mean/valid={stats.get('mean_time_per_valid_s')}s")
        print(f"  Salvo: data/samples/{out.name}\n")

        # checkpoint de segurança
        mean_time = stats.get("mean_time_per_valid_s") or 999
        if mean_time > args.timeout_per_sample_s:
            print(f"  ALERTA: tempo médio por amostra ({mean_time}s) excede "
                  f"timeout. Abortando.")
            break

    total_elapsed = time.perf_counter() - total_t
    log_path = SAMPLES / "run_log.json"
    with open(log_path, "w") as f:
        json.dump({
            "args": vars(args),
            "total_elapsed_s": round(total_elapsed, 3),
            "batches": all_stats,
        }, f, indent=2)
    print(f"\nTotal: {sum(s['n_valid'] for s in all_stats)} amostras válidas "
          f"em {total_elapsed:.1f}s")
    print(f"Log:  {log_path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
