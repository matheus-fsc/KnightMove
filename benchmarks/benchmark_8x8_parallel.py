"""benchmark_8x8_parallel.py — Benchmark paralelo no tabuleiro 8×8.

Roda TODOS os métodos simultaneamente, cada um fixado em 1 núcleo de CPU
(os.sched_setaffinity), pelo MESMO tempo de parede T. Cada método tenta achar
o máximo de tours (ciclos hamiltonianos fechados) possível em T segundos;
a métrica de comparação é tours/s.

Métodos (todos genéricos para n=8):
  1. flagship    — knight_tours.py (cycle-space + propagação + UF incremental)
  2. minimal_v2  — incremental_subtour/scaling_minimal_v2.MinimalV2 (pressão+UF)
  3. theory      — local_phase_heuristic theory (f∞ por nível)
  4. uniforme    — theory com f=0.25 ∀L
  5. naive       — backtracking ingênuo (sem otimização alguma)
  6. warnsdorff  — backtracking clássico com ordenação de Warnsdorff
  7. z3_puro     — Z3 SAT (grau-2) + corte de sub-ciclos
  8. z3_mandatory— Z3 + arestas obrigatórias dos cantos

Uso:
    .venv/bin/python benchmark_8x8_parallel.py [T_segundos]

Default T = 300s. Cada worker escreve seu resultado em
data/parallel_8x8/<metodo>.json ao terminar. Ao final o processo principal
agrega e imprime/salva o resumo em data/parallel_8x8/summary.json.
"""
from __future__ import annotations

# --- bootstrap: torna core/ importavel a partir de qualquer subpasta ---
import sys as _sys, pathlib as _pathlib
for _anc in _pathlib.Path(__file__).resolve().parents:
    if (_anc / "core").is_dir() and (_anc / "pyproject.toml").exists():
        _sys.path.insert(0, str(_anc / "core"))
        break
# --- fim do bootstrap ---

import json
import multiprocessing as mp
import os
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUTDIR = ROOT / "data" / "parallel_8x8"

N = 8
SEED = 0


# ---------------------------------------------------------------------------
# infra de núcleo / escrita
# ---------------------------------------------------------------------------

def _pin(core: int):
    try:
        os.sched_setaffinity(0, {core})
    except Exception:
        pass


def _write(method: str, payload: dict):
    OUTDIR.mkdir(parents=True, exist_ok=True)
    tmp = OUTDIR / f"{method}.json.tmp"
    final = OUTDIR / f"{method}.json"
    tmp.write_text(json.dumps(payload, indent=2, default=str))
    tmp.replace(final)


# ---------------------------------------------------------------------------
# métodos — cada um devolve dict com n_tours, n_nos (ou None), tempo_s
# ---------------------------------------------------------------------------

def run_flagship(timeout: float) -> dict:
    import numpy as np
    from knight_tours import (
        build_graph, State, _propagate_initial, _choose_next_edge,
        fix_and_propagate, _is_complete_tour,
        OK, COMPLETE_TOUR, CONTRADICTION, SUBTOUR,
    )
    ctx = build_graph(N)
    rng = np.random.default_rng(SEED)
    state = State(ctx['V'], ctx['E'])
    status = _propagate_initial(state, ctx)
    cnt = {'nodes': 0, 'tours': 0}
    deadline = time.perf_counter() + timeout
    stop = {'v': False}

    def bt():
        cnt['nodes'] += 1
        if (cnt['nodes'] & 0x3FFF) == 0 and time.perf_counter() > deadline:
            stop['v'] = True
            return True
        e, first_val = _choose_next_edge(state, ctx, rng)
        if e == -1:
            if _is_complete_tour(state, ctx):
                cnt['tours'] += 1
            return False
        for val in (first_val, 1 - first_val):
            snap = state.snapshot()
            st = fix_and_propagate(state, ctx, e, val)
            if st == OK:
                if bt():
                    state.restore(snap)
                    return True
            elif st == COMPLETE_TOUR:
                if _is_complete_tour(state, ctx):
                    cnt['tours'] += 1
            state.restore(snap)
            if stop['v']:
                return True
        return False

    t0 = time.perf_counter()
    if status not in (CONTRADICTION, SUBTOUR):
        if status == COMPLETE_TOUR:
            if _is_complete_tour(state, ctx):
                cnt['tours'] += 1
        else:
            bt()
    dt = time.perf_counter() - t0
    return {"tempo_s": dt, "n_tours": cnt['tours'], "n_nos": cnt['nodes'],
            "nos_por_tour": cnt['nodes'] / max(1, cnt['tours'])}


def run_minimal_v2(timeout: float) -> dict:
    sys.path.insert(0, str(ROOT / "incremental_subtour"))
    from scaling_minimal_v2 import MinimalV2
    bt = MinimalV2(n=N, alvo=10**9, timeout=timeout)
    r = bt.executar()
    return {"tempo_s": r["tempo_s"], "n_tours": r["n_tours"],
            "n_nos": r["n_nos"], "nos_por_tour": r["nos_por_tour"]}


def run_theory(timeout: float) -> dict:
    sys.path.insert(0, str(ROOT / "incremental_subtour"))
    sys.path.insert(0, str(ROOT / "local_phase_heuristic"))
    from theory_heuristic import TheoryHeuristicV2
    bt = TheoryHeuristicV2(n=N, alvo=10**9, timeout=timeout)
    r = bt.executar()
    return {"tempo_s": r["tempo_s"], "n_tours": r["n_tours"],
            "n_nos": r["n_nos"], "nos_por_tour": r["nos_por_tour"]}


def run_uniforme(timeout: float) -> dict:
    sys.path.insert(0, str(ROOT / "incremental_subtour"))
    sys.path.insert(0, str(ROOT / "local_phase_heuristic"))
    from theory_heuristic import TheoryHeuristicV2
    f_unif = {L: 0.25 for L in range(N)}
    bt = TheoryHeuristicV2(n=N, alvo=10**9, timeout=timeout, f_inf=f_unif)
    r = bt.executar()
    return {"tempo_s": r["tempo_s"], "n_tours": r["n_tours"],
            "n_nos": r["n_nos"], "nos_por_tour": r["nos_por_tour"]}


def run_naive(timeout: float) -> dict:
    from benchmark_tabela9_completo import NaiveBacktracker
    bt = NaiveBacktracker(n=N, alvo=10**9, timeout=timeout)
    r = bt.executar()
    return {"tempo_s": r["tempo_s"], "n_tours": r["n_tours"],
            "n_nos": r["n_nos"], "nos_por_tour": r["nos_por_tour"]}


def run_warnsdorff(timeout: float) -> dict:
    from benchmark_8x8 import warnsdorff_backtracking
    tours, elapsed = warnsdorff_backtracking(N, K=10**9, seed=SEED,
                                             timeout=timeout)
    return {"tempo_s": elapsed, "n_tours": len(tours),
            "n_nos": None, "nos_por_tour": None}


def _run_z3(timeout: float, com_mandatory: bool) -> dict:
    sys.path.insert(0, str(ROOT / "incremental_subtour"))
    from scaling_minimal_v2 import build_knight_graph, mandatory_corner_edges
    from z3 import Bool, Not, Or, PbEq, Solver, Sum, sat

    edges_uv, V, v2e = build_knight_graph(N)
    E = len(edges_uv)
    x = [Bool(f"x_{i}") for i in range(E)]
    s = Solver()
    for vtx, es in v2e.items():
        s.add(PbEq([(x[e], 1) for e in es], 2))
    if com_mandatory:
        for e in mandatory_corner_edges(N, v2e):
            s.add(x[e])

    t0 = time.perf_counter()
    tours = []
    tentativas = 0
    while True:
        if time.perf_counter() - t0 > timeout:
            break
        tentativas += 1
        if s.check() != sat:
            break
        m = s.model()
        ativos = [e for e in range(E)
                  if m[x[e]] is not None and bool(m[x[e]])]
        adj = [[] for _ in range(V)]
        for e in ativos:
            u, v = edges_uv[e]
            adj[u].append(v)
            adj[v].append(u)
        bfs = [0]; seen = {0}; qi = 0
        while qi < len(bfs):
            node = bfs[qi]; qi += 1
            for nb in adj[node]:
                if nb not in seen:
                    seen.add(nb)
                    bfs.append(nb)
        if len(seen) == V:
            tours.append(ativos)
            s.add(Sum([x[e] if e in ativos else Not(x[e])
                       for e in range(E)]) < E)
        else:
            ac = [e for e in range(E)
                  if edges_uv[e][0] in seen and edges_uv[e][1] in seen
                  and m[x[e]] is not None and bool(m[x[e]])]
            if ac:
                s.add(Or([Not(x[e]) for e in ac]))
    dt = time.perf_counter() - t0
    return {"tempo_s": dt, "n_tours": len(tours), "tentativas": tentativas,
            "n_nos": None, "nos_por_tour": None}


def run_z3_puro(timeout: float) -> dict:
    return _run_z3(timeout, com_mandatory=False)


def run_z3_mandatory(timeout: float) -> dict:
    return _run_z3(timeout, com_mandatory=True)


# ---------------------------------------------------------------------------
# orquestração
# ---------------------------------------------------------------------------

METHODS = [
    ("flagship",     "knight_tours flagship (cycle-space+UF)", run_flagship),
    ("minimal_v2",   "incremental_subtour MinimalV2 (pressão+UF)", run_minimal_v2),
    ("theory",       "theory f∞ por nível",                  run_theory),
    ("uniforme",     "theory uniforme f=0.25",               run_uniforme),
    ("naive",        "backtracking ingênuo (sem otimização)", run_naive),
    ("warnsdorff",   "Warnsdorff BT clássico",               run_warnsdorff),
    ("z3_puro",      "Z3 puro (grau-2 + corte sub-ciclo)",   run_z3_puro),
    ("z3_mandatory", "Z3 + cantos obrigatórios",             run_z3_mandatory),
]


def _worker(method: str, label: str, fn, core: int, timeout: float):
    _pin(core)
    _write(method, {"method": method, "label": label, "core": core,
                    "status": "running", "timeout_s": timeout,
                    "started_at": time.time()})
    t0 = time.perf_counter()
    try:
        r = fn(timeout)
        r.update({"method": method, "label": label, "core": core,
                  "status": "done", "timeout_s": timeout,
                  "wall_s": time.perf_counter() - t0})
        if r.get("tempo_s"):
            r["tours_por_s"] = r["n_tours"] / r["tempo_s"]
        _write(method, r)
    except Exception as exc:
        import traceback
        _write(method, {"method": method, "label": label, "core": core,
                        "status": "error", "error": str(exc),
                        "trace": traceback.format_exc()})


def main():
    timeout = float(sys.argv[1]) if len(sys.argv) > 1 else 300.0
    OUTDIR.mkdir(parents=True, exist_ok=True)

    ncpu = os.cpu_count() or 8
    print(f">>> BENCHMARK 8×8 PARALELO — {len(METHODS)} métodos, "
          f"{timeout:.0f}s cada, 1 núcleo/método (de {ncpu} CPUs) <<<")
    for i, (m, label, _) in enumerate(METHODS):
        print(f"  core {i:2d}: {m:<13} — {label}")
    print(f"\n  saída: {OUTDIR}/<metodo>.json")
    print(f"  ETA  : ~{timeout:.0f}s\n", flush=True)

    ctx = mp.get_context("fork")
    procs = []
    for i, (method, label, fn) in enumerate(METHODS):
        core = i % ncpu
        p = ctx.Process(target=_worker, args=(method, label, fn, core, timeout))
        p.start()
        procs.append((method, p))

    t_start = time.perf_counter()
    for method, p in procs:
        p.join()
    wall = time.perf_counter() - t_start

    # agregação
    print("\n" + "=" * 92)
    print(f"RESUMO — benchmark 8×8 paralelo  (T={timeout:.0f}s, wall={wall:.1f}s)")
    print("=" * 92)
    print(f"{'Método':<34}{'t(s)':>8}{'tours':>10}{'tours/s':>11}"
          f"{'speedup':>10}{'nós/tour':>10}")
    print("-" * 92)

    results = {}
    for method, _ in procs:
        f = OUTDIR / f"{method}.json"
        results[method] = json.loads(f.read_text()) if f.exists() else {}

    # baseline para speedup = naive (se válido), senão menor tours/s>0
    def rate(r):
        t = r.get("tempo_s") or 0
        return (r.get("n_tours", 0) / t) if t > 0 else 0.0

    base = rate(results.get("naive", {})) or 1e-9

    rows = []
    for method, label, _ in METHODS:
        r = results.get(method, {})
        if r.get("status") == "error":
            print(f"{label:<34}{'ERRO':>8}  {r.get('error','')[:40]}")
            rows.append((method, label, None))
            continue
        t = r.get("tempo_s") or 0
        tours = r.get("n_tours", 0)
        rt = rate(r)
        sp = rt / base if base > 0 else 0
        npt = r.get("nos_por_tour")
        npt_s = f"{npt:.1f}" if isinstance(npt, (int, float)) and npt else "-"
        print(f"{label:<34}{t:>8.1f}{tours:>10}{rt:>11.1f}{sp:>9.1f}×{npt_s:>10}")
        rows.append((method, label, rt))

    print("-" * 92)
    print("  speedup relativo ao backtracking ingênuo (naive).")

    summary = {
        "n": N, "timeout_s": timeout, "wall_s": wall,
        "results": results,
        "ranking_tours_por_s": sorted(
            [(m, rate(results.get(m, {}))) for m, _, _ in METHODS],
            key=lambda x: -x[1]),
    }
    (OUTDIR / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print(f"\n  → resumo salvo em {OUTDIR / 'summary.json'}")


if __name__ == "__main__":
    main()
