#!/usr/bin/env python3
"""
engine_6x6_sampler.py
=====================
Engine de amostragem de caminhos hamiltonianos no tabuleiro 6×6 do cavalo.
Réplica do engine_8x8.py com BOARD=6, usada para a Fase A de validação:
queremos confirmar que a amostragem Z3 reproduz os invariantes que o
catalogador exaustivo (cavalo_loop_destruicao_6x6.py) já estabeleceu.

Diferenças em relação ao engine_8x8.py:
  - BOARD = 6
  - break_symmetry default = False (queremos comparar com GT sem viés)
  - notação chess: ranks 1..6 (em vez de 1..8)
"""

from time import time
from z3 import Solver, Bool, Sum, If, Or, sat, is_true

BOARD = 6
KNIGHT_DELTAS = [(1,2),(2,1),(-1,2),(-2,1),(1,-2),(2,-1),(-1,-2),(-2,-1)]


# ── construção do grafo ───────────────────────────────────────────────

def knight_nbrs(x, y):
    out = []
    for dx, dy in KNIGHT_DELTAS:
        nx, ny = x + dx, y + dy
        if 0 <= nx < BOARD and 0 <= ny < BOARD:
            out.append((nx, ny))
    def deg(n):
        return sum(1 for dx, dy in KNIGHT_DELTAS
                   if 0 <= n[0]+dx < BOARD and 0 <= n[1]+dy < BOARD)
    return sorted(out, key=deg)

def build_graph():
    nodes = [(x, y) for x in range(BOARD) for y in range(BOARD)]
    adj = {n: knight_nbrs(*n) for n in nodes}
    edges = set()
    for n, nbrs in adj.items():
        for m in nbrs:
            edges.add(tuple(sorted([n, m])))
    return nodes, adj, edges

def get_ciclos_de_4(nodes, adj):
    ciclos = set()
    for u in nodes:
        for v in adj[u]:
            for w in adj[v]:
                if w != u:
                    for x in adj[w]:
                        if x != v and u in adj[x]:
                            ciclos.add(tuple(sorted([u, v, w, x])))
    return ciclos

def spanning_tree(adj, root):
    parent = {root: None}
    stack = [root]
    visited = {root}
    te = set()
    while stack:
        u = stack.pop()
        for v in reversed(adj[u]):
            if v not in visited:
                visited.add(v)
                parent[v] = u
                te.add(tuple(sorted([u, v])))
                stack.append(v)
    return parent, te

def extract_loops(parent, te, all_edges):
    loops = []
    for e in sorted(all_edges):
        if e in te:
            continue
        u, v = e
        pu = []
        c = u
        while c is not None:
            pu.append(c)
            c = parent[c]
        pv = []
        c = v
        while c is not None:
            pv.append(c)
            c = parent[c]
        loops.append({
            "c1": [list(p) for p in pu],
            "c2": [list(p) for p in pv],
            "col": [list(u), list(v)],
        })
    return loops

def get_loop_edges(lp):
    es = set()
    def tog(u, v):
        e = tuple(sorted([tuple(u), tuple(v)]))
        if e in es: es.remove(e)
        else: es.add(e)
    for i in range(len(lp["c1"]) - 1): tog(lp["c1"][i], lp["c1"][i+1])
    for i in range(len(lp["c2"]) - 1): tog(lp["c2"][i], lp["c2"][i+1])
    tog(lp["col"][0], lp["col"][1])
    return es


# ── notação de chess para arestas ─────────────────────────────────────

def node_label(n):
    """(row, col) -> 'A6' ... 'F1' (rank 6 = top)."""
    r, c = n
    return chr(ord('A') + c) + str(BOARD - r)

def edge_label(e):
    u, v = e
    return f"{node_label(u)}-{node_label(v)}"


# ── API de amostragem ─────────────────────────────────────────────────

def sample_signatures(start, end, n_target,
                      max_attempts=None,
                      timeout_s=None,
                      break_symmetry=False,
                      verbose=False,
                      log=print):
    """
    Amostra caminhos hamiltonianos de `start` a `end` no grafo do cavalo 6×6.

    Padrão `break_symmetry=False` (para validação): queremos comparar com
    o ground truth exaustivo sem nenhum viés adicional do solver.
    Para o run de A4 (quantificação do viés), passar break_symmetry=True.

    Demais parâmetros e semântica: idênticos ao engine_8x8.sample_signatures.
    """
    t0 = time()
    if max_attempts is None:
        max_attempts = 5 * n_target

    nodes, adj, edges_orig = build_graph()
    virtual = tuple(sorted([start, end]))
    edges_aug = set(edges_orig) | {virtual}

    adj_aug = {n: list(adj[n]) for n in nodes}
    adj_aug[start].append(end)
    adj_aug[end].append(start)

    par, te = spanning_tree(adj_aug, nodes[0])
    loops = extract_loops(par, te, edges_aug)
    if verbose:
        log(f"  H1_aug = {len(loops)}  (loops fundamentais no grafo aumentado)")

    solver = Solver()
    lvars = [Bool(f"L{i}") for i in range(len(loops))]
    e2l = {e: [] for e in edges_aug}
    for i, lp in enumerate(loops):
        for e in get_loop_edges(lp):
            e2l[e].append(i)

    def ea(e):
        idxs = e2l.get(e, [])
        if not idxs:
            return False
        return (Sum([If(lvars[i], 1, 0) for i in idxs]) % 2 == 1)

    for n in nodes:
        inc = [e for e in edges_aug if n in e]
        solver.add(Sum([If(ea(e), 1, 0) for e in inc]) == 2)

    solver.add(ea(virtual) == True)

    sym_applied = False
    if break_symmetry and start == (0, 0):
        edge_sym = tuple(sorted([(0, 0), (1, 2)]))
        if edge_sym in edges_aug:
            solver.add(ea(edge_sym) == True)
            sym_applied = True

    c4 = get_ciclos_de_4(nodes, adj)
    for ciclo in c4:
        es4 = [
            tuple(sorted([ciclo[0], ciclo[1]])),
            tuple(sorted([ciclo[1], ciclo[2]])),
            tuple(sorted([ciclo[2], ciclo[3]])),
            tuple(sorted([ciclo[3], ciclo[0]])),
        ]
        if all(e in edges_aug for e in es4):
            solver.add(Sum([If(ea(e), 1, 0) for e in es4]) <= 3)

    edges_s = sorted(edges_orig)
    sigs = []
    attempts = 0
    exhausted = False
    timed_out = False

    while len(sigs) < n_target:
        if attempts >= max_attempts:
            break
        if timeout_s is not None and (time() - t0) > timeout_s:
            timed_out = True
            break

        attempts += 1
        if solver.check() != sat:
            exhausted = True
            break

        m = solver.model()
        active = [e for e in edges_orig if is_true(m.evaluate(ea(e)))]
        adj_r = {n: [] for n in nodes}
        for u, v in active:
            adj_r[u].append(v)
            adj_r[v].append(u)

        vis = {start}
        q = [start]
        while q:
            c = q.pop(0)
            for nb in adj_r[c]:
                if nb not in vis:
                    vis.add(nb)
                    q.append(nb)

        if len(vis) == len(nodes):
            sigs.append([is_true(m.evaluate(ea(e))) for e in edges_s])
            if verbose and len(sigs) % 100 == 0:
                log(f"  -> {len(sigs)}/{n_target}  "
                    f"(attempts={attempts}, elapsed={time()-t0:.1f}s)")

        solver.add(Or([lvars[i] != is_true(m.evaluate(lvars[i]))
                       for i in range(len(lvars))]))

        if len(vis) < len(nodes):
            seen = set(vis)
            for seed in [n for n in nodes if n not in seen]:
                if seed in seen:
                    continue
                comp = {seed}
                q2 = [seed]
                while q2:
                    c = q2.pop(0)
                    for nb in adj_r[c]:
                        if nb not in comp:
                            comp.add(nb)
                            q2.append(nb)
                in_e = [e for e in edges_aug if e[0] in comp and e[1] in comp]
                solver.add(Sum([If(ea(e), 1, 0) for e in in_e]) <= len(comp) - 1)
                seen |= comp

    elapsed = time() - t0
    edge_strs = [edge_label(e) for e in edges_s]

    return {
        "start": list(start),
        "end": list(end),
        "n_target": n_target,
        "n_found": len(sigs),
        "n_attempts": attempts,
        "exhausted": exhausted,
        "timed_out": timed_out,
        "elapsed_s": round(elapsed, 2),
        "break_symmetry_applied": sym_applied,
        "edges": edge_strs,
        "signatures": sigs,
    }


# ── smoke test ────────────────────────────────────────────────────────

if __name__ == "__main__":
    import json
    print(f"Grafo {BOARD}×{BOARD}: V=36, E=80, H1=45")
    print(f"Smoke test: amostrando 20 caminhos de A6 (=(0,0)) → F6 (=(0,5))...\n")
    out = sample_signatures(
        start=(0, 0), end=(0, 5),
        n_target=20,
        verbose=True,
        break_symmetry=False,
    )
    summary = {k: v for k, v in out.items() if k != "signatures"}
    summary["edges"] = f"<{len(out['edges'])} edge labels>"
    print(json.dumps(summary, indent=2))
    if out["signatures"]:
        print(f"\nPrimeira assinatura tem {sum(out['signatures'][0])} arestas ativas "
              f"(esperado 35 = caminho hamiltoniano de 36 vértices)")
