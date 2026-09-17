#!/usr/bin/env python3
"""
knight_8x8_allsat_async.py
==========================================================================
AVISO DE ESCALA: O tabuleiro 8x8 possui H1=105 loops fundamentais (arestas
não-árvore no grafo do cavalo). O número total de ciclos hamiltonianos
fechados é da ordem de 10^13, tornando a enumeração completa completamente
inviável em qualquer hardware atual (anos de CPU-time).

Este script opera em modo AMOSTRAL: coleta MAX_SOLUTIONS (padrão: 50.000)
representantes canônicos — suficiente para verificar estatisticamente se o
grupo diedral D4 age livremente sobre o espaço de soluções.

Objetivo: verificar se total_soluções % 8 == 0, confirmando (ou refutando)
que cada órbita D4 tem exatamente 8 elementos, i.e., nenhuma solução possui
simetria interna.

Hardware alvo: RTX 4060 8GB VRAM | 16GB RAM | i7-13ª geração
Paralelismo: 6 workers em ProcessPoolExecutor (2 threads livres para o SO)
==========================================================================
"""

import asyncio
import json
import os
import sys
import time
from collections import deque
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp

# ==========================================================================
# CONFIGURAÇÃO
# ==========================================================================
BOARD = 8
START_NODE = (0, 0)
MAX_SOLUTIONS = int(os.getenv("MAX_SOLUTIONS", "50000"))
MAX_WORKERS = 6
NUM_PARTITIONS = 8          # 2^3: fixa L_0, L_1, L_2 para criar 8 regiões disjuntas
CHECKPOINT_FILE = "checkpoint_8x8.json"
OUTPUT_FILE = "assinaturas_8x8.json"
CHECKPOINT_INTERVAL = 1000  # salva checkpoint a cada N soluções globais
DEDUP_REPORT_INTERVAL = 500 # loga estatísticas de dedup a cada N soluções
BATCH_SIZE = 50             # soluções por batch enviadas ao queue

# ==========================================================================
# GRUPO DIEDRAL D4 — 7 transformações não-identidade para n×n
# ==========================================================================
def rot90(x, y, n=BOARD):  return (y, n - 1 - x)
def rot180(x, y, n=BOARD): return (n - 1 - x, n - 1 - y)
def rot270(x, y, n=BOARD): return (n - 1 - y, x)
def refX(x, y, n=BOARD):   return (n - 1 - x, y)
def refY(x, y, n=BOARD):   return (x, n - 1 - y)
def refD1(x, y, n=BOARD):  return (y, x)
def refD2(x, y, n=BOARD):  return (n - 1 - y, n - 1 - x)

TRANSFORMACOES = [rot90, rot180, rot270, refX, refY, refD1, refD2]


def transpor_aresta(u, v, transformacao, n=BOARD):
    """Aplica rotação/espelhamento D4 aos dois vértices de uma aresta."""
    nu = transformacao(*u, n)
    nv = transformacao(*v, n)
    return tuple(sorted([nu, nv]))


# ==========================================================================
# GRAFO DO CAVALO 8×8
# ==========================================================================
KNIGHT_DELTAS = [(1, 2), (2, 1), (-1, 2), (-2, 1),
                 (1, -2), (2, -1), (-1, -2), (-2, -1)]


def neighbors(x, y, n=BOARD):
    out = []
    for dx, dy in KNIGHT_DELTAS:
        nx, ny = x + dx, y + dy
        if 0 <= nx < n and 0 <= ny < n:
            out.append((nx, ny))
    def degree(p):
        return sum(1 for dx, dy in KNIGHT_DELTAS
                   if 0 <= p[0] + dx < n and 0 <= p[1] + dy < n)
    return sorted(out, key=degree)


def build_graph(n=BOARD):
    nodes = [(x, y) for x in range(n) for y in range(n)]
    adj = {node: neighbors(*node, n) for node in nodes}
    edges = set()
    for node, nbrs in adj.items():
        for m in nbrs:
            edges.add(tuple(sorted([node, m])))
    return nodes, adj, edges


# ==========================================================================
# ÁRVORE GERADORA DFS E LOOPS FUNDAMENTAIS (CICLOS DE HOMOLOGIA H1)
# ==========================================================================
def compute_spanning_tree(adj, start=START_NODE):
    parent = {start: None}
    level = {start: 0}
    order = []
    stack = [start]
    tree_edges = set()
    visited = {start}
    while stack:
        u = stack.pop()
        order.append(u)
        for v in reversed(adj[u]):
            if v not in visited:
                visited.add(v)
                parent[v] = u
                level[v] = level[u] + 1
                tree_edges.add(tuple(sorted([u, v])))
                stack.append(v)
    return parent, level, order, tree_edges


def extract_loops(parent, tree_edges, edges):
    """Cada aresta não-árvore define um loop fundamental na base H1."""
    loops = []
    for e in sorted(edges):
        if e in tree_edges:
            continue
        u, v = e
        path_u = []
        curr = u
        while curr is not None:
            path_u.append(curr)
            curr = parent[curr]
        path_v = []
        curr = v
        while curr is not None:
            path_v.append(curr)
            curr = parent[curr]
        loops.append({
            "caminho1": [list(p) for p in path_u],
            "caminho2": [list(p) for p in path_v],
            "colisao": [list(u), list(v)],
            "tamanho_ciclo": len(path_u) + len(path_v),
        })
    return loops


# ==========================================================================
# WORKER — executa em processo separado via ProcessPoolExecutor
# ==========================================================================
def worker_search(partition_id, partition_bits, result_queue, stop_event,
                  max_solutions_per_worker):
    """
    All-SAT Z3 dentro de uma partição disjunta do espaço de busca.

    A partição é definida fixando os primeiros 3 loop vars (L_0, L_1, L_2)
    aos bits do partition_id, criando 8 regiões mutuamente exclusivas.
    Soluções são enviadas ao result_queue em batches de BATCH_SIZE.
    """
    from z3 import Solver, Bool, Sum, If, Not, Or, sat, is_true

    nodes, adj, edges = build_graph()
    parent, _level, _order, tree_edges = compute_spanning_tree(adj)
    loops = extract_loops(parent, tree_edges, edges)

    solver = Solver()
    loop_vars = [Bool(f"L_{i}") for i in range(len(loops))]

    # --- Mapa aresta → loops que a contêm (via XOR em GF(2)) ---
    def get_loop_edge_set(loop):
        edges_set = set()
        def toggle(u, v):
            e = tuple(sorted([tuple(u), tuple(v)]))
            if e in edges_set:
                edges_set.remove(e)
            else:
                edges_set.add(e)
        c1, c2 = loop['caminho1'], loop['caminho2']
        for i in range(len(c1) - 1):
            toggle(c1[i], c1[i + 1])
        for i in range(len(c2) - 1):
            toggle(c2[i], c2[i + 1])
        toggle(loop['colisao'][0], loop['colisao'][1])
        return edges_set

    edge_to_loops = {e: [] for e in edges}
    for i, lp in enumerate(loops):
        for e in get_loop_edge_set(lp):
            edge_to_loops[e].append(i)

    # --- Ativação de aresta via GF(2): aresta ativa ⟺ soma ímpar dos loops ---
    edge_active = {}
    for e in edges:
        if not edge_to_loops[e]:
            edge_active[e] = False
        else:
            edge_active[e] = (
                Sum([If(loop_vars[i], 1, 0) for i in edge_to_loops[e]]) % 2 == 1
            )

    # --- Restrição de grau = 2 em todos os vértices ---
    for n in nodes:
        incident = [e for e in edges if n in e]
        solver.add(Sum([If(edge_active[e], 1, 0) for e in incident]) == 2)

    # --- Symmetry Breaking Predicates (SBP) do grupo D4 ---
    arestas_ordenadas = sorted(list(edges))

    def lex_leq(vec_a, vec_b):
        """vec_a <=_lex vec_b como vetores booleanos Z3."""
        cond = True
        for a, b in reversed(list(zip(vec_a, vec_b))):
            cond = If(a == b, cond, b)
        return cond

    vetor_original = [edge_active[e] for e in arestas_ordenadas]
    for transformacao in TRANSFORMACOES:
        vetor_transformado = [
            edge_active[transpor_aresta(e[0], e[1], transformacao)]
            for e in arestas_ordenadas
        ]
        solver.add(lex_leq(vetor_original, vetor_transformado))

    # --- Restrições de partição: fixa L_0, L_1, L_2 ---
    for idx, val in enumerate(partition_bits):
        solver.add(loop_vars[idx] if val else Not(loop_vars[idx]))

    found = 0
    batch = []

    while not stop_event.is_set() and found < max_solutions_per_worker:
        from z3 import sat
        if solver.check() != sat:
            break

        m = solver.model()
        active_edges = [
            e for e in edges
            if edge_to_loops[e] and is_true(m.evaluate(edge_active[e]))
        ]

        # --- Verificação de conectividade hamiltoniana (BFS) ---
        adj_res = {nd: [] for nd in nodes}
        for u, v in active_edges:
            adj_res[u].append(v)
            adj_res[v].append(u)

        visited = set()
        components = []
        for node in nodes:
            if node not in visited:
                comp, q = [], [node]
                visited.add(node)
                while q:
                    curr = q.pop(0)
                    comp.append(curr)
                    for nbr in adj_res[curr]:
                        if nbr not in visited:
                            visited.add(nbr)
                            q.append(nbr)
                components.append(comp)

        assinatura = tuple(is_true(m.evaluate(v)) for v in loop_vars)

        if len(components) == 1:
            found += 1
            batch.append(assinatura)

            if len(batch) >= BATCH_SIZE:
                result_queue.put(("batch", partition_id, batch))
                batch = []

            # Cláusula de bloqueio: impede Z3 de reutilizar esta assinatura
            clausula_bloqueio = Or([
                Not(loop_vars[i]) if assinatura[i] else loop_vars[i]
                for i in range(len(loops))
            ])
            solver.add(clausula_bloqueio)
        else:
            # Eliminação de sub-tours: impõe espanamento de cada componente
            for comp in components:
                comp_set = set(comp)
                in_edges = [e for e in edges
                            if e[0] in comp_set and e[1] in comp_set]
                solver.add(
                    Sum([If(edge_active[e], 1, 0) for e in in_edges])
                    <= len(comp) - 1
                )

    # Flush do batch final
    if batch:
        result_queue.put(("batch", partition_id, batch))

    result_queue.put(("done", partition_id, found))
    return found


# ==========================================================================
# CHECKPOINT
# ==========================================================================
def load_checkpoint():
    if os.path.exists(CHECKPOINT_FILE):
        with open(CHECKPOINT_FILE) as f:
            data = json.load(f)
        sigs = set(tuple(s) for s in data["assinaturas"])
        total = data["total"]
        print(f"[Checkpoint] Retomando de {data['timestamp']} — {total} soluções já coletadas.")
        return sigs, total
    return set(), 0


def save_checkpoint(assinaturas, total):
    data = {
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "total": total,
        "assinaturas": [list(sig) for sig in assinaturas],
    }
    with open(CHECKPOINT_FILE, "w") as f:
        json.dump(data, f)
    print(f"  [Checkpoint] Salvo: {total} soluções → {CHECKPOINT_FILE}")


# ==========================================================================
# ANÁLISE FINAL — Lema de Burnside
# ==========================================================================
def analyze_results(total):
    print()
    print("=" * 62)
    print("ANÁLISE FINAL — AÇÃO DO GRUPO D4")
    print("=" * 62)
    print(f"  Total de representantes canônicos: {total}")
    print(f"  Total % 8 = {total % 8}")
    print()

    if total % 8 == 0:
        orbitas = total // 8
        print(f"  Número de órbitas D4: {orbitas}")
        print()
        print("  D4 age livremente: nenhuma solução tem simetria interna")
    else:
        r = total % 8
        print(f"  ATENÇÃO: {r} soluções não se encaixam em órbitas completas de tamanho 8.")
        print()
        print("  Lema de Burnside — análise de órbitas de tamanho < 8:")
        print("  (Tamanhos possíveis: 1, 2, 4 — divisores próprios de 8)")
        print()
        # Órbita de tamanho 4: fixadas por rot180 (única rotação de ordem 2)
        # Órbita de tamanho 2: fixadas por rot180 E por uma reflexão
        # Órbita de tamanho 1: fixadas por todos os 8 elementos de D4
        if r % 4 == 0:
            print(f"  Possível: {r // 4} órbitas de tamanho 4 (fixadas por rot180)")
        if r % 2 == 0:
            print(f"  Possível: {r // 2} órbitas de tamanho 2 (fixadas por rot180 + reflexão)")
        print()
        print("  Nota: Para determinar o tipo exato, execute com MAX_SOLUTIONS maior")
        print("  e verifique quais soluções são invariantes sob cada gerador de D4.")

    print("=" * 62)
    print()


# ==========================================================================
# HELPER — leitura não-bloqueante do queue (roda em thread pool)
# ==========================================================================
def _queue_get_timeout(q, timeout):
    try:
        return q.get(timeout=timeout)
    except Exception:
        return None


# ==========================================================================
# COORDINATOR — coroutine asyncio que agrega resultados dos workers
# ==========================================================================
async def coordinator(result_queue, all_signatures, initial_total):
    """
    Drena o result_queue assincronamente.
    - Deduplica entre workers a cada DEDUP_REPORT_INTERVAL soluções.
    - Salva checkpoint a cada CHECKPOINT_INTERVAL soluções.
    - Imprime ETA baseado na taxa das últimas 500 soluções encontradas.
    """
    loop = asyncio.get_running_loop()
    total = initial_total
    duplicatas_rejeitadas = 0
    workers_done = 0
    last_checkpoint = total
    last_dedup_report = total
    recent_times = deque(maxlen=500)  # timestamps das últimas 500 soluções novas

    print(f"[Coordinator] Iniciado. Coletando até {MAX_SOLUTIONS} soluções...\n")

    while workers_done < NUM_PARTITIONS:
        msg = await loop.run_in_executor(None, _queue_get_timeout, result_queue, 0.15)

        if msg is None:
            await asyncio.sleep(0.05)
            continue

        kind = msg[0]

        if kind == "batch":
            _, partition_id, batch = msg
            novas = 0
            for sig in batch:
                if sig not in all_signatures:
                    all_signatures.add(sig)
                    total += 1
                    novas += 1
                    recent_times.append(time.monotonic())
                else:
                    duplicatas_rejeitadas += 1

            if novas:
                print(f"  [P{partition_id}] +{novas:3d} novas | Total: {total}/{MAX_SOLUTIONS}")

            # Relatório de dedup entre workers
            if total - last_dedup_report >= DEDUP_REPORT_INTERVAL:
                print(f"  [Dedup] {duplicatas_rejeitadas} duplicatas rejeitadas até agora.")
                last_dedup_report = total

            # Checkpoint + ETA
            if total - last_checkpoint >= CHECKPOINT_INTERVAL:
                save_checkpoint(all_signatures, total)
                last_checkpoint = total

                if len(recent_times) >= 2:
                    elapsed = recent_times[-1] - recent_times[0] + 1e-9
                    rate = (len(recent_times) - 1) / elapsed  # sol/s
                    restante_s = (MAX_SOLUTIONS - total) / max(rate, 1e-9)
                    print(f"  [ETA] {restante_s / 60:.1f} min restantes "
                          f"| Taxa: {rate * 60:.0f} sol/min")

            if total >= MAX_SOLUTIONS:
                print(f"\n[Coordinator] Limite MAX_SOLUTIONS={MAX_SOLUTIONS} atingido.")
                break

        elif kind == "done":
            _, partition_id, count = msg
            workers_done += 1
            print(f"  [P{partition_id}] Partição concluída com {count} soluções locais. "
                  f"({workers_done}/{NUM_PARTITIONS} partições finalizadas)")

    return total, all_signatures


# ==========================================================================
# MAIN ASYNC
# ==========================================================================
async def async_main():
    print()
    print("=" * 62)
    print(" Knight Tour 8×8 — All-SAT Assíncrono com D4 SBP")
    print(f" MAX_SOLUTIONS={MAX_SOLUTIONS} | WORKERS={MAX_WORKERS} | PARTIÇÕES={NUM_PARTITIONS}")
    print("=" * 62)
    print()

    # Carrega checkpoint existente
    all_signatures, initial_total = load_checkpoint()

    if initial_total >= MAX_SOLUTIONS:
        print(f"Checkpoint já atingiu o limite ({initial_total} soluções). Pulando busca.")
        analyze_results(initial_total)
        _export(all_signatures, initial_total)
        return

    # 8 partições: fixa L_0, L_1, L_2 a cada combinação binária de 3 bits
    partitions = [
        [(i >> bit) & 1 == 1 for bit in range(3)]
        for i in range(NUM_PARTITIONS)
    ]

    # Cada worker pode encontrar no máximo esta fração do restante
    remaining = MAX_SOLUTIONS - initial_total
    max_per_worker = (remaining // NUM_PARTITIONS) + CHECKPOINT_INTERVAL

    # Manager para queue e evento compartilhados entre processos
    manager = mp.Manager()
    result_queue = manager.Queue(maxsize=2000)
    stop_event = manager.Event()

    print(f"Submetendo {NUM_PARTITIONS} partições ao pool de {MAX_WORKERS} workers...\n")

    loop = asyncio.get_running_loop()

    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = [
            loop.run_in_executor(
                executor,
                worker_search,
                i,
                partitions[i],
                result_queue,
                stop_event,
                max_per_worker,
            )
            for i in range(NUM_PARTITIONS)
        ]

        # Coordinator corre concorrentemente no mesmo event loop
        total, all_signatures = await coordinator(result_queue, all_signatures, initial_total)

        # Sinaliza workers para encerrarem
        stop_event.set()

        # Aguarda todos os workers (eles verificam stop_event entre iterações)
        await asyncio.gather(*futures, return_exceptions=True)

    # Checkpoint e análise final
    save_checkpoint(all_signatures, total)
    analyze_results(total)
    _export(all_signatures, total)


def _export(all_signatures, total):
    """Exporta assinaturas_8x8.json com mesma estrutura do 6×6."""
    nodes, adj, edges = build_graph()
    parent, _level, _order, tree_edges = compute_spanning_tree(adj)
    loops = extract_loops(parent, tree_edges, edges)

    data = {
        "descricao": {
            "tabuleiro": "8x8",
            "V": len(nodes),
            "E": len(edges),
            "H1": len(loops),
            "amostras": total,
            "max_solutions_configurado": MAX_SOLUTIONS,
            "modo": "amostral",
            "nota": (
                "Enumeração completa inviável (~10^13 ciclos). "
                "Este arquivo contém representantes canônicos D4."
            ),
        },
        "assinaturas_booleanas": [list(sig) for sig in all_signatures],
    }
    with open(OUTPUT_FILE, "w") as f:
        json.dump(data, f)
    print(f"Exportado: {OUTPUT_FILE} ({total} assinaturas, H1={len(loops)}).")


# ==========================================================================
# ENTRY POINT
# ==========================================================================
def main():
    # Garante que subprocessos usem 'fork' no Linux para herdar globals Z3
    try:
        mp.set_start_method("fork")
    except RuntimeError:
        pass  # já definido
    asyncio.run(async_main())


if __name__ == "__main__":
    main()
