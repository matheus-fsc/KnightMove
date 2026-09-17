#!/usr/bin/env python3
"""
cavalo_loop_destruicao_6x6.py
==============================
Adaptado de cavalo_engine_8x8_path_v3.py

O QUE FAZ:
  - Enumera TODOS os caminhos hamiltonianos fechados no tabuleiro 6x6
    via backtracking exaustivo (viável: ~9.862 soluções)
  - Para cada solução, registra a sequência de destruição de loops
    fundamentais a cada passo (usando o número ciclomático H1)
  - Cataloga por casa: quantos loops são destruídos quando aquela
    casa é visitada, em qual passo, e em qual contexto topológico
  - Identifica loops determinísticos (frequência 100%) vs livres
  - Detecta pares de exclusão mútua (correlação negativa forte)

SAÍDA:
  - destruction_catalogue.json   : catálogo completo por solução
  - destruction_map.json         : mapa por casa (destruição média, max, por passo)
  - loop_correlations.json       : pares de loops com correlação negativa forte
  - summary_report.txt           : relatório legível
"""

import json, math, time, collections
from itertools import combinations

# ── configuração ──────────────────────────────────────────────────────────────

BOARD = 6
TOTAL = BOARD * BOARD
MOVES = [(2,1),(2,-1),(-2,1),(-2,-1),(1,2),(1,-2),(-1,2),(-1,-2)]

# loops obrigatórios conhecidos do paper (para validação cruzada)
KNOWN_MANDATORY = [7, 34, 44]


# ── utilidades de coordenadas ─────────────────────────────────────────────────

def vid(r, c):
    """Vértice como inteiro único."""
    return r * BOARD + c

def vrc(v):
    """Inteiro -> (row, col)."""
    return divmod(v, BOARD)

def label(v):
    """Notação xadrez: A6, B3, etc."""
    r, c = vrc(v)
    return chr(65 + c) + str(BOARD - r)

def coord_str(v):
    r, c = vrc(v)
    return f"({r},{c})"


# ── construção do grafo ───────────────────────────────────────────────────────

def build_graph():
    """
    Retorna:
      ADJ  : lista de adjacência [lista de vértices vizinhos]
      EDGES: conjunto de arestas como pares (u,v) com u < v
    """
    ADJ = [[] for _ in range(TOTAL)]
    EDGES = set()
    for r in range(BOARD):
        for c in range(BOARD):
            v = vid(r, c)
            for dr, dc in MOVES:
                nr, nc = r + dr, c + dc
                if 0 <= nr < BOARD and 0 <= nc < BOARD:
                    u = vid(nr, nc)
                    ADJ[v].append(u)
                    if u > v:
                        EDGES.add((v, u))
    return ADJ, EDGES

ADJ, EDGES = build_graph()
EDGES_LIST = sorted(EDGES)


# ── número ciclomático H1 ─────────────────────────────────────────────────────

def count_components(alive_set):
    """Número de componentes conexas no subgrafo induzido por alive_set."""
    if not alive_set:
        return 0
    visited = set()
    comps = 0
    for start in alive_set:
        if start in visited:
            continue
        comps += 1
        stack = [start]
        while stack:
            v = stack.pop()
            if v in visited:
                continue
            visited.add(v)
            for u in ADJ[v]:
                if u in alive_set and u not in visited:
                    stack.append(u)
    return comps

def cyclomatic_dim(alive_set):
    """
    H1(G) = |E| - |V| + componentes_conexas
    Dimensão do espaço de ciclos fundamentais no subgrafo residual.
    """
    V = len(alive_set)
    if V == 0:
        return 0
    E = sum(
        1 for v in alive_set
        for u in ADJ[v]
        if u in alive_set and u > v
    )
    C = count_components(alive_set)
    return max(0, E - V + C)

def loops_destroyed(v, alive_set):
    """
    Quantos loops fundamentais são destruídos ao remover v de alive_set.
    ΔH1 = H1(G) - H1(G \\ {v}) = deg(v) - 1 + ΔC
    onde ΔC é a variação no número de componentes.
    """
    before = cyclomatic_dim(alive_set)
    next_alive = alive_set - {v}
    after = cyclomatic_dim(next_alive)
    return before - after


# ── busca exaustiva de caminhos hamiltonianos fechados ────────────────────────

def find_all_closed_tours(start, progress_every=500000):
    """
    Backtracking completo: encontra todos os ciclos hamiltonianos
    que partem e voltam a 'start'.

    Otimizações:
      - degree tracking (deg[v] = nº de vizinhos vivos de v) para Warnsdorff O(1)
      - poda: candidato com deg=0 fora do passo de fechamento é beco sem saída
      - bytearray de alive em vez de set (lookup mais rápido)
      - relatório periódico para evitar percepção de "loop"

    Retorna lista de tuplas (sequência de vértices).
    """
    solutions = []
    path = [start]
    alive = bytearray([1] * TOTAL)
    alive[start] = 0

    deg = [sum(1 for u in ADJ[v] if alive[u]) for v in range(TOTAL)]

    counter = [0]
    t_start = time.time()

    def backtrack():
        counter[0] += 1
        if counter[0] % progress_every == 0:
            elapsed = time.time() - t_start
            rate = counter[0] / elapsed if elapsed > 0 else 0
            print(f"    [progresso] backtracks={counter[0]:>11,}  "
                  f"sols={len(solutions):>6,}  "
                  f"depth={len(path):2d}/{TOTAL}  "
                  f"taxa={rate:>9,.0f}/s  "
                  f"elapsed={elapsed:5.0f}s")

        n = len(path)
        current = path[-1]

        if n == TOTAL:
            if start in ADJ[current]:
                solutions.append(tuple(path))
            return

        candidates = [u for u in ADJ[current] if alive[u]]
        if not candidates:
            return

        # Warnsdorff: ordenar por grau no subgrafo vivo (ascendente)
        candidates.sort(key=lambda u: deg[u])

        # Poda: candidato com deg=0 só é viável no passo de fechamento.
        # Caso contrário, ao visitá-lo ficaríamos sem saída.
        if deg[candidates[0]] == 0 and n < TOTAL - 1:
            return

        for nxt in candidates:
            alive[nxt] = 0
            for w in ADJ[nxt]:
                deg[w] -= 1
            path.append(nxt)
            backtrack()
            path.pop()
            alive[nxt] = 1
            for w in ADJ[nxt]:
                deg[w] += 1

    backtrack()
    return solutions


# ── catalogação de destruição por solução ─────────────────────────────────────

def catalogue_solution(solution):
    """
    Para uma solução (sequência de vértices), retorna:
      destruction_seq : lista de int, destruição em cada passo
      dim_seq         : H1 residual após cada passo
      critical_step   : passo onde a destruição é máxima
      dead_end_risk   : passo onde H1 cai pela primeira vez a 0
    """
    alive = set(range(TOTAL))
    destruction_seq = []
    dim_seq = []

    for v in solution:
        killed = loops_destroyed(v, alive)
        alive.discard(v)
        dim_after = cyclomatic_dim(alive)
        destruction_seq.append(killed)
        dim_seq.append(dim_after)

    max_kill = max(destruction_seq)
    critical_step = destruction_seq.index(max_kill)

    # primeiro passo onde H1 chega a 0 (ou nunca)
    dead_end_risk = next(
        (i for i, d in enumerate(dim_seq) if d == 0), None
    )

    return {
        "destruction_seq": destruction_seq,
        "dim_seq": dim_seq,
        "total_destroyed": sum(destruction_seq),
        "max_single_kill": max_kill,
        "critical_step": critical_step,
        "critical_vertex": label(solution[critical_step]),
        "dead_end_risk_step": dead_end_risk,
    }


# ── análise estatística ───────────────────────────────────────────────────────

def compute_loop_signatures(all_solutions):
    """
    Para análise de correlação, representamos cada solução como
    vetor booleano sobre as EDGES (aresta ativa ou não).
    Retorna matriz n_soluções × n_arestas.
    """
    edge_idx = {e: i for i, e in enumerate(EDGES_LIST)}
    n_edges = len(EDGES_LIST)
    sigs = []
    for sol in all_solutions:
        sig = [False] * n_edges
        for i in range(len(sol)):
            u = sol[i]
            v = sol[(i + 1) % len(sol)]
            e = (min(u, v), max(u, v))
            if e in edge_idx:
                sig[edge_idx[e]] = True
        sigs.append(sig)
    return sigs

def compute_correlations(sigs, top_n=20):
    """
    Calcula correlação de Pearson entre todas as dimensões (arestas).
    Retorna os top_n pares com correlação mais negativa.
    """
    n = len(sigs)
    nd = len(sigs[0]) if sigs else 0
    freq = [sum(s[i] for s in sigs) / n for i in range(nd)]

    # só analisar arestas "livres" (nem sempre ativa, nem nunca)
    live_dims = [i for i, f in enumerate(freq) if 0 < f < 1]

    corrs = []
    for i, j in combinations(live_dims, 2):
        fi, fj = freq[i], freq[j]
        cov = sum((sigs[k][i] - fi) * (sigs[k][j] - fj) for k in range(n)) / n
        si = math.sqrt(fi * (1 - fi))
        sj = math.sqrt(fj * (1 - fj))
        if si > 0 and sj > 0:
            corrs.append((round(cov / (si * sj), 4), i, j))

    corrs.sort()
    return corrs[:top_n], freq

def find_mandatory_and_impossible(freq):
    """
    Obrigatórios: freq = 1.0
    Impossíveis: freq = 0.0
    """
    mandatory = [i for i, f in enumerate(freq) if f == 1.0]
    impossible = [i for i, f in enumerate(freq) if f == 0.0]
    free = [i for i, f in enumerate(freq) if 0 < f < 1]
    return mandatory, impossible, free


# ── relatório de texto ────────────────────────────────────────────────────────

def write_report(summary, destruction_map, correlations, freq, path):
    lines = []
    w = lines.append

    w("=" * 65)
    w("CATALOGAÇÃO DE DESTRUIÇÃO DE LOOPS — PASSEIO DO CAVALO 6×6")
    w("=" * 65)
    w("")
    w(f"Tabuleiro           : {BOARD}×{BOARD}")
    w(f"H1 inicial          : {summary['H1_initial']}  (dim. espaço de ciclos)")
    w(f"Total de soluções   : {summary['total_solutions']}")
    w("")

    w("── SOLUÇÕES POR CASA DE PARTIDA ──────────────────────────")
    for lbl, cnt in sorted(summary['solutions_by_start'].items()):
        bar = "█" * (cnt * 30 // max(summary['solutions_by_start'].values()))
        w(f"  {lbl}: {cnt:4d}  {bar}")
    w("")

    w("── DESTRUIÇÃO GLOBAL ──────────────────────────────────────")
    ds = summary['destruction_stats']
    w(f"  Média total destruída por solução : {ds['mean_total_destroyed']}")
    w(f"  Média do kill máximo por solução  : {ds['mean_max_kill']}")
    w(f"  Maior kill único (ever)           : {ds['max_single_kill_ever']}")
    w("")

    w("── PASSO DE PICO DE DESTRUIÇÃO ────────────────────────────")
    peak = summary['peak_step_distribution']
    max_cnt = max(peak.values()) if peak else 1
    for step in sorted(peak.keys()):
        cnt = peak[step]
        bar = "█" * (cnt * 40 // max_cnt)
        w(f"  passo {int(step):2d}: {cnt:5d}  {bar}")
    w("")

    w("── ARESTAS OBRIGATÓRIAS (freq=100%) ───────────────────────")
    mandatory, impossible, free = find_mandatory_and_impossible(freq)
    w(f"  {len(mandatory)} arestas obrigatórias")
    for i in mandatory:
        u, v = EDGES_LIST[i]
        w(f"    aresta {i:3d}: {label(u)}-{label(v)}")
    w("")

    w("── ARESTAS IMPOSSÍVEIS (freq=0%) ──────────────────────────")
    w(f"  {len(impossible)} arestas impossíveis")
    for i in impossible:
        u, v = EDGES_LIST[i]
        w(f"    aresta {i:3d}: {label(u)}-{label(v)}")
    w("")

    w("── TOP 10 CORRELAÇÕES NEGATIVAS ───────────────────────────")
    w("  (pares de arestas com exclusão mútua)")
    for r, i, j in correlations[:10]:
        ui, vi = EDGES_LIST[i]
        uj, vj = EDGES_LIST[j]
        w(f"  r={r:7.4f}  aresta{i:3d}({label(ui)}-{label(vi)}) ↔ aresta{j:3d}({label(uj)}-{label(vj)})")
    w("")

    w("── CASAS MAIS PERIGOSAS (destruição média alta) ───────────")
    danger = sorted(
        destruction_map.items(),
        key=lambda x: -x[1]['mean_destruction']
    )
    for lbl, data in danger[:10]:
        w(f"  {lbl}  grau={data['degree']:2d}  "
          f"média={data['mean_destruction']:.3f}  "
          f"max={data['max_destruction']}")
    w("")

    w("── CASAS MAIS SEGURAS (destruição média baixa) ────────────")
    for lbl, data in danger[-10:]:
        w(f"  {lbl}  grau={data['degree']:2d}  "
          f"média={data['mean_destruction']:.3f}  "
          f"max={data['max_destruction']}")
    w("")

    w("=" * 65)

    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))

    # também imprime no terminal
    print("\n".join(lines))


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    t0 = time.time()

    print("=" * 65)
    print("CATALOGADOR DE DESTRUIÇÃO DE LOOPS — CAVALO 6×6")
    print("=" * 65)

    full_alive = set(range(TOTAL))
    H1_initial = cyclomatic_dim(full_alive)
    print(f"V={TOTAL}  E={len(EDGES)}  H1={H1_initial}")
    print()

    # ── 1. ENCONTRAR TODAS AS SOLUÇÕES ────────────────────────────────────────
    print("Fase 1: Busca exaustiva de ciclos hamiltonianos...")
    print("(Busca única ancorada em v=0; demais casas via rotação)")
    print("(Warnsdorff + degree tracking + poda de dead-ends)\n")

    # Todo ciclo hamiltoniano fechado passa por todos os 36 vértices, então
    # basta procurar a partir de v=0. Cada ciclo não-direcionado aparece duas
    # vezes (uma por direção). As variantes por casa de partida são geradas
    # depois por rotação cíclica.
    t_s = time.time()
    directed_cycles_at_0 = find_all_closed_tours(0)
    elapsed = time.time() - t_s
    print(f"\n  v=0 ({label(0)}): {len(directed_cycles_at_0)} ciclos direcionados  [{elapsed:.1f}s]")
    print(f"  → {len(directed_cycles_at_0)//2} ciclos não-direcionados (esperado ≈9862)")

    # Gera as variantes por casa de partida via rotação.
    # Cada ciclo (v_0, v_1, ..., v_{N-1}) rotacionado a partir do índice k
    # vira (v_k, v_{k+1}, ..., v_{N-1}, v_0, ..., v_{k-1}), o que equivale
    # a "começar pela casa v_k" sem refazer o backtracking.
    print("\n  Gerando rotações (uma por casa de partida)...")
    all_solutions = []
    solutions_by_start = {label(v): 0 for v in range(TOTAL)}
    for cyc in directed_cycles_at_0:
        for k in range(TOTAL):
            rotated = cyc[k:] + cyc[:k]
            all_solutions.append(rotated)
            solutions_by_start[label(rotated[0])] += 1

    print(f"\nTotal: {len(all_solutions)} sequências (todas rotações)  |  "
          f"Tempo fase 1: {time.time()-t0:.1f}s\n")

    # ── 2. CATALOGAR DESTRUIÇÕES ──────────────────────────────────────────────
    print("Fase 2: Catalogando destruição de loops...")

    catalogue = []
    por_casa = {v: collections.defaultdict(list) for v in range(TOTAL)}

    for idx, sol in enumerate(all_solutions):
        data = catalogue_solution(sol)
        record = {
            "id": idx,
            "start": label(sol[0]),
            "sequence": [label(v) for v in sol],
            **data,
        }
        catalogue.append(record)

        # acumular por casa e por passo
        for step, v in enumerate(sol):
            por_casa[v][step].append(data["destruction_seq"][step])

        if (idx + 1) % 1000 == 0:
            print(f"  {idx+1}/{len(all_solutions)}...")

    print(f"  Catalogação completa: {len(catalogue)} registros\n")

    # ── 3. MAPA DE DESTRUIÇÃO POR CASA ────────────────────────────────────────
    print("Fase 3: Construindo mapa de destruição por casa...")

    destruction_map = {}
    for v in range(TOTAL):
        all_kills = []
        by_step = {}
        for step, kills in por_casa[v].items():
            by_step[str(step)] = {
                "n_appearances": len(kills),
                "mean": round(sum(kills) / len(kills), 4),
                "max": max(kills),
                "min": min(kills),
                "distribution": {
                    str(k): cnt
                    for k, cnt in sorted(
                        collections.Counter(kills).items()
                    )
                },
            }
            all_kills.extend(kills)

        if all_kills:
            destruction_map[label(v)] = {
                "vertex_id": v,
                "label": label(v),
                "row_col": list(vrc(v)),
                "degree": len(ADJ[v]),
                "total_appearances": len(all_kills),
                "mean_destruction": round(sum(all_kills) / len(all_kills), 4),
                "max_destruction": max(all_kills),
                "min_destruction": min(all_kills),
                "destruction_distribution": {
                    str(k): cnt
                    for k, cnt in sorted(
                        collections.Counter(all_kills).items()
                    )
                },
                "by_step": by_step,
            }

    # ── 4. CORRELAÇÕES ENTRE LOOPS (ARESTAS) ──────────────────────────────────
    print("Fase 4: Calculando correlações entre arestas (loops)...")

    sigs = compute_loop_signatures(all_solutions)
    correlations, freq = compute_correlations(sigs, top_n=50)

    # pares com r ≈ -0.771 (invariante do paper)
    target = -0.771
    invariant_pairs = [
        (r, i, j) for r, i, j in correlations
        if abs(r - target) < 0.05
    ]
    print(f"  Pares com r ≈ {target}: {len(invariant_pairs)}")

    # ── 5. ESTATÍSTICAS GLOBAIS ────────────────────────────────────────────────
    all_total = [r["total_destroyed"] for r in catalogue]
    all_max = [r["max_single_kill"] for r in catalogue]
    all_peak_step = [r["critical_step"] for r in catalogue]

    peak_dist = dict(collections.Counter(all_peak_step))

    summary = {
        "board": f"{BOARD}x{BOARD}",
        "H1_initial": H1_initial,
        "total_edges": len(EDGES),
        "total_solutions": len(all_solutions),
        "solutions_by_start": solutions_by_start,
        "destruction_stats": {
            "mean_total_destroyed": round(sum(all_total) / len(all_total), 3),
            "mean_max_kill": round(sum(all_max) / len(all_max), 3),
            "max_single_kill_ever": max(all_max),
        },
        "peak_step_distribution": {str(k): v for k, v in sorted(peak_dist.items())},
        "mandatory_edges": [
            {"idx": i, "edge": [label(u) for u in EDGES_LIST[i]]}
            for i, f in enumerate(freq) if f == 1.0
        ],
        "impossible_edges": [
            {"idx": i, "edge": [label(u) for u in EDGES_LIST[i]]}
            for i, f in enumerate(freq) if f == 0.0
        ],
        "invariant_pairs_r771": [
            {
                "r": r,
                "edge_a": {"idx": i, "label": f"{label(EDGES_LIST[i][0])}-{label(EDGES_LIST[i][1])}"},
                "edge_b": {"idx": j, "label": f"{label(EDGES_LIST[j][0])}-{label(EDGES_LIST[j][1])}"},
            }
            for r, i, j in invariant_pairs
        ],
    }

    # ── 6. SALVAR RESULTADOS ───────────────────────────────────────────────────
    print("\nSalvando resultados...")

    with open("destruction_catalogue.json", "w") as f:
        json.dump(catalogue, f, separators=(',', ':'))
    print("  destruction_catalogue.json")

    with open("destruction_map.json", "w") as f:
        json.dump(destruction_map, f, indent=2, ensure_ascii=False)
    print("  destruction_map.json")

    corr_out = {
        "edge_labels": [f"{label(u)}-{label(v)}" for u, v in EDGES_LIST],
        "top_negative_correlations": [
            {
                "r": r,
                "dim_i": i, "edge_i": f"{label(EDGES_LIST[i][0])}-{label(EDGES_LIST[i][1])}",
                "dim_j": j, "edge_j": f"{label(EDGES_LIST[j][0])}-{label(EDGES_LIST[j][1])}",
            }
            for r, i, j in correlations
        ],
        "invariant_pairs_r771": summary["invariant_pairs_r771"],
    }
    with open("loop_correlations.json", "w") as f:
        json.dump(corr_out, f, indent=2, ensure_ascii=False)
    print("  loop_correlations.json")

    with open("summary.json", "w") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)
    print("  summary.json")

    # ── 7. RELATÓRIO LEGÍVEL ───────────────────────────────────────────────────
    write_report(summary, destruction_map, correlations, freq, "summary_report.txt")
    print("  summary_report.txt")

    print(f"\nTempo total: {time.time()-t0:.1f}s")
    print("Pronto.")


if __name__ == "__main__":
    main()
