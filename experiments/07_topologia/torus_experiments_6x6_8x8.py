#!/usr/bin/env python3
"""
Torus Knight's Tour Topology — Experiments A, B, and C
======================================================
Investigating the relation between local degree-2 corners and global topological obstruction Q.
"""

import numpy as np
from itertools import combinations
import math

# ============================================================
# UTILITIES FOR GRAPH AND GF2 RANK (General for N x N)
# ============================================================

KNIGHT_MOVES = [(1,2),(1,-2),(-1,2),(-1,-2),(2,1),(2,-1),(-2,1),(-2,-1)]

class TorusExperimenter:
    def __init__(self, n):
        self.N = n
        self.V = n * n
        self.build_graph()

    def vid(self, r, c):
        return r * self.N + c

    def build_graph(self):
        # 1. Torus edge set (allowing multiple moves between same vertices mod N if N is small, though N>=6 has no double edges)
        torus_edge_set = {}
        for r in range(self.N):
            for c in range(self.N):
                for dr, dc in KNIGHT_MOVES:
                    r2_raw, c2_raw = r + dr, c + dc
                    r2, c2 = r2_raw % self.N, c2_raw % self.N
                    wraps_r = (r2_raw != r2)
                    wraps_c = (c2_raw != c2)
                    v1, v2 = self.vid(r, c), self.vid(r2, c2)
                    if v1 == v2:
                        continue
                    key = (min(v1,v2), max(v1,v2))
                    if key not in torus_edge_set:
                        torus_edge_set[key] = []
                    torus_edge_set[key].append((wraps_r, wraps_c, dr, dc, r, c, r2, c2))

        # 2. Plane edge set
        plane_edge_set = set()
        for r in range(self.N):
            for c in range(self.N):
                for dr, dc in KNIGHT_MOVES:
                    r2, c2 = r + dr, c + dc
                    if 0 <= r2 < self.N and 0 <= c2 < self.N:
                        v1, v2 = self.vid(r, c), self.vid(r2, c2)
                        plane_edge_set.add((min(v1,v2), max(v1,v2)))

        edges_to_remove = set(torus_edge_set.keys()) - plane_edge_set

        self.all_edges = []
        for (v1, v2), moves in sorted(torus_edge_set.items()):
            r1, c1, r2, c2 = v1 // self.N, v1 % self.N, v2 // self.N, v2 % self.N
            if (v1, v2) in edges_to_remove:
                wrap_types = set()
                for wr, wc, *_ in moves:
                    if wr and wc: wrap_types.add('diag')
                    elif wr and not wc: wrap_types.add('y-only')
                    elif wc and not wr: wrap_types.add('x-only')
                cat = 'x-only' if 'x-only' in wrap_types else ('y-only' if 'y-only' in wrap_types else 'diag')
            else:
                cat = 'internal'
            self.all_edges.append(((r1,c1), (r2,c2), v1, v2, cat))

        # Build wrap indices
        self.wrap_indices = []
        self.wrap_cats = []
        for i, e in enumerate(self.all_edges):
            if e[4] != 'internal':
                self.wrap_indices.append(i)
                self.wrap_cats.append(e[4])

    def gf2_rank(self, M):
        if M.size == 0 or M.shape[0] == 0 or M.shape[1] == 0:
            return 0
        M = M.copy().astype(np.int64) % 2
        rows, cols = M.shape
        pivot_row = 0
        for col in range(cols):
            found = -1
            for row in range(pivot_row, rows):
                if M[row, col] & 1:
                    found = row
                    break
            if found == -1:
                continue
            M[[pivot_row, found]] = M[[found, pivot_row]]
            for row in range(rows):
                if row != pivot_row and M[row, col] & 1:
                    M[row] ^= M[pivot_row]
            pivot_row += 1
        return pivot_row

    def compute_Q(self, active_mask):
        active_idx = np.where(active_mask)[0]
        E = len(active_idx)
        if E == 0:
            return 0, [], {}
        
        amap = {old: new for new, old in enumerate(active_idx)}
        degree = np.zeros(self.V, dtype=int)
        adj = {v: [] for v in range(self.V)}
        for idx in active_idx:
            e = self.all_edges[idx]
            v1, v2 = e[2], e[3]
            degree[v1] += 1
            degree[v2] += 1
            adj[v1].append(idx)
            adj[v2].append(idx)
        
        deg2 = [v for v in range(self.V) if degree[v] == 2]
        
        inc = np.zeros((self.V, E), dtype=np.int64)
        for new, old in enumerate(active_idx):
            e = self.all_edges[old]
            inc[e[2], new] = 1
            inc[e[3], new] = 1
        rank_inc = self.gf2_rank(inc)
        
        if len(deg2) < 2:
            return 0, [(v // self.N, v % self.N) for v in deg2], {v: int(degree[v]) for v in range(self.V) if degree[v] > 0}
        
        forced = [adj[v][0] for v in deg2]
        xvecs = []
        for i in range(len(forced)):
            for j in range(i+1, len(forced)):
                vec = np.zeros(E, dtype=np.int64)
                vec[amap[forced[i]]] = 1
                vec[amap[forced[j]]] = 1
                xvecs.append(vec)
        
        if not xvecs:
            return 0, [(v // self.N, v % self.N) for v in deg2], {v: int(degree[v]) for v in range(self.V) if degree[v] > 0}
        
        combined = np.vstack([inc, np.array(xvecs, dtype=np.int64)])
        Q = self.gf2_rank(combined) - rank_inc
        return Q, [(v // self.N, v % self.N) for v in deg2], {v: int(degree[v]) for v in range(self.V) if degree[v] > 0}

    def make_mask(self, removed_wrap_idx_list):
        """removed_wrap_idx_list: indices into self.wrap_indices."""
        mask = np.ones(len(self.all_edges), dtype=bool)
        for wi in removed_wrap_idx_list:
            mask[self.wrap_indices[wi]] = False
        return mask


# ============================================================
# EXPERIMENT A: CORNER SURGERY (N=6)
# ============================================================

def run_experiment_a():
    print("=" * 80)
    print(" EXPERIMENTO A: CIRURGIA NOS CANTOS (TORO 6×6)")
    print("=" * 80)
    
    exp6 = TorusExperimenter(6)
    corners = [(0,0), (0,5), (5,0), (5,5)]
    
    print("\n1. Análise de Arestas Incidentes aos Cantos:")
    corner_wraps_dict = {c: [] for c in corners}
    
    for r, c in corners:
        v = exp6.vid(r, c)
        incident_edges = []
        for i, e in enumerate(exp6.all_edges):
            if e[2] == v or e[3] == v:
                incident_edges.append((i, e))
        
        print(f"\nCanto ({r},{c}):")
        wrap_count = 0
        internal_count = 0
        for i, e in incident_edges:
            other_v = e[3] if e[2] == v else e[2]
            other_r, other_c = other_v // 6, other_v % 6
            cat = e[4]
            if cat == 'internal':
                internal_count += 1
            else:
                wrap_count += 1
                # Find the index of this wrap in wrap_indices
                wi = exp6.wrap_indices.index(i)
                corner_wraps_dict[(r,c)].append(wi)
            print(f"  -> Aresta para ({other_r},{other_c}): {cat}")
            
        print(f"  Resumo: Total = {len(incident_edges)} edges (Internas = {internal_count}, Wraps = {wrap_count})")
        print(f"  Grau no Toro = {len(incident_edges)} (esperado: 8), Grau no Plano = {internal_count} (esperado: 2)")

    print("\n2. Verificação de Compartilhamento de Wraps entre Cantos:")
    sharing_found = False
    all_pairs = list(combinations(corners, 2))
    for c1, c2 in all_pairs:
        set1 = set(corner_wraps_dict[c1])
        set2 = set(corner_wraps_dict[c2])
        shared = set1.intersection(set2)
        if len(shared) > 0:
            print(f"  AVISO: Canto {c1} e Canto {c2} compartilham {len(shared)} wrap(s): {shared}")
            sharing_found = True
            
    if not sharing_found:
        print("  ✓ NENHUMA wrap-edge é compartilhada entre os cantos!")
        print("  Conclusão: O conjunto de wraps dos cantos é perfeitamente disjunto.")
        
    # Gather all corner wraps
    all_corner_wraps = []
    for c in corners:
        all_corner_wraps.extend(corner_wraps_dict[c])
    all_corner_wraps = sorted(list(set(all_corner_wraps)))
    
    print(f"\n3. EXPERIMENTO CHAVE: Remoção de todas as wraps dos 4 cantos")
    print(f"  Total de wraps incidentes aos cantos: {len(all_corner_wraps)} (esperado: 24)")
    
    mask_24 = exp6.make_mask(all_corner_wraps)
    Q, deg2, degs = exp6.compute_Q(mask_24)
    
    print(f"  Q calculado = {Q}")
    print(f"  deg2_count = {len(deg2)}")
    print(f"  Vértices de grau 2: {deg2}")
    
    print(f"\nPERGUNTA: Q = 0, 1, 2, ou 3?")
    print(f"RESPOSTA: Q = {Q}")
    
    return exp6, all_corner_wraps, Q, deg2


# ============================================================
# EXPERIMENT B: WHY 41 AND NOT 24? (N=6)
# ============================================================

def run_experiment_b(exp6, all_corner_wraps, Q_init, deg2_init):
    print("\n" + "=" * 80)
    print(" EXPERIMENTO B: POR QUE 41 E NÃO 24? (TORO 6×6)")
    print("=" * 80)
    
    print(f"\nInício: 24 wraps removidas (Q={Q_init}, deg2={len(deg2_init)} {deg2_init})")
    
    # Greedy search from the 24 wraps
    current_set = list(all_corner_wraps)
    remaining = [w for w in range(len(exp6.wrap_indices)) if w not in current_set]
    
    greedy_steps = []
    Q_curr = Q_init
    step = 0
    
    while remaining and Q_curr < 3:
        best_wi = None
        best_Q = -1
        best_d2 = None
        
        for wi in remaining:
            trial = current_set + [wi]
            Q, d2, _ = exp6.compute_Q(exp6.make_mask(trial))
            if Q > best_Q or (Q == best_Q and best_d2 is not None and len(d2) > len(best_d2)):
                best_Q = Q
                best_wi = wi
                best_d2 = d2
        
        current_set.append(best_wi)
        remaining.remove(best_wi)
        Q_curr = best_Q
        step += 1
        
        ei = exp6.wrap_indices[best_wi]
        e = exp6.all_edges[ei]
        cat = e[4]
        
        # Check incident vertices of the added wrap edge
        (r1, c1), (r2, c2) = e[0], e[1]
        
        greedy_steps.append((step, cat, (r1,c1), (r2,c2), best_Q, len(best_d2), best_d2))
        print(f"  Passo {step:2d}: +{cat:7s} {e[0]}→{e[1]}  → Q={best_Q}, deg2={len(best_d2)} {best_d2}, total_removidas={len(current_set)}")

    print(f"\nResumo da Busca Gulosa:")
    print(f"  - Adições necessárias além das 24 iniciais para obter Q=3: {step}")
    print(f"  - Total de wraps removidas: {len(current_set)}")
    
    print(f"\n1. Análise das Arestas Adicionadas:")
    for step_id, cat, src, dst, Q_val, d2c, d2v in greedy_steps:
        # Determine if incident to internal or border vertices
        # In a 6x6 grid, borders are rows/cols 0 or 5. Internals are rows/cols 1, 2, 3, 4.
        def is_border(r, c):
            return r == 0 or r == 5 or c == 0 or c == 5
        def classify_v(r, c):
            if (r,c) in [(0,0),(0,5),(5,0),(5,5)]:
                return "Canto"
            elif is_border(r, c):
                return "Borda"
            else:
                return "Interno"
                
        t1, t2 = classify_v(src[0], src[1]), classify_v(dst[0], dst[1])
        print(f"  Aresta {src}→{dst} ({cat}): conecta [{t1}] a [{t2}]")
        
    print(f"\n2. Verificação do deg2_count durante a busca:")
    deg2_changes = [deg2_init] + [s[6] for s in greedy_steps]
    deg2_counts = [len(deg2_init)] + [s[5] for s in greedy_steps]
    print(f"  Evolução de deg2_count: {deg2_counts}")
    
    print(f"\nINTERPRETAÇÃO CONCEITUAL:")
    print("  As wraps extras de fato alteram o espaço quociente F2^E / rowspace(∂1)?")
    print(f"  -> Q passa de {Q_init} para 3, enquanto o deg2_count evolui de {len(deg2_init)} para {deg2_counts[-1]}.")
    if deg2_counts == [4] * len(deg2_counts):
        print("  -> SIM! O deg2_count permaneceu EXATAMENTE 4 (apenas os 4 cantos) durante todo o processo.")
        print("  Isso prova que as 32 wraps adicionais removidas NÃO criaram novos vértices de grau 2,")
        print("  mas foram matematicamente necessárias para desacoplar as restrições globais no quociente homológico!")
    else:
        print(f"  -> O deg2_count variou: {deg2_counts}")
        print("  -> Analisando os conjuntos de grau-2 em cada passo:")
        for s_idx, d2_set in enumerate(deg2_changes):
            print(f"     Passo {s_idx}: {d2_set}")


# ============================================================
# EXPERIMENT C: SCALING — TORO 8×8
# ============================================================

def run_experiment_c():
    print("\n" + "=" * 80)
    print(" EXPERIMENTO C: SCALING — TORO 8×8")
    print("=" * 80)
    
    exp8 = TorusExperimenter(8)
    
    print(f"\n1. Classificação das Arestas no Toro 8×8:")
    print(f"  Total de arestas no Toro 8×8: {len(exp8.all_edges)}")
    
    cats = [e[4] for e in exp8.all_edges]
    print(f"  - Internas: {cats.count('internal')}")
    print(f"  - Wraps: {len(exp8.wrap_indices)} no total:")
    print(f"    * x-only: {cats.count('x-only')}")
    print(f"    * y-only: {cats.count('y-only')}")
    print(f"    * diag: {cats.count('diag')}")
    
    # 2. Incident wraps for 8x8 corners
    corners8 = [(0,0), (0,7), (7,0), (7,7)]
    print("\n2. Wraps incidentes aos 4 Cantos no Toro 8×8:")
    corner8_wraps_dict = {c: [] for c in corners8}
    
    for r, c in corners8:
        v = exp8.vid(r, c)
        incident = []
        for i, e in enumerate(exp8.all_edges):
            if e[2] == v or e[3] == v:
                incident.append((i, e))
                
        print(f"\nCanto ({r},{c}):")
        w_counts = {'x-only': 0, 'y-only': 0, 'diag': 0, 'internal': 0}
        for i, e in incident:
            w_counts[e[4]] += 1
            if e[4] != 'internal':
                wi = exp8.wrap_indices.index(i)
                corner8_wraps_dict[(r,c)].append(wi)
        print(f"  Incidentes: total={len(incident)} (internas={w_counts['internal']}, x-only={w_counts['x-only']}, y-only={w_counts['y-only']}, diag={w_counts['diag']})")

    # Check sharing in 8x8
    sharing8 = False
    all_pairs8 = list(combinations(corners8, 2))
    for c1, c2 in all_pairs8:
        shared = set(corner8_wraps_dict[c1]).intersection(set(corner8_wraps_dict[c2]))
        if len(shared) > 0:
            sharing8 = True
            print(f"  AVISO: Canto {c1} e Canto {c2} compartilham {len(shared)} wraps!")
    if not sharing8:
        print("  ✓ Cantos do toro 8×8 também possuem conjuntos de wraps perfeitamente disjuntos.")

    # 3. Compute Q for different 8x8 configurations
    print("\n3. Computação de Q no Toro 8×8:")
    
    # a) Toro puro
    mask_pure = exp8.make_mask([])
    Q_pure, d2_pure, _ = exp8.compute_Q(mask_pure)
    print(f"  a) Toro Puro: Q = {Q_pure}, deg2 = {len(d2_pure)}")
    
    # b) Plano 8x8 (remove todas as wraps)
    mask_plane = exp8.make_mask(range(len(exp8.wrap_indices)))
    Q_plane, d2_plane, _ = exp8.compute_Q(mask_plane)
    print(f"  b) Plano 8×8 (todas as wraps removidas): Q = {Q_plane}, deg2 = {len(d2_plane)} {d2_plane}")
    
    # c) Remove exatamente as wraps dos 4 cantos (4 x 6 = 24 wraps)
    all_corner8_wraps = []
    for c in corners8:
        all_corner8_wraps.extend(corner8_wraps_dict[c])
    all_corner8_wraps = sorted(list(set(all_corner8_wraps)))
    
    print(f"  c) Remoção EXCLUSIVA das {len(all_corner8_wraps)} wraps dos 4 cantos:")
    mask_c = exp8.make_mask(all_corner8_wraps)
    Q_c, d2_c, _ = exp8.compute_Q(mask_c)
    print(f"     Q = {Q_c}")
    print(f"     deg2_count = {len(d2_c)}")
    print(f"     Vértices de grau 2: {d2_c}")


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":
    exp6, all_corner_wraps, Q_init, deg2_init = run_experiment_a()
    run_experiment_b(exp6, all_corner_wraps, Q_init, deg2_init)
    run_experiment_c()
