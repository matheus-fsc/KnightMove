#!/usr/bin/env python3
"""
Generate two figures illustrating GF(2) cycles on the knight graph.

Figure 1: DFS tree — two paths converge on the same vertex → fundamental cycle
Figure 2: XOR of two cycles → shared edges cancel, producing a new cycle
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
import numpy as np
from collections import defaultdict, deque

# ──────────────────────────────────────────
# Knight graph on a small board
# ──────────────────────────────────────────
OFFSETS = [(1,2),(2,1),(2,-1),(1,-2),(-1,-2),(-2,-1),(-2,1),(-1,2)]

def knight_graph(n):
    adj = defaultdict(list)
    for r in range(n):
        for c in range(n):
            for dr, dc in OFFSETS:
                nr, nc = r+dr, c+dc
                if 0 <= nr < n and 0 <= nc < n:
                    adj[(r,c)].append((nr,nc))
    return adj

def board_pos(r, c, n):
    """(row, col) → plot (x, y), row 0 at top."""
    return (c * 1.1, (n - 1 - r) * 1.1)

# ═══════════════════════════════════════════
# FIGURE 1: DFS tree + fundamental cycle
# ═══════════════════════════════════════════
def gen_figure1():
    n = 5
    adj = knight_graph(n)

    # Run DFS from (0,0) — iterative, recording discovery order
    start = (2, 2)  # center of board → richer DFS tree, better back-edges
    visited = {}  # vertex → discovery index
    parent = {}
    tree_edges = []
    discovery = 0

    def dfs(v):
        nonlocal discovery
        visited[v] = discovery
        discovery += 1
        for w in sorted(adj[v]):
            if w not in visited:
                parent[w] = v
                tree_edges.append((v, w))
                dfs(w)

    parent[start] = None
    dfs(start)

    # Classify non-tree edges
    tree_set = set()
    for u, v in tree_edges:
        tree_set.add((u,v))
        tree_set.add((v,u))

    back_edges = []
    for v in adj:
        for w in adj[v]:
            if (v,w) not in tree_set and v < w:
                back_edges.append((v, w))

    # Find a back-edge producing a cycle where BOTH paths to the LCA are ≥ 2
    best = None
    for u, w in back_edges:
        # Trace ancestors
        anc_u = []
        node = u
        while node is not None:
            anc_u.append(node)
            node = parent.get(node)

        anc_w = []
        node = w
        while node is not None:
            anc_w.append(node)
            node = parent.get(node)

        set_u = set(anc_u)
        lca = None
        for x in anc_w:
            if x in set_u:
                lca = x
                break

        # Path from u → LCA
        path_a = []
        node = u
        while node != lca:
            path_a.append(node)
            node = parent[node]
        path_a.append(lca)

        # Path from w → LCA
        path_b = []
        node = w
        while node != lca:
            path_b.append(node)
            node = parent[node]
        path_b.append(lca)

        la, lb = len(path_a), len(path_b)
        cycle_len = la + lb - 1  # LCA counted once

        # We want both paths ≥ 2 edges for visual clarity
        if la >= 3 and lb >= 3 and 6 <= cycle_len <= 12:
            score = -abs(la - lb)  # prefer balanced
            if best is None or score > best[0]:
                best = (score, u, w, lca, path_a, path_b)

    if best is None:
        # Fallback: pick anything with both ≥ 2
        for u, w in back_edges:
            anc_u, anc_w = [], []
            node = u
            while node is not None: anc_u.append(node); node = parent.get(node)
            node = w
            while node is not None: anc_w.append(node); node = parent.get(node)
            set_u = set(anc_u)
            for x in anc_w:
                if x in set_u: lca = x; break
            pa, pb = [], []
            node = u
            while node != lca: pa.append(node); node = parent[node]
            pa.append(lca)
            node = w
            while node != lca: pb.append(node); node = parent[node]
            pb.append(lca)
            if len(pa) >= 2 and len(pb) >= 2:
                best = (0, u, w, lca, pa, pb)
                break

    if best is None:
        # Last resort
        u, w = back_edges[0]
        anc_u = []
        node = u
        while node is not None: anc_u.append(node); node = parent.get(node)
        anc_w = []
        node = w
        while node is not None: anc_w.append(node); node = parent.get(node)
        set_u = set(anc_u)
        for x in anc_w: 
            if x in set_u: lca = x; break
        pa, pb = [], []
        node = u
        while node != lca: pa.append(node); node = parent[node]
        pa.append(lca)
        node = w
        while node != lca: pb.append(node); node = parent[node]
        pb.append(lca)
        best = (0, u, w, lca, pa, pb)

    _, back_u, back_w, lca, path_a, path_b = best

    # Build the cycle: path_a (u→LCA) + reversed path_b[:-1] (LCA→w) 
    cycle = path_a + list(reversed(path_b[:-1]))

    # Edge sets for coloring
    path_a_edges = set()
    for i in range(len(path_a)-1):
        path_a_edges.add(tuple(sorted([path_a[i], path_a[i+1]])))

    path_b_edges = set()
    for i in range(len(path_b)-1):
        path_b_edges.add(tuple(sorted([path_b[i], path_b[i+1]])))

    back_edge = tuple(sorted([back_u, back_w]))

    print(f"  Ciclo: {cycle}  ({len(cycle)} vértices)")
    print(f"  Back-edge: {back_u} ↔ {back_w}")
    print(f"  LCA: {lca}")
    print(f"  Caminho A ({len(path_a)-1} arestas): {' → '.join(str(v) for v in path_a)}")
    print(f"  Caminho B ({len(path_b)-1} arestas): {' → '.join(str(v) for v in path_b)}")

    # ── Plot ──
    fig, axes = plt.subplots(1, 3, figsize=(21, 7.5))
    fig.patch.set_facecolor('#0d1117')

    titles = [
        '① Árvore DFS\n(arestas de árvore em cinza)',
        '② Dois caminhos → mesma casa\nBack-edge fecha o ciclo',
        '③ Ciclo fundamental em GF(2)\n$c \\in \\ker(\\partial_1) \\subset \\mathbb{F}_2^{|E|}$'
    ]
    title_colors = ['#58a6ff', '#f0883e', '#3fb950']

    all_verts = sorted(adj.keys())

    for ax_idx, ax in enumerate(axes):
        ax.set_facecolor('#0d1117')
        ax.set_aspect('equal')
        xmin, xmax = -0.7, (n-1)*1.1 + 0.7
        ymin, ymax = -0.7, (n-1)*1.1 + 0.7
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.axis('off')
        ax.set_title(titles[ax_idx], color=title_colors[ax_idx],
                     fontsize=13, fontweight='bold', pad=14, linespacing=1.3)

        # Chess board
        for r in range(n):
            for c in range(n):
                x, y = board_pos(r, c, n)
                clr = '#1a2030' if (r+c)%2==0 else '#151b25'
                rect = plt.Rectangle((x-0.48, y-0.48), 0.96, 0.96,
                                     facecolor=clr, edgecolor='#2a3040',
                                     linewidth=0.5, zorder=0)
                ax.add_patch(rect)

        # ─── Panel 1: Full DFS tree ───
        if ax_idx == 0:
            # All tree edges
            for u, v in tree_edges:
                x1,y1 = board_pos(*u, n)
                x2,y2 = board_pos(*v, n)
                e = tuple(sorted([u,v]))
                if e in path_a_edges:
                    c_, lw_, a_ = '#f97583', 3.0, 0.85
                elif e in path_b_edges:
                    c_, lw_, a_ = '#79c0ff', 3.0, 0.85
                else:
                    c_, lw_, a_ = '#3d4450', 1.0, 0.5
                ax.plot([x1,x2],[y1,y2], color=c_, lw=lw_, alpha=a_,
                        zorder=1, solid_capstyle='round')

            # Back edge (subtle dashed)
            x1,y1 = board_pos(*back_u, n)
            x2,y2 = board_pos(*back_w, n)
            ax.plot([x1,x2],[y1,y2], color='#d2a8ff', lw=2, ls='--',
                    alpha=0.4, zorder=1)

            # Vertices with DFS order
            for v in all_verts:
                x, y = board_pos(*v, n)
                in_pa = v in set(path_a)
                in_pb = v in set(path_b)
                is_lca = v == lca
                is_be = v in (back_u, back_w)

                if is_lca:
                    fc, sz, ec_ = '#f0883e', 220, '#ffffff'
                elif is_be:
                    fc, sz, ec_ = '#d2a8ff', 180, '#ffffff'
                elif in_pa:
                    fc, sz, ec_ = '#f97583', 150, '#8b949e'
                elif in_pb:
                    fc, sz, ec_ = '#79c0ff', 150, '#8b949e'
                else:
                    fc, sz, ec_ = '#252b35', 55, '#3d4450'

                ax.scatter(x, y, s=sz, c=fc, edgecolors=ec_, linewidths=1, zorder=3)

                # DFS number
                idx = visited.get(v, '')
                if v in set(path_a) or v in set(path_b) or idx != '' and idx < 12:
                    ax.text(x, y, str(idx), ha='center', va='center',
                            fontsize=6, color='white', fontweight='bold', zorder=4,
                            path_effects=[pe.withStroke(linewidth=2, foreground='#0d1117')])

                if v == start:
                    ax.annotate('raiz\n(0,0)', (x,y), textcoords="offset points",
                               xytext=(-20,-18), ha='center', fontsize=8,
                               color='#8b949e', fontweight='bold')

        # ─── Panel 2: Two paths highlighted + back-edge ───
        elif ax_idx == 1:
            # Dim tree
            for u, v in tree_edges:
                x1,y1 = board_pos(*u, n)
                x2,y2 = board_pos(*v, n)
                e = tuple(sorted([u,v]))
                if e in path_a_edges:
                    ax.plot([x1,x2],[y1,y2], color='#f97583', lw=4, alpha=0.9,
                            zorder=2, solid_capstyle='round')
                elif e in path_b_edges:
                    ax.plot([x1,x2],[y1,y2], color='#79c0ff', lw=4, alpha=0.9,
                            zorder=2, solid_capstyle='round')
                else:
                    ax.plot([x1,x2],[y1,y2], color='#3d4450', lw=0.6, alpha=0.2, zorder=0)

            # Back-edge
            x1,y1 = board_pos(*back_u, n)
            x2,y2 = board_pos(*back_w, n)
            ax.plot([x1,x2],[y1,y2], color='#d2a8ff', lw=4, ls='--',
                    alpha=0.95, zorder=2)

            # Label back-edge
            mx, my = (x1+x2)/2, (y1+y2)/2
            ax.annotate('back-edge', (mx, my), textcoords="offset points",
                        xytext=(14, 10), ha='left', fontsize=10,
                        color='#d2a8ff', fontweight='bold',
                        arrowprops=dict(arrowstyle='->', color='#d2a8ff', lw=1.5))

            # Vertices
            for v in all_verts:
                x, y = board_pos(*v, n)
                in_pa = v in set(path_a)
                in_pb = v in set(path_b)
                is_lca = v == lca
                is_be = v in (back_u, back_w)

                if is_lca:
                    fc, sz, ec_ = '#f0883e', 240, '#ffffff'
                elif is_be:
                    fc, sz, ec_ = '#d2a8ff', 200, '#ffffff'
                elif in_pa and in_pb:
                    fc, sz, ec_ = '#f0883e', 180, '#ffffff'
                elif in_pa:
                    fc, sz, ec_ = '#f97583', 160, '#ffffff'
                elif in_pb:
                    fc, sz, ec_ = '#79c0ff', 160, '#ffffff'
                else:
                    fc, sz, ec_ = '#181d25', 25, '#2a3040'

                ax.scatter(x, y, s=sz, c=fc, edgecolors=ec_, linewidths=1.2, zorder=3)

                # Label cycle vertices
                if in_pa or in_pb:
                    label = f'({v[0]},{v[1]})'
                    ax.text(x, y, label, ha='center', va='center',
                            fontsize=6, color='white', fontweight='bold', zorder=4,
                            path_effects=[pe.withStroke(linewidth=2, foreground='#0d1117')])

            # Path labels with arrows
            if len(path_a) >= 3:
                mid = len(path_a)//2
                mx, my = board_pos(*path_a[mid], n)
                ax.annotate(f'Caminho A\n({len(path_a)-1} arestas)',
                            (mx, my), textcoords="offset points",
                            xytext=(-35, 18), ha='center', fontsize=10,
                            color='#f97583', fontweight='bold',
                            bbox=dict(boxstyle='round,pad=0.3', fc='#0d1117', ec='#f97583', alpha=0.8))
            if len(path_b) >= 3:
                mid = len(path_b)//2
                mx, my = board_pos(*path_b[mid], n)
                ax.annotate(f'Caminho B\n({len(path_b)-1} arestas)',
                            (mx, my), textcoords="offset points",
                            xytext=(35, -18), ha='center', fontsize=10,
                            color='#79c0ff', fontweight='bold',
                            bbox=dict(boxstyle='round,pad=0.3', fc='#0d1117', ec='#79c0ff', alpha=0.8))

            # LCA label
            lx, ly = board_pos(*lca, n)
            ax.annotate('LCA', (lx, ly), textcoords="offset points",
                        xytext=(0, 18), ha='center', fontsize=11,
                        color='#f0883e', fontweight='bold',
                        bbox=dict(boxstyle='round,pad=0.2', fc='#0d1117', ec='#f0883e', alpha=0.8))

        # ─── Panel 3: The fundamental cycle ───
        elif ax_idx == 2:
            # Dim everything
            for u, v in tree_edges:
                x1,y1 = board_pos(*u, n)
                x2,y2 = board_pos(*v, n)
                ax.plot([x1,x2],[y1,y2], color='#3d4450', lw=0.4, alpha=0.12, zorder=0)

            # Draw cycle with direction arrows
            for i in range(len(cycle)):
                v1 = cycle[i]
                v2 = cycle[(i+1) % len(cycle)]
                x1,y1 = board_pos(*v1, n)
                x2,y2 = board_pos(*v2, n)

                e = tuple(sorted([v1, v2]))
                if e == back_edge:
                    color_ = '#d2a8ff'
                    ls_ = '--'
                else:
                    color_ = '#3fb950'
                    ls_ = '-'

                ax.plot([x1,x2],[y1,y2], color=color_, lw=4.5, ls=ls_,
                        alpha=0.9, zorder=2, solid_capstyle='round')

                # Direction arrow at midpoint
                mx, my = (x1+x2)/2, (y1+y2)/2
                dx, dy = (x2-x1)*0.12, (y2-y1)*0.12
                ax.annotate('', xy=(mx+dx, my+dy), xytext=(mx-dx, my-dy),
                            arrowprops=dict(arrowstyle='->', color=color_,
                                           lw=2.5, mutation_scale=18),
                            zorder=5)

            # Vertices
            for v in all_verts:
                x, y = board_pos(*v, n)
                in_cyc = v in set(cycle)
                if in_cyc:
                    ax.scatter(x, y, s=200, c='#3fb950', edgecolors='#ffffff',
                               linewidths=1.5, zorder=3)
                    label = f'({v[0]},{v[1]})'
                    ax.text(x, y, label, ha='center', va='center',
                            fontsize=6, color='white', fontweight='bold', zorder=4,
                            path_effects=[pe.withStroke(linewidth=2, foreground='#0d1117')])
                else:
                    ax.scatter(x, y, s=18, c='#181d25', edgecolors='#2a3040',
                               linewidths=0.5, zorder=1)

            # GF(2) annotation
            gf2_lines = [
                'Vetor indicador:  $\\mathbf{c}[e] = 1$ se $e \\in$ ciclo, $0$ c.c.',
                'Todo vértice tem grau par  →  $\\mathbf{c} \\in \\ker(\\partial_1)$'
            ]
            for idx_, line_ in enumerate(gf2_lines):
                ax.text(0.5, -0.06 - idx_*0.06, line_, transform=ax.transAxes,
                        ha='center', fontsize=9.5, color='#8b949e')

    # Legend
    legend_elements = [
        mpatches.Patch(facecolor='#f97583', label='Caminho A (árvore → LCA)'),
        mpatches.Patch(facecolor='#79c0ff', label='Caminho B (árvore → LCA)'),
        mpatches.Patch(facecolor='#d2a8ff', label='Back-edge (não-árvore)'),
        mpatches.Patch(facecolor='#3fb950', label='Ciclo fundamental GF(2)'),
        mpatches.Patch(facecolor='#f0883e', label='LCA (ancestral comum)'),
    ]
    fig.legend(handles=legend_elements, loc='lower center', ncol=5,
               fontsize=10, frameon=False, labelcolor='#c9d1d9')

    fig.suptitle('Ciclos Fundamentais GF(2) no Grafo do Cavalo  —  Tabuleiro 5×5',
                 color='#f0f6fc', fontsize=17, fontweight='bold', y=0.98)

    plt.tight_layout(rect=[0, 0.07, 1, 0.93])
    fig.savefig('gf2_cycle_dfs_tree.png', dpi=200, facecolor='#0d1117',
                bbox_inches='tight', pad_inches=0.3)
    plt.close()
    print("✓ Figura 1 salva: gf2_cycle_dfs_tree.png")


# ═══════════════════════════════════════════
# FIGURE 2: XOR of two cycles
# ═══════════════════════════════════════════
def gen_figure2():
    n = 5
    adj = knight_graph(n)

    # Find short cycles via BFS
    def find_short_cycles(adj, max_len=8, max_count=200):
        cycles_found = []
        seen_edge_sets = set()
        for start in sorted(adj.keys()):
            queue = deque([(start, [start])])
            while queue and len(cycles_found) < max_count:
                v, path = queue.popleft()
                if len(path) > max_len:
                    continue
                for w in adj[v]:
                    if w == start and len(path) >= 4:
                        edge_set = frozenset(
                            tuple(sorted([path[i], path[(i+1)%len(path)]]))
                            for i in range(len(path))
                        )
                        if edge_set not in seen_edge_sets:
                            seen_edge_sets.add(edge_set)
                            cycles_found.append(path[:])
                    elif w not in set(path) and len(path) < max_len:
                        queue.append((w, path + [w]))
        return cycles_found

    cycles = find_short_cycles(adj)
    print(f"  Encontrados {len(cycles)} ciclos curtos")

    # Find pair with 1–2 shared edges and nice visual result
    best = None
    for i in range(len(cycles)):
        ei = set(tuple(sorted([cycles[i][k], cycles[i][(k+1)%len(cycles[i])]])) 
                 for k in range(len(cycles[i])))
        for j in range(i+1, len(cycles)):
            ej = set(tuple(sorted([cycles[j][k], cycles[j][(k+1)%len(cycles[j])]])) 
                     for k in range(len(cycles[j])))
            shared = ei & ej
            xor = ei ^ ej
            if 1 <= len(shared) <= 2 and 5 <= len(xor) <= 10:
                # Check result is a single cycle (all degrees even)
                deg = defaultdict(int)
                for e in xor:
                    deg[e[0]] += 1; deg[e[1]] += 1
                is_cycle = all(d == 2 for d in deg.values())
                if is_cycle:
                    score = 100 - abs(len(ei) - len(ej)) - len(xor)
                    if best is None or score > best[0]:
                        best = (score, i, j, shared, xor, ei, ej)

    if best is None:
        # Accept any pair with shared edges
        for i in range(len(cycles)):
            ei = set(tuple(sorted([cycles[i][k], cycles[i][(k+1)%len(cycles[i])]])) 
                     for k in range(len(cycles[i])))
            for j in range(i+1, len(cycles)):
                ej = set(tuple(sorted([cycles[j][k], cycles[j][(k+1)%len(cycles[j])]])) 
                         for k in range(len(cycles[j])))
                shared = ei & ej
                xor = ei ^ ej
                if 1 <= len(shared) <= 3 and len(xor) >= 4:
                    best = (0, i, j, shared, xor, ei, ej)
                    break
            if best: break

    _, ci, cj, shared, xor_edges, edges_c1, edges_c2 = best
    c1, c2 = cycles[ci], cycles[cj]
    only_c1 = edges_c1 - edges_c2
    only_c2 = edges_c2 - edges_c1

    # Check if XOR result is a cycle
    deg = defaultdict(int)
    for e in xor_edges:
        deg[e[0]] += 1; deg[e[1]] += 1
    xor_is_cycle = all(d == 2 for d in deg.values())

    print(f"  Ciclo c₁: {c1}  ({len(edges_c1)} arestas)")
    print(f"  Ciclo c₂: {c2}  ({len(edges_c2)} arestas)")
    print(f"  Compartilhadas: {len(shared)}")
    print(f"  c₁ ⊕ c₂: {len(xor_edges)} arestas, {'ciclo' if xor_is_cycle else 'não-ciclo'}")

    # ── Plot ──
    fig, axes = plt.subplots(1, 3, figsize=(21, 7.5))
    fig.patch.set_facecolor('#0d1117')

    titles = [
        f'Ciclo $c_1$ — {len(edges_c1)} arestas',
        f'Ciclo $c_2$ — {len(edges_c2)} arestas\n({len(shared)} aresta{"s" if len(shared)>1 else ""} em comum)',
        f'$c_1 \\oplus c_2$  em GF(2) — {len(xor_edges)} arestas\n{"Ciclo" if xor_is_cycle else "Caminho"} resultante!'
    ]
    tcolors = ['#f97583', '#79c0ff', '#d2a8ff']

    all_graph_edges = set()
    for v in adj:
        for w in adj[v]:
            all_graph_edges.add(tuple(sorted([v,w])))

    for ax_idx, ax in enumerate(axes):
        ax.set_facecolor('#0d1117')
        ax.set_aspect('equal')
        ax.set_xlim(-0.7, (n-1)*1.1+0.7)
        ax.set_ylim(-0.7, (n-1)*1.1+0.7)
        ax.axis('off')
        ax.set_title(titles[ax_idx], color=tcolors[ax_idx],
                     fontsize=13, fontweight='bold', pad=14, linespacing=1.3)

        # Board
        for r in range(n):
            for c in range(n):
                x, y = board_pos(r, c, n)
                clr = '#1a2030' if (r+c)%2==0 else '#151b25'
                rect = plt.Rectangle((x-0.48, y-0.48), 0.96, 0.96,
                                     facecolor=clr, edgecolor='#2a3040',
                                     linewidth=0.5, zorder=0)
                ax.add_patch(rect)

        # Dim all knight edges
        for e in all_graph_edges:
            u, v = e
            x1,y1 = board_pos(*u, n)
            x2,y2 = board_pos(*v, n)
            ax.plot([x1,x2],[y1,y2], color='#3d4450', lw=0.25, alpha=0.08, zorder=0)

        if ax_idx == 0:  # ── c₁ ──
            for e in edges_c1:
                u, v = e
                x1,y1 = board_pos(*u, n); x2,y2 = board_pos(*v, n)
                if e in shared:
                    ax.plot([x1,x2],[y1,y2], color='#f0883e', lw=5, alpha=0.95,
                            zorder=2, solid_capstyle='round')
                    # Highlight label
                    mx, my = (x1+x2)/2, (y1+y2)/2
                    ax.plot(mx, my, 's', color='#f0883e', ms=6, zorder=5)
                else:
                    ax.plot([x1,x2],[y1,y2], color='#f97583', lw=4, alpha=0.9,
                            zorder=2, solid_capstyle='round')
            # Vertices
            for v in sorted(adj.keys()):
                x, y = board_pos(*v, n)
                if v in set(c1):
                    ax.scatter(x, y, s=160, c='#f97583', edgecolors='#ffffff',
                               linewidths=1.2, zorder=3)
                    ax.text(x, y, f'({v[0]},{v[1]})', ha='center', va='center',
                            fontsize=5.5, color='white', fontweight='bold', zorder=4,
                            path_effects=[pe.withStroke(linewidth=1.5, foreground='#0d1117')])
                else:
                    ax.scatter(x, y, s=15, c='#181d25', edgecolors='#2a3040',
                               linewidths=0.5, zorder=1)

        elif ax_idx == 1:  # ── c₂ ──
            for e in edges_c2:
                u, v = e
                x1,y1 = board_pos(*u, n); x2,y2 = board_pos(*v, n)
                if e in shared:
                    ax.plot([x1,x2],[y1,y2], color='#f0883e', lw=5, alpha=0.95,
                            zorder=2, solid_capstyle='round')
                    mx, my = (x1+x2)/2, (y1+y2)/2
                    ax.plot(mx, my, 's', color='#f0883e', ms=6, zorder=5)
                else:
                    ax.plot([x1,x2],[y1,y2], color='#79c0ff', lw=4, alpha=0.9,
                            zorder=2, solid_capstyle='round')
            for v in sorted(adj.keys()):
                x, y = board_pos(*v, n)
                if v in set(c2):
                    ax.scatter(x, y, s=160, c='#79c0ff', edgecolors='#ffffff',
                               linewidths=1.2, zorder=3)
                    ax.text(x, y, f'({v[0]},{v[1]})', ha='center', va='center',
                            fontsize=5.5, color='white', fontweight='bold', zorder=4,
                            path_effects=[pe.withStroke(linewidth=1.5, foreground='#0d1117')])
                else:
                    ax.scatter(x, y, s=15, c='#181d25', edgecolors='#2a3040',
                               linewidths=0.5, zorder=1)

        elif ax_idx == 2:  # ── XOR ──
            # Draw cancelled edges (ghosted + X)
            for e in shared:
                u, v = e
                x1,y1 = board_pos(*u, n); x2,y2 = board_pos(*v, n)
                ax.plot([x1,x2],[y1,y2], color='#f0883e', lw=3, alpha=0.25,
                        ls=':', zorder=1)
                mx, my = (x1+x2)/2, (y1+y2)/2
                ax.plot(mx, my, 'X', color='#f0883e', ms=16, markeredgewidth=3, zorder=6)

            # Draw XOR result edges
            for e in xor_edges:
                u, v = e
                x1,y1 = board_pos(*u, n); x2,y2 = board_pos(*v, n)
                # Color by origin
                if e in only_c1:
                    ax.plot([x1,x2],[y1,y2], color='#d2a8ff', lw=4.5, alpha=0.9,
                            zorder=2, solid_capstyle='round')
                else:
                    ax.plot([x1,x2],[y1,y2], color='#bc8cff', lw=4.5, alpha=0.9,
                            zorder=2, solid_capstyle='round')

            # Vertices
            xor_verts = set()
            for e in xor_edges:
                xor_verts.add(e[0]); xor_verts.add(e[1])

            for v in sorted(adj.keys()):
                x, y = board_pos(*v, n)
                if v in xor_verts:
                    ax.scatter(x, y, s=180, c='#d2a8ff', edgecolors='#ffffff',
                               linewidths=1.2, zorder=3)
                    ax.text(x, y, f'({v[0]},{v[1]})', ha='center', va='center',
                            fontsize=5.5, color='white', fontweight='bold', zorder=4,
                            path_effects=[pe.withStroke(linewidth=1.5, foreground='#0d1117')])
                else:
                    ax.scatter(x, y, s=15, c='#181d25', edgecolors='#2a3040',
                               linewidths=0.5, zorder=1)

            # Bottom annotations
            ax.text(0.5, -0.06,
                    '$1 \\oplus 1 = 0$ em $\\mathbb{F}_2$  →  arestas em comum cancelam',
                    transform=ax.transAxes, ha='center', fontsize=10, color='#8b949e')
            ax.text(0.5, -0.12,
                    '$c_1 \\oplus c_2 \\in \\ker(\\partial_1)$  (todo vértice mantém grau par)',
                    transform=ax.transAxes, ha='center', fontsize=10, color='#3fb950')

    # Legend
    legend_elements = [
        mpatches.Patch(facecolor='#f97583', label='Arestas de $c_1$ apenas'),
        mpatches.Patch(facecolor='#79c0ff', label='Arestas de $c_2$ apenas'),
        mpatches.Patch(facecolor='#f0883e', label='Arestas compartilhadas (cancelam em GF(2))'),
        mpatches.Patch(facecolor='#d2a8ff', label='$c_1 \\oplus c_2$ — resultado'),
    ]
    fig.legend(handles=legend_elements, loc='lower center', ncol=4,
               fontsize=10, frameon=False, labelcolor='#c9d1d9')

    fig.suptitle('Operação XOR de Ciclos em GF(2) no Grafo do Cavalo',
                 color='#f0f6fc', fontsize=17, fontweight='bold', y=0.98)

    plt.tight_layout(rect=[0, 0.07, 1, 0.93])
    fig.savefig('gf2_xor_cycles.png', dpi=200, facecolor='#0d1117',
                bbox_inches='tight', pad_inches=0.3)
    plt.close()
    print("✓ Figura 2 salva: gf2_xor_cycles.png")


if __name__ == '__main__':
    print("Gerando figuras GF(2)...\n")
    print("─── Figura 1: Árvore DFS + Ciclo Fundamental ───")
    gen_figure1()
    print()
    print("─── Figura 2: XOR de Ciclos ───")
    gen_figure2()
    print("\nConcluído!")
