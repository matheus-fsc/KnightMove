#!/usr/bin/env python3
"""
Generate two paper-quality figures (white background) for the GF(2) section.

Figure 1 (gf2_xor_cycles.png):
  3 panels (large): c₁, c₂, c₁⊕c₂ on the 6×6 board

Figure 2 (gf2_cycle_dfs_tree.png):
  2×2 grid:
    (0,0) Binary tree showing two branches converging → cycle
    (0,1) That cycle inscribed on the board
    (1,0) Second cycle on its own board
    (1,1) Third cycle on its own board
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
import numpy as np
from collections import defaultdict, deque
from matplotlib.lines import Line2D

OFFSETS = [(1,2),(2,1),(2,-1),(1,-2),(-1,-2),(-2,-1),(-2,1),(-1,2)]
N = 6

BOARD_LIGHT = '#f0f0f0'
BOARD_DARK  = '#dcdcdc'
BOARD_EDGE  = '#c0c0c0'
BG          = '#ffffff'
TEXT_DARK   = '#222222'
TEXT_MID    = '#555555'
NODE_IDLE   = '#d0d0d0'

CYCLE_COLORS = ['#d62728', '#1f77b4', '#2ca02c']
C1_COLOR  = '#d62728'
C2_COLOR  = '#1f77b4'
XOR_COLOR = '#7b2d8e'
SHARED_COLOR = '#ff7f0e'
TREE_COLOR = '#bbbbbb'
PATH_A_COLOR = '#d62728'
PATH_B_COLOR = '#1f77b4'
LCA_COLOR    = '#ff7f0e'
COLLISION_COLOR = '#7b2d8e'


def knight_graph(n):
    adj = defaultdict(list)
    for r in range(n):
        for c in range(n):
            for dr, dc in OFFSETS:
                nr, nc = r+dr, c+dc
                if 0 <= nr < n and 0 <= nc < n:
                    adj[(r,c)].append((nr,nc))
    return adj


def board_pos(r, c, n=N):
    return (c, n - 1 - r)


def draw_board(ax, n=N):
    for r in range(n):
        for c in range(n):
            x, y = board_pos(r, c, n)
            clr = BOARD_LIGHT if (r+c) % 2 == 0 else BOARD_DARK
            rect = plt.Rectangle((x - 0.45, y - 0.45), 0.9, 0.9,
                                 facecolor=clr, edgecolor=BOARD_EDGE,
                                 linewidth=0.6, zorder=0)
            ax.add_patch(rect)
    ax.set_xlim(-0.6, n - 0.4)
    ax.set_ylim(-0.6, n - 0.4)
    ax.set_aspect('equal')
    ax.axis('off')


def draw_idle_nodes(ax, adj, active_verts=set(), n=N):
    for v in sorted(adj.keys()):
        if v not in active_verts:
            x, y = board_pos(*v, n)
            ax.plot(x, y, 'o', color=NODE_IDLE, ms=4, mec=BOARD_EDGE,
                    mew=0.4, zorder=1)


def draw_edge(ax, u, v, color, lw=2.5, ls='-', alpha=0.9, zorder=2, n=N):
    x1, y1 = board_pos(*u, n)
    x2, y2 = board_pos(*v, n)
    ax.plot([x1, x2], [y1, y2], color=color, lw=lw, ls=ls,
            alpha=alpha, zorder=zorder, solid_capstyle='round')


def draw_node(ax, v, color, sz=7, ec='white', ew=1.2, zorder=4, n=N):
    x, y = board_pos(*v, n)
    ax.plot(x, y, 'o', color=color, ms=sz, mec=ec, mew=ew, zorder=zorder)


def label_node(ax, v, text, color=TEXT_DARK, n=N, fontsize=7):
    x, y = board_pos(*v, n)
    ax.text(x, y, text, ha='center', va='center', fontsize=fontsize,
            color='white', fontweight='bold', zorder=5,
            path_effects=[pe.withStroke(linewidth=2, foreground=color)])


def draw_cycle_on_board(ax, adj, cycle, color, lw=3, node_sz=7):
    verts = set(cycle)
    draw_idle_nodes(ax, adj, verts)
    for i in range(len(cycle)):
        u = cycle[i]
        v = cycle[(i + 1) % len(cycle)]
        draw_edge(ax, u, v, color, lw=lw)
    for v in cycle:
        draw_node(ax, v, color, sz=node_sz, ec='white')
    for v in cycle:
        x, y = board_pos(*v)
        ax.text(x, y, f'{v[0]},{v[1]}', ha='center', va='center',
                fontsize=5.5, color='white', fontweight='bold', zorder=5,
                path_effects=[pe.withStroke(linewidth=1.5, foreground=color)])


# ═══════════════════════════════════════════
# FIGURE 2: Binary tree + cycles on boards
# ═══════════════════════════════════════════
def gen_figure_dfs():
    adj = knight_graph(N)

    # Deterministic cycles for Figure 2 (dfs_tree) on 6x6 board
    tree_cycle = [(0,0), (1,2), (2,4), (4,5), (3,3), (2,1)]
    tree_root = (0,0)
    tree_collision = (4,5)
    tree_pa = [(0,0), (1,2), (2,4), (4,5)]
    tree_pb = [(0,0), (2,1), (3,3), (4,5)]

    # Other fundamental cycles: one of length 4, one of length 6
    cycle_4 = [(0,1), (2,2), (4,1), (2,0)]
    cycle_6_other = [(0,2), (1,4), (2,2), (4,3), (3,1), (1,0)]

    others = [cycle_4, cycle_6_other]
    all_cycles = [tree_cycle] + others

    print(f"  Ciclo da árvore: {len(tree_cycle)} vértices, "
          f"pa={len(tree_pa)-1} arestas, pb={len(tree_pb)-1} arestas")
    print(f"    Raiz={tree_root}, Colisão={tree_collision}")
    for i, c in enumerate(others):
        print(f"  Ciclo extra {i+1}: {len(c)} vértices")

    # ── Layout: 2×2 ──
    fig, axes = plt.subplots(2, 2, figsize=(12, 11))
    fig.patch.set_facecolor(BG)

    # ═══ Panel (0,0): Binary tree ═══
    ax = axes[0, 0]
    ax.set_facecolor(BG)
    ax.axis('off')

    lca = tree_root
    pa = tree_pa
    pb = tree_pb
    collision = tree_collision
    k_a = len(pa) - 1
    k_b = len(pb) - 1
    max_k = max(k_a, k_b)

    spread = 1.5
    y_top = max_k
    tree_positions = {}

    tree_positions[('root', 0)] = (0, y_top)
    for i in range(1, k_a):
        tree_positions[('a', i)] = (-spread, y_top - i * (y_top / k_a))
    for i in range(1, k_b):
        tree_positions[('b', i)] = (spread, y_top - i * (y_top / k_b))
    tree_positions[('collision', 0)] = (0, 0)

    def draw_tree_edge(ax, pos1, pos2, color, lw=3):
        ax.plot([pos1[0], pos2[0]], [pos1[1], pos2[1]],
                color=color, lw=lw, alpha=0.85, zorder=1, solid_capstyle='round')
        mx = (pos1[0] + pos2[0]) / 2
        my = (pos1[1] + pos2[1]) / 2
        dx = (pos2[0] - pos1[0]) * 0.08
        dy = (pos2[1] - pos1[1]) * 0.08
        ax.annotate('', xy=(mx+dx, my+dy), xytext=(mx-dx, my-dy),
                    arrowprops=dict(arrowstyle='->', color=color, lw=2,
                                   mutation_scale=14), zorder=2)

    def draw_tree_node(ax, pos, label, color, sz=350):
        ax.scatter(*pos, s=sz, c=color, edgecolors='white', linewidths=1.5, zorder=3)
        ax.text(pos[0], pos[1], label, ha='center', va='center',
                fontsize=8, color='white', fontweight='bold', zorder=4,
                path_effects=[pe.withStroke(linewidth=2, foreground=color)])

    prev_a = tree_positions[('root', 0)]
    for i in range(1, k_a):
        cur = tree_positions[('a', i)]
        draw_tree_edge(ax, prev_a, cur, PATH_A_COLOR)
        prev_a = cur
    draw_tree_edge(ax, prev_a, tree_positions[('collision', 0)], PATH_A_COLOR)

    prev_b = tree_positions[('root', 0)]
    for i in range(1, k_b):
        cur = tree_positions[('b', i)]
        draw_tree_edge(ax, prev_b, cur, PATH_B_COLOR)
        prev_b = cur
    draw_tree_edge(ax, prev_b, tree_positions[('collision', 0)], PATH_B_COLOR)

    draw_tree_node(ax, tree_positions[('root', 0)],
                   f'{lca[0]},{lca[1]}', LCA_COLOR, sz=400)
    ax.text(tree_positions[('root', 0)][0], tree_positions[('root', 0)][1] + 0.45,
            'Raiz (LCA)', ha='center', fontsize=9, color=LCA_COLOR, fontweight='bold')

    for i in range(1, k_a):
        v = pa[i]
        draw_tree_node(ax, tree_positions[('a', i)],
                       f'{v[0]},{v[1]}', PATH_A_COLOR, sz=320)
    for i in range(1, k_b):
        v = pb[i]
        draw_tree_node(ax, tree_positions[('b', i)],
                       f'{v[0]},{v[1]}', PATH_B_COLOR, sz=320)

    draw_tree_node(ax, tree_positions[('collision', 0)],
                   f'{collision[0]},{collision[1]}', COLLISION_COLOR, sz=400)
    ax.text(tree_positions[('collision', 0)][0],
            tree_positions[('collision', 0)][1] - 0.45,
            'Colisão (back-edge)', ha='center', fontsize=9,
            color=COLLISION_COLOR, fontweight='bold')

    ax.text(-spread - 0.3, y_top / 2, 'Caminho A', ha='right', va='center',
            fontsize=10, color=PATH_A_COLOR, fontweight='bold', rotation=90)
    ax.text(spread + 0.3, y_top / 2, 'Caminho B', ha='left', va='center',
            fontsize=10, color=PATH_B_COLOR, fontweight='bold', rotation=-90)

    ax.set_xlim(-spread - 1.0, spread + 1.0)
    ax.set_ylim(-0.8, y_top + 0.8)
    ax.set_aspect('equal')
    ax.set_title('Expansão em árvore binária\n(dois ramos → mesma casa = ciclo)',
                 fontsize=11, color=TEXT_DARK, fontweight='bold', pad=10,
                 linespacing=1.3)

    # ═══ Panels with cycles on boards ═══
    board_panels = [(0, 1), (1, 0), (1, 1)]
    panel_labels = [
        f'Ciclo da árvore — {len(all_cycles[0])} arestas',
        f'Outro ciclo fundamental — {len(all_cycles[1])} arestas',
        f'Outro ciclo fundamental — {len(all_cycles[2])} arestas',
    ]

    for idx, (row, col) in enumerate(board_panels):
        ax = axes[row, col]
        ax.set_facecolor(BG)
        draw_board(ax)

        cycle = all_cycles[idx]
        color = CYCLE_COLORS[idx % len(CYCLE_COLORS)]

        draw_cycle_on_board(ax, adj, cycle, color)
        ax.set_title(panel_labels[idx], fontsize=11, color=color,
                     fontweight='bold', pad=10)

    legend_elements = [
        Line2D([0], [0], color=LCA_COLOR, lw=0, marker='o', ms=8,
               mec='white', label='Raiz / LCA'),
        Line2D([0], [0], color=PATH_A_COLOR, lw=2.5, label='Caminho A'),
        Line2D([0], [0], color=PATH_B_COLOR, lw=2.5, label='Caminho B'),
        Line2D([0], [0], color=COLLISION_COLOR, lw=0, marker='o', ms=8,
               mec='white', label='Colisão (back-edge)'),
    ]
    fig.legend(handles=legend_elements, loc='lower center', ncol=4,
               fontsize=9, frameon=True, fancybox=True,
               edgecolor='#cccccc', facecolor='#fafafa',
               labelcolor=TEXT_DARK)

    fig.suptitle('Ciclos Fundamentais GF(2) no Grafo do Cavalo — Tabuleiro 6×6',
                 color=TEXT_DARK, fontsize=14, fontweight='bold', y=0.98)

    plt.tight_layout(rect=[0, 0.05, 1, 0.95])
    fig.savefig('gf2_cycle_dfs_tree.png', dpi=200, facecolor=BG,
                bbox_inches='tight', pad_inches=0.2)
    plt.close()
    print("  ✓ gf2_cycle_dfs_tree.png salva")


# ═══════════════════════════════════════════
# FIGURE 1: XOR of two cycles (larger)
# ═══════════════════════════════════════════
def gen_figure_xor():
    adj = knight_graph(N)

    def find_short_cycles(max_len=8, max_count=300):
        cycles_found = []
        seen = set()
        for start in sorted(adj.keys()):
            queue = deque([(start, [start])])
            while queue and len(cycles_found) < max_count:
                v, path = queue.popleft()
                if len(path) > max_len:
                    continue
                for w in adj[v]:
                    if w == start and len(path) >= 4:
                        es = frozenset(
                            tuple(sorted([path[i], path[(i+1) % len(path)]]))
                            for i in range(len(path))
                        )
                        if es not in seen:
                            seen.add(es)
                            cycles_found.append(path[:])
                    elif w not in set(path) and len(path) < max_len:
                        queue.append((w, path + [w]))
        return cycles_found

    cycles = find_short_cycles()
    print(f"  Encontrados {len(cycles)} ciclos curtos")

    best = None
    for i in range(len(cycles)):
        ei = set(tuple(sorted([cycles[i][k], cycles[i][(k+1) % len(cycles[i])]]))
                 for k in range(len(cycles[i])))
        for j in range(i+1, len(cycles)):
            ej = set(tuple(sorted([cycles[j][k], cycles[j][(k+1) % len(cycles[j])]]))
                     for k in range(len(cycles[j])))
            shared = ei & ej
            xor = ei ^ ej
            if 1 <= len(shared) <= 2 and 5 <= len(xor) <= 10:
                deg = defaultdict(int)
                for e in xor:
                    deg[e[0]] += 1
                    deg[e[1]] += 1
                if all(d == 2 for d in deg.values()):
                    score = 100 - abs(len(ei) - len(ej)) - len(xor)
                    if best is None or score > best[0]:
                        best = (score, i, j, shared, xor, ei, ej)

    if best is None:
        for i in range(len(cycles)):
            ei = set(tuple(sorted([cycles[i][k], cycles[i][(k+1) % len(cycles[i])]]))
                     for k in range(len(cycles[i])))
            for j in range(i+1, len(cycles)):
                ej = set(tuple(sorted([cycles[j][k], cycles[j][(k+1) % len(cycles[j])]]))
                         for k in range(len(cycles[j])))
                shared = ei & ej
                xor = ei ^ ej
                if 1 <= len(shared) <= 3 and len(xor) >= 4:
                    best = (0, i, j, shared, xor, ei, ej)
                    break
            if best:
                break

    _, ci, cj, shared, xor_edges, edges_c1, edges_c2 = best
    c1, c2 = cycles[ci], cycles[cj]

    print(f"  c₁: {c1} ({len(edges_c1)} arestas)")
    print(f"  c₂: {c2} ({len(edges_c2)} arestas)")
    print(f"  Compartilhadas: {len(shared)}")
    print(f"  c₁⊕c₂: {len(xor_edges)} arestas")

    fig, axes = plt.subplots(1, 3, figsize=(16, 8))
    fig.patch.set_facecolor(BG)

    titles = [
        f'Ciclo $c_1$ — {len(edges_c1)} arestas',
        f'Ciclo $c_2$ — {len(edges_c2)} arestas\n({len(shared)} em comum)',
        f'$c_1 \\oplus c_2$ — {len(xor_edges)} arestas\n(arestas comuns cancelam)',
    ]
    tcolors = [C1_COLOR, C2_COLOR, XOR_COLOR]

    for ax_idx, ax in enumerate(axes):
        ax.set_facecolor(BG)
        draw_board(ax)

        verts_here = set()

        if ax_idx == 0:
            for e in edges_c1:
                u, v = e
                if e in shared:
                    draw_edge(ax, u, v, SHARED_COLOR, lw=4.5, alpha=0.9, zorder=2)
                else:
                    draw_edge(ax, u, v, C1_COLOR, lw=3.5, alpha=0.9, zorder=2)
            verts_here = set(c1)
            for v in verts_here:
                draw_node(ax, v, C1_COLOR, sz=8, ec='white')

        elif ax_idx == 1:
            for e in edges_c2:
                u, v = e
                if e in shared:
                    draw_edge(ax, u, v, SHARED_COLOR, lw=4.5, alpha=0.9, zorder=2)
                    x1, y1 = board_pos(*u)
                    x2, y2 = board_pos(*v)
                    mx, my = (x1+x2)/2, (y1+y2)/2
                    ax.plot(mx, my, 's', color=SHARED_COLOR, ms=6, zorder=5)
                else:
                    draw_edge(ax, u, v, C2_COLOR, lw=3.5, alpha=0.9, zorder=2)
            verts_here = set(c2)
            for v in verts_here:
                draw_node(ax, v, C2_COLOR, sz=8, ec='white')

        elif ax_idx == 2:
            for e in shared:
                u, v = e
                draw_edge(ax, u, v, SHARED_COLOR, lw=3, ls=':', alpha=0.3, zorder=1)
                x1, y1 = board_pos(*u)
                x2, y2 = board_pos(*v)
                mx, my = (x1+x2)/2, (y1+y2)/2
                ax.plot(mx, my, 'X', color=SHARED_COLOR, ms=14,
                        markeredgewidth=2.5, zorder=6)

            for e in xor_edges:
                u, v = e
                draw_edge(ax, u, v, XOR_COLOR, lw=4, alpha=0.9, zorder=2)

            xor_verts = set()
            for e in xor_edges:
                xor_verts.add(e[0])
                xor_verts.add(e[1])
            verts_here = xor_verts
            for v in verts_here:
                draw_node(ax, v, XOR_COLOR, sz=8, ec='white')

        draw_idle_nodes(ax, adj, verts_here)

        for v in verts_here:
            x, y = board_pos(*v)
            ax.text(x, y, f'{v[0]},{v[1]}', ha='center', va='center',
                    fontsize=6, color='white', fontweight='bold', zorder=5,
                    path_effects=[pe.withStroke(linewidth=1.5,
                                               foreground=tcolors[ax_idx])])

        ax.set_title(titles[ax_idx], color=tcolors[ax_idx],
                     fontsize=13, fontweight='bold', pad=12, linespacing=1.3)

    axes[2].text(0.5, -0.07,
                 r'$1 \oplus 1 = 0$ em $\mathbb{F}_2$'
                 r' $\rightarrow$ arestas comuns cancelam',
                 transform=axes[2].transAxes, ha='center', fontsize=10,
                 color=TEXT_MID)
    axes[2].text(0.5, -0.13,
                 r'$c_1 \oplus c_2 \in \ker(\partial_1)$'
                 r' (grau par preservado)',
                 transform=axes[2].transAxes, ha='center', fontsize=10,
                 color=XOR_COLOR, fontweight='bold')

    legend_elements = [
        mpatches.Patch(facecolor=C1_COLOR, label='Arestas de $c_1$'),
        mpatches.Patch(facecolor=C2_COLOR, label='Arestas de $c_2$'),
        mpatches.Patch(facecolor=SHARED_COLOR, label='Compartilhadas (cancelam)'),
        mpatches.Patch(facecolor=XOR_COLOR, label='$c_1 \\oplus c_2$'),
    ]
    fig.legend(handles=legend_elements, loc='lower center', ncol=4,
               fontsize=10, frameon=True, fancybox=True,
               edgecolor='#cccccc', facecolor='#fafafa',
               labelcolor=TEXT_DARK)

    fig.suptitle('Operação XOR em GF(2) — Tabuleiro 6×6',
                 color=TEXT_DARK, fontsize=15, fontweight='bold', y=0.98)

    plt.tight_layout(rect=[0, 0.08, 1, 0.93])
    fig.savefig('gf2_xor_cycles.png', dpi=200, facecolor=BG,
                bbox_inches='tight', pad_inches=0.2)
    plt.close()
    print("  ✓ gf2_xor_cycles.png salva")


if __name__ == '__main__':
    print("Gerando figuras GF(2) (fundo branco)...\n")
    print("─── Figura XOR ───")
    gen_figure_xor()
    print()
    print("─── Figura DFS / Árvore Binária ───")
    gen_figure_dfs()
    print("\nConcluído!")
