"""print_backtrack_tree.py — imprime a ÁRVORE DE BACKTRACKING de knight_tours.py.

Reimplementa `_backtrack` com instrumentação (não altera o módulo original):
cada nó da árvore = uma decisão (aresta, valor) e o status devolvido pela
propagação R2 + Union-Find.

Uso:
    python print_backtrack_tree.py                 # 6x6, K=1, seed=0
    python print_backtrack_tree.py -n 6 -K 3 --seed 0
    python print_backtrack_tree.py --dot arvore.dot --json arvore.json
    python print_backtrack_tree.py --max-nodes 500

Legenda dos status (folhas/ramos):
    OK        propagação consistente, desce mais um nível
    CONTRA    contradição (grau > 2 ou vértice sem 2 arestas disponíveis)
    SUBTOUR   fechou um ciclo curto (Union-Find detectou componente < n²)
    TOUR      fechou o ciclo hamiltoniano  ✔
    DEADEND   sem arestas livres mas não é tour
"""

# --- bootstrap: torna core/ importavel a partir de qualquer subpasta ---
import sys as _sys, pathlib as _pathlib
for _anc in _pathlib.Path(__file__).resolve().parents:
    if (_anc / "core").is_dir() and (_anc / "pyproject.toml").exists():
        _sys.path.insert(0, str(_anc / "core"))
        break
# --- fim do bootstrap ---

import argparse
import json
import sys

import numpy as np

import knight_tours as kt

STATUS_NAME = {kt.OK: 'OK', kt.CONTRADICTION: 'CONTRA',
               kt.SUBTOUR: 'SUBTOUR', kt.COMPLETE_TOUR: 'TOUR'}


def sq(v, n):
    """Nome algébrico da casa (a1..f6) para v = row*n + col."""
    r, c = divmod(int(v), n)
    return f"{chr(ord('a') + c)}{r + 1}"


class Tree:
    def __init__(self, n, max_nodes):
        self.n = n
        self.nodes = []          # dicts
        self.children = {}       # id -> [ids]
        self.max_nodes = max_nodes
        self.truncated = False

    def add(self, parent, edge_uv, val, status, depth, n_free):
        if len(self.nodes) >= self.max_nodes:
            self.truncated = True
            return None
        nid = len(self.nodes)
        self.nodes.append({
            'id': nid, 'parent': parent, 'depth': depth,
            'u': edge_uv[0], 'v': edge_uv[1], 'val': val,
            'status': status, 'n_free': n_free,
            'subtree': 0, 'found_tour': False,
        })
        self.children.setdefault(parent, []).append(nid)
        self.children.setdefault(nid, [])
        return nid

    def finish(self):
        """Calcula tamanho de subárvore e propaga 'found_tour' para a raiz."""
        for nid in range(len(self.nodes) - 1, -1, -1):
            nd = self.nodes[nid]
            tot = 1
            for c in self.children.get(nid, ()):
                tot += self.nodes[c]['subtree']
                if self.nodes[c]['found_tour']:
                    nd['found_tour'] = True
            if nd['status'] == 'TOUR':
                nd['found_tour'] = True
            nd['subtree'] = tot


# --------------------------------------------------------------------------
# Backtracking instrumentado (espelha kt._backtrack)
# --------------------------------------------------------------------------

def backtrack_traced(state, ctx, rng, tours, K, tree, parent, depth):
    if len(tours) >= K:
        return True

    e, first_val = kt._choose_next_edge(state, ctx, rng)

    if e == -1:
        if kt._is_complete_tour(state, ctx):
            tours.append(kt._extract_tour(state, ctx))
            if parent >= 0:
                tree.nodes[parent]['status'] = 'TOUR'
        elif parent >= 0:
            tree.nodes[parent]['status'] = 'DEADEND'
        return len(tours) >= K

    u = int(ctx['edge_endpoints'][e, 0])
    v = int(ctx['edge_endpoints'][e, 1])

    for val in (first_val, 1 - first_val):
        snap = state.snapshot()
        status = kt.fix_and_propagate(state, ctx, e, val)
        nid = tree.add(parent, (u, v), val, STATUS_NAME[status],
                       depth, int(state.n_free))
        if status == kt.OK:
            if nid is not None:
                if backtrack_traced(state, ctx, rng, tours, K, tree, nid,
                                    depth + 1):
                    state.restore(snap)
                    return True
            else:
                if kt._backtrack(state, ctx, rng, tours, K):
                    state.restore(snap)
                    return True
        elif status == kt.COMPLETE_TOUR:
            if kt._is_complete_tour(state, ctx):
                tours.append(kt._extract_tour(state, ctx))
            if len(tours) >= K:
                state.restore(snap)
                return True
        state.restore(snap)

    return len(tours) >= K


def trace(n, K, seed, max_nodes):
    ctx = kt.build_graph(n)
    rng = np.random.default_rng(seed)
    state = kt.State(ctx['V'], ctx['E'])
    tree = Tree(n, max_nodes)

    n_free_before = state.n_free
    status = kt._propagate_initial(state, ctx)
    forced = n_free_before - state.n_free

    tours = []
    if status not in (kt.CONTRADICTION, kt.SUBTOUR):
        if status == kt.COMPLETE_TOUR and kt._is_complete_tour(state, ctx):
            tours.append(kt._extract_tour(state, ctx))
        else:
            backtrack_traced(state, ctx, rng, tours, K, tree, -1, 0)
    tree.finish()
    return tree, tours, ctx, forced, STATUS_NAME.get(status, 'OK')


# --------------------------------------------------------------------------
# Impressão ASCII
# --------------------------------------------------------------------------

MARK = {'TOUR': '✔ TOUR', 'CONTRA': '✘ contradição',
        'SUBTOUR': '✘ subciclo', 'DEADEND': '✘ beco',
        'OK': ''}


def print_ascii(tree, out=sys.stdout):
    n = tree.n

    def render(nid, prefix, is_last):
        nd = tree.nodes[nid]
        conn = '└── ' if is_last else '├── '
        rel = '=1' if nd['val'] == 1 else '=0'
        lbl = f"{sq(nd['u'], n)}-{sq(nd['v'], n)}{rel}"
        tag = MARK[nd['status']] if nd['status'] in MARK else nd['status']
        star = ' ★' if nd['found_tour'] and nd['status'] == 'OK' else ''
        extra = f"livres={nd['n_free']}"
        line = f"{prefix}{conn}[{nd['id']}] {lbl}  {extra}"
        if tag:
            line += f"  {tag}"
        line += star
        print(line, file=out)
        kids = tree.children.get(nid, [])
        newpref = prefix + ('    ' if is_last else '│   ')
        for i, c in enumerate(kids):
            render(c, newpref, i == len(kids) - 1)

    roots = tree.children.get(-1, [])
    print('raiz (após propagação R2 inicial)', file=out)
    for i, r in enumerate(roots):
        render(r, '', i == len(roots) - 1)


def to_dot(tree):
    n = tree.n
    color = {'TOUR': 'green', 'CONTRA': 'red', 'SUBTOUR': 'orange',
             'DEADEND': 'gray', 'OK': 'black'}
    L = ['digraph backtrack {', '  node [shape=box, fontname="monospace"];',
         '  root [label="raiz", shape=ellipse];']
    for nd in tree.nodes:
        lbl = (f"{sq(nd['u'], n)}-{sq(nd['v'], n)}="
               f"{nd['val']}\\n{nd['status']}")
        L.append(f'  n{nd["id"]} [label="{lbl}", '
                 f'color={color.get(nd["status"], "black")}];')
    for nd in tree.nodes:
        p = 'root' if nd['parent'] < 0 else f'n{nd["parent"]}'
        L.append(f'  {p} -> n{nd["id"]};')
    L.append('}')
    return '\n'.join(L)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('-n', type=int, default=6, help='lado do tabuleiro')
    ap.add_argument('-K', type=int, default=1, help='tours a encontrar')
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--max-nodes', type=int, default=20000)
    ap.add_argument('--dot', help='grava a árvore em formato Graphviz')
    ap.add_argument('--json', help='grava a árvore em JSON')
    ap.add_argument('--quiet', action='store_true', help='só estatísticas')
    args = ap.parse_args()

    tree, tours, ctx, forced, init_status = trace(
        args.n, args.K, args.seed, args.max_nodes)

    print(f"tabuleiro {args.n}x{args.n}: V={ctx['V']} E={ctx['E']}")
    print(f"propagação R2 inicial: {forced} arestas fixadas "
          f"(status {init_status})")
    print(f"nós na árvore: {len(tree.nodes)}"
          + ("  [TRUNCADA]" if tree.truncated else ""))
    from collections import Counter
    cnt = Counter(nd['status'] for nd in tree.nodes)
    print("folhas/ramos por status:", dict(cnt))
    print(f"tours encontrados: {len(tours)}"
          f"  (todos válidos: {all(kt.verify_tour(t, args.n) for t in tours)})")
    print()

    if not args.quiet:
        print_ascii(tree)

    if args.dot:
        with open(args.dot, 'w') as f:
            f.write(to_dot(tree))
        print(f"\n[dot] {args.dot}  (render: dot -Tpng {args.dot} -o arvore.png)")
    if args.json:
        with open(args.json, 'w') as f:
            json.dump({'n': args.n, 'K': args.K, 'seed': args.seed,
                       'nodes': tree.nodes,
                       'children': {str(k): v for k, v in tree.children.items()}},
                      f, indent=1)
        print(f"[json] {args.json}")


if __name__ == '__main__':
    main()
