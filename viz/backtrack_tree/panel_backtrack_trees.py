"""panel_backtrack_trees.py — painel comparativo de árvores de backtracking.

Roda a busca de knight_tours.py com várias seeds (cada uma produz uma raiz de
ramificação diferente) e gera um painel HTML com as árvores lado a lado,
para comparar estados de busca: corredor R2, leque de ramificação, subciclos.

Uso:
    python panel_backtrack_trees.py                        # seeds 8,0,7,10
    python panel_backtrack_trees.py -K 3 --seeds 0 1 2 3
    python panel_backtrack_trees.py -n 8 -K 2 --seeds 0 1
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
import os
from collections import Counter

import knight_tours as kt
import print_backtrack_tree as ptree
import render_backtrack_tree as rtree


def stats(tree, tours, n):
    c = Counter(nd['status'] for nd in tree.nodes)
    depths = [nd['depth'] for nd in tree.nodes]
    first_branch = next((nd['depth'] for nd in tree.nodes
                         if len(tree.children.get(nd['id'], [])) > 1), -1)
    return {
        'nos': len(tree.nodes),
        'ok': c['OK'], 'subtour': c['SUBTOUR'],
        'contra': c['CONTRA'], 'deadend': c['DEADEND'], 'tour': c['TOUR'],
        'prof': max(depths) if depths else 0,
        'primeira_ramif': first_branch,
        'corredor': first_branch,
        'nos_por_tour': round(len(tree.nodes) / max(len(tours), 1), 1),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('-n', type=int, default=6)
    ap.add_argument('-K', type=int, default=3)
    ap.add_argument('--seeds', type=int, nargs='+', default=[8, 0, 7, 10])
    ap.add_argument('--outdir', default='backtrack_tree_example')
    ap.add_argument('--max-nodes', type=int, default=20000)
    a = ap.parse_args()

    os.makedirs(a.outdir, exist_ok=True)
    panels = []

    for seed in a.seeds:
        tree, tours, ctx, forced, init = ptree.trace(a.n, a.K, seed, a.max_nodes)
        st = stats(tree, tours, a.n)
        ok = all(kt.verify_tour(t, a.n) for t in tours)

        data = {'n': a.n, 'K': a.K, 'seed': seed, 'nodes': tree.nodes,
                'children': {str(k): v for k, v in tree.children.items()}}
        base = f"{a.outdir}/arvore_{a.n}x{a.n}_seed{seed}"
        with open(base + '.json', 'w') as f:
            json.dump(data, f)
        svg = rtree.render(data, 'LR')
        with open(base + '.svg', 'w') as f:
            f.write(svg)

        panels.append((seed, st, svg, len(tours), ok, forced))
        print(f"seed {seed:3d}: {st['nos']:4d} nós  "
              f"corredor R2 até prof {st['corredor']}  "
              f"subciclos {st['subtour']:3d}  tours {len(tours)} (ok={ok})")

    html = [
        '<title>Árvores de backtracking — cavalo 6×6</title>',
        '<style>',
        ':root{--bg:#fbfbfd;--fg:#1a1c22;--mut:#5b6270;--line:#dfe3ea;--card:#fff}',
        ':root:not([data-theme="light"]){}',
        '@media(prefers-color-scheme:dark){:root:not([data-theme="light"])',
        '{--bg:#14161a;--fg:#e8eaef;--mut:#9aa2b1;--line:#2a2f38;--card:#1b1e24}}',
        ':root[data-theme="dark"]{--bg:#14161a;--fg:#e8eaef;--mut:#9aa2b1;',
        '--line:#2a2f38;--card:#1b1e24}',
        'body{background:var(--bg);color:var(--fg);',
        'font:14px/1.5 ui-sans-serif,system-ui,sans-serif;margin:0;padding:28px}',
        'h1{font-size:20px;margin:0 0 4px}',
        'p.sub{color:var(--mut);margin:0 0 26px}',
        '.card{background:var(--card);border:1px solid var(--line);',
        'border-radius:10px;margin-bottom:22px;overflow:hidden}',
        '.hd{padding:12px 16px;border-bottom:1px solid var(--line);',
        'display:flex;gap:18px;align-items:baseline;flex-wrap:wrap}',
        '.hd b{font-size:15px}',
        '.m{color:var(--mut);font:12px ui-monospace,Menlo,monospace}',
        '.m em{color:var(--fg);font-style:normal;font-weight:600}',
        '.scroll{overflow-x:auto;padding:10px}',
        '.scroll svg{display:block}',
        '</style>',
        f'<h1>Árvores de backtracking — cavalo {a.n}×{a.n}</h1>',
        f'<p class="sub">Mesma busca (K={a.K} tours), quatro seeds. A seed muda '
        'apenas o desempate na escolha de aresta — e com isso a raiz da região '
        'de ramificação.</p>',
    ]
    for seed, st, svg, ntours, ok, forced in panels:
        html += [
            '<div class="card"><div class="hd">',
            f'<b>seed {seed}</b>',
            f'<span class="m"><em>{st["nos"]}</em> nós</span>',
            f'<span class="m">corredor R2 até prof <em>{st["corredor"]}</em></span>',
            f'<span class="m"><em>{st["subtour"]}</em> subciclos</span>',
            f'<span class="m"><em>{st["contra"]}</em> contradições</span>',
            f'<span class="m"><em>{ntours}</em> tours · <em>{st["nos_por_tour"]}</em> nós/tour</span>',
            '</div><div class="scroll">', svg, '</div></div>',
        ]
    out = f"{a.outdir}/painel_{a.n}x{a.n}.html"
    with open(out, 'w') as f:
        f.write('\n'.join(html))
    print(f"\npainel: {out}")


if __name__ == '__main__':
    main()
