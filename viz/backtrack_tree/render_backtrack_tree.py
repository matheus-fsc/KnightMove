"""render_backtrack_tree.py — desenha a árvore de backtracking em SVG puro.

Não depende de graphviz nem de nenhuma biblioteca externa: lê o JSON
produzido por print_backtrack_tree.py --json e emite um SVG.

Uso:
    python render_backtrack_tree.py arvore.json -o arvore.svg
    python render_backtrack_tree.py arvore.json -o arvore.svg --orient TB
"""

import argparse
import json

COLOR = {
    'TOUR':    ('#1a7f37', '#d7f5dd'),
    'CONTRA':  ('#b42318', '#fde3e1'),
    'SUBTOUR': ('#b25a00', '#fdeccc'),
    'DEADEND': ('#666666', '#e8e8e8'),
    'OK':      ('#24406e', '#e4ecf9'),
}


def sq(v, n):
    r, c = divmod(int(v), n)
    return f"{chr(ord('a') + c)}{r + 1}"


def layout(nodes, children, orient):
    """Tidy-tree simples: folhas em ordem, internos na média dos filhos."""
    pos = {}
    counter = [0]

    def walk(nid):
        kids = children.get(str(nid), children.get(nid, []))
        if not kids:
            pos[nid] = counter[0]
            counter[0] += 1
        else:
            for c in kids:
                walk(c)
            pos[nid] = (pos[kids[0]] + pos[kids[-1]]) / 2.0

    for r in children.get('-1', children.get(-1, [])):
        walk(r)
    return pos


def render(data, orient='LR'):
    n = data['n']
    nodes = {nd['id']: nd for nd in data['nodes']}
    children = data['children']
    cross = layout(nodes, children, orient)

    NW, NH = 150, 26          # caixa
    GAP_D, GAP_C = 172, 34    # espaçamento por profundidade / entre irmãos
    PAD = 30

    def xy(nid):
        nd = nodes[nid]
        if orient == 'LR':
            return PAD + nd['depth'] * GAP_D, PAD + cross[nid] * GAP_C
        return PAD + cross[nid] * (NW + 12), PAD + nd['depth'] * GAP_C * 1.6

    maxd = max(nd['depth'] for nd in nodes.values())
    nleaf = max(cross.values()) + 1
    if orient == 'LR':
        W = PAD * 2 + maxd * GAP_D + NW
        H = PAD * 2 + nleaf * GAP_C + NH
    else:
        W = PAD * 2 + nleaf * (NW + 12)
        H = PAD * 2 + maxd * GAP_C * 1.6 + NH

    o = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W:.0f}" '
         f'height="{H:.0f}" viewBox="0 0 {W:.0f} {H:.0f}" '
         f'font-family="ui-monospace,Menlo,Consolas,monospace">',
         f'<rect width="{W:.0f}" height="{H:.0f}" fill="#fbfbfd"/>']

    # arestas primeiro (ficam atrás)
    for nid, nd in nodes.items():
        p = nd['parent']
        if p < 0:
            continue
        x1, y1 = xy(p)
        x2, y2 = xy(nid)
        if orient == 'LR':
            xa, ya = x1 + NW, y1 + NH / 2
            xb, yb = x2, y2 + NH / 2
            mid = (xa + xb) / 2
            d = f"M{xa:.1f},{ya:.1f} C{mid:.1f},{ya:.1f} {mid:.1f},{yb:.1f} {xb:.1f},{yb:.1f}"
        else:
            xa, ya = x1 + NW / 2, y1 + NH
            xb, yb = x2 + NW / 2, y2
            mid = (ya + yb) / 2
            d = f"M{xa:.1f},{ya:.1f} C{xa:.1f},{mid:.1f} {xb:.1f},{mid:.1f} {xb:.1f},{yb:.1f}"
        w = 2.0 if nd['found_tour'] else 1.0
        col = '#1a7f37' if nd['found_tour'] else '#b9c0cc'
        o.append(f'<path d="{d}" fill="none" stroke="{col}" stroke-width="{w}"/>')

    # nós
    for nid, nd in nodes.items():
        x, y = xy(nid)
        stroke, fill = COLOR.get(nd['status'], COLOR['OK'])
        sw = 2.0 if nd['found_tour'] else 1.0
        lbl = (f"{sq(nd['u'], n)}-{sq(nd['v'], n)}="
               f"{nd['val']}  liv:{nd['n_free']}")
        o.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{NW}" height="{NH}" '
                 f'rx="5" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"/>')
        o.append(f'<text x="{x + 8:.1f}" y="{y + 17:.1f}" font-size="11.5" '
                 f'fill="#111">{lbl}</text>')
        if nd['status'] != 'OK':
            mk = {'TOUR': '✔', 'SUBTOUR': '↺', 'CONTRA': '✘', 'DEADEND': '·'}
            o.append(f'<text x="{x + NW - 14:.1f}" y="{y + 17:.1f}" '
                     f'font-size="13" fill="{stroke}">{mk[nd["status"]]}</text>')

    # legenda
    lx, ly = PAD, H - 14
    items = [('OK', 'propaga'), ('SUBTOUR', 'subciclo (↺)'),
             ('CONTRA', 'contradição (✘)'), ('TOUR', 'tour completo (✔)')]
    for i, (k, txt) in enumerate(items):
        s, f = COLOR[k]
        cx = lx + i * 170
        o.append(f'<rect x="{cx}" y="{ly - 10}" width="12" height="12" rx="3" '
                 f'fill="{f}" stroke="{s}"/>')
        o.append(f'<text x="{cx + 18}" y="{ly}" font-size="11" fill="#333">{txt}</text>')

    o.append('</svg>')
    return '\n'.join(o)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('json_file')
    ap.add_argument('-o', '--out', default='arvore.svg')
    ap.add_argument('--orient', choices=['LR', 'TB'], default='LR')
    a = ap.parse_args()
    with open(a.json_file) as f:
        data = json.load(f)
    svg = render(data, a.orient)
    with open(a.out, 'w') as f:
        f.write(svg)
    print(f"{a.out}  ({len(data['nodes'])} nós, {len(svg) // 1024} KB)")


if __name__ == '__main__':
    main()
