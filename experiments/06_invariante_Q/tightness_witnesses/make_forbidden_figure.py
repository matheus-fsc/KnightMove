"""Figura dos 3 ciclos proibidos do 6x6, no estilo do paper.

Le forbidden_cycles_6x6/data/results/forbidden_cycles.json (rotulos de aresta)
e emite TikZ com 3 tabuleiros. Nao recalcula nada.
"""
import json, argparse, re

def sq(lbl):                      # 'A6' -> (col, row) 0-indexado, A=0, linha 1 embaixo
    c = ord(lbl[0]) - ord('A'); r = int(lbl[1:]) - 1
    return c, r

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--src', default='../forbidden_cycles_6x6/data/results/forbidden_cycles.json')
    ap.add_argument('--out', default='../paper/fig_forbidden.tex')
    a = ap.parse_args(); n = 6
    d = json.load(open(a.src))
    cycles = [(f['index'], f['support_size'], f['edge_labels']) for f in d['forbidden']]

    def board(labels, capt, xshift):
        L = [r"\begin{scope}[xshift=%.2fcm]" % xshift,
             r"\draw[step=1,black!12,very thin] (0,0) grid (%d,%d);" % (n, n),
             r"\draw[black!35] (0,0) rectangle (%d,%d);" % (n, n)]
        verts = set()
        for lab in labels:
            u, v = lab.split('-'); (x1, y1), (x2, y2) = sq(u), sq(v)
            verts |= {(x1, y1), (x2, y2)}
            L.append(r"\draw[fcyc] (%.1f,%.1f) -- (%.1f,%.1f);" % (x1+.5, y1+.5, x2+.5, y2+.5))
        for (x, y) in sorted(verts):
            L.append(r"\fill[fvert] (%.1f,%.1f) circle (2.6pt);" % (x+.5, y+.5))
        L += [r"\node[anchor=north,font=\small] at (%.1f,-0.30) {%s};" % (n/2, capt),
              r"\end{scope}"]
        return "\n".join(L)

    T = [r"% gerado por tightness_witnesses/make_forbidden_figure.py -- NAO EDITAR A MAO",
         r"\begin{figure}[!ht]", r"\centering",
         r"\tikzset{fcyc/.style={line width=1.5pt,red!72!black},",
         r"         fvert/.style={fill=black!80}}",
         r"\begin{tikzpicture}[scale=0.46]"]
    for k, (idx, size, labs) in enumerate(cycles):
        T.append(board(labs, r"(%s) $f_%d$, $|S|=%d$" % ("abc"[k], idx, size), k*(n+2.2)))
    T += [r"\end{tikzpicture}",
          r"\caption{\textbf{Os tr\^es ciclos proibidos do $6\times6$.} Cada $f_i$ est\'a em",
          r"  $Z_1(G_6;\Fdois)$ mas \emph{n\~ao} em $\Span(\Ham(6))$: juntos representam as",
          r"  tr\^es dimens\~oes de deficit ($\beta_1 = 45$, $\rk(\Ham) = 42$), com suportes",
          r"  de tamanho $6$, $6$ e $8$. S\~ao \emph{representantes de cosets}, obtidos como",
          r"  ciclos fundamentais de uma \'arvore geradora BFS; a geometria de cada um",
          r"  depende dessa escolha e n\~ao \'e can\^onica --- o que \'e intr\'inseco \'e o",
          r"  subespa\c{c}o que eles geram no quociente. Fonte: \texttt{forbidden\_cycles\_6x6/}.}",
          r"\label{fig:proibidos}", r"\end{figure}"]
    open(a.out, 'w').write("\n".join(T))
    print("figura ->", a.out, "|", len(cycles), "ciclos")

if __name__ == "__main__":
    main()
