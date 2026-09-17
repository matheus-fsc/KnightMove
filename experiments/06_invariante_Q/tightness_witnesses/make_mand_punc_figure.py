"""Figura de dois paineis: Mand(n) e Punc(n) = G_n menos as arestas obrigatorias.

(a) os 4 cantos de grau 2 e suas 8 arestas obrigatorias
(b) o grafo perfurado: bulk conexo + 4 cantos isolados = 5 componentes

Uso: python make_mand_punc_figure.py --n 6 --out ../paper/fig_mand_punc.tex
"""
from __future__ import annotations
import argparse

MOVES = [(1,2),(2,1),(-1,2),(-2,1),(1,-2),(2,-1),(-1,-2),(-2,-1)]

def build(n):
    V = [(r,c) for r in range(n) for c in range(n)]
    S = set(V); E = []
    for (r,c) in V:
        for dr,dc in MOVES:
            b = (r+dr, c+dc)
            if b in S and (r,c) < b: E.append(((r,c), b))
    return V, sorted(E)

def components(V, E, skip):
    adj = {v:set() for v in V}
    for a,b in E:
        if frozenset((a,b)) in skip: continue
        adj[a].add(b); adj[b].add(a)
    seen=set(); comps=[]
    for s in V:
        if s in seen: continue
        comp={s}; seen.add(s); st=[s]
        while st:
            u=st.pop()
            for w in adj[u]:
                if w not in seen: seen.add(w); comp.add(w); st.append(w)
        comps.append(comp)
    return comps

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=6)
    ap.add_argument("--out", default="../paper/fig_mand_punc.tex")
    a = ap.parse_args(); n = a.n
    V, E = build(n)
    corners = [(0,0),(0,n-1),(n-1,0),(n-1,n-1)]
    mand = [e for e in E if e[0] in corners or e[1] in corners]
    mandset = {frozenset(e) for e in mand}
    comps = components(V, E, mandset)
    assert len(mand) == 8, f"esperava 8 obrigatorias, achei {len(mand)}"
    assert len(comps) == 5, f"esperava 5 componentes, achei {len(comps)}"
    beta1 = len(E) - len(V) + 1
    print(f"n={n}: |V|={len(V)} |E|={len(E)} beta1={beta1} |Mand|={len(mand)} componentes(Punc)={len(comps)}")

    def xy(u): return (u[1], n-1-u[0])

    def panel(edges, hl, capt, xshift, isolate):
        L = [r"\begin{scope}[xshift=%.2fcm]" % xshift,
             r"\draw[step=1,black!10,very thin] (0,0) grid (%d,%d);" % (n,n),
             r"\draw[black!35] (0,0) rectangle (%d,%d);" % (n,n)]
        for e in edges:
            (x1,y1),(x2,y2) = xy(e[0]), xy(e[1])
            st = "mand" if frozenset(e) in hl else "bulk"
            L.append(r"\draw[%s] (%.1f,%.1f) -- (%.1f,%.1f);" % (st,x1+.5,y1+.5,x2+.5,y2+.5))
        for u in V:
            x,y = xy(u)
            if u in corners:
                st = "cisol" if isolate else "cnr"
                L.append(r"\fill[%s] (%.1f,%.1f) circle (3.4pt);" % (st,x+.5,y+.5))
            else:
                L.append(r"\fill[blk] (%.1f,%.1f) circle (1.5pt);" % (x+.5,y+.5))
        L += [r"\node[anchor=north,font=\small] at (%.1f,-0.30) {%s};" % (n/2,capt),
              r"\end{scope}"]
        return "\n".join(L)

    keep = [e for e in E if frozenset(e) not in mandset]
    T = [r"% gerado por tightness_witnesses/make_mand_punc_figure.py -- NAO EDITAR A MAO",
         r"\begin{figure}[!ht]", r"\centering",
         r"\tikzset{bulk/.style={line width=0.4pt,black!18},",
         r"  mand/.style={line width=1.7pt,red!75!black},",
         r"  cnr/.style={fill=red!75!black}, cisol/.style={fill=red!75!black},",
         r"  blk/.style={fill=black!55}}",
         r"\begin{tikzpicture}[scale=0.46]"]
    T.append(panel(E, mandset, r"(a) 8 obrigat\'orias", 0, False))
    T.append(panel(keep, set(), r"(b) 5 componentes", n+5.0, True))
    T += [r"\end{tikzpicture}",
          r"\caption{\textbf{Os cantos e o grafo perfurado} ($n = %d$)." % n,
          r"  \textbf{(a)} Cada canto tem grau exatamente~2, logo suas duas arestas",
          r"  pertencem a \emph{todo} tour: s\~ao as $|\Mand(n)| = 8$ arestas",
          r"  obrigat\'orias (vermelho). Um elemento de $Z_1$ usa as duas ou nenhuma,",
          r"  o que define os quatro funcionais de canto.",
          r"  \textbf{(b)} Removendo-as, o grafo se parte em exatamente",
          r"  \textbf{cinco} componentes: o \emph{bulk} --- conexo, o \'unico",
          r"  ingrediente geom\'etrico do Teorema~\ref{thm:Qn3} --- e os quatro",
          r"  cantos, agora isolados. Da\'i",
          r"  $\dim Z_{\mathrm{bulk}} = \beta_1 - 4$ (Lema~\ref{lem:dim-zbulk}),",
          r"  formalizado como \texttt{card\_connectedComponent\_punc}",
          r"  (Se\c{c}\~ao~\ref{sec:lean}).}",
          r"\label{fig:mand-punc}", r"\end{figure}"]
    open(a.out,"w").write("\n".join(T))
    print("figura ->", a.out)

if __name__ == "__main__":
    main()
