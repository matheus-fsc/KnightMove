"""Gera a figura do certificado: H_A, H_B e H_A xor H_B = hexagono.

Produz TikZ com 3 tabuleiros. Usa um certificado REAL, auditado.
Uso:  python make_certificate_figure.py --n 8 --out ../paper/fig_certificate.tex
"""
from __future__ import annotations
import argparse, itertools
from typing import List, Set, Tuple
from switcher_search import build, build_switchers, ham_path

def hexes_in_bulk(adj, n, corner_set, limit=4000):
    """6-ciclos simples que evitam cantos."""
    out=[]
    for a in range(n*n):
        if a in corner_set: continue
        for path in dfs_cycles(adj, a, corner_set, 6):
            h=tuple(path)
            if min(h)==a and h[1]<h[-1]: out.append(h)
            if len(out)>=limit: return out
    return out

def dfs_cycles(adj, start, corner_set, L):
    stack=[(start,[start],{start})]
    while stack:
        v,path,seen=stack.pop()
        if len(path)==L:
            if start in adj[v]: yield path
            continue
        for w in adj[v]:
            if w in corner_set or w in seen or w<start: continue
            stack.append((w,path+[w],seen|{w}))

def edges(seq, closed=False):
    p=list(zip(seq,seq[1:]))+([(seq[-1],seq[0])] if closed else [])
    return {frozenset(e) for e in p}

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--n',type=int,default=8)
    ap.add_argument('--out',default='../paper/fig_certificate.tex')
    a=ap.parse_args(); n=a.n
    adj,_=build(n) if isinstance(build(n),tuple) else (build(n),None)
    corner_set={0,n-1,n*(n-1),n*n-1}
    print(f"n={n}, procurando hexagono no bulk...")
    cands=hexes_in_bulk(adj,n,corner_set)
    ctr=(n-1)/2
    def centrality(h):
        return max(max(abs(u%n-ctr),abs(u//n-ctr)) for u in h)
    cands.sort(key=centrality)
    for h in cands:
        for rot in range(3):
            for w in build_switchers(adj,h,rot,max_internal=3,max_w=200):
                v=w["cycle"]; i2,i3=w["P2_internal"],w["P3_internal"]
                VW=set(v)|set(i2)|set(i3)
                allowed=set(range(n*n))-(VW-{v[0],v[3]})
                hp,exh,_=ham_path(adj,allowed,v[0],v[3])
                if hp is None: continue
                pathA=[v[0],v[1]]+i2+[v[5],v[4]]+i3[::-1]+[v[2],v[3]]
                pathB=[v[0],v[5]]+i2[::-1]+[v[1],v[2]]+i3+[v[4],v[3]]
                HA=edges(pathA)|edges(hp); HB=edges(pathB)|edges(hp)
                C=edges(list(v),closed=True)
                assert HA^HB==C, "XOR != hexagono"
                assert all(len({x for e in H for x in e})==n*n for H in (HA,HB))
                print(f"  certificado OK | hex={v} | P2={i2} P3={i3} | |V(W)|={len(VW)}")
                emit(a.out,n,HA,HB,C,v,VW,i2,i3)
                return
    print("nenhum certificado encontrado")

def emit(path,n,HA,HB,C,v,VW,i2,i3):
    def xy(u): return (u%n, n-1-u//n)
    shared=HA&HB; onlyA=HA-shared; onlyB=HB-shared
    def board(H,capt,mark_a=frozenset(),mark_b=frozenset(),dots=True):
        L=[r"\draw[step=1,black!12,very thin] (0,0) grid (%d,%d);"%(n,n),
           r"\draw[black!35] (0,0) rectangle (%d,%d);"%(n,n)]
        for e in sorted(H,key=lambda t:sorted(t)):
            p_,q_=sorted(e); (x1,y1),(x2,y2)=xy(p_),xy(q_)
            st = "cedA" if e in mark_a else ("cedB" if e in mark_b else "shared")
            L.append(r"\draw[%s] (%.1f,%.1f) -- (%.1f,%.1f);"%(st,x1+.5,y1+.5,x2+.5,y2+.5))
        if dots:
            for u in set(i2)|set(i3):
                x,y=xy(u); L.append(r"\fill[black!45] (%.1f,%.1f) circle (3pt);"%(x+.5,y+.5))
            for u in v:
                x,y=xy(u); L.append(r"\fill[hexv] (%.1f,%.1f) circle (4pt);"%(x+.5,y+.5))
        L.append(r"\node[anchor=north,font=\small] at (%.1f,-0.35) {%s};"%(n/2,capt))
        return "\n".join(L)
    T=[r"% gerado por tightness_witnesses/make_certificate_figure.py -- NAO EDITAR A MAO",
       r"\begin{figure}[!ht]",r"\centering",
       r"\tikzset{shared/.style={line width=0.45pt,black!20},",
       r"  cedA/.style={line width=1.6pt,red!75!black},",
       r"  cedB/.style={line width=1.6pt,blue!70!black},",
       r"  hexv/.style={fill=black!80}}",
       r"\begin{tikzpicture}[scale=0.34]"]
    T.append(board(HA,r"(a) $H_A$",mark_a=onlyA))
    T.append(r"\begin{scope}[xshift=%dcm]"%(n+3)); T.append(board(HB,r"(b) $H_B$",mark_b=onlyB)); T.append(r"\end{scope}")
    T.append(r"\begin{scope}[xshift=%dcm]"%(2*(n+3))); T.append(board(C,r"(c) $H_A \oplus H_B = C$",mark_a=onlyA,mark_b=onlyB)); T.append(r"\end{scope}")
    T+=[r"\end{tikzpicture}",
        r"\caption{\textbf{Um certificado.} Dois tours fechados do cavalo no $%d\times%d$"%(n,n),
        r"  \emph{que coincidem fora de um gadget de oito v\'ertices}: em cinza as arestas",
        r"  comuns; em vermelho as %d que s\'o est\~ao em $H_A$; em azul as %d que s\'o est\~ao"%(len(onlyA),len(onlyB)),
        r"  em $H_B$. A diferen\c{c}a sim\'etrica (c) \'e o hex\'agono $C$ --- todo o resto",
        r"  cancela ---, exibindo $C \in \Span(\Ham(n))$ por constru\c{c}\~ao.}",
        r"\label{fig:certificado}",r"\end{figure}"]
    open(path,'w').write("\n".join(T))
    print("  figura ->",path)

if __name__=="__main__": main()
