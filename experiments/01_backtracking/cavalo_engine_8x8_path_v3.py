#!/usr/bin/env python3
"""
cavalo_engine_8x8_path.py  v3
INCLUI:
  1. Aresta virtual (S,E) para GF(2)
  2. Quebra de Simetria Diedral (Symmetry Breaking)
  3. Poda A Priori de Subciclos de 4 nós (Anti-Short-Circuiting)
"""
import json, os, math
from itertools import combinations
from z3 import Solver, Bool, Sum, If, sat, is_true, Or

BOARD = 8
AMOSTRAS_POR_PAR = int(os.getenv("AMOSTRAS_POR_PAR", "500"))
PARES_CANONICOS = [
    ((0,0),(0,7)),  # canto->canto mesma borda   0->1
    ((0,0),(7,6)),  # canto->borda diagonal       0->1
    ((0,0),(6,5)),  # canto->interior profundo    0->1
    ((0,0),(4,3)),  # canto->quase centro         0->1
    ((0,1),(6,4)),  # borda->interior medio       1->0
    ((1,0),(6,3)),  # borda->interior oposto      1->0
]
KNIGHT_DELTAS = [(1,2),(2,1),(-1,2),(-2,1),(1,-2),(2,-1),(-1,-2),(-2,-1)]

def knight_nbrs(x, y):
    out = []
    for dx,dy in KNIGHT_DELTAS:
        nx,ny = x+dx, y+dy
        if 0<=nx<BOARD and 0<=ny<BOARD: out.append((nx,ny))
    def deg(n): return sum(1 for dx,dy in KNIGHT_DELTAS if 0<=n[0]+dx<BOARD and 0<=n[1]+dy<BOARD)
    return sorted(out, key=deg)

def build_graph():
    nodes = [(x,y) for x in range(BOARD) for y in range(BOARD)]
    adj   = {n: knight_nbrs(*n) for n in nodes}
    edges = set()
    for n,nbrs in adj.items():
        for m in nbrs: edges.add(tuple(sorted([n,m])))
    return nodes, adj, edges

def get_ciclos_de_4(nodes, adj):
    """Mapeia todos os ciclos de 4 nós usando apenas arestas reais."""
    ciclos = set()
    for u in nodes:
        for v in adj[u]:
            for w in adj[v]:
                if w != u:
                    for x in adj[w]:
                        if x != v and u in adj[x]:
                            ciclo = tuple(sorted([u, v, w, x]))
                            ciclos.add(ciclo)
    return ciclos

def spanning_tree(adj, root):
    parent={root:None}; stack=[root]; visited={root}; te=set()
    while stack:
        u=stack.pop()
        for v in reversed(adj[u]):
            if v not in visited:
                visited.add(v); parent[v]=u; te.add(tuple(sorted([u,v]))); stack.append(v)
    return parent, te

def extract_loops(parent, te, all_edges):
    loops=[]
    for e in sorted(all_edges):
        if e in te: continue
        u,v=e
        pu=[]; c=u
        while c is not None: pu.append(c); c=parent[c]
        pv=[]; c=v
        while c is not None: pv.append(c); c=parent[c]
        loops.append({"c1":[list(p) for p in pu],"c2":[list(p) for p in pv],"col":[list(u),list(v)]})
    return loops

def get_loop_edges(lp):
    es=set()
    def tog(u,v):
        e=tuple(sorted([tuple(u),tuple(v)]))
        if e in es: es.remove(e)
        else: es.add(e)
    for i in range(len(lp["c1"])-1): tog(lp["c1"][i],lp["c1"][i+1])
    for i in range(len(lp["c2"])-1): tog(lp["c2"][i],lp["c2"][i+1])
    tog(lp["col"][0],lp["col"][1])
    return es

def solve_path(nodes, adj, edges_orig, start, end, n_amostras):
    virtual = tuple(sorted([start, end]))
    edges_aug = set(edges_orig) | {virtual}
    adj_aug = {n:[] for n in nodes}
    for u,v in edges_orig: adj_aug[u].append(v); adj_aug[v].append(u)
    adj_aug[start].append(end); adj_aug[end].append(start)
    par, te = spanning_tree(adj_aug, nodes[0])
    loops   = extract_loops(par, te, edges_aug)
    print(f"  H1_aug={len(loops)}")

    solver = Solver()
    lvars  = [Bool(f"L{i}") for i in range(len(loops))]
    e2l = {e:[] for e in edges_aug}
    for i,lp in enumerate(loops):
        for e in get_loop_edges(lp): e2l[e].append(i)

    def ea(e):
        idxs = e2l.get(e,[])
        if not idxs: return False
        return (Sum([If(lvars[i],1,0) for i in idxs]) % 2 == 1)

    # Grau 2 para todos (ciclo no grafo aumentado)
    for n in nodes:
        inc = [e for e in edges_aug if n in e]
        solver.add(Sum([If(ea(e),1,0) for e in inc]) == 2)

    # Aresta virtual SEMPRE ativa
    solver.add(ea(virtual) == True)

    # -----------------------------------------------------------------
    # INJEÇÃO 1: QUEBRA DE SIMETRIA (DIHEDRAL SYMMETRY BREAKING)
    # -----------------------------------------------------------------
    if start == (0, 0):
        edge_sym = tuple(sorted([(0, 0), (1, 2)]))
        if edge_sym in edges_aug:
            solver.add(ea(edge_sym) == True)

    # -----------------------------------------------------------------
    # INJEÇÃO 2: PODA INVERSA DE SUBCICLOS DE 4 NÓS
    # -----------------------------------------------------------------
    c4 = get_ciclos_de_4(nodes, adj)
    for ciclo in c4:
        e1 = tuple(sorted([ciclo[0], ciclo[1]]))
        e2 = tuple(sorted([ciclo[1], ciclo[2]]))
        e3 = tuple(sorted([ciclo[2], ciclo[3]]))
        e4 = tuple(sorted([ciclo[3], ciclo[0]]))

        if all(e in edges_aug for e in [e1, e2, e3, e4]):
            solver.add(Sum(If(ea(e1),1,0), If(ea(e2),1,0), If(ea(e3),1,0), If(ea(e4),1,0)) <= 3)

    edges_s = sorted(edges_orig)
    sigs = []

    # -----------------------------------------------------------------
    # INJEÇÃO 3: TESTE DE RUPTURA TOPOLÓGICA (PROVA DE INVARIANTE)
    # -----------------------------------------------------------------
    # Escolha duas dimensões que apresentaram r = -1.0 no seu JSON
    # Por exemplo, para o par (0,0)->(0,7), as dimensões 47 e 71.
    
    idx_a = 47
    idx_b = 71
    
    edge_a = edges_s[idx_a]
    edge_b = edges_s[idx_b]
    
    # Forçamos o solver a passar pelas duas arestas simultaneamente
    solver.add(ea(edge_a) == True)
    solver.add(ea(edge_b) == True)
    
    print(f"\n  [TESTE DE RUPTURA] Tentando forçar passagem por dim{idx_a} e dim{idx_b} simultaneamente...")
    
    # Se a correlação -1.0 for real, o Z3 não conseguirá encontrar NENHUMA 
    # solução e retornará 'unsat' imediatamente.
    if solver.check() != sat:
        print("  [RESULTADO] UNSAT! Ruptura confirmada. É matematicamente impossível!")
        # return [], edges_s # Aborta a função, pois provamos o ponto
    else:
        print("  [RESULTADO] SAT? A correlação não era perfeita!")
    # -----------------------------------------------------------------

    print(f"  Amostrando {n_amostras} caminhos...")
    while len(sigs) < n_amostras:
        if solver.check() != sat:
            print(f"  Espaco esgotado: {len(sigs)} amostras."); break
        m = solver.model()
        active = [e for e in edges_orig if is_true(m.evaluate(ea(e)))]
        adj_r={n:[] for n in nodes}
        for u,v in active: adj_r[u].append(v); adj_r[v].append(u)

        vis={start}; q=[start]
        while q:
            c=q.pop(0)
            for nb in adj_r[c]:
                if nb not in vis: vis.add(nb); q.append(nb)

        if len(vis)==len(nodes):
            sigs.append([is_true(m.evaluate(ea(e))) for e in edges_s])
            if len(sigs)%100==0: print(f"  -> {len(sigs)}")

        # Força o Z3 a buscar uma nova solução
        solver.add(Or([lvars[i] != is_true(m.evaluate(lvars[i])) for i in range(len(lvars))]))

        # Lazy Constraint: Se gerou subciclos maiores, adiciona restrição dinâmica
        if len(vis)<len(nodes):
            seen=set(vis)
            for seed in [n for n in nodes if n not in seen]:
                if seed not in seen:
                    comp={seed}; q2=[seed]
                    while q2:
                        c=q2.pop(0)
                        for nb in adj_r[c]:
                            if nb not in comp: comp.add(nb); q2.append(nb)
                    in_e=[e for e in edges_aug if e[0] in comp and e[1] in comp]
                    solver.add(Sum([If(ea(e),1,0) for e in in_e]) <= len(comp)-1)
                    seen|=comp

    return sigs, edges_s

def analisar(sigs, edges_s, label):
    n=len(sigs); nd=len(sigs[0]) if sigs else 0
    if n<10: print(f"  [{label}] insuficiente ({n})."); return {"label":label,"n":n}
    freq=[sum(s[i] for s in sigs)/n for i in range(nd)]
    pesos=[sum(s) for s in sigs]
    med=sum(pesos)/n; std=math.sqrt(sum((p-med)**2 for p in pesos)/n)
    obr=[i for i,f in enumerate(freq) if f==1.0]
    imp=[i for i,f in enumerate(freq) if f==0.0]
    liv=[i for i,f in enumerate(freq) if 0<f<1]

    corrs=[]
    for i,j in combinations(liv,2):
        fi,fj=freq[i],freq[j]
        cov=sum((sigs[k][i]-fi)*(sigs[k][j]-fj) for k in range(n))/n
        si=math.sqrt(fi*(1-fi)); sj=math.sqrt(fj*(1-fj))
        if si>0 and sj>0: corrs.append((cov/(si*sj),i,j))

    corrs.sort()
    alvo,thr=-0.771,0.05
    persiste=[(r,i,j) for r,i,j in corrs if abs(r-alvo)<thr]

    print(f"\n  [{label}] n={n}")
    print(f"  Hamming: media={med:.1f}  std={std:.2f}  esperado~{nd//2}")
    print(f"  Obrig={len(obr)} Impos={len(imp)} Livres={len(liv)}")
    print(f"  Top correlacoes negativas:")
    for r,i,j in corrs[:8]: print(f"    dim{i}<->dim{j}  r={r:.4f}")

    if persiste:
        print(f"  *** INVARIANTE r~{alvo} detectado: {len(persiste)} pares ***")
        for r,i,j in persiste[:5]: print(f"      dim{i}<->dim{j}  r={r:.4f}")
    else:
        print(f"  Invariante r~{alvo} nao detectado.")

    return {"label":label,"n":n,"nd":nd,"med":med,"std":std,
            "obr":len(obr),"imp":len(imp),"liv":len(liv),
            "top_neg":corrs[:10],"persiste_0771":persiste}

    def testar_ruptura(nodes, adj, edges_orig, start, end, e1, e2):
        """
        Modo Juiz: Instancia um solver limpo, aplica as regras da física do tabuleiro
        e tenta forçar a passagem pelas arestas e1 e e2 simultaneamente.
        Retorna True se for UNSAT (Lei Provada), False se for SAT (Falso Positivo).
        """
        from z3 import Solver, Bool, Sum, If, sat
        virtual = tuple(sorted([start, end]))
        edges_aug = set(edges_orig) | {virtual}
        
        adj_aug = {n:[] for n in nodes}
        for u,v in edges_orig: adj_aug[u].append(v); adj_aug[v].append(u)
        adj_aug[start].append(end); adj_aug[end].append(start)
        
        par, te = spanning_tree(adj_aug, nodes[0])
        loops = extract_loops(par, te, edges_aug)
        
        solver = Solver()
        lvars = [Bool(f"J_{i}") for i in range(len(loops))]
        e2l = {e:[] for e in edges_aug}
        for i,lp in enumerate(loops):
            for e in get_loop_edges(lp): e2l[e].append(i)
            
        def ea(e):
            idxs = e2l.get(e,[])
            if not idxs: return False
            return (Sum([If(lvars[i],1,0) for i in idxs]) % 2 == 1)

        # 1. Regra do Grau 2
        for n in nodes:
            inc = [e for e in edges_aug if n in e]
            solver.add(Sum([If(ea(e),1,0) for e in inc]) == 2)

        # 2. Aresta Virtual Ativa
        solver.add(ea(virtual) == True)
        
        # 3. Quebra de Simetria Diedral
        if start == (0, 0):
            edge_sym = tuple(sorted([(0, 0), (1, 2)]))
            if edge_sym in edges_aug:
                solver.add(ea(edge_sym) == True) 
                
        # 4. Poda a Priori de Ciclos de 4
        c4 = get_ciclos_de_4(nodes, adj)
        for ciclo in c4:
            e_c1 = tuple(sorted([ciclo[0], ciclo[1]]))
            e_c2 = tuple(sorted([ciclo[1], ciclo[2]]))
            e_c3 = tuple(sorted([ciclo[2], ciclo[3]]))
            e_c4 = tuple(sorted([ciclo[3], ciclo[0]]))
            if all(e in edges_aug for e in [e_c1, e_c2, e_c3, e_c4]):
                solver.add(Sum(If(ea(e_c1),1,0), If(ea(e_c2),1,0), If(ea(e_c3),1,0), If(ea(e_c4),1,0)) <= 3)

        # --- O GOLPE FINAL: FORÇAR A RUPTURA ---
        solver.add(ea(e1) == True)
        solver.add(ea(e2) == True)
        
        # Se retornar diferente de sat (UNSAT), a lei é verdadeira!
        return solver.check() != sat

def main():
    print(f"Grafo {BOARD}x{BOARD}")
    nodes,adj,edges=build_graph()
    print(f"  V={len(nodes)} E={len(edges)} H1={len(edges)-len(nodes)+1}")
    pares=[]

    print("\nParidade:")
    for s,e in PARES_CANONICOS:
        cs,ce=(s[0]+s[1])%2,(e[0]+e[1])%2
        ok="OK" if cs!=ce else "SKIP"
        print(f"  {ok} {s}->{e}")
        if cs!=ce: pares.append((s,e))

    todos=[]
    for start,end in pares:
        print(f"\n{'='*55}\nPar: {start} -> {end}\n{'='*55}")
        label=f"{start[0]}{start[1]}_{end[0]}{end[1]}"
        # Adicionado o parâmetro 'adj' na chamada da função
        sigs,es=solve_path(nodes,adj,edges,start,end,AMOSTRAS_POR_PAR)
        res=analisar(sigs,es,label)
        res["sigs"]=sigs
        todos.append(res)

    out={"board":BOARD,"amostras":AMOSTRAS_POR_PAR,
         "resultados":[{k:v for k,v in r.items() if k!="sigs"} for r in todos],
         "assinaturas":{r["label"]:r["sigs"] for r in todos if "sigs" in r}}

    with open("resultados_8x8_v3.json","w") as f: json.dump(out,f,indent=2)
    print("\nSalvo: resultados_8x8_v3.json")

if __name__=="__main__":
    main()
