"""
Análise do espaço residual após propagação em cascata.

Como propagation_engine concluiu que nenhuma das regras R2..R6 fixa
nada novo a partir das 8 obrigatórias dos cantos, o espaço residual
tem n_free = 72.

Esta análise descreve a estrutura desse residual:
  - rótulo, posição e grau no subgrafo induzido pelas livres
  - grafo de dependência: nó = aresta livre, ligação =
    coocorrer em algum par/tripla/quadra/XOR ainda ativa
  - componentes conexas
  - distribuição condicional P(x_e=1 | x_{e'}=1) nos tours reais

Saída: data/residual_variables.json
"""

from __future__ import annotations

import json
from collections import defaultdict
from itertools import combinations
from pathlib import Path

import numpy as np

from propagation_engine import (
    DATA_6x6,
    DATA_XOR,
    FREE,
    PropagationEngine,
    ROOT,
    carregar_dados,
    rotulo_aresta,
)


def construir_residual(dados):
    eng = PropagationEngine(
        edges_uv=dados["edges_uv"],
        freq=dados["freq"],
        pares=dados["pares"],
        triplas=dados["triplas"],
        quadras=dados["quadras"],
        xor=dados["xor"],
        xor_paridade=dados["xor_paridade_alvo"],
    )
    eng.run_to_fixpoint([1, 2, 3, 4, 5, 6])
    return eng


def grafo_dependencia(eng, dados, livres_set):
    """Constrói grafo (nós = arestas livres) das dependências residuais.

    Restrição residual ativa = constraint cuja avaliação ainda depende
    de >1 variável livre (suas demais variáveis já fixadas ou ela é
    intrinsecamente multivariada). Aqui adicionamos arestas par a par
    para cada constraint em que ambas variáveis livres aparecem.
    """
    arestas = defaultdict(set)  # e -> conjunto de e' co-restritos

    def adicionar(constraint_vars, fonte):
        livres = [e for e in constraint_vars if eng.estado[e] == FREE]
        if len(livres) < 2:
            return
        for a, b in combinations(livres, 2):
            arestas[a].add(b)
            arestas[b].add(a)
        return fonte

    for p in dados["pares"]:
        adicionar(p, "par")
    for t in dados["triplas"]:
        adicionar(t, "tripla")
    for q in dados["quadras"]:
        adicionar(q, "quadra")
    for s in dados["xor"]:
        adicionar(s, "xor")
    # estrutura de incidência (grau-2 acopla arestas no mesmo vértice)
    for v, es in eng.v2e.items():
        adicionar(tuple(es), "vertice")

    # componentes conexas via BFS
    componentes = []
    visit = set()
    for n in livres_set:
        if n in visit:
            continue
        comp = []
        stack = [n]
        while stack:
            x = stack.pop()
            if x in visit:
                continue
            visit.add(x)
            comp.append(x)
            for y in arestas[x]:
                if y not in visit and y in livres_set:
                    stack.append(y)
        componentes.append(sorted(comp))
    componentes.sort(key=len, reverse=True)
    return arestas, componentes


def distribuicao_condicional(inc, livres, freq, amostragem_max=300):
    """Calcula coeficiente de correlação (não normalizado) e P(e|e').

    Retorna estatísticas agregadas — não a matriz inteira para evitar
    JSON gigante.
    """
    livres = list(livres)
    inc_l = inc[:, livres].astype(np.float64)
    f = inc_l.mean(axis=0)
    cov = (inc_l.T @ inc_l) / inc_l.shape[0] - np.outer(f, f)
    std = np.sqrt(np.clip(f * (1 - f), 1e-12, None))
    corr = cov / np.outer(std, std)
    np.fill_diagonal(corr, 0.0)

    abs_corr = np.abs(corr)
    triu = np.triu_indices_from(corr, k=1)
    flat = abs_corr[triu]
    stats = {
        "corr_abs_max": float(flat.max()),
        "corr_abs_p99": float(np.quantile(flat, 0.99)),
        "corr_abs_p90": float(np.quantile(flat, 0.90)),
        "corr_abs_p50": float(np.quantile(flat, 0.50)),
        "corr_abs_mean": float(flat.mean()),
        "frac_corr_abs_gt_005": float((flat > 0.05).mean()),
        "frac_corr_abs_gt_010": float((flat > 0.10).mean()),
        "frac_corr_abs_gt_025": float((flat > 0.25).mean()),
    }

    # Top 20 pares por |corr|
    n_l = len(livres)
    i_idx, j_idx = triu
    top = np.argsort(flat)[-20:][::-1]
    pares_top = []
    for k in top:
        i, j = int(i_idx[k]), int(j_idx[k])
        pares_top.append(
            {
                "e1": int(livres[i]),
                "e2": int(livres[j]),
                "corr": float(corr[i, j]),
                "p_e1": float(f[i]),
                "p_e2": float(f[j]),
                "p_coocorr": float((inc_l[:, i] * inc_l[:, j]).mean()),
            }
        )

    return stats, pares_top


def main() -> None:
    dados = carregar_dados()
    eng = construir_residual(dados)
    snap = eng.snapshot()

    livres = snap["free"]
    livres_set = set(livres)
    print(f"Pós-propagação: fixadas_1={snap['n_fixadas_1']}, "
          f"fixadas_0={snap['n_fixadas_0']}, FREE={snap['n_free']}")

    # 2.1 Estrutura das livres
    edges_uv = dados["edges_uv"]
    freq = dados["freq"]
    estrutura = []
    cls_borda = {0, 5, 30, 35}
    for e in livres:
        u, v = edges_uv[e]
        ru, cu = u // 6, u % 6
        rv, cv = v // 6, v % 6
        # classe simples: interior se ambos extremos não estão na borda externa
        def cls(idx):
            r, c = idx // 6, idx % 6
            if (r in (0, 5)) and (c in (0, 5)):
                return "canto"
            if r in (0, 5) or c in (0, 5):
                return "borda"
            return "interior"
        estrutura.append(
            {
                "idx": int(e),
                "label": rotulo_aresta((u, v)),
                "u": int(u),
                "v": int(v),
                "freq": float(freq[e]),
                "cls_u": cls(u),
                "cls_v": cls(v),
            }
        )

    # 2.2 Grafo de dependência residual
    _, componentes = grafo_dependencia(eng, dados, livres_set)
    print(f"Componentes do grafo residual: {len(componentes)}")
    print("Tamanhos: " + ", ".join(str(len(c)) for c in componentes[:10]))

    # 2.3 Distribuição condicional
    inc = dados["incidencia"]
    stats, top_pares = distribuicao_condicional(inc, livres, freq)
    print("Stats de correlação (entre livres):")
    for k, v in stats.items():
        print(f"  {k}: {v:.4f}")

    # 2.4 Bound do espaço residual
    n_free = snap["n_free"]
    bound_2n = 2 ** n_free
    print(f"\nBound 2^n_free = 2^{n_free} = {bound_2n:.3e}")
    print(f"Ground truth tours = 9862")
    print(f"Densidade = 9862 / 2^{n_free} = {9862 / bound_2n:.3e}")

    # Bound refinado por componente: enumerar atribuições consistentes
    # com as restrições residuais (pares, triplas, quadras, XOR, grau-2)
    # — só faz sentido se a componente for pequena.
    comp_info = []
    for ci, comp in enumerate(componentes):
        info = {"componente": ci, "tamanho": len(comp), "arestas": comp[:50]}
        comp_info.append(info)

    saida = {
        "n_fixadas_1": snap["n_fixadas_1"],
        "fixadas_1": snap["fixadas_1"],
        "n_fixadas_0": snap["n_fixadas_0"],
        "fixadas_0": snap["fixadas_0"],
        "n_free": n_free,
        "free_edges": livres,
        "estrutura_livres": estrutura,
        "n_componentes": len(componentes),
        "tamanhos_componentes": [len(c) for c in componentes],
        "componentes": comp_info,
        "correlacoes_stats": stats,
        "top_pares_corr": top_pares,
        "upper_bound_residual_2n": bound_2n,
        "ground_truth_tours": 9862,
        "densidade_bruta": 9862 / bound_2n,
    }

    out = ROOT / "data" / "residual_variables.json"
    out.write_text(json.dumps(saida, indent=2))
    print(f"\nSalvo em {out}")


if __name__ == "__main__":
    main()
