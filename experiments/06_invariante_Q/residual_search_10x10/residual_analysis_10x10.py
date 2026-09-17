"""Análise do espaço residual no 10×10."""

from __future__ import annotations

import json
from collections import defaultdict
from itertools import combinations
from pathlib import Path

import numpy as np

from propagation_engine_10x10 import (
    FREE,
    PropagationEngine10,
    ROOT,
    carregar_dados,
)


def grafo_dependencia(eng, dados, livres_set):
    arestas = defaultdict(set)

    def adicionar(vars_, _fonte):
        livres = [e for e in vars_ if eng.estado[e] == FREE]
        if len(livres) < 2:
            return
        for a, b in combinations(livres, 2):
            arestas[a].add(b)
            arestas[b].add(a)

    for p in dados["pares"]:
        adicionar(p, "par")
    for s in dados["xor"]:
        adicionar(s, "xor")
    for v, es in eng.v2e.items():
        adicionar(tuple(es), "vertice")

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


def distribuicao_condicional(inc, livres, freq):
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
    i_idx, j_idx = triu
    top = np.argsort(flat)[-15:][::-1]
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
    dados = carregar_dados(strict_pair_mode="todos")
    eng = PropagationEngine10(
        edges_uv=dados["edges_uv"],
        freq=dados["freq"],
        pares=dados["pares"],
        triplas=dados["triplas"],
        quadras=dados["quadras"],
        xor=dados["xor"],
        xor_paridade=dados["xor_paridade_alvo"],
    )
    eng.run_to_fixpoint([1, 2, 3, 6])
    snap = eng.snapshot()
    livres = snap["free"]
    livres_set = set(livres)
    print(
        f"Pós-propagação: fixadas_1={snap['n_fixadas_1']}, "
        f"fixadas_0={snap['n_fixadas_0']}, FREE={snap['n_free']}"
    )

    # Estrutura espacial das livres
    def cls(idx):
        r, c = idx // 10, idx % 10
        if (r in (0, 9)) and (c in (0, 9)):
            return "canto"
        if r in (0, 9) or c in (0, 9):
            return "borda"
        return "interior"

    edges_uv = dados["edges_uv"]
    estrutura_resumo = defaultdict(int)
    for e in livres:
        u, v = edges_uv[e]
        estrutura_resumo[(cls(u), cls(v))] += 1
    estrutura_resumo = {f"{a}-{b}": int(v) for (a, b), v in estrutura_resumo.items()}

    # Histograma de freq nas livres
    freq = dados["freq"]
    fl = freq[livres]
    bins = np.linspace(0, 1, 11)
    hist, _ = np.histogram(fl, bins=bins)
    print(f"freq das livres: min={fl.min():.3f}, max={fl.max():.3f}, mean={fl.mean():.3f}")

    # Grafo de dependência
    _, componentes = grafo_dependencia(eng, dados, livres_set)
    print(f"Componentes do grafo residual: {len(componentes)}")
    print("Tamanhos: " + ", ".join(str(len(c)) for c in componentes[:8]))

    # Correlações
    stats, top_pares = distribuicao_condicional(dados["amostras"], livres, freq)
    print("Stats de correlação (entre livres):")
    for k, v in stats.items():
        print(f"  {k}: {v:.4f}")

    n_free = snap["n_free"]
    bound = 2.0 ** n_free
    print(f"\nBound ingênuo 2^{n_free} ≈ {bound:.3e}")
    print(f"2-fatores estimados (~10^15) → razão {1e15 / bound:.3e}")

    saida = {
        "n_fixadas_1": snap["n_fixadas_1"],
        "fixadas_1": snap["fixadas_1"],
        "n_fixadas_0": snap["n_fixadas_0"],
        "fixadas_0": snap["fixadas_0"],
        "n_free": n_free,
        "free_edges": livres,
        "estrutura_resumo_classes": estrutura_resumo,
        "freq_livres_hist": hist.tolist(),
        "freq_livres_min": float(fl.min()),
        "freq_livres_max": float(fl.max()),
        "freq_livres_mean": float(fl.mean()),
        "n_componentes": len(componentes),
        "tamanhos_componentes": [len(c) for c in componentes],
        "correlacoes_stats": stats,
        "top_pares_corr": top_pares,
        "upper_bound_2n": bound,
    }
    out = ROOT / "data" / "residual_variables.json"
    out.write_text(json.dumps(saida, indent=2))
    print(f"\nSalvo em {out}")


if __name__ == "__main__":
    main()
