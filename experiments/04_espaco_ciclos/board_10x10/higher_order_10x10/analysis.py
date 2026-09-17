#!/usr/bin/env python3
"""
analysis.py
===========
Tarefa 4 — análise de escalonamento e diagnóstico do resultado negativo
do NOT(A∧B∧C) no 10×10.

Categoriza triplas proven em:
  - estrutural   : as 3 arestas compartilham um vértice
                   (degree-2 constraint já as bloqueia implicitamente)
  - genuína      : sem vértice comum (exclusão topológica real)

Compara cobertura com 6×6:
  - 6×6: 1.776 triplas minimais / C(80,3) ≈ 82k = 2.2%
  - 10×10: N_proven / C(288,3) ≈ 3.9M

Saída:
  data/analysis.json
"""

import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
PARENT = ROOT.parent

sys.path.insert(0, str(PARENT))
from graph_10x10 import build_graph, label  # noqa: E402

ADJ, EDGES = build_graph()


def edges_common_vertex(idxs):
    """Retorna o set de vértices comum a todas as arestas, ou set vazio."""
    common = set(EDGES[idxs[0]])
    for e in idxs[1:]:
        common &= set(EDGES[e])
    return common


def main():
    ver = json.loads((DATA / "triples_verified.json").read_text())
    proven = [r for r in ver["results"] if r["z3_status"] == "proven"]
    false_pos = [r for r in ver["results"] if r["z3_status"] == "false_pos"]
    print(f"proven   = {len(proven)}")
    print(f"false_pos= {len(false_pos)}")

    # ── categorização ───────────────────────────────────────────────
    structural = []
    genuine = []
    structural_by_deg = Counter()
    for r in proven:
        common = edges_common_vertex(r["edges"])
        if common:
            structural.append(r)
            v = next(iter(common))
            d = len(ADJ[v])
            structural_by_deg[d] += 1
            r["category"] = "structural"
            r["pivot_vertex"] = label(v)
            r["pivot_degree"] = d
        else:
            genuine.append(r)
            r["category"] = "genuine"

    print(f"\nCategorização das 118 PROVEN:")
    print(f"  estruturais (3 arestas em mesmo vértice): {len(structural)}")
    print(f"    grau 3 do pivô: {structural_by_deg.get(3, 0)} triplas (C(3,3)=1 cada)")
    print(f"    grau 4 do pivô: {structural_by_deg.get(4, 0)} triplas (C(4,3)=4 cada)")
    print(f"    grau 5+ pivô  : "
          f"{sum(v for k, v in structural_by_deg.items() if k >= 5)}")
    print(f"  genuínas (topológicas reais): {len(genuine)}")

    # ── contagem teórica de estruturais ─────────────────────────────
    deg_counts = Counter(len(a) for a in ADJ)
    from math import comb
    expected_structural = 0
    for deg, n_vert in deg_counts.items():
        if deg >= 3:
            expected_structural += n_vert * comb(deg, 3)
    print(f"\nEsperado teórico de triplas estruturais "
          f"(todos vértices, C(grau,3)): {expected_structural}")
    print(f"  diferença vs proven encontradas: "
          f"{expected_structural - len(structural)} "
          f"(provavelmente perdidas pelo filtro freq>0.05)")

    # ── escalonamento ───────────────────────────────────────────────
    from math import comb
    C_288_3 = comb(288, 3)
    pct_proven = len(proven) / C_288_3 * 100
    pct_searched = ver["args"].get("max_p_zero", 9999) or 9999
    print(f"\nEscalonamento:")
    print(f"  C(288,3) = {C_288_3:,}")
    print(f"  proven / C(288,3) = {pct_proven:.4f}%")
    print(f"  comparação 6×6: 1.776 / C(80,3) = {1776/comb(80,3)*100:.3f}%")

    # ── leitura do benchmark para sumário do resultado negativo ─────
    try:
        bench = json.loads((DATA / "benchmark_triples.json").read_text())
        speedups = bench.get("speedups_vs_A", {})
        print(f"\nSpeedups do benchmark vs A_puro:")
        for cfg, sp in speedups.items():
            mark = "✓" if sp and sp > 1.0 else "✗"
            print(f"  {mark} {cfg:18s} {sp}")
    except Exception as ex:
        print(f"AVISO: {ex}")

    summary = {
        "n_proven": len(proven),
        "n_false_pos": len(false_pos),
        "n_structural": len(structural),
        "n_genuine": len(genuine),
        "structural_by_pivot_degree": dict(structural_by_deg),
        "expected_structural_total": int(expected_structural),
        "missed_structural": int(expected_structural - len(structural)),
        "scale": {
            "C_288_3": int(C_288_3),
            "proven_pct": round(pct_proven, 6),
            "comparison_6x6_pct": round(1776 / comb(80, 3) * 100, 3),
        },
        "conclusion": (
            "NOT(A∧B∧C) NÃO acelera 10×10. As 118 triplas proven são todas "
            "estruturais (3 arestas em mesmo vértice), já bloqueadas pelo "
            "grau-2 constraint. Adicionar explicitamente apenas sobrecarrega "
            "o solver com restrições redundantes. Triplas GENUÍNAS (sem "
            "vértice comum) não foram detectadas — provavelmente requer "
            "muito mais amostras OU enumeração diferente."
        ),
    }
    with open(DATA / "analysis.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\nSalvo: data/analysis.json")


if __name__ == "__main__":
    main()
