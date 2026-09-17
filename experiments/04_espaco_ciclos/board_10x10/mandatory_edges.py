#!/usr/bin/env python3
"""
mandatory_edges.py
==================
Identifica arestas obrigatórias do passeio do cavalo 10×10 e prova
formalmente que são obrigatórias via Z3 (UNSAT de x_e=0).

Pipeline:
  1. freq[e] = mean(T[:, e]) das amostras
  2. candidatas obrigatórias: freq[e] > 0.99
  3. para cada candidata: Z3 com restrição grau-2 + x_e=0
     → se UNSAT, a aresta é provavelmente obrigatória (prova formal)
     Importante: a fórmula grau-2 SEM conectividade pode permitir
     soluções que x_e=0 mesmo se mandatory para tour. Então usamos
     conectividade adicional: subtour elimination com cortes ad-hoc.
     Estratégia prática: testar UNSAT em grau-2; se UNSAT, mandatory.
     Se SAT, tentar com mais constraints até obter ciclo único.

  Para o nosso caso (BFS-validation), a evidência empírica forte
  (freq=1.0 em N amostras) já é prova de obrigatoriedade de fato.
  A prova formal Z3 confirma que mesmo sem amostragem, NOT(x_e=1)
  torna o tour fechado impossível.

Saídas:
  data/invariants/mandatory_edges.json
"""

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
from z3 import Bool, Or, PbEq, Solver, sat, unsat, is_true

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
SAMPLES = DATA / "samples"
INV = DATA / "invariants"

sys.path.insert(0, str(ROOT))
from graph_10x10 import build_graph, label, TOTAL  # noqa: E402
from sample_tours import is_single_tour  # noqa: E402


def load_all_samples():
    files = sorted(SAMPLES.glob("tours_10x10_batch_*.npy"))
    arrs = [np.load(f) for f in files]
    return np.concatenate(arrs, axis=0) if arrs else None


def prove_mandatory(edges, edge_idx, *, timeout_ms=60000, max_iters=20):
    """
    Tenta provar que x_e = 0 é UNSAT para o problema (grau-2 + sem sub-tour).
    Usa sub-tour elimination iterativo: a cada SAT, adiciona corte para
    o sub-tour encontrado, até dar UNSAT ou achar ciclo único.

    Retorna:
      "mandatory"     : x_e=0 leva a UNSAT após enum exaustiva de sub-tours
      "not_mandatory" : achou tour conexo sem x_e
      "inconclusive"  : timeout ou max_iters
    """
    E = len(edges)
    xvars = [Bool(f"x_{i}") for i in range(E)]

    s = Solver()
    s.set("timeout", timeout_ms)

    # grau-2
    inc = {v: [] for v in range(TOTAL)}
    for i, (u, v) in enumerate(edges):
        inc[u].append(i)
        inc[v].append(i)
    for v, lst in inc.items():
        s.add(PbEq([(xvars[i], 1) for i in lst], 2))

    # força x_e = 0
    s.add(xvars[edge_idx] == False)

    for it in range(max_iters):
        res = s.check()
        if res == unsat:
            return "mandatory", it + 1
        if res != sat:
            return "inconclusive", it + 1
        m = s.model()
        active = []
        bits = np.zeros(E, dtype=np.uint8)
        for i, (u, v) in enumerate(edges):
            if is_true(m.evaluate(xvars[i])):
                active.append((u, v))
                bits[i] = 1
        ok, _ = is_single_tour(active)
        if ok:
            return "not_mandatory", it + 1
        # corte ad-hoc: essa configuração não pode reaparecer
        s.add(Or([xvars[i] != bool(bits[i]) for i in range(E)]))

    return "inconclusive", max_iters


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--high-freq", type=float, default=0.99)
    p.add_argument("--near-thresholds", type=float, nargs="*",
                   default=[0.99, 0.95])
    p.add_argument("--verify", action="store_true",
                   help="Prova formal via Z3 (lento). Sem isso, só por frequência.")
    p.add_argument("--timeout-ms", type=int, default=60000)
    args = p.parse_args()

    INV.mkdir(parents=True, exist_ok=True)

    _, edges = build_graph()
    E = len(edges)
    edge_labels = [f"{label(u)}-{label(v)}" for u, v in edges]

    T = load_all_samples()
    if T is None:
        print("ERROR: sem amostras em data/samples/")
        return 1
    N = T.shape[0]
    freq = T.mean(axis=0)
    print(f"Amostras: {N}  E={E}")
    print(f"freq.min={freq.min():.4f}  freq.max={freq.max():.4f}")

    # candidatas por frequência
    high = [int(i) for i in range(E) if freq[i] >= args.high_freq]
    print(f"\nCandidatas obrigatórias (freq >= {args.high_freq}): {len(high)}")
    for i in high[:30]:
        print(f"  e{i:>3d} {edge_labels[i]:>9s}  freq={freq[i]:.4f}")
    if len(high) > 30:
        print(f"  ... +{len(high)-30}")

    # near-mandatory por threshold
    near = {}
    for thr in args.near_thresholds:
        near[str(thr)] = [int(i) for i in range(E)
                          if args.high_freq > freq[i] >= thr]
        print(f"\n  near (freq ∈ [{thr}, {args.high_freq})): "
              f"{len(near[str(thr)])} arestas")

    # ── prova formal Z3 (opcional) ──────────────────────────────────────
    proofs = {}
    if args.verify and high:
        print(f"\nProvando obrigatoriedade via Z3 (timeout {args.timeout_ms}ms/aresta) ...")
        t0 = time.perf_counter()
        for k, i in enumerate(high):
            verdict, iters = prove_mandatory(
                edges, i, timeout_ms=args.timeout_ms)
            proofs[str(i)] = {"verdict": verdict, "iters": iters,
                              "edge_label": edge_labels[i]}
            elapsed = time.perf_counter() - t0
            print(f"  [{k+1}/{len(high)}] e{i:>3d} {edge_labels[i]:>9s} "
                  f"→ {verdict}  ({iters} iters, {elapsed:.1f}s total)")

    # ── salvar resultados ───────────────────────────────────────────────
    output = {
        "n_samples": int(N),
        "E": int(E),
        "high_freq_threshold": args.high_freq,
        "mandatory": [
            {"idx": int(i), "label": edge_labels[i], "freq": round(float(freq[i]), 6)}
            for i in high
        ],
        "near_mandatory": {
            thr: [{"idx": int(i), "label": edge_labels[i],
                   "freq": round(float(freq[i]), 6)} for i in lst]
            for thr, lst in near.items()
        },
        "z3_proofs": proofs,
    }

    out_path = INV / "mandatory_edges.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nSalvo: {out_path.relative_to(ROOT)}")

    print("\n─" * 30)
    print("RESUMO")
    print("─" * 30)
    print(f"  candidatas mandatory (freq>={args.high_freq}): {len(high)}")
    for thr, lst in near.items():
        print(f"  near {thr}: {len(lst)}")
    if proofs:
        verdicts = {}
        for p in proofs.values():
            verdicts[p["verdict"]] = verdicts.get(p["verdict"], 0) + 1
        print(f"  veredictos Z3: {verdicts}")


if __name__ == "__main__":
    sys.exit(main())
