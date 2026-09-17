#!/usr/bin/env python3
"""
kl_divergence.py
================
Mede o viés de amostragem do Z3 vs ground truth exaustivo no 6×6.

Importante:
  O sampler Z3 produz CAMINHOS hamiltonianos start→end (via aresta virtual
  no espaço de ciclos GF(2)). O catálogo destruction_catalogue.json contém
  CICLOS fechados, que em geral não correspondem aos caminhos do Z3 (a
  aresta {start,end} não é uma aresta do cavalo para a maioria dos pares).

  Portanto o ground truth para a KL é gerado pelo benchmark/gt_paths_6x6.py
  (enumeração BT específica de caminhos start→end), e este script consome
  esses arquivos gt_paths_6x6_<start>_<end>.json.

KL:
  P_z3(x)      = freq da assinatura x nas amostras Z3
  P_uniform(x) = 1/|S_gt| se x ∈ S_gt, 0 senão
  KL_partial   = sum_{x ∈ S_gt} P_z3(x) * log(P_z3(x) / P_uniform(x))

  Reportamos também:
    - fração das amostras Z3 dentro do suporte do GT
    - cobertura: quantas assinaturas distintas do GT o Z3 visitou

Uso:
  # 1. enumerar GT para cada par usado
  python benchmark/gt_paths_6x6.py --start 0 0 --end 0 5

  # 2. rodar KL contra Z3 fresco
  python benchmark/kl_divergence.py --source benchmark/results/benchmark_z3_fresh.json

Saída:
  benchmark/results/kl_divergence_6x6.json
"""

import argparse
import json
import math
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from cavalo_loop_destruicao_6x6 import vid, label, BOARD  # noqa: E402

Z3_FRESH = ROOT / "benchmark" / "results" / "benchmark_z3_fresh.json"
GT_DIR = ROOT / "benchmark" / "results"
OUT_PATH = ROOT / "benchmark" / "results" / "kl_divergence_6x6.json"


def load_gt_for_pair(start_rc, end_rc):
    """Carrega o arquivo gt_paths_6x6_<start>_<end>.json se existir."""
    s_lbl = label(vid(*start_rc))
    e_lbl = label(vid(*end_rc))
    path = GT_DIR / f"gt_paths_6x6_{s_lbl}_{e_lbl}.json"
    if not path.exists():
        return None, path
    data = json.loads(path.read_text())
    sig_counter = Counter()
    for entry in data["signature_counts"]:
        sig_counter[tuple(bool(b) for b in entry["sig"])] = entry["count"]
    return {
        "n_paths_total": data["n_paths_total"],
        "n_unique_signatures": data["n_unique_signatures"],
        "sig_counter": sig_counter,
        "source_file": path.name,
    }, path


def compute_kl(z3_signatures, gt_sig_counter):
    """KL(P_z3 || P_gt_uniform). Assinaturas fora do suporte do GT contadas à parte."""
    z3_counts = Counter(z3_signatures)
    n_z3 = sum(z3_counts.values())
    gt_support = set(gt_sig_counter.keys())
    n_gt_unique = len(gt_support)

    in_support = 0
    out_support = 0
    kl_partial = 0.0
    for sig, c in z3_counts.items():
        p_z3 = c / n_z3
        if sig in gt_support:
            in_support += c
            p_gt = 1.0 / n_gt_unique
            kl_partial += p_z3 * math.log(p_z3 / p_gt)
        else:
            out_support += c

    visited_unique = len(set(z3_counts.keys()) & gt_support)
    return {
        "n_z3_samples": n_z3,
        "n_z3_unique": len(z3_counts),
        "n_gt_unique": n_gt_unique,
        "in_support": in_support,
        "out_support": out_support,
        "fraction_in_support": round(in_support / n_z3, 6) if n_z3 else None,
        "kl_partial_nats": round(kl_partial, 6),
        "kl_partial_bits": round(kl_partial / math.log(2), 6),
        "coverage_gt": round(visited_unique / n_gt_unique, 6) if n_gt_unique else None,
    }


def analyze_run(run, gt_sig_counter):
    sigs = run.get("signatures", [])
    if not sigs:
        return None
    z3_sigs = [tuple(bool(b) for b in s) for s in sigs]
    out = compute_kl(z3_sigs, gt_sig_counter)
    out["start"] = run["start"]
    out["end"] = run["end"]
    out["target_k"] = run["target_k"]
    out["found_k"] = run["found_k"]
    out["method"] = run["method"]
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--source", action="append",
                   help="JSON com runs Z3 (pode repetir; default fresco+invariantes)")
    p.add_argument("--min-k", type=int, default=1000,
                   help="Só analisar runs com found_k >= min-k")
    args = p.parse_args()

    if not args.source:
        args.source = [str(Z3_FRESH)]
        inv_path = ROOT / "benchmark" / "results" / "benchmark_z3_fresh_invariants.json"
        if inv_path.exists():
            args.source.append(str(inv_path))

    all_results = []
    gt_cache = {}

    for src_str in args.source:
        src = Path(src_str)
        if not src.exists():
            print(f"  AVISO: {src} não existe — pulando")
            continue
        print(f"\nCarregando runs Z3 de {src.relative_to(ROOT)} ...")
        z3_data = json.loads(src.read_text())

        for pair_block in z3_data["pairs"]:
            key = (tuple(pair_block["start"]), tuple(pair_block["end"]))
            for run in pair_block["runs"]:
                if run.get("found_k", 0) < args.min_k:
                    continue
                if key not in gt_cache:
                    gt, gt_path = load_gt_for_pair(*key)
                    if gt is None:
                        print(f"  AVISO: GT ausente para {key} "
                              f"(rode: python benchmark/gt_paths_6x6.py "
                              f"--start {key[0][0]} {key[0][1]} "
                              f"--end {key[1][0]} {key[1][1]})")
                        gt_cache[key] = None
                    else:
                        print(f"  GT {key}: {gt['n_paths_total']:,} caminhos, "
                              f"{gt['n_unique_signatures']:,} assinaturas distintas "
                              f"({gt['source_file']})")
                        gt_cache[key] = gt
                gt = gt_cache[key]
                if gt is None:
                    continue
                r = analyze_run(run, gt["sig_counter"])
                if r is None:
                    continue
                r["n_gt_unique"] = gt["n_unique_signatures"]
                r["n_gt_paths"] = gt["n_paths_total"]
                all_results.append(r)
                print(f"    K={r['target_k']:>5} {r['method']:25s}  "
                      f"KL={r['kl_partial_bits']:.4f} bits  "
                      f"cobertura={r['coverage_gt']:.3%}  "
                      f"fora_sup={r['out_support']}/{r['n_z3_samples']}")

    if not all_results:
        print("\nNenhum resultado calculado. "
              "Verifique que gt_paths_6x6_*.json existem e que Z3 fresco salvou as assinaturas.")
        return 1

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_PATH, "w") as f:
        json.dump({
            "sources": args.source,
            "min_k_filter": args.min_k,
            "results": all_results,
        }, f, indent=2)
    print(f"\nSalvo: {OUT_PATH.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
