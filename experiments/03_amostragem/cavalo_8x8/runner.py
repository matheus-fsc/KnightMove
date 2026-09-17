#!/usr/bin/env python3
"""
runner.py
=========
Driver serial de amostragem 8×8.

Para cada par (start, end) e índice de lote, chama
`engine_8x8.sample_signatures` e salva o resultado em
data/samples/batch_<task_id>.json. Mantém checkpoint atômico
em data/checkpoint.json para retomar após interrupção.

Uso típico:

  # primeira execução (defaults: 2 pares × 2 lotes × 50 amostras = 200)
  python3 runner.py

  # retomar do checkpoint
  python3 runner.py --resume

  # rodar mais agressivo
  python3 runner.py --samples-per-batch 200 --n-batches-per-pair 25 --resume

  # config próprio
  python3 runner.py --config meu_config.json
"""

import argparse
import json
import os
import sys

from engine_8x8 import sample_signatures
from checkpoint import save_checkpoint, load_checkpoint, init_state


# ── defaults pequenos para validação rápida da Fase 1 ────────────────
# Para Fase 3 (full sweep), passar --samples-per-batch e --n-batches-per-pair
# pela linha de comando, ou trocar o config persistido.
DEFAULT_CONFIG = {
    "board": 8,
    "samples_per_batch": 50,
    "n_batches_per_pair": 2,
    "max_attempts_factor": 5,        # max_attempts = factor * samples_per_batch
    "timeout_s_per_batch": 600,
    "break_symmetry": True,
    "pairs": [
        # Pares com paridade oposta (caminho hamiltoniano possível).
        # Lista pequena para Fase 1; cobertura D₄ canônica vem na Fase 3.
        [[0, 0], [0, 7]],   # canto A8 → canto H8
        [[0, 1], [6, 4]],   # borda B8 → interior E2
    ],
}


def parity_ok(start, end):
    return ((start[0] + start[1]) % 2) != ((end[0] + end[1]) % 2)


def task_id(pair_idx, batch_idx):
    return f"p{pair_idx:02d}_b{batch_idx:03d}"


def write_batch(out_dir, tid, payload):
    path = os.path.join(out_dir, "samples", f"batch_{tid}.json")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f)
    return path


def merge_config(base, overrides):
    cfg = dict(base)
    for k, v in overrides.items():
        if v is not None:
            cfg[k] = v
    return cfg


def main():
    ap = argparse.ArgumentParser(
        description="Driver serial de amostragem hamiltoniana 8×8 com checkpoint."
    )
    ap.add_argument("--out-dir", default="data",
                    help="diretório raiz de saída (default: data)")
    ap.add_argument("--config", default=None,
                    help="caminho para config.json customizado")
    ap.add_argument("--resume", action="store_true",
                    help="retomar a partir do checkpoint existente")
    ap.add_argument("--samples-per-batch", type=int, default=None)
    ap.add_argument("--n-batches-per-pair", type=int, default=None)
    ap.add_argument("--timeout-s-per-batch", type=int, default=None)
    ap.add_argument("--verbose", action="store_true",
                    help="repassa verbose=True ao engine")
    args = ap.parse_args()

    out_dir = os.path.abspath(args.out_dir)
    cfg_path = os.path.join(out_dir, "config.json")
    chk_path = os.path.join(out_dir, "checkpoint.json")
    os.makedirs(out_dir, exist_ok=True)

    # Carrega config: explícito > checkpoint existente > default
    if args.config:
        with open(args.config, "r", encoding="utf-8") as f:
            base_config = json.load(f)
    elif args.resume and os.path.exists(cfg_path):
        with open(cfg_path, "r", encoding="utf-8") as f:
            base_config = json.load(f)
    else:
        base_config = dict(DEFAULT_CONFIG)

    config = merge_config(base_config, {
        "samples_per_batch": args.samples_per_batch,
        "n_batches_per_pair": args.n_batches_per_pair,
        "timeout_s_per_batch": args.timeout_s_per_batch,
    })

    # Persiste config canônico (espelha o que de fato foi usado)
    with open(cfg_path, "w", encoding="utf-8") as f:
        json.dump(config, f, indent=2, ensure_ascii=False)

    # Validação de paridade dos pares
    raw_pairs = [(tuple(s), tuple(e)) for s, e in config["pairs"]]
    pairs = []
    print("=" * 70)
    print("RUNNER 8×8 SERIAL")
    print("=" * 70)
    print()
    print("Validação de paridade dos pares (caminho hamiltoniano exige paridade oposta):")
    for i, (s, e) in enumerate(raw_pairs):
        if parity_ok(s, e):
            print(f"  OK    par {i}: {s} → {e}")
            pairs.append((s, e))
        else:
            print(f"  SKIP  par {i}: {s} → {e}  (mesma paridade — descartado)")
    if not pairs:
        print("\nNenhum par válido. Abortando.")
        return 2

    # Carrega ou inicializa checkpoint
    state = load_checkpoint(chk_path) if args.resume else None
    if state is None:
        state = init_state(config)
    completed = set(state["completed_batches"])

    # Constrói lista de tarefas pendentes
    tasks = []
    for pidx, (s, e) in enumerate(pairs):
        for bidx in range(config["n_batches_per_pair"]):
            tid = task_id(pidx, bidx)
            if tid not in completed:
                tasks.append((tid, pidx, bidx, s, e))

    n_total = len(pairs) * config["n_batches_per_pair"]
    print()
    print(f"Configuração da sessão:")
    print(f"  out_dir            : {out_dir}")
    print(f"  pares válidos      : {len(pairs)}")
    print(f"  lotes por par      : {config['n_batches_per_pair']}")
    print(f"  amostras por lote  : {config['samples_per_batch']}")
    print(f"  total de lotes     : {n_total}")
    print(f"  já completos       : {len(completed)}")
    print(f"  pendentes          : {len(tasks)}")
    print(f"  break_symmetry     : {config.get('break_symmetry', True)}")
    print(f"  timeout/batch      : {config.get('timeout_s_per_batch')}s")
    print()

    if not tasks:
        print("Nada a fazer (sessão completa).")
        return 0

    sps = config["samples_per_batch"]
    factor = config.get("max_attempts_factor", 5)
    timeout = config.get("timeout_s_per_batch", None)
    bsym = config.get("break_symmetry", True)

    for k, (tid, pidx, bidx, s, e) in enumerate(tasks):
        print(f"[{k+1}/{len(tasks)}] {tid}  start={s} end={e}  alvo={sps}")
        result = sample_signatures(
            start=s, end=e, n_target=sps,
            max_attempts=factor * sps,
            timeout_s=timeout,
            break_symmetry=bsym,
            verbose=args.verbose,
        )

        payload = {
            "task_id": tid,
            "pair_idx": pidx,
            "batch_idx": bidx,
            **result,
        }
        path = write_batch(out_dir, tid, payload)

        completed.add(tid)
        state["completed_batches"] = sorted(completed)
        state["totals"]["n_attempts"] += result["n_attempts"]
        state["totals"]["n_found"] += result["n_found"]
        state["totals"]["elapsed_s"] = round(
            state["totals"]["elapsed_s"] + result["elapsed_s"], 2
        )
        save_checkpoint(state, chk_path)

        flags = []
        if result["exhausted"]: flags.append("EXHAUSTED")
        if result["timed_out"]: flags.append("TIMED_OUT")
        flag_str = " " + " ".join(f"[{f}]" for f in flags) if flags else ""
        print(f"   → {result['n_found']}/{sps} amostras  "
              f"em {result['elapsed_s']}s  "
              f"(attempts={result['n_attempts']}){flag_str}")
        print(f"   ↳ {os.path.relpath(path, out_dir)}")

    print()
    print("=" * 70)
    print("CONCLUÍDO")
    print(f"  amostras coletadas (sessão) : {state['totals']['n_found']}")
    print(f"  tentativas (sessão)         : {state['totals']['n_attempts']}")
    print(f"  tempo acumulado (sessão)    : {state['totals']['elapsed_s']:.1f}s")
    print(f"  checkpoint                  : {os.path.relpath(chk_path, out_dir)}")
    print("=" * 70)
    return 0


if __name__ == "__main__":
    sys.exit(main())
