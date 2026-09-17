#!/usr/bin/env python3
"""
parallel_runner.py
==================
Driver paralelo de amostragem 8×8 via multiprocessing.Pool.

Características:
  - N workers Z3 independentes (default 8)
  - Master único faz I/O do checkpoint (atomic, sem race condition)
  - Log persistente: data_parallel/logs/run_YYYYMMDD_HHMMSS.log
    com symlink data_parallel/logs/current.log apontando p/ o último.
    Linhas line-buffered, flushadas a cada task → tail -f friendly.
  - SIGINT (Ctrl+C) finaliza pool limpo; checkpoint preservado.
  - Resume: pula task_ids já no checkpoint.

Uso:
  # default: 8 workers, 8 pares canônicos, 100 lotes × 200 amostras por par
  nohup ../.venv/bin/python -u parallel_runner.py > /tmp/8x8_pr.out 2>&1 &
  disown
  tail -f data_parallel/logs/current.log

  # retomar
  python parallel_runner.py --resume

  # ajustar volume
  python parallel_runner.py --n-batches-per-pair 500 --samples-per-batch 200
"""

import argparse
import json
import os
import signal
import sys
import time
from datetime import datetime
from multiprocessing import Pool

from checkpoint import save_checkpoint, load_checkpoint, init_state


# 8 pares canônicos D₄-distintos, paridade validada.
DEFAULT_CONFIG = {
    "board": 8,
    "samples_per_batch": 200,
    "n_batches_per_pair": 100,
    "max_attempts_factor": 5,
    "timeout_s_per_batch": 1800,
    "break_symmetry": True,
    "pairs": [
        [[0, 0], [0, 7]],   # 0: canto → canto (mesma borda)
        [[0, 0], [7, 6]],   # 1: canto → canto (borda diagonal)
        [[0, 0], [4, 3]],   # 2: canto → quase-centro
        [[0, 0], [2, 5]],   # 3: canto → interior próximo
        [[0, 1], [6, 4]],   # 4: borda → interior
        [[0, 1], [7, 5]],   # 5: borda → borda lateral oposta
        [[0, 1], [4, 4]],   # 6: borda → centro
        [[1, 2], [7, 5]],   # 7: interior → borda
    ],
}


# ── utilidades ───────────────────────────────────────────────────────

def parity_ok(s, e):
    return ((s[0] + s[1]) % 2) != ((e[0] + e[1]) % 2)

def task_id(pair_idx, batch_idx):
    return f"p{pair_idx:02d}_b{batch_idx:04d}"

def write_batch(out_dir, tid, payload):
    path = os.path.join(out_dir, "samples", f"batch_{tid}.json")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f)
    return path


# ── worker ───────────────────────────────────────────────────────────

def worker_run(args):
    """
    Executado em processo separado (Pool). Importa Z3/engine só aqui
    para que o master não carregue Z3, evitando interações com fork.
    Cada worker recebe seed determinístico para diversificar amostras.
    """
    pair_idx, batch_idx, start, end, sps, max_att, timeout, bsym, board = args
    tid = task_id(pair_idx, batch_idx)
    try:
        from z3 import set_param
        seed = (pair_idx * 1009 + batch_idx) % (2 ** 31)
        set_param("smt.random_seed", seed)
        set_param("sat.random_seed", seed)

        if board == 6:
            from engine_6x6_sampler import sample_signatures
        elif board == 10:
            from engine_10x10 import sample_signatures
        else:
            from engine_8x8 import sample_signatures
        result = sample_signatures(
            start=tuple(start), end=tuple(end), n_target=sps,
            max_attempts=max_att, timeout_s=timeout,
            break_symmetry=bsym, verbose=False,
        )
        result["task_id"] = tid
        result["pair_idx"] = pair_idx
        result["batch_idx"] = batch_idx
        result["seed"] = seed
        return result
    except Exception as exc:
        return {
            "task_id": tid, "pair_idx": pair_idx, "batch_idx": batch_idx,
            "error": repr(exc),
            "n_found": 0, "n_attempts": 0, "elapsed_s": 0.0,
            "exhausted": False, "timed_out": False,
        }


# ── logger persistente ──────────────────────────────────────────────

class TimestampedLogger:
    def __init__(self, log_path, current_link):
        os.makedirs(os.path.dirname(log_path), exist_ok=True)
        self.path = log_path
        # line-buffered (buffering=1) para que tail -f veja em tempo real
        self.f = open(log_path, "a", buffering=1, encoding="utf-8")
        try:
            if os.path.lexists(current_link):
                os.unlink(current_link)
            os.symlink(os.path.basename(log_path), current_link)
        except OSError:
            pass

    def log(self, msg):
        ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        line = f"[{ts}] {msg}"
        print(line, flush=True)
        self.f.write(line + "\n")

    def close(self):
        try:
            self.f.close()
        except Exception:
            pass


def fmt_eta(seconds):
    if seconds < 60:    return f"{seconds:.0f}s"
    if seconds < 3600:  return f"{seconds/60:.1f}m"
    if seconds < 86400: return f"{seconds/3600:.1f}h"
    return f"{seconds/86400:.1f}d"


# ── main ─────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description="Driver paralelo (multiprocessing) de amostragem 8×8."
    )
    ap.add_argument("--out-dir", default="data_parallel",
                    help="diretório raiz de saída (default: data_parallel)")
    ap.add_argument("--config", default=None,
                    help="caminho para config.json customizado")
    ap.add_argument("--workers", type=int, default=8,
                    help="número de processos paralelos (default: 8)")
    ap.add_argument("--resume", action="store_true",
                    help="retomar do checkpoint existente")
    ap.add_argument("--samples-per-batch", type=int, default=None)
    ap.add_argument("--n-batches-per-pair", type=int, default=None)
    ap.add_argument("--pairs-subset", default=None,
                    help="índices dos pares a processar, ex: '0,1' (default: todos)")
    args = ap.parse_args()

    out_dir = os.path.abspath(args.out_dir)
    cfg_path = os.path.join(out_dir, "config.json")
    chk_path = os.path.join(out_dir, "checkpoint.json")
    logs_dir = os.path.join(out_dir, "logs")
    os.makedirs(logs_dir, exist_ok=True)

    log_name = f"run_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    logger = TimestampedLogger(
        os.path.join(logs_dir, log_name),
        os.path.join(logs_dir, "current.log"),
    )

    # Config: explícito > resume > default
    if args.config:
        with open(args.config, "r", encoding="utf-8") as f:
            base_config = json.load(f)
    elif args.resume and os.path.exists(cfg_path):
        with open(cfg_path, "r", encoding="utf-8") as f:
            base_config = json.load(f)
    else:
        base_config = dict(DEFAULT_CONFIG)

    config = dict(base_config)
    if args.samples_per_batch is not None:
        config["samples_per_batch"] = args.samples_per_batch
    if args.n_batches_per_pair is not None:
        config["n_batches_per_pair"] = args.n_batches_per_pair

    with open(cfg_path, "w", encoding="utf-8") as f:
        json.dump(config, f, indent=2, ensure_ascii=False)

    # Cabeçalho
    logger.log("=" * 70)
    logger.log("PARALLEL_RUNNER 8×8")
    logger.log("=" * 70)
    logger.log(f"out_dir            : {out_dir}")
    logger.log(f"workers            : {args.workers}")
    logger.log(f"samples_per_batch  : {config['samples_per_batch']}")
    logger.log(f"n_batches_per_pair : {config['n_batches_per_pair']}")
    logger.log(f"timeout_s_per_batch: {config.get('timeout_s_per_batch')}")
    logger.log(f"break_symmetry     : {config.get('break_symmetry', True)}")
    logger.log(f"resume             : {args.resume}")
    logger.log("")
    subset = None
    if args.pairs_subset:
        subset = set(int(x) for x in args.pairs_subset.split(","))
        logger.log(f"pairs_subset      : {sorted(subset)}")

    logger.log("Pares canônicos:")
    raw_pairs = [(tuple(s), tuple(e)) for s, e in config["pairs"]]
    pairs = []
    for i, (s, e) in enumerate(raw_pairs):
        if subset is not None and i not in subset:
            logger.log(f"  par {i}: {s} → {e}  SKIP (fora do subset)")
            continue
        if parity_ok(s, e):
            logger.log(f"  par {i}: {s} → {e}  OK")
            pairs.append((i, s, e))
        else:
            logger.log(f"  par {i}: {s} → {e}  PARIDADE FAIL — pulando")
    if not pairs:
        logger.log("Nenhum par válido. Abortando.")
        logger.close()
        return 2

    # Checkpoint
    state = load_checkpoint(chk_path) if args.resume else None
    if state is None:
        state = init_state(config)
    completed = set(state["completed_batches"])

    # Tasks
    tasks = []
    for pidx, s, e in pairs:
        for bidx in range(config["n_batches_per_pair"]):
            tid = task_id(pidx, bidx)
            if tid in completed:
                continue
            tasks.append((
                pidx, bidx, list(s), list(e),
                config["samples_per_batch"],
                config["samples_per_batch"] * config.get("max_attempts_factor", 5),
                config.get("timeout_s_per_batch"),
                config.get("break_symmetry", True),
                config.get("board", 8),
            ))

    n_total = len(pairs) * config["n_batches_per_pair"]
    logger.log("")
    logger.log(f"Lotes totais: {n_total}  (completos: {len(completed)}, pendentes: {len(tasks)})")
    if not tasks:
        logger.log("Nada a fazer (sessão completa).")
        logger.close()
        return 0

    # Estimativa wall-clock: ~22s/lote em 1 core p/ 100 amostras
    # Escala linear com samples_per_batch e inversamente com workers
    est_per_task = 22.0 * (config["samples_per_batch"] / 100.0)
    est_total = est_per_task * len(tasks) / args.workers
    logger.log(f"Estimativa inicial : ~{fmt_eta(est_total)}  "
               f"({est_per_task:.0f}s/lote × {len(tasks)} lotes / {args.workers} workers)")
    logger.log("=" * 70)

    t0 = time.time()
    n_done = 0
    pool = None

    def graceful_shutdown(signum, frame):
        if pool is not None:
            logger.log(f"Sinal {signum} recebido. Finalizando pool...")
            pool.terminate()

    signal.signal(signal.SIGTERM, graceful_shutdown)

    try:
        with Pool(processes=args.workers) as pool:
            for result in pool.imap_unordered(worker_run, tasks):
                n_done += 1
                tid = result["task_id"]

                if "error" in result:
                    logger.log(f"[{n_done}/{len(tasks)}] {tid}  ERRO: {result['error']}")
                    continue

                payload = {
                    **result,
                    "start": list(result["start"]),
                    "end": list(result["end"]),
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

                wall = time.time() - t0
                rate_h = (n_done / wall) * 3600
                eta = (len(tasks) - n_done) / max(n_done / wall, 1e-9)
                flags = []
                if result["exhausted"]: flags.append("EXHAUSTED")
                if result["timed_out"]: flags.append("TIMED_OUT")
                flag_str = " " + " ".join(f"[{f}]" for f in flags) if flags else ""
                logger.log(
                    f"[{n_done}/{len(tasks)}] {tid}  "
                    f"start={result['start']} end={result['end']}  "
                    f"{result['n_found']}/{result['n_target']}  "
                    f"{result['elapsed_s']:.1f}s  "
                    f"(rate={rate_h:.0f}/h, eta={fmt_eta(eta)}){flag_str}"
                )
    except KeyboardInterrupt:
        logger.log("Ctrl+C recebido. Checkpoint preservado, retomar com --resume.")
    finally:
        elapsed = time.time() - t0
        logger.log("=" * 70)
        logger.log(
            f"FIM  | concluídos nesta sessão: {n_done} | "
            f"wall={fmt_eta(elapsed)} | "
            f"amostras_acumuladas={state['totals']['n_found']:,} | "
            f"checkpoint={os.path.relpath(chk_path, out_dir)}"
        )
        logger.log("=" * 70)
        logger.close()

    return 0


if __name__ == "__main__":
    sys.exit(main())
