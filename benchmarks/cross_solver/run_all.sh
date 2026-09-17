#!/usr/bin/env bash
# Pipeline de benchmark Z3 vs Backtracking 6×6.
# Roda em série (não paralelo) para que cada medição seja em CPU isolada.
set -e

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PY="${ROOT}/.venv/bin/python"
cd "${ROOT}"

echo "═════════════════════════════════════════════════════════════════"
echo "  BENCHMARK 6×6  —  Z3 GF(2)  vs  BACKTRACKING"
echo "═════════════════════════════════════════════════════════════════"
echo

echo "[1/6] BT baseline 6×6 instrumentado ..."
"${PY}" benchmark/baseline_bt_6x6.py 2>&1 | tee benchmark/results/baseline_bt_6x6_run.log
echo

echo "[2/6] Z3 agregação dos batches existentes ..."
"${PY}" benchmark/benchmark_z3.py --aggregate-only
echo

echo "[3/6] Z3 fresco (sym=False, K∈{10,100,1000}) ..."
"${PY}" benchmark/benchmark_z3.py --fresh-only --ks 10 100 1000 --n-pairs 3 \
    2>&1 | tee benchmark/results/z3_fresh_run.log
echo

echo "[4/6] Z3 + invariantes (H4) ..."
"${PY}" benchmark/benchmark_z3.py --fresh-only --invariants \
    --ks 10 100 1000 --n-pairs 3 \
    2>&1 | tee benchmark/results/z3_fresh_invariants_run.log
echo

echo "[5/6] Enumeração GT de caminhos hamiltonianos (3 primeiros pares) ..."
for spec in "0 0 0 5" "0 0 5 4" "0 0 3 2"; do
    set -- $spec
    "${PY}" benchmark/gt_paths_6x6.py --start $1 $2 --end $3 $4 \
        2>&1 | tee -a benchmark/results/gt_paths_run.log
done
echo

echo "[6/6] KL-divergência Z3 vs GT (K≥1000) ..."
"${PY}" benchmark/kl_divergence.py \
    --source benchmark/results/benchmark_z3_fresh.json \
    --source benchmark/results/benchmark_z3_fresh_invariants.json
echo

echo "── relatório final ──"
"${PY}" benchmark/report.py
