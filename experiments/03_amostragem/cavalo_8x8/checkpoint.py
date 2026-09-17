"""
checkpoint.py
=============
Estado persistente da sessão de amostragem.

Estrutura do checkpoint:
{
  "config": {...},                  # snapshot do config no início
  "started_at":  ISO 8601 UTC,
  "updated_at":  ISO 8601 UTC,
  "completed_batches": ["p00_b000", ...],   # task_ids concluídos
  "totals": {
    "n_attempts": int,              # somatórios para diagnóstico
    "n_found":    int,
    "elapsed_s":  float,
  },
}

Escrita atômica: write to .tmp + os.replace, evita corrupção em interrupção.
"""

import json
import os
import tempfile
from datetime import datetime, timezone


def now_iso():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def save_checkpoint(state, path):
    state = dict(state)
    state["updated_at"] = now_iso()
    parent = os.path.dirname(path) or "."
    os.makedirs(parent, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=parent, prefix=".checkpoint_", suffix=".tmp")
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as f:
            json.dump(state, f, indent=2, ensure_ascii=False)
        os.replace(tmp, path)
    except Exception:
        if os.path.exists(tmp):
            os.unlink(tmp)
        raise


def load_checkpoint(path):
    if not os.path.exists(path):
        return None
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def init_state(config):
    ts = now_iso()
    return {
        "config": config,
        "started_at": ts,
        "updated_at": ts,
        "completed_batches": [],
        "totals": {"n_attempts": 0, "n_found": 0, "elapsed_s": 0.0},
    }
