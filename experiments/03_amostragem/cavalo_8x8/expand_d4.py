"""
expand_d4.py
============
Expansão D₄ das amostras existentes — sem rodar Z3 de novo.

A ideia:
  Toda transformação T ∈ D₄ preserva o grafo do cavalo (knight moves
  são invariantes sob rotações e reflexões do tabuleiro). Logo, se P é
  um caminho hamiltoniano de start → end, então T(P) é um caminho
  hamiltoniano de T(start) → T(end), com as mesmas propriedades
  topológicas, apenas com vértices/arestas reetiquetados por T.

  Aplicando todas as 8 transformações de D₄ a cada amostra original,
  geramos 8N amostras que cobrem toda a órbita D₄ do par (start, end).
  Quando agrupadas no "canonical pool" (lex-min do par), as estatísticas
  do pool ficam D₄-invariantes por construção — exatamente o que era
  natural no 6×6 fechado.

API:
  transform_batch(batch, t)            — gera o batch após D₄ transform t
  canonical_pair(start, end)           — lex-min representante D₄ do par
  expand_with_d4(batches)              — gera 8 imagens por batch
  expand_and_pool_by_canonical(batches)— agrupa por pool canônico
"""

import glob
import json
import os
from collections import defaultdict

from d4_orbits import (
    label_to_edge, edge_to_label,
    d4_apply_edge, d4_apply_node,
)


def transform_batch(batch, t):
    """
    Aplica a t-ésima transformação D₄ a um batch.
    Reetiqueta arestas (signature) e (start, end), preservando a ordem
    canônica das arestas (sorted) usada pelo engine.
    """
    edge_labels = batch["edges"]
    edge_tuples = [label_to_edge(s) for s in edge_labels]
    label_to_idx = {edge_to_label(e): k for k, e in enumerate(edge_tuples)}

    # remap[i] = j  onde  edges[j] = D₄_t(edges[i])
    remap = [None] * len(edge_tuples)
    for i, e in enumerate(edge_tuples):
        e_t = d4_apply_edge(e, t)
        remap[i] = label_to_idx[edge_to_label(e_t)]

    new_signatures = []
    for sig in batch["signatures"]:
        new_sig = [False] * len(sig)
        for i, j in enumerate(remap):
            new_sig[j] = sig[i]
        new_signatures.append(new_sig)

    new_start = list(d4_apply_node(tuple(batch["start"]), t))
    new_end   = list(d4_apply_node(tuple(batch["end"]), t))

    out = {k: v for k, v in batch.items()
           if k not in {"signatures", "start", "end"}}
    out["start"] = new_start
    out["end"] = new_end
    out["signatures"] = new_signatures
    out["d4_t"] = t
    out["task_id"] = batch.get("task_id", "?") + f"_t{t}"
    return out


def canonical_pair(start, end):
    """Lex-min do par não-ordenado {start, end} sobre as 8 imagens D₄."""
    images = set()
    for t in range(8):
        s_t = d4_apply_node(tuple(start), t)
        e_t = d4_apply_node(tuple(end), t)
        images.add(tuple(sorted([s_t, e_t])))
    return min(images)


def expand_with_d4(batches):
    """Lista expandida: 8 × len(batches) com todas as transformações D₄."""
    out = []
    for b in batches:
        for t in range(8):
            out.append(transform_batch(b, t))
    return out


def expand_and_pool_by_canonical(batches):
    """
    Para cada batch, gera as 8 imagens D₄ e agrupa pelo par canônico.
    Retorna dict { canonical_pair_tuple : list[batch] }.
    Cada pool canônico é D₄-invariante por construção.
    """
    pools = defaultdict(list)
    for b in batches:
        canon = canonical_pair(b["start"], b["end"])
        for t in range(8):
            pools[canon].append(transform_batch(b, t))
    return pools


# ── smoke test ───────────────────────────────────────────────────────

if __name__ == "__main__":
    paths = sorted(glob.glob("data/samples/batch_*.json"))
    if not paths:
        print("Sem batches em data/samples/. Rode o runner primeiro.")
        raise SystemExit(1)

    batches = []
    for p in paths:
        with open(p, "r", encoding="utf-8") as f:
            batches.append(json.load(f))

    n_orig = sum(b["n_found"] for b in batches)
    print(f"Lotes carregados: {len(batches)}")
    print(f"Amostras originais: {n_orig:,}")

    pools = expand_and_pool_by_canonical(batches)
    print(f"\nPools canônicos: {len(pools)}")
    for canon, plist in pools.items():
        n = sum(b["n_found"] for b in plist)
        s, e = canon
        s_lbl = chr(ord('A') + s[1]) + str(8 - s[0])
        e_lbl = chr(ord('A') + e[1]) + str(8 - e[0])
        print(f"  canon {s_lbl}→{e_lbl}: {len(plist)} lotes expandidos, {n:,} amostras")

    # Validação: D4 transform preserva nº de arestas ativas em cada signature
    sample = batches[0]["signatures"][0]
    n_active_orig = sum(sample)
    for t in range(8):
        tb = transform_batch(batches[0], t)
        n_active_t = sum(tb["signatures"][0])
        assert n_active_t == n_active_orig, (
            f"D₄_{t} não preservou nº de arestas: {n_active_t} ≠ {n_active_orig}"
        )
    print(f"\n✓ Sanity: 8 transformações preservam {n_active_orig} arestas ativas por signature")
