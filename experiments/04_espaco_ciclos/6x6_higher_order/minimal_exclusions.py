#!/usr/bin/env python3
"""
minimal_exclusions.py
=====================
Busca exclusões minimais de triplas e quadras no espaço de tours 6×6.

Uma tupla (e1, ..., ek) é EXCLUÍDA se P(x_{e1}=1 ∧ ... ∧ x_{ek}=1) = 0.
É MINIMAL se nenhuma subtupla própria já é excluída — isto é, o
"motivo" da exclusão emerge só com a tupla completa.

Para triplas: usa o tensor 3-way P(e1,e2,e3) = (T^T T) ⊗ T / N,
calculado em chunks por e3.

Para quadras: enumera combinações usando o tensor de pares já
calculado para podar (pular se algum subpar é excluído).

Saídas:
  data/exclusions_triples.json
  data/exclusions_quads.json
"""

import argparse
import json
import sys
import time
from pathlib import Path
from itertools import combinations

import numpy as np

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"

sys.path.insert(0, str(ROOT.parent))
from cavalo_loop_destruicao_6x6 import EDGES_LIST, label  # noqa: E402


def compute_triple_tensor(T):
    """
    Tensor P[e1,e2,e3] = freq de coexistência das 3 arestas.
    Implementação: P[:,:,e3] = (T * T[:,e3:e3+1])^T @ (T * T[:,e3:e3+1]) / N
    """
    N, E = T.shape
    T32 = T.astype(np.int32)
    P3 = np.zeros((E, E, E), dtype=np.int32)
    for e3 in range(E):
        mask = T32 * T32[:, e3:e3 + 1]  # (N, E); row 1 only se tour i tem e3
        P3[:, :, e3] = mask.T @ mask
    return P3.astype(np.int64)  # contagens inteiras


def find_minimal_triples(T, freq, P2_count, P3_count, *, freq_threshold=0):
    """
    Triplas com P(∧)=0 e cada subpar com P>0.
    Retorna lista de dicts.
    """
    N, E = T.shape
    live = [e for e in range(E) if freq[e] > freq_threshold]

    out = []
    n_excluded = 0
    n_minimal = 0
    for i in range(len(live)):
        e1 = live[i]
        for j in range(i + 1, len(live)):
            e2 = live[j]
            if P2_count[e1, e2] == 0:
                continue  # já é par excluído — tripla não-minimal
            for k in range(j + 1, len(live)):
                e3 = live[k]
                if P2_count[e1, e3] == 0 or P2_count[e2, e3] == 0:
                    continue  # tripla envolve par excluído
                if P3_count[e1, e2, e3] == 0:
                    n_excluded += 1
                    n_minimal += 1
                    out.append((e1, e2, e3,
                                int(P2_count[e1, e2]),
                                int(P2_count[e1, e3]),
                                int(P2_count[e2, e3])))
    return out, n_excluded, n_minimal


def find_minimal_quads(T, freq, P2_count, P3_count, *, freq_threshold=0,
                       triples_minimal_set=None,
                       sample_if_above=5_000_000,
                       rng_seed=2026):
    """
    Quadras com P(∧)=0 e nenhum subtripla/subpar excluído.
    Implementação vetorizada com bit-packing: T comprimida para uint64
    blocks (N/64 + remainder), AND + popcount = contagem da quadra.
    """
    N, E = T.shape
    live = [e for e in range(E) if freq[e] > freq_threshold]
    n_live = len(live)
    from math import comb
    total_quads = comb(n_live, 4)

    # bit-packing: cada coluna de T vira um vetor de uint64 de tamanho ceil(N/64)
    Nb = (N + 63) // 64
    T_bits = np.zeros((E, Nb), dtype=np.uint64)
    for e in range(E):
        bits = np.packbits(T[:, e][::-1], bitorder='little')
        # packbits dá uint8 com N+7//8; precisamos uint64 com Nb*8 bytes
        padded = np.zeros(Nb * 8, dtype=np.uint8)
        # converte T[:, e] (bool/uint8) → uint64 packed bit-order pequeno-endian
        col = T[:, e].astype(np.uint8)
        for k in range(N):
            if col[k]:
                padded[k // 8] |= np.uint8(1 << (k % 8))
        T_bits[e] = padded.view(np.uint64)

    sampled = False
    if total_quads > sample_if_above:
        sampled = True
        rng = np.random.default_rng(rng_seed)
        all_pairs = []
        seen = set()
        while len(seen) < sample_if_above:
            qs = rng.choice(n_live, size=4, replace=False)
            qs.sort()
            t = tuple(qs.tolist())
            if t not in seen:
                seen.add(t)
        quads_idx = np.array(list(seen), dtype=np.int32)
    else:
        quads_idx = np.array(list(combinations(range(n_live), 4)),
                             dtype=np.int32)
    quads_global = np.array([(live[a], live[b], live[c], live[d])
                              for a, b, c, d in quads_idx], dtype=np.int32)
    M = quads_global.shape[0]

    out = []
    n_excluded = 0
    n_minimal = 0
    n_pair_pruned = 0
    n_triple_pruned = 0
    t0 = time.perf_counter()

    # processamento em chunks vetorizado
    CHUNK = 100_000
    for start in range(0, M, CHUNK):
        chunk = quads_global[start:start + CHUNK]
        # poda por par excluído (vetorizada)
        e1, e2, e3, e4 = chunk[:, 0], chunk[:, 1], chunk[:, 2], chunk[:, 3]
        pair_ok = (
            (P2_count[e1, e2] > 0) &
            (P2_count[e1, e3] > 0) &
            (P2_count[e1, e4] > 0) &
            (P2_count[e2, e3] > 0) &
            (P2_count[e2, e4] > 0) &
            (P2_count[e3, e4] > 0)
        )
        n_pair_pruned += int((~pair_ok).sum())
        kept = chunk[pair_ok]
        if len(kept) == 0:
            continue

        # poda por tripla com P=0 (não-minimal se contém)
        e1, e2, e3, e4 = kept[:, 0], kept[:, 1], kept[:, 2], kept[:, 3]
        triple_ok = (
            (P3_count[e1, e2, e3] > 0) &
            (P3_count[e1, e2, e4] > 0) &
            (P3_count[e1, e3, e4] > 0) &
            (P3_count[e2, e3, e4] > 0)
        )
        n_triple_pruned += int((~triple_ok).sum())
        candidates = kept[triple_ok]

        # avaliar P(quad)=0 nas candidatas
        for q in candidates:
            a, b, c, d = q
            ab = T_bits[a] & T_bits[b]
            cd = T_bits[c] & T_bits[d]
            combined = ab & cd
            # popcount via np.unpackbits trick ou bit_count
            cnt = sum(int(x).bit_count() for x in combined.tolist())
            if cnt == 0:
                n_minimal += 1
                n_excluded += 1
                out.append((int(a), int(b), int(c), int(d)))

        elapsed = time.perf_counter() - t0
        print(f"    [chunk {start//CHUNK+1}/{(M+CHUNK-1)//CHUNK}] "
              f"processed={min(start+CHUNK, M):>9,}/{M:,}  "
              f"pair_pruned={n_pair_pruned:>8,}  "
              f"triple_pruned={n_triple_pruned:>8,}  "
              f"minimal={n_minimal}  elapsed={elapsed:.1f}s",
              flush=True)

    elapsed = time.perf_counter() - t0
    return out, n_excluded, n_minimal, total_quads, sampled, M, elapsed


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--max-triples", type=int, default=None,
                   help="Limita output (default: salva todas)")
    p.add_argument("--skip-quads", action="store_true")
    p.add_argument("--sample-quads-above", type=int, default=5_000_000)
    args = p.parse_args()

    DATA.mkdir(parents=True, exist_ok=True)

    print("Carregando T e frequências ...")
    T = np.load(DATA / "incidence_matrix_6x6.npy")
    freq = np.load(DATA / "freq_singles.npy")
    N, E = T.shape
    print(f"  T shape={T.shape}  E={E}")

    # contagens inteiras (mais robusto que float pra teste ==0)
    print("\nContagens de pares (int)...")
    P2_count = (T.astype(np.int32).T @ T.astype(np.int32))

    edges = list(EDGES_LIST)
    edge_labels = [f"{label(u)}-{label(v)}" for u, v in edges]

    # ── TRIPLAS ──────────────────────────────────────────────────────
    print("\nConstruindo tensor 3-way P3 ...")
    t0 = time.perf_counter()
    P3_count = compute_triple_tensor(T)
    print(f"  shape={P3_count.shape}  tempo={time.perf_counter()-t0:.1f}s")

    print("\nBuscando triplas minimais (P=0, nenhum subpar excluído) ...")
    t0 = time.perf_counter()
    triples, n_excluded_t, n_minimal_t = find_minimal_triples(
        T, freq, P2_count, P3_count)
    print(f"  triplas minimais excluídas: {n_minimal_t}")
    print(f"  tempo: {time.perf_counter()-t0:.1f}s")

    # também contagem total de triplas excluídas (incluindo não-minimais)
    # quaisquer 3 vivas com P3=0
    n_excl_all_t = 0
    n_total_t = 0
    live = [e for e in range(E) if freq[e] > 0]
    for i in range(len(live)):
        for j in range(i + 1, len(live)):
            for k in range(j + 1, len(live)):
                e1, e2, e3 = live[i], live[j], live[k]
                n_total_t += 1
                if P3_count[e1, e2, e3] == 0:
                    n_excl_all_t += 1

    n_nonmin_t = n_excl_all_t - n_minimal_t
    print(f"\nResumo triplas:")
    print(f"  total enumeradas       : {n_total_t}")
    print(f"  excluídas (P=0)        : {n_excl_all_t}")
    print(f"  minimais               : {n_minimal_t}")
    print(f"  reduzem a subpar excl. : {n_nonmin_t}")

    # serializa triplas minimais
    out_triples = []
    for e1, e2, e3, c12, c13, c23 in triples:
        out_triples.append({
            "edges": [int(e1), int(e2), int(e3)],
            "edge_names": [edge_labels[e1], edge_labels[e2], edge_labels[e3]],
            "freq_individual": [round(float(freq[e1]), 4),
                                round(float(freq[e2]), 4),
                                round(float(freq[e3]), 4)],
            "count_pairs": [int(c12), int(c13), int(c23)],
            "freq_pairs": [round(c12 / N, 4), round(c13 / N, 4),
                           round(c23 / N, 4)],
            "minimal": True,
        })

    with open(DATA / "exclusions_triples.json", "w") as f:
        json.dump({
            "n_tours": N,
            "n_edges": E,
            "n_triples_enum": n_total_t,
            "n_triples_excluded": n_excl_all_t,
            "n_triples_minimal": n_minimal_t,
            "n_triples_redundant": n_nonmin_t,
            "minimal_triples": out_triples,
        }, f, indent=2)
    print(f"  Salvo: data/exclusions_triples.json")

    if out_triples[:10]:
        print("\nTop 10 triplas minimais:")
        for t in out_triples[:10]:
            print(f"  {t['edge_names']}  freq={t['freq_individual']}  "
                  f"par_freqs={t['freq_pairs']}")

    if args.skip_quads:
        print("\n[--skip-quads] pulando quadras.")
        return

    # ── QUADRAS ──────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("Buscando quadras minimais ...")
    triples_set = {(t[0], t[1], t[2]) for t in triples}
    quads_res = find_minimal_quads(
        T, freq, P2_count, P3_count,
        triples_minimal_set=triples_set,
        sample_if_above=args.sample_quads_above,
    )
    quads, n_ex_q, n_min_q, total_q, sampled, processed, t_q = quads_res
    print(f"  total quadras (C(E,4))  : {total_q}")
    print(f"  amostragem ativa?       : {sampled}")
    print(f"  quadras processadas     : {processed}")
    print(f"  excluídas (entre as proc.): {n_ex_q}")
    print(f"  minimais (P=0, sem sub-excl.): {n_min_q}")
    print(f"  tempo                   : {t_q:.1f}s")

    out_quads = []
    for q in quads:
        e1, e2, e3, e4 = q
        out_quads.append({
            "edges": [int(e1), int(e2), int(e3), int(e4)],
            "edge_names": [edge_labels[e1], edge_labels[e2],
                           edge_labels[e3], edge_labels[e4]],
            "freq_individual": [round(float(freq[e]), 4) for e in q],
            "minimal": True,
        })

    with open(DATA / "exclusions_quads.json", "w") as f:
        json.dump({
            "n_tours": N,
            "n_edges": E,
            "total_quads_C_E_4": total_q,
            "sampled": bool(sampled),
            "n_processed": processed,
            "n_excluded_processed": n_ex_q,
            "n_minimal": n_min_q,
            "minimal_quads": out_quads,
        }, f, indent=2)
    print(f"  Salvo: data/exclusions_quads.json")

    if out_quads[:5]:
        print("\nTop 5 quadras minimais:")
        for q in out_quads[:5]:
            print(f"  {q['edge_names']}  freq={q['freq_individual']}")


if __name__ == "__main__":
    main()
