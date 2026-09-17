#!/usr/bin/env python3
"""rank_ham_torus.py — Caracterização do rank(Ham(G_T(n)))."""

# --- bootstrap: torna core/ importavel a partir de qualquer subpasta ---
import sys as _sys, pathlib as _pathlib
for _anc in _pathlib.Path(__file__).resolve().parents:
    if (_anc / "core").is_dir() and (_anc / "pyproject.toml").exists():
        _sys.path.insert(0, str(_anc / "core"))
        break
# --- fim do bootstrap ---

import sys
import time
from collections import defaultdict

import numpy as np

sys.path.insert(0, "/home/math/Dev/knight_tour")
from knight_tours_torus import knight_tours, verify_tour, build_graph


# ── infraestrutura ────────────────────────────────────────────────────────────

def get_context(n):
    """build_graph + edge_to_idx + orbit_id (Z_n × Z_n)."""
    ctx = build_graph(n, n)
    V = ctx["V"]
    E = ctx["E"]
    edges = ctx["edge_endpoints"]

    edge_list = [(int(edges[ei, 0]), int(edges[ei, 1])) for ei in range(E)]
    edge_to_idx = {e: ei for ei, e in enumerate(edge_list)}

    # Órbitas sob translação Z_n × Z_n: indexadas por (Δr, Δc) mod n
    # com canonicalização para (Δr, Δc) ~ (-Δr, -Δc) (mesma aresta)
    orbit_id = np.full(E, -1, dtype=np.int32)
    canon_to_oid = {}
    for ei, (u, v) in enumerate(edge_list):
        ru, cu = divmod(u, n)
        rv, cv = divmod(v, n)
        dr1, dc1 = (rv - ru) % n, (cv - cu) % n
        dr2, dc2 = (n - dr1) % n, (n - dc1) % n
        canon = min((dr1, dc1), (dr2, dc2))
        if canon not in canon_to_oid:
            canon_to_oid[canon] = len(canon_to_oid)
        orbit_id[ei] = canon_to_oid[canon]

    return {
        "n": n, "V": V, "E": E, "beta1": E - V + 1,
        "edges": edge_list, "edge_to_idx": edge_to_idx,
        "orbit_id": orbit_id, "n_orbits": len(canon_to_oid),
        "orbit_canon": {v: k for k, v in canon_to_oid.items()},
        "ctx": ctx,
    }


def tour_to_signature(tour, edge_to_idx, E, V):
    sig = np.zeros(E, dtype=np.uint8)
    for i in range(V):
        u = int(tour[i])
        v = int(tour[(i + 1) % V])
        sig[edge_to_idx[(min(u, v), max(u, v))]] = 1
    return sig


class IncrementalRank:
    """Mantém RREF parcial em GF(2): para cada base[i] = (pivot_col, vec),
    pivot_col é a primeira coordenada não-zero, e vetor é reduzido pelos
    pivôs anteriores. Não exige sort por pivot."""

    def __init__(self, E):
        self.E = E
        self.pivot_cols = []
        self.basis = []

    def add(self, v):
        v = v.copy()
        for pcol, bvec in zip(self.pivot_cols, self.basis):
            if v[pcol]:
                v ^= bvec
        nz = np.flatnonzero(v)
        if len(nz) > 0:
            self.pivot_cols.append(int(nz[0]))
            self.basis.append(v)
            return True
        return False

    @property
    def rank(self):
        return len(self.basis)


def collect_tours(n, K_target, batch_per_seed=None, seeds_max=2000, log=True):
    """Coleta K_target tours diversos. Usa MUITAS seeds com batches pequenos."""
    if batch_per_seed is None:
        # batch grande o suficiente para amortizar overhead mas não dominar
        batch_per_seed = max(50, K_target // 100)
    tours_all = []
    seen = set()
    t0 = time.time()
    for seed in range(seeds_max):
        remaining = K_target - len(tours_all)
        if remaining <= 0:
            break
        batch_size = min(batch_per_seed, remaining)
        batch = knight_tours(n, n, batch_size, seed=seed)
        new = 0
        for t in batch:
            key = tuple(int(x) for x in t)
            if key not in seen:
                seen.add(key)
                tours_all.append(t)
                new += 1
        if log and (seed % 20 == 0 or len(tours_all) >= K_target):
            elapsed = time.time() - t0
            rate = len(tours_all) / max(elapsed, 1e-6)
            print(f"    seed={seed:>4}  total={len(tours_all):>6}  "
                  f"+{new:>4} novos  t={elapsed:.1f}s  ({rate:.0f}/s)")
    return tours_all[:K_target]


# ── Experimento A ─────────────────────────────────────────────────────────────

def experiment_A(n=6, K=10000, log_every=100):
    print(f"\n{'═'*72}")
    print(f"EXPERIMENTO A — rank real (n={n}, K={K})")
    print(f"{'═'*72}")

    info = get_context(n)
    V, E = info["V"], info["E"]
    print(f"V={V}  E={E}  β₁={info['beta1']}  #órbitas Z_{n}×Z_{n}={info['n_orbits']}")
    orbit_sizes = [int((info["orbit_id"] == o).sum())
                   for o in range(info["n_orbits"])]
    print(f"Tamanhos das órbitas: {orbit_sizes}")
    print(f"Canônicos: {info['orbit_canon']}")

    print(f"\nColetando até {K} tours diversos...")
    t0 = time.time()
    tours = collect_tours(n, K, log=True)
    t_collect = time.time() - t0
    print(f"  {len(tours)} tours em {t_collect:.1f}s")

    print(f"\nVerificando tours (sample)...")
    sample_idx = list(range(0, len(tours), max(1, len(tours) // 20)))[:20]
    for i in sample_idx:
        assert verify_tour(tours[i], n, n), f"tour {i} inválido"
    print(f"  {len(sample_idx)} samples OK")

    print(f"\nRank incremental...")
    t0 = time.time()
    rank_eng = IncrementalRank(E)
    rank_log = []
    for k, tour in enumerate(tours):
        sig = tour_to_signature(tour, info["edge_to_idx"], E, V)
        rank_eng.add(sig)
        if (k + 1) % log_every == 0:
            rank_log.append((k + 1, rank_eng.rank))
    t_rank = time.time() - t0

    rank_final = rank_eng.rank
    print(f"  rank final: {rank_final} / β₁={info['beta1']}")
    print(f"  deficit:    {info['beta1'] - rank_final}")
    print(f"  rank/β₁:    {rank_final/info['beta1']:.3f}")
    print(f"  tempo:      {t_rank:.1f}s")

    print(f"\nEvolução do rank:")
    n_logs = len(rank_log)
    indices = sorted(set(list(range(0, n_logs, max(1, n_logs // 25))) + [n_logs - 1]))
    for i in indices:
        k, r = rank_log[i]
        print(f"  k={k:>6}: rank={r}")

    # Convergência: últimas 10 amostras devem ser iguais
    last_vals = [r for _, r in rank_log[-10:]]
    converged = len(set(last_vals)) == 1
    print(f"\nConvergiu nas últimas 10 amostras? {converged}")
    print(f"  últimos valores: {last_vals}")

    return info, tours, rank_eng, rank_log


# ── Experimento B ─────────────────────────────────────────────────────────────

def translate_edge(edge_uv, n, da, db):
    """Translação (da, db): vértice w = r*n + c → (r+da, c+db) mod n."""
    u, v = edge_uv
    ru, cu = divmod(u, n)
    rv, cv = divmod(v, n)
    u2 = ((ru + da) % n) * n + ((cu + db) % n)
    v2 = ((rv + da) % n) * n + ((cv + db) % n)
    return (min(u2, v2), max(u2, v2))


def build_translation_perm(n, da, db, edge_to_idx):
    """Permutação σ: e → σ(e) como np.array de ints (índice de aresta)."""
    E = len(edge_to_idx)
    perm = np.zeros(E, dtype=np.int32)
    for e, ei in edge_to_idx.items():
        perm[ei] = edge_to_idx[translate_edge(e, n, da, db)]
    return perm


def experiment_B(info, rank_eng):
    print(f"\n{'═'*72}")
    print(f"EXPERIMENTO B — estrutura dos vetores da base (n={info['n']})")
    print(f"{'═'*72}")

    n = info["n"]
    V, E = info["V"], info["E"]
    edge_to_idx = info["edge_to_idx"]

    basis = rank_eng.basis
    print(f"Base com {len(basis)} vetores.")

    # ── B.1: Invariância por translação ────────────────────────────────────
    print(f"\n[B.1] Invariância de cada vetor da base por translação:")
    print(f"  Para cada vetor b, conta quantas translações σ ∈ Z_n×Z_n fixam b")
    print(f"  (estabilizador). Reporta histograma do tamanho do estabilizador.")

    # Compute permutations once
    perms = {}
    for da in range(n):
        for db in range(n):
            perms[(da, db)] = build_translation_perm(n, da, db, edge_to_idx)

    stab_sizes = []
    for b in basis:
        stab = 0
        for (da, db), perm in perms.items():
            if np.array_equal(b, b[perm]):
                stab += 1
        stab_sizes.append(stab)

    hist = defaultdict(int)
    for s in stab_sizes:
        hist[s] += 1
    print(f"  Histograma |estab| → contagem:")
    for s in sorted(hist):
        print(f"    |estab|={s:>3} ({n*n//s:>3} órbitas-de-orbita): {hist[s]} vetores")

    # ── B.2: Uniformidade nas órbitas de arestas ────────────────────────────
    print(f"\n[B.2] Cada vetor da base é uniforme nas {info['n_orbits']} órbitas?")
    orbit_id = info["orbit_id"]
    n_orbits = info["n_orbits"]

    uniform_count = 0
    distrib_summary = []
    for bi, b in enumerate(basis):
        # Conta 1's em cada órbita
        per_orbit = []
        for o in range(n_orbits):
            mask = (orbit_id == o)
            n_ones = int(b[mask].sum())
            n_total = int(mask.sum())
            per_orbit.append((n_ones, n_total))
        # "Uniforme" = todas as órbitas têm 0 ou todas têm o seu tamanho
        is_uniform = all((n_ones == 0 or n_ones == n_total) for n_ones, n_total in per_orbit)
        if is_uniform:
            uniform_count += 1
        distrib_summary.append((bi, per_orbit, is_uniform))

    print(f"  Vetores uniformes nas órbitas: {uniform_count} / {len(basis)}")
    if uniform_count > 0:
        print(f"  (primeiros uniformes:)")
        n_shown = 0
        for bi, per_orbit, is_u in distrib_summary:
            if is_u and n_shown < 5:
                pattern = "/".join(f"{n_ones}/{n_total}" for n_ones, n_total in per_orbit)
                print(f"    b[{bi}]: {pattern}")
                n_shown += 1

    # ── B.3: distribuição de 1's por órbita para vetores típicos ────────────
    print(f"\n[B.3] Distribuição de 1's por órbita para os primeiros 5 vetores:")
    for bi, per_orbit, _ in distrib_summary[:5]:
        s = " | ".join(f"o{o}: {n_ones:>3}/{n_total}"
                      for o, (n_ones, n_total) in enumerate(per_orbit))
        print(f"    b[{bi}]: {s}")

    # ── B.4: dimensão do subespaço Z_n × Z_n - invariante do span ──────────
    print(f"\n[B.4] dim(invariante)(Span Ham) sob Z_n × Z_n:")
    print(f"  Calcula a projeção orbital: x_orbital ∈ GF(2)^{n_orbits}")
    print(f"  com (x_orbital)_o = Σ_{{e ∈ órbita o}} x[e]  (mod 2)")
    # Construir matriz orbital
    M_orbital = np.zeros((len(basis), n_orbits), dtype=np.uint8)
    for bi, b in enumerate(basis):
        for o in range(n_orbits):
            mask = (orbit_id == o)
            M_orbital[bi, o] = int(b[mask].sum()) % 2

    # Rank GF(2) da matriz orbital
    def gf2_rank_matrix(M):
        A = M.copy()
        rows, cols = A.shape
        r = 0
        for c in range(cols):
            if r >= rows:
                break
            piv = next((i for i in range(r, rows) if A[i, c]), None)
            if piv is None:
                continue
            A[[r, piv]] = A[[piv, r]]
            for i in range(rows):
                if i != r and A[i, c]:
                    A[i] ^= A[r]
            r += 1
        return r

    rank_orbital = gf2_rank_matrix(M_orbital)
    print(f"  dim Span(projeções orbitais) = {rank_orbital} (≤ {n_orbits})")

    # ── B.5: paridade Σ x[e] = n² mod 2 ────────────────────────────────────
    print(f"\n[B.5] paridade Σx[e] mod 2 dos vetores da base:")
    parities = [int(b.sum()) % 2 for b in basis]
    n0 = parities.count(0)
    n1 = parities.count(1)
    print(f"  Σ=par: {n0},  Σ=ímpar: {n1}")
    print(f"  (Cada tour tem Σ = V = {V}, mod 2 = {V % 2})")


# ── Experimento C ─────────────────────────────────────────────────────────────

def experiment_C(ns=(4, 6, 8), K_per_n=None):
    print(f"\n{'═'*72}")
    print(f"EXPERIMENTO C — scaling do rank em n ∈ {list(ns)}")
    print(f"{'═'*72}")

    if K_per_n is None:
        K_per_n = {4: 5000, 6: 10000, 8: 3000}

    results = []
    for n in ns:
        K = K_per_n.get(n, 5000)
        print(f"\n── n={n} (K={K}) ──────────────────────────")
        info = get_context(n)
        V, E = info["V"], info["E"]
        print(f"  V={V} E={E} β₁={info['beta1']} #órbitas={info['n_orbits']}")
        orbit_sizes = [int((info["orbit_id"] == o).sum())
                       for o in range(info["n_orbits"])]
        print(f"  Órbitas: {orbit_sizes}")

        t0 = time.time()
        tours = collect_tours(n, K, log=False)
        t_collect = time.time() - t0
        print(f"  {len(tours)} tours em {t_collect:.1f}s")

        t0 = time.time()
        rank_eng = IncrementalRank(E)
        rank_log = []
        for k, tour in enumerate(tours):
            sig = tour_to_signature(tour, info["edge_to_idx"], E, V)
            rank_eng.add(sig)
            log_every = max(50, K // 50)
            if (k + 1) % log_every == 0:
                rank_log.append((k + 1, rank_eng.rank))
        t_rank = time.time() - t0

        # Convergence
        last_vals = [r for _, r in rank_log[-5:]] if len(rank_log) >= 5 else []
        converged = len(set(last_vals)) == 1 if last_vals else False

        results.append({
            "n": n, "V": V, "E": E, "beta1": info["beta1"],
            "n_orbits": info["n_orbits"], "orbit_sizes": orbit_sizes,
            "K": len(tours), "rank": rank_eng.rank,
            "deficit": info["beta1"] - rank_eng.rank,
            "rank_over_beta1": rank_eng.rank / max(info["beta1"], 1),
            "converged": converged, "last_vals": last_vals,
            "t_collect": t_collect, "t_rank": t_rank,
        })
        print(f"  rank={rank_eng.rank} / β₁={info['beta1']}  "
              f"deficit={info['beta1']-rank_eng.rank}  "
              f"conv={converged} (últimos: {last_vals})")

    # ── Tabela final ────────────────────────────────────────────────────────
    print(f"\n{'─'*78}")
    print(f"{'n':>3} | {'β₁':>5} | {'rank':>5} | {'def':>5} | "
          f"{'r/β₁':>6} | {'orb':>4} | {'K':>6} | {'conv?':>6}")
    print("─" * 78)
    for r in results:
        print(f"{r['n']:>3} | {r['beta1']:>5} | {r['rank']:>5} | "
              f"{r['deficit']:>5} | {r['rank_over_beta1']:>6.3f} | "
              f"{r['n_orbits']:>4} | {r['K']:>6} | {str(r['converged']):>6}")
    print("─" * 78)

    return results


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    print("rank_ham_torus.py — 2026-05-22")
    print("=" * 72)

    t0 = time.time()

    # Exp A: n=6, K=10000 (ajustar se demorar muito)
    info, tours, rank_eng, rank_log = experiment_A(n=6, K=10000, log_every=200)

    # Exp B usa o resultado de A
    experiment_B(info, rank_eng)

    # Exp C: n=4, 6, 8
    experiment_C(ns=(4, 6, 8))

    print(f"\nTempo total: {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
