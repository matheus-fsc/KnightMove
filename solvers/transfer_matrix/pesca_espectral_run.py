"""
Pesca espectral na T_bulk(6) REAL — busca de invariantes de transferência.
==========================================================================

T_bulk(6): transfer matrix coluna-a-coluna do passeio do cavalo 6×6
(transfer_build.py). Estado = (deg_c, deg_{c+1}) ∈ {0,1,2}^6 × {0,1,2}^6.
K = 33.346 estados, nnz = 1.966.428.

ATENÇÃO (PASSO 0): este é um objeto DIFERENTE do knight_transfer.cpp.
Aquele é um broken-profile DP (estado PState{lo,hi}, 4 bits/célula) que
apenas CONTA tours (saída escalar 9862) e nunca materializa uma matriz.
A matriz de 33.346 estados das âncoras é a T_bulk coluna-a-coluna, salva
em data/transfer_n6_bulk.npz. Exportamos os .txt a partir dela.

INVARIANTE DE TRANSFERÊNCIA (alvo aqui) = operador S com S·T = T·S, que
bloco-diagonaliza T e corta o espaço de estados. NÃO confundir com a
obstrução Q=3 (restrição estática sobre F_2 em Z_1, já saturada — outro
objeto, ver cycle_space_hyperplanes/).
"""
from __future__ import annotations
import json
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

THIS = Path(__file__).parent
DATA = THIS / "data"
N = 6

# âncoras (Tabela 17 / eigenvalues_n6.json)
ANCHOR = dict(K=33346, nnz=1966428, two_factors=36236,
              lam1=70.4766, lam2=-39.0867, gap=31.39, ratio=0.5546)


def parse_state(key: str):
    a, b = key.split("|")
    at = tuple(int(x) for x in a.strip("()").split(",") if x.strip())
    bt = tuple(int(x) for x in b.strip("()").split(",") if x.strip())
    return (at, bt)


def load():
    T = sp.load_npz(DATA / f"transfer_n{N}_bulk.npz").astype(np.float64)
    Tnc2 = sp.load_npz(DATA / f"transfer_n{N}_nocp2.npz").astype(np.float64)
    idx_raw = json.load(open(DATA / f"transfer_n{N}_idx.json"))
    idx = {parse_state(k): v for k, v in idx_raw.items()}
    K = T.shape[0]
    state_by_i = [None] * len(idx)
    label_by_i = [None] * len(idx)
    for k, v in idx_raw.items():
        if v < K:                       # idx pode ser estendido (idx_ext)
            label_by_i[v] = k
    for s, v in idx.items():
        if v < K:
            state_by_i[v] = s
    return T, Tnc2, idx, state_by_i, label_by_i


# ======================================================================
# PASSO 1 — Export
# ======================================================================

def passo1_export(T, label_by_i):
    print("=" * 70)
    print("PASSO 1 — Export da matriz (a partir da T_bulk real .npz)")
    print("=" * 70)
    K = T.shape[0]
    states_path = DATA / "Tbulk6_states.txt"
    triples_path = DATA / "Tbulk6_triples.txt"

    with open(states_path, "w") as f:
        for i in range(K):
            f.write(label_by_i[i] + "\n")

    coo = T.tocoo()
    order = np.lexsort((coo.col, coo.row))
    rows, cols, vals = coo.row[order], coo.col[order], coo.data[order].astype(int)
    # escrita eficiente
    buf = np.column_stack([rows, cols, vals])
    np.savetxt(triples_path, buf, fmt="%d")

    n_states = sum(1 for _ in open(states_path))
    n_triples = sum(1 for _ in open(triples_path))
    print(f"  {states_path.name}: {n_states} estados")
    print(f"  {triples_path.name}: {n_triples} triplas (nnz)")
    ok = (n_states == ANCHOR["K"] and n_triples == ANCHOR["nnz"])
    print(f"  esperado: K={ANCHOR['K']}, nnz={ANCHOR['nnz']}  ->  "
          f"{'OK' if ok else 'DIVERGE!'}")
    assert ok, "contagens divergem — matriz exportada não é a do paper"
    return n_states, n_triples


# ======================================================================
# PASSO 2 — Fidelidade (gate obrigatório)
# ======================================================================

def passo2_fidelidade(T, Tnc2, idx):
    print("\n" + "=" * 70)
    print("PASSO 2 — Fidelidade ANTES do espectro (gate obrigatório)")
    print("=" * 70)
    fail = []
    K, nnz = T.shape[0], T.nnz
    print(f"  K   = {K:>8}  (esperado {ANCHOR['K']})")
    print(f"  nnz = {nnz:>8}  (esperado {ANCHOR['nnz']})")
    if K != ANCHOR["K"]:
        fail.append("K")
    if nnz != ANCHOR["nnz"]:
        fail.append("nnz")

    # --- 2-fatores = 36.236 via contração de fronteira (NÃO é o traço!) ---
    s0 = ((0,) * N, (0,) * N)
    s_target = ((2,) * N, (0,) * N)
    v = np.zeros(K, dtype=np.float64)
    v[idx[s0]] = 1.0
    for _ in range(N - 2):           # n-2 passos bulk
        v = T.T @ v
    v = Tnc2.T @ v                   # 1 passo no_cp2 (fecha últimas colunas)
    two_factors = int(round(v[idx[s_target]]))
    print(f"  2-fatores (contração e_s0·T^(n-2)·T_nc2·e_target) = {two_factors}"
          f"  (esperado {ANCHOR['two_factors']})")
    if two_factors != ANCHOR["two_factors"]:
        fail.append("two_factors")

    # --- N(6) = 9862 tours (CONEXOS) ---
    print("  N(6)=9862 tours CONEXOS: NÃO extraível desta matriz.")
    print("    razão: o estado é só o perfil de GRAUS de coluna (deg_c,deg_cp1);")
    print("    ele é cego à conectividade — conta TODOS os 2-fatores (36.236),")
    print("    não distingue o sub-tour conexo. (Consistente com a tese do projeto.)")

    # --- espectro: λ1, λ2, gap, ρ ---
    t0 = time.time()
    vals = spla.eigs(T, k=12, which="LM", return_eigenvectors=False)
    vals = vals[np.argsort(-np.abs(vals))]
    lam1 = vals[0].real
    lam2 = vals[1].real
    gap = abs(lam1) - abs(lam2)
    ratio = abs(lam2) / abs(lam1)
    print(f"  λ1 = {lam1:+.4f}  (esperado {ANCHOR['lam1']:+.4f})")
    print(f"  λ2 = {lam2:+.4f}  (esperado {ANCHOR['lam2']:+.4f})")
    print(f"  gap = {gap:.4f}  (esperado ~{ANCHOR['gap']})")
    print(f"  ρ = |λ2|/|λ1| = {ratio:.4f}  (esperado ~{ANCHOR['ratio']})")
    print(f"  (eigs em {time.time()-t0:.2f}s)")
    if abs(lam1 - ANCHOR["lam1"]) > 0.05:
        fail.append("lam1")
    if abs(lam2 - ANCHOR["lam2"]) > 0.05:
        fail.append("lam2")
    if abs(ratio - ANCHOR["ratio"]) > 0.01:
        fail.append("ratio")

    if fail:
        print(f"\n  !!! FIDELIDADE FALHOU em: {fail} — ABORTANDO PASSO 3")
        raise SystemExit(1)
    print("\n  >> FIDELIDADE CONFIRMADA — espectro confiável, pesca liberada")
    return dict(lam1=lam1, lam2=lam2, gap=gap, ratio=ratio, two_factors=two_factors)


# ======================================================================
# PASSO 3 — A pesca
# ======================================================================

def passo3_pesca(T, idx, state_by_i):
    print("\n" + "=" * 70)
    print("PASSO 3 — A pesca (degenerescências + operadores comutantes)")
    print("=" * 70)
    K = T.shape[0]
    report = {}

    # --- espectro grande para procurar degenerescências ---
    kk = 60
    vals, vecs = spla.eigs(T, k=kk, which="LM")
    order = np.argsort(-np.abs(vals))
    vals, vecs = vals[order], vecs[:, order]

    print(f"\n  Top-{kk} autovalores (|λ| decrescente), procura de degenerescência:")
    # degenerescência = MESMO valor complexo repetido (não par conjugado!)
    tol = 1e-6
    degen = []
    used = [False] * kk
    for i in range(kk):
        if used[i]:
            continue
        group = [i]
        for j in range(i + 1, kk):
            if not used[j] and abs(vals[i] - vals[j]) < tol:
                group.append(j)
                used[j] = True
        if len(group) > 1:
            degen.append((complex(vals[i]), len(group)))
    # pares conjugados (mesmo |λ|, valores distintos) — ESPERADO p/ matriz real
    conj_pairs = 0
    for i in range(kk):
        for j in range(i + 1, kk):
            if abs(vals[i] - np.conj(vals[j])) < 1e-6 and abs(vals[i].imag) > 1e-6:
                conj_pairs += 1
    for r in range(min(10, kk)):
        z = vals[r]
        print(f"    λ_{r:<2} = {z.real:+10.4f} {z.imag:+10.4f}i   "
              f"|λ|={abs(z):8.4f}")
    print(f"  degenerescências (valor complexo repetido, mult>1, tol={tol}): "
          f"{len(degen)}")
    for z, m in degen:
        print(f"    {z:.4f}  multiplicidade {m}")
    print(f"  pares conjugados λ,λ̄ (|λ| igual, valores distintos): {conj_pairs}")
    print("    -> pares conjugados são ARTEFATO de matriz real, NÃO simetria.")
    report["degeneracies"] = [(str(z), m) for z, m in degen]
    report["conjugate_pairs"] = conj_pairs

    # --- (a) operador de paridade diagonal -------------------------------
    print("\n  [a] paridade de (r+c) como operador diagonal S=diag(±1):")
    # candidato natural: σ(s) = (-1)^{Σ deg_c[r]}
    sigma = np.array([(-1.0) ** (sum(state_by_i[i][0]) % 2) for i in range(K)])
    S = sp.diags(sigma)
    comm = (S @ T) - (T @ S)
    nrm = abs(comm).max() if comm.nnz else 0.0
    print(f"      ‖S·T − T·S‖_max = {nrm:.3e}")
    print("      ARGUMENTO GERAL: T_bulk é irredutível (Perron λ1 real simples),")
    print("      logo o grafo de transição é fortemente conexo. Para S=diag(σ),")
    print("      S·T=T·S exige σ(s)=σ(s') em TODA aresta ⇒ σ constante. Portanto")
    print("      NENHUM operador de paridade diagonal não-trivial comuta com T.")
    print("      A 'bipartição (r+c)' não é simetria estática: ela é uma graduação")
    print("      Z2 que depende do índice de coluna — é a própria T que a desloca")
    print("      (por isso λ2<0 sem espectro simétrico λ→−λ).")
    report["parity_diag_commutator_norm"] = float(nrm)
    report["parity_diag_commutes"] = bool(nrm < 1e-9)

    # --- (b)/(c) reflexão de fronteira r -> n-1-r (única D4 que preserva sweep)
    print("\n  [b/c] reflexão vertical R_v: r -> n-1-r (reverte deg-tuples):")
    def reflect(s):
        a, b = s
        return (a[::-1], b[::-1])
    perm = np.empty(K, dtype=np.int64)
    bijection = True
    for i in range(K):
        si = state_by_i[i]
        ri = reflect(si)
        j = idx.get(ri)
        if j is None or j >= K:
            bijection = False
            break
        perm[i] = j
    print(f"      involução bem-definida em todos os {K} estados? {bijection}")
    if bijection:
        assert np.array_equal(perm[perm], np.arange(K)), "perm não é involução"
        # comutação: T[σ(i),σ(j)] == T[i,j]  <=>  T[perm][:,perm] == T
        Tp = T[perm][:, perm]
        diff = (Tp - T)
        nrm_b = abs(diff).max() if diff.nnz else 0.0
        print(f"      ‖P·T − T·P‖_max = {nrm_b:.3e}  "
              f"({'COMUTA' if nrm_b < 1e-9 else 'NÃO comuta'})")
        fix = int((perm == np.arange(K)).sum())     # estados palíndromos
        dim_sym = (K + fix) // 2
        dim_anti = (K - fix) // 2
        print(f"      estados fixos (deg-tuples palíndromos): {fix}")
        print(f"      setor simétrico  dim = (K+fix)/2 = {dim_sym}")
        print(f"      setor antissim.  dim = (K−fix)/2 = {dim_anti}")
        print(f"      maior bloco = {dim_sym};  fator de redução = "
              f"K/{dim_sym} = {K/dim_sym:.3f}")
        report["Rv_commutes"] = bool(nrm_b < 1e-9)
        report["Rv_commutator_norm"] = float(nrm_b)
        report["Rv_fixed_states"] = fix
        report["Rv_block_sym"] = dim_sym
        report["Rv_block_anti"] = dim_anti
        report["Rv_reduction_factor"] = K / dim_sym

        # --- Perron invariante por R_v? ---
        # Perron = autovetor de λ1 (real positivo, simples)
        lam_vals, lam_vecs = spla.eigs(T, k=1, which="LR")
        perron = np.real(lam_vecs[:, 0])
        perron = perron / np.linalg.norm(perron)
        Pperron = perron[perm]
        # sinal ambíguo; alinhar
        if np.dot(Pperron, perron) < 0:
            Pperron = -Pperron
        inv_err = np.linalg.norm(Pperron - perron)
        print(f"      Perron invariante por R_v? ‖P·v1 − v1‖ = {inv_err:.3e}  "
              f"({'SIM' if inv_err < 1e-3 else 'NÃO — ERRO!'})")
        report["perron_Rv_invariant"] = bool(inv_err < 1e-3)
        report["perron_Rv_error"] = float(inv_err)

    # --- D4: quais elementos preservam o sweep de coluna? ---
    print("\n  D₄ vs sweep de coluna (quais comutam com o transfer-por-coluna):")
    print("    identidade ......... preserva colunas .......... comuta (trivial)")
    print("    reflexão vertical R_v (r→n-1-r) .. preserva colunas .. COMUTA  ← (b)")
    print("    reflexão horizontal (c→n-1-c) .... INVERTE o sweep ... mapeia T→Tᵀ")
    print("    rotação 180° (=R_v∘R_h) .......... inverte o sweep ... mapeia T→Tᵀ")
    print("    reflexões diagonais / rot 90° .... trocam linha↔coluna  quebram sweep")
    print("    => só {id, R_v} ⊂ D₄ qualificam: subgrupo Z₂, fator 2 no máximo.")
    report["D4_sweep_preserving"] = ["identity", "R_v (row reflection)"]
    return report


# ======================================================================
def main():
    T, Tnc2, idx, state_by_i, label_by_i = load()
    passo1_export(T, label_by_i)
    fid = passo2_fidelidade(T, Tnc2, idx)
    rep = passo3_pesca(T, idx, state_by_i)
    rep["fidelity"] = fid
    rep["K"] = T.shape[0]
    rep["nnz"] = int(T.nnz)
    json.dump(rep, open(DATA / "pesca_espectral_n6.json", "w"), indent=2)
    print("\n  resultado salvo em data/pesca_espectral_n6.json")


if __name__ == "__main__":
    main()
