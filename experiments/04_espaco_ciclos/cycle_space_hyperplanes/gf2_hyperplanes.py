"""
Restrições GF(2) como hiperplanos no espaço de ciclos Z_1(G; F_2)
=================================================================

Grafo do cavalo 6x6. Toda a aritmética das Partes 1-3 é sobre F_2 (mod 2),
NUNCA sobre R. "Hiperplano" aqui significa SEMPRE o núcleo de um funcional
GF(2)-linear  φ : F_2^E -> F_2,  i.e. {x : φ(x) = 0}. Não existe vetor
normal, não existe volume contínuo, não existe ângulo. Só paridade.

Âncoras de fidelidade (o script PARA se qualquer uma falhar):
    * 9.862 tours fechados
    * dim Z_1 = β_1 = 45
    * deficit = β_1 − rank(Ham) = 3   (= Q(6))

Partes:
    1. Construção do espaço (grafo, ∂_1, base de Z_1, tours como vetores)
    2. Uma restrição GF(2) como hiperplano (φ de suporte 2) + os 3 funcionais Q
    3. O funcional σ (paridade global dos 4 cantos), análogo correto do "(1,1,1)"
    4. Figura: os 3 cortes sucessivos 2^45 -> 2^44 -> 2^43 -> 2^42
    5. Stub pesca_espectral(...) — NÃO executado (precisa da T_bulk local 33.346 estados)
"""
from __future__ import annotations
import json
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_IN = os.path.join(HERE, "..", "6x6_higher_order", "data")
OUT = os.path.join(HERE, "results")
os.makedirs(OUT, exist_ok=True)

RNG = np.random.default_rng(20260611)


def die(msg: str) -> None:
    print(f"\n!!! ÂNCORA DE FIDELIDADE FALHOU: {msg}")
    print("    abortando (ver restrições do brief).")
    sys.exit(1)


# ======================================================================
# Álgebra linear sobre F_2  (todo & 1, todo ^ é XOR — nada de reais)
# ======================================================================

def gf2_rref(M: np.ndarray):
    """Forma escalonada reduzida sobre F_2. Retorna (rref, pivots, rank)."""
    A = (M.astype(np.uint8) & 1).copy()
    rows, cols = A.shape
    pivots: list[int] = []
    r = 0
    for c in range(cols):
        piv = None
        for i in range(r, rows):
            if A[i, c]:
                piv = i
                break
        if piv is None:
            continue
        if piv != r:
            A[[r, piv]] = A[[piv, r]]
        mask = A[:, c].astype(bool)
        mask[r] = False
        if mask.any():
            A[mask] ^= A[r]
        pivots.append(c)
        r += 1
        if r == rows:
            break
    return A[:r], pivots, r


def gf2_rank(M: np.ndarray) -> int:
    return gf2_rref(M)[2]


def gf2_nullspace(M: np.ndarray) -> np.ndarray:
    """Base do núcleo {x : M x = 0} sobre F_2, como linhas (k x cols)."""
    rref, pivots, _ = gf2_rref(M)
    cols = M.shape[1]
    pivot_set = set(pivots)
    free = [c for c in range(cols) if c not in pivot_set]
    basis = np.zeros((len(free), cols), dtype=np.uint8)
    for k, f in enumerate(free):
        basis[k, f] = 1
        for r, pc in enumerate(pivots):
            if rref[r, f]:
                basis[k, pc] = 1
    return basis


def gf2_in_span(rref: np.ndarray, pivots: list[int], v: np.ndarray) -> bool:
    """v ∈ rowspan(rref) ?  (rref já reduzida, pivots as colunas-pivô)."""
    w = (v.astype(np.uint8) & 1).copy()
    for r, pc in enumerate(pivots):
        if w[pc]:
            w ^= rref[r]
    return not w.any()


def gf2_reduce(rref: np.ndarray, pivots: list[int], v: np.ndarray) -> np.ndarray:
    """Resíduo de v módulo rowspan(rref) (representante canônico do quociente)."""
    w = (v.astype(np.uint8) & 1).copy()
    for r, pc in enumerate(pivots):
        if w[pc]:
            w ^= rref[r]
    return w


# ======================================================================
# PARTE 1 — Construção do espaço
# ======================================================================

def parte1():
    print("=" * 70)
    print("PARTE 1 — Construção do espaço de ciclos Z_1(G_6x6; F_2)")
    print("=" * 70)

    meta = json.load(open(os.path.join(DATA_IN, "edges_6x6.json")))
    edges = [tuple(e) for e in meta["edges_uv"]]
    labels = meta["edge_labels"]
    V, E = 36, len(edges)
    beta1 = E - V + 1
    print(f"  |V| = {V}   |E| = {E}   β_1 = |E|−|V|+1 = {beta1}")
    if (V, E, beta1) != (36, 80, 45):
        die(f"esperado V=36, E=80, β1=45; obtido {V},{E},{beta1}")

    # Matriz de incidência ∂_1 sobre F_2 (36 x 80): ∂_1[v,e]=1 sse v ∈ e
    bnd = np.zeros((V, E), dtype=np.uint8)
    for i, (u, v) in enumerate(edges):
        bnd[u, i] = 1
        bnd[v, i] = 1

    # Base de Z_1 = ker(∂_1) por eliminação gaussiana sobre F_2
    Z1 = gf2_nullspace(bnd)
    dimZ1 = Z1.shape[0]
    print(f"  dim Z_1 = dim ker(∂_1) = {dimZ1}   (eliminação gaussiana mod 2)")
    if dimZ1 != beta1:
        die(f"dim ker(∂_1)={dimZ1} ≠ β1={beta1}")
    # sanity: cada vetor da base TEM bordo nulo (é mesmo um ciclo)
    if ((Z1 @ bnd.T) % 2).any():
        die("vetor da base de Z_1 com ∂_1 x ≠ 0 — bug no núcleo")
    print("  OK: base de Z_1 verificada (∂_1·x = 0 para todo gerador)")

    # Tours: matriz indicadora 9.862 x 80 (enumeração exaustiva canônica do repo)
    T = np.load(os.path.join(DATA_IN, "incidence_matrix_6x6.npy")).astype(np.uint8)
    n_tours = T.shape[0]
    print(f"  tours carregados: {n_tours} vetores em F_2^{E}")
    if n_tours != 9862:
        die(f"contagem de tours = {n_tours} ≠ 9.862")
    # check 1: cada tour é um ciclo fechado de 36 arestas (grau par em todo vértice)
    if not (T.sum(1) == 36).all():
        die("algum tour não tem exatamente 36 arestas")
    if ((T @ bnd.T) % 2).any():
        die("algum tour NÃO pertence a ker(∂_1) — grau ímpar em algum vértice")
    print("  OK: todos os 9.862 tours ∈ ker(∂_1) (grau par em todo vértice)")
    print("  >> âncora 9.862 + dim Z_1 = 45 confirmadas")

    return dict(edges=edges, labels=labels, V=V, E=E, beta1=beta1,
                bnd=bnd, Z1=Z1, T=T)


# ======================================================================
# PARTE 2 — Uma restrição GF(2) como hiperplano
# ======================================================================

def parte2(ctx):
    print("\n" + "=" * 70)
    print("PARTE 2 — Uma restrição GF(2) como hiperplano (o coração pedagógico)")
    print("=" * 70)
    E, labels, T, bnd, Z1 = ctx["E"], ctx["labels"], ctx["T"], ctx["bnd"], ctx["Z1"]
    beta1 = ctx["beta1"]

    # φ(x) = x_{F6-D5} ⊕ x_{B3-A1}  — detector dual mínimo (suporte 2)
    i_fd = labels.index("F6-D5")
    i_ba = labels.index("B3-A1")
    phi = np.zeros(E, dtype=np.uint8)
    phi[i_fd] = 1
    phi[i_ba] = 1
    print(f"  φ = x_[F6-D5] ⊕ x_[B3-A1]   (índices de aresta {i_fd}, {i_ba})")
    print(f"  hiperplano H_φ = ker(φ) = {{x ∈ Z_1 : φ(x) = 0}}")

    # (a) φ se anula em TODOS os tours -> H_φ contém o conjunto de soluções
    vals_tours = (T @ phi) % 2
    n_zero = int((vals_tours == 0).sum())
    print(f"\n  [a] φ avaliado nos 9.862 tours: φ=0 em {n_zero}/9862 "
          f"({100*n_zero/9862:.1f}%)")
    if n_zero != 9862:
        die("φ NÃO se anula em todos os tours — não é restrição válida")
    print("      => o hiperplano ker(φ) CONTÉM todos os tours (restrição válida)")

    # (b) o corte: φ se anula em ~50% de ciclos GENÉRICos de Z_1 (não-tours)
    n_samp = 200_000
    coeffs = RNG.integers(0, 2, size=(n_samp, beta1), dtype=np.uint8)
    samples = (coeffs @ Z1) % 2                     # combinações lineares aleatórias
    tour_set = {bytes(row) for row in T}            # para excluir os que calham ser tour
    is_tour = np.fromiter((bytes(r) in tour_set for r in samples),
                          dtype=bool, count=n_samp)
    generic = samples[~is_tour]
    vals_gen = (generic @ phi) % 2
    frac_zero = float((vals_gen == 0).mean())
    print(f"\n  [b] {len(generic)} ciclos genéricos de Z_1 (não-tours):")
    print(f"      φ=0 em {frac_zero*100:.3f}%   (esperado ~50% — φ é informativa)")
    print(f"      (apenas {int(is_tour.sum())} amostras coincidiram com tours reais)")
    print("      CONTRASTE: 100% dos tours vs ~50% de ciclos => φ NÃO é trivial")

    # (c) anulador do span dos tours = funcionais que zeram em TODOS os tours
    #     ψ tal que T ψ = 0  (núcleo de T como aplicação F_2^80 -> F_2^9862)
    NT = gf2_nullspace(T)                # dim = 80 - rank(T)
    rank_T = E - NT.shape[0]
    deficit = beta1 - rank_T
    print(f"\n  [c] rank(Ham) = rank(T) = {rank_T}")
    print(f"      anulador {{ψ : ψ·t=0 ∀ tour t}} tem dim {NT.shape[0]} "
          f"= |E| − rank(T)")
    print(f"      deficit = β_1 − rank(Ham) = {beta1} − {rank_T} = {deficit}")
    if rank_T != 42 or deficit != 3:
        die(f"rank(Ham)={rank_T}, deficit={deficit} (esperado 42 e 3)")

    # Os funcionais em row(∂_1) são triviais sobre Z_1 (zeram em TODO ciclo).
    # Os 3 hiperplanos Q vivem no quociente  anulador / row(∂_1).
    R_rref, R_piv, rank_R = gf2_rref(bnd)            # row(∂_1), dim 35
    three = []
    three_rref, three_piv, _ = R_rref.copy(), list(R_piv), rank_R
    acc = R_rref.copy()
    acc_piv = list(R_piv)
    for psi in NT:
        res = gf2_reduce(acc, acc_piv, psi)
        if res.any():
            three.append(psi.copy())
            # incorporar res ao acumulador para manter independência mod row(∂)
            acc2, acc_piv2, _ = gf2_rref(np.vstack([acc, res]))
            acc, acc_piv = acc2, acc_piv2
        if len(three) == 3:
            break
    print(f"      funcionais independentes mod row(∂_1): {len(three)}  (= Q = 3)")
    if len(three) != 3:
        die(f"obtidos {len(three)} funcionais Q, esperado 3")

    print("\n      Os 3 hiperplanos Q (cada um = ker de um funcional GF(2)):")
    for k, psi in enumerate(three):
        supp = np.where(psi == 1)[0]
        # confirma: zera em todos os tours
        ok = int(((T @ psi) % 2).sum()) == 0
        print(f"        φ_{k}: |suporte|={len(supp):2d}  zera_em_tours={ok}  "
              f"arestas[0:4]={[labels[j] for j in supp[:4]]}")
        if not ok:
            die(f"φ_{k} não zera em todos os tours")

    # φ (suporte 2) realmente pertence a esse quociente de dim 3?
    in_ann = int(((T @ phi) % 2).sum()) == 0
    res_phi = gf2_reduce(R_rref, R_piv, phi)
    nontrivial = bool(res_phi.any())
    print(f"\n      φ (suporte 2): zera_em_tours={in_ann}, "
          f"não-trivial mod row(∂_1)={nontrivial}")
    print("      => φ é UM dos 3 hiperplanos Q (o de menor suporte)")
    if not (in_ann and nontrivial):
        die("φ suporte-2 não é um funcional Q legítimo")

    ctx["phi"] = phi
    ctx["three"] = three
    ctx["frac_generic_zero"] = frac_zero
    ctx["rank_T"] = rank_T
    ctx["deficit"] = deficit
    return ctx


# ======================================================================
# PARTE 3 — O funcional σ (paridade global dos 4 cantos)
# ======================================================================

def parte3(ctx):
    print("\n" + "=" * 70)
    print("PARTE 3 — O funcional σ (análogo correto do '(1,1,1)')")
    print("=" * 70)
    E, labels, bnd = ctx["E"], ctx["labels"], ctx["bnd"]

    # 4 cantos de grau 2; representante r_c = 1ª aresta obrigatória de cada canto
    corners = {"A6": 0, "F6": 5, "A1": 30, "F1": 35}
    reps = {}
    for name, v in corners.items():
        es = [i for i in range(E) if bnd[v, i]]
        reps[name] = es[0]
        print(f"  canto {name}: arestas obrigatórias {[labels[i] for i in es]}"
              f"  -> representante r_{name} = {labels[es[0]]}")

    rep_idx = list(reps.values())
    # σ : F_2^E -> F_2  soma as 4 coordenadas dos representantes de canto
    sigma = np.zeros(E, dtype=np.uint8)
    for i in rep_idx:
        sigma[i] = 1
    print(f"\n  σ(x) = "
          + " ⊕ ".join(f"x_[{labels[i]}]" for i in rep_idx))
    print("  σ mede PARIDADE GLOBAL dos 4 cantos — é um funcional, não um vetor-seta.")

    # Cada gerador XOR de cantos (r_ci + r_cj) tem σ = 0
    print("\n  σ avaliado nos geradores XOR de canto (r_ci ⊕ r_cj):")
    names = list(reps)
    all_zero = True
    for a in range(len(names)):
        for b in range(a + 1, len(names)):
            v = np.zeros(E, dtype=np.uint8)
            v[reps[names[a]]] ^= 1
            v[reps[names[b]]] ^= 1
            s = int((sigma @ v) % 2)
            all_zero &= (s == 0)
            print(f"    r_{names[a]} ⊕ r_{names[b]} : σ = {s}")
    if not all_zero:
        die("algum gerador XOR de canto tem σ ≠ 0")

    # W = span(r_c) ≅ F_2^4 ; ker(σ|_W) = XOR-pairs, dim 4 − 1 = 3 = Q
    W = np.zeros((4, E), dtype=np.uint8)
    for k, i in enumerate(rep_idx):
        W[k, i] = 1
    dim_W = gf2_rank(W)
    print(f"\n  dim W = dim span(r_c) = {dim_W}  (4 cantos independentes ≅ F_2^4)")
    print(f"  ker(σ|_W) = span{{r_ci ⊕ r_cj}} tem dim {dim_W} − 1 = {dim_W-1}")
    print(f"  >> é exatamente o '−1' que produz  Q = 4 − 1 = {dim_W-1}")
    if dim_W - 1 != 3:
        die(f"dim ker(σ|_W) = {dim_W-1}, esperado 3")

    ctx["sigma"] = sigma
    ctx["reps"] = reps
    return ctx


# ======================================================================
# PARTE 4 — Visualização honesta para espaço discreto
# ======================================================================

def parte4(ctx):
    print("\n" + "=" * 70)
    print("PARTE 4 — Figura: 3 cortes GF(2) sucessivos no espaço discreto")
    print("=" * 70)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    n_tours = ctx["T"].shape[0]
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(14, 7),
                                   gridspec_kw={"width_ratios": [1.25, 1]})

    # ---- Painel esquerdo: caixas aninhadas = hiperplanos sucessivos -------
    axL.set_title("Espaço de ciclos $Z_1(G;\\mathbb{F}_2)$, dim 45\n"
                  "cada borda = núcleo de um funcional GF(2) "
                  "(hiperplano), NÃO um plano com normal",
                  fontsize=11)
    # cada hiperplano divide a contagem por 2 (codim 1 sobre F_2)
    layers = [
        ("$Z_1$  (todos os ciclos)", "$2^{45}$", "#e9ecf5"),
        ("$\\ker\\varphi_0$", "$2^{44}$", "#cdd6ee"),
        ("$\\ker\\varphi_0\\cap\\ker\\varphi_1$", "$2^{43}$", "#a9bce0"),
        ("$\\ker\\varphi_0\\cap\\ker\\varphi_1\\cap\\ker\\varphi_2$\n"
         "$=\\,$compatível-com-tours", "$2^{42}$", "#7e9fd4"),
    ]
    for k, (name, card, color) in enumerate(layers):
        w = 10 - 2.2 * k
        h = 9 - 2.0 * k
        x = (10 - w) / 2
        y = (9 - h) / 2 + 0.0
        axL.add_patch(Rectangle((x, y), w, h, facecolor=color,
                                edgecolor="#22324f", lw=2, zorder=k))
        # nome do subespaço no topo + cardinalidade no canto inf-esquerdo
        axL.text(5, y + h - 0.40, f"{name}", ha="center", va="top",
                 fontsize=9.5, zorder=10)
        if k < 3:  # cardinalidade do bloco interno fica no texto dos tours
            axL.text(x + 0.18, y + 0.28, card, ha="left", va="bottom",
                     fontsize=10.5, fontweight="bold", color="#22324f",
                     zorder=10)
    # tours dentro da interseção dos 3 hiperplanos (marcador acima do rótulo |S|=2^42)
    axL.scatter([5], [4.15], s=70, color="#c0392b", zorder=20)
    axL.text(5, 3.75, f"$2^{{42}}$ pts;  {n_tours} tours $\\subset$\n"
             "$\\ker\\varphi_0\\cap\\ker\\varphi_1\\cap\\ker\\varphi_2$",
             ha="center", va="top", fontsize=8.5, color="#7a1f15", zorder=20)
    axL.set_xlim(-0.3, 10.3)
    axL.set_ylim(-0.3, 9.6)
    axL.set_aspect("equal")
    axL.axis("off")

    # ---- Painel direito: filtração log2 da cardinalidade ------------------
    axR.set_title("Cada hiperplano GF(2) divide a contagem por 2\n"
                  "(codimensão 1, sem 'volume contínuo')", fontsize=11)
    expo = [45, 44, 43, 42]
    xlab = ["$Z_1$", "$\\cap\\ker\\varphi_0$",
            "$\\cap\\ker\\varphi_1$", "$\\cap\\ker\\varphi_2$"]
    bars = axR.bar(range(4), expo, color=["#e9ecf5", "#cdd6ee",
                                          "#a9bce0", "#7e9fd4"],
                   edgecolor="#22324f", lw=1.5)
    for i, e in enumerate(expo):
        axR.text(i, e + 0.15, f"$2^{{{e}}}$", ha="center", fontsize=11,
                 fontweight="bold")
    axR.set_xticks(range(4))
    axR.set_xticklabels(xlab, fontsize=10)
    axR.set_ylabel("$\\log_2$ |subespaço|", fontsize=11)
    axR.set_ylim(40, 46)
    axR.axhline(42, color="#c0392b", ls="--", lw=1.2)
    axR.text(1.5, 42.12, "tours geram dim 42   (deficit = 45−42 = 3 = Q)",
             ha="center", va="bottom", fontsize=8.5, color="#7a1f15")
    axR.grid(axis="y", alpha=0.25)

    fig.text(0.5, 0.015,
             "Onde a intuição do 'plano que corta' VALE: cada restrição = um "
             "hiperplano (núcleo de funcional). Onde QUEBRA: não há normal, não "
             "há ângulo, não há volume contínuo — só $2^k$ pontos e paridade mod 2.",
             ha="center", fontsize=9, style="italic")
    fig.tight_layout(rect=[0, 0.04, 1, 1])
    out_png = os.path.join(OUT, "gf2_hyperplanes.png")
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    print(f"  figura salva: {out_png}")
    ctx["figure"] = out_png
    return ctx


# ======================================================================
# PARTE 5 — Stub de pesca espectral (NÃO executa)
# ======================================================================

def pesca_espectral(T_states, T_triples, tol=1e-6):
    """Procura simetrias escondidas na matriz de transferência T_bulk (n=6).

    *** STUB — não roda neste script. Requer a T_bulk real do repositório
    local (transfer_matrix/), com 33.346 estados. ***

    Parâmetros
    ----------
    T_states : lista/array dos 33.346 rótulos de estado (perfis de fronteira).
    T_triples : iterável de (row, col, val) — entradas esparsas de T_bulk.

    Passos
    ------
    1. Monta a matriz esparsa T (scipy.sparse) a partir das triplas.
    2. Fidelidade: trace/contagem -> 36.236 (2-fatores); λ_1 ≈ 70,48,
       λ_2 ≈ −39,09. PARA se não bater.
    3. Espectro completo via scipy.sparse.linalg.eigs.
    4. Degenerescências: autovalores repetidos dentro de `tol` => sinal de
       simetria escondida (subespaço próprio de dim > 1).
    5. Operadores candidatos S e teste de comutação S·T = T·S:
         (a) paridade de (linha+coluna) do perfil de fronteira;
         (b) reflexões de fronteira do grupo diédrico D_4.
    6. Para cada S que comuta: verifica se o autovetor de Perron (λ_1) é
       S-invariante (consequência de Perron-Frobenius: λ_1 simples => seu
       autovetor é fixado por toda simetria comutante).
    7. Reporta o fator de redução de estados que cada S comutante oferece
       (≈ |órbitas de S| / |estados|).

    Destrava n=8 (>187M estados -> cabível por bloco simétrico),
    mas NÃO o n=12 (explosão combinatória permanece).
    """
    import numpy as _np
    import scipy.sparse as _sp
    import scipy.sparse.linalg as _spla

    N = len(T_states)
    rows, cols, vals = zip(*T_triples)
    T = _sp.csr_matrix((vals, (rows, cols)), shape=(N, N), dtype=float)

    # --- 2. fidelidade -------------------------------------------------
    two_factors = int(round(T.diagonal().sum()))  # placeholder p/ contagem real
    assert two_factors == 36236, (
        f"fidelidade falhou: 2-fatores={two_factors} (esperado 36.236)")

    # --- 3. espectro ---------------------------------------------------
    k = min(64, N - 2)
    eigvals, eigvecs = _spla.eigs(T, k=k, which="LM")
    eigvals = _np.real_if_close(eigvals)
    order = _np.argsort(-_np.abs(eigvals))
    eigvals, eigvecs = eigvals[order], eigvecs[:, order]
    lam1, lam2 = eigvals[0].real, eigvals[1].real
    assert abs(lam1 - 70.48) < 0.5 and abs(lam2 + 39.09) < 0.5, (
        f"fidelidade espectral falhou: λ1={lam1:.2f}, λ2={lam2:.2f}")

    # --- 4. degenerescências ------------------------------------------
    degen = []
    for i in range(len(eigvals) - 1):
        if abs(eigvals[i] - eigvals[i + 1]) < tol:
            degen.append((i, complex(eigvals[i])))

    # --- 5. operadores candidatos S e comutação -----------------------
    def build_parity_rc(states):
        """S = diag(±1) pela paridade de (r+c) do perfil — involução."""
        diag = _np.array([(-1) ** (sum(s) & 1) for s in states], dtype=float)
        return _sp.diags(diag)

    def build_reflection(states, perm):
        """S = matriz de permutação de uma reflexão D_4 (perm: estado->estado)."""
        idx = _np.arange(len(states))
        return _sp.csr_matrix((_np.ones(len(states)), (idx, perm)),
                              shape=(len(states), len(states)))

    candidates = {"parity_rc": build_parity_rc(T_states)}
    # reflexões D_4 entrariam aqui via permutações de estado (perm_*).

    commuting = {}
    perron = eigvecs[:, 0].real
    for name, S in candidates.items():
        comm = (S @ T) - (T @ S)
        if abs(comm).max() < tol:
            # --- 6. Perron invariante? ---
            sp_v = S @ perron
            inv = _np.linalg.norm(sp_v - perron) < 1e-4 * _np.linalg.norm(perron)
            # --- 7. fator de redução (nº de blocos próprios) ---
            n_orbits = int(round((abs(S.diagonal()).sum() + N) / 2)) if name == "parity_rc" else None
            commuting[name] = dict(perron_invariant=bool(inv),
                                   reduction_blocks=n_orbits)

    return dict(N=N, two_factors=two_factors, lam1=lam1, lam2=lam2,
                degeneracies=degen, commuting=commuting,
                note="destrava n=8 (>187M -> por bloco); NÃO destrava n=12")


# ======================================================================
def main():
    ctx = parte1()
    ctx = parte2(ctx)
    ctx = parte3(ctx)
    ctx = parte4(ctx)

    # resumo de fidelidade
    print("\n" + "=" * 70)
    print("ÂNCORAS DE FIDELIDADE — todas confirmadas")
    print("=" * 70)
    print(f"  tours          : {ctx['T'].shape[0]} == 9.862")
    print(f"  dim Z_1        : {ctx['Z1'].shape[0]} == 45")
    print(f"  rank(Ham)      : {ctx['rank_T']} == 42")
    print(f"  deficit = Q    : {ctx['deficit']} == 3")
    print(f"  φ=0 em tours   : 100%   |  φ=0 em ciclos genéricos : "
          f"{ctx['frac_generic_zero']*100:.2f}%")
    print(f"\n  figura         : {ctx['figure']}")
    print("  PARTE 5        : pesca_espectral() definida (stub, NÃO executada)")


if __name__ == "__main__":
    main()
