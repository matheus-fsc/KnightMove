#!/usr/bin/env python3
"""
coset_tests.py
==============
Validacao da reformulacao por *coset* (COSET_PROMPT.md): os quatro
funcionais de canto y_i : Z_1 -> F_2, a fibra y = 1 onde vivem os tours,
a cota rk(Ham) <= beta_1 - 3 e a classificacao Z_1/Span(Ham).

Aritmetica exata sobre GF(2), vetores como inteiros de Python (bit k =
aresta k). Enumeracao EXAUSTIVA de tours onde ela existe -- nada de
amostragem (ranks amostrais ja' produziram um falso deficit neste
projeto).

CLI:
  python coset_tests.py --scan            # varredura c(Punc) (T2 neg.)
  python coset_tests.py --test T1,T2,T3   # testes rapidos
  python coset_tests.py --test T4,T5,T6,T7 --boards 6x6,5x8
  python coset_tests.py --all             # tudo menos 6x7 (caro)
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
MOVES = [(1, 2), (2, 1), (-1, 2), (-2, 1), (1, -2), (2, -1), (-1, -2), (-2, -1)]


# ---------------------------------------------------------------- grafo

class Board:
    """Grafo do cavalo num tabuleiro R x C, com indexacao de arestas fixa."""

    def __init__(self, R: int, C: int):
        self.R, self.C = R, C
        self.V = R * C
        adj: List[List[int]] = [[] for _ in range(self.V)]
        edges: List[Tuple[int, int]] = []
        for r in range(R):
            for c in range(C):
                i = r * C + c
                for dr, dc in MOVES:
                    rr, cc = r + dr, c + dc
                    if 0 <= rr < R and 0 <= cc < C:
                        j = rr * C + cc
                        adj[i].append(j)
                        if i < j:
                            edges.append((i, j))
        for a in adj:
            a.sort()
        edges.sort()
        self.adj = adj
        self.edges = edges
        self.E = len(edges)
        self.eidx = {e: k for k, e in enumerate(edges)}
        self.beta1 = self.E - self.V + 1

    # -- rotulos (mesma convencao do resto do projeto: coluna 'A'+c, linha R-r)
    def label(self, v: int) -> str:
        r, c = divmod(v, self.C)
        return f"{chr(ord('A') + c)}{self.R - r}"

    def elabel(self, k: int) -> str:
        u, v = self.edges[k]
        return f"{self.label(u)}-{self.label(v)}"

    def corners(self) -> List[int]:
        R, C = self.R, self.C
        return [0, C - 1, (R - 1) * C, (R - 1) * C + C - 1]

    def mandatory(self) -> List[List[int]]:
        """Para cada canto, os indices das suas arestas incidentes."""
        out = []
        for c in self.corners():
            out.append(sorted(self.eidx[(min(c, u), max(c, u))]
                              for u in self.adj[c]))
        return out

    def boundary_rows(self) -> List[int]:
        """Linhas de d_1: uma por vertice, bits = arestas incidentes."""
        rows = [0] * self.V
        for k, (u, v) in enumerate(self.edges):
            rows[u] |= 1 << k
            rows[v] |= 1 << k
        return rows

    def components(self) -> int:
        """Componentes de G. (No 3x3 o centro e' isolado: c(G) = 2.)"""
        seen = [False] * self.V
        comps = 0
        for s in range(self.V):
            if seen[s]:
                continue
            comps += 1
            seen[s] = True
            stack = [s]
            while stack:
                x = stack.pop()
                for u in self.adj[x]:
                    if not seen[u]:
                        seen[u] = True
                        stack.append(u)
        return comps

    def punc_components(self) -> int:
        """Componentes de Punc = G menos as arestas incidentes a cantos."""
        corner = set(self.corners())
        seen = [False] * self.V
        comps = 0
        for s in range(self.V):
            if seen[s]:
                continue
            comps += 1
            seen[s] = True
            stack = [s]
            while stack:
                x = stack.pop()
                if x in corner:
                    continue
                for u in self.adj[x]:
                    if not seen[u] and u not in corner:
                        seen[u] = True
                        stack.append(u)
        return comps


# ------------------------------------------------------------- GF(2)

class Basis:
    """Base incremental em forma escalonada (pivo = bit mais alto)."""

    def __init__(self, vecs: Sequence[int] = ()):
        self.piv: Dict[int, int] = {}
        for v in vecs:
            self.add(v)

    def reduce(self, v: int) -> int:
        while v:
            p = v.bit_length() - 1
            r = self.piv.get(p)
            if r is None:
                return v
            v ^= r
        return 0

    def add(self, v: int) -> bool:
        v = self.reduce(v)
        if v == 0:
            return False
        self.piv[v.bit_length() - 1] = v
        return True

    def contains(self, v: int) -> bool:
        return self.reduce(v) == 0

    @property
    def rank(self) -> int:
        return len(self.piv)

    def vectors(self) -> List[int]:
        return list(self.piv.values())


def gf2_rank(vecs: Sequence[int]) -> int:
    return Basis(vecs).rank


def nullspace(rows: Sequence[int], ncols: int) -> List[int]:
    """Base do nucleo (a direita) da matriz dada pelas linhas `rows`."""
    # RREF com pivo = coluna MAIS BAIXA presente (ordem crescente de coluna)
    mat = [r for r in rows if r]
    pivots: List[int] = []
    pr: List[int] = []
    used = 0
    for c in range(ncols):
        pick = None
        for i in range(used, len(mat)):
            if (mat[i] >> c) & 1:
                pick = i
                break
        if pick is None:
            continue
        mat[used], mat[pick] = mat[pick], mat[used]
        for i in range(len(mat)):
            if i != used and (mat[i] >> c) & 1:
                mat[i] ^= mat[used]
        pivots.append(c)
        pr.append(used)
        used += 1
    pivset = set(pivots)
    free = [c for c in range(ncols) if c not in pivset]
    out: List[int] = []
    for f in free:
        v = 1 << f
        for c, i in zip(pivots, pr):
            if (mat[i] >> f) & 1:
                v |= 1 << c
        out.append(v)
    return out


# --------------------------------------------------- enumeracao de tours

def enumerate_tours(b: Board) -> List[int]:
    """Todos os tours fechados como mascaras de aresta, cada um UMA vez.

    Comeca em 0 (canto, grau 2) e fixa o segundo vertice no menor vizinho,
    o que quebra a simetria ida/volta sem perder nenhum ciclo.
    """
    total, adj, eidx = b.V, b.adj, b.eidx
    used = [False] * total
    used[0] = True
    path = [0]
    out: List[int] = []
    first = min(adj[0])

    def viable() -> bool:
        cur = path[-1]
        free = [v for v in range(total) if not used[v]]
        if not free:
            return True
        for v in free:
            k = 0
            for u in adj[v]:
                if (not used[u]) or u == cur or u == 0:
                    k += 1
                    if k >= 2:
                        break
            if k < 2:
                return False
        seen = {cur}
        stack = [cur]
        cnt = 0
        while stack:
            x = stack.pop()
            for u in adj[x]:
                if not used[u] and u not in seen:
                    seen.add(u)
                    stack.append(u)
                    cnt += 1
        return cnt == len(free)

    def rec():
        cur = path[-1]
        if len(path) == total:
            if 0 in adj[cur]:
                m = 0
                for a, c in zip(path, path[1:] + path[:1]):
                    m ^= 1 << eidx[(min(a, c), max(a, c))]
                out.append(m)
            return
        for u in adj[cur]:
            if used[u]:
                continue
            if len(path) == 1 and u != first:
                continue
            used[u] = True
            path.append(u)
            if viable():
                rec()
            path.pop()
            used[u] = False

    rec()
    return out


# ------------------------------------------------------------- contexto

class Ctx:
    """Tudo que os testes compartilham para um tabuleiro."""

    def __init__(self, R: int, C: int, with_tours: bool = False):
        self.b = b = Board(R, C)
        self.name = f"{R}x{C}"
        self.mand = b.mandatory()
        self.Z1 = nullspace(b.boundary_rows(), b.E)
        self.c_punc = b.punc_components()
        self.c_g = b.components()
        self.tours: List[int] | None = None
        if with_tours:
            t0 = time.time()
            self.tours = enumerate_tours(b)
            self.t_enum = time.time() - t0

    # y_i(z) = coeficiente de z na PRIMEIRA aresta do canto i
    def y(self, z: int) -> int:
        v = 0
        for i, es in enumerate(self.mand):
            v |= ((z >> es[0]) & 1) << i
        return v

    def y_welldefined(self, z: int) -> bool:
        return all(((z >> es[0]) & 1) == ((z >> es[1]) & 1) for es in self.mand)


# --------------------------------------------------------------- testes

def T1(ctx: Ctx) -> dict:
    """y_i bem-definido: sobre uma BASE de Z_1 (basta, por linearidade:
    'coef(e1) = coef(e2)' e' condicao linear)."""
    degs = [len(ctx.b.adj[c]) for c in ctx.b.corners()]
    bad = [i for i, z in enumerate(ctx.Z1) if not ctx.y_welldefined(z)]
    return {
        "board": ctx.name, "dim_Z1": len(ctx.Z1), "beta1": ctx.b.beta1,
        "corner_degrees": degs,
        "mand_edges_distinct": len({e for es in ctx.mand for e in es}) == 8,
        "basis_violations": len(bad),
        "verdict": "PASS" if not bad and degs == [2, 2, 2, 2] else "FAIL",
    }


def T2(ctx: Ctx) -> dict:
    """Posto dos 4 funcionais em Z_1, e dim da fibra y = 0."""
    ys = [ctx.y(z) for z in ctx.Z1]
    # posto do mapa y : Z_1 -> F_2^4  ==  posto da matriz 4 x dim(Z1)
    cols = [0] * 4
    for j, v in enumerate(ys):
        for i in range(4):
            if (v >> i) & 1:
                cols[i] |= 1 << j
    r = gf2_rank(cols)
    dim_ker = len(ctx.Z1) - r
    # identidade estrutural  r = 8 + c(G) - c(Punc); com G conexo, r = 9 - c(Punc)
    pred_r = 8 + ctx.c_g - ctx.c_punc
    dim_zbulk = (ctx.b.E - 8) - ctx.b.V + ctx.c_punc
    return {
        "board": ctx.name, "beta1": ctx.b.beta1, "c_G": ctx.c_g,
        "rank_y": r, "c_Punc": ctx.c_punc, "predicted_rank_9_minus_cPunc": pred_r,
        "dim_fiber_y0": dim_ker, "dim_ZBulk_formula": dim_zbulk,
        "beta1_minus_4": ctx.b.beta1 - 4,
        "bulk_connected": ctx.c_punc == 5,
        "identity_holds": r == pred_r and dim_ker == dim_zbulk,
        "verdict": "PASS" if (r == pred_r and dim_ker == dim_zbulk and
                              (r == 4) == (ctx.c_punc == 5)) else "FAIL",
    }


def Q_paper(b: Board) -> int:
    """Q pela Definicao 2.2 do paper: dim pi(Span{v_ij}), pi mod row(d_1)."""
    rows = b.boundary_rows()
    base = gf2_rank(rows)
    mand = sorted({e for es in b.mandatory() for e in es})
    v = [(1 << mand[0]) | (1 << e) for e in mand[1:]]   # gera todos os v_ij
    return gf2_rank(list(rows) + v) - base


def T3(ctx: Ctx) -> dict:
    """Cota rk(Ham) <= beta_1 - 3, e comparacao com Q do paper."""
    q = Q_paper(ctx.b)
    ys = [ctx.y(z) for z in ctx.Z1]
    r = gf2_rank([sum(((v >> i) & 1) << j for j, v in enumerate(ys))
                  for i in range(4)])
    out = {
        "board": ctx.name, "beta1": ctx.b.beta1, "rank_y": r,
        "Q_paper_def": q, "bound_beta1_minus_3": ctx.b.beta1 - 3,
        "coset_bound": ctx.b.beta1 - (r - 1) if r else None,
    }
    if ctx.tours is not None:
        rk = gf2_rank(ctx.tours)
        out["n_tours"] = len(ctx.tours)
        out["rank_Ham"] = rk
        out["deficit"] = ctx.b.beta1 - rk
        out["bound_respected"] = rk <= ctx.b.beta1 - 3
    out["Q_equals_rank_y_minus_1"] = (q == r - 1) if r else None
    out["verdict"] = "PASS" if (out.get("bound_respected", True) and
                                (q == max(r - 1, 0))) else "FAIL"
    return out


def T4(ctx: Ctx) -> dict:
    """Todo tour esta' na fibra y = 1."""
    assert ctx.tours is not None
    bad = sum(1 for t in ctx.tours if ctx.y(t) != 0b1111)
    illdef = sum(1 for t in ctx.tours if not ctx.y_welldefined(t))
    return {
        "board": ctx.name, "n_tours": len(ctx.tours),
        "enum_seconds": round(ctx.t_enum, 1),
        "tours_not_in_fiber_one": bad, "tours_y_illdefined": illdef,
        "verdict": "PASS" if bad == 0 and illdef == 0 else "FAIL",
    }


def zbulk_basis(ctx: Ctx) -> List[int]:
    """Base do NUCLEO de y em Z_1, i.e. de Z_bulk.

    Nao basta filtrar a base de Z_1 por y(z)=0 -- a fibra e' um subespaco,
    nao um subconjunto da base escolhida. Eliminacao simultanea sobre F_2^4
    (valor de y) e F_2^E (o proprio ciclo).
    """
    piv: Dict[int, Tuple[int, int]] = {}
    ker: List[int] = []
    for z in ctx.Z1:
        yv, zz = ctx.y(z), z
        while yv:
            p = yv.bit_length() - 1
            if p not in piv:
                break
            pv, pz = piv[p]
            yv ^= pv
            zz ^= pz
        if yv == 0:
            ker.append(zz)
        else:
            piv[yv.bit_length() - 1] = (yv, zz)
    return ker


def T5(ctx: Ctx) -> dict:
    """rk(Ham) = beta_1 - 3  <=>  span{tau - tau'} = Z_bulk."""
    assert ctx.tours is not None
    ts = ctx.tours
    rk = gf2_rank(ts)
    t0 = ts[0]
    Bd = Basis(t ^ t0 for t in ts[1:])            # span{tau + tau'}
    Zb = Basis(zbulk_basis(ctx))
    dim_zb = Zb.rank
    contained = all(Zb.contains(v) for v in Bd.vectors())
    equal = contained and Bd.rank == dim_zb
    lhs = (rk == ctx.b.beta1 - 3)
    return {
        "board": ctx.name, "beta1": ctx.b.beta1, "rank_Ham": rk,
        "deficit": ctx.b.beta1 - rk,
        "dim_span_diffs": Bd.rank, "dim_ZBulk": dim_zb,
        "diffs_subset_ZBulk": contained, "diffs_equal_ZBulk": equal,
        "lhs_rank_eq_beta1_minus_3": lhs,
        "equivalence_holds": lhs == equal,
        "verdict": "PASS" if lhs == equal else "FAIL",
    }


def forbidden_patterns(ctx: Ctx) -> dict:
    """Os 3 ciclos proibidos do 6x6, lidos por rotulo de aresta."""
    p = ROOT.parent / "forbidden_cycles_6x6/data/results/forbidden_cycles.json"
    if not p.exists() or ctx.name != "6x6":
        return {}
    lab2k = {ctx.b.elabel(k): k for k in range(ctx.b.E)}
    lab2k.update({"-".join(reversed(s.split("-"))): k for s, k in list(lab2k.items())})
    out = {}
    for fc in json.load(p.open())["forbidden"]:
        m = 0
        for s in fc["edge_labels"]:
            m ^= 1 << lab2k[s]
        assert m in (m,) and gf2_rank([m]) == 1
        # confere que e' mesmo um ciclo (todo grau par)
        deg = {}
        for s in fc["edge_labels"]:
            for v in s.split("-"):
                deg[v] = deg.get(v, 0) + 1
        out[fc["source_edge_label"]] = {
            "mask": m, "is_cycle": all(d % 2 == 0 for d in deg.values()),
            "y": format(ctx.y(m), "04b")[::-1],
        }
    return out


def T6(ctx: Ctx) -> dict:
    """A classificacao: Z_1/Span(Ham) e o mapa induzido por y."""
    assert ctx.tours is not None
    H = Basis(ctx.tours)
    dim_quot = len(ctx.Z1) - H.rank
    # y induzido no quociente: bem definido pois y(Span Ham) c {0,1}
    yspan = {ctx.y(v) for v in H.vectors()}
    ok_welldef = yspan <= {0b0000, 0b1111}
    # imagem de y em Z_1, modulo {0,1}
    img = {ctx.y(z) for z in ctx.Z1}
    imgspan = Basis(list(img))
    # pre-imagem de {0,1}: Z_bulk mais UM representante da fibra y=1 (um tour).
    # O nucleo do mapa induzido e' essa pre-imagem modulo Span(Ham).
    ker_basis = zbulk_basis(ctx)
    extra = Basis(ker_basis + [ctx.tours[0]])
    ker_dim = extra.rank - H.rank
    H_in_extra = all(extra.contains(v) for v in H.vectors())
    res = {
        "board": ctx.name, "dim_Z1": len(ctx.Z1), "rank_Ham": H.rank,
        "dim_quotient": dim_quot,
        "y_of_Ham_span": sorted(format(v, "04b")[::-1] for v in yspan),
        "y_induced_welldefined": ok_welldef,
        "rank_image_y": imgspan.rank,
        "dim_target_F2_4_mod_one": imgspan.rank - 1,
        "dim_preimage_of_{0,1}": extra.rank,
        "Ham_inside_preimage": H_in_extra,
        "dim_kernel_of_induced_y": ker_dim,
        "is_isomorphism": bool(ok_welldef and H_in_extra and ker_dim == 0
                               and dim_quot == imgspan.rank - 1),
    }
    # as 6 (ou k) dimensoes extras tem padrao de canto trivial?
    if ker_dim:
        Zb = Basis(ker_basis)
        # completa Span(Ham) ate' a pre-imagem: os geradores acrescentados sao
        # exatamente as dimensoes que y NAO ve.
        HB = Basis(ctx.tours)
        witnesses = [z for z in ker_basis if HB.add(z)]
        res["extra_dims"] = ker_dim
        res["n_witnesses"] = len(witnesses)
        res["extra_all_in_ZBulk"] = all(Zb.contains(w) for w in witnesses)
        res["extra_y_patterns"] = sorted({format(ctx.y(w), "04b")[::-1]
                                          for w in witnesses})
    fp = forbidden_patterns(ctx)
    if fp:
        res["forbidden_cycles"] = {k: {"y": v["y"], "is_cycle": v["is_cycle"]}
                                   for k, v in fp.items()}
        masks = [v["mask"] for v in fp.values()]
        # independencia modulo {0,1}: os padroes y devem gerar F_2^4/{1}
        pats = Basis([ctx.y(m) for m in masks] + [0b1111])
        res["forbidden_y_rank_with_one"] = pats.rank
        res["forbidden_independent_mod_one"] = pats.rank == 4
        res["forbidden_independent_in_quotient"] = (
            gf2_rank(ctx.tours + masks) == H.rank + 3)
    res["verdict"] = "PASS" if res["is_isomorphism"] or ker_dim else "FAIL"
    return res


def T7(ctx: Ctx) -> dict:
    """Q(R,C) direto + conexidade do bulk."""
    q = Q_paper(ctx.b)
    d = {
        "board": ctx.name, "Q_paper_def": q, "c_Punc": ctx.c_punc,
        "bulk_connected": ctx.c_punc == 5,
    }
    if ctx.tours is not None:
        rk = gf2_rank(ctx.tours)
        d.update(n_tours=len(ctx.tours), rank_Ham=rk,
                 beta1=ctx.b.beta1, deficit=ctx.b.beta1 - rk,
                 tight=(ctx.b.beta1 - rk == q))
    d["verdict"] = "PASS" if q == 3 and ctx.c_punc == 5 else "SEE"
    return d


# ---------------------------------------------------------------- scan

def scan(maxdim: int = 12) -> List[dict]:
    """Varredura de c(Punc): procura o CONTROLE NEGATIVO (bulk desconexo)."""
    out = []
    for R in range(3, maxdim + 1):
        for C in range(R, maxdim + 1):
            b = Board(R, C)
            degs = [len(b.adj[c]) for c in b.corners()]
            if degs != [2, 2, 2, 2]:
                out.append({"board": f"{R}x{C}", "corner_degrees": degs,
                            "c_Punc": None, "note": "canto sem grau 2"})
                continue
            c = b.punc_components()
            out.append({"board": f"{R}x{C}", "V": b.V, "E": b.E,
                        "beta1": b.beta1, "c_Punc": c,
                        "bulk_connected": c == 5,
                        "rank_y_pred": 9 - c, "Q_paper": Q_paper(b)})
    return out


# ----------------------------------------------------------------- main

def parse_board(s: str) -> Tuple[int, int]:
    a, b = s.lower().split("x")
    return int(a), int(b)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--test", default="")
    ap.add_argument("--boards", default="6x6,8x8,6x7,5x8,5x6")
    ap.add_argument("--scan", action="store_true")
    ap.add_argument("--tours", action="store_true",
                    help="enumera tours mesmo se so' T1/T2/T3/T7 forem pedidos")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    results: dict = {}

    if args.scan or args.all:
        s = scan()
        results["scan"] = s
        print("=" * 72)
        print(" VARREDURA c(Punc)  -- controle negativo de T2")
        print("=" * 72)
        print(f" {'board':8s} {'V':>4s} {'E':>4s} {'b1':>4s} {'cPunc':>6s} "
              f"{'rk(y)':>6s} {'Q':>3s}  bulk")
        for row in s:
            if row.get("c_Punc") is None:
                print(f" {row['board']:8s}  (cantos com grau "
                      f"{row['corner_degrees']}) -- fora do escopo")
                continue
            print(f" {row['board']:8s} {row['V']:4d} {row['E']:4d} "
                  f"{row['beta1']:4d} {row['c_Punc']:6d} "
                  f"{row['rank_y_pred']:6d} {row['Q_paper']:3d}  "
                  f"{'conexo' if row['bulk_connected'] else '** DESCONEXO **'}")
        print()

    tests = [t.strip() for t in args.test.split(",") if t.strip()]
    if args.all:
        tests = ["T1", "T2", "T3", "T4", "T5", "T6", "T7"]
    if not tests:
        if args.out:
            Path(args.out).write_text(json.dumps(results, indent=2))
        return

    boards = [parse_board(x) for x in args.boards.split(",")]
    needs_tours = bool({"T4", "T5", "T6"} & set(tests)) or args.tours

    for (R, C) in boards:
        t0 = time.time()
        ctx = Ctx(R, C, with_tours=needs_tours)
        hdr = f" TABULEIRO {R}x{C}   (V={ctx.b.V} E={ctx.b.E} beta1={ctx.b.beta1}"
        if ctx.tours is not None:
            hdr += f" tours={len(ctx.tours)}"
        print("=" * 72)
        print(hdr + ")")
        print("=" * 72)
        for t in tests:
            fn = globals()[t]
            if t in ("T4", "T5", "T6") and not ctx.tours:
                print(f" [{t}] SKIP (sem tours)")
                continue
            r = fn(ctx)
            results.setdefault(t, {})[ctx.name] = r
            print(f" [{t}] {r.get('verdict')}")
            for k, v in r.items():
                if k in ("board", "verdict"):
                    continue
                print(f"      {k:34s}: {v}")
        print(f" ({time.time() - t0:.1f}s)\n")

    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(json.dumps(results, indent=2, default=str))
        print(f"-> {args.out}")


if __name__ == "__main__":
    main()
