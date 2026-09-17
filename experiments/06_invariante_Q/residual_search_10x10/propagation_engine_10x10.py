"""
Motor de propagação em cascata para o passeio fechado do cavalo 10×10.

Adaptado de residual_search/propagation_engine.py. Regras ativas:
  R1 obrigatórias diretas (freq=1.0)
  R2 propagação de grau-2 (principal)
  R3 pares estritos: P(A∧B)=0 nas 5000 amostras (84 pares no 10×10)
  R5 quadras estritas (se houver) — não presentes
  R6 cláusulas XOR (3 chosen_pairs com paridade alvo 0)

R4 (triplas) não é adicionada porque o experimento higher-order do
10×10 mostrou que todas as 118 triplas proven são estruturais
(3 arestas no mesmo vértice → redundante com R2).
"""

from __future__ import annotations

import glob
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

import numpy as np

FREE = -1
ROOT = Path(__file__).resolve().parent
DATA_10x10 = ROOT.parent / "board_10x10" / "data"
DATA_XOR = ROOT.parent / "xor_clauses_benchmark" / "data"


def carregar_dados(strict_pair_mode: str = "todos"):
    """strict_pair_mode = 'todos' usa os 84 pares com P(A∧B)=0 nas amostras;
    'candidatos' usa apenas pares dentre os 17 candidatos r<-0.5."""
    edges = json.load(open(DATA_10x10 / "edges_10x10.json"))
    edges_uv = edges["edges_uv"]
    edge_labels = edges["edge_labels"]

    files = sorted(
        glob.glob(str(DATA_10x10 / "samples" / "tours_10x10_batch_*.npy"))
    )
    T = np.concatenate([np.load(f) for f in files], axis=0).astype(np.uint8)
    freq = T.mean(axis=0)
    n_arestas = T.shape[1]
    n_amostras = T.shape[0]

    mand_json = json.load(open(DATA_10x10 / "invariants" / "mandatory_edges.json"))
    mandatory_idx = [m["idx"] for m in mand_json["mandatory"]]

    # Pares estritos
    if strict_pair_mode == "candidatos":
        c = json.load(open(DATA_10x10 / "invariants" / "candidate_clauses.json"))[
            "candidates"
        ]
        pares = []
        for r in c:
            a, b = r["edge_a_idx"], r["edge_b_idx"]
            if (T[:, a] * T[:, b]).mean() == 0.0:
                pares.append((a, b))
    else:
        # todos pares com P=0 nas amostras (entre arestas com freq>0)
        T32 = T.astype(np.int32)
        coex = T32.T @ T32  # (E,E)
        np.fill_diagonal(coex, 1)  # evitar i==j
        iu, ju = np.triu_indices(n_arestas, k=1)
        mask = coex[iu, ju] == 0
        # ambos com freq>0 (caso contrário trivial)
        m_freq = (freq[iu] > 0) & (freq[ju] > 0)
        idx = np.where(mask & m_freq)[0]
        pares = [(int(iu[i]), int(ju[i])) for i in idx]

    # XOR
    xor_meta = json.load(open(DATA_XOR / "xor_10x10_clauses.json"))
    xor = [tuple(p) for p in xor_meta["chosen_pairs"]]

    return {
        "edges_uv": edges_uv,
        "edge_labels": edge_labels,
        "freq": freq,
        "amostras": T,
        "n_amostras": n_amostras,
        "n_arestas": n_arestas,
        "mandatory_idx": mandatory_idx,
        "pares": pares,
        "triplas": [],  # R4 ≡ R2 no 10×10 (118 triplas estruturais)
        "quadras": [],  # não construídas; potencialmente redundantes
        "xor": xor,
        "xor_paridade_alvo": [0 for _ in xor],
    }


def rotulo_aresta_10x10(uv, edge_labels) -> str:
    return edge_labels[uv] if isinstance(uv, int) else f"e({uv[0]},{uv[1]})"


# ---------------------------------------------------------------------------
# Motor (mesmo da versão 6×6, ajustes para 100 vértices)
# ---------------------------------------------------------------------------

@dataclass
class PropagationEngine10:
    edges_uv: list[list[int]]
    freq: np.ndarray
    pares: list[tuple[int, int]] = field(default_factory=list)
    triplas: list[tuple[int, int, int]] = field(default_factory=list)
    quadras: list[tuple[int, int, int, int]] = field(default_factory=list)
    xor: list[tuple[int, ...]] = field(default_factory=list)
    xor_paridade: list[int] = field(default_factory=list)

    def __post_init__(self) -> None:
        self.n_arestas = len(self.edges_uv)
        self.estado = np.full(self.n_arestas, FREE, dtype=np.int8)

        self.v2e: dict[int, list[int]] = {}
        for ei, (u, v) in enumerate(self.edges_uv):
            self.v2e.setdefault(u, []).append(ei)
            self.v2e.setdefault(v, []).append(ei)
        self.n_vertices = max(self.v2e) + 1
        self.grau = {v: len(es) for v, es in self.v2e.items()}

        self.pares_arr = (
            np.array(self.pares, dtype=np.int32) if self.pares else np.zeros((0, 2), dtype=np.int32)
        )
        self.triplas_arr = (
            np.array(self.triplas, dtype=np.int32) if self.triplas else np.zeros((0, 3), dtype=np.int32)
        )
        self.quadras_arr = (
            np.array(self.quadras, dtype=np.int32) if self.quadras else np.zeros((0, 4), dtype=np.int32)
        )

        max_grau = max(self.grau.values())
        self.v2e_mat = np.full((self.n_vertices, max_grau), -1, dtype=np.int32)
        self.v2e_len = np.zeros(self.n_vertices, dtype=np.int32)
        for v, es in self.v2e.items():
            self.v2e_mat[v, : len(es)] = es
            self.v2e_len[v] = len(es)
        self.grau_arr = np.array(
            [self.grau.get(v, 0) for v in range(self.n_vertices)], dtype=np.int32
        )

        self.log: list[tuple[int, int, str]] = []
        self.gravar_log = True

    def fix(self, e: int, valor: int, motivo: str) -> bool:
        atual = int(self.estado[e])
        if atual == valor:
            return False
        if atual != FREE and atual != valor:
            raise ValueError(
                f"contradição em e={e}: já fixado {atual} vs novo {valor}, motivo={motivo}"
            )
        self.estado[e] = valor
        if self.gravar_log:
            self.log.append((e, valor, motivo))
        return True

    def propagar_R1(self) -> int:
        n = 0
        for e in range(self.n_arestas):
            f = float(self.freq[e])
            if f == 1.0:
                n += int(self.fix(e, 1, "R1:freq=1"))
            elif f == 0.0:
                n += int(self.fix(e, 0, "R1:freq=0"))
        return n

    def propagar_R2(self) -> int:
        mat = self.v2e_mat
        vmask = mat >= 0
        idx = np.where(vmask, mat, 0)
        vals = self.estado[idx]
        vals = np.where(vmask, vals, 99)
        n_um = (vals == 1).sum(axis=1)
        n_zero = (vals == 0).sum(axis=1)

        if (n_um > 2).any():
            v = int(np.argmax(n_um > 2))
            raise ValueError(f"R2: vértice {v} com {n_um[v]} uns")
        if (self.grau_arr - n_zero < 2).any():
            v = int(np.argmax(self.grau_arr - n_zero < 2))
            raise ValueError(
                f"R2: vértice {v} com só {self.grau_arr[v]-n_zero[v]} não-zero"
            )

        n = 0
        for v in np.where(n_um == 2)[0]:
            row = mat[v, : self.v2e_len[v]]
            for e in row:
                if self.estado[e] == FREE:
                    n += int(self.fix(int(e), 0, f"R2:v={v} já 2 uns"))
        for v in np.where(self.grau_arr - n_zero == 2)[0]:
            row = mat[v, : self.v2e_len[v]]
            for e in row:
                if self.estado[e] == FREE:
                    n += int(self.fix(int(e), 1, f"R2:v={v} grau efetivo 2"))
        return n

    def propagar_R3(self) -> int:
        if self.pares_arr.size == 0:
            return 0
        vals = self.estado[self.pares_arr]
        n_um = (vals == 1).sum(axis=1)
        if (n_um == 2).any():
            raise ValueError("R3: par estrito com ambos em 1")
        ativas = np.where(n_um == 1)[0]
        n = 0
        for r in ativas:
            a, b = int(self.pares_arr[r, 0]), int(self.pares_arr[r, 1])
            if self.estado[a] == FREE:
                n += int(self.fix(a, 0, f"R3:par[{r}]"))
            if self.estado[b] == FREE:
                n += int(self.fix(b, 0, f"R3:par[{r}]"))
        return n

    def propagar_R5(self) -> int:
        return 0  # sem quadras carregadas

    def propagar_R6(self) -> int:
        n = 0
        for k, sup in enumerate(self.xor):
            vals = self.estado[list(sup)]
            livres = [e for e, v in zip(sup, vals) if v == FREE]
            soma_fix = int((vals == 1).sum())
            alvo = self.xor_paridade[k]
            if not livres:
                if soma_fix % 2 != alvo:
                    raise ValueError(f"R6: XOR {k} viola paridade")
                continue
            if len(livres) == 1:
                e_ult = livres[0]
                val = (alvo - soma_fix % 2) % 2
                n += int(self.fix(e_ult, val, f"R6:xor[{k}]"))
        return n

    def run_to_fixpoint(self, regras: Iterable[int]) -> int:
        regras = set(regras)
        iters = 0
        while True:
            iters += 1
            mudanca = 0
            if 1 in regras:
                mudanca += self.propagar_R1()
            if 2 in regras:
                mudanca += self.propagar_R2()
            if 3 in regras:
                mudanca += self.propagar_R3()
            if 5 in regras:
                mudanca += self.propagar_R5()
            if 6 in regras:
                mudanca += self.propagar_R6()
            if mudanca == 0:
                return iters

    def snapshot(self) -> dict:
        s = self.estado
        return {
            "n_fixadas_1": int((s == 1).sum()),
            "n_fixadas_0": int((s == 0).sum()),
            "n_free": int((s == FREE).sum()),
            "fixadas_1": [int(i) for i in np.where(s == 1)[0]],
            "fixadas_0": [int(i) for i in np.where(s == 0)[0]],
            "free": [int(i) for i in np.where(s == FREE)[0]],
        }


# ---------------------------------------------------------------------------
# Execução por níveis
# ---------------------------------------------------------------------------

def executar_niveis(dados) -> dict:
    niveis = {
        0: ([1], "R1: obrigatórias diretas"),
        1: ([1, 2], "R1+R2: obrigatórias + grau-2"),
        2: ([1, 2, 3], "+ R3: pares estritos"),
        3: ([1, 2, 3, 6], "+ R6: XOR clauses"),
    }
    relatorio = {
        "niveis": [],
        "n_arestas": dados["n_arestas"],
        "n_mandatorias_iniciais": len(dados["mandatory_idx"]),
        "n_pares_estritos": len(dados["pares"]),
        "n_xor": len(dados["xor"]),
        "n_amostras": dados["n_amostras"],
    }
    inc = dados["amostras"]
    freq = dados["freq"]
    for nivel, (regras, descricao) in niveis.items():
        eng = PropagationEngine10(
            edges_uv=dados["edges_uv"],
            freq=freq,
            pares=dados["pares"],
            triplas=dados["triplas"],
            quadras=dados["quadras"],
            xor=dados["xor"],
            xor_paridade=dados["xor_paridade_alvo"],
        )
        try:
            iters = eng.run_to_fixpoint(regras)
        except ValueError as exc:
            relatorio["niveis"].append(
                {"nivel": nivel, "descricao": descricao, "erro": str(exc)}
            )
            return relatorio

        snap = eng.snapshot()
        viola = []
        for e in snap["fixadas_1"]:
            if not np.isclose(freq[e], 1.0):
                viola.append({"aresta": e, "freq": float(freq[e]), "fixada": 1})
        for e in snap["fixadas_0"]:
            if int(inc[:, e].sum()) > 0:
                viola.append({"aresta": e, "freq": float(freq[e]), "fixada": 0})
        relatorio["niveis"].append(
            {
                "nivel": nivel,
                "descricao": descricao,
                "regras": regras,
                "iteracoes": iters,
                "n_fixadas_1": snap["n_fixadas_1"],
                "n_fixadas_0": snap["n_fixadas_0"],
                "n_free": snap["n_free"],
                "fixadas_1": snap["fixadas_1"],
                "fixadas_0": snap["fixadas_0"],
                "free": snap["free"],
                "violacoes": viola,
            }
        )
    return relatorio


def main() -> None:
    dados = carregar_dados(strict_pair_mode="todos")
    print(
        f"Carregado: {dados['n_arestas']} arestas, {dados['n_amostras']} amostras, "
        f"{len(dados['mandatory_idx'])} obrigatórias, "
        f"{len(dados['pares'])} pares estritos, "
        f"{len(dados['xor'])} XORs"
    )
    rel = executar_niveis(dados)
    out = ROOT / "data" / "propagation_log.json"
    out.write_text(json.dumps(rel, indent=2))
    print(f"\nRelatório salvo em {out}\n")
    for n in rel["niveis"]:
        if "erro" in n:
            print(f"Nível {n['nivel']}: ERRO {n['erro']}")
            continue
        print(
            f"Nível {n['nivel']} ({n['descricao']}): "
            f"fixadas_1={n['n_fixadas_1']:3d}  "
            f"fixadas_0={n['n_fixadas_0']:3d}  "
            f"FREE={n['n_free']:3d}  "
            f"iters={n['iteracoes']}  "
            f"violações={len(n['violacoes'])}"
        )


if __name__ == "__main__":
    main()
