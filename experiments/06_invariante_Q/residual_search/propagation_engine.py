"""
Motor de propagação em cascata para o passeio fechado do cavalo 6x6.

Estado por aresta x_e em {0, 1, FREE}. As regras R1..R6 são aplicadas
até ponto fixo. A execução por níveis cumulativos (0..5) registra
quantas variáveis cada nível adicional consegue fixar.

Regras:
  R1  freq[e] == 1.0          → x_e := 1
      freq[e] == 0.0          → x_e := 0
  R2  v com 2 arestas em 1    → as demais incidentes em v := 0
      v com (deg-2) arestas 0 → as 2 restantes := 1
  R3  par (e1,e2) excluído e x_e1=1 → x_e2 := 0
  R4  tripla (e1,e2,e3) excluída e x_e1=x_e2=1 → x_e3 := 0
  R5  quadra (e1..e4) excluída e três delas =1 → a quarta := 0
  R6  cláusula XOR S com paridade alvo: se sobra apenas uma livre,
      seu valor é determinado pela paridade.

O suporte das XOR utilizadas tem paridade alvo 0 (verificado contra
os 9.862 tours fechados em complex_orbit/data/tours_closed_6x6.npy).
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

import numpy as np

FREE = -1
ROOT = Path(__file__).resolve().parent
DATA_6x6 = ROOT.parent / "6x6_higher_order" / "data"
DATA_XOR = ROOT.parent / "forbidden_cycles_6x6" / "data" / "results"
TOURS_PATH = ROOT.parent / "complex_orbit" / "data" / "tours_closed_6x6.npy"


# ---------------------------------------------------------------------------
# Carregamento de dados
# ---------------------------------------------------------------------------

def carregar_dados():
    edges_uv = json.load(open(DATA_6x6 / "edges_6x6.json"))["edges_uv"]
    freq = np.load(DATA_6x6 / "freq_singles.npy")
    inc = np.load(DATA_6x6 / "incidence_matrix_6x6.npy").astype(np.int8)

    pares = json.load(open(DATA_6x6 / "exclusions_pairs.json"))["exclusions"]
    triplas = json.load(open(DATA_6x6 / "exclusions_triples.json"))["minimal_triples"]
    quadras = json.load(open(DATA_6x6 / "exclusions_quads.json"))["minimal_quads"]
    xor = json.load(open(DATA_XOR / "dual_detectors.json"))

    return {
        "edges_uv": edges_uv,
        "freq": freq,
        "incidencia": inc,
        "pares": [tuple(p["edges"]) for p in pares],
        "triplas": [tuple(t["edges"]) for t in triplas],
        "quadras": [tuple(q["edges"]) for q in quadras],
        "xor": [tuple(c["edge_indices"]) for c in xor],
        "xor_paridade_alvo": [0 for _ in xor],
    }


def rotulo_aresta(uv: tuple[int, int]) -> str:
    cols = "ABCDEF"
    u, v = uv
    return f"{cols[u % 6]}{6 - u // 6}-{cols[v % 6]}{6 - v // 6}"


# ---------------------------------------------------------------------------
# Motor
# ---------------------------------------------------------------------------

@dataclass
class PropagationEngine:
    edges_uv: list[tuple[int, int]]
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

        # Estruturas NumPy para propagação vetorizada
        self.pares_arr = (
            np.array(self.pares, dtype=np.int32) if self.pares else np.zeros((0, 2), dtype=np.int32)
        )
        self.triplas_arr = (
            np.array(self.triplas, dtype=np.int32) if self.triplas else np.zeros((0, 3), dtype=np.int32)
        )
        self.quadras_arr = (
            np.array(self.quadras, dtype=np.int32) if self.quadras else np.zeros((0, 4), dtype=np.int32)
        )

        # v2e como matriz padded: linhas = vértices, colunas = índices de aresta (-1 = vazio)
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
        self.contradicao: tuple[int, int, str] | None = None

    # --- operações primitivas ---------------------------------------------
    def fix(self, e: int, valor: int, motivo: str) -> bool:
        atual = int(self.estado[e])
        if atual == valor:
            return False
        if atual != FREE and atual != valor:
            self.contradicao = (e, valor, motivo)
            raise ValueError(
                f"contradição em e={e} ({rotulo_aresta(self.edges_uv[e])}): "
                f"já fixado {atual}, motivo do conflito={motivo}"
            )
        self.estado[e] = valor
        if self.gravar_log:
            self.log.append((e, valor, motivo))
        return True

    # --- regras ------------------------------------------------------------
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
        # Conta uns/zeros por vértice usando a matriz v2e padded
        mat = self.v2e_mat  # (V, max_grau), -1 onde vazio
        vmask = mat >= 0
        # estado nas posições válidas
        idx = np.where(vmask, mat, 0)
        vals = self.estado[idx]  # (V, max_grau)
        vals = np.where(vmask, vals, 99)  # 99 = "fora"
        n_um = (vals == 1).sum(axis=1)
        n_zero = (vals == 0).sum(axis=1)

        if (n_um > 2).any():
            v = int(np.argmax(n_um > 2))
            raise ValueError(f"R2: vértice {v} tem {n_um[v]} arestas em 1")
        if (self.grau_arr - n_zero < 2).any():
            v = int(np.argmax(self.grau_arr - n_zero < 2))
            raise ValueError(
                f"R2: vértice {v} tem só {self.grau_arr[v]-n_zero[v]} não-zero"
            )

        n = 0
        # Caso A: n_um == 2 → demais livres viram 0
        verts_full = np.where(n_um == 2)[0]
        for v in verts_full:
            row = mat[v, : self.v2e_len[v]]
            for e in row:
                if self.estado[e] == FREE:
                    n += int(self.fix(int(e), 0, f"R2:v={v} já tem 2 uns"))
        # Caso B: deg - n_zero == 2 e há livres → todas livres viram 1
        verts_tight = np.where((self.grau_arr - n_zero) == 2)[0]
        for v in verts_tight:
            row = mat[v, : self.v2e_len[v]]
            for e in row:
                if self.estado[e] == FREE:
                    n += int(self.fix(int(e), 1, f"R2:v={v} grau efetivo 2"))
        return n

    def propagar_R3(self) -> int:
        if self.pares_arr.size == 0:
            return 0
        arr = self.pares_arr  # (P, 2)
        vals = self.estado[arr]
        n_um = (vals == 1).sum(axis=1)
        if (n_um == 2).any():
            raise ValueError("R3: par excluído com ambos em 1")
        # linhas com 1 um e 1 livre: a livre vira 0
        ativas = np.where(n_um == 1)[0]
        n = 0
        for r in ativas:
            a, b = int(arr[r, 0]), int(arr[r, 1])
            if self.estado[a] == FREE:
                n += int(self.fix(a, 0, f"R3:par[{r}]"))
            if self.estado[b] == FREE:
                n += int(self.fix(b, 0, f"R3:par[{r}]"))
        return n

    def propagar_R4(self) -> int:
        if self.triplas_arr.size == 0:
            return 0
        arr = self.triplas_arr  # (T, 3)
        vals = self.estado[arr]
        n_um = (vals == 1).sum(axis=1)
        if (n_um == 3).any():
            raise ValueError("R4: tripla excluída com todas em 1")
        ativas = np.where(n_um == 2)[0]
        n = 0
        for r in ativas:
            row = arr[r]
            for e_ in row:
                if self.estado[e_] == FREE:
                    n += int(self.fix(int(e_), 0, f"R4:tripla[{r}]"))
        return n

    def propagar_R5(self) -> int:
        if self.quadras_arr.size == 0:
            return 0
        arr = self.quadras_arr  # (Q, 4)
        vals = self.estado[arr]
        n_um = (vals == 1).sum(axis=1)
        if (n_um == 4).any():
            raise ValueError("R5: quadra excluída com todas em 1")
        ativas = np.where(n_um == 3)[0]
        n = 0
        for r in ativas:
            row = arr[r]
            for e_ in row:
                if self.estado[e_] == FREE:
                    n += int(self.fix(int(e_), 0, f"R5:quadra[{r}]"))
        return n

    def propagar_R6(self) -> int:
        n = 0
        for k, sup in enumerate(self.xor):
            vals = self.estado[list(sup)]
            livres = [e for e, v in zip(sup, vals) if v == FREE]
            soma_fix = int((vals == 1).sum())
            alvo = self.xor_paridade[k]
            if not livres:
                if soma_fix % 2 != alvo:
                    raise ValueError(
                        f"R6: XOR {k} com paridade {soma_fix%2} ≠ alvo {alvo}"
                    )
                continue
            if len(livres) == 1:
                e_ult = livres[0]
                resto = soma_fix % 2
                val = (alvo - resto) % 2
                n += int(self.fix(e_ult, val, f"R6:xor[{k}]"))
        return n

    # --- ciclo até ponto fixo ---------------------------------------------
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
            if 4 in regras:
                mudanca += self.propagar_R4()
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
# Execução em níveis e verificação
# ---------------------------------------------------------------------------

def executar_niveis(dados) -> dict:
    """Roda 6 níveis cumulativos R1..R{n+1} e devolve relatório completo."""
    niveis = {
        0: ([1], "R1: obrigatórias diretas"),
        1: ([1, 2], "R1+R2: obrigatórias + grau-2"),
        2: ([1, 2, 3], "+ R3: pares excluídos"),
        3: ([1, 2, 3, 4], "+ R4: triplas excluídas"),
        4: ([1, 2, 3, 4, 5], "+ R5: quadras excluídas"),
        5: ([1, 2, 3, 4, 5, 6], "+ R6: cláusulas XOR"),
    }

    relatorio = {"niveis": []}
    inc = dados["incidencia"]
    freq = dados["freq"]

    for nivel, (regras, descricao) in niveis.items():
        eng = PropagationEngine(
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
                {
                    "nivel": nivel,
                    "descricao": descricao,
                    "regras": regras,
                    "erro": str(exc),
                }
            )
            return relatorio

        snap = eng.snapshot()
        # Verificação: nenhum fixadas_1 com freq<1 e nenhum fixadas_0 com freq>0
        viola = []
        for e in snap["fixadas_1"]:
            if not np.isclose(freq[e], 1.0):
                viola.append({"aresta": e, "freq": float(freq[e]), "fixada": 1})
        for e in snap["fixadas_0"]:
            if freq[e] > 0:
                # confirmar nos tours reais: nenhum tour deve ter x_e=1
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
                "log_size": len(eng.log),
                "violacoes": viola,
            }
        )
    relatorio["n_arestas"] = len(dados["edges_uv"])
    return relatorio


def main() -> None:
    dados = carregar_dados()
    rel = executar_niveis(dados)
    saida = ROOT / "data" / "propagation_log.json"
    saida.write_text(json.dumps(rel, indent=2))
    print(f"Relatório salvo em {saida}\n")
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
