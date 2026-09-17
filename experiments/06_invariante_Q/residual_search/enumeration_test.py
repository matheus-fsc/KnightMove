"""
Enumeração do residual via backtracking com propagação.

Como n_free = 72, 2^72 ≈ 4.7e21 inviabiliza enumeração ingênua.
Estratégia: backtracking com propagação a cada decisão, filtrando
em folha pela conectividade do subgrafo de arestas em 1.

Métricas registradas:
  - nós explorados
  - folhas chegadas
  - 2-fatores válidos (todas regras satisfeitas, sem checar conectividade)
  - tours conexos (é o ground truth: 9862)

Saída: data/residual_enumeration.json
"""

from __future__ import annotations

import json
import time
from collections import deque
from pathlib import Path

import numpy as np

from propagation_engine import (
    FREE,
    PropagationEngine,
    ROOT,
    carregar_dados,
)


class Enumerador:
    def __init__(self, dados, ordem_variavel: str = "freq_dist"):
        self.dados = dados
        self.eng = PropagationEngine(
            edges_uv=dados["edges_uv"],
            freq=dados["freq"],
            pares=dados["pares"],
            triplas=dados["triplas"],
            quadras=dados["quadras"],
            xor=dados["xor"],
            xor_paridade=dados["xor_paridade_alvo"],
        )
        # Propagar nível 5 antes de qualquer branch
        self.eng.run_to_fixpoint([1, 2, 3, 4, 5, 6])
        self.eng.gravar_log = False  # backtracking não precisa de log
        self.eng.log = []

        self.n_arestas = self.eng.n_arestas
        self.ordem_variavel = ordem_variavel

        # Pré-computa ordem global de prioridade
        freq = self.eng.freq
        if ordem_variavel == "freq_dist":
            # menor distância a 0 ou 1 = mais decidido = priorizar
            score = np.minimum(freq, 1 - freq)
            self.prioridade = np.argsort(score).tolist()
        else:
            self.prioridade = list(range(self.n_arestas))

        self.nos = 0
        self.folhas = 0
        self.fatores2 = 0
        self.tours = 0
        self.t0 = 0.0

    # ---------- viabilidade ----------------------------------------------
    def viavel(self) -> bool:
        """Confere R2: nenhum vértice com >2 uns ou com (deg-n_zero) < 2."""
        s = self.eng.estado
        for v, es in self.eng.v2e.items():
            vals = s[es]
            n_um = int((vals == 1).sum())
            n_zero = int((vals == 0).sum())
            if n_um > 2:
                return False
            if self.eng.grau[v] - n_zero < 2:
                return False
        return True

    def escolher_variavel(self) -> int:
        for e in self.prioridade:
            if self.eng.estado[e] == FREE:
                return e
        return -1

    # ---------- conectividade --------------------------------------------
    def conexo_em_um(self) -> bool:
        """Verifica se o subgrafo das arestas em 1 é um único ciclo de 36 vértices."""
        s = self.eng.estado
        ativos = [e for e, v in enumerate(s) if v == 1]
        if len(ativos) != 36:
            return False
        # adjacência
        adj = [[] for _ in range(self.eng.n_vertices)]
        for e in ativos:
            u, v = self.eng.edges_uv[e]
            adj[u].append(v)
            adj[v].append(u)
        # cada vértice tem grau 2 (R2 garante); BFS a partir de 0
        bfs = deque([0])
        seen = {0}
        while bfs:
            x = bfs.popleft()
            for y in adj[x]:
                if y not in seen:
                    seen.add(y)
                    bfs.append(y)
        return len(seen) == 36

    # ---------- branch ---------------------------------------------------
    def branch(self) -> None:
        self.nos += 1

        # Propaga; se cair em contradição, backtrack
        snapshot_estado = self.eng.estado.copy()
        try:
            self.eng.run_to_fixpoint([2, 3, 4, 5, 6])
        except ValueError:
            self.eng.estado = snapshot_estado
            return

        if not self.viavel():
            self.eng.estado = snapshot_estado
            return

        e = self.escolher_variavel()
        if e == -1:
            # todas fixadas — é folha
            self.folhas += 1
            if (self.eng.estado == 1).sum() == 36 and self.viavel():
                self.fatores2 += 1
                if self.conexo_em_um():
                    self.tours += 1
            self.eng.estado = snapshot_estado
            return

        # Tenta valor mais provável primeiro
        f = float(self.eng.freq[e])
        ordem_val = (1, 0) if f >= 0.5 else (0, 1)
        for val in ordem_val:
            inner_estado = self.eng.estado.copy()
            self.eng.estado[e] = val
            self.branch()
            self.eng.estado = inner_estado

        self.eng.estado = snapshot_estado

    def executar(self) -> dict:
        self.t0 = time.perf_counter()
        self.branch()
        dt = time.perf_counter() - self.t0
        return {
            "tempo_s": dt,
            "nos_explorados": self.nos,
            "folhas": self.folhas,
            "n_2fatores": self.fatores2,
            "n_tours_conexos": self.tours,
        }


def main() -> None:
    dados = carregar_dados()
    enumerador = Enumerador(dados, ordem_variavel="freq_dist")
    print(
        f"Iniciando enumeração: n_arestas={enumerador.n_arestas}, "
        f"livres_inicio={(enumerador.eng.estado == FREE).sum()}"
    )
    res = enumerador.executar()
    print("Resultado:")
    for k, v in res.items():
        print(f"  {k}: {v}")
    out = ROOT / "data" / "residual_enumeration.json"
    out.write_text(json.dumps(res, indent=2))
    print(f"\nSalvo em {out}")
    print(
        "\nGround truth: 9.862 tours fechados — "
        f"{'OK' if res['n_tours_conexos'] == 9862 else 'DIVERGE'}"
    )


if __name__ == "__main__":
    main()
