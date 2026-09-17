"""
Convergência de freq[e] em vértices interiores para n ∈ {6, 8, 10, 12, 14}.

Hipótese: para arestas cujos dois vértices têm interior_level ≥ d, freq[e]
converge a um valor universal quando n cresce (i.e., a "fase local" do
cavalo existe empiricamente — motiva prova via G_∞ no plano infinito).

interior_level(v) = min(r, c, n-1-r, n-1-c)
  - 0 = borda
  - 1, 2, … = anéis sucessivos para dentro
  - n=6 → max 2;  n=8 → max 3;  n=10 → max 4;  n=12 → max 5;  n=14 → max 6.

interior_level(e) = min(level(u), level(v)).

Amostragem: SamplingV2 = MinimalV2 com tiebreaks aleatorizados:
  - entre vértices de pressão máxima, escolha uniforme
  - entre arestas livres do vértice escolhido, escolha uniforme
  - ordem de val (0,1) shuffled
Resultado: K=500 tours diversos por n, em segundos.

Saídas: data/edge_freq_scaling.json + data/plots/edge_freq_convergence.png
"""

from __future__ import annotations

import json
import sys
import time
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from scaling_minimal_v2 import MinimalV2, FREE  # noqa: E402
from subtour_detector import OK  # noqa: E402


class SamplingV2(MinimalV2):
    def __init__(self, n: int, alvo: int, timeout: float, seed: int = 0):
        super().__init__(n=n, alvo=alvo, timeout=timeout)
        self.rng = np.random.default_rng(seed)
        self.tour_records: list[np.ndarray] = []

    def escolher_var(self) -> int:
        s = self.estado
        mat = self.v2e_mat
        vmask = mat >= 0
        idx = np.where(vmask, mat, 0)
        vals = s[idx]
        vals = np.where(vmask, vals, 99)
        n_um = (vals == 1).sum(axis=1)
        n_zero = (vals == 0).sum(axis=1)
        n_free_v = (vals == -1).sum(axis=1)
        pressao = n_um.astype(np.int64) * 100 + (self.grau_arr - n_zero - 2)
        pressao = np.where(n_free_v > 0, pressao, -10**9)
        max_p = pressao.max()
        if max_p <= -10**8:
            return -1
        cands_v = np.where(pressao == max_p)[0]
        v_alvo = int(self.rng.choice(cands_v))
        row = mat[v_alvo, : self.v2e_len[v_alvo]]
        livres = [int(e) for e in row if self.estado[e] == FREE]
        if not livres:
            return -1
        return int(self.rng.choice(livres))

    def branch(self):
        if self.parar:
            return
        self.nos += 1
        if (self.nos & 0xFFFF) == 0:
            if time.perf_counter() - self.t0 > self.timeout:
                self.parar = True
                return

        snap_est = self.estado.copy()
        snap_seen = self.last_seen.copy()
        cp = self.detector.checkpoint()

        try:
            self.propagar()
        except ValueError:
            self.podas["R2"] += 1
            self.estado = snap_est
            self.last_seen = snap_seen
            self.detector.rollback(cp)
            return

        if self.sync_det() != OK:
            self.estado = snap_est
            self.last_seen = snap_seen
            self.detector.rollback(cp)
            return

        e = self.escolher_var()
        if e == -1:
            self.folhas += 1
            if (self.estado == 1).sum() == self.V:
                self.fatores2 += 1
                if self.conexo():
                    self.tours += 1
                    self.tour_records.append((self.estado == 1).astype(np.uint8))
                    if self.tours >= self.alvo:
                        self.parar = True
            self.estado = snap_est
            self.last_seen = snap_seen
            self.detector.rollback(cp)
            return

        ordem = [1, 0] if self.rng.random() < 0.5 else [0, 1]
        for val in ordem:
            if self.parar:
                break
            inner_est = self.estado.copy()
            inner_seen = self.last_seen.copy()
            inner_cp = self.detector.checkpoint()
            self.estado[e] = val
            self.branch()
            self.estado = inner_est
            self.last_seen = inner_seen
            self.detector.rollback(inner_cp)

        self.estado = snap_est
        self.last_seen = snap_seen
        self.detector.rollback(cp)


def interior_level_vertice(i: int, n: int) -> int:
    r, c = i // n, i % n
    return min(r, c, n - 1 - r, n - 1 - c)


def interior_level_aresta(edges_uv, n: int, e: int) -> int:
    u, v = edges_uv[e]
    return min(interior_level_vertice(u, n), interior_level_vertice(v, n))


def coletar(n: int, K: int, seed: int = 42, timeout: float = 300.0,
            n_restarts: int = 1):
    """Coleta K_total = K * n_restarts tours, com n_restarts seeds
    independentes para diversificar. Cada run interno acumula K tours."""
    print(f"\n=== n={n}  K_por_run={K}  n_restarts={n_restarts} ===", flush=True)
    todos: list[np.ndarray] = []
    nodes_total = 0
    tempo_total = 0.0
    base = SamplingV2(n=n, alvo=1, timeout=timeout, seed=seed)
    E_local = base.E
    edges_uv_local = base.edges_uv
    V_local = base.V
    del base
    for k in range(n_restarts):
        s = SamplingV2(n=n, alvo=K, timeout=timeout, seed=seed + 1000 * k + 7)
        r = s.executar()
        if "erro" in r:
            raise RuntimeError(f"erro em n={n}: {r['erro']}")
        todos.extend(s.tour_records)
        nodes_total += s.nos
        tempo_total += r["tempo_s"]
        if (k + 1) % 5 == 0 or k == n_restarts - 1:
            print(f"  restart {k+1:2d}/{n_restarts}: acum {len(todos)} tours, "
                  f"{tempo_total:.2f}s", flush=True)
    T = np.stack(todos, axis=0)
    freq = T.mean(axis=0)
    s = type("S", (), {"E": E_local, "edges_uv": edges_uv_local,
                       "V": V_local, "tours": len(todos)})()
    r = {"tempo_s": tempo_total, "nos_por_tour": nodes_total / max(1, len(todos))}
    levels = np.array(
        [interior_level_aresta(s.edges_uv, n, e) for e in range(s.E)],
        dtype=np.int32,
    )

    # agregação por nível
    por_nivel = {}
    for L in sorted(set(levels.tolist())):
        mask = levels == L
        por_nivel[int(L)] = {
            "n_arestas": int(mask.sum()),
            "freq_mean": float(freq[mask].mean()),
            "freq_std": float(freq[mask].std()),
            "freq_min": float(freq[mask].min()),
            "freq_max": float(freq[mask].max()),
        }

    return {
        "n": n,
        "V": s.V,
        "E": s.E,
        "K": int(s.tours),
        "tempo_s": r["tempo_s"],
        "nodes_per_tour": r["nos_por_tour"],
        "por_nivel": por_nivel,
        "freq_all": freq.tolist(),
        "levels_all": levels.tolist(),
    }


def main():
    seeds = {6: 11, 8: 13, 10: 17, 12: 19, 14: 23}
    # K_por_run × n_restarts → tours diversos por n
    plano = {
        6:  (50, 40),
        8:  (50, 40),
        10: (50, 40),
        12: (50, 30),
        14: (25, 20),  # n=14 mais caro; menor amostra
    }
    resultados = []
    for n in [6, 8, 10, 12, 14]:
        K, n_restarts = plano[n]
        r = coletar(n, K, seed=seeds[n], n_restarts=n_restarts)
        print(f"  por nível:")
        for L, st in r["por_nivel"].items():
            print(
                f"    L={L}: n_arestas={st['n_arestas']:>4} "
                f"freq_mean={st['freq_mean']:.4f} ± {st['freq_std']:.4f}  "
                f"(min={st['freq_min']:.3f}, max={st['freq_max']:.3f})"
            )
        resultados.append(r)

    # salvar
    out_dir = ROOT / "data"
    out_dir.mkdir(exist_ok=True)
    saida = {"K": K, "seeds": seeds, "resultados": resultados}
    (out_dir / "edge_freq_scaling.json").write_text(json.dumps(saida, indent=2))
    print(f"\nSalvo em {out_dir / 'edge_freq_scaling.json'}")

    # plot
    plot(resultados, out_dir / "plots" / "edge_freq_convergence.png")


def plot(resultados, out_path: Path):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    ns = [r["n"] for r in resultados]
    max_L = max(L for r in resultados for L in r["por_nivel"])

    fig, axes = plt.subplots(1, 2, figsize=(13, 7))

    # Painel 1: freq média por nível (linhas = nível, x = n)
    ax = axes[0]
    cmap = plt.cm.viridis
    for L in range(max_L + 1):
        xs, ys, errs = [], [], []
        for r in resultados:
            if L in r["por_nivel"]:
                xs.append(r["n"])
                ys.append(r["por_nivel"][L]["freq_mean"])
                errs.append(r["por_nivel"][L]["freq_std"])
        if not xs:
            continue
        cor = cmap(L / max(1, max_L))
        ax.errorbar(
            xs, ys, yerr=errs, marker="o", lw=1.8, capsize=3,
            color=cor, label=f"L={L} ({'borda' if L==0 else 'interior+'+str(L)})",
        )
    ax.axhline(0.25, color="black", linestyle="--", alpha=0.5,
               label="2/8 = 0.25 (esperado uniform)")
    ax.set_xticks(ns)
    ax.set_xlabel("n (lado do tabuleiro)")
    ax.set_ylabel("freq média da aresta")
    ax.set_title("Convergência de freq por nível de interioridade")
    ax.set_ylim(0, 0.55)
    ax.grid(linestyle=":", alpha=0.4)
    ax.legend(loc="upper right", fontsize=8, ncol=2)

    # Painel 2: freq[interior_max] vs n com referência
    ax2 = axes[1]
    # para cada n: freq média no MAIOR nível disponível (mais interior)
    max_L_n = [(r["n"], max(r["por_nivel"]),
                r["por_nivel"][max(r["por_nivel"])]) for r in resultados]
    xs = [x[0] for x in max_L_n]
    ys = [x[2]["freq_mean"] for x in max_L_n]
    es = [x[2]["freq_std"] for x in max_L_n]
    labels = [f"L_max={x[1]}, n_e={x[2]['n_arestas']}" for x in max_L_n]
    ax2.errorbar(xs, ys, yerr=es, marker="s", lw=2, capsize=4,
                 color="#1f4e79", label="freq em L_max(n)")
    for x, y, lab in zip(xs, ys, labels):
        ax2.annotate(lab, (x, y), xytext=(5, 7),
                     textcoords="offset points", fontsize=8)
    ax2.axhline(0.25, color="black", linestyle="--", alpha=0.5,
                label="0.25 (limite simétrico esperado)")
    ax2.set_xticks(xs)
    ax2.set_xlabel("n")
    ax2.set_ylabel("freq média no anel mais interior")
    ax2.set_title("Convergência: maior anel interior por n")
    ax2.set_ylim(0.15, 0.4)
    ax2.grid(linestyle=":", alpha=0.4)
    ax2.legend(loc="lower right", fontsize=9)

    fig.suptitle(
        "Existência empírica da fase local — freq[e] em vértices interiores",
        fontsize=12, y=1.02,
    )
    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"plot salvo: {out_path}")


if __name__ == "__main__":
    main()
