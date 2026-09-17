"""
Medição do TAMANHO do espaço de estados da T_bulk de graus (column-transfer).
=============================================================================

SÓ ENUMERA ESTADOS. Não monta matriz, não calcula espectro, não conta nada.
Objetivo: saber quantos estados a T_bulk de graus tem em n=8, para decidir
viabilidade de hardware.

Estado = (deg_c, deg_{c+1}), deg ∈ {0,1,2}^n. Reusa `transitions` /
`forward_options` de transfer_build.py (a MESMA BFS que produziu 33.346 em n=6).
Estado guardado BITPACKED em int (2 bits/célula): n=8 -> 32 bits.

Lembrete: este é o objeto CEGO À CONECTIVIDADE (cresce mais devagar). NÃO é o
broken-profile de knight_transfer.cpp (que estoura >187M em n=8).
"""
from __future__ import annotations
import resource
import sys
import time
from collections import deque
from pathlib import Path

THIS = Path(__file__).parent
sys.path.insert(0, str(THIS))
from transfer_build import transitions  # a mesma usada em build_T


def pack(s, n):
    deg_c, deg_cp1 = s
    x = 0
    for d in deg_c:
        x = (x << 2) | d
    for d in deg_cp1:
        x = (x << 2) | d
    return x


def unpack(x, n):
    cells = []
    for _ in range(2 * n):
        cells.append(x & 3)
        x >>= 2
    cells.reverse()
    return (tuple(cells[:n]), tuple(cells[n:]))


def rss_gb():
    # ru_maxrss em KB no Linux
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024 * 1024)


def enumerate_states(n, mem_cap_gb=14.0, time_cap_s=1500.0,
                     report_every=1_000_000, verbose=True):
    """BFS de estados alcançáveis a partir de s0=(0,0). Retorna dict de stats.

    Aborta graciosamente se RSS > mem_cap_gb ou tempo > time_cap_s.
    """
    s0 = ((0,) * n, (0,) * n)
    p0 = pack(s0, n)
    visited = {p0}
    queue = deque([p0])

    t0 = time.time()
    last_report = 0
    pops = 0
    curve = []  # (n_states, elapsed_s, rss_gb)
    status = "complete"

    while queue:
        p = queue.popleft()
        s = unpack(p, n)
        for s_next in transitions(s, n, forbid_cp2=False):
            pn = pack(s_next, n)
            if pn not in visited:
                visited.add(pn)
                queue.append(pn)
        pops += 1

        if pops % 20000 == 0:  # checar limites em intervalo fixo de pops
            el = time.time() - t0
            mem = rss_gb()
            if mem > mem_cap_gb:
                status = "mem_cap"
                break
            if el > time_cap_s:
                status = "time_cap"
                break

        if len(visited) - last_report >= report_every:
            last_report = len(visited)
            el = time.time() - t0
            mem = rss_gb()
            curve.append((len(visited), el, mem))
            if verbose:
                print(f"    estados={len(visited):>12,}  fila={len(queue):>12,}  "
                      f"processados={pops:>12,}  t={el:7.1f}s  rss={mem:5.2f}GB",
                      flush=True)

    el = time.time() - t0
    mem = rss_gb()
    curve.append((len(visited), el, mem))
    return dict(n=n, n_states=len(visited), processed=pops,
                queue_left=len(queue), elapsed_s=el, rss_gb=mem,
                status=status, curve=curve)


def main():
    print("=" * 70)
    print("PASSO 1 — Gate de fidelidade: reproduzir 33.346 estados em n=6")
    print("=" * 70)
    r6 = enumerate_states(6, verbose=False)
    print(f"  n=6: estados enumerados = {r6['n_states']:,}  "
          f"(esperado 33.346)  status={r6['status']}")
    print(f"       tempo={r6['elapsed_s']:.2f}s  rss={r6['rss_gb']:.2f}GB")
    if r6["n_states"] != 33346:
        print("  !!! GATE FALHOU — enumeração não confiável, abortando.")
        sys.exit(1)
    print("  >> GATE OK\n")

    print("=" * 70)
    print("PASSO 2 (intermediário, âncora) — n=7  [deixar terminar]")
    print("=" * 70)
    r7 = enumerate_states(7, mem_cap_gb=14.0, time_cap_s=1000.0,
                          report_every=200_000, verbose=True)
    print(f"  n=7: estados={r7['n_states']:,}  status={r7['status']}  "
          f"t={r7['elapsed_s']:.1f}s  rss={r7['rss_gb']:.2f}GB  "
          f"fila={r7['queue_left']:,}\n")

    print("=" * 70)
    print("PASSO 2 — n=8 (o número que importa)  [cap 14GB, ~12min]")
    print("=" * 70)
    r8 = enumerate_states(8, mem_cap_gb=14.0, time_cap_s=720.0,
                          report_every=500_000, verbose=True)
    print(f"\n  n=8: estados={r8['n_states']:,}  status={r8['status']}")
    print(f"       processados={r8['processed']:,}  fila_restante={r8['queue_left']:,}")
    print(f"       tempo={r8['elapsed_s']:.1f}s  rss={r8['rss_gb']:.2f}GB")

    import json
    out = dict(n4_5_6_ref={4: 342, 5: 7322, 6: 33346},
               n6=r6, n7=r7, n8=r8)
    json.dump(out, open(THIS / "data" / "measure_states.json", "w"), indent=2)
    print("\n  salvo em data/measure_states.json")


if __name__ == "__main__":
    main()
