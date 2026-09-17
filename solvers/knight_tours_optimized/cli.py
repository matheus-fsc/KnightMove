"""cli.py — CLI interativa para o pacote knight_tours_optimized.

Uso:
    python -m knight_tours_optimized.cli

Funcionalidades:
  1) Gerar K tours fechados (amostragem)
  2) Enumerar TODOS os tours fechados (exaustivo)
  3) Caminho Hamiltoniano aberto start → end
  4) Decomposição D₄ (canônicos × expansão)
  5) Visualizar tour/caminho em ASCII
  6) Benchmark sequencial × paralelo
  7) Estimar tempo para enumerar TODOS (n > 6)
  0) Sair
"""

import json
import os
import sys
import time
from pathlib import Path

import numpy as np

from .tours import knight_tours, verify_tour
from .paths import knight_path, verify_path, _vertex_color
from .symmetry import (
    d4_symmetries, canonical_form, expand_d4, knight_tours_canonical,
)
from .parallel import knight_tours_parallel

try:
    from .core_numba import knight_tours_numba, NUMBA_AVAILABLE
    from .parallel_prefix import knight_tours_prefix_parallel
    OPTIMIZED_AVAILABLE = True
except ImportError:
    OPTIMIZED_AVAILABLE = False
    NUMBA_AVAILABLE = False


# Contagens conhecidas (literatura) para enumeração exaustiva
KNOWN_CLOSED_TOURS = {
    6: 9_862,
    8: 13_267_364_410_532,   # McKay 1997
    10: None,                 # desconhecido (estimado >> 10^22)
}

DATA_DIR = Path(__file__).parent / "data"


# ---------------------------------------------------------------------------
# Utilidades
# ---------------------------------------------------------------------------

def banner():
    print()
    print("╔" + "═" * 60 + "╗")
    print("║" + "  Knight's Tour Optimized — CLI Interativa".ljust(60) + "║")
    print("║" + "  D₄ · multiprocessing · open paths".ljust(60) + "║")
    print("╚" + "═" * 60 + "╝")


def ask(prompt, default=None, cast=str, validator=None):
    while True:
        suf = f" [{default}]" if default is not None else ""
        raw = input(f"{prompt}{suf}: ").strip()
        if not raw and default is not None:
            return default
        try:
            val = cast(raw)
        except (ValueError, TypeError):
            print(f"  valor inválido para {cast.__name__}")
            continue
        if validator and not validator(val):
            print("  valor fora do intervalo permitido")
            continue
        return val


def ask_yes_no(prompt, default=False):
    d = "S/n" if default else "s/N"
    raw = input(f"{prompt} [{d}]: ").strip().lower()
    if not raw:
        return default
    return raw in ("s", "sim", "y", "yes")


def fmt_int(x):
    return f"{x:,}".replace(",", ".")


def render_board(seq, n, closed=True):
    """Renderiza tabuleiro n×n com a ordem do tour em cada casa."""
    seq = [int(v) for v in seq]
    order = {v: i for i, v in enumerate(seq)}
    width = len(str(n * n)) + 1
    out = []
    for r in range(n):
        row = []
        for c in range(n):
            v = r * n + c
            row.append(str(order.get(v, "·")).rjust(width))
        out.append(" ".join(row))
    if closed:
        out.append(f"  (fechado: {seq[-1]} → {seq[0]})")
    else:
        out.append(f"  (aberto: {seq[0]} → {seq[-1]})")
    return "\n".join(out)


# ---------------------------------------------------------------------------
# Operações
# ---------------------------------------------------------------------------

def op_sample_tours(n):
    K = ask("Quantos tours (K)", default=100, cast=int, validator=lambda x: x > 0)
    seed_s = input("Seed (Enter = aleatório): ").strip()
    seed = int(seed_s) if seed_s else None

    print()
    print("Estratégias:")
    print("  1) NumPy sequencial")
    if NUMBA_AVAILABLE:
        print("  2) Numba sequencial (JIT, ~2× vs numpy)")
    print("  3) Paralelo por seeds (mais rápido p/ amostragem em n grande)")
    if OPTIMIZED_AVAILABLE:
        print("  4) Paralelo por prefixo (sem duplicatas, cap)")
        print("  5) Paralelo por prefixo + Numba (cap)")
        print("  6) Paralelo por prefixo + Numba EXHAUST (enumeração)")
    strategy = ask("Estratégia", default=1, cast=int,
                   validator=lambda x: x in (1, 2, 3, 4, 5, 6))

    n_workers = 1
    if strategy in (3, 4, 5, 6):
        n_workers = ask("Workers", default=8, cast=int,
                        validator=lambda x: x >= 1)

    print(f"\nGerando {K} tours em n={n}...")
    if strategy == 2 and NUMBA_AVAILABLE:
        knight_tours_numba(n, 5, seed=0)  # warmup
    t0 = time.perf_counter()
    if strategy == 1:
        tours = knight_tours(n, K, seed=seed)
    elif strategy == 2:
        tours = knight_tours_numba(n, K, seed=seed)
    elif strategy == 3:
        tours = knight_tours_parallel(n, K, n_workers=n_workers, seed=seed)
    elif strategy == 4:
        tours = knight_tours_prefix_parallel(n, K, depth=4,
                                              n_workers=n_workers, seed=seed,
                                              use_numba=False)
    elif strategy == 5:
        tours = knight_tours_prefix_parallel(n, K, depth=4,
                                              n_workers=n_workers, seed=seed,
                                              use_numba=True)
    elif strategy == 6:
        tours = knight_tours_prefix_parallel(n, K, depth=4,
                                              n_workers=n_workers, seed=seed,
                                              use_numba=True, exhaust=True)
    dt = time.perf_counter() - t0
    unique = len(set(tuple(int(x) for x in t) for t in tours))
    cov = unique / max(len(tours), 1) * 100
    print(f"  cobertura: {unique}/{len(tours)} únicos ({cov:.1f}%)")

    print(f"\n✓ {len(tours)} tours em {dt:.3f}s ({len(tours)/max(dt,1e-9):.0f} tours/s)")
    ok = all(verify_tour(t, n) for t in tours)
    print(f"  verify_tour em todos: {'OK ✓' if ok else 'FALHOU ✗'}")

    if tours and ask_yes_no("Visualizar um tour", default=False):
        idx = ask(f"Qual índice (0..{len(tours)-1})", default=0, cast=int,
                  validator=lambda x: 0 <= x < len(tours))
        print()
        print(render_board(tours[idx], n, closed=True))

    if tours and ask_yes_no("Salvar em JSON", default=False):
        DATA_DIR.mkdir(exist_ok=True)
        fname = DATA_DIR / f"tours_n{n}_K{len(tours)}.json"
        with open(fname, "w") as f:
            json.dump([t.tolist() for t in tours], f)
        print(f"  salvo em {fname}")


def op_enumerate_all(n):
    if n > 6:
        kn = KNOWN_CLOSED_TOURS.get(n)
        print(f"\n⚠  Enumeração exaustiva em n={n} pode demorar MUITO.")
        if kn:
            print(f"   Tours conhecidos: {fmt_int(kn)}")
        if not ask_yes_no("Continuar mesmo assim", default=False):
            return

    # Cap proporcional para garantir exaustão
    K_cap = 30_000 if n == 6 else 200_000
    K_cap = ask("Cap K (alto = exaustivo)", default=K_cap, cast=int,
                validator=lambda x: x > 0)
    print(f"\nEnumerando até esgotar (cap K={fmt_int(K_cap)})...")
    t0 = time.perf_counter()
    tours = knight_tours(n, K_cap, seed=0)
    dt = time.perf_counter() - t0

    unique = set(tuple(int(x) for x in t) for t in tours)
    print(f"\n  {fmt_int(len(tours))} tours em {dt:.3f}s")
    print(f"  únicos: {fmt_int(len(unique))}")
    kn = KNOWN_CLOSED_TOURS.get(n)
    if kn is not None:
        status = "✓ COMPLETO" if len(unique) == kn else f"✗ esperava {fmt_int(kn)}"
        print(f"  vs literatura: {status}")

    if unique and ask_yes_no("Computar decomposição D₄", default=True):
        syms = d4_symmetries(n)
        canon = set()
        for t in unique:
            canon.add(canonical_form(np.asarray(t, dtype=np.int32), syms))
        print(f"  órbitas D₄:  {fmt_int(len(canon))}")
        print(f"  razão total/canônicos: {len(unique)/max(len(canon),1):.3f}")


def op_open_path(n):
    V = n * n
    print(f"\n  vértices: 0..{V-1}; layout v = row*n + col")
    start = ask("start", default=0, cast=int,
                validator=lambda x: 0 <= x < V)
    end = ask("end", default=V - 1, cast=int,
              validator=lambda x: 0 <= x < V)
    K = ask("K (quantos caminhos)", default=1, cast=int,
            validator=lambda x: x > 0)
    seed_s = input("Seed (Enter = aleatório): ").strip()
    seed = int(seed_s) if seed_s else None

    cs, ce = _vertex_color(start, n), _vertex_color(end, n)
    print(f"\n  cor(start)={cs}, cor(end)={ce}, V={V} ({'par' if V%2==0 else 'ímpar'})")

    t0 = time.perf_counter()
    paths = knight_path(n, start, end, K=K, seed=seed)
    dt = time.perf_counter() - t0

    print(f"\n  {len(paths)} caminho(s) em {dt:.3f}s")
    if paths:
        ok = all(verify_path(p, n, start, end) for p in paths)
        print(f"  verify_path: {'OK ✓' if ok else 'FALHOU ✗'}")
        if ask_yes_no("Visualizar 1 caminho", default=True):
            print()
            print(render_board(paths[0], n, closed=False))


def op_d4_decomposition(n):
    if n > 8:
        print(f"⚠  D₄ exaustivo em n={n} requer enumerar TODOS os tours.")
        if not ask_yes_no("Continuar", default=False):
            return

    K_cap = 30_000 if n == 6 else 5_000
    K_cap = ask("Cap K (amostra/enumeração)", default=K_cap, cast=int,
                validator=lambda x: x > 0)

    print(f"\nGerando até {fmt_int(K_cap)} tours...")
    t0 = time.perf_counter()
    tours = knight_tours(n, K_cap, seed=0)
    t_gen = time.perf_counter() - t0
    unique = set(tuple(int(x) for x in t) for t in tours)

    t0 = time.perf_counter()
    syms = d4_symmetries(n)
    canon = set()
    for t in unique:
        canon.add(canonical_form(np.asarray(t, dtype=np.int32), syms))
    t_canon = time.perf_counter() - t0

    t0 = time.perf_counter()
    expanded = expand_d4(canon, n)
    t_expand = time.perf_counter() - t0

    print(f"\n  Geração ({fmt_int(len(unique))} tours únicos):  {t_gen:.3f}s")
    print(f"  Canonicalização ({fmt_int(len(canon))} órbitas): {t_canon:.3f}s")
    print(f"  Expansão D₄ ({fmt_int(len(expanded))} reproduzidos): {t_expand:.3f}s")
    print(f"  Razão: {len(unique)/max(len(canon),1):.3f}× (máx 8)")
    print(f"  Reconstrução consistente: {'✓' if expanded == unique else '✗'}")

    # Distribuição de tamanhos de órbitas
    orbit_size = {}
    for ct in canon:
        ct_arr = np.asarray(ct, dtype=np.int32)
        images = set()
        for perm in syms:
            transformed = perm[ct_arr]
            from .symmetry import _normalize_closed_tour
            images.add(_normalize_closed_tour(transformed))
        orbit_size[len(images)] = orbit_size.get(len(images), 0) + 1
    print(f"\n  Distribuição de tamanhos de órbita:")
    for sz in sorted(orbit_size):
        print(f"    tamanho {sz}: {fmt_int(orbit_size[sz])} órbitas")


def op_visualize_tour(n):
    print("\n1) Gerar tour e visualizar")
    print("2) Visualizar caminho aberto start→end")
    op = ask("Escolha", default=1, cast=int, validator=lambda x: x in (1, 2))
    seed_s = input("Seed (Enter = aleatório): ").strip()
    seed = int(seed_s) if seed_s else None

    if op == 1:
        tours = knight_tours(n, 1, seed=seed)
        if not tours:
            print("  ✗ nenhum tour")
            return
        print()
        print(render_board(tours[0], n, closed=True))
    else:
        V = n * n
        start = ask("start", default=0, cast=int,
                    validator=lambda x: 0 <= x < V)
        end = ask("end", default=V - 1, cast=int,
                  validator=lambda x: 0 <= x < V)
        paths = knight_path(n, start, end, K=1, seed=seed)
        if not paths:
            print("  ✗ nenhum caminho")
            return
        print()
        print(render_board(paths[0], n, closed=False))


def op_benchmark(n):
    K = ask("K para benchmark", default=5000, cast=int,
            validator=lambda x: x > 0)
    workers = ask("Workers paralelos", default=8, cast=int,
                  validator=lambda x: x >= 1)
    if OPTIMIZED_AVAILABLE:
        from .parallel_prefix import full_benchmark_optimized
        full_benchmark_optimized(n=n, K=K, n_workers=workers)
        return

    print(f"\nSequencial (n={n}, K={K})...")
    t0 = time.perf_counter()
    knight_tours(n, K, seed=0)
    t_seq = time.perf_counter() - t0
    print(f"  {t_seq:.3f}s  ({K/max(t_seq,1e-9):.0f} tours/s)")
    print(f"Paralelo ({workers} workers)...")
    t0 = time.perf_counter()
    knight_tours_parallel(n, K, n_workers=workers, seed=0)
    t_par = time.perf_counter() - t0
    print(f"  {t_par:.3f}s  ({K/max(t_par,1e-9):.0f} tours/s)")
    print(f"\n  Speedup: {t_seq/max(t_par,1e-9):.2f}×")
    print(f"  Eficiência: {(t_seq/max(t_par,1e-9))/workers*100:.0f}%")


def op_estimate_total(n):
    """Estima tempo para enumerar TODOS os tours, medindo cada estratégia."""
    K = 20000 if n >= 8 else 5000
    print(f"\nMedindo throughput em n={n} (K={K})...\n")

    rates = {}

    # NumPy seq
    t0 = time.perf_counter()
    knight_tours(n, K, seed=0)
    t = time.perf_counter() - t0
    rates['NumPy seq'] = K / max(t, 1e-9)

    # Numba seq
    if NUMBA_AVAILABLE:
        knight_tours_numba(n, 5, seed=0)  # warmup
        t0 = time.perf_counter()
        knight_tours_numba(n, K, seed=0)
        t = time.perf_counter() - t0
        rates['Numba seq'] = K / max(t, 1e-9)

    # Par seeds 8w
    t0 = time.perf_counter()
    knight_tours_parallel(n, K, n_workers=8, seed=0)
    t = time.perf_counter() - t0
    rates['Par seeds 8w (amostra)'] = K / max(t, 1e-9)

    # Par prefix 8w + Numba (cobertura única)
    if OPTIMIZED_AVAILABLE:
        t0 = time.perf_counter()
        tours = knight_tours_prefix_parallel(n, K, depth=4, n_workers=8,
                                              seed=0, use_numba=True)
        t = time.perf_counter() - t0
        unique = len(set(tuple(int(x) for x in t) for t in tours))
        rates['Par prefix+Numba 8w (únicos)'] = unique / max(t, 1e-9)

    for label, r in rates.items():
        print(f"  {label:35s} {r:>10.0f} tours/s")

    print()
    print("Cenários de extrapolação:")
    for tour_label, total in [
        ("6×6 (9.862)", 9_862),
        ("8×8 ~ 3 trilhões (estim. usuário)", 3_000_000_000_000),
        ("8×8 ~ 1.3 × 10¹³ (McKay 1997)", 13_267_364_410_532),
    ]:
        print(f"\n  {tour_label}:")
        for label, r in rates.items():
            t_total = total / r
            print(f"    {label:35s} {fmt_time(t_total)}")


def fmt_time(seconds):
    if seconds < 60:
        return f"{seconds:.1f}s"
    if seconds < 3600:
        return f"{seconds/60:.1f}min"
    if seconds < 86400:
        return f"{seconds/3600:.1f}h"
    if seconds < 86400 * 365:
        return f"{seconds/86400:.1f} dias"
    return f"{seconds/86400/365:.2f} anos"


# ---------------------------------------------------------------------------
# Loop principal
# ---------------------------------------------------------------------------

MENU = """
Operações:
  1) Gerar K tours fechados (amostragem)
  2) Enumerar TODOS os tours fechados
  3) Caminho Hamiltoniano aberto start→end
  4) Decomposição D₄ (canônicos × expansão)
  5) Visualizar tour/caminho em ASCII
  6) Benchmark sequencial × paralelo
  7) Estimar tempo de enumeração completa
  9) Trocar tamanho n
  0) Sair
"""

OPS = {
    1: op_sample_tours,
    2: op_enumerate_all,
    3: op_open_path,
    4: op_d4_decomposition,
    5: op_visualize_tour,
    6: op_benchmark,
    7: op_estimate_total,
}


def main():
    banner()
    n = ask("Tamanho do tabuleiro n", default=6, cast=int,
            validator=lambda x: x >= 4)
    while True:
        print(MENU)
        op = ask("Escolha", default=1, cast=int,
                 validator=lambda x: x in (0, 1, 2, 3, 4, 5, 6, 7, 9))
        if op == 0:
            print("Tchau!")
            return
        if op == 9:
            n = ask("Novo n", default=n, cast=int, validator=lambda x: x >= 4)
            continue
        fn = OPS.get(op)
        if fn is None:
            print("inválido")
            continue
        try:
            fn(n)
        except KeyboardInterrupt:
            print("\n[interrompido]")
        except Exception as exc:
            print(f"\n✗ erro: {exc}")
        if not ask_yes_no("\nContinuar", default=True):
            return


if __name__ == "__main__":
    try:
        main()
    except (EOFError, KeyboardInterrupt):
        print("\nTchau!")
