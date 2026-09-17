"""sheared_winding.py — Decodificador de windings no toro com cisalhamento.

Para cada tour fechado válido, decodifica os winding numbers (w_x, w_y) e
paridades. O cisalhamento é compensado durante a decodificação.

Predição teórica:
  - Se s = n/2 (n par): Q > 0 (classe de paridade forçada)
  - Se gcd(s, n) = 1: Q = 0 (ambas paridades realizadas)

Interface pública:
    decode_sheared_winding(tour, n, m, s) -> (w_x, w_y, ambig_count)
    medir_sheared(n, m, s, K_por_seed, n_seeds, label="")
    batch_simulation(board_configs)
"""

# --- bootstrap: torna core/ importavel a partir de qualquer subpasta ---
import sys as _sys, pathlib as _pathlib
for _anc in _pathlib.Path(__file__).resolve().parents:
    if (_anc / "core").is_dir() and (_anc / "pyproject.toml").exists():
        _sys.path.insert(0, str(_anc / "core"))
        break
# --- fim do bootstrap ---

import time
from collections import Counter
from math import gcd

import numpy as np

from knight_tours_sheared import knight_tours, verify_tour, sheared_step, _KNIGHT_MOVES


def _decode_dx(dx_aparente, m):
    """Decodifica dx_aparente ∈ ℤ/m para dx_real ∈ {-2, -1, 1, 2}.
    
    Se múltiplas realizações, retorna (dx_real, ambig=True).
    Se nenhuma, retorna (None, ambig=True) indicando erro.
    """
    candidates = []
    for dx_real in (-2, -1, 1, 2):
        if (dx_real % m) == (dx_aparente % m):
            candidates.append(dx_real)
    
    if not candidates:
        return None, True  # Ambiguous/error case
    
    if len(candidates) > 1:
        return candidates[0], True
    return candidates[0], False


def _decode_dy(dy_base, n):
    """Decodifica dy_base ∈ ℤ/n para dy_real ∈ {-2, -1, 1, 2}.
    
    Se múltiplas realizações, retorna (dy_real, ambig=True).
    Se nenhuma, retorna (None, ambig=True) indicando erro.
    """
    candidates = []
    for dy_real in (-2, -1, 1, 2):
        if (dy_real % n) == (dy_base % n):
            candidates.append(dy_real)
    
    if not candidates:
        return None, True  # Ambiguous/error case
    
    if len(candidates) > 1:
        return candidates[0], True
    return candidates[0], False


def decode_sheared_winding(tour, n, m, s):
    """Decodifica winding numbers do tour no toro com cisalhamento s.
    
    Algoritmo:
      Para cada passo v_i -> v_{i+1}:
        1. Compute dx_aparente = x_{i+1} - x_i (mod m)
        2. Decodifique dx_real, compute wrap_x = (x_i + dx_real) // m
        3. Compute dy_base = (y_{i+1} - y_i - wrap_x * s) mod n
        4. Decodifique dy_real
        5. Accumule somas
      6. Retorne (w_x, w_y) = (sum_dx // m, sum_dy // n)
    
    Parâmetros:
        tour: lista/array de vértices
        n: número de linhas
        m: número de colunas
        s: parâmetro de cisalhamento
    
    Retorna:
        (w_x, w_y, ambig_count): winding numbers e número de passos ambíguos
    
    Levanta:
        RuntimeError se sanity checks falham
        ValueError se decodificação ambígua/inválida
    """
    V = n * m
    sum_dx = 0
    sum_dy = 0
    ambig_count = 0
    
    for i in range(V):
        u = int(tour[i])
        v = int(tour[(i + 1) % V])
        
        uy, ux = divmod(u, m)
        vy, vx = divmod(v, m)
        
        # Step 1: Decodifique dx
        dx_aparente = (vx - ux) % m
        dx_real, dx_ambig = _decode_dx(dx_aparente, m)
        if dx_real is None:
            raise ValueError(f"passo {i}: impossível decodificar dx_aparente={dx_aparente} em mod {m}")
        if dx_ambig:
            ambig_count += 1
        
        # Step 2: Compute wrap_x
        wrap_x = (ux + dx_real) // m
        
        # Step 3: Decodifique dy (com compensação de cisalhamento)
        dy_base = (vy - uy - wrap_x * s) % n
        dy_real, dy_ambig = _decode_dy(dy_base, n)
        if dy_real is None:
            raise ValueError(f"passo {i}: impossível decodificar dy_base={dy_base} em mod {n}")
        if dy_ambig:
            ambig_count += 1
        
        sum_dx += dx_real
        sum_dy += dy_real
    
    # Sanity checks
    if sum_dx % m != 0:
        raise RuntimeError(
            f"sanity falhou: Σdx={sum_dx} não é múltiplo de m={m}")
    adjusted_y = sum_dy + s * (sum_dx // m)
    if adjusted_y % n != 0:
        raise RuntimeError(
            f"sanity falhou: Σdy ajustado={adjusted_y} não é múltiplo de n={n}")
    
    w_x = sum_dx // m
    w_y = adjusted_y // n
    
    return w_x, w_y, ambig_count



def dedup_tour_key(t, n, m):
    """Gera chave para deduplicação de tours (normalização por simetria de arestas)."""
    V = n * m
    ar = sorted(tuple(sorted([int(t[i]), int(t[(i + 1) % V])]))
                for i in range(V))
    return tuple(ar)


def medir_sheared(n, m, s, K_por_seed, n_seeds, label=""):
    """Mede distribuição de paridades em tours do toro com cisalhamento.
    
    Parâmetros:
        n: número de linhas
        m: número de colunas
        s: parâmetro de cisalhamento
        K_por_seed: número de tours por seed
        n_seeds: número de seeds
        label: rótulo para impressão
    
    Retorna:
        (par_conj, w_dist, todos_tours): contador de paridades, winding dist, lista de tours
    """
    chaves = set()
    par_conj = Counter()
    w_dist = Counter()
    todos_tours = []
    
    t0 = time.time()
    bruto = 0
    erros = 0
    
    for seed_idx in range(n_seeds):
        tours = knight_tours(n, m, s, K_por_seed, seed=seed_idx)
        bruto += len(tours)
        
        for tour in tours:
            assert verify_tour(tour, n, m, s), f"tour inválido (seed={seed_idx})"
            
            k = dedup_tour_key(tour, n, m)
            if k in chaves:
                continue
            chaves.add(k)
            
            try:
                w_y, w_x, ambig = decode_sheared_winding(tour, n, m, s)
                todos_tours.append(tour)
                par_conj[(abs(w_y) % 2, abs(w_x) % 2)] += 1
                w_dist[(w_y, w_x)] += 1
            except (ValueError, RuntimeError) as e:
                erros += 1
                if erros <= 3:  # Log primeiros 3 erros
                    print(f"  ⚠ Erro ao decodificar tour: {e}")
    
    t1 = time.time()
    
    total = sum(par_conj.values())
    
    # Computar Q (número de classes de paridade realizadas)
    Q = len([c for c in par_conj if par_conj[c] > 0])
    
    # Computar gcd(s, n) para predição teórica
    g = gcd(s, n)
    s_pred = "Q=0 (teórico)" if g == 1 else "Q>0 (teórico)" if s == n // 2 else "indefinido"
    
    print(f"\n{'='*68}")
    print(f"  {label}  [toro {n}x{m}, s={s}]   {n_seeds} seeds × K={K_por_seed}")
    print(f"  gcd(s, n) = {g}  →  {s_pred}")
    print(f"{'='*68}")
    print(f"  tempo: {t1-t0:.2f}s  brutos={bruto}  únicos={total}  erros={erros}")
    print(f"  Q (classes de paridade): {Q}/4")
    print(f"  distribuição conjunta (par_y, par_x):")
    
    classes = [(0, 0), (0, 1), (1, 0), (1, 1)]
    for c in classes:
        cnt = par_conj.get(c, 0)
        pct = 100 * cnt / total if total else 0
        print(f"    (y={c[0]}, x={c[1]}) : {cnt:6d}   {pct:5.2f}%")
    
    if total > 0:
        sigma = (total * 0.25 * 0.75) ** 0.5
        print(f"  σ esperado por classe sob H0=25%: ±{sigma:.1f}")
        for c in classes:
            cnt = par_conj.get(c, 0)
            z = (cnt - total / 4) / sigma if sigma > 0 else 0
            print(f"    z(y={c[0]},x={c[1]}) = {z:+.2f}")
    
    # Estatísticas de winding
    w_x_vals = sorted(set(w for w, _ in w_dist))
    w_y_vals = sorted(set(w for _, w in w_dist))
    print(f"  w_x range: [{min(w_x_vals)}, {max(w_x_vals)}]" if w_x_vals else "  (sem dados)")
    print(f"  w_y range: [{min(w_y_vals)}, {max(w_y_vals)}]" if w_y_vals else "  (sem dados)")
    
    return par_conj, w_dist, todos_tours


def batch_simulation(board_configs):
    """Executa simulação em lote com múltiplas configurações.
    
    Parâmetro:
        board_configs: lista de tuplas (n, m, s_values, K_por_seed, n_seeds)
                      onde s_values é lista de valores de s para testar
    
    Exemplo:
        configs = [
            (6, 6, range(1, 6), 10, 2),      # 6x6, s∈{1,2,3,4,5}, 10 tours, 2 seeds
            (8, 6, [1, 2, 3, 4, 5, 6, 7], 8, 2),
        ]
        batch_simulation(configs)
    """
    print("\n" + "="*68)
    print("  BATCH SIMULATION: TORO COM CISALHAMENTO")
    print("="*68)
    
    for n, m, s_values, K_por_seed, n_seeds in board_configs:
        print(f"\nTabuleiro {n}×{m}:")
        for s in s_values:
            s_norm = s % n
            label = f"Cisalhamento s={s_norm}"
            try:
                par_conj, w_dist, _ = medir_sheared(
                    n, m, s_norm, K_por_seed, n_seeds, label=label
                )
            except Exception as e:
                print(f"  ✗ Erro na simulação: {e}")
                import traceback
                traceback.print_exc()


if __name__ == "__main__":
    # Exemplo: simular 6x6 e 8x6 com diversos s
    configs = [
        (6, 6, range(1, 6), 10, 2),
        (8, 6, range(1, 6), 8, 2),
    ]
    batch_simulation(configs)
