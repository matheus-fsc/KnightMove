
# --- bootstrap: torna core/ importavel a partir de qualquer subpasta ---
import sys as _sys, pathlib as _pathlib
for _anc in _pathlib.Path(__file__).resolve().parents:
    if (_anc / "core").is_dir() and (_anc / "pyproject.toml").exists():
        _sys.path.insert(0, str(_anc / "core"))
        break
# --- fim do bootstrap ---
# bias_check_6x6.py
import random
from collections import Counter
from scipy.stats import chisquare
import knight_tours_sheared as kt
import sheared_winding as sw

random.seed(42)

N, M = 6, 6
TARGET = 10_000

def collect_tours(s, target):
    """Coleta 'target' tours únicos para o shear s."""
    tours = set()
    attempts = 0
    while len(tours) < target:
        tour = kt.generate_tour(s, N, M)
        if tour is not None:
            # Converte para tupla de tuplas (imutável) para o set
            tour_tuple = tuple(tour)
            if tour_tuple not in tours:
                tours.add(tour_tuple)
        attempts += 1
        if attempts % 1000 == 0:
            print(f"  [s={s}] {len(tours)} tours coletados ({attempts} tentativas)...")
    print(f"  [s={s}] Coleta concluída: {len(tours)} tours em {attempts} tentativas.")
    return list(tours)

def compute_parity_counts(tours, s):
    """Retorna dicionário com contagens das 4 classes de paridade."""
    counts = {(0,0):0, (0,1):0, (1,0):0, (1,1):0}
    for tour in tours:
        # Fecha o ciclo para cálculo do winding
        closed = list(tour) + [tour[0]]
        wx_mod, wy_mod = sw.get_winding_parity(closed, N, M, s)
        counts[(wx_mod, wy_mod)] += 1
    return counts

def print_table(s, total, counts):
    print(f"\n--- Resultados para s = {s} ---")
    print(f"Total de tours únicos: {total}")
    print(f"{'Classe':>8} | {'Contagem':>8} | {'Percentagem':>11}")
    print("-" * 38)
    for (wx, wy), cnt in sorted(counts.items()):
        pct = 100.0 * cnt / total
        print(f"({wx},{wy})   | {cnt:8d} | {pct:10.2f}%")

    # Determina classes não‑nulas e faz teste χ²
    non_zero = {k: v for k, v in counts.items() if v > 0}
    if len(non_zero) > 1:
        obs = list(non_zero.values())
        expected = [total / len(non_zero)] * len(non_zero)
        chi2, p = chisquare(f_obs=obs, f_exp=expected)
        print(f"\nTeste χ² para uniformidade nas classes não‑nulas:")
        print(f"  χ² = {chi2:.4f}, valor‑p = {p:.4f}")
        if p < 0.05:
            print("  → VIÉS SIGNIFICATIVO (rejeita‑se H0 a 5%)")
        else:
            print("  → Distribuição compatível com uniformidade.")
    else:
        print("\n  Apenas uma classe não‑nula – teste χ² não aplicável.")

# Cenário A: s = 1
print("Coletando tours para s = 1 (desacoplamento, Q=0 esperado)...")
tours_s1 = collect_tours(s=1, target=TARGET)
counts_s1 = compute_parity_counts(tours_s1, s=1)
print_table(s=1, total=len(tours_s1), counts=counts_s1)

# Cenário B: s = 3
print("\nColetando tours para s = 3 (colapso perfeito, Q>0 esperado)...")
tours_s3 = collect_tours(s=3, target=TARGET)
counts_s3 = compute_parity_counts(tours_s3, s=3)
print_table(s=3, total=len(tours_s3), counts=counts_s3)

# Verificação explícita de que classes com wy ímpar são zero para s=3
print("\nVerificação para s=3:")
if counts_s3[(0,1)] == 0 and counts_s3[(1,1)] == 0:
    print("  ✓ Classes com wy ímpar são EXATAMENTE zero (colapso topológico confirmado).")
else:
    print("  ✗ ERRO: tours com wy ímpar encontrados – desvio da topologia esperada.")
