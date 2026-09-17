
# --- bootstrap: torna core/ importavel a partir de qualquer subpasta ---
import sys as _sys, pathlib as _pathlib
for _anc in _pathlib.Path(__file__).resolve().parents:
    if (_anc / "core").is_dir() and (_anc / "pyproject.toml").exists():
        _sys.path.insert(0, str(_anc / "core"))
        break
# --- fim do bootstrap ---
import time
from knight_tours import knight_tours as kt_orig
from knight_tours_hybrid import knight_tours as kt_hybrid

def run_exhaustive_benchmark():
    print("=========================================================")
    print(" BENCHMARK EXAUSTIVO 6x6: Original vs Híbrida (LUT)")
    print("=========================================================\n")
    
    # 6x6 tem 9862 tours estruturais únicos fechados.
    # Passando um K muito grande forçamos o Backtracking a varrer 100% da árvore
    n = 6
    K_INF = 100000 
    
    print(f"Iniciando busca exaustiva em {n}x{n}...\n")
    
    # 1. Testar Híbrida
    print("Executando: Heurística Híbrida (com LUT)...")
    t0 = time.time()
    tours_hyb, nodes_hyb = kt_hybrid(n, K_INF, seed=42)
    t1 = time.time()
    time_hyb = t1 - t0
    
    print(f" -> Tours encontrados: {len(tours_hyb)}")
    print(f" -> Tempo decorrido:   {time_hyb:.4f}s")
    print(f" -> Nós Explorados:    {nodes_hyb}\n")
    
    # 2. Testar Original
    print("Executando: Heurística Original (Estática)...")
    t0 = time.time()
    tours_orig = kt_orig(n, K_INF, seed=42)
    t1 = time.time()
    time_orig = t1 - t0
    
    print(f" -> Tours encontrados: {len(tours_orig)}")
    print(f" -> Tempo decorrido:   {time_orig:.4f}s\n")
    
    # Comparação
    speedup = time_orig / time_hyb if time_hyb > 0 else float('inf')
    print("---------------------------------------------------------")
    print(f" GANHO DE VELOCIDADE (SPEEDUP): {speedup:.2f}x")
    print("---------------------------------------------------------")

if __name__ == '__main__':
    run_exhaustive_benchmark()
