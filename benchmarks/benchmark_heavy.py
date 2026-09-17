
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

def run_heavy_benchmark():
    print("=========================================================")
    print(" BENCHMARK PESADO: 5000 Tours por Tabuleiro")
    print("=========================================================\n")
    
    K = 5000
    tabuleiros = [6, 8, 10, 12]
    
    for n in tabuleiros:
        # Pular 6x6 com K=5000 se só tiver 9862 e já testamos exaustivo? Não, 5000 é ok.
        print(f"--- Tabuleiro {n}x{n} (Buscando {K} tours) ---")
        
        # Testar Híbrida
        t0 = time.time()
        tours_hyb, nodes_hyb = kt_hybrid(n, K, seed=99)
        t1 = time.time()
        time_hyb = t1 - t0
        
        # Testar Original
        t0 = time.time()
        tours_orig = kt_orig(n, K, seed=99)
        t1 = time.time()
        time_orig = t1 - t0
        
        speedup = time_orig / time_hyb if time_hyb > 0 else float('inf')
        
        print(f"[Original] Tempo: {time_orig:.4f}s")
        print(f"[Híbrida]  Tempo: {time_hyb:.4f}s | Nós Explorados: {nodes_hyb}")
        print(f"-> Speedup: {speedup:.2f}x\n")

if __name__ == '__main__':
    run_heavy_benchmark()
