
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

def run_benchmark():
    print("=========================================================")
    print(" BENCHMARK: Heurística Original vs Híbrida Contínua")
    print("=========================================================\n")
    
    # K fixo para teste
    K = 100
    tabuleiros = [6, 8, 10, 12]
    
    for n in tabuleiros:
        print(f"--- Tabuleiro {n}x{n} (Buscando {K} tours) ---")
        
        # Testar Híbrida primeiro
        t0 = time.time()
        # knight_tours_hybrid retorna (tours, nodes)
        tours_hyb, nodes_hyb = kt_hybrid(n, K, seed=42)
        t1 = time.time()
        time_hyb = t1 - t0
        
        # Testar Original
        t0 = time.time()
        # knight_tours original não retorna os nós, precisamos medir apenas o tempo
        tours_orig = kt_orig(n, K, seed=42)
        t1 = time.time()
        time_orig = t1 - t0
        
        # Como o código original não exportava os nós, podemos avaliar a aceleração (speedup) de tempo
        speedup = time_orig / time_hyb if time_hyb > 0 else float('inf')
        
        print(f"[Original] Tempo: {time_orig:.4f}s")
        print(f"[Híbrida]  Tempo: {time_hyb:.4f}s | Nós Explorados: {nodes_hyb}")
        print(f"-> Speedup: {speedup:.2f}x\n")

if __name__ == '__main__':
    run_benchmark()
