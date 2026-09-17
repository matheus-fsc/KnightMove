import json
import time
from klein_full_search import build_graph, State, Metrics, search

def run_search(n, m, node_limit, enable_topo_prune):
    t0 = time.perf_counter()
    ctx = build_graph(n, m)
    ctx["enable_topo_prune"] = enable_topo_prune
    state = State(int(ctx["V"]), int(ctx["E"]))
    metrics = Metrics()
    canonical_seen = set()

    search(state, ctx, metrics, canonical_seen, 0, None, node_limit)
    
    elapsed = time.perf_counter() - t0
    return elapsed, metrics

def main():
    n, m = 6, 6
    node_limit = 500_000
    
    print("============================================================")
    print("BENCHMARK DE RUPTURA TOPOLÓGICA: GARRAFA DE KLEIN (6x6)")
    print("============================================================")
    print("Executando Rodada A (Com TOPO_PRUNE)...")
    time_a, metrics_a = run_search(n, m, node_limit, True)
    
    print("Executando Rodada B (Sem TOPO_PRUNE)...")
    time_b, metrics_b = run_search(n, m, node_limit, False)
    
    speedup = time_b / time_a if time_a > 0 else float('inf')
    
    print("\n============================================================")
    print("Métrica               | Com TOPO_PRUNE | Sem TOPO_PRUNE | Delta / Speedup")
    print("------------------------------------------------------------")
    print(f"Nós Visitados         | {metrics_a.nodes:<14} | {metrics_b.nodes:<14} | -")
    print(f"Tempo Total (s)       | {time_a:<14.3f} | {time_b:<14.3f} | {speedup:.2f}x Speedup")
    print(f"Tours Válidos / nós   | {metrics_a.canonical_tours:<14} | {metrics_b.canonical_tours:<14} | -")
    print(f"Podas Topológicas     | {metrics_a.topo_prunes:<14} | 0 (Desativado) | -")
    print("============================================================")
    
    results = {
        "board": f"{n}x{m}",
        "node_limit": node_limit,
        "round_a_topo_on": {
            "time_s": time_a,
            "nodes": metrics_a.nodes,
            "tours": metrics_a.canonical_tours,
            "topo_prunes": metrics_a.topo_prunes,
            "subtour_prunes": metrics_a.subtour_prunes
        },
        "round_b_topo_off": {
            "time_s": time_b,
            "nodes": metrics_b.nodes,
            "tours": metrics_b.canonical_tours,
            "topo_prunes": 0,
            "subtour_prunes": metrics_b.subtour_prunes
        },
        "speedup": speedup
    }
    
    with open("klein_speedup_results.json", "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
        
    print("\nResultados salvos em 'klein_speedup_results.json'.")

if __name__ == "__main__":
    main()
