import json
import time
from klein_predictive_search import build_graph, State, Metrics, search

def run_search(n, m, node_limit, enable_lookahead, enable_heuristic, use_strict_limit=False):
    t0 = time.perf_counter()
    ctx = build_graph(n, m)
    ctx["enable_lookahead"] = enable_lookahead
    ctx["enable_heuristic"] = enable_heuristic
    ctx["use_strict_limit"] = use_strict_limit
    
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
    print("BENCHMARK PREDICTIVE & HEURISTIC KLEIN (6x6)")
    print("============================================================")
    
    print("Executando Rodada A (Clássica: Sem Look-ahead, Sem Heurística)...")
    time_a, metrics_a = run_search(n, m, node_limit, False, False)
    
    print("Executando Rodada B (Otimizada: Com Look-ahead, Com Heurística)...")
    time_b, metrics_b = run_search(n, m, node_limit, True, True)
    
    print("\n============================================================")
    print("Métrica                         | Rodada A (Clássico) | Rodada B (Otimizado)")
    print("------------------------------------------------------------")
    print(f"Nós Visitados                   | {metrics_a.nodes:<19} | {metrics_b.nodes:<19}")
    print(f"Tempo Total (s)                 | {time_a:<19.3f} | {time_b:<19.3f}")
    print(f"Tours Válidos (parcial/limite)  | {metrics_a.canonical_tours:<19} | {metrics_b.canonical_tours:<19}")
    print(f"Poda Topológica Clássica (Ciclo)| {metrics_a.topo_prunes:<19} | {metrics_b.topo_prunes:<19}")
    print(f"Poda Preditiva (Look-ahead)     | {metrics_a.predictive_prunes:<19} | {metrics_b.predictive_prunes:<19}")
    print(f"Subciclos rejeitados pelo UF    | {metrics_a.subtour_prunes:<19} | {metrics_b.subtour_prunes:<19}")
    print(f"Contradições locais             | {metrics_a.contradictions:<19} | {metrics_b.contradictions:<19}")
    
    tour_rate_a = metrics_a.canonical_tours / time_a if time_a > 0 else 0
    tour_rate_b = metrics_b.canonical_tours / time_b if time_b > 0 else 0
    heuristic_efficiency = tour_rate_b / tour_rate_a if tour_rate_a > 0 else 0
    
    print("------------------------------------------------------------")
    print(f"Taxa de Soluções (tours/seg)   | {tour_rate_a:<19.1f} | {tour_rate_b:<19.1f}")
    print(f"Eficiência da Heurística (Speedup na busca de soluções)      : {heuristic_efficiency:.2f}x")
    print("============================================================")
    
    results = {
        "board": f"{n}x{m}",
        "node_limit": node_limit,
        "round_a_classical": {
            "time_s": time_a,
            "nodes": metrics_a.nodes,
            "tours": metrics_a.canonical_tours,
            "topo_prunes": metrics_a.topo_prunes,
            "predictive_prunes": metrics_a.predictive_prunes,
            "subtour_prunes": metrics_a.subtour_prunes,
            "contradictions": metrics_a.contradictions,
            "tour_rate": tour_rate_a
        },
        "round_b_optimized": {
            "time_s": time_b,
            "nodes": metrics_b.nodes,
            "tours": metrics_b.canonical_tours,
            "topo_prunes": metrics_b.topo_prunes,
            "predictive_prunes": metrics_b.predictive_prunes,
            "subtour_prunes": metrics_b.subtour_prunes,
            "contradictions": metrics_b.contradictions,
            "tour_rate": tour_rate_b
        },
        "heuristic_efficiency": heuristic_efficiency
    }
    
    with open("benchmark_predictive_results.json", "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
        
    print("\nResultados do benchmark salvos em 'benchmark_predictive_results.json'.")

if __name__ == "__main__":
    main()
