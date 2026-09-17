import sys
sys.path.append("/home/math/Dev/knight_tour")

from klein_predictive_search import build_graph, State, Metrics, search

def test_board(n, m, node_limit):
    ctx = build_graph(n, m)
    ctx["enable_lookahead"] = True
    ctx["enable_heuristic"] = True
    ctx["use_strict_limit"] = False # Test standard V - d first
    
    state = State(int(ctx["V"]), int(ctx["E"]))
    metrics = Metrics()
    canonical_seen = set()

    search(state, ctx, metrics, canonical_seen, 0, None, node_limit)
    return metrics

def main():
    boards = [
        (6, 8),
        (6, 10),
        (6, 12),
        (6, 14)
    ]
    node_limit = 200_000
    
    print("========================================================================")
    print("EXPERIMENTO: O DESPERTAR DA PODA PREDITIVA EM GRELAH RECTANGULAR")
    print("========================================================================")
    for n, m in boards:
        print(f"\nTestando tabuleiro {n}x{m} (V={n*m}) com limite de {node_limit} nós...")
        try:
            metrics = test_board(n, m, node_limit)
            print(f"  Nós Visitados                : {metrics.nodes}")
            print(f"  Poda Preditiva Ativa         : {metrics.predictive_prunes}")
            print(f"  Poda Clássica (Ciclo)        : {metrics.topo_prunes}")
            print(f"  Tours Encontrados (Canônicos): {metrics.canonical_tours}")
        except Exception as e:
            print(f"  Erro ao testar {n}x{m}: {e}")

if __name__ == "__main__":
    main()
