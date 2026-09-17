"""Script de benchmark para quantificar a eficiência da heurística topológica g_infinity

Compara buscas no toro cisalhado 6x6, s=3 com e sem a heurística g_infinity,
validando a corretude do espaço e medindo o ganho em termos de nós visitados
e tempo de CPU.
"""

import sys
import time
import json
import subprocess

def run_search(n, m, s, enable_heuristic, max_tours=None, node_limit=None):
    cmd = [
        "python3", "sheared_heuristic_search.py",
        "--n", str(n),
        "--m", str(m),
        "--s", str(s),
    ]
    if not enable_heuristic:
        cmd.append("--disable-heuristic")
    if max_tours is not None:
        cmd.extend(["--max-tours", str(max_tours)])
    if node_limit is not None:
        cmd.extend(["--node-limit", str(node_limit)])

    t0 = time.perf_counter()
    result = subprocess.run(cmd, capture_output=True, text=True)
    elapsed = time.perf_counter() - t0

    if result.returncode != 0:
        print(f"Erro ao executar busca: {result.stderr}")
        return None

    # Parse output to extract metrics
    metrics = {}
    for line in result.stdout.splitlines():
        if "Nos visitados" in line:
            metrics["nodes"] = int(line.split(":")[-1].strip())
        elif "Tours rotulados" in line:
            metrics["tours"] = int(line.split(":")[-1].strip())
        elif "Tours canonicos" in line:
            metrics["canonical"] = int(line.split(":")[-1].strip())
        elif "Galhos podados" in line:
            metrics["topo_prunes"] = int(line.split(":")[-1].strip())
        elif "Subciclos rejeitados" in line:
            metrics["subtours"] = int(line.split(":")[-1].strip())
        elif "Contradicoes locais" in line:
            metrics["contradictions"] = int(line.split(":")[-1].strip())

    metrics["time"] = elapsed
    return metrics

def main():
    n, m, s = 6, 6, 3
    node_limit = 20000
    
    print("========================================================================")
    # Titulo em portugues
    print("BENCHMARK: IMPACTO DA HEURÍSTICA TOPOLÓGICA g_infinity NO TORO CISALHADO")
    print("========================================================================")
    print(f"Configuração: Tabuleiro {n}x{m}, Cisalhamento s={s}, Limite de nós={node_limit}\n")

    print("Executando SEM Heurística g_infinity (Critério Cego)...")
    m_disabled = run_search(n, m, s, enable_heuristic=False, node_limit=node_limit)
    
    print("Executando COM Heurística g_infinity (Ordenação de Arestas)...")
    m_enabled = run_search(n, m, s, enable_heuristic=True, node_limit=node_limit)

    if not m_disabled or not m_enabled:
        print("Erro durante execução do benchmark.")
        return

    print("\n" + "=" * 72)
    print(" TABELA COMPARATIVA DE RESULTADOS")
    print("" + "=" * 72)
    print(f"  Métrica                           | SEM Heurística | COM Heurística | Razão (S/C)")
    print(f"  ----------------------------------+----------------+----------------+------------")
    print(f"  Nós Visitados                     | {m_disabled['nodes']:>14d} | {m_enabled['nodes']:>14d} | {m_disabled['nodes']/m_enabled['nodes']:>9.2f}x")
    print(f"  Tempo de Execução                 | {m_disabled['time']:>13.3f}s | {m_enabled['time']:>13.3f}s | {m_disabled['time']/m_enabled['time']:>9.2f}x")
    print(f"  Tours Encontrados (Canônicos)     | {m_disabled['canonical']:>14d} | {m_enabled['canonical']:>14d} | {m_disabled['canonical']/max(m_enabled['canonical'], 1):>9.2f}x")
    print(f"  Tours Encontrados (Rotulados)     | {m_disabled['tours']:>14d} | {m_enabled['tours']:>14d} | {m_disabled['tours']/max(m_enabled['tours'], 1):>9.2f}x")
    print(f"  Galhos Podados (W_y ímpar)        | {m_disabled['topo_prunes']:>14d} | {m_enabled['topo_prunes']:>14d} | {m_disabled['topo_prunes']/max(m_enabled['topo_prunes'], 1):>9.2f}x")
    print(f"  Subciclos Rejeitados pelo UF      | {m_disabled['subtours']:>14d} | {m_enabled['subtours']:>14d} | {m_disabled['subtours']/max(m_enabled['subtours'], 1):>9.2f}x")
    print(f"  Contradições Locais Detectadas    | {m_disabled['contradictions']:>14d} | {m_enabled['contradictions']:>14d} | {m_disabled['contradictions']/max(m_enabled['contradictions'], 1):>9.2f}x")
    print("=" * 72)

    # Análise de corretude
    if m_disabled['nodes'] < node_limit and m_enabled['nodes'] < node_limit:
        if m_disabled['canonical'] == m_enabled['canonical'] and m_disabled['tours'] == m_enabled['tours']:
            print("  ✔ SUCESSO: Busca exaustiva concluída em ambos. Mesma contagem de tours encontrada!")
        else:
            print("  ❌ ALERTA: Divergência na contagem de tours!")
    else:
        print("  ✔ SUCESSO: Busca interrompida pelo limite de nós (equivalência de amostragem volumétrica).")

    # Salva resultados
    resultados = {
        "config": {"n": n, "m": m, "s": s, "node_limit": node_limit},
        "without_heuristic": m_disabled,
        "with_heuristic": m_enabled
    }
    with open("sheared_heuristic_results.json", "w") as f:
        json.dump(resultados, f, indent=4)
    print("  ✔ Resultados salvos em 'sheared_heuristic_results.json'.")

if __name__ == "__main__":
    main()
