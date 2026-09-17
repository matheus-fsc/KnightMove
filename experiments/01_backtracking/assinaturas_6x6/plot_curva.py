#!/usr/bin/env python3
import json
import pandas as pd
import matplotlib.pyplot as plt


def main():
    print("Carregando assinaturas_6x6.json...")
    try:
        with open("assinaturas_6x6.json", "r") as f:
            data = json.load(f)
    except FileNotFoundError:
        print("Erro: Arquivo assinaturas_6x6.json não encontrado.")
        return

    assinaturas = data["assinaturas_booleanas"]
    h1 = data["descricao"]["H1"]

    colunas = [f"L_{i}" for i in range(h1)]
    df = pd.DataFrame(assinaturas, columns=colunas)

    # Quantidade de True por solução (peso da assinatura).
    df["qtd_true"] = df.sum(axis=1)

    # Curva: para cada peso, quantas soluções existem.
    curva = df["qtd_true"].value_counts().sort_index()
    print("\nDistribuição (qtd_true -> número de soluções):")
    print(curva.to_string())

    plt.style.use("seaborn-v0_8-whitegrid")
    plt.figure(figsize=(9, 5))
    plt.plot(curva.index, curva.values, marker="o", linewidth=2)
    plt.xlabel("Quantidade de True na assinatura")
    plt.ylabel("Número de soluções")
    plt.title("Curva de distribuição dos pesos das assinaturas (6x6)")
    plt.tight_layout()

    file_out = "curva_pesos_assinaturas.png"
    plt.savefig(file_out, dpi=300, bbox_inches="tight")
    print(f"\nCurva exportada com sucesso para: {file_out}")

    if "agg" not in plt.get_backend().lower():
        plt.show()


if __name__ == "__main__":
    main()
