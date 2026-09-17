#!/usr/bin/env python3
import json
import numpy as np
import pandas as pd
import seaborn as sns
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
    H1 = data["descricao"]["H1"]

    # 1. Criar o DataFrame
    colunas = [f"L_{i}" for i in range(H1)]
    df = pd.DataFrame(assinaturas, columns=colunas)

    # 2. Isolar APENAS a "Zona de Comutação Livre" (Os 14 loops azuis)
    # Por que fazemos isso? Porque a correlação de Pearson (df.corr) mede como a
    # variação de A afeta B. Loops 100% ou 0% não têm variação (variância = 0), 
    # o que geraria um erro matemático (divisão por zero / valores NaN).
    frequencia = df.mean() * 100
    loops_livres = frequencia[(frequencia > 0.0) & (frequencia < 100.0)].index.tolist()

    print(f"\nExtraindo a matriz de correlação para {len(loops_livres)} Loops Livres:")
    print(loops_livres)

    # Filtramos o DataFrame para conter apenas esses loops dinâmicos
    df_livres = df[loops_livres]
    if df_livres.shape[1] < 2:
        print("Poucos loops livres para calcular correlação.")
        return

    # 3. Calcular a Matriz de Correlação
    corr_matrix = df_livres.corr()

    # Reordena para agrupar loops mais conectados visualmente.
    ordem = corr_matrix.abs().mean().sort_values(ascending=False).index
    corr_matrix = corr_matrix.loc[ordem, ordem]

    # Esconde metade superior para remover informação duplicada.
    mask = np.triu(np.ones_like(corr_matrix, dtype=bool), k=1)

    # Mostra rótulo numérico apenas para correlações mais relevantes.
    limiar_anotacao = 0.35
    labels = corr_matrix.apply(
        lambda col: col.map(lambda v: f"{v:.2f}" if abs(v) >= limiar_anotacao else "")
    )

    # 4. Configurar e Plotar o Heatmap
    sns.set_theme(style="white", context="talk")
    n = len(corr_matrix)
    fig_w = max(9, n * 0.8)
    fig_h = max(8, n * 0.75)
    plt.figure(figsize=(fig_w, fig_h))

    # A mágica do Seaborn
    # vmin=-1, vmax=1 garante que a escala de cores vá do Vermelho (-1) ao Azul (+1)
    sns.heatmap(
        corr_matrix, 
        mask=mask,
        annot=labels,            # Mostra só os valores fortes para evitar poluição visual
        fmt="",
        cmap="RdBu_r",           # Paleta de cores (Vermelho a Azul)
        center=0,                # Branco fica no 0 (Sem correlação)
        vmin=-1, 
        vmax=1, 
        square=True, 
        linewidths=.35,
        linecolor="#e8e8e8",
        cbar_kws={"shrink": .82, "label": "Correlação de Pearson"},
        annot_kws={"fontsize": 9}
    )

    plt.title("Heatmap de Correlação entre Loops Livres (6x6)", fontsize=16, pad=14)
    plt.xticks(rotation=45, ha="right")
    plt.yticks(rotation=0)
    plt.tight_layout()

    # 5. Salvar a imagem
    file_out = "heatmap_correlacao.png"
    plt.savefig(file_out, dpi=300, bbox_inches='tight')
    print(f"\nHeatmap exportado com sucesso para: {file_out}")
    
    if "agg" not in plt.get_backend().lower():
        plt.show()

if __name__ == "__main__":
    main()
