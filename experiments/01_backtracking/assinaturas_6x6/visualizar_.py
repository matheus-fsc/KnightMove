#!/usr/bin/env python3
import json
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def main():
    # 1. Carregar os Dados extraídos pelo Z3
    print("Carregando assinaturas_6x6.json...")
    try:
        with open("assinaturas_6x6.json", "r") as f:
            data = json.load(f)
    except FileNotFoundError:
        print("Erro: Arquivo assinaturas_6x6.json não encontrado. Rode o extrator primeiro.")
        return

    assinaturas = data["assinaturas_booleanas"]
    H1 = data["descricao"]["H1"]
    num_amostras = len(assinaturas)

    # 2. Criar o DataFrame do Pandas
    # Linhas = Soluções (Amostras), Colunas = Loops L0 a L44
    colunas = [f"L_{i}" for i in range(H1)]
    df = pd.DataFrame(assinaturas, columns=colunas)

    # 3. Calcular a Média de Ativação (Frequência em %)
    # Como True=1 e False=0, a média * 100 dá a exata porcentagem de ativação!
    frequencia = df.mean() * 100

    # 4. Classificar os Loops
    deterministicos = frequencia[frequencia == 100.0]
    mortos = frequencia[frequencia == 0.0]
    livres = frequencia[(frequencia > 0.0) & (frequencia < 100.0)]

    print(f"\n--- RELATÓRIO DA ASSINATURA TOPOLÓGICA ---")
    print(f"Total de Soluções Analisadas: {num_amostras}")
    print(f"Loops Determinísticos (Pilares 100%): {len(deterministicos)}")
    print(f"Loops Mortos (Inviáveis 0%): {len(mortos)}")
    print(f"Loops Livres (Zona de Comutação): {len(livres)}")

    # 5. Configurar o Estilo do Gráfico (Dark Mode)
    plt.style.use('dark_background')
    fig, ax = plt.subplots(figsize=(14, 6))

    # Desenhar as barras
    bars = ax.bar(frequencia.index, frequencia.values, color='gray')

    # Colorir de acordo com o papel estrutural do loop
    for i, bar in enumerate(bars):
        valor = frequencia.iloc[i]
        if valor == 100.0:
            bar.set_color('#5ad06b')  # Verde Neon (Pilares Obrigatórios)
        elif valor == 0.0:
            bar.set_color('#ff5577')  # Vermelho (Nunca ativados)
        else:
            bar.set_color('#66c8ff')  # Azul (Loops Livres - Emaranhado)

    # Estética do Gráfico
    ax.set_title('Assinatura do Caminho Hamiltoniano (Frequência de Ativação dos Loops)', fontsize=16, pad=20)
    ax.set_ylabel('Taxa de Ativação (%)', fontsize=12)
    ax.set_xlabel('IDs dos Loops Fundamentais da Árvore DFS', fontsize=12)

    # Linhas de referência
    ax.axhline(y=100, color='#5ad06b', linestyle='--', alpha=0.3)
    ax.axhline(y=50, color='white', linestyle=':', alpha=0.3)
    ax.axhline(y=0, color='#ff5577', linestyle='--', alpha=0.3)

    plt.xticks(rotation=90, fontsize=9)
    plt.ylim(-5, 110)
    plt.tight_layout()

 # 6. Salvar e Mostrar
    file_out = "assinatura_grafico.png"
    plt.savefig(file_out, dpi=300, bbox_inches='tight') # <-- O parêntese fecha aqui!
    print(f"\nGráfico exportado com sucesso para: {file_out}")
    plt.show()

if __name__ == "__main__":
    main()
