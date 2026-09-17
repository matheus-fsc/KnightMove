import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os

# Dados empíricos
L = np.array([0, 1, 2, 3, 4, 5])
f_obs = np.array([0.528, 0.193, 0.198, 0.294, 0.247, 0.261])

# O limite assintótico teórico (toro infinito)
ASYMPTOTIC_LIMIT = 0.25

# Função teórica: f(L) = 0.25 + A * exp(-lambda * L)
def model_func(L, A, lambda_):
    return ASYMPTOTIC_LIMIT + A * np.exp(-lambda_ * L)

def fit_and_evaluate(L_data, f_data, label):
    # Estimativa inicial: A tenta ser a diferença do primeiro ponto pro limite, lambda razoável
    p0 = [f_data[0] - ASYMPTOTIC_LIMIT, 1.0]
    
    popt, _ = curve_fit(model_func, L_data, f_data, p0=p0)
    A, lambda_ = popt
    
    # Métricas de avaliação
    f_pred = model_func(L_data, A, lambda_)
    mse = np.mean((f_data - f_pred)**2)
    ss_res = np.sum((f_data - f_pred)**2)
    ss_tot = np.sum((f_data - np.mean(f_data))**2)
    r2 = 1 - (ss_res / ss_tot)
    
    print(f"--- {label} ---")
    print(f"A      = {A:.4f}")
    print(f"lambda = {lambda_:.4f}")
    print(f"MSE    = {mse:.6f}")
    print(f"R²     = {r2:.4f}\n")
    
    return A, lambda_

if __name__ == '__main__':
    print("Iniciando a Suíte de Análise da Heurística de Fase Local...\n")

    # 1. Cenário A: Usando todos os dados
    A_all, l_all = fit_and_evaluate(L, f_obs, "Cenário A (Com Borda: L=0 a 5)")

    # 2. Cenário B: Excluindo a borda
    mask = L >= 1
    A_bulk, l_bulk = fit_and_evaluate(L[mask], f_obs[mask], "Cenário B (Apenas Bulk: L>=1)")

    # 3. Plotagem para visualização da aderência
    L_plot = np.linspace(0, 6, 100)
    plt.figure(figsize=(9, 6))

    # Plot dos pontos reais
    plt.scatter(L, f_obs, color='black', zorder=5, label='Dados Empíricos', s=60)
    plt.scatter(L[0], f_obs[0], color='red', zorder=6, label='Borda (L=0)', s=80, marker='s')

    # Plot das curvas ajustadas
    plt.plot(L_plot, model_func(L_plot, A_all, l_all), 'r--', alpha=0.7, 
             label=f'Ajuste A (Com Borda)\n$0.25 + {A_all:.3f} e^{{-{l_all:.2f}L}}$')
    plt.plot(L_plot, model_func(L_plot, A_bulk, l_bulk), 'b-', linewidth=2, 
             label=f'Ajuste B (Apenas Bulk)\n$0.25 {A_bulk:+.3f} e^{{-{l_bulk:.2f}L}}$')

    # Assíntota teórica
    plt.axhline(ASYMPTOTIC_LIMIT, color='gray', linestyle=':', label='Limite Toro Infinito (0.25)')

    plt.title("Ajuste da Heurística de Fase Local $f_\\infty(L)$", fontsize=14)
    plt.xlabel("Nível de Interioridade (L)", fontsize=12)
    plt.ylabel("Frequência de Aresta $f(L)$", fontsize=12)
    plt.legend(fontsize=11)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()

    # Salvar e mostrar o plot
    output_path = os.path.join(os.getcwd(), "fit_heuristica_fase.png")
    plt.savefig(output_path, dpi=300)
    print(f"Gráfico salvo com sucesso em: {output_path}")
