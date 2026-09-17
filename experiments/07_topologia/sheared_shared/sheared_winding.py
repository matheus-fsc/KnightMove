# sheared_winding.py
def get_winding_parity(closed_tour, N, M, s):
    """
    Calcula a paridade (mod 2) dos números de enrolamento (wx, wy)
    para um tour fechado no toro cisalhado.

    closed_tour : lista de coordenadas (x,y) começando e terminando em (0,0)
    N, M        : dimensões do tabuleiro
    s           : parâmetro de cisalhamento

    Retorna (wx_mod2, wy_mod2)
    """
    total_wrap_x = 0
    total_wrap_y = 0
    for i in range(len(closed_tour) - 1):
        x1, y1 = closed_tour[i]
        x2, y2 = closed_tour[i+1]
        dx = x2 - x1
        dy = y2 - y1
        raw_x = x1 + dx
        wrap_x = raw_x // N
        new_x = raw_x - wrap_x * N   # não usado diretamente, só para validação
        raw_y = y1 + dy + s * wrap_x
        wrap_y = raw_y // M
        total_wrap_x += wrap_x
        total_wrap_y += wrap_y
    return total_wrap_x % 2, total_wrap_y % 2
