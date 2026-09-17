# knight_tours_sheared.py
import random
import sys

sys.setrecursionlimit(10000)

KNIGHT_MOVES = [
    (1, 2), (2, 1), (-1, 2), (-2, 1),
    (1, -2), (2, -1), (-1, -2), (-2, -1)
]

def sheared_move(x, y, dx, dy, N, M, s):
    """Retorna (novo_x, novo_y) no toro cisalhado e os wraps em x e y."""
    raw_x = x + dx
    wrap_x = raw_x // N
    new_x = raw_x - wrap_x * N
    raw_y = y + dy + s * wrap_x
    wrap_y = raw_y // M
    new_y = raw_y - wrap_y * M
    return new_x, new_y, wrap_x, wrap_y

def get_neighbors(x, y, N, M, s):
    """Retorna lista de vizinhos (nx, ny) a partir de (x,y)."""
    neighbors = []
    for dx, dy in KNIGHT_MOVES:
        nx, ny, _, _ = sheared_move(x, y, dx, dy, N, M, s)
        neighbors.append((nx, ny))
    return neighbors

def generate_tour(s, N=6, M=6, max_attempts=50):
    """
    Gera um tour fechado do cavalo no toro cisalhado N×M com shear s.
    O tour inicia em (0,0) e retorna a (0,0) após N*M passos.
    A primeira aresta é sorteada uniformemente entre os 8 movimentos.
    """
    start = (0, 0)

    for attempt in range(max_attempts):
        # Escolhe primeiro movimento uniformemente
        first_moves = []
        for dx, dy in KNIGHT_MOVES:
            nx, ny, _, _ = sheared_move(0, 0, dx, dy, N, M, s)
            if (nx, ny) != start:
                first_moves.append(((nx, ny), (dx, dy)))
        random.shuffle(first_moves)

        for first_square, (dx0, dy0) in first_moves:
            path = [start, first_square]
            visited = {start, first_square}

            def dfs(current):
                if len(path) == N * M:
                    # Verifica se o último vértice fecha o ciclo com o início
                    for dx, dy in KNIGHT_MOVES:
                        nx, ny, _, _ = sheared_move(current[0], current[1], dx, dy, N, M, s)
                        if (nx, ny) == start:
                            return True
                    return False

                # Obtém vizinhos não visitados com grau (heurística de Warnsdorff)
                candidates = []
                for nx, ny in get_neighbors(current[0], current[1], N, M, s):
                    if (nx, ny) not in visited:
                        deg = len([v for v in get_neighbors(nx, ny, N, M, s) if v not in visited])
                        candidates.append(((nx, ny), deg))

                # Ordena por grau (Warnsdorff) com desempate aleatório
                random.shuffle(candidates)
                candidates.sort(key=lambda x: x[1])

                for (nx, ny), _ in candidates:
                    path.append((nx, ny))
                    visited.add((nx, ny))
                    if dfs((nx, ny)):
                        return True
                    visited.remove((nx, ny))
                    path.pop()
                return False

            if dfs(first_square):
                return path

    return None
