#include "knight.h"
#include <set>
#include <tuple>

bool segmentosColidem(int,int,int,int,int,int,int,int);

Knight::Knight(int startX, int startY, Tabuleiro* t) {
    x = startX;
    y = startY;
    this->board = t;
    movGeometrico = nullptr;
    numMoves = 0;

    for (int i = 0; i < board->getLinhas(); i++)
        for (int j = 0; j < board->getColunas(); j++)
            board->at(i, j) = -1;

    board->at(y, x) = 0; // marca pos inicial
    gerarMovimentos(); 
}

void Knight::gerarMovimentos() {
    delete[] movGeometrico;
    movGeometrico = new int[16]; // 8 movimentos * 2 inteiros (dx, dy)

    numMoves = 0;
    for (int dx = -2; dx <= 2; dx++) {
        for (int dy = -2; dy <= 2; dy++) {
            if (dx == 0 && dy == 0) continue;
            if (dx * dx + dy * dy == 5) {
                movGeometrico[numMoves * 2] = dx;
                movGeometrico[numMoves * 2 + 1] = dy;
                numMoves++;
            }
        }
    }
}


bool Knight::validMove(int newX, int newY) const {
    return (newX >= 0 && newX < board->getColunas() &&
        newY >= 0 && newY < board->getLinhas() &&
        board->at(newY, newX) == -1);
}

void Knight::move(int newX, int newY) {
    if (validMove(newX, newY)) {
        x = newX;
        y = newY;
        std::cout << "Movimento para (" << x << ", " << y << ")\n";
    }
    else {
        std::cout << "Movimento inválido.\n";
    }
}

// Novo backtrack: armazena todos os caminhos
bool Knight::backtrack(int step, std::vector<std::pair<int, int>>& caminhoAtual) {
    if (step == board->getLinhas() * board->getColunas()) {
        caminhos.push_back(caminhoAtual);
        return false; // Não para na primeira solução
    }

    for (int i = 0; i < numMoves; i++) {
        int dx = movGeometrico[i * 2];
        int dy = movGeometrico[i * 2 + 1];
        int nx = x + dx;
        int ny = y + dy;

        if (validMove(nx, ny)) {
            board->at(ny, nx) = step;
            int oldX = x, oldY = y;
            x = nx; y = ny;
            caminhoAtual.push_back({nx, ny});

            backtrack(step + 1, caminhoAtual);

            x = oldX; y = oldY;
            board->at(ny, nx) = -1;
            caminhoAtual.pop_back();
        }
    }
    return false;
}

void Knight::caminhada() {
    caminhos.clear();
    std::vector<std::pair<int, int>> caminhoAtual;
    caminhoAtual.push_back({x, y});
    backtrack(1, caminhoAtual);
    if (!caminhos.empty()) {
        std::cout << "Total de caminhos encontrados: " << caminhos.size() << "\n";
        // Exemplo: imprime o primeiro caminho
        std::cout << "Primeiro caminho encontrado:" << std::endl;
        for (const auto& pos : caminhos[0]) {
            std::cout << "(" << pos.first << ", " << pos.second << ") ";
        }
        std::cout << std::endl;

        // Exibe o número de colisões geométricas
        int colisoes = contarColisoesGeometricas();
        std::cout << "Colisões geométricas entre caminhos: " << colisoes << std::endl;
    } else {
        std::cout << "Não existe passeio do cavalo para esse tabuleiro.\n";
    }
}

// Função auxiliar: conta colisões entre um caminho e todos os outros
int Knight::contarColisoesComOutros(const std::vector<std::pair<int, int>>& caminho, size_t idxCaminho) {
    int colisoes = 0;
    std::set<std::tuple<int,int,int,int,int,int,int,int>> colisoesRegistradas;
    for (size_t a = 1; a < caminho.size(); ++a) {
        int x1 = caminho[a-1].first, y1 = caminho[a-1].second;
        int x2 = caminho[a].first,   y2 = caminho[a].second;
        for (size_t j = 0; j < caminhos.size(); ++j) {
            if (j == idxCaminho) continue;
            const auto& outro = caminhos[j];
            for (size_t b = 1; b < outro.size(); ++b) {
                int x3 = outro[b-1].first, y3 = outro[b-1].second;
                int x4 = outro[b].first,   y4 = outro[b].second;
                if (segmentosColidem(x1, y1, x2, y2, x3, y3, x4, y4)) {
                    auto key = std::make_tuple(x1,y1,x2,y2,x3,y3,x4,y4);
                    if (colisoesRegistradas.insert(key).second) {
                        ++colisoes;
                    }
                }
            }
        }
    }
    return colisoes;
}

void Knight::exportarCaminho(const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Erro ao abrir arquivo para escrita.\n";
        return;
    }

    int idx = 1;
    for (size_t i = 0; i < caminhos.size(); ++i) {
        const auto& caminho = caminhos[i];
        file << "Caminho " << idx++ << ":\n";
        for (const auto& pos : caminho) {
            file << pos.first << "," << pos.second << "\n";
        }
        int colisoes = contarColisoesComOutros(caminho, i);
        file << "Colisoes: " << colisoes << "\n\n";
    }
    file.close();
    std::cout << "Caminhos exportados para: " << filename << "\n";
}

int Knight::contarColisoesGeometricas() {
    int colisoes = 0;
    std::set<std::tuple<int,int,int,int,int,int,int,int>> colisoesRegistradas;

    for (size_t i = 0; i < caminhos.size(); ++i) {
        const auto& caminhoA = caminhos[i];
        for (size_t j = i + 1; j < caminhos.size(); ++j) {
            const auto& caminhoB = caminhos[j];
            // Para cada segmento de A
            for (size_t a = 1; a < caminhoA.size(); ++a) {
                int x1 = caminhoA[a-1].first, y1 = caminhoA[a-1].second;
                int x2 = caminhoA[a].first,   y2 = caminhoA[a].second;
                // Para cada segmento de B
                for (size_t b = 1; b < caminhoB.size(); ++b) {
                    int x3 = caminhoB[b-1].first, y3 = caminhoB[b-1].second;
                    int x4 = caminhoB[b].first,   y4 = caminhoB[b].second;
                    if (segmentosColidem(x1, y1, x2, y2, x3, y3, x4, y4)) {
                        // Evita contar colisões duplicadas
                        auto key = std::make_tuple(x1,y1,x2,y2,x3,y3,x4,y4);
                        if (colisoesRegistradas.insert(key).second) {
                            ++colisoes;
                        }
                    }
                }
            }
        }
    }
    return colisoes;
}

// Função utilitária para checar interseção de segmentos
bool segmentosColidem(
    int x1, int y1, int x2, int y2,
    int x3, int y3, int x4, int y4)
{
    auto ccw = [](int ax, int ay, int bx, int by, int cx, int cy) {
        return (cy - ay) * (bx - ax) > (by - ay) * (cx - ax);
        };
    return (ccw(x1, y1, x3, y3, x4, y4) != ccw(x2, y2, x3, y3, x4, y4)) &&
        (ccw(x1, y1, x2, y2, x3, y3) != ccw(x1, y1, x2, y2, x4, y4));
}

Knight::~Knight() {
    delete[] movGeometrico;
}

int Knight::contarColisoesNoCaminho(const std::vector<std::pair<int, int>>& caminho) {
    int colisoes = 0;
    std::set<std::tuple<int,int,int,int>> colisoesRegistradas;
    for (size_t i = 1; i < caminho.size(); ++i) {
        int x1 = caminho[i-1].first, y1 = caminho[i-1].second;
        int x2 = caminho[i].first,   y2 = caminho[i].second;
        for (size_t j = i+1; j < caminho.size(); ++j) {
            int x3 = caminho[j-1].first, y3 = caminho[j-1].second;
            int x4 = caminho[j].first,   y4 = caminho[j].second;
            if (segmentosColidem(x1, y1, x2, y2, x3, y3, x4, y4)) {
                // Evita contar colisões duplicadas
                auto key = std::make_tuple(x1,y1,x2,y2);
                if (colisoesRegistradas.insert(key).second) {
                    ++colisoes;
                }
            }
        }
    }
    return colisoes;
}