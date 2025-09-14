#include "board.h"

Tabuleiro::Tabuleiro(int &linhas, int &colunas) {
    this->x= linhas;
    this->y= colunas;

    this->mapa = new int* [linhas];

    for (int i = 0; i < linhas; i++) {
        this->mapa[i] = new int[colunas];
        for (int j = 0; j < colunas; j++) {
            this->mapa[i][j] = 0;
        }
    }
}

void Tabuleiro::print() const {
    for (int i = 0; i < this->x; i++) {
        for (int j = 0; j < this->y; j++) {
            std::cout << this->mapa[i][j] << "\t";
        }
        std::cout << std::endl;
    }
}

Tabuleiro::~Tabuleiro() {
    for (int i = 0; i < x; i++) {
        delete[] this->mapa[i];
    }
    delete[] this->mapa;
}


int& Tabuleiro::at(int linha, int coluna) {
   return this->mapa[linha][coluna];
}