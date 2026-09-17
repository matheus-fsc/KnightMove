#ifndef KNIGHT_H
#define KNIGHT_H

#include <cmath>
#include <iostream>
#include <vector>
#include <utility>
#include "board.h"  
#include <fstream>

class Knight {
private:
    int x, y;
    int* movGeometrico;
    int numMoves;        
    // Novo: armazena todos os caminhos encontrados
    std::vector<std::vector<std::pair<int, int>>> caminhos;
    int colisao;
    Tabuleiro* board;

public:
    Knight(int, int, Tabuleiro*);
    ~Knight(void);
    bool validMove(int, int) const;
    void move(int, int);
    void caminhada(void);
    void getPosition(int& outX, int& outY) const { outX = x; outY = y;}
    void gerarMovimentos(void);
    int contarColisoesGeometricas(void);
	int contarColisoesNoCaminho(const std::vector<std::pair<int, int>>& caminho);
    void exportarCaminho(const std::string& filename);
	int contarColisoesComOutros(const std::vector<std::pair<int, int>>& caminho, size_t idxCaminho);
    bool backtrack(int, std::vector<std::pair<int, int>>& caminhoAtual);
    // Novo: acessar todos os caminhos
    const std::vector<std::vector<std::pair<int, int>>>& getCaminhos() const { return caminhos; }
};

#endif