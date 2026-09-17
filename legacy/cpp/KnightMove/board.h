#ifndef BOARD_H
#define BOARD_H

#include <iostream>

class Tabuleiro {
private:
    int x, y;
    int** mapa;
public:
    int& at(int, int);
    void print() const;
	int getLinhas() const { return x; }
	int getColunas() const { return y; }
    Tabuleiro(int&, int&);
	~Tabuleiro();
};

#endif