#include "knight.h"
#include "board.h"

int main(int argc, char* argv[]) {
    int x = 0, y = 0;
	int knightX = 0, knightY = 0;

    if (argc == 3) {
		x = std::atoi(argv[1]);
		y = std::atoi(argv[2]);
	}
	else if (argc == 5) {
		x = std::atoi(argv[1]);
		y = std::atoi(argv[2]);
		knightX = std::atoi(argv[3]);
		knightY = std::atoi(argv[4]);
	}
	else {
		std::cout << "Iniciando com valores padrão (5x5 tabuleiro, cavalo em 0,0)\n";
		x = 5;
		y = 5;
	}
	
    Tabuleiro t(x, y);
    Knight k(knightX,knightY, &t);

    k.caminhada();
	k.exportarCaminho("caminho.txt");
    return 0;
}
