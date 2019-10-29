#include "GA.h"
#include <iostream>

int main(int argc, char ** argv) {
	
	uint32_t seed = 20u;
	
	if ( argc > 1 ) { seed = std::stoi(argv[1]); }
	
	std::cout << "Genetic Algorithm for finding maximum of function: "
		<< "f(x) = x sin( 10PIx) + 1,0 for all x c [-1, 2] \n";
		
	std::cout << "Starting with seed=" << seed << '\n';

	GeneticAlgorithm GA{ seed };
	GA.evolve();
}

