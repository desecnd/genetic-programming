// SimpleGeneticAlgorithm.cpp : This file contains the 'main' function. Program execution begins and ends there.

#include "SGA.h"

int main(int argc, char ** argv) {
	std::ios_base::sync_with_stdio(false);

	unsigned seed = 5489u;
	
	if(argc == 2) {
		seed = std::stoi(argv[1]);	
	}
	
	SGA sga{ seed, 10000, SelectMethod::ROULETTE, ScalingType::LINEAR, [&](double x) { return x * x; } };
	sga.evolution();

	return 0;
}