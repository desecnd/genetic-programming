#pragma once

#include <bitset>
#include <vector>
#include <array>
#include <iostream>
#include <random>
#include <iomanip>
#include <memory>
#include <algorithm> 

template<typename T>
const T pi = std::acos(-T(1));

class GeneticAlgorithm {

	static constexpr int chromLen{ 22 };
	static constexpr int popSize{ 50 };
	static constexpr double pCross{ 0.25 };
	static constexpr double pMut{ 0.01 };
	const unsigned int seed{ 20 };
	int generation{ 0 };
	int mutations{ 0 };
	int crossings{ 0 };

	using chromosome_t = std::bitset<chromLen>; 
	using fitness_t = double;

	struct Individual {
		chromosome_t chrom;
		fitness_t fitness;

		Individual();
	};

	struct Population {
		std::array<Individual, popSize> population;
		std::array<fitness_t, popSize> prefixSum;
		fitness_t sum;
		fitness_t avg;
		fitness_t max;
		fitness_t min;

		Population();
		// calls decodePop(), and updates population stats
		auto calcStats() -> void;

		// uses eval function on every individual
		auto decodePop() -> void;

		auto getBestIndividual() -> const Individual&;
	};

	// population.back() means last population
	std::vector< Population > populations;

	std::mt19937 gen;
	std::uniform_real_distribution<double> fractDis;
	std::uniform_int_distribution<size_t> crossPointDis;
	std::bernoulli_distribution mutate{ pMut };
	std::bernoulli_distribution cross{ pCross };

	// converts binary to int value
	static auto decodeChrom(const chromosome_t & chrom) -> int;

	// casts chromosome value to range [-1, 2]
	static auto castChrom(const chromosome_t & chrom)->fitness_t;

	// evaluates casted and decoded chromosome value 
	// based on given fitness function
	static auto evalChrom(const chromosome_t & chrom)->fitness_t;

	// returns reference to chosen Individual in last population
	auto selection(const Population & pop)->const Individual&; 

	// takes population and crosses individuals with each other
	auto crossing(Individual & p1, Individual & p2, size_t point)->void;

	// takes chromosome and applies mutation
	auto mutation(Population & pop)->void;

	// Creates new population based on last population
	// Selection, Crossing, Mutation
	auto evolvePopulation(Population & pop) -> void;

	auto debug() -> void;

public:
	GeneticAlgorithm(unsigned int seed = 20u);

	// Performs whole evolution of GA
	auto evolve(int generations = 150) -> void;

	friend std::ostream& operator<<(std::ostream & stream, const Individual & ind);
	friend std::ostream& operator<<(std::ostream & stream, const Population & pop);
};

