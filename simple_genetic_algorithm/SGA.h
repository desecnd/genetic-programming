#pragma once
#include <array>
#include <iostream>
#include <random>
#include <utility>
#include <iomanip>
#include <functional>
#include <numeric>
#include <algorithm>

enum class ScalingType {
	NONE, LINEAR,
};

enum class SelectMethod {
	ROULETTE, DETERMINISTIC, TRUNCATION,
};

class SGA {

	static constexpr int popSize		{ 30 };		// even number
	static constexpr int chromosomeLen	{ 10 };
	static constexpr double mutProb		{ 0.03 };
	static constexpr double crossProb	{ 0.6 };

	const unsigned int inputSeed		{ 0u };
	const int maxGenerations			{ 100 };
	
	using allele_t = bool;
	using fitness_t = double;
	using chromosome_t = std::array<allele_t, chromosomeLen>;

	struct Individual {
		chromosome_t genotype;
		fitness_t fitness;
		fitness_t expectedCopies;    //  individual fitness / avg population fitness
		int id;

		void decodeGenotype();	// decodes genotype and updates fitness value
		Individual();
		Individual(chromosome_t chrom);
		~Individual();


		friend std::ostream & operator<<(std::ostream & stream, const SGA::Individual & ind) {
			stream << "ID = " << std::setw(3) << ind.id;
			stream << ", Genotype: ";
			for (auto & al : ind.genotype)
				stream << al;
			stream << ", Fitness = " << ind.fitness;
			return stream;
		}
	};

	struct Population {
		std::array<Individual, popSize> population;

		fitness_t sumFitness;
		fitness_t avgFitness;
		fitness_t maxFitness;
		fitness_t minFitness;

		Population();
		~Population();

		void resetIDs();
		void calculateStatistics();
		void decodeIndividuals(std::function<fitness_t(fitness_t)>& fitFunc);
		void updateIndividuals();
		void scaleLinear();
		void scaleNone();
	};

	SelectMethod selectMethod;
	ScalingType scalingType;
	std::function<fitness_t(fitness_t)> fitnessFunction;

	Population lastPop;
	Population currPop;

	// Counters
	int mutCnt;
	int crossCnt;
	int generations;

	// Random generators
	std::random_device rd;
	std::mt19937_64 gen;
	std::uniform_real_distribution<fitness_t> fractionDis;
	std::uniform_int_distribution<int> crossPointDis;
	std::bernoulli_distribution mutFlip;
	std::bernoulli_distribution crossFlip;

	// Scaling functions, update Individuals "fitness" fields in pop
	
	inline auto getRandomFraction() -> fitness_t { return fractionDis(gen); }

	auto scalePopulation(Population & pop) ->void;
	auto selectPopulation(Population & last, Population & curr) -> void;
	auto deterministicSelection(const Population & last, Population & curr) -> void;
	auto rouletteSelection(const Population & last, Population & curr) -> void;
	auto truncationSelection(const Population & last, Population & curr) -> void; 
	auto crossPopulation(Population & curr) -> void;
	auto calculatePopulation(Population & pop) ->void;
	auto simpleCrossing(const Individual & parent1, const Individual & parent2) -> std::pair<Individual, Individual>;
	auto mutatePopulation(Population & pop) -> void;
	auto evolvePopulation(Population & last, Population & curr) ->void;
	auto debug(std::ostream & file, Population & pop) -> void;

public:
	SGA(unsigned int seed = 0u, int maxGen = 100, SelectMethod selectM = SelectMethod::ROULETTE,
		ScalingType scaleType = ScalingType::NONE,
		std::function<fitness_t(fitness_t)> fitFunc = [&](fitness_t x) { return x; });
	~SGA();

	void evolution() {
		
		std::cout << "Starting with seed: " << inputSeed << '\n';
		std::cout << "MaxGenerations = " << maxGenerations << '\n';
		debug(std::cout, lastPop);

		for(int i = 1; i <= maxGenerations; i++)
			evolvePopulation(lastPop, currPop);
		debug(std::cout, currPop);
	}

};

