/* Genetic Player Class for Prisoners Dilemma Trust Game
* 
* Genetic Algorithm finds optimal strategy based on last 3 games 
* GA Player makes first 3 moves, then depending on other player moves
* he chooses his optimal 
* 
*/

#pragma once
#include <utility>
#include <array>
#include <vector>
#include <random>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <string>
#include <sstream>
#include <bitset>

namespace gap {
class GAP {
	static constexpr size_t gameRounds{ 150 };
	static constexpr size_t chromLen{ 70 };
	static constexpr size_t premoveIndex{ 64 };
	static constexpr size_t popSize{ 30 };
	static constexpr double pMut{ 0.01 };
	static constexpr double pCross{ 0.25 };

	const unsigned int seed{ 20u };
	int mutations{ 0 };
	int crossings{ 0 };
	
	using move_t = bool;

	// rounds go from left to right =>  (3 ago) -> (2ago) -> (last) 
	// first move is always mine in players perspective
	using round_t = std::pair<move_t, move_t>;
	using fitness_t = double; 

	static const move_t cooperate{ true };
	static const move_t deceive{ false };

	using chromosome_t = std::array<move_t, chromLen>;

	struct GeneticPlayer {
		chromosome_t strategy;
		fitness_t fitness;

		GeneticPlayer();

		// Sets player pre moves (6 last bits in chromosome) 
		auto getPreMoves(std::vector<round_t>& preMoves) -> void;
	};

	struct Population {
		std::array<GeneticPlayer, popSize> pop;
		fitness_t sum, avg;
		fitness_t max, min;

		auto calcStats() -> void;
		Population();
	};

	std::vector<Population> populations;

	std::mt19937 gen;
	std::uniform_real_distribution<fitness_t> fractDis;
	std::uniform_int_distribution<size_t> crossPointDis{ 1, chromLen - 1 };
	std::bernoulli_distribution mutate{ pMut };
	std::bernoulli_distribution cross{ pCross };

	auto currPop() -> Population& { return populations.back(); }

	inline static auto moveSymbol(move_t move) -> char {
		if (move == cooperate) return 'C';
		else return 'D';
	}


	// get index, based on 6-bit value, return round format (me,him)(me,him)(me,him)
	auto getRoundFormat(unsigned int index)->std::string;

	// returns genetic player strategy based on 3 last moves
	// order of moves is first, second third, so third is last move in game
	auto getStrategyIndex(const round_t & first,const round_t & second, const round_t & thrid)->size_t;

	// returns both players fitnesses
	auto payoff(const round_t & round)-> std::pair<fitness_t, fitness_t>;

	// plays one game for every player, updates fitness sum
	auto playGame(GeneticPlayer & p1, GeneticPlayer & p2) -> void;

	// play games among every player and sets their fitness scores
	auto tournament(Population & pop) -> void;

	auto selection(Population & next) -> void;

	auto crossing(Population & next) -> void;

	auto mutation(Population & next) -> void;

	auto exportPlayer(const GeneticPlayer& gp) -> void;

	auto getBestPlayer(Population & pop) -> const GeneticPlayer&;

	auto debug() -> void;

public:
	//auto evolveStrategy() -> void;

	auto evolve(int generations = 50) -> void;

	GAP(unsigned int seed = 20u);

	friend std::ostream& operator<<(std::ostream & stream, const GeneticPlayer& gp) {
		stream << "Strategy: ";
		for (const move_t & move : gp.strategy)
			stream << moveSymbol(move);
		stream << "    fitness: " << gp.fitness;
		return stream;
	}
	friend std::ostream& operator<<(std::ostream & stream, const Population &pop) {
		stream << "PopulationStats: "
			<< "    sum: " << pop.sum
			<< "    avg: " << pop.avg
			<< "    min: " << pop.min
			<< "    max: " << pop.max << '\n';
		int index{ };
		for (const GeneticPlayer & gp : pop.pop) {
			stream << std::setw(3) << index++ << "--> ";
			stream << gp << '\n';
		}
		return stream;
	}
};

};