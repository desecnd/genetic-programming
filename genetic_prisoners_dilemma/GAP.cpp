#include "GAP.h"

using namespace gap;

GAP::GAP(unsigned int seed) : seed(seed), populations(1),
	gen(seed) {
	std::bernoulli_distribution flip{ 0.5 };
	for (GeneticPlayer & gp : populations.back().pop) {
		for (move_t & move : gp.strategy)
			if (flip(gen)) move = cooperate;
	}
}

auto GAP::payoff(const round_t & round) -> std::pair<fitness_t, fitness_t> {
	std::pair<fitness_t, fitness_t> payoff;
	if (round.first == cooperate) {
		if (round.second == cooperate)
			payoff = std::make_pair(3, 3);
		else
			payoff = std::make_pair(0, 5);
	} else {
		if (round.second == deceive)
			payoff = std::make_pair(1, 1);
		else
			payoff = std::make_pair(5, 0);
	}
	return payoff;
}

auto GAP::getStrategyIndex(const round_t & threeAgo, const round_t & twoAgo, const round_t & last)->size_t {
	size_t index{ 0 };

	// index is 6 bit number, 2^6 => [0, 63] value
	// binary works like arrays, 
	//	> 2 least significant bits -> last
	//  > 2 most significant bits -> 3 ago rounds


	// example:			3 ago						2 ago					 last
	// (11)(00)(00) -> (I Cooperate, he Cooperates) ( I Deceive, he Deceives) ( I Deceive, he Deceives)
	index += (static_cast<size_t>(threeAgo.first)	<< 5);
	index += (static_cast<size_t>(threeAgo.second)	<< 4);
	index += (static_cast<size_t>(twoAgo.first)		<< 3);
	index += (static_cast<size_t>(twoAgo.second)	<< 2);
	index += (static_cast<size_t>(last.first)		<< 1);
	index += (static_cast<size_t>(last.second)		<< 0);

	return index;
}

auto GAP::playGame(GeneticPlayer & p1, GeneticPlayer & p2) -> void {
	std::vector<round_t> p1rounds;
	std::vector<round_t> p2rounds;

	fitness_t p1fitness{ 0.0 };
	fitness_t p2fitness{ 0.0 };

	p1rounds.reserve(gameRounds + 10);
	p2rounds.reserve(gameRounds + 10);

	p1.getPreMoves(p1rounds);
	p2.getPreMoves(p2rounds);

	for (size_t i = 0; i < gameRounds; i++) {
		move_t p1move{ p1.strategy[getStrategyIndex(
			p1rounds[i], p1rounds[i + 1], p1rounds[i + 2])] };
		move_t p2move{ p2.strategy[getStrategyIndex(
			p2rounds[i], p2rounds[i + 1], p2rounds[i + 2])] };
	
		// for every player, first move is his own move
		p1rounds.push_back(std::make_pair(p1move, p2move));
		p2rounds.push_back(std::make_pair(p2move, p1move));

		std::pair<fitness_t, fitness_t> payoffs{ payoff(std::make_pair(p1move, p2move)) };

		p1fitness += payoffs.first;
		p2fitness += payoffs.second;
	}

	p1.fitness += p1fitness;
	p2.fitness += p2fitness;
}

auto GAP::tournament(Population & pop) -> void {
	for (size_t i = 0; i < popSize; i++) {
		for (size_t j = i + 1; j < popSize; j++) {
			playGame(pop.pop[i], pop.pop[j]);
		}
	}

	// fitness is avg player's match score
	for (GeneticPlayer & gp : pop.pop) {
		gp.fitness /= popSize - 1;
	}
}

auto GAP::evolve(int generations) -> void {
	tournament(currPop());
	currPop().calcStats();

	for (int i = 1; i <= generations; i++) {
		std::cerr << "generation: #" << i << '\n';
		Population next;

		selection(next);
		crossing(next);
		mutation(next);

		tournament(next);
		next.calcStats();

		populations.push_back(next);
	}

	debug();
	exportPlayer(getBestPlayer(currPop()));
}

auto GAP::selection(Population & next) -> void {
	std::array<fitness_t, popSize> prefSum{ };
	fitness_t currSum{0.0};
	size_t index{ 0 };

	for (const GeneticPlayer & gp : populations.back().pop) {
		currSum += gp.fitness;
		prefSum[index++] = currSum;
	}

	for (GeneticPlayer & gp : next.pop) {
		fitness_t choice{ fractDis(gen) * populations.back().sum };
		auto it{ std::lower_bound(prefSum.begin(), prefSum.end(), choice) };
		uint32_t x{ std::distance(prefSum.begin(), it) };
		
		gp = populations.back().pop[x];
	}
}

auto GAP::crossing(Population & next) -> void {
	for (int i = 0; i < popSize; i += 2) {
		if (!cross(gen)) continue;
		
		crossings++;
		size_t crossPoint{ crossPointDis(gen) };

		chromosome_t & chrom1{ next.pop[i].strategy };
		chromosome_t & chrom2{ next.pop[i + 1].strategy };

		std::swap_ranges(chrom1.begin(), chrom1.begin() + crossPoint,
			chrom2.begin());
	}
}

auto GAP::mutation(Population & next) -> void {
	for (GeneticPlayer & gp : next.pop)
	for (move_t & move : gp.strategy) {
		if (!mutate(gen)) continue;

		mutations++;
		if (move == cooperate) move = deceive;
		else move = cooperate;
	}
}

auto GAP::debug() -> void {
	std::cerr << currPop() << "\n";
}

auto GAP::getRoundFormat(unsigned int index) -> std::string {
	std::stringstream ss;

	std::string x{ std::bitset<6>(index).to_string() };
	std::array<move_t, 6> seq;

	// transition from reverse binary to indexed array 
	for (size_t i = 0; i < 6; i++) 
		seq[i] = static_cast<move_t>(x[i] - '0');

	for (size_t i = 0; i < 6; i += 2)
		ss << "(" << moveSymbol(seq[i]) << "," << moveSymbol(seq[i + 1]) << ")";

	return ss.str();
}

auto GAP::exportPlayer(const GeneticPlayer& gp) -> void {
	std::cerr << "Premoves: ";
	for (size_t i = premoveIndex; i < chromLen; i+=2) 
		std::cerr << "(" << moveSymbol(gp.strategy[i]) << "," << moveSymbol(gp.strategy[i + 1]) << ")";

	std::cerr << '\n';

	for (size_t i = 0; i < 64; i++) {
		std::cerr << "[" << i << "]";
		std::cerr << getRoundFormat(i) << " ";
		std::cerr << "  ----> " << moveSymbol(gp.strategy[i]) << '\n';
	}
}

auto GAP::getBestPlayer(Population & pop) -> const GeneticPlayer& {
	fitness_t currMax{ std::numeric_limits<fitness_t>::min() };
	size_t index{ 0 };
	for (size_t i = 0; i < popSize; i++) {
		if (pop.pop[i].fitness > currMax) {
			index = i;
			currMax = pop.pop[i].fitness;
		}
	}
	return pop.pop[index];
}

GAP::GeneticPlayer::GeneticPlayer() : strategy(), fitness() {} 

auto GAP::GeneticPlayer::getPreMoves(std::vector<round_t>& preMoves) -> void {
	// index:		0		1		2		3		4		5	
	//			(my3ago, him3ago)(my2ago, him2ago)(mylast, himlast)
	// example: 101100 -> (I Cooperate, he Deceives)(I Cooperate, he Cooperates)(I Deceive, He Deceives)
	preMoves.push_back(std::make_pair(strategy[64], strategy[65]));
	preMoves.push_back(std::make_pair(strategy[66], strategy[67]));
	preMoves.push_back(std::make_pair(strategy[68], strategy[69]));
}

GAP::Population::Population() : pop(), sum(), avg(), max(), min() {}

auto GAP::Population::calcStats() -> void {
	sum = 0;
	max = min = pop[0].fitness;
	for (const GeneticPlayer & gp : pop) {
		sum += gp.fitness;
		max = std::max(max, gp.fitness);
		min = std::min(min, gp.fitness);
	}
	avg = sum / popSize;
}


