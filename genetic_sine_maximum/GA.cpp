#include "GA.h"

GeneticAlgorithm::GeneticAlgorithm(unsigned int seed) :
	seed(seed), populations(1), gen(seed), fractDis(0.0, 1.0), 
	crossPointDis(1, chromLen-1) {

	std::bernoulli_distribution flip{ 0.5 };

	// generate random population
	for (Individual & ind : populations.back().population)
		for (int i = 0; i < chromLen; i++) 
			if (flip(gen)) ind.chrom.set(i);
		
}

auto GeneticAlgorithm::evolve(int generations) -> void {
	populations.back().calcStats();

	for (int i = 1; i <= generations; i++) {
		generation++;
		Population newPop{};
		evolvePopulation(newPop);
		newPop.calcStats();
		populations.push_back(newPop);	
	}

	//debug();
	std::cout << "Evolved: " << generations << " generations.\n";
	std::cout << "Mutations:" << mutations << '\n';
	std::cout << "Crossings:" << crossings << '\n';
	std::cout << "Best Individual:{\n    " << populations.back().getBestIndividual() << "\n}\n";

}

auto GeneticAlgorithm::debug() -> void {
	std::cerr << "Generation: " << populations.size()-1 << "\n";
	std::cerr << "Crossings: " << crossings << '\n';
	std::cerr << "Mutations: " << mutations << '\n';
	std::cerr << populations.back();
}

auto GeneticAlgorithm::decodeChrom(const chromosome_t & chrom) -> int {
	return static_cast<int>(chrom.to_ulong());
}

auto GeneticAlgorithm::castChrom(const chromosome_t & chrom) -> fitness_t {
	fitness_t bin_val{ static_cast<fitness_t>(decodeChrom(chrom)) };
	fitness_t intervalBegin{ -1.0 };
	fitness_t intervalEnd{ 2.0 };
	fitness_t maxVal{ static_cast<fitness_t>(1 << chromLen) };
	fitness_t nIntervals{ intervalEnd - intervalBegin };

	return intervalBegin + (bin_val*nIntervals) / maxVal;
}

auto GeneticAlgorithm::evalChrom(const chromosome_t & chrom) -> fitness_t {
	fitness_t x{ castChrom(chrom) };
	return x * std::sin(10 * pi<fitness_t> * x) + 1.0;
}

auto GeneticAlgorithm::selection(const Population & last) -> const Individual& {
	fitness_t choice{ fractDis(gen) * last.sum };
	auto it = std::upper_bound(last.prefixSum.begin(), last.prefixSum.end(), choice);
	return last.population[std::distance(last.prefixSum.begin(), it)];
}

auto GeneticAlgorithm::crossing(Individual & p1, Individual & p2, size_t point) -> void {
	chromosome_t ch1{ p1.chrom };
	chromosome_t ch2 { p2.chrom};

	// bitset counts from left [n, n-1, ... 0]
	for (int i = chromLen-1; i >= chromLen-point; i--) {
		if (p1.chrom.test(i)) ch2.set(i);
		else ch2.reset(i);

		if (p2.chrom.test(i)) ch1.set(i);
		else ch1.reset(i);
	}

	p1.chrom = ch1;
	p2.chrom = ch2;
}

auto GeneticAlgorithm::mutation(Population & pop)-> void {
	for (Individual & ind : pop.population) 
	for (int i = 0; i < chromLen; i++) {
		if (mutate(gen)) {
			ind.chrom.flip(i);
			mutations++;
		}
	}
}

auto GeneticAlgorithm::evolvePopulation(Population & pop) -> void{
	for (int i = 0; i < popSize; i += 2) {
		pop.population[i] = selection(populations.back());
		pop.population[i+1] = selection(populations.back());

		if (cross(gen)) {
			size_t point{ crossPointDis(gen) };
			crossing(pop.population[i], pop.population[i + 1], point);
			crossings++;
		}
	}
	mutation(pop);
}

GeneticAlgorithm::Population::Population() : population(), prefixSum(),
sum(), avg(), max(), min() {}

auto GeneticAlgorithm::Population::calcStats() -> void{
	decodePop();

	fitness_t cSum{ 0.0 };
	max = min = population.front().fitness;

	int index{ };
	for (const Individual & ind : population) {
		fitness_t curr{ ind.fitness };
		cSum += curr;

		prefixSum[index++] = cSum;
		max = std::max(max, curr);
		min = std::min(min, curr);
	}
	sum = cSum;
	avg = cSum / popSize;
}

auto GeneticAlgorithm::Population::decodePop() -> void {
	for (Individual & ind : population) {
		ind.fitness = GeneticAlgorithm::evalChrom(ind.chrom);
	}
}

auto GeneticAlgorithm::Population::getBestIndividual() -> const Individual& {
	fitness_t currMax{ population[0].fitness };
	int best{ 0 };
	int index{ };
	for(const Individual & ind : population ) {
		if (ind.fitness > currMax) {
			best = index;
			currMax = ind.fitness;
		}
			
		index++;
	}
	return population[best];
}

auto operator<<(std::ostream & stream, const GeneticAlgorithm::Population & pop) -> std::ostream&{
	stream << "Fitness Stats:   ";
	stream << "sum: " << std::setw(7) << pop.sum;
	stream << "    avg: " << std::setw(7) << pop.avg;
	stream << "    max: " << std::setw(7) << pop.max;
	stream << "    min: " << pop.min << "\n";
	stream << "---------|\n";
	for (int i = 0; i < pop.population.size(); i++) {
		stream << "Ind#" << std::setw(3) << i << " -- " << pop.population[i] << '\n';
		//stream << "PrefSum: " << pop.prefixSum[i] << '\n';
	}
	return stream;
}

GeneticAlgorithm::Individual::Individual() : chrom(), fitness() {}

auto operator<<(std::ostream & stream, const GeneticAlgorithm::Individual & ind) ->std::ostream& {
	stream << "genotype:" << ind.chrom;
	stream << "    decoded:" << std::setw(7) << GeneticAlgorithm::decodeChrom(ind.chrom);
	stream << "    casted(-1,2):" << std::setw(7) << GeneticAlgorithm::castChrom(ind.chrom);
	stream << "    fitness:" << std::setw(7) << ind.fitness;;
	return stream;
}
