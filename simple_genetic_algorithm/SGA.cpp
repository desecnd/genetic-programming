#include "SGA.h"

SGA::SGA(unsigned int seed, int maxGen, SelectMethod selectM, ScalingType scaleType, std::function<fitness_t(fitness_t)> fitFunc) :
	inputSeed(seed), maxGenerations(maxGen), selectMethod(selectM), scalingType(scaleType), 
	fitnessFunction(fitFunc), lastPop(), currPop(), mutCnt(0), 
	crossCnt(0), generations(), rd(), gen(( seed? seed : rd() )),
	fractionDis(0, std::nextafter(1, std::numeric_limits<fitness_t>::max())),
	crossPointDis(1, chromosomeLen-1), mutFlip(mutProb), crossFlip(crossProb)
{
	// Generating random population
	// Setting it as "lastPop" 

	std::bernoulli_distribution bitFlip(0.5);
	int id{ };
	for (Individual & ind : lastPop.population) {
		ind.id = id++;
		for (allele_t & al : ind.genotype)
			al = bitFlip(gen);
	}

	calculatePopulation(lastPop);
}

SGA::~SGA() { }

void SGA::debug(std::ostream & file,  Population & pop) {

	file << "Generation #" << generations << '\n';
	file << "Stats: {\n";
	file << "SumF = " << pop.sumFitness << '\n';
	file << "AvgF = " << pop.avgFitness << '\n';
	file << "MaxF = " << pop.maxFitness << '\n';
	file << "MinF = " << pop.minFitness << '\n';
	file << "mutN = " << mutCnt << '\n';
	file << "croN = " << crossCnt << '\n';
	file << "} Population: \n";

	
	for (const Individual & ind : pop.population) {
		file << ind << '\n';
	}

	file << "|-----------------------\n\n";
}

void SGA::scalePopulation(Population & pop) {

	switch (this->scalingType) {
	case ScalingType::NONE: pop.scaleNone(); break;
	case ScalingType::LINEAR: pop.scaleLinear(); break;
	default: std::cerr << "Scaling Error! Unknown Scaling Type\n";
	}

	pop.calculateStatistics();
	pop.updateIndividuals();
	pop.resetIDs();
}

void SGA::evolvePopulation(Population & last, Population & curr) {
	generations++;

	// Gets last population, and based on choosen selection method 
	// fills current population genotypes 
	selectPopulation(last, curr);	

	// Chooses crossing, based on crossing type
	// Cross adjecent individuals in current population
	crossPopulation(curr);				

	// Iterates over alleles of individuals in population, applying mutations
	mutatePopulation(curr);

	// decodes new population (calc individual fitness, based on fitness function)
	// calculates stats (sum, avg, min, max)
	// resets population IDs
	// calculates other individuals stats (expected Copies)
	calculatePopulation(curr);	


	// dont scale if last generation
	if (generations != maxGenerations) 
	{
		// scale population based on scalingType
		scalePopulation(curr);
	}

	//debug(std::cerr, curr);			// debug
	last = curr;						// change generations
}

void SGA::crossPopulation(Population & curr) {
	for (int i = 0; i < popSize; i += 2) {

		std::pair<Individual, Individual> childsPair{ 
			simpleCrossing(curr.population[i], curr.population[i+1]) 
		};

		curr.population[i] = childsPair.first;
		curr.population[i+1] = childsPair.second;
	}
}

void SGA::calculatePopulation(Population & pop) {
	pop.decodeIndividuals(fitnessFunction);
	pop.calculateStatistics();
	pop.updateIndividuals();
	pop.resetIDs();
}

void SGA::selectPopulation(Population & last, Population & curr) {

	switch (this->selectMethod) {
	case SelectMethod::ROULETTE: rouletteSelection(last, curr); break;
	case SelectMethod::DETERMINISTIC: deterministicSelection(last, curr); break;
	case SelectMethod::TRUNCATION: truncationSelection(last, curr); break;
	default: std::cerr << "Selecting Error! Unknown selecting type\n";
	}

}

void SGA::rouletteSelection(const Population & last, Population & curr) {
	for (Individual & ind : curr.population) {

		fitness_t partialSum{ };
		fitness_t randSum{ SGA::getRandomFraction() * last.sumFitness };

		auto it{ last.population.begin() };

		while (partialSum <= randSum) {
			partialSum += it->fitness;
			it++;
		}
		ind = *(--it);
	}
}

void SGA::deterministicSelection(const Population & last, Population & curr) {

	std::array<fitness_t, popSize> fractionPart{ };
	std::array<int, popSize> selectPerm{ };

	std::iota(selectPerm.begin(), selectPerm.end(), 0);

	int nextID{ };
	int spaceLeft{ popSize };

	for (const Individual & ind : last.population) {
		int sureCopies{ static_cast<int>( ind.expectedCopies) };

		fractionPart[ind.id] = ind.expectedCopies - static_cast<fitness_t>(sureCopies);

		// Copy sure amount of individuals into next gen
		while (sureCopies--) {
			curr.population[nextID++] = ind;
			spaceLeft--;
		}
	}

	// sort Individuals (ids) accordingly to their expected copies fractions 
	std::sort(selectPerm.begin(), selectPerm.end(), [&](int a, int b) {
		return fractionPart[a] > fractionPart[b];
	});

	// Based on sorted order, copy individuals
	for (int i = 0; spaceLeft && i < popSize; i++) {
		curr.population[nextID++] = last.population[selectPerm[i]];
		spaceLeft--;
	}
}

void SGA::truncationSelection(const Population & last, Population & curr) {

	std::array<fitness_t, popSize> fractionPart{ };
	std::array<int, popSize> selectPerm{ };

	std::iota(selectPerm.begin(), selectPerm.end(), 0);

	int nextID{ };
	int spaceLeft{ popSize };

	for (const Individual & ind : last.population) {
		int sureCopies{ static_cast<int>(ind.expectedCopies) };

		fractionPart[ind.id] = ind.expectedCopies - static_cast<fitness_t>(sureCopies);

		// Copy sure amount of individuals into next gen
		while (sureCopies--) {
			curr.population[nextID++] = ind;
			spaceLeft--;
		}
	}

	// sort Individuals (ids) accordingly to their expected copies fractions 
	std::sort(selectPerm.begin(), selectPerm.end(), [&](int a, int b) {
		return fractionPart[a] > fractionPart[b];
	});

	// Based on sorted order, try bernoully's distribution
	// untill you fill whole population
	for (int i = 0; spaceLeft; i++) {
		int t{ selectPerm[i%popSize] };  // sorted order;

		std::bernoulli_distribution copyFlip(fractionPart[t]);
		if (copyFlip(gen)) {
			curr.population[nextID++] = last.population[t];
			spaceLeft--;
		}
	}
}

std::pair<SGA::Individual, SGA::Individual> 
	SGA::simpleCrossing(const Individual & parent1, const Individual & parent2) {

	// if not flipped crossing, just return selected parents
	if (!crossFlip(gen))	
		return std::make_pair(parent1, parent2);
	
	// select crossing points
	int cp{ crossPointDis(gen) };

	// Increase counter
	crossCnt++;

	chromosome_t ch1{ parent1.genotype };
	chromosome_t ch2{ parent2.genotype };

	// simple genotype crossing
	std::swap_ranges(ch1.begin(), ch1.begin() + cp, ch2.begin());
	
	// create individuals based on their genotype
	Individual child1{ ch1 };
	Individual child2{ ch2 };

	// return crossed pair of children
	return std::make_pair(child1, child2);
}

void SGA::mutatePopulation(Population & pop) {

	for (Individual & ind : pop.population) {
		for (allele_t & all : ind.genotype) {
			if (this->mutFlip(gen)) {
				all = !all;
				this->mutCnt++;
			}
		}
	}
}

SGA::Population::Population() : 
	population(), sumFitness(), avgFitness(), maxFitness(), minFitness() {
}

SGA::Population::~Population() {}

void SGA::Population::scaleLinear(){
	fitness_t multiplier{ 2.0 };
	fitness_t delta{ };
	fitness_t a{ 1.0 };
	fitness_t b{ 0.0 };

	if (minFitness != maxFitness && multiplier != 1.0) {
		if (minFitness > (multiplier * avgFitness - maxFitness) / (multiplier - 1.0)) {
			delta = maxFitness - avgFitness;
			a = (multiplier - 1.0) * avgFitness / delta;
			b = avgFitness * (maxFitness - multiplier * avgFitness) / delta;
		}
		else {
			delta = avgFitness - minFitness;
			a = avgFitness / delta;
			b = -minFitness * avgFitness / delta;
		}
	}

	for (Individual & ind : population) {
		ind.fitness = fabs(a * ind.fitness + b);
	}
}

void SGA::Population::scaleNone(){

}

void SGA::Population::resetIDs() {
	int IDcnt{ };
	for (Individual & ind : this->population) {
		ind.id = IDcnt++;
	}
}

void SGA::Population::decodeIndividuals(std::function<fitness_t(fitness_t)>& fitFunc) {
	for (Individual & ind : this->population) {
		ind.decodeGenotype();
		ind.fitness = fitFunc(ind.fitness);
	}
}

void SGA::Population::updateIndividuals() {
	for (Individual & ind : this->population) {
		ind.expectedCopies = ind.fitness / this->avgFitness;
	}
}

void SGA::Population::calculateStatistics() {
	sumFitness = avgFitness = 0;
	maxFitness = std::numeric_limits<fitness_t>::min();
	minFitness = std::numeric_limits<fitness_t>::max();

	for (const Individual & ind : this->population) {
		sumFitness += ind.fitness;
		maxFitness = std::max(maxFitness, ind.fitness);
		minFitness = std::min(minFitness, ind.fitness);
	}

	avgFitness = sumFitness / popSize;
}

SGA::Individual::Individual() { }

SGA::Individual::Individual(chromosome_t chrom) : 
	genotype(chrom), fitness( 0.0 ), expectedCopies(0), id(-1) {
}

SGA::Individual::~Individual() {}

void SGA::Individual::decodeGenotype() {
	fitness_t x{ };

	int bitCnt{ chromosomeLen - 1 };
	for (allele_t & al : this->genotype) {
		if (al) x += std::pow(2.0, bitCnt);
		bitCnt--;
	}

	this->fitness = x;
}


