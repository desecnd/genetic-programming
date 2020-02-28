#include "pathfinder.h"

using namespace geo;

Pathfinder::Pathfinder(unsigned int seed = 44) 
    : rng{seed}, xDistr{MIN_X, MAX_X}, yDistr{MIN_Y, MAX_Y}, fraction{0, 1} {}

std::vector<Point> Pathfinder::findBestPath(Circle& queen, Point destination, std::vector<Circle>& sites,  int nOfGenerations) {
    robot = queen;
    dest = destination;
    obstacles = sites;
    nOfGenerations += nOfGen();

    if ( nOfGen() == 0 ) {
        populations.emplace_back(popSize);
        randomize(populations.back());
        evaluate(populations.back());
        calcStats(populations.back());
    }

    while ( nOfGen()  < nOfGenerations ) {
        populations.emplace_back(popSize);
        inherit( populations[nOfGen() - 1], populations[nOfGen() - 2]); 
        evaluate(populations.back());
        calcStats(populations.back());
    }

    Individual best { populations.back().getBest() };

    print(populations.back());

    std::vector<Point> ans;
    for (auto p : best.chrom) 
        ans.push_back(p.first);

    return ans;
}

void Pathfinder::print(Population& pop) {
    std::cerr << "Population of size: " << pop.size << "\n\n";
    for (int i = 0; i < pop.size; i++) {
        std::cerr << i + 1 << ". - ";
        pop.individuals[i].debug();
        std::cerr << "\n";
    }
    std::cerr << "-----------------------------------\n\n";
}

Pathfinder::fitness_t Pathfinder::calcGoodCost(chrom_t& chrom) {
    return wdi * distance(chrom) + wsm * smooth(chrom) + wcl * clear(chrom);
}

Pathfinder::fitness_t Pathfinder::calcBadCost(chrom_t& chrom, fitness_t maxCost) {
    int intersections = 0;
    for (int i = 0; i < chrom.size(); i++) 
        if ( !chrom[i].second ) intersections += 2;

    return double(intersections) + 2.0 + maxCost;
}

bool Pathfinder::markWrong(chrom_t& chrom) {
    bool allGood = true;
    chrom.back().second = true;

    for (int i = 0; i < chrom.size() - 1; i++) {
        chrom[i].second = true;

        Point a{chrom[i].first}, b{chrom[i + 1].first};

        for ( Circle& c : obstacles ) {
            Point p { c.o };
            long double border = c.r + robot.r; 

            if ( segPoint(a, b, p) < border ) {
                chrom[i].second = false;
                allGood = false;
            }
        }
    } 

    return allGood;
}

double Pathfinder::distance(chrom_t& chrom) {
    double sum = 0;
    for (int i = 0; i < chrom.size() - 1; i++) 
        sum += abs(chrom[i].first - chrom[i + 1].first);
    return sum;
}

double Pathfinder::smooth(chrom_t& chrom) {
   double maxS { 0.0 };
   for (int i = 1; i < chrom.size() - 1; i++) {
       Point a{ chrom[i].first }; 
       Point b{ chrom[i - 1].first };
       Point c{ chrom[i + 1].first };

       double s { angle(b - a, c - a) / std::min(abs(b - a), abs(c - a)) };
       maxS = std::max(maxS, s);
   } 
   return maxS;
}

double Pathfinder::clear(chrom_t& chrom) {
    double maxC { 0.0 };
    for (int i = 0; i < chrom.size() - 1; i++) {

        Point a{chrom[i].first}, b{chrom[i + 1].first};

        for ( Circle& c : obstacles ) {
            Point p { c.o };
            double cl { segPoint(a, b, p) - c.r };

            if ( cl < 0 ) cl *= -clearParam;

            maxC = std::max(maxC, cl);
        }
    } 

    return maxC;
}

void Pathfinder::randomize(Population &pop) {
    for (Individual& ind : pop.individuals ) {
        int nodes = std::uniform_int_distribution<int>(2, getMaxChromLen())(rng);
        nodes = std::max(0, nodes - 2);

        ind.chrom.emplace_back(robot.o, true);
        for (int i = 0; i < nodes; i++) ind.chrom.emplace_back(getRandomPoint(), true);
        ind.chrom.emplace_back(dest, true);
    }
}

Pathfinder::Population::Population(size_t n) : size(n), 
    individuals(size), prefixSum(size), sum{0.0}, avg{0.0}, max{0.0}, min{0.0} {
}

void Pathfinder::calcStats(Population& pop) {
    pop.sum = 0.0;
    pop.min = pop.max = pop.individuals.front().fitness;
    int index { 0 };

    for (Individual& ind : pop.individuals) {
       fitness_t curr { ind.fitness };
       pop.sum += curr;

       pop.prefixSum[index++] = pop.sum; 
       pop.min = std::min(pop.min, curr);
       pop.max = std::max(pop.max, curr);
    }
    pop.avg = pop.sum / pop.size;
}

const Pathfinder::Individual& Pathfinder::select(Population& pop) {
    /// Wheel Selection
    fitness_t choice{ fraction(rng) * pop.sum };
    auto it { std::upper_bound(pop.prefixSum.begin(), pop.prefixSum.end(), choice) };
    return pop.individuals[ std::distance(pop.prefixSum.begin(), it) ];
}

const Pathfinder::Individual& Pathfinder::Population::getBest() { 
    return *max_element(individuals.begin(), individuals.end()); 
}

void Pathfinder::evaluate(Population& pop) {

    fitness_t maxCost = 0.0;
    for ( Individual &ind : pop.individuals ) {
        remove(ind.chrom); insert(ind.chrom);
        swap(ind.chrom); 
        smallMutate(ind.chrom);
        largeMutate(ind.chrom);

        ind.valid = markWrong(ind.chrom);
        if ( ind.valid ) ind.cost = calcGoodCost(ind.chrom);

        maxCost = std::max(ind.cost, maxCost);
    } 

    for ( Individual &ind : pop.individuals ) {
        if ( !ind.valid ) ind.cost = calcBadCost(ind.chrom, maxCost);
        ind.fitness = std::max(costBorder - ind.cost, 0.0);
    }

}

void Pathfinder::inherit(Population& curr, Population& last) {
    for (int i = 0; i < curr.size; i += 2) {
        curr.individuals[i] = select(last);
        curr.individuals[i+1] = select(last);
        cross(curr.individuals[i].chrom, curr.individuals[i + 1].chrom);
    }
}


void Pathfinder::cross(chrom_t& chrom1, chrom_t& chrom2) {
    if ( !crossRoll(rng) ) return;

    int c1 { -1 }, c2 { -1 };
    for (int i = 0; i < chrom1.size() - 1; i++) 
        if ( !chrom1[i].second ) c1 = i;

    for (int i = 0; i < chrom2.size() - 1; i++) 
        if ( !chrom2[i].second ) c2 = i;

    if ( c1 == -1 ) c1 = std::uniform_int_distribution<int>(0, chrom1.size() - 2)(rng);
    if ( c2 == -1 ) c2 = std::uniform_int_distribution<int>(0, chrom2.size() - 2)(rng);

    chrom_t newChrom1, newChrom2;
    for (int i = 0; i <= c1; i++) newChrom2.push_back(chrom1[i]);
    for (int i = 0; i <= c2; i++) newChrom1.push_back(chrom2[i]);
    for (int i = c1 + 1; i < chrom1.size(); i++) newChrom2.push_back(chrom1[i]);
    for (int i = c2 + 1; i < chrom2.size(); i++) newChrom1.push_back(chrom2[i]);

    chrom1 = newChrom1;
    chrom2 = newChrom2;
}

void Pathfinder::swap(chrom_t& chrom) {
    if ( !swapRoll(rng) ) return;
    else if ( chrom.size() <= 3 ) return;

    size_t n { chrom.size() };
    size_t p { std::uniform_int_distribution<size_t>(1, n - 2)(rng) };
    chrom_t newChrom;

    newChrom.push_back(chrom.front());
    for (size_t i = p + 1; i < n - 1; i++) newChrom.push_back(chrom[i]); 
    for (size_t i = 1; i <= p; i++) newChrom.push_back(chrom[i]); 
    newChrom.push_back(chrom.back());
    chrom = newChrom;
}

void Pathfinder::insert(chrom_t& chrom) {
    size_t n { chrom.size() };
    for (int i = 1; i < n; i++) {
        if ( n >= Pathfinder::getMaxChromLen() ) return;
        else if ( !insertRoll(rng) ) continue;

        chrom.insert(chrom.begin() + i++, { getRandomPoint(), false});
        n++;
    }
}

void Pathfinder::remove(chrom_t& chrom) {
    size_t n { chrom.size() };
    for (int i = 1; i < n - 1; i++) {
        if ( n <= 2 ) return;
        else if ( !removeRoll(rng) ) continue;

        chrom.erase(chrom.begin() + i--);
        n--;
    }
} 

size_t Pathfinder::smallDelta(size_t z) {
    return std::uniform_int_distribution<size_t>(0, z)(rng);
}

size_t Pathfinder::largeDelta(size_t z) {
    return std::uniform_int_distribution<size_t>(0, z)(rng);
}

void Pathfinder::smallMutate(chrom_t& chrom) {
    size_t n { chrom.size() };
    for (int i = 1; i < n - 1; i++) {
        if ( !smallMutateRoll(rng) ) continue;

        if ( roll(rng) ) chrom[i].first.x -= smallDelta(chrom[i].first.x - MIN_X); 
        else chrom[i].first.x += smallDelta(MAX_X - chrom[i].first.x);

        if ( roll(rng) ) chrom[i].first.y -= smallDelta(chrom[i].first.y - MIN_Y); 
        else chrom[i].first.y += smallDelta(MAX_Y - chrom[i].first.y);
    }
} 

void Pathfinder::largeMutate(chrom_t& chrom) {
    size_t n { chrom.size() };
    for (int i = 1; i < n - 1; i++) {
        if ( !largeMutateRoll(rng) ) continue;

        if ( roll(rng) ) chrom[i].first.x -= largeDelta(chrom[i].first.x - MIN_X); 
        else chrom[i].first.x += largeDelta(MAX_X - chrom[i].first.x);

        if ( roll(rng) ) chrom[i].first.y -= largeDelta(chrom[i].first.y - MIN_Y); 
        else chrom[i].first.y += largeDelta(MAX_Y - chrom[i].first.y);
    }
}

void Pathfinder::Individual::debug() {
    std::cerr << "valid? " << std::boolalpha << valid << " - ";
    for (auto p : chrom ) {
        std::cerr << "{ " << p.first << "," << std::boolalpha << p.second << " }";
    }
    std::cerr << ", cost: " << cost << ", fitness: " << fitness << "\n";
}
