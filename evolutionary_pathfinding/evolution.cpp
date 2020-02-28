#include "evolution.h"
using namespace gen;

vector<Point> evolveBestPath(Circle& queen, Point destination, vector<Circle>& sites,  int generations) {
    robot = queen;
    dest = destination;
    obstacles = sites;

    populations.emplace_back(popSize, true);
    populations.back().calcStats();

    for (int g = 2; g <= generations; g++) {
        generation++;
        populations.emplace_back(popSize, false);
        populations.back().evolvePopulation( populations[populations.size() - 2 ] ); 
        populations.back().calcStats();
    }

    Individual best { populations.back().getBest() };

    vector<Point> ans;
    for (auto p : best.chrom) 
        ans.push_back(p.first);

    return ans;
}

fitness_t gen::calcGoodCost(chrom_t& chrom) {
    return wdi * distance(chrom) + wsm * smooth(chrom) + wcl * clear(chrom);
}

fitness_t gen::calcBadCost(chrom_t& chrom, fitness_t maxCost) {
    int intersections = 0;
    for (int i = 0; i < chrom.size(); i++) 
        if ( !chrom[i].second ) intersections += 2;

    return double(intersections) + 2.0 + maxCost;
}

bool gen::markWrong(chrom_t& chrom) {
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

double gen::distance(chrom_t& chrom) {
    double sum = 0;
    for (int i = 0; i < chrom.size() - 1; i++) 
        sum += abs(chrom[i].first - chrom[i + 1].first);
    return sum;
}

double gen::smooth(chrom_t& chrom) {
   double maxS { 0.0 };
   for (int i = 1; i < chrom.size() - 1; i++) {
       Point a{ chrom[i].first }; 
       Point b{ chrom[i - 1].first };
       Point c{ chrom[i + 1].first };

       double s { angle(b - a, c - a) / min(abs(b - a), abs(c - a)) };
       maxS = max(maxS, s);
   } 
   return maxS;
}

double gen::clear(chrom_t& chrom) {
    double maxC { 0.0 };
    for (int i = 0; i < chrom.size() - 1; i++) {

        Point a{chrom[i].first}, b{chrom[i + 1].first};

        for ( Circle& c : obstacles ) {
            Point p { c.o };
            double cl { segPoint(a, b, p) - c.r };

            if ( cl < 0 ) cl *= -clearParam;

            maxC = max(maxC, cl);
        }
    } 

    return maxC;
}



gen::Population::Population(int n, bool random = false) : size(n), 
    population(size), prefixSum(size), sum{}, avg{}, max{}, min{} {

    if ( random ) {
        for (Individual& ind : population ) {
            int nodes = uniform_int_distribution<int>(2, getMaxChromLen())(rng);
            nodes = std::max(0, nodes - 2);

            ind.chrom.emplace_back(robot.o, true);
            for (int i = 0; i < nodes; i++) ind.chrom.emplace_back(getRandomPoint(), true);
            ind.chrom.emplace_back(dest, true);
        }
    }
}

void gen::Population::calcStats() {
    sum = 0.0;
    min = max = population.front().fitness;
    int index { 0 };

    for (Individual& ind : population) {
       fitness_t curr { ind.fitness };
       sum += curr;

       prefixSum[index++] = sum; 
       min = std::min(min, curr);
       max = std::max(max, curr);
    }
    avg = sum / size;
}

const Individual& gen::Population::select() {
    /// Wheel Selection
    fitness_t choice{ fraction(rng) * sum };
    auto it { std::upper_bound(prefixSum.begin(), prefixSum.end(), choice) };
    return population[ std::distance(prefixSum.begin(), it) ];
}

const Individual& gen::Population::getBest() { 
    return *max_element(population.begin(), population.end()); 
}

void gen::Population::evolvePopulation(Population& last) {
    for (int i = 0; i < size; i += 2) {
        population[i] = last.select();
        population[i+1] = last.select();
        cross(population[i].chrom, population[i + 1].chrom);
    }

    fitness_t maxCost = 0.0;
    for ( Individual &ind : population ) {
        remove(ind.chrom); insert(ind.chrom);
        swap(ind.chrom); 
        smallMutate(ind.chrom);
        largeMutate(ind.chrom);

        ind.valid = markWrong(ind.chrom);
        if ( ind.valid ) ind.cost = calcGoodCost(ind.chrom);

        maxCost = std::max(ind.cost, maxCost);
    } 

    for ( Individual &ind : population ) {
        if ( !ind.valid ) ind.cost = calcBadCost(ind.chrom, maxCost);
        ind.fitness = std::max(costBorder - ind.cost, 0.0);
    }
}


void gen::cross(chrom_t& chrom1, chrom_t& chrom2) {
    if ( !crossRoll(rng) ) return;

    int c1 { -1 }, c2 { -1 };
    for (int i = 0; i < chrom1.size() - 1; i++) 
        if ( !chrom1[i].second ) c1 = i;

    for (int i = 0; i < chrom2.size() - 1; i++) 
        if ( !chrom2[i].second ) c2 = i;

    if ( c1 == -1 ) c1 = uniform_int_distribution<int>(0, chrom1.size() - 2)(rng);
    if ( c2 == -1 ) c2 = uniform_int_distribution<int>(0, chrom2.size() - 2)(rng);

    chrom_t newChrom1, newChrom2;
    for (int i = 0; i <= c1; i++) newChrom2.push_back(chrom1[i]);
    for (int i = 0; i <= c2; i++) newChrom1.push_back(chrom2[i]);
    for (int i = c1 + 1; i < chrom1.size(); i++) newChrom2.push_back(chrom1[i]);
    for (int i = c2 + 1; i < chrom2.size(); i++) newChrom1.push_back(chrom2[i]);

    chrom1 = newChrom1;
    chrom2 = newChrom2;
}

void gen::swap(chrom_t& chrom) {
    if ( !swapRoll(rng) ) return;
    else if ( chrom.size() <= 3 ) return;

    size_t n { chrom.size() };
    size_t p { uniform_int_distribution<size_t>(1, n - 2)(rng) };
    chrom_t newChrom;

    newChrom.push_back(chrom.front());
    for (size_t i = p + 1; i < n - 1; i++) newChrom.push_back(chrom[i]); 
    for (size_t i = 1; i <= p; i++) newChrom.push_back(chrom[i]); 
    newChrom.push_back(chrom.back());
    chrom = newChrom;
}

void gen::insert(chrom_t& chrom) {
    size_t n { chrom.size() };
    for (int i = 1; i < n; i++) {
        if ( n >= gen::getMaxChromLen() ) return;
        else if ( !insertRoll(rng) ) continue;

        chrom.insert(chrom.begin() + i++, { getRandomPoint(), false});
        n++;
    }
}

void gen::remove(chrom_t& chrom) {
    size_t n { chrom.size() };
    for (int i = 1; i < n - 1; i++) {
        if ( n <= 2 ) return;
        else if ( !removeRoll(rng) ) continue;

        chrom.erase(chrom.begin() + i--);
        n--;
    }
} 

size_t gen::smallDelta(size_t z) {
    return uniform_int_distribution<size_t>(0, z)(rng);
}

size_t gen::largeDelta(size_t z) {
    return uniform_int_distribution<size_t>(0, z)(rng);
}

void gen::smallMutate(chrom_t& chrom) {
    size_t n { chrom.size() };
    for (int i = 1; i < n - 1; i++) {
        if ( !smallMutateRoll(rng) ) continue;

        if ( roll(rng) ) chrom[i].first.x += smallDelta(chrom[i].first.x - MIN_X); 
        else chrom[i].first.x -= smallDelta(MAX_X - chrom[i].first.x);

        if ( roll(rng) ) chrom[i].first.y += smallDelta(chrom[i].first.y - MIN_Y); 
        else chrom[i].first.y -= smallDelta(MAX_Y - chrom[i].first.y);
    }
} 

void gen::largeMutate(chrom_t& chrom) {
    size_t n { chrom.size() };
    for (int i = 1; i < n - 1; i++) {
        if ( !largeMutateRoll(rng) ) continue;

        if ( roll(rng) ) chrom[i].first.x += largeDelta(chrom[i].first.x - MIN_X); 
        else chrom[i].first.x -= largeDelta(MAX_X - chrom[i].first.x);

        if ( roll(rng) ) chrom[i].first.y += largeDelta(chrom[i].first.y - MIN_Y); 
        else chrom[i].first.y -= largeDelta(MAX_Y - chrom[i].first.y);
    }
}