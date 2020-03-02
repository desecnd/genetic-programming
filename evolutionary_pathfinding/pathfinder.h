/*
* author: pavveu
* structure for evolutionary pathfinding
* enviroment is rectangle of constant size 
* every obstacle (including robot) is a circle
*/

#ifndef PATHFINDER_213888_H
#define PATHFINDER_213888_H

#include "geometry.h"

#include <SFML/Window.hpp>
#include <SFML/Graphics.hpp>

#include <iostream>
#include <iomanip>
#include <sstream>
#include <chrono>
#include <random>


class Pathfinder { 
private:

// DEBUG SFML
sf::RenderWindow window;
sf::View view;
sf::Font font;


// USED TYPES
using chrom_t = std::vector<std::pair<geo::Point, bool>>;
using fitness_t = double;

// MEMBERS 
bool local = false;
const size_t popSize = 100;

double clearParam { 1 };
double wdi { 2 };
double wsm { 2000 };
double wcl { 2 };
fitness_t costBorder { 10000.0 };

const int MAX_X = 1920; 
const int MAX_Y = 1000;
const int MIN_X = 0; 
const int MIN_Y = 0;

std::mt19937 rng;
std::uniform_int_distribution<int> xDistr;
std::uniform_int_distribution<int> yDistr;
std::uniform_real_distribution<double> fraction;

std::bernoulli_distribution crossRoll { 0.3 };
std::bernoulli_distribution swapRoll { 0.01 };
std::bernoulli_distribution insertRoll { 0.05 };
std::bernoulli_distribution removeRoll { 0.05 };
std::bernoulli_distribution smallMutateRoll { 0.05 };
std::bernoulli_distribution largeMutateRoll { 0.05 };
std::bernoulli_distribution smootheRoll { 0.1 };
std::bernoulli_distribution roll { 0.5 };

std::vector<geo::Circle> obstacles;
geo::Circle robot;
geo::Point dest;

/// TODO -------------------------------------------------
// solve geometric problems connected to double precision
/// ------------------------------------------------------
    struct Individual {
        /// Individual is a simple structure containing one chromosome and fitness values
        chrom_t chrom; bool valid;
        fitness_t cost, fitness;
        Individual() : chrom{}, valid{ true }, cost{}, fitness{} {}
        bool operator<(const Individual& ind) { return fitness < ind.fitness; }
        void debug();
    };

    struct Population {
        size_t size;
        std::vector<Individual> individuals;
        std::vector<fitness_t> prefixSum;
		fitness_t sum, avg, max,min;

		Population(size_t n);
        const Individual& getBest();
    };

    std::vector<Population> populations;

    // Inline functions
    size_t nOfGen() { return populations.size(); }
    size_t getMaxChromLen() { return ( local? nOfGen() : obstacles.size() ); }
public:
    geo::Point getRandomPoint() { return geo::Point(xDistr(rng), yDistr(rng)); }
private:

    // COST METHODS
    fitness_t calcGoodCost(chrom_t& chrom);
    fitness_t calcBadCost(chrom_t& chrom, fitness_t maxCost);

    bool markWrong(chrom_t& chrom);
    double distance(chrom_t& chrom);
    double smooth(chrom_t& chrom);
    double clear(chrom_t& chrom);

    // OPERATORS HELPER FUNCTIONS
    size_t smallDelta(size_t z);
    size_t largeDelta(size_t z);

    // CHROMOSOME OPERATORS
    void cross(chrom_t& chrom1, chrom_t& chrom2);
    void swap(chrom_t& chrom);
    void insert(chrom_t& chrom); 
    void remove(chrom_t& chrom);
    // TODO - choose better distributions, higher gen -> chances of 0 increases
    void smallMutate(chrom_t& chrom);
    void largeMutate(chrom_t& chrom);
    // --------------

    // POPULATION OPERATORS
    const Individual& select(Population& pop);
	void inherit(Population& curr, Population& last);
	void evaluate(Population& pop);
    void randomize(Population& pop);
    void calcStats(Population& pop);
    void print(Population& pop);
    void draw(Population pop, int x, bool pauseAfter);

    Pathfinder();
    Pathfinder(const Pathfinder& );

public:
    std::vector<geo::Point> findBestPath(geo::Circle& queen, geo::Point destination, std::vector<geo::Circle>& sites, int nOfGenerations);

    void test(geo::Circle& queen, std::vector<geo::Circle>& sites);

    static Pathfinder& getPathfinder() {
        static Pathfinder pathfinder;
        return pathfinder;
    }
};

#endif