/*
* author: pavveu
* structure for evolutionary pathfinding
* enviroment is rectangle of constant size 
* every obstacle (including robot) is a circle
*/

#pragma once 

#include "geometry.h"
#include <iostream>
#include <random>

namespace gen { 

using chrom_t = vector<pair<Point, bool>>;
using fitness_t = double;

bool local = false;
int generation = 0;
size_t popSize = 100;

double clearParam { 0.5 };
double wdi { 0.5 };
double wsm { 1.5 };
double wcl { 2.0 };
fitness_t costBorder { 1000 };

mt19937 rng(35);
const int MAX_X = 1920; 
const int MAX_Y = 1000;
const int MIN_X = 0; 
const int MIN_Y = 0;
uniform_int_distribution<int> xDistr(MIN_X, MAX_X);
uniform_int_distribution<int> yDistr(MIN_Y, MAX_Y);
uniform_real_distribution<double> fraction(0, 1);

bernoulli_distribution crossRoll { 1 };
bernoulli_distribution swapRoll { 0.1 };
bernoulli_distribution insertRoll { 0.1 };
bernoulli_distribution removeRoll { 0.1 };
bernoulli_distribution smallMutateRoll { 0.1 };
bernoulli_distribution largeMutateRoll { 0.1 };
bernoulli_distribution smootheRoll { 0.1 };
bernoulli_distribution roll { 0.5 };

vector<Circle> obstacles;
Circle robot;
Point dest;

/// TODO -------------------------------------------------
// solve geometric problems connected to double precision
/// ------------------------------------------------------

    vector<Point> evolveBestPath(Circle& queen, Point destination, vector<Circle>& sites, int generations);

    fitness_t calcGoodCost(chrom_t& chrom);
    fitness_t calcBadCost(chrom_t& chrom, fitness_t maxCost);

    bool markWrong(chrom_t& chrom);
    double distance(chrom_t& chrom);
    double smooth(chrom_t& chrom);
    double clear(chrom_t& chrom);

    struct Individual {
        chrom_t chrom; bool valid;
        fitness_t cost, fitness;
        Individual() : chrom{}, valid{ true }, cost{}, fitness{} {}
        bool operator<(Individual& ind) { return fitness < ind.fitness; }
    };

    struct Population {
        size_t size;
        vector<Individual> population;
        vector<fitness_t> prefixSum;
		fitness_t sum, avg, max,min;

        // TODO random generation of population
		Population(int n, bool random);
        void calcStats();
        const Individual& select();
        const Individual& getBest();
		void evolvePopulation(Population& last);
    };

    vector<Population> populations;

    // --- OPERATORS
    void cross(chrom_t& chrom1, chrom_t& chrom2);
    void swap(chrom_t& chrom);
    void insert(chrom_t& chrom); 
    void remove(chrom_t& chrom);
    // TODO - choose better distributions, higher gen -> chances of 0 increases
    size_t smallDelta(size_t z);
    size_t largeDelta(size_t z);
    // ODOT 
    void smallMutate(chrom_t& chrom);
    void largeMutate(chrom_t& chrom);
    // --------------

    Point inline getRandomPoint() { return Point(xDistr(rng), yDistr(rng)); }
    size_t inline getMaxChromLen() { return ( local? generation : obstacles.size() ); }


    void test();
};

