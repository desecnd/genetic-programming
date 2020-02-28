/*
* Evolutionary Pathfinding 
*/

#include <iostream>
#include "evolution.h"
#include "geometry.h"
using namespace std;

int main() {
    ios::sync_with_stdio(0);
    int nSite; cin >> nSite;
    cerr << nSite << "\n";
    
    vector<Circle> sites(nSite);
    for (Circle& s : sites) 
        cin >> s.o.x >> s.o.y >> s.r;

    Circle queen({219,219}, 30.0);
    Point dest { 810, 460 };


    vector<Point> path { gen::evolveBestPath(queen, dest, sites, 150) };

    for (Point p : path) cout << p << " ";
    cout << endl;

    return 0;
}
