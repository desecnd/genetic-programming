/*
* Evolutionary Pathfinding 
*/
#include "pathfinder.h"
#include "geometry.h"

#include <iostream>
using namespace std;
using namespace geo;

int main() {
    ios::sync_with_stdio(0);
    int nSite; cin >> nSite;
    cerr << nSite << "\n";
    
    vector<Circle> sites(nSite);
    for (Circle& s : sites) 
        cin >> s.o.x >> s.o.y >> s.r;

    Circle queen({266,266}, 30.0);
    Point dest { 790, 135 };


    vector<Point> path { Pathfinder::getPathfinder().findBestPath(queen, dest, sites, 20) };

    for (Point p : path) cout << p << " ";
    cout << endl;

    return 0;
}
