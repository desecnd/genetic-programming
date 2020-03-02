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
    cerr << nSite << std::endl;
    
    vector<Circle> sites(nSite);
    for ( int i = 0; i < nSite; i++) {
        int index; cin >> index;
        cin >> sites[index].o.x >> sites[index].o.y >> sites[index].r;
    }

    Circle queen({417,750}, 30.0);
    Pathfinder::getPathfinder().test(queen, sites);


    return 0;
}
