/* Game of Trust 
*
* Game is played by 2 players, in each turn both of them 
* have to decide whether they cooperate or deceive second player
* there are 4 paybacks, each players get score according to payback table
* its says what score each player gets on every one of 4 different situation
* (C, C), (C, D), (D, C), (D, D)
* Players should earn maximum possible score over x turns
*
* Strategy may vary depending on Payback Table and n of turns
*
* Genetic ALgorithm Player simulates game among population and chooses
* best strategy (that depends on 3 last turns) after x generations
*/

#include <iostream>
#include "GAP.h"

int main(int argc, char ** argv) {
	std::ios::sync_with_stdio(false);
	
	uint32_t seed = 20;
	if ( argc > 1 ) {
		seed = std::stoi(argv[1]);
	}
	
	gap::GAP GAPlayer(seed);
	GAPlayer.evolve(50);
}
