#include <cstdlib>
#include <iostream>
// Define printf()
#if 1  // Choose between Mersenne Twister and SFMT generator
  #include "randomc.h"
  #define  STOC_BASE CRandomMersenne     // define random number generator base class
#else
  #include "sfmt.h"
  #define  STOC_BASE CRandomSFMT
#endif

int main() {

	int seed_a = 0; //time
	int seed_b; //threadnr
	
	// Array of seeds
	const int seeds[] = 
		{
			seed_a,
			seed_b
		};
	STOC_BASE RanObject(0);

	// Generate random seeds
	RanObject.RandomInitByArray(seeds, 2); // Initialize RNG
	
	// Get a random number
	std::cout << RanObject.Random() << std::endl;
	std::cout << RanObject.Random() << std::endl;
	std::cout << RanObject.Random() << std::endl;
	std::cout << RanObject.Random() << std::endl;
	std::cout << RanObject.IRandomX(0,1) << std::endl;
	std::cout << RanObject.IRandomX(0,1) << std::endl;
	std::cout << RanObject.IRandomX(0,1) << std::endl;

	

	return 0;
}

