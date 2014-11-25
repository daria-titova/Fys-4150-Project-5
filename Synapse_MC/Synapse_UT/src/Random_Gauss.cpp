#include <Random_Gauss.h>
#include <iostream>
#include <math.h>
#include <armadillo>

#include <chrono>
#include <random>

using namespace arma;
using namespace std;

rand_gauss::rand_gauss()
{}

double rand_gauss::rand(){
    // normal (Gaussian) distribution

      // construct a trivial random generator engine from a time-based seed:
      unsigned seed = chrono::system_clock::now().time_since_epoch().count();

      //unsigned seed = -1;
      default_random_engine generator(seed);

      normal_distribution<double> distribution(0.0, 1.0/sqrt(2.0));

      double random_number=distribution(generator);

      return random_number;
      //return 0;
    }

