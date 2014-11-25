#ifndef CRANIC_H
#define CRANIC_H
#include <armadillo>
using namespace arma;

class Crank_Nicolson
{public:

    //constructor
    Crank_Nicolson();
    //functions
    void Crank_Nicolson_Scheme (vec &V, double alpha, int n, int m, double dx);

};

#endif // CRANIC_H
