#ifndef EXPLICIT_H
#define EXPLICIT_H
#include <armadillo>
using namespace arma;


class Explicit
{public:

    //constructor
    Explicit();

    //functions
    void Explicit_Scheme (vec &V, double alpha, int n, int m, double dx);

};

#endif // EXPLICIT_H
