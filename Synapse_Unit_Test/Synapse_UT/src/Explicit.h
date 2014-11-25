#ifndef EXPLICIT_H
#define EXPLICIT_H
#include <armadillo>
using namespace arma;


class Explicit
{public:

    //constructor
    Explicit();

    //functions
    void Explicit_Scheme (mat &U, double alpha, int n, int m, double dx, double t);

};

#endif // EXPLICIT_H
