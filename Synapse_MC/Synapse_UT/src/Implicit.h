#ifndef IMPLICIT_H
#define IMPLICIT_H
#include <armadillo>
using namespace arma;


class Implicit
{public:

    //constructor
    Implicit();

    //functions
    void Implicit_Scheme (vec &V, double alpha, int n, int m, double dx);

};

#endif // IMPLICIT_H
