#ifndef IMPLICIT_H
#define IMPLICIT_H
#include <armadillo>
using namespace arma;


class Implicit
{public:

    //constructor
    Implicit();

    //functions
    void Implicit_Scheme (mat &U, double alpha, int n, int m, double dx, double t);

};

#endif // IMPLICIT_H
