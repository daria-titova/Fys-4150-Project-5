#ifndef TRIDIAG_H
#define TRIDIAG_H
#include <armadillo>
using namespace arma;

class tridiag
{
public:

    tridiag ();
    void tridiag_solver(mat &A, vec &V, int m);
};

#endif // TRIDIAG_H
