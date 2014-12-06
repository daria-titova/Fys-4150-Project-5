#ifndef INITIALIZE_H
#define INITIALIZE_H
#include <armadillo>
using namespace arma;

class Initialize {
public:
    Initialize();

    void insert(int &n, int &m, double &t, double &dt, double &dx, double &alpha);
    void print_out (char *file_name, mat &V, int &m, double &dx);
};

#endif // INITIALIZE_H
