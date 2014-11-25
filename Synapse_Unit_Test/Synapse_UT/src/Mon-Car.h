#ifndef MONCAR_H
#define MONCAR_H
#include <armadillo>
using namespace arma;


class Monte_Carlo
{public:

    //constructor
    Monte_Carlo();

    //functions
    void Monte_Carlo_boxes (int n, double dx, double dt, double t);

};

#endif // MONCAR_H
