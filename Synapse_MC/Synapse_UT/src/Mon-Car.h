#ifndef MONCAR_H
#define MONCAR_H
#include <armadillo>
using namespace arma;

class Monte_Carlo
{public:

    //constructor
    Monte_Carlo();
    //functions
    void Monte_Carlo_boxes (int n, double dt);
    void Monte_Carlo_vector (int n, double dt);
    void Monte_Carlo_Gauss_boxes (int n, double dt);
    void Monte_Carlo_Gauss_vector (int n, double dt);

};
#endif // MONCAR_H
