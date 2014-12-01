#include <iostream>
#include <armadillo>
#include <Explicit.h>
#include <Implicit.h>
#include <Closed_form.h>
#include <Cra-Nic.h>
#include <Mon-Car.h>
#include <initialize.h>
using namespace arma;
using namespace std;

int main()
{
    int n, m;
    double t_final, dt, dx, alpha;

    Initialize variables;
    variables.insert(n, m, t_final, dt, dx, alpha);

    vec V_Ex(m), V_Im(m), V_CN(m), V_MC(m);
    V_Ex.zeros();    //initial condition, t=0;
    V_Im.zeros();
    V_CN.zeros();
    V_MC.zeros();

    Explicit method;
   // method.Explicit_Scheme(V_Ex, alpha, n, m, dx);

    Implicit solve;
   // solve.Implicit_Scheme(V_Im, alpha, n, m, dx);

    Crank_Nicolson result;
   // result.Crank_Nicolson_Scheme(V_CN, alpha, n, m, dx);

    Monte_Carlo flip;
    flip.Monte_Carlo_Gauss_vector(n, dt);

    Closed_form plot;
   // plot.Closed_form_solution(n, m, t_final, dx);

    return 0;
}
