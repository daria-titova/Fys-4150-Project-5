#include <iostream>
#include <armadillo>
#include <Explicit.h>
#include <Implicit.h>
#include <Closed_form.h>
#include <Mon-Car.h>
#include <Initialize.h>
using namespace arma;
using namespace std;

int main()
{
    int n, m;
    double t_final, dx, dt, alpha;

    Initialize variables;
    variables.insert(n, m, t_final, dt, dx, alpha);

    cout<<"alpha="<<alpha<<endl;
    cout<<"dt="<<dt<<endl;
    cout<<"dx=dy="<<dx<<endl;

    //create array U
    mat U(m,m);

    //initial conditions
   /* for (int i=0; i<m; i++) {
        for (int j=0; j<m; j++) {
            U(i,j)=(1-j*dx)*exp(i*dx);}}*/


    mat U_Ex(m,m), U_Im(m,m), U_MC(m,m);
    U.zeros();//initial condition, t=0;
    U_Ex=U;
    U_Im=U;
    U_MC=U;

    Explicit method;
  //  method.Explicit_Scheme(U_Ex, alpha, n, m, dx, t_final);

    Implicit solve;
   // solve.Implicit_Scheme(alpha, m, dx, t_final);

    Monte_Carlo flip;
    flip.Monte_Carlo_boxes_Gauss(n, dt);

    Closed_form test;
  //  test.Closed_form_solution(n, m, t_final, dx);


    return 0;
}
