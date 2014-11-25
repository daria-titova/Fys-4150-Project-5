#include <iostream>
#include <armadillo>
#include <Explicit.h>
#include <Implicit.h>
#include <Closed_form.h>
#include <Cra-Nic.h>
#include <Mon-Car.h>

using namespace arma;
using namespace std;

int main()
{
    int n, m;
    double a=0.0, b=1.0, t_begin=0.0, t_final;
    cout<<"Insert the number of steps in time"<<endl<<"n=";
    cin>>n;
    cout<<"Insert the end of the interval in time"<<endl<<"t_final=";
    cin>>t_final;
    cout<<"Insert the number of steps for position"<<endl<<"m=";
    cin>>m;

    //define steps for t and x:
    double dt=(t_final-t_begin)/(n-1);
    double dx=(b-a)/(m-1);
    double alpha=dt/(dx*dx);
    cout<<"alpha="<<alpha<<endl;
    cout<<"dt="<<dt<<endl;
    cout<<"dx=dy="<<dx<<endl;

    //create array U
    mat U(n,m);

    //boundary conditions
  /*  U.col(0).fill();
    U.col(m-1).zeros();
    U.row(0).fill(exp());
    U.row(m-1).fill(0.0);*/

    //initial conditions
    for (int i=0; i<m; i++) {
        for (int j=0; j<m; j++) {
            U(i,j)=(1-j*dx)*exp(i*dx);}}

    mat U_Ex(m,m), U_Im(m,m), U_CN(m,m), U_MC(m,m);
    U_Ex=U;    //initial condition, t=0;
    U_Im=U;
    U_MC=U;

    Explicit method;
   // method.Explicit_Scheme(U_Ex, alpha, n, m, dx, t_final);

    Implicit solve;
    solve.Implicit_Scheme(U_Im, alpha, n, m, dx, t_final);

    Crank_Nicolson result;
   // result.Crank_Nicolson_Scheme(V_CN, alpha, n, m, dx);


    Monte_Carlo flip;
   // flip.Monte_Carlo_boxes(n, dx, dt, t_final);

    Closed_form test;
    test.Closed_form_solution(n, m, t_final/20, dx);

    cout<<"The end"<<endl;
    cout<<"End MC"<<endl;

    return 0;
}
