#include <Implicit.h>
#include <iostream>
#include <math.h>
#include <armadillo>
#include <tridiag.h>
#include <initialize.h>
using namespace arma;
using namespace std;

Implicit::Implicit()
{
}

void Implicit::Implicit_Scheme (vec &V, double alpha, int n, int m, double dx)
{
    mat A(m,m);
    A.zeros();
    A.diag().fill(1+2*alpha);
    A.diag(1).fill(-alpha);
    A.diag(-1).fill(-alpha);
   // A.print();

    int j=1;
    while (j<n){

          V(0)=1.0;     //boundary condition
          vec U(m-2);
          for (int k=0; k<m-2; k++)
                   {U(k)=V(k+1);}

                   U(0)+=alpha;
                   tridiag solver;
                   solver.tridiag_solver(A, U, m-2);

                   for (int k=0; k<m-2; k++)
                   {V(k+1)=U(k);}

           V(m-1)=0.0;   //boundary condition


        j++;}

    Initialize data;
    data.print_out("Implicit.txt", V, m, dx);

   return;
}

