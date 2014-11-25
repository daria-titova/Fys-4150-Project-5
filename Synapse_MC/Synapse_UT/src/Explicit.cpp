#include <Explicit.h>
#include <iostream>
#include <math.h>
#include <armadillo>
using namespace arma;
using namespace std;

Explicit::Explicit()
{
}

void Explicit::Explicit_Scheme (vec &V, double alpha, int n, int m, double dx)
{
    mat A(m,m);
    A.zeros();
    A.diag().fill(1-2*alpha);
    A.diag(1).fill(alpha);
    A.diag(-1).fill(alpha);
    //A.print();

    int j=1;
    vec V_new(m);
    while (j<n){
    for (int i=0; i<m; i++)
    {if (i==0) V_new(i)=1.0;      //boundary condition
     else{ if (i<(m-1)) V_new(i)=alpha*V(i-1) + (1-2*alpha)*V(i) + alpha*V(i+1);
           else V_new(i)=0.0;}    //boundary condition
        } V=V_new;
    j++;}

   // V.print("Explicit=");

    ofstream myfile;
    myfile.open ("Explicit.txt");
    for (int i=0; i<m; i++)
       myfile <<i*dx<<" "<<V(i)<<endl;
       myfile.close();

   return;
}

