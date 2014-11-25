#include <tridiag.h>
#include <iostream>
#include <math.h>
#include <armadillo>
using namespace std;
using namespace arma;

tridiag::tridiag(){}

void tridiag::tridiag_solver(mat &A, vec &V, int m)
{ vec a(m);
    for (int i=0; i<m; i++){
            if (i==0) a(i)=0;
            else a(i)=A(i,i-1);}
  vec b(m);
  b=A.diag();

  vec c(m);
    for (int i=0; i<m; i++){
        if (i==m-1) c(i)=0;
        else c(i)=A(i,i+1);}

    //forward substitution
  vec c1(m);
    for (int i=0; i<m-1; i++){
        if (i==0) c1(i)=c(i)/b(i);
        else c1(i)=c(i)/(b(i)-a(i)*c1(i-1));}

  vec f1(m);
    for (int i=0; i<m; i++){
            if (i==0) f1(i)=V(i)/b(i);
            else f1(i)=(V(i)-a(i)*f1(i-1))/(b(i)-a(i)*c1(i-1));}

    //backward substitution
        for (int i=m-1; i>=0; i--){
            if (i==m-1) V(i)=f1(i);
                else V(i)=f1(i)-c1(i)*V(i+1);}

 return;
}

