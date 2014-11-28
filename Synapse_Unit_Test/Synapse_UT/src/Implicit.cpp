#include <Implicit.h>
#include <iostream>
#include <math.h>
#include <armadillo>
#include <tridiag.h>
#include <omp.h>
using namespace arma;
using namespace std;

Implicit::Implicit()
{
}

void Implicit::Implicit_Scheme (mat &U, double alpha, int n, int m, double dx, double t)
{

#pragma omp parallel
    {
        int k=0;
        mat U_new(m,m);
        U_new.zeros();
        double diff;

        U_new=U;

        /*    while (k<1000){

        diff=0.0;

        for (int i=0; i<m; i++) {
             for (int j=0; j<m; j++)
             {if (i==0) U_new(i,j)=(1-j*dx)*exp(t);      //boundary condition
              else{ if (i==(m-1)) U_new(i,j)=(1-j*dx)*exp(1+t);  //boundary condition
                else {if (j==0) U_new(i,j)=exp(i*dx+t);
                    else {if (j==m-1) U_new(i,j)=0.0;
                        else U_new(i,j)=(1/(1+4*alpha))*(alpha*(U_new(i-1,j)+U_new(i+1,j)+U_new(i,j-1)+U_new(i,j+1))+U(i,j));
                             diff+=fabs(U(i,j)-U_new(i,j));}}}}}
        U=U_new;
        k++;
        diff/=pow((m),2.0);}*/


        while (k<=1000){

            diff=0.0;

            for (int i=0; i<m-2; i++) {
                for (int j=0; j<m-2; j++)
                {U_new(i+1,j+1)=(1/(1+4*alpha))*(alpha*(U_new(i,j+1)+U_new(i+2,j+1)+U_new(i+1,j)+U_new(i+1,j+2))+U(i+1,j+1));
                    diff+=fabs(U(i+1,j+1)-U_new(i+1,j+1));}}
            U=U_new;
            k++;
            diff/=pow((m),2.0);}



        // U.print("Implicit=");


        ofstream myfile;
        myfile.open ("Implicit.txt");
        for (int i=0; i<m; i++) {
            for (int j=0; j<m; j++)
                myfile <<i*dx<<" "<<j*dx<<" "<<U(i,j)<<endl;
            myfile<<endl; }
        myfile.close();
    }

   return;
}

