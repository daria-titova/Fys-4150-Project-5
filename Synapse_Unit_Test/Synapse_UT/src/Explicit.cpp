#include <Explicit.h>
#include <iostream>
#include <math.h>
#include <armadillo>
#include <Initialize.h>
using namespace arma;
using namespace std;

Explicit::Explicit()
{
}

void Explicit::Explicit_Scheme (mat &U, double alpha, int n, int m, double dx, double t)
{
    int k=1;
    mat U_new(m,m);
    U_new.zeros();

    /*   while (k<n){
    for (int i=0; i<m; i++) {
      for (int j=0; j<m; j++)
           {
          if (i==0) U_new(i,j)=(1-j*dx)*exp(t);      //boundary condition
          else{ if (i==(m-1)) U_new(i,j)=(1-j*dx)*exp(1+t);  //boundary condition
              else {if (j==0) //U_new(i,j)=U(i,j) + alpha*(U(i+1,j) + U(i-1,j) + U(i,j+1) - 4*U(i,j));
                      U_new(i,j)=exp(i*dx+t);
                  else {if (j==m-1) //U_new(i,j)=U(i,j) + alpha*(U(i+1,j) + U(i-1,j) + U(i,j-1) - 4*U(i,j));
                          U_new(i,j)=0.0;
                      else U_new(i,j)=U(i,j) + alpha*(U(i+1,j) + U(i-1,j) + U(i,j+1) + U(i,j-1) - 4*U(i,j));}
                    }
              }
            }} U=U_new; k++;}*/

    while (k<n){
        for (int i=0; i<m; i++)
        {for (int j=0; j<m; j++)
            {   if (i==0) U_new(i,j)=1.0;      //boundary condition
                else{ if (i==(m-1)) U_new(i,j)=0.0;  //boundary condition
                    else {if (j==0) U_new(i,j)=U(i,j) + alpha*(U(i+1,j) + U(i-1,j) + 2*U(i,j+1) - 4*U(i,j));
                        else {if (j==m-1) U_new(i,j)=U(i,j) + alpha*(U(i+1,j) + U(i-1,j) + 2*U(i,j-1) - 4*U(i,j));
                            else U_new(i,j)=U(i,j) + alpha*(U(i+1,j) + U(i-1,j) + U(i,j+1) + U(i,j-1) - 4*U(i,j));}
                    }}}}   U=U_new; k++;}

    // V.print("Explicit=");

    ofstream myfile;
    myfile.open ("Explicit.txt");
    for (int i=0; i<m; i++) {
        for (int j=0; j<m; j++)
            myfile <<i*dx<<" "<<j*dx<<" "<<U(i,j)<<endl;
        myfile<<endl; }
    myfile.close();

    Initialize data;
    data.print_out("Explicit.txt", U, m, dx);


    return;
}

