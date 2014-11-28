#include <Closed_form.h>
#include <iostream>
#include <math.h>
#include <armadillo>
using namespace arma;
using namespace std;

Closed_form::Closed_form()
{
}

void Closed_form::Closed_form_solution(int n, int m, double t, double dx)
{
    mat U(m,m);
    int k;

    //we choose k to be like this in order to avoid the case when
    //when the exponent in the closed-form solution is not decreasing
    //when t<<1,
    //viz. to kill the exponent in order to get steady-state
    //
    if (t<1) k=1/t;
    else k=n;



   /* for (int i=0; i<m; i++){
        for (int j=0; j<m; j++) {
        double sum=0.0;
        for (int l=0; l<k-1; l++) {
            for (int p=0; p<k-1; p++)
                sum=sum+(2/((l+1)*M_PI))*sin((l+1)*M_PI*i*dx)*exp(-M_PI*M_PI*((l+1)*(l+1)+(p+1)*(p+1))*t);}
        U(i,j)=1.0-i*dx-sum;
        }}

    U(m-1,m-1)=0.0;     //for the sace of printing
   // U.print("Closed_form=");*/

    for (int i=0; i<m; i++){
            for (int j=0; j<m; j++) {
                U(i,j)=(1-j*dx)*exp(i*dx+t);
            }}


    ofstream myfile;
    myfile.open ("Closed.txt");
    for (int i=0; i<m; i++){
        for (int j=0; j<m; j++)
           myfile <<i*dx<<" "<<j*dx<<" "<<U(i,j)<<endl;
        myfile<<endl;   }
       myfile.close();

}

