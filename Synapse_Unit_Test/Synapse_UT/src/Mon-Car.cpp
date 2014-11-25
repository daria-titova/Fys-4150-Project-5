#include <iostream>
#include <math.h>
#include <Mon-Car.h>
#include <armadillo>
using namespace std;
using namespace arma;

Monte_Carlo::Monte_Carlo(){}

void Monte_Carlo::Monte_Carlo_boxes(int n, double dx, double dt, double t)
{
    cout<<"Start MC"<<endl;

    int N=1000; //the number of particles in (x=0; y=0)
    double D=1.0;
    double lo=sqrt(2*D*dt);
    int m=1/lo;
    cout<<"lo="<<lo<<endl;
    cout<<"m="<<m<<endl;


    mat u(m,m);
    u.zeros();
    u(0,0)=N;


    int time=0;
    while (time < n) {


        for(int i = 0; i < m; i++){
            for (int j=0; j<m; j++)

            { int l=u(i,j);
               for (int k=0; k<l; k++)
               {  double direction = (double)((double)rand() / (RAND_MAX));

                   if (direction <=0.5) //then jump on x (=j)

                   {double epsilon = (double)((double)rand() / (RAND_MAX));
                       if (epsilon <= 0.5)
                       { if (j==0) { u(i,j)+=0;}
                         else { u(i,j)-=1;
                                u(i,j-1)+=1;

                                if ((j-1) == 0) u(i,j-1)-=1;

                                l=u(i,j);}}

                       else { if (j+1 >= m) u(i,j)+=0;
                            else { if (j==0) { u(i,j+1)+=1;
                                               l-=1; }
                                  else { u(i,j)-=1;
                                         u(i,j+1)+=1;
                                         l=u(i,j); }}}}

                       else             //jump on y (=j)
                   {double epsilon = (double)((double)rand() / (RAND_MAX));
                       if (epsilon <= 0.5) {  if (i==0) {u(i,j)-=0;
                                                         l=u(i,j);}
                                             else
                           {
                               u(i,j)-=1;
                               u(i-1,j)+=1;
                               l=u(i,j);
                           } }

                       else { if (i==m-1) {
                               u(i,j)-=1;
                               l=u(i,j);
                           } else {
                               u(i,j)-=1;
                               u(i+1,j)+=1;
                               l=u(i,j);}} }

                 }   //thus we've made some moves
                 }}
        time++;// cout<<"Next time step"<<time<<endl;
    }

    cout<<"End MC"<<endl;


       ofstream myfile;
       myfile.open ("Mon-Car.txt");
       for (int i=0; i<m; i++){
           for (int j=0; j<m; j++)
           myfile <<i*dx<<" "<<j*dx<<" "<<u(i,j)<<endl;}
          myfile.close();
return;}
