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

    int N=10; //the number of particles in (x=0; y=0)
    double D=1.0;
    double lo=sqrt(2*D*dt);
    double m=1/lo;
    cout<<"lo="<<lo<<endl;
    cout<<"m="<<m<<endl;


    mat u(m,m);
    u.zeros();
    u(0,0)=N;


    int time=0;
    while (time < n) {


        for(int i = 0; i < u.size(); ++i)
        { int l=u[i];
            for (int j=0; j<l; j++)
           {double epsilon = (double)((double)rand() / (RAND_MAX));
            if (epsilon <= 0.5) {

                  if (i==0) {
                      u[i]+=0;
                  } else {
                      u[i]-=1;
                      u[i-1]+=1;
                      if ((i-1) == 0) u[i-1]-=1;
                      l=u[i];
                  }

            } else {
                if (i+1 >= 1.0/lo)
                    u[i]+=0;
                else {
                    if (i==0) {
                        u[i+1]+=1;
                        l-=1;
                    }
                    else
                    {
                        u[i]-=1;
                        u[i+1]+=1;
                        l=u[i];
                    }
                }}  }   //thus we've made some moves
             }time++;// cout<<"Next time step"<<time<<endl;
    }

    cout<<"End MC"<<endl;


       ofstream myfile;
       myfile.open ("Mon-Car.txt");
       for (int i=0; i<u.size(); i++)
          myfile <<i*lo<<" "<<((double) u[i]/N)<<endl;
          myfile.close();
return;}
