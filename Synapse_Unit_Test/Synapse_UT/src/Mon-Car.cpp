#include <iostream>
#include <math.h>
#include <Mon-Car.h>
#include <armadillo>
#include <Random.h>
#include <omp.h>
using namespace std;
using namespace arma;

Monte_Carlo::Monte_Carlo(){}

void Monte_Carlo::Monte_Carlo_boxes(int n, double dt)
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
                         {if (j==0) { u(i,j)+=0;
                                     l-=1;}
                         else { if ((j-1) == 0) {u(i,j)-=1;
                                                if (u(i,j-1)>=N) u(i,j-1)+=0;
                                                else u(i,j-1)+=1;}
                                else { u(i,j)-=1;
                                       u(i,j-1)+=1;}
                               l=u(i,j);}}

                      else { if (j==0) {u(i,j+1)+=1;
                                        l-=1;}
                             else { if (j==m-1) { u(i,j)-=1;}
                                    else { u(i,j)-=1;
                                         u(i,j+1)+=1;}
                                   l=u(i,j);}}}

                  else             //jump on y (=j)
                   {double epsilon = (double)((double)rand() / (RAND_MAX));
                     if (epsilon <= 0.5)
                       {if (i==0) { if (j<=m-2) {u(i,j)-=0;      //was j==0
                                               l-=1;}
                                    else { u(i,j)-=1;
                                          l=u(i,j);}}

                           else {if (j==0) { {if (i-1==0) u(i,j)-=0;
                                             else {if (u(i-1,j)>=N) u(i-1,j)+=0;
                                                  else u(i-1,j)+=1;}}
                                            l-=1;}

                                 else {u(i,j)-=1;
                                       u(i-1,j)+=1;
                                       l=u(i,j);} } }

                     else { if (i==0) { if (j==0) { {if (u(i+1,j)>=N) u(i+1,j)+=0;
                                                     else u(i+1,j)+=1;}
                                                   l-=1;}
                                        else {u(i,j)-=1;
                                              u(i+1,j)+=1;
                                              l=u(i,j);}}
                            else
                            {  if (i==m-1) {u(i,j)-=1;
                                          l=u(i,j);
                                          if (j<=m-2) u(i,j)+=1;     //was j==0
                                          l-=1;}
                             else { if (j==0) { {if (u(i+1,j)>=N) u(i+1,j)+=0;
                                                else u(i+1,j)+=1;}
                                               l-=1;}
                                    else {u(i,j)-=1;
                                          u(i+1,j)+=1;
                                          l=u(i,j);} } } } } } } }
        time++; }

       ofstream myfile;
       myfile.open ("Mon-Car.txt");
       for (int i=0; i<m; i++){
           for (int j=0; j<m; j++)
           myfile <<i*lo<<" "<<j*lo<<" "<<u(i,j)<<endl;
           myfile<<endl;}
          myfile.close();
return;}




void Monte_Carlo::Monte_Carlo_boxes_Gauss(int n, double dt)

{

    cout<<"Start MC"<<endl;

    int N=500; //the number of particles in (x=0; y=0)
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

                     rand_gauss rn;
                     double l=lo*rn.rand();
                     int p = fabs(l/lo);
                   //  int p = (l/lo);
                   //  if (((double)l/lo-p)>=0.5) p+=1;

                      if ((epsilon <= 0.5 && l>=0.0) || (epsilon > 0.5 && l<0.0))
                         {if (j==0) { u(i,j)+=0;
                                     l-=1;}
                         else { if ((j-p) == 0) {u(i,j)-=1;
                                                if (u(i,j-p)>=N) u(i,j-p)+=0;
                                                else u(i,j-p)+=1;}
                                else { u(i,j)-=1;
                                  if (j-p<0) u(i,j)+=0;
                                  else
                                       u(i,j-p)+=1;}
                               l=u(i,j);}}

                      else { if (j==0 && j+p>0 && j+p<m) {u(i,j+p)+=1;
                                        l-=1;}
                          else { if (p==0 || j+p>m-1) u(i,j)+=0;

                                 else{ if (j==m-1 || j+p>=m-1) { u(i,j)-=1;}
                                    else { u(i,j)-=1;
                                         u(i,j+p)+=1;}
                                   l=u(i,j);}}}}

                  else             //jump on y (=j)
                   {double epsilon = (double)((double)rand() / (RAND_MAX));

                     rand_gauss rn;
                     double l_new=lo*rn.rand();
                     int p = fabs(l_new/lo);
                    // if (((double)l_new/lo-p)>=0.5) p+=1;

                     if ((epsilon <= 0.5 && l>=0.0) || (epsilon > 0.5 && l<0.0))
                       {if (i==0) { if (j<=m-2) {u(i,j)-=0;  //was j==0
                                               l-=1;}
                                    else {u(i,j)-=1;
                                          l=u(i,j);}}

                           else {if (j==0) { {if (i-p<=0) u(i,j)-=0;
                                             else {if (u(i-p,j)>=N) u(i-p,j)+=0;
                                                  else u(i-p,j)+=1;}}
                                            l-=1;}

                                 else {if (i-p<=0) u(i,j)-=0;
                                 else {u(i,j)-=1;
                                       u(i-p,j)+=1;
                                       l=u(i,j);}} } }

                     else { if (i==0) { if (j==0 && i+p>0 && i+p<m) { {if (u(i+p,j)>=N) u(i+p,j)+=0;
                                                     else u(i+p,j)+=1;}
                                                   l-=1;}
                                    else { if (p==0 || i+p>m-1) u(i,j)+=0;
                                          else        { u(i,j)-=1;
                                                      u(i+p,j)+=1;
                                                      l=u(i,j);}}}
                            else
                            {  if (i==m-1 || i+p>m-1) {u(i,j)-=1;
                                          l=u(i,j);
                                          if (j<=m-2) u(i,j)+=1;   //was j==0
                                          l-=1;}
                             else { if (j==0) { if (i+p>m-1) u(i,j)+=0;
                                             else
                                              {if (i+p>m-1 || u(i+p,j)>=N) u(i,j)+=0;
                                                else u(i+p,j)+=1;}
                                               l-=1;}
                                    else {u(i,j)-=1;
                                             if (i+p>m-1) u(i,j)+=0;
                                             else
                                          u(i+p,j)+=1;
                                          l=u(i,j);} } } } } } } }
        time++; }

       ofstream myfile;
       myfile.open ("Mon-Car.txt");
       for (int i=0; i<m; i++){
           for (int j=0; j<m; j++)
           myfile <<i*lo<<" "<<j*lo<<" "<<u(i,j)<<endl;
       myfile<<endl;}
          myfile.close();


return;}

