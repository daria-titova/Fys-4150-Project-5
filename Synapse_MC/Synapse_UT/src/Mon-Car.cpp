#include <iostream>
#include <math.h>
#include <Mon-Car.h>
#include <armadillo>
#include <vector>
#include <Random_Gauss.h>
#include <initialize.h>
#include <omp.h>
using std::vector;
using namespace std;
using namespace arma;

Monte_Carlo::Monte_Carlo(){}

void Monte_Carlo::Monte_Carlo_vector(int n, double dt)
{
    int N=1000; //the number of particles in x_0
    double D=1.0;
    double lo=sqrt(2*D*dt);
    cout<<"lo="<<lo<<endl;
    cout<<"m="<<1/lo<<endl;  //the number of boxes for the histogram

    vector<double> u(N,0);

    int time=0;
    while (time < n) {

//#pragma omp parallel num_threads(2)

        for(int i = 0; i < u.size(); ++i)
        { //epsilon is random from [0;1]
            double epsilon = (double)((double)rand() / (RAND_MAX));
            if (epsilon <= 0.5) {

                Monte_Carlo make;
                make.negative_move(u, i, lo);

            } else {

                Monte_Carlo jump;
                jump.positive_move(u, i, lo);

            }  }
             time++; }

       ofstream myfile;
       myfile.open ("Mon-Car.txt");
       for (int i=0; i<u.size(); i++)
          { if (u[i]==0) u[i]+=0;
            else  myfile <<u[i]<<endl; }
       myfile.close();
return;}


void Monte_Carlo::Monte_Carlo_Gauss_vector(int n, double dt)
{
    int N=100; //the number of particles in x_0
    double D=1.0;
    double lo=sqrt(2*D*dt);
    cout<<"lo="<<lo<<endl;
    cout<<"m="<<1/lo<<endl;

    vector<double> u(N,0);

    int time=0;
    while (time < n)
    { for(int i = 0; i < u.size(); ++i)
        { double epsilon = (double)((double)rand() / (RAND_MAX));

            rand_gauss rn;
            double r=rn.rand();
            double l=lo*r;

            if ((epsilon <= 0.5 && l>=0.0) || (epsilon > 0.5 && l<0.0)){

               /* l=fabs(l);
                Monte_Carlo make;
                make.negative_move(u, i, l);*/

                  if (u[i]==0.0) {
                      u[i]+=0;
                  } else {
                      if ((u[i]-fabs(l))<=0.0)
                      {
                          u.erase(u.begin()+(i));
                          i=i-1;
                      }
                      else
                      {   u.insert(u.begin(), u[i]-fabs(l));
                          u.erase(u.begin()+(i+1));

                           }
                  } //stop

            } else {

              /* l=fabs(l);
                Monte_Carlo jump1;
                jump1.positive_move(u, i, l);*/

                if (u[i]==0.0){ if (fabs(l)>=1.0) u[i]+=0;

                    else {//u.push_back(u[i]+fabs(l));
                          u.insert(u.begin(), u[i]+fabs(l));
                          i+=1;
                    }
                    }
                else {

                    if ((u[i]+fabs(l))>=1.0)
                    {
                        u.erase(u.begin()+(i));
                        i-=1;
                    }
                     else {
                        u.insert(u.begin(), u[i]+fabs(l));
                        u.erase(u.begin()+(i+1));

                        }} //stop
            }
        }
             time++; }

       ofstream myfile;
       myfile.open ("Mon-Car.txt");
       for (int i=0; i<u.size(); i++)
       {if (u[i]==0) u[i]+=0;
         else myfile <<u[i]<<endl;}
       myfile.close();

return;}


void Monte_Carlo::negative_move(vector<double> &u, int &i, double &l)
{if (u[i]==0.0) {
        u[i]+=0;
    } else {
        if ((u[i]-l)<=0.0)
        {   u.erase(u.begin()+(i));
            i=i-1;
        }
        else
        {   u.insert(u.begin(), u[i]-l);
            u.erase(u.begin()+(i+1)); } }
    return;
}


void Monte_Carlo::positive_move(vector<double> &u, int &i, double &l)
{
    if (u[i]==0.0) { if (l>=1.0) u[i]+=0;
                    else { u.insert(u.begin(), u[i]+l);
                          i+=1;} }
    else {
        if ((u[i]+l)>=1.0)    //was >=
        {   u.erase(u.begin()+(i));
            i-=1;}
        else {
            u.insert(u.begin(), u[i]+l);
            u.erase(u.begin()+(i+1));
            }}
return;}











void Monte_Carlo::Monte_Carlo_boxes(int n, double dt)
{

/*#pragma omp parallel num_threads(1)
    {cout<<"Lev durak, ne hochet pomogat' kisushe"<<endl;}*/

    int N=1000; //the number of particles in x_0
    double D=1.0;
    double lo=sqrt(2*D*dt);
    int m=1/lo;
    cout<<"lo="<<lo<<endl;
    cout<<"m="<<m<<endl;   //the number of boxes


    vector<int> u(m,0);    //the vector of boxes
    u[0]=N;                //N walkers in the first box (initial condition)


    int time=0;
    while (time < n) {

        for(int i = 0; i < u.size(); i++)
        { int l=u[i];
            for (int j=0; j<l; j++)
           {// epsilon is random between [0;1]
            double epsilon = (double)((double)rand() / (RAND_MAX));

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
                if (i+1 >= m)
                    u[i]-=1;
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
                }}  }
             }time++;}

    ofstream myfile;
       myfile.open ("Mon-Car.txt");
       for (int i=0; i<u.size(); i++)
          myfile <<i*lo<<" "<<((double) u[i]/N)<<endl;
          myfile.close();
return;}



void Monte_Carlo::Monte_Carlo_Gauss_boxes(int n, double dt)
{
    int N=1000; //the number of particles in x_0
    double D=1.0;
    double lo=sqrt(2*D*dt);
    cout<<"lo="<<lo<<endl;
    cout<<"1/lo="<<1/lo<<endl;

    int m;
    m=1/lo+1;
    cout<<"m="<<m<<endl;
    double dm=1/m;

    vector<int> u(m,0);
    u[0]=N;

    int time=0;
    while (time < n) {


        for(int i = 0; i < u.size(); ++i)
        { int l=u[i];
            for (int j=0; j<l; j++)
            { rand_gauss rn;

              double l=lo*rn.rand();
              int p = fabs(l/lo);
              double epsilon = (double)((double)rand() / (RAND_MAX));

                if ((epsilon <= 0.5 && l>=0.0) || (epsilon > 0.5 && l<0.0)) {

                  if (i==0) {
                      u[i]+=0;
                  } else {
                      u[i]-=1;
                      u[i-p]+=1;
                      if ((i-p) <= 0) u[i-p]-=1;
                  }

            } else {
                if (i+p >= u.size())
                    u[i]+=0;
                else {
                    if (i==0) {
                        if (p==0) u[i+p]+=0;
                        else {
                            u[i+p]+=1;
                            l-=1;
                        }
                    }
                    else
                    {
                        u[i]-=1;
                        u[i+p]+=1;
                    }
                }}  }
             }
        time++; }

       ofstream myfile;
       myfile.open ("Mon-Car_Gauss.txt");
       for (int i=0; i<u.size(); i++)
           myfile <<(i)*lo<<" "<<(u[i])<<endl;
        myfile.close();

return;}

