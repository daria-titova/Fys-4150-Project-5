#include <iostream>
#include <math.h>
#include <Mon-Car.h>
#include <armadillo>
#include <vector>
#include <random.h>
#include <Random_Gauss.h>
using std::vector;
using namespace std;
using namespace arma;

Monte_Carlo::Monte_Carlo(){}

void Monte_Carlo::Monte_Carlo_vector(int n, double dt)
{
    cout<<"Start MC"<<endl;

    int N=1000; //the number of particles in x_0
    double D=1.0;
    double lo=sqrt(2*D*dt);
    cout<<"lo="<<lo<<endl;
    cout<<"m="<<1/lo<<endl;


    vector<double> u(N,0);

    int time=0;
    while (time < n) {


        for(int i = 0; i < u.size(); ++i)
        { double epsilon = (double)((double)rand() / (RAND_MAX));
            if (epsilon <= 0.5) {

                  if (u[i]==0.0) {
                      u[i]+=0;
                  } else {
                      if ((u[i]-lo)<=0.0)
                      {
                          //u[i]+=0;
                          u.erase(u.begin()+(i));
                          i=i-1;
                      }
                      else
                      {   u.insert(u.begin()+i, u[i]-lo);
                          u.erase(u.begin()+(i+1));

                           }
                  }

            } else {
                if ((u[i]+lo)>=1.0)
                {
                    u.erase(u.begin()+(i));
                    i-=1;
                }
                else {
                    if (u[i]==0.0) {
                       // u.insert(i, 0.0);
                        u.push_back(u[i]+lo);
                        }
                    else
                    {

                        u.push_back(u[i]+lo);
                        u.erase(u.begin()+(i));
                        i-=1;

                        }
                }}  }   //thus we've made some moves
             time++;// cout<<"Next time step"<<time<<endl;
    }

    cout<<"End MC"<<endl;


        // for(int i = 0; i < u.size(); ++i)
        // {cout << (double) u[i] <<endl;}

       ofstream myfile;
       myfile.open ("Mon-Car.txt");
       for (int i=0; i<u.size(); i++)
       {
           if (u[i]==0) u[i]+=0;
           else
               myfile <<u[i]<<endl;
       }
       myfile.close();
return;}









void Monte_Carlo::Monte_Carlo_Algo(int n, double dt)
{
    cout<<"Start MC"<<endl;

    int N=1000; //the number of particles in x_0
    double D=1.0;
    double lo=sqrt(2*D*dt);
    cout<<"lo="<<lo<<endl;
    cout<<"m="<<1/lo<<endl;


    vector<int> u(1/lo,0);
    u[0]=N;


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


    //   for(int i = 0; i < u.size(); ++i)
    //   {cout << (double) u[i]/N <<endl;}

       ofstream myfile;
       myfile.open ("Mon-Car.txt");
       for (int i=0; i<u.size(); i++)
          myfile <<i*lo<<" "<<((double) u[i])<<endl;
          myfile.close();
return;}






void Monte_Carlo::Monte_Carlo_Gauss_boxes(int n, double dt)
{
    cout<<"Start MC"<<endl;

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

       cout<<"End MC"<<endl;


       ofstream myfile;
       myfile.open ("Mon-Car_Gauss.txt");
       for (int i=0; i<u.size(); i++)
           myfile <<(i)*lo<<" "<<(u[i])<<endl;
        myfile.close();

return;}







void Monte_Carlo::Monte_Carlo_Gauss_vector(int n, double dt)
{
    cout<<"Start MC"<<endl;

    int N=1000; //the number of particles in x_0
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

                  if (u[i]==0.0) {
                      u[i]+=0;
                  } else {
                      if ((u[i]-fabs(l))<=0.0)
                      {
                          u.erase(u.begin()+(i));
                          i=i-1;
                      }
                      else
                      {   u.insert(u.begin()+i, u[i]-fabs(l));
                          u.erase(u.begin()+(i+1));

                           }
                  }

            } else { if ((u[i]+l)>=1.0)
                {
                    u.erase(u.begin()+(i));
                    i-=1;
                }
                else {
                    if (u[i]==0.0) {
                        u.push_back(u[i]+fabs(l));
                        }
                    else
                    {

                        u.push_back(u[i]+fabs(l));
                        u.erase(u.begin()+(i));
                        i-=1;

                        }
                }}  }   //thus we've made some moves
             time++;// cout<<"Next time step"<<time<<endl;
    }

    cout<<"End MC"<<endl;


        // for(int i = 0; i < u.size(); ++i)
        // {cout << (double) u[i] <<endl;}

       ofstream myfile;
       myfile.open ("Mon-Car_Gauss_vec.txt");
       for (int i=0; i<u.size(); i++)
       {if (u[i]==0) u[i]+=0;
         else myfile <<u[i]<<endl;}
        //  myfile <<u[i]<<endl;
       myfile.close();
return;}




