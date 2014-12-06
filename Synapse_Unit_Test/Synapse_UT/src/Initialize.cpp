#include <Initialize.h>
#include <iostream>
using namespace std;

Initialize::Initialize() {}

void Initialize::insert(int &n, int &m, double &t, double &dt, double &dx, double &alpha)
{
    double a=0.0, b=1.0, t_begin=0.0;
    cout<<"Insert the number of steps in time"<<endl<<"n=";
    cin>>n;
    cout<<"Insert the end of the interval in time"<<endl<<"t_final=";
    cin>>t;
    cout<<"Insert the number of steps for position"<<endl<<"m=";
    cin>>m;

    //define steps for t and x,y:
    dt=(t-t_begin)/(n-1);
    dx=(b-a)/(m-1);
    alpha=dt/(dx*dx);
    cout<<"alpha="<<alpha<<endl;
    cout<<"dt="<<dt<<endl;
    cout<<"dx=dy="<<dx<<endl;

}
void Initialize::print_out(char *file_name, mat &V, int &m, double &dx)
{
    ofstream myfile;
    myfile.open (file_name);
    for (int i=0; i<m; i++)
    {for (int j=0; j<m; j++)
        myfile <<i*dx<<" "<<j*dx<<" "<<V(i,j)<<endl;
    myfile<<endl;}
    myfile.close();}

