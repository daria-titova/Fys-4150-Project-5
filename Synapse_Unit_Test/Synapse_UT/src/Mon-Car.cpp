#include <iostream>
#include <math.h>
#include <Mon-Car.h>
#include <armadillo>
#include <Random.h>
#include <omp.h>
using namespace std;
using namespace arma;

Monte_Carlo::Monte_Carlo(){}

void Monte_Carlo::Monte_Carlo_vector_Gauss(int n, double dt)
{
    int N=100; //the number of particles in x_0
    double D=1.0;
    double lo=sqrt(2*D*dt);
    int m=1/lo;
    cout<<"lo="<<lo<<endl;
    cout<<"m="<<m<<endl; //the number of boxes for the histogram

    vector<double> x(N,0);
    vector<double> y(N,0);
    vec test_u(m);
    test_u.zeros();


    double x_old, y_old;
    double null=1e-15;
    int time=0;

    while (time < n) {
        //#pragma omp parallel num_threads(2)
        for(int i = 0; i < x.size(); i++)
        { //epsilon is random from [0;1]
            double epsilon_x = (double)((double)rand() / (RAND_MAX));
            double epsilon_y = (double)((double)rand() / (RAND_MAX));

            //rand_gauss rn;
            rand_gauss rnx;
            rand_gauss rny;
            double rx=rnx.rand();
            double ry=rny.rand();
            double lx=lo*rx;
            double ly=lo*ry;

            x_old=x[i];
            y_old=y[i];

            if ((epsilon_x <= 0.5 && lx>=0.0) || (epsilon_x > 0.5 && lx<0.0)) {
                x.insert(x.begin(), x[i]-fabs(lx));
                x.erase(x.begin()+(i+1));
            } else {
                x.insert(x.begin(), x[i]+fabs(lx));
                x.erase(x.begin()+(i+1));
            }
            if ((epsilon_y <= 0.5 && ly>=0.0) || (epsilon_y > 0.5 && ly<0.0)) {
                y.insert(y.begin(), y[i]-fabs(ly));
                y.erase(y.begin()+(i+1));
            } else {
                y.insert(y.begin(), y[i]+fabs(ly));
                y.erase(y.begin()+(i+1));
            }

            double x_new =x[0];
            double y_new =y[0];

            int k=0.0;
            if (y_new>null && y_new<1.0) //if y=(0;1)
                k=y_new*m;


            if ((x_new<(1/(double)m) && x_new > 0.0) && (y_new>0.0 && y_new<1.0))
            {   if (test_u(k)<N)  //preserv N at (x=0; y)
                    test_u(k)++;
                else {x.erase(x.begin());
                    y.erase(y.begin());} }


            if  (x_new<0.0 || x_new>=1.0      //if x<0 or x>=1, the boundary condiitons
                 || y_new<0.0 || y_new>=1.0    //if y<0 or y>=1, the boundary condiitons
                 ||((x_new<=null && x_new>=0.0) && (y_new<=null && y_new>=0.0)))   // if the move was back to (x=0; y=0)
            {   x.erase(x.begin());
                y.erase(y.begin());
                i-=1;}


            if (x_old<(1/(double)m) && x_old>0.0)
            {x.insert(x.begin(), x_old);
                y.insert(y.begin(), y_old);  //increment test_u if not (0,0)
                i+=1;}


            if  (  (x_old<=null && x_old>=0.0)        //if the move was from x=0
                   )
            {   x.insert(x.begin(), x_old);
                y.insert(y.begin(), y_old);  //increment test_u if not (0,0)
                i+=1;}

            if  (  (y_new<0.0)     //if the move is to y<0 and..
                   && ((x_old>=(1/(double)m) || (x_old>null && x_old<(1/(double)m) && test_u((int)y_old*m)<N))
                       && (y_old>=null)) )   //and x_old=(0;1) and there are not N elements yet in the k_th box
            {   //x.insert(x.begin(), x_new);
                //y.insert(y.begin(), y_new - int(y_new/1) +1.0);
                x.insert(x.begin(), x_old);
                y.insert(y.begin(), y_old);
                i+=1;}

            if  (  (y_new>1.0)   //if the move is to y>1.0 and..
                   // && (x_new>null && x_new<1.0)
                   //and x_old=(0;1) and there are not N elements yet in the k_th box
                   && (x_old>=(1/(double)m) || (x_old>null && x_old<lo && test_u((int)y_old*m)<N)))
            {  // x.insert(x.begin(), x_new);
                // y.insert(y.begin(), y_new-int(y_new/1));
                x.insert(x.begin(), x_old);
                y.insert(y.begin(), y_old);
                i+=1;}


        }
        time++; }

    test_u.print();

    ofstream myfile;
    myfile.open ("Mon-Car.txt");
    for (int i=0; i<x.size(); i++)
    { if (x[i]<=null || y[i]<=null) x[i]+=0;
        else
            myfile <<x[i]<<" "<<y[i]<<endl;}
    myfile.close();
    return;}





void Monte_Carlo::Monte_Carlo_vector(int n, double dt)
{
    int N=500; //the number of particles in x_0
    double D=1.0;
    double lo=sqrt(2*D*dt);
    int m=1/lo;
    cout<<"lo="<<lo<<endl;
    cout<<"m="<<m<<endl; //the number of boxes for the histogram

    vector<double> x(N,0);
    vector<double> y(N,0);
    vec test_u(m);


    double x_old, y_old;
    double null=1e-10;
    int time=0;

    while (time < n) {
        //#pragma omp parallel num_threads(2)
        for(int i = 0; i < x.size(); i++)
        { //epsilon is random from [0;1]
            double epsilon_x = (double)((double)rand() / (RAND_MAX));
            double epsilon_y = (double)((double)rand() / (RAND_MAX));

            x_old=x[i];
            y_old=y[i];

            if (epsilon_x <= 0.5) {
                x.insert(x.begin(), x[i]-lo);
                x.erase(x.begin()+(i+1));
            } else {
                x.insert(x.begin(), x[i]+lo);
                x.erase(x.begin()+(i+1));
            }
            if (epsilon_y <= 0.5) {
                y.insert(y.begin(), y[i]-lo);
                y.erase(y.begin()+(i+1));
            } else {
                y.insert(y.begin(), y[i]+lo);
                y.erase(y.begin()+(i+1));
            }

            double x_new =x[0];
            double y_new =y[0];


            if ((x_new<=null && x_new>=0.0) && (y_new>null && y_new<1.0))
            {int k=y_new*m;

                if (test_u(k)<N)  //preserv N at (x=0; y)
                    test_u(k)+=1;
                else {x.erase(x.begin());
                    y.erase(y.begin());} }


            if ((x_new<=null && x_new!=0.0) || x_new<0.0 || y_new<0.0 || (y_new<=null && y_new!=0.0)
                    || x_new>=1.0 || y_new>=1.0 || ((x_new<=null && x_new>=0.0) && (y_new<=null && y_new>=0.0)) )
            {   x.erase(x.begin());
                y.erase(y.begin());
                i-=1;

            }

            if ((x_old<=null && x_old>=0.0))
            {   x.insert(x.begin(), x_old);
                y.insert(y.begin(), y_old);
                i+=1;}

            if ((y_new<0.0 || (y_new<=null && y_new!=0.0)) && (x_new>null && x_new<1.0) )
            {   x.insert(x.begin(), x_new);
                y.insert(y.begin(), y_new+1.0);
                i+=1;}

            if ((y_new>=1.0) && (x_new>null && x_new<1.0) )
            {   x.insert(x.begin(), x_new);
                y.insert(y.begin(), y_new-1.0);
                i+=1;}

        }
        time++; }
    ofstream myfile;
    myfile.open ("Mon-Car.txt");
    for (int i=0; i<x.size(); i++)
    { if (x[i]<=null || y[i]<=null) x[i]+=0;
        else
            myfile <<x[i]<<" "<<y[i]<<endl;}
    myfile.close();
    return;}





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
                { double direction = (double)((double)rand() / (RAND_MAX));
                    if (direction <=0.5)
                        //then jump on x (=j)
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

                    else //jump on y (=j)
                    {double epsilon = (double)((double)rand() / (RAND_MAX));
                        if (epsilon <= 0.5)
                        {if (i==0) { if (j<=m-2) {u(i,j)-=0; //was j==0
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
                            { if (i==m-1) {u(i,j)-=1;
                                    l=u(i,j);
                                    if (j<=m-2) u(i,j)+=1; //was j==0
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
    int N=100; //the number of particles in (x=0; y=0)
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
            { int ll=u(i,j);
                for (int k=0; k<ll; k++)
                { double direction = (double)((double)rand() / (RAND_MAX));
                    if (direction <=0.5)
                        //then jump on x (=j)
                    {double epsilon = (double)((double)rand() / (RAND_MAX));
                        rand_gauss rn;
                        double l=lo*rn.rand();
                        int p = fabs(l/lo);

                        if ((epsilon <= 0.5 && l>=0.0) || (epsilon > 0.5 && l<0.0))
                        {if (j==0) { u(i,j)+=0;
                                ll-=1;}
                            else { if ((j-p) == 0) {u(i,j)-=1;
                                    if (u(i,j-p)>=N) u(i,j-p)+=0;
                                    else u(i,j-p)+=1;}
                                else { u(i,j)-=1;
                                    if (j-p<0) u(i,j)+=0;
                                    else
                                        u(i,j-p)+=1;}
                                ll=u(i,j);}}
                        else { if (j==0 && j+p>0 && j+p<m) {u(i,j+p)+=1;
                                ll-=1;}
                            else { if (p==0 || j+p>m-1) u(i,j)+=0;
                                else{ if (j==m-1 || j+p>=m-1) { u(i,j)-=1;}
                                    else { u(i,j)-=1;
                                        u(i,j+p)+=1;}
                                    ll=u(i,j);}}}}
                    else //jump on y (=j)
                    {double epsilon = (double)((double)rand() / (RAND_MAX));
                        rand_gauss rn;
                        double l=lo*rn.rand();   //l_new
                        int p = fabs(l/lo);     //l_new
                        // if (((double)l_new/lo-p)>=0.5) p+=1;

                        if ((epsilon <= 0.5 && l>=0.0) || (epsilon > 0.5 && l<0.0))
                        {if (i==0) { if (j<=m-2) {u(i,j)-=0; //was j==0
                                    ll-=1;}
                                else {u(i,j)-=1;
                                    ll=u(i,j);}}
                            else {if (j==0) { {if (i-p<=0) u(i,j)-=0;
                                        else {if (u(i-p,j)>=N) u(i-p,j)+=0;
                                            else u(i-p,j)+=1;}}
                                    ll-=1;}
                                else {if (i-p<=0) u(i,j)-=0;
                                    else {u(i,j)-=1;
                                        u(i-p,j)+=1;
                                        ll=u(i,j);}} } }

                        else { if (i==0) { if (j==0 && i+p>0 && i+p<m) { {if (u(i+p,j)>=N) u(i+p,j)+=0;
                                        else u(i+p,j)+=1;}
                                    ll-=1;}
                                else { if (p==0 || i+p>m-1) u(i,j)+=0;
                                    else { u(i,j)-=1;
                                        u(i+p,j)+=1;
                                        ll=u(i,j);}}}
                            else
                            { if (i==m-1 || i+p>m-1) {u(i,j)-=1;
                                    ll=u(i,j);
                                    if (j<=m-2) u(i,j)+=1; //was j==0
                                    ll-=1;}
                                else { if (j==0) { if (i+p>m-1) u(i,j)+=0;
                                        else
                                        {if (i+p>m-1 || u(i+p,j)>=N) u(i,j)+=0;
                                            else u(i+p,j)+=1;}
                                        ll-=1;}
                                    else {u(i,j)-=1;
                                        if (i+p>m-1) u(i,j)+=0;
                                        else
                                            u(i+p,j)+=1;
                                        ll=u(i,j);} } } } } } } }
        time++; }
    ofstream myfile;
    myfile.open ("Mon-Car.txt");
    for (int i=0; i<m; i++){
        for (int j=0; j<m; j++)
            myfile <<i*lo<<" "<<j*lo<<" "<<u(i,j)<<endl;
        myfile<<endl;}
    myfile.close();
    return;}




