#include <unittest++/UnitTest++.h>
#include <tridiag.h>
#include <Explicit.h>
#include <Implicit.h>
#include <Cra-Nic.h>
#include <iostream>
#include <Mon-Car.h>
#include <vector>
using std::vector;

TEST(tridiagonal_solver) {
    mat a;
    a.zeros(4,4);
    a.diag().fill(2);
    a.diag(1).fill(1);
    a.diag(-1).fill(1);

    vec b(4);
    b.fill(2);

    vec x(4);
    x(0)=0.8;
    x(1)=0.4;
    x(2)=0.4;
    x(3)=0.8;

tridiag my;
my.tridiag_solver(a, b, 4);

vec dif=abs(b-x);
CHECK(dif.max() < 1e-5);
}


TEST(Explicit_solver){
    vec v(4);
    v.zeros();

    vec result(4);
    result(0)=1.0;
    result(1)=0.42;
    result(2)=0.09;
    result(3)=0.0;

Explicit my;
my.Explicit_Scheme(v, 0.3, 4, 4, 1.0/3.0);
vec dif=abs(result-v);

    CHECK(dif.max() < 1e-4);}


TEST(Implicit_solver){
    vec v(4);
    v.zeros();

    vec result(4);
    result(0)=1.0;
    result(1)=0.414785;
    result(2)=0.130049;
    result(3)=0.0;

Implicit my;
my.Implicit_Scheme(v, 0.3, 4, 4, 1.0/3.0);

vec dif=abs(result-v);

CHECK(dif.max() < 1e-3);}


TEST(Crank_Nicolson_solver){

    vec v(4);
    v.zeros();

    vec result(4);
    result(0)=1.0;
    result(1)=0.412601;
    result(2)=0.112343;
    result(3)=0.0;

Crank_Nicolson my;
my.Crank_Nicolson_Scheme(v, 0.3, 4, 4, 1.0/3.0);

vec dif=abs(result-v);

CHECK(dif.max() < 1e-3);}



TEST(positiv_move){

    double first[] = {0.0, 0.0, 0.0, 0.2, 0.3, 0.6, 0.8, 0.7, 1.0, 1.0, 0.5};
    vector<double> v(first, first + sizeof(first) / sizeof(double) );

    double second[] = {0.8, 0.7, 0.5, 0.5, 0.5, 0.0, 0.0, 0.0};
    vector<double> result(second, second + sizeof(second) / sizeof(double) );

Monte_Carlo solve;
double l=0.5;
for (int i=0; i<v.size(); i++)
{ solve.positive_move(v, i, l);}

vector <double> diff(result.size());
double sum=0;
for (int i=0; i<result.size(); i++)
{diff[i]=abs(result[i]-v[i]);
    sum+=diff[i];}

CHECK(sum < 1e-3);}



TEST(negativ_move){

    double first[] = {0.0, 0.0, 0.0, 0.2, 0.3, 0.6, 0.8, 0.7, 1.0, 1.0, 0.5};
    vector<double> v(first, first + sizeof(first) / sizeof(double) );

    double second[] = {0.5, 0.5, 0.2, 0.3, 0.1, 0.0, 0.0, 0.0};
    vector<double> result(second, second + sizeof(second) / sizeof(double) );

Monte_Carlo solve;
double l=0.5;
for (int i=0; i<v.size(); i++)
{ solve.negative_move(v, i, l);}

vector <double> diff(result.size());
double sum=0;
for (int i=0; i<result.size(); i++)
{diff[i]=abs(result[i]-v[i]);
    sum+=diff[i];}

CHECK(sum < 1e-3);}


int main()
{
return UnitTest::RunAllTests();
}

