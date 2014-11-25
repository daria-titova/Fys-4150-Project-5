#include <unittest++/UnitTest++.h>
#include <tridiag.h>
#include <Explicit.h>
#include <Implicit.h>
#include <Cra-Nic.h>
#include <iostream>

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



int main()
{
return UnitTest::RunAllTests();
}

