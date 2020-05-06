//
// Created by shaharnik on 24/04/2020.
//

#include "solver.hpp"
#include <complex>
#include <iostream>

using namespace std;
using namespace solver;

//Constructor
RealVariable::RealVariable(double a ,double b, double c)
{
    this->a = a;
    this->b = b;
    this->c = c;
}

//Operators

//plus
RealVariable solver::operator+(const RealVariable& x1, const RealVariable& x2) // (a,b,c)+(d,e,f) = (a+d,b+e,c+f)
{ 
    return RealVariable(x1.a + x2.a, x1.b + x2.b, x1.c + x2.c);
}

RealVariable solver::operator+(const RealVariable& x, const double& y) 
{
    return RealVariable(x.a, x.b, x.c + y);
}

RealVariable solver::operator+(const double& y, const RealVariable& x) 
{
    return x+y;
}


//subtraction 
RealVariable solver::operator-(const RealVariable &x, const RealVariable& y) 
{
    return RealVariable(x.a - y.a, x.b - y.b, x.c - y.c);
}

RealVariable solver::operator- (const RealVariable & x, const double y)
{
    return RealVariable(x.a, x.b, x.c - y);
}

RealVariable solver::operator-(const double& y, const RealVariable& x)
{
    return RealVariable(x.a, x.b, y - x.c);
}


//multiplication
RealVariable solver::operator*(const RealVariable& x, const RealVariable& y)
{
    return RealVariable(x.a * y.a, x.b * y.b, x.c * y.c);
}

RealVariable solver::operator*(const RealVariable& x, const double y)
{
    return RealVariable(x.a, x.b, x.c*y);
}

RealVariable solver::operator*(const double y, const RealVariable& x)
{
    return RealVariable(x.a*y, x.b*y, x.c*y);
}

//division
RealVariable solver::operator/ (const RealVariable& x, const double y)
{
    return RealVariable(x.a / y, x.b / y, x.c / y);
}

//power 
RealVariable solver::operator^(const RealVariable& x, const double y) 
{
    if(y == 0) 
    {
        return RealVariable(0, 0, 1);
    }

    else if(y == 2) 
    {
        return RealVariable(1, 0, 0);
    }

    if(y==2 && x.c != 0 &&x.b != 0) return RealVariable(pow(x.b, y), x.b * x.c*y, pow(x.c, y));
    if(y==2 && x.c==0 && x.b!=0) return RealVariable(pow(x.b, y), 0, 0);
}

//equal
RealVariable solver::operator==(const RealVariable& x, const RealVariable& y)
{
    return x - y;
}

RealVariable solver::operator==(const RealVariable& x, const double y)
{
    return x - y;
}

RealVariable solver::operator==(const double y, const RealVariable& x)
{
    return y - x;
}


//-------------------------------------------------------------------------------------------------


//Constractor

ComplexVariable::ComplexVariable(double _re, double _im)
{
    this->_re = _re;
    this->_im = _im;
}

//Operations

//plus
ComplexVariable operator+(const ComplexVariable& x1, const ComplexVariable& x2) //(a+bi) + (c+di) = (a+c) + (b+d)i
{
    return ComplexVariable(x1._re + x2._re, x1._im + x2._im);
}

ComplexVariable operator+(const ComplexVariable& x1,const double k) // (a+bi) + k = (a+k) + bi
{
    return ComplexVariable(x1._re + k, x1._im);
}

ComplexVariable operator+(const double k,const ComplexVariable& x1) // k + (a+bi) = (k+a) + bi 
{
    return ComplexVariable(k+x1._re, x1._im);
}

//Subtraction
ComplexVariable operator-(const ComplexVariable& x1,const ComplexVariable& x2)//(a+bi) - (c+di) = (a-c) + (b-d)i
{
    return ComplexVariable(x1._re - x2._re, x1._im - x2._im);
}

ComplexVariable operator-(const ComplexVariable& x1,const double k)// (a+bi) - k = (a-k) + bi 
{
    return ComplexVariable(x1._re - k, x1._im);
}

ComplexVariable operator-(const double k,const ComplexVariable& x1) // k - (a+bi) = (k-a) + bi
{
    return ComplexVariable(k-x1._re, x1._im);
}


//multiplication
ComplexVariable operator*(const ComplexVariable& x1,const ComplexVariable& x2) // (a+bi)*(c+di) = (ac-bd)+(ad+cb)i
{
    double new_re = x1._re * x2._re - x1._im * x2._im;
    double new_im = x1._re * x2._im + x1._im * x2._re;
    return ComplexVariable(new_re,new_im);
}

ComplexVariable operator*(const ComplexVariable& x1, const double k) // (a+bi)*k = ak + bki
{
    return ComplexVariable(x1._re * k, x1._im * k);
}

ComplexVariable operator*(const double k,const ComplexVariable& x1)// k*(a+bi) = ka + kbi
{
    return ComplexVariable(x1._re * k, x1._im * k);
}

//division
ComplexVariable operator/(const ComplexVariable& x1,const ComplexVariable& x2) // (a+bi)/(c+di) = (ac+bd)/(c^2+d^2)+(-ad+bc)/(c^2+d^2)
{
    double new_re = x1._re * x2._re + x1._im * x2._im;
    double new_im = -x1._re * x2._im + x1._im * x1._re;
    double d = x2._re * x2._re + x2._im * x2._im;
    return ComplexVariable((new_re)/(d),(new_im)/(d));
}

ComplexVariable operator/(const ComplexVariable& x,const double k) //(a+bi)/k = a/k+bi/k
{
    return ComplexVariable(x._re/k,x._im/k);
}

ComplexVariable operator/(const double k,const ComplexVariable& x) // k/(a+bi) = k/a + k/bi 
{
    return ComplexVariable(k/x._re,k/x._im);
}

//power
ComplexVariable operator^(const ComplexVariable& x,const double k) // (a+bi)^k
{
    if(k == 0) //k = 0
    {
        return ComplexVariable(1.0,0.0);
    }
    else if(k == 1)
    {
       return ComplexVariable(x._re,x._im);
    }
    else if(k == 2) // k=2 ==> (a^2-b^2)+2abi
    {
        double new_re = x._re * x._re - x._im * x._im;
        double new_im = 2 * x._re * x._im;
        return ComplexVariable(new_re,new_im);
    }
}


//equal 
ComplexVariable operator==(const ComplexVariable& x1 ,const ComplexVariable& x2) //(a+bi) = (c+di) ==> (a-c) - (b-d)i = 0  
{
    return ComplexVariable(x1._re - x2._re, x1._im - x2._im);
}

ComplexVariable operator==(const ComplexVariable& x,const double k) // (a+bi) == k ==> (a-k)+bi == 0
{
    return ComplexVariable(x._re-k,x._im);
}

ComplexVariable operator==(const double k,const ComplexVariable& x) // k == (a+bi) ==> (k-a)-bi == 0
{
    return ComplexVariable(k-x._re, -x._im);
}

//-------------------------------------------------------------------------------------------------

double solver::solve(const RealVariable& x) 
{
    cout<<"a.x:"<<x.a<<" x.b:"<<x.b<<" x.c: "<<x.c<<endl;
    if(x.a == 0) 
    {
        if(x.b == 0 && x.c!=0)
            throw std::out_of_range {"There is no Real solution! "};
        else return x.c/-x.b;
    } 
    else if(x.a == 0) 
    {
        return abs(x.c / x.b); 
    } 
    else if(x.c < 0) 
    {
        double x1 = (-x.b + sqrt(x.b * x.b - 4 * x.a * x.c)) / (2 * x.a);
        return x1;
    } 
    else 
    {
        throw runtime_error {"There is no Real solution! "};
    }


    // double solver::solve(const ComplexVariable& x){

    // } 
}
