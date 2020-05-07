//
// Created by shaharnik on 24/04/2020.
//

#include "solver.hpp"
#include <complex>
#include <iostream>
#include <exception>

using namespace std;
using namespace solver;


//------------------------------------------RealVariable------------------------------------------

//Operators

double solver::solve(const RealVariable& x) 
{
    //cout << "a.x:" << x.a << "x.b:" << x.b << " x.c: " << x.c <<endl;
    if(x.a == 0.0) // bx+c 
    {
        if(x.b == 0.0 && x.c != 0.0) 
            throw std::out_of_range ("There is no Real solution! ");
        else return x.c / -x.b; // x = b/-c 
    } 
    else if(x.a == 0.0)
        return abs(x.c / x.b); 
     
    else if(x.c < 0.0) 
    {
        double x1 = (-x.b + sqrt(x.b * x.b - 4 * x.a * x.c)) / (2.0 * x.a);
        return x1;
    } 
    else 
        throw runtime_error ("There is no Real solution! ");
}

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
RealVariable solver::operator^(const RealVariable& x, const int y) 
{
    if(y == 0) 
    {
        return RealVariable(0.0, 0.0, 1.0);
    }

    else if(y == 2) 
    {
        return RealVariable(1.0, 0.0, 0.0);
    }

    if(y == 2.0 && x.c != 0.0 && x.b != 0.0) return RealVariable(pow(x.b, y), x.b * x.c*y, pow(x.c, y));
    if(y == 2.0 && x.c==0.0 && x.b!=0.0) return RealVariable(pow(x.b, y), 0.0, 0.0);
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

//------------------------------------------ComplexVariable------------------------------------------

//Operations

complex<double> solver::solve(const ComplexVariable& x)
{
    complex<double> ans;
    if(x._a == complex<double>(0,0))
    {
        if(x._b == complex<double>(0,0) && x._c != complex<double>(0,0)) 
            throw std::out_of_range ("There is no solution! ");
        else return x._c / -x._b; // x = b/-c 
    } 

    else if(x._a == complex<double>(0,0))
        return (-x._c / x._b); 
    
    else
    {
        complex<double> ans = ( (complex<double> (-1) * x._b ) + sqrt( (x._b * x._b) - 
        (complex<double> (4, 0) * x._a * x._c)) ) / ( complex<double>(2, 0) * x._a );
        return ans;
    }     
}

//plus
ComplexVariable solver::operator+(const ComplexVariable& x1, const ComplexVariable& x2) //(a+bi) + (c+di) = (a+c) + (b+d)i
{
    return ComplexVariable(x1._a + x2._a, x1._b + x2._b, x1._c + x2._c);
}

ComplexVariable solver::operator+(const ComplexVariable& x1,const double k) // (a+bi) + k = (a+k) + bi
{
    return ComplexVariable(x1._a, x1._b, x1._c+k);
}

ComplexVariable solver::operator+(const double k,const ComplexVariable& x1) // k + (a+bi) = (k+a) + bi 
{
    return ComplexVariable(x1._a, x1._b, x1._c+k);
}

ComplexVariable solver::operator+(const complex<double> k,const ComplexVariable& x1) // k + (a+bi) = (k+a) + bi 
{
    return ComplexVariable(x1._a, x1._b, x1._c+k);
}

ComplexVariable solver::operator+(const ComplexVariable& x1,const complex<double> k) // k + (a+bi) = (k+a) + bi 
{
    return ComplexVariable(x1._a, x1._b, x1._c+k);
}


//Subtraction
ComplexVariable solver::operator-(const ComplexVariable& x1,const ComplexVariable& x2)//(a+bi) - (c+di) = (a-c) + (b-d)i
{
    return ComplexVariable(x1._a - x2._a, x1._b - x2._b, x1._c - x2._c);
}

ComplexVariable solver::operator-(const ComplexVariable& x1,const double k)// (a+bi) - k = (a-k) + bi 
{
    return ComplexVariable(x1._a, x1._b, x1._c-k);
}

ComplexVariable solver::operator-(const double k,const ComplexVariable& x1) // k - (a+bi) = (k-a) + bi
{
    return ComplexVariable( -x1._a, -x1._b, k - x1._c);
}

ComplexVariable solver::operator-(const ComplexVariable& x1,const complex<double> k)// (a+bi) - k = (a-k) + bi 
{
    return ComplexVariable(x1._a, x1._b, x1._c-k);
}

ComplexVariable solver::operator-(const complex<double> k,const ComplexVariable& x1) // k - (a+bi) = (k-a) + bi
{
    return ComplexVariable( -x1._a, -x1._b, k - x1._c);
}


//multiplication
ComplexVariable solver::operator*(const ComplexVariable& x1,const ComplexVariable& x2) // (a+bi)*(c+di) = (ac-bd)+(ad+cb)i
{
    complex<double> new_a = x1._a * x2._c + x1._b * x2._b + x1._c * x2._a;
    complex<double> new_b = x1._b * x2._c + x1._c * x2._b; 
    complex<double> new_c = x1._c * x2._c;
    return ComplexVariable(new_a, new_b, new_c);
}

ComplexVariable solver::operator*(const ComplexVariable& x1, const double k) // (a+bi)*k = ak + bki
{
    return ComplexVariable(x1._a * k, x1._b * k, x1._c * k);
}

ComplexVariable solver::operator*(const double k,const ComplexVariable& x1)// k*(a+bi) = ka + kbi
{
    return ComplexVariable(x1._a * k, x1._b * k, x1._c * k);
}

ComplexVariable solver::operator*(const ComplexVariable& x1, const complex<double> k) 
{
    return ComplexVariable(x1._a * k, x1._b * k, x1._c * k);
}

ComplexVariable solver::operator*(const complex<double> k,const ComplexVariable& x1)
{
    return ComplexVariable(x1._a * k, x1._b * k, x1._c * k);
}


//division
ComplexVariable solver::operator/(const ComplexVariable& x1,const ComplexVariable& x2) // (a+bi)/(c+di) = (ac+bd)/(c^2+d^2)+(-ad+bc)/(c^2+d^2)
{
    complex<double> new_a = x1._a,new_b = x1._b, new_c = x1._c;
    if(x1._a != complex<double>(0,0) && x2._a != complex<double>(0,0)){
        new_a = 0;
        new_c = x1._c + x1._a/x2._a;
         
    }

    if(x1._a != complex<double>(0) && x2._b != complex<double>(0)){
       new_a = 0;
       new_b = x1._b + x1._a/x2._b;
    }

    if(x1._b != complex<double>(0) && x2._b != complex<double>(0)){
       new_c = x1._c + x1._b / x2._b; 
       new_b = 0;
    }
    return ComplexVariable(new_a, new_b, new_c);
}

ComplexVariable solver::operator/(const ComplexVariable& x,const double k) //(a+bi)/k = a/k+bi/k
{
    if(k == 0.0)
        throw runtime_error("invalid");
    else 
        return ComplexVariable(x._a/k,x._b/k, x._c/k);
}

ComplexVariable solver::operator/(const double k,const ComplexVariable& x) // k/(a+bi) = k/a + k/bi 
{
    return ComplexVariable(k/x._a,k/x._b, k/x._c);
}

ComplexVariable solver::operator/(const ComplexVariable& x,const complex<double> k) 
{
        return ComplexVariable(x._a/k,x._b/k, x._c/k);
}

ComplexVariable solver::operator/(const complex<double> k,const ComplexVariable& x)  
{
    return ComplexVariable(k/x._a,k/x._b, k/x._c);
}


//power
ComplexVariable solver::operator^(const ComplexVariable& x,const int k) // (a+bi)^k
{
    complex<double> new_a,new_b,new_c;
    if(k == 0) 
        return ComplexVariable(0, 0, 1);
    
    if(k == 1) 
        return ComplexVariable(x._a,x._b,x._c);
    
    if(k == 2)
        return ComplexVariable(x._a + (x._b * x._b), 0, x._c);
    
    else throw runtime_error("invalid"); 
      
}


//equal 
ComplexVariable solver::operator==(const ComplexVariable& x1 ,const ComplexVariable& x2) //(a+bi) = (c+di) ==> (a-c) - (b-d)i = 0  
{
    return ComplexVariable(x1._a - x2._a, x1._b - x2._b, x1._c - x2._c);
}

ComplexVariable solver::operator==(const ComplexVariable& x,const double k) // (a+bi) == k ==> (a-k)+bi == 0
{
    complex<double> new_c =  x._c + complex<double>(-1)*k; 
    return ComplexVariable(x._a,x._b, new_c);
}

ComplexVariable solver::operator==(const double k,const ComplexVariable& x) // k == (a+bi) ==> (k-a)-bi == 0
{
    complex<double> new_c =  x._c + complex<double>(-1)*k; 
    return ComplexVariable(x._a,x._b, new_c);
}

ComplexVariable solver::operator==(const ComplexVariable& x,const complex<double> k) 
{
    complex<double> new_c =  x._c + complex<double>(-1)*k; 
    return ComplexVariable(x._a,x._b, new_c);
}

ComplexVariable solver::operator==(const complex<double> k,const ComplexVariable& x) 
{
    complex<double> new_c =  x._c + complex<double>(-1)*k; 
    return ComplexVariable(x._a,x._b, new_c);
}

//-------------------------------------------------------------------------------------------------

ostream &operator<<(ostream &o, complex<double> x) 
{
    if(x.imag() >= 0)
        return (o << x.real() << "+" << x.imag() << "i");
     
    else 
        return (o << x.real() << x.imag() << "i");
}
