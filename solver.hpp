//
// Created by shaharnik on 24/04/2020.
//

#pragma once
#include <stdlib.h>
#include <iostream>
#include <complex>

using namespace std;

namespace solver
{

    class RealVariable 
    {
        public:
            double a,b,c;

        //Constructor
        RealVariable(): a(0), b(1), c(0){};
        RealVariable(double a ,double b, double c):a(a),b(b),c(c){};
        ~RealVariable();


        // RealVariable::RealVariable(double a ,double b, double c){
        //     this->a = a;
        //     this->b = b;
        //     this->c = c;
        // }

        //Operators

        //plus
        friend RealVariable operator+(const RealVariable&, const RealVariable&);
        friend RealVariable operator+(const RealVariable&, const double&);
        friend RealVariable operator+(const double&, const RealVariable&);

        //subtraction
        friend RealVariable operator-(const RealVariable&, const RealVariable&);
        friend RealVariable operator-(const RealVariable&, const double);
        friend RealVariable operator-(const double& y, const RealVariable& x);

        //multiplication
        friend RealVariable operator*(const RealVariable &, const RealVariable &);
        friend RealVariable operator*(const RealVariable &,const double);
        friend RealVariable operator*(const double, const RealVariable &);
        
        //division
        friend RealVariable operator/(const RealVariable&, const RealVariable&);
        friend RealVariable operator/(const RealVariable&, const double);
        friend RealVariable operator/(const double, const RealVariable&);

        //power
        friend RealVariable operator^(const RealVariable&, const int);

        //equal
        friend RealVariable operator ==(const RealVariable&, const RealVariable&);
        friend RealVariable operator ==(const RealVariable&, const double);
        friend RealVariable operator ==(const double, const RealVariable&);

    }; // end real


    class ComplexVariable
    {
        public:
            complex<double> _a,_b,_c;

    
            ComplexVariable(): _a(0),_b(1),_c(0){};
            ComplexVariable(const complex<double> _a, const complex<double> _b, const complex<double> _c): _a(_a), _b(_b), _c(_c){};
            ~ComplexVariable();

        // ComplexVariable::ComplexVariable(const complex<double> _a, const complex<double> _b, const complex<double> _c)
        // {
        //     this->_a;
        //     this->_b;
        //     this->_c;
        // }


        // ///Operators

        //plus 
        friend ComplexVariable operator+(const ComplexVariable&,const ComplexVariable&);
        friend ComplexVariable operator+(const ComplexVariable&,const double);
        friend ComplexVariable operator+(const double,const ComplexVariable&);
        friend ComplexVariable operator+(const ComplexVariable&,const complex<double>);
        friend ComplexVariable operator+(const complex<double>,const ComplexVariable&);


        //subtraction 
        friend ComplexVariable operator-(const ComplexVariable&,const ComplexVariable&);
        friend ComplexVariable operator-(const ComplexVariable&,const double);
        friend ComplexVariable operator-(const double,const ComplexVariable&);
        friend ComplexVariable operator-(const ComplexVariable&,const complex<double>);
        friend ComplexVariable operator-(const complex<double>,const ComplexVariable&);
      

        //multiplication
        friend ComplexVariable operator*(const ComplexVariable&,const ComplexVariable&);
        friend ComplexVariable operator*(const ComplexVariable&,const double);
        friend ComplexVariable operator*(const double,const ComplexVariable&);
        friend ComplexVariable operator*(const ComplexVariable&,const complex<double>);
        friend ComplexVariable operator*(const complex<double>,const ComplexVariable&);
    

        //division
        friend ComplexVariable operator/(const ComplexVariable&,const ComplexVariable&);
        friend ComplexVariable operator/(const ComplexVariable&,const double);
        friend ComplexVariable operator/(const double,const ComplexVariable&);
        friend ComplexVariable operator/(const ComplexVariable&,const complex<double>);
        friend ComplexVariable operator/(const complex<double>,const ComplexVariable&);
  

        //power 
        friend ComplexVariable operator^(const ComplexVariable&,const int);


        //equal
        friend ComplexVariable operator==(const ComplexVariable&,const ComplexVariable&);
        friend ComplexVariable operator==(const ComplexVariable&,const double);
        friend ComplexVariable operator==(const double,const ComplexVariable&);
        friend ComplexVariable operator==(const ComplexVariable&,const complex<double>);
        friend ComplexVariable operator==(const complex<double>,const ComplexVariable&);

    }; // end complex 

    double solve(const RealVariable&);
    complex<double> solve(const ComplexVariable&);
    ostream &operator<<(ostream&, complex<double>);

};
