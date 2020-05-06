//
// Created by shaharnik on 24/04/2020.
//

#pragma once
#include <complex>

namespace solver
{

    class RealVariable 
    {
        public:
            double a;
            double b;
            double c;

        //Constructor
        RealVariable(double a ,double b, double c);
        RealVariable(): a(0), b(1), c(0){}


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
        friend RealVariable operator^(const RealVariable&, const double);

        //equal
        friend RealVariable operator ==(const RealVariable&, const RealVariable&);
        friend RealVariable operator ==(const RealVariable&, const double);
        friend RealVariable operator ==(const double, const RealVariable&);

    }; // end real


    class ComplexVariable
    {
        public: 
           double _re;
           double _im;

        //Constructor
        ComplexVariable(double _re, double _im);
        ComplexVariable(): _re(0.0), _im(0.0){}

        ///Operators

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
        friend ComplexVariable operator^(const ComplexVariable&,const double);

        //equal
        friend ComplexVariable operator==(const ComplexVariable&,const ComplexVariable&);
        friend ComplexVariable operator==(const ComplexVariable&,const double);
        friend ComplexVariable operator==(const double,const ComplexVariable&);
        friend ComplexVariable operator==(const ComplexVariable&,const complex<double>);
        friend ComplexVariable operator==(const complex<double>,const ComplexVariable&);

    }; // end complex 

    double solve(const RealVariable&);
    double solve(const ComplexVariable&);
}
