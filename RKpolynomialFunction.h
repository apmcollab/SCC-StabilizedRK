//
// RKsteadyStateCoeff.cpp 
//
// Author: Chris Anderson  
// (C) UCLA 2012 
//
//
// ############################################################
// RKsteadyStateCoeff.cpp 
//
// Version : May 26, 2004 
//
// Chris Anderson UCLA 
//###############################################################
//

#include <iostream>
#include <iomanip>

#ifndef RK_POLYNOMIAL_FUNCTION_
#define RK_POLYNOMIAL_FUNCTION_

class RKpolynomialFunction
{

public:

RKpolynomialFunction()
{
    degree = -1;
    a      =  nullptr;
}

RKpolynomialFunction(const RKpolynomialFunction& P)
{
    degree = P.degree;                             // copy degree
    a = new double[degree + 1];                    // allocate space for new coefficients
    for(int i=0; i < degree+1; i++) a[i] = P.a[i]; // copy over coefficients
}


RKpolynomialFunction(int deg)
{
    degree = deg;
    a      = new double[degree+1];   // create space for coefficients
    int i;
    for(i = 0; i < degree+1; i++)
    {
    a[i] = 0.0;
    }
}

RKpolynomialFunction(int n, double coefficients[])
{
    degree = n;
    a      = new double[degree+1];   // create space for coefficients
    int i;
    for(i = 0; i < degree+1; i++)
    {
    a[i] = coefficients[i];          // assign values
    }
}

~RKpolynomialFunction()
{
    if(a != nullptr) delete [] a;   // remove space for coefficients
}

RKpolynomialFunction differentiate()
{
    int i;
    long rDegree;
    
    if(degree > 0)
    {rDegree = degree - 1;}
    else
    {rDegree = 0;}

    RKpolynomialFunction R(rDegree);

    if(degree > 0)
    {
    for(i = 0; i <= degree - 1; i++)
    {R.a[i] =  a[i+1]*double(i+1);}
    }
    else
    {
    R.a[0] = 0.0;
    }

    return R;
}

RKpolynomialFunction operator=(const RKpolynomialFunction& P)
{
    if(this == &P) return *this;                   // return if already =
    if(a != nullptr) delete[] a;                         // remove "this" coefficients
    degree = P.degree;                             // copy degree
    a = new double[degree + 1];                    // allocate space for new coefficients
    for(int i=0; i < degree+1; i++) a[i] = P.a[i]; // copy over coefficients
    return *this;
}

RKpolynomialFunction operator+(const RKpolynomialFunction& P)
{
    RKpolynomialFunction  R;
    int i;

    if(degree >= P.degree)
    {
      R = *this;
      for(i = 0; i <= P.degree; i++)
      {R.a[i] += P.a[i];}
    }
    else
    {
       R = P;
       for(i = 0; i <= degree; i++)
       {R.a[i] += a[i];}
    }
    return R;
}

RKpolynomialFunction operator-(const RKpolynomialFunction& P)
{
    RKpolynomialFunction  R;
    int i;

    if(degree >= P.degree)
    {
      R = *this;
      for(i = 0; i <= P.degree; i++)
      {R.a[i] -= P.a[i];}
    }
    else
    {
       R = -P;
       for(i = 0; i <= degree; i++)
       {R.a[i] += a[i];}
    }
    return R;
}
RKpolynomialFunction operator*(const RKpolynomialFunction& P)
{
    int i; int k;
    RKpolynomialFunction R(degree+P.degree);

    for(i = 0; i <= degree; i++)
    {
    for(k = i; k <= i+P.degree; k++)
    {R.a[k] += a[i]*P.a[k-i];}
    }

    return R;
}

RKpolynomialFunction operator-() const
{
    RKpolynomialFunction  R(*this);
    int i;
    for(i = 0; i <= degree; i++)
    {R.a[i] = -a[i];}
    return R;
}

double operator()(double x)
{
//  Evaluate using Horner's method
//
    double result = 0.0;
    int i;
    for(i = degree; i>=0; i--)
    {
    result = x*result + a[i];
    }
    return result;
}

const double&  operator[](long i) const
{
   if(i > degree)
   {return zero;}

   return *(a+i);
}

double&  operator[](long i)
{
   if(i > degree)
   {
   resizeTo(i);
   }
    return *(a+i);
}

RKpolynomialFunction operator*(double alpha)
{
    RKpolynomialFunction R(*this);
   long i;
   for(i =0; i <= degree; i++) R[i] = alpha*R[i];
   return R;

}
RKpolynomialFunction operator/(double alpha)
{
    RKpolynomialFunction R(*this);
   long i;
   for(i =0; i <= degree; i++) R[i] = R[i]/alpha;
   return R;

}
friend RKpolynomialFunction operator*(double alpha, const RKpolynomialFunction& P)
{
   RKpolynomialFunction R(P);
   long i;
   for(i =0; i <= P.degree; i++) R[i] = alpha*R[i];
   return R;
}

RKpolynomialFunction shift(double p)
{
     RKpolynomialFunction R(degree);

     RKpolynomialFunction S(1);
     S[0] = -p; S[1] = 1.0;      // S = (x - p)
     RKpolynomialFunction Q(S);

     R[0] = a[0];
     long i;
     for(i = 1; i <= degree; i++)
     {
     R = R + a[i]*Q;
     Q = Q*S;
     }
     return R;
}

RKpolynomialFunction scale(double alpha)
{
     RKpolynomialFunction R(*this);
     long i;

     double s = 1.0;
     for(i = 1; i <= R.degree; i++)
     {
     s = s*alpha;
     R.a[i] = R.a[i]*s;
     }
     return R;
}

friend std::ostream& operator <<(std::ostream& outStream, const RKpolynomialFunction& P)
{
    outStream.setf(std::ios::fixed);
    outStream.precision(2);
    if(P.degree >= 0)
    {outStream << std::setw(2) << P.a[0];}
    int i;
    if(P.degree >= 1)
    {
    for(i = 1; i <= P.degree; i++)
    outStream << " + " << std::setw(2) << P.a[i] << "x^";
    outStream.setf(std::ios::left);
    outStream << std::setw(2) << i;
    outStream.setf(std::ios::right);
    }
    return outStream;
}
//
//    Initialize
//
void initialize()
{
    if(a != nullptr) delete [] a;
    degree = -1;
    a      =  nullptr;
}

void initialize(const RKpolynomialFunction& P)
{
    if(a != nullptr) delete [] a;
    degree = P.degree;                             // copy degree
    a = new double[degree + 1];                    // allocate space for new coefficients
    for(int i=0; i < degree+1; i++) a[i] = P.a[i]; // copy over coefficients
}

void initialize(int deg)
{
    if(a != nullptr) delete [] a;
    degree = deg;
    a      = new double[degree+1];
    int i;
    for(i = 0; i < degree+1; i++)
    {
    a[i] =0.0;
    }
}

void initialize(int n, double coefficients[])
{
    if(a != nullptr) delete [] a;
    degree = n;
    a      = new double[degree+1];
    int i;
    for(i = 0; i < degree+1; i++)
    {
    a[i] = coefficients[i];
    }
}

void resizeTo(int n)
{
    double* newCoefficients = new double[n+1];
    long i;
    if(n <= degree)
    {
    for(i = 0; i <= n; i++)newCoefficients[i] = a[i];
    }
    else
    {
    for(i = 0; i <= degree; i++) newCoefficients[i] = a[i];
    for(i = degree + 1; i <= n; i++) newCoefficients[i] = 0.0;
    }
    if(a != nullptr) delete [] a;
    a       = newCoefficients;
    degree = n;
}

int getDegree() const
{
     long i;
     int deg = 0;
     for(i = degree; i >= 0; i--)
     {
     if((deg == 0)&&(a[i] != 0.0)) deg = i;
     }
     return deg;
}

double* getCoefficients() const
{
     return a;
}

    public  :

    double*           a;           // ptr to coefficient array
    int               degree;

    constexpr static double zero = 0.0;
};

#endif

