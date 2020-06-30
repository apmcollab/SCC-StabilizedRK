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

#include "RKsteadyStateCoeff.h"

#include <iostream>
#include <iomanip>
using namespace std;

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void RKsteadyStateCoeff::getTchebyShiftFactors(double M, double gamma,
double& delta, double& beta, RKpolynomialFunction& ChebyPoly)
{
    RKpolynomialFunction P = ChebyPoly;
    
    RKpolynomialFunction Pprime;
    Pprime = P.differentiate();

    PprimeCondition Pdelta;
    Pdelta.initialize(P,Pprime,M);

    double a = 0.0; 
    double b;
    
    double tolerance = 1.0e-12;
    
    if(fabs(gamma-2.0) > tolerance) 
    {
    if(gamma <  1.0) b = 2.0/gamma;
    if(gamma >= 1.0) b = gamma*2.0;
    if(gamma >    M) b = M-M/10.0;
   
	delta = Pdelta.getRoot(&a, &b,&tolerance);
    }
    else
    {
    delta = 0.0;
    }
    

    P = P.scale(2.0/(M-delta));
    P = P.shift(-((2.0*M)/(M-delta) - 1)*((M-delta)/2.0));
    beta = P(0.0);
}

//
// Given an array of RKpolynomials of length maxDegree+1, this routine
// computes the first maxDegree + 1 Chebyshev polynomials
// using the recurrance 
//
// P_0     = 1
// P_1     = x
// P_(n+1) = (2*x)*P_(n) - P_(n-1)
//
void createChebyshevPolynomials(RKpolynomialFunction* P, long maxDegree)
{
   long i;

   RKpolynomialFunction Xpoly(1);
   Xpoly[0] = 0.0; Xpoly[1] = 1.0;

   P[0].initialize(0);
   P[0][0] = 1.0;

   P[1].initialize(1);
   P[1][0] = 0.0; P[1][1] = 1.0;

   for(i = 1; i < maxDegree; i++)
   {
   P[i+1] =2.0*Xpoly*P[i] - P[i-1];
   }
}

RKcoefficients::RKcoefficients()
{
    alphaCoeff = 0;
    stageOrder = 0;
    gamma      = 0;
}

RKcoefficients::RKcoefficients(long stageOrder, double gamma)
{ 
    this->stageOrder = stageOrder;
    this->gamma      = gamma;

    long i;
    alphaCoeff   = new double*[stageOrder];
    for(i = 0; i < stageOrder; i++) alphaCoeff[i] = new double[stageOrder];

    RKsteadyStateCoeff::getRKcoefficients(stageOrder, gamma, alphaCoeff);
}

RKcoefficients::~RKcoefficients()
{
    long i;

    if(alphaCoeff != 0) 
    {
    for(i = 0; i < stageOrder; i++) delete [] alphaCoeff[i];
    delete [] alphaCoeff;
    }
}

//
// RKsteadyStateCoeff Static initializer
//

long             RKsteadyStateCoeff::cacheStorageSize       = 100;
RKcoefficients** RKsteadyStateCoeff::coeffCache             = new RKcoefficients*[100];
long             RKsteadyStateCoeff::cacheSize              = 0;
long             RKsteadyStateCoeff::cacheStorageIncrement  = 100;


double** RKsteadyStateCoeff::getRKcoefficientsPtr(long stageOrder, 
double gamma)
{
    long i;
    //
    // Check cache for coefficients
    //
    for(i = 0; i < cacheSize; i++)
    {
    if((stageOrder == coeffCache[i]->stageOrder)
     &&(gamma      == coeffCache[i]->gamma))
    {return coeffCache[i]->alphaCoeff;}
    }
    //
    // Create new cache entry and 
    // then return the coefficients
    //
    if(cacheSize + 1 > cacheStorageSize) expandCache(100);

    coeffCache[cacheSize] = new RKcoefficients(stageOrder,gamma);
    cacheSize++;

    /*
    long j;
    for(j = 0; j < stageOrder; j++)
    {
    for(i = 0; i < stageOrder; i++)
    {
    printf("%10.7e  ",(coeffCache[cacheSize-1]->alphaCoeff)[i][j]);
    }
    printf("\n \n");
    }
    */
    return coeffCache[cacheSize-1]->alphaCoeff;
}

void RKsteadyStateCoeff::expandCache(long storageIncrement)
{
    long i;

    // allocate new cache pointer array 

    long newSize              = cacheStorageSize + storageIncrement;
    RKcoefficients** cachePtr = new RKcoefficients*[newSize];

    // copy over elements
    for(i = 0; i < cacheSize; i++) 
    {cachePtr[i] = coeffCache[i];}

    // remove existing elements
    delete [] coeffCache;

    // assign new pointer to old pointer
    coeffCache      = cachePtr;

    cacheStorageSize = newSize;
}



void RKsteadyStateCoeff::getRKcoefficients(long stageOrder, double gamma,
double** alphaCoefficients) 
{
    long i;  long j; long k; 

    double delta;
    double beta;
    double M;

    double* coeffPtr;
    double* deltaK = new double[stageOrder];
    double* MK     = new double[stageOrder];
    double* betaK  = new double[stageOrder];

    long degree = stageOrder;
//
//  Create Chebyshev polynomials
//
    RKpolynomialFunction* T = new RKpolynomialFunction[degree+1];
    createChebyshevPolynomials(T,degree);
//
// Step #1 compute shift factors 
//
    deltaK[0] =    0.0;
    MK[0]     =   gamma;

    double deg;

    for(k = 2; k <= degree; k++)
    {
    deg = k;
    M   = gamma*deg*deg;

    getTchebyShiftFactors(M, gamma,delta, beta,T[k]);
    MK[k-1]     = gamma*deg*deg;
    deltaK[k-1] = delta;
    betaK[k-1]  = beta;
    }
//
//  Create the coefficients matrix
//
    RKpolynomialFunction P;

    double** alphaMatrix = 0;
    double** alphaRHS    = 0;

    create2Darray(alphaMatrix,degree,degree);
    create2Darray(alphaRHS,degree,degree);

    for(i = 0; i < degree; i++)
    {
    for(j = 0; j < degree; j++)
    {
		alphaCoefficients[i][j] = 0.0;
	}}


    for(k = 0; k < degree; k++)
    {
    alphaMatrix[0][k] = 1.0;
    }

    if(degree == 1)
    {
    alphaRHS[0][0]   = gamma/MK[degree-1];
    }
    else
    {
    alphaMatrix[1][1] = gamma/MK[degree-1];
    alphaRHS[0][0]    = alphaMatrix[1][1];
    }


    double mFactor;
    double mFactorProd;

    for(k = 2; k <= degree-1; k++)
    {

    P      = T[k];
    M      = MK[k-1];
    delta  = deltaK[k-1];
    beta   = betaK[k-1];

    P = P.scale(2.0/(M-delta));
    P = P.shift(-((2.0*M)/(M-delta) - 1)*((M-delta)/2.0));
    P = P/beta;

    coeffPtr = P.getCoefficients();

    mFactor     = MK[k-1]/MK[degree-1];

    mFactorProd = mFactor;
    for(i = 1; i <= k; i++)
    {
    alphaMatrix[i][k]   = coeffPtr[i]*mFactorProd;
    alphaRHS[i-1][k-1]  = coeffPtr[i]*mFactorProd;
    mFactorProd *= mFactor;
    }

    }
//
//  Final column of right hand side matrix
//
    if(degree > 1)
    {
    k = degree;

    P      = T[k];
    M      = MK[k-1];
    delta  = deltaK[k-1];
    beta   = betaK[k-1];

    P = P.scale(2.0/(M-delta));
    P = P.shift(-((2.0*M)/(M-delta) - 1)*((M-delta)/2.0));
    P = P/beta;

    coeffPtr = P.getCoefficients();
    mFactor  = MK[k-1]/MK[degree-1];

    mFactorProd = mFactor;
    for(i = 1; i <= k; i++)
    {
    alphaRHS[i-1][degree-1] = coeffPtr[i]*mFactorProd;
    mFactorProd *= mFactor;
    }
    }
//
//
//  Create coefficients by backsolving

    double*  c = new double[degree];

    for(k = 1; k <= degree; k++)
    {

    for(i = 0; i < k; i++)
    {
    c[i] = alphaRHS[i][k-1];
    }

    backSolve(k, c, alphaMatrix);

    for(i = 0; i < k; i++)
    {
    alphaCoefficients[i][k-1] = c[i];
    }

    }
//
//   clean up 
//
	delete [] deltaK;
    delete [] MK;
    delete [] betaK;
    delete [] c;

    delete2Darray(alphaMatrix);
    delete2Darray(alphaRHS);

    delete [] T;
}

void RKsteadyStateCoeff::backSolve(long N, double* c, double** A)
{
    long i; long j;
    double     sum;

    c[N-1] = c[N-1]/A[N-1][N-1];

    for(i = N-2; i >= 0; i--)
    {
    sum = 0.0;
    for(j = N-1; j >= i+1; j--)
    {
    sum += A[i][j]*c[j];
    }
    c[i] = (c[i] - sum)/A[i][i];
    }
}

void RKsteadyStateCoeff::create2Darray(double** &a, long m, long n)
{
    double* dataPtr = new double[m*n];
	a               = new double*[m];

    long i;

    for(i = 0; i < m; i++) 
    {
    a[i] =  &(dataPtr[i*n]); 
    }

    for(i = 0; i < m*n; i++) 
    {
    dataPtr[i] = 0.0;
    }
}

void RKsteadyStateCoeff::delete2Darray(double** &a)
{
	double* dataPtr = a[0];
    delete [] dataPtr;
    delete [] a;
}



#include <stdlib.h>
#include <float.h> 
#include <math.h> 
/* 
   d1machForZeroin is a localized version of d1mach -- d1mach
   translated by f2c (version 19980403) and then updated to
   C++ syntax. 
*/
double PprimeCondition::d1machForZeroin_(long *i) 
{ 
	switch(*i){ 
 	  case 1: return DBL_MIN; 
 	  case 2: return DBL_MAX; 
	  case 3: return DBL_EPSILON/FLT_RADIX; 
 	  case 4: return DBL_EPSILON; 
 	  case 5: return log10(double(FLT_RADIX)); 
	  } 
 	fprintf(stderr, "invalid argument: d1mach(%ld)\n", *i); 
 	exit(1); return 0; /*/+ for compilers that complain of missing return values +/ */
} 
//
//  getRoot is a localized version of Zeroin.f -- 
//  zeroin.f translated by f2c (version 19980403)and then 
//  updated to C++ syntax. 
//
//  Original Comments :
//
//  zeroin : a zero of the function  f(x) is computed in the interval ax,bx.  
//
//  Input..  
//
//  ax     left endpoint of initial interval  
//  bx     right endpoint of initial interval  
//  f      function subprogram which evaluates f(x) for any x in  
//         the interval  ax,bx  
//  tol    desired length of the interval of uncertainty of the  
//         final result (.ge.0.)  
//
//  Output..  
//
//  zeroin abscissa approximating a zero of  f  in the interval ax,bx  
//
//  It is assumed  that   f(ax)   and   f(bx)   have  opposite  signs  
//  this is checked, and an error message is printed if this is not  
//  satisfied.   zeroin  returns a zero  x  in the given interval  
//  ax,bx  to within a tolerance  4*macheps*abs(x)+tol, where macheps  is 
//  the  relative machine precision defined as the smallest representable 
//  number such that  1.+macheps .gt. 1.  
//
//  This function subprogram is a slightly  modified  translation  of  
//  the algol 60 procedure  zero  given in  richard brent, algorithms for 
//  minimization without derivatives, prentice-hall, inc. (1973).  
//
double PprimeCondition::getRoot(double* ax, double* bx, double *tol)
{

    static long c__4 = 4;
    double ret_val, d__1, d__2;
    static double a, b, c__, d__, e, p, q, r__, s;
    static double fa, fb, fc, xm, eps, tol1;

    eps = d1machForZeroin_(&c__4);
    tol1 = eps + 1.;

    a = *ax;
    b = *bx;
    fa = this->operator()(a);
    fb = this->operator()(b);
    
    if(fabs(fa) < *tol) return a;
    if(fabs(fb) < *tol) return b;
/*     check that f(ax) and f(bx) have different signs */
    if (fa == 0. || fb == 0.) {
	goto L20;
    }
    if (fa * (fb / fabs(fb)) <= 0.) {
	goto L20;
    }

    printf("##########  Root Finder Error  ############# \n");
    printf("The function values at the initial interval \n");
    printf("endpoionts are not of opposite sign. \n");
    printf("The root value returned is inaccurate.\n\n");

    ret_val = 0;
    return ret_val;
L20:
    c__ = a;
    fc = fa;
    d__ = b - a;
    e = d__;
L30:
    if (fabs(fc) >= fabs(fb)) {
	goto L40;
    }
    a = b;
    b = c__;
    c__ = a;
    fa = fb;
    fb = fc;
    fc = fa;
L40:
    tol1 = eps * 2. * fabs(b) + *tol * .5;
    xm = (c__ - b) * .5;
    if (fabs(xm) <= tol1 || fb == 0.) {
	goto L150;
    }

/* see if a bisection is forced */

    if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
	goto L50;
    }
    d__ = xm;
    e = d__;
    goto L110;
L50:
    s = fb / fa;
    if (a != c__) {
	goto L60;
    }

/* linear interpolation */

    p = xm * 2. * s;
    q = 1. - s;
    goto L70;

/* inverse quadratic interpolation */

L60:
    q = fa / fc;
    r__ = fb / fc;
    p = s * (xm * 2. * q * (q - r__) - (b - a) * (r__ - 1.));
    q = (q - 1.) * (r__ - 1.) * (s - 1.);
L70:
    if (p <= 0.) {
	goto L80;
    }
    q = -q;
    goto L90;
L80:
    p = -p;
L90:
    s = e;
    e = d__;
    if (p * 2. >= xm * 3. * q - (d__1 = tol1 * q, fabs(d__1)) || p >= (d__2 = 
	    s * .5 * q, fabs(d__2))) {
	goto L100;
    }
    d__ = p / q;
    goto L110;
L100:
    d__ = xm;
    e = d__;
L110:
    a = b;
    fa = fb;
    if (fabs(d__) <= tol1) {
	goto L120;
    }
    b += d__;
    goto L140;
L120:
    if (xm <= 0.) {
	goto L130;
    }
    b += tol1;
    goto L140;
L130:
    b -= tol1;
L140:
    fb = this->operator()(b);
    if (fb * (fc / fabs(fc)) > 0.) {
	goto L20;
    }
    goto L30;
L150:
    ret_val = b;
    return ret_val;
} /* zeroin_ */


RKpolynomialFunction::RKpolynomialFunction()
{
    degree = -1;
    a      =  0;
}

RKpolynomialFunction::RKpolynomialFunction(const RKpolynomialFunction& P)
{
    degree = P.degree;                             // copy degree
    a = new double[degree + 1];                    // allocate space for new coefficients
    for(int i=0; i < degree+1; i++) a[i] = P.a[i]; // copy over coefficients
}


RKpolynomialFunction::RKpolynomialFunction(int deg)
{
    degree = deg;
    a      = new double[degree+1];   // create space for coefficients
    int i;
    for(i = 0; i < degree+1; i++)
    {
    a[i] = 0.0;
    }
}

RKpolynomialFunction::RKpolynomialFunction(int n, double coefficients[])
{
    degree = n;
    a      = new double[degree+1];   // create space for coefficients
    int i;
    for(i = 0; i < degree+1; i++)
    {
    a[i] = coefficients[i];          // assign values
    }
}

RKpolynomialFunction::~RKpolynomialFunction()
{
    if(a != 0) delete [] a;   // remove space for coefficients
}

RKpolynomialFunction RKpolynomialFunction::differentiate()
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

RKpolynomialFunction RKpolynomialFunction::operator=(const RKpolynomialFunction& P)
{
    if(this == &P) return *this;                   // return if already =
    if(a != 0) delete[] a;                         // remove "this" coefficients
    degree = P.degree;                             // copy degree
    a = new double[degree + 1];                    // allocate space for new coefficients
    for(int i=0; i < degree+1; i++) a[i] = P.a[i]; // copy over coefficients
    return *this;
}

RKpolynomialFunction RKpolynomialFunction::operator+(const RKpolynomialFunction& P)
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

RKpolynomialFunction RKpolynomialFunction::operator-(const RKpolynomialFunction& P)
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
RKpolynomialFunction RKpolynomialFunction::operator*(const RKpolynomialFunction& P)
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

RKpolynomialFunction RKpolynomialFunction::operator-() const
{
    RKpolynomialFunction  R(*this);
    int i;
    for(i = 0; i <= degree; i++)
    {R.a[i] = -a[i];}
    return R;
}

double RKpolynomialFunction::operator()(double x)
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

const double&  RKpolynomialFunction::operator[](long i) const
{
   if(i > degree)
   {return zero;}
	return *(a+i);
}

double&  RKpolynomialFunction::operator[](long i)
{
   if(i > degree)
   {
   resizeTo(i);
   }
	return *(a+i);
}
RKpolynomialFunction RKpolynomialFunction::operator*(double alpha)
{
	RKpolynomialFunction R(*this);
   long i;
   for(i =0; i <= degree; i++) R[i] = alpha*R[i];
   return R;

}
RKpolynomialFunction RKpolynomialFunction::operator/(double alpha)
{
	RKpolynomialFunction R(*this);
   long i;
   for(i =0; i <= degree; i++) R[i] = R[i]/alpha;
   return R;

}
RKpolynomialFunction operator*(double alpha, const RKpolynomialFunction& P)
{
  	RKpolynomialFunction R(P);
   long i;
   for(i =0; i <= P.degree; i++) R[i] = alpha*R[i];
   return R;
}

RKpolynomialFunction RKpolynomialFunction::shift(double p)
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

RKpolynomialFunction RKpolynomialFunction::scale(double alpha)
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

ostream& operator <<(ostream& outStream, const RKpolynomialFunction& P)
{
    outStream.setf(ios::fixed);
    outStream.precision(2);
    if(P.degree >= 0)
    {outStream << setw(2) << P.a[0];}
    int i;
    if(P.degree >= 1)
    {
    for(i = 1; i <= P.degree; i++)
    {
    outStream << " + " << setw(2) << P.a[i] << "x^";
    outStream.setf(ios::left);
    outStream << setw(2) << i;
    outStream.setf(ios::right);
    }
    }
    return outStream;
}
//
//    Initialize
//
void RKpolynomialFunction::initialize()
{
    if(a != 0) delete [] a;
    degree = -1;
    a      =  0;
}

void RKpolynomialFunction::initialize(const RKpolynomialFunction& P)
{
    if(a != 0) delete [] a;
    degree = P.degree;                             // copy degree
    a = new double[degree + 1];                    // allocate space for new coefficients
    for(int i=0; i < degree+1; i++) a[i] = P.a[i]; // copy over coefficients
}

void RKpolynomialFunction::initialize(int deg)
{
    if(a != 0) delete [] a;
    degree = deg;
    a      = new double[degree+1];
    int i;
    for(i = 0; i < degree+1; i++)
    {
    a[i] =0.0;
    }
}

void RKpolynomialFunction::initialize(int n, double coefficients[])
{
    if(a != 0) delete [] a;
    degree = n;
    a      = new double[degree+1];
    int i;
    for(i = 0; i < degree+1; i++)
    {
    a[i] = coefficients[i];
    }
}

void RKpolynomialFunction::resizeTo(int n)
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
    if(a != 0) delete [] a;
    a       = newCoefficients;
    degree = n;
}

int RKpolynomialFunction::getDegree() const
{
     long i;
     int deg = 0;
     for(i = degree; i >= 0; i--)
     {
     if((deg == 0)&&(a[i] != 0.0)) deg = i;
     }
     return deg;
}

double* RKpolynomialFunction::getCoefficients() const
{
     return a;
}

const double RKpolynomialFunction::zero = 0.0; 
