//
// RKsteadyStateCoeff.h 
//
// Author: Chris Anderson  
// (C) UCLA 2012 
//
//
// ############################################################
// RKsteadyStateCoeff.h 
//
// Version : May 26, 2004 
//
// Chris Anderson UCLA 
//###############################################################
//

#include <iostream>
#include <iomanip>
#include <vector>

#include "RKpolynomialFunction.h"
#include "PprimeCondition.h"

#ifndef RK_STEADY_STATE_COEFF_
#define RK_STEADY_STATE_COEFF_



//
// class RKsteadyStateCoeff is a class for constructing coefficients for
// stabilized Runge-Kutta methods useful for evolving ODE systems to
// steady state.
//
// The determination of the coefficients is described in
// CAM report 04-26, Christopher R. Anderson and Christopher Elion,
// "Accelerated Solutions of Nonlinear Equations Using Stabilized Runge-Kutta Methods"
// April 2004 (ftp://ftp.math.ucla.edu/pub/camreport/cam04-26.pdf)
//
//
// stageOrder = the number of stages used in the Runga-Kutta method.
//              If the problem is not very stiff, then low order
//              (1-3) methods work best, while high order methods
//              (4-10) work better for mildly stiff, or stiff equations.
//
// RKgamma    = The factor that determines the size of the stability
//              region and the magnitude of the damping associated with of the
//              stabilized method. (0 < RKgamma < 2.0) The stability interval is
//              [-RKgamma*stageOrder*stageOrder, 0].
//
//              Larger RKgamma gives methods with larger stability regions, but
//              smaller damping; while smaller RKgamma gives smaller stability
//              regions with more damping.
//
//              Good values for RKgamma are 1.5 <= RKgamma <= 1.75.
//
//
// If the s is the stage order of the method then
// the coefficients of the RK method are assumed to be in an s X s array,
// with the  upper left (s-1)X(s-1) block being the coefficients used to
// determine the stages K[1] through K[s-1].
//
// The last column of the coefficient array stores the coefficients
// of the linear combination of the stages used to construct the solution.
//
//
//    | alpha[0][0] | alpha[0][1] | alpha[0][2] | .....   | alpha[0][s-2]   | alpha[0][s -1]  |
//    |     0       | alpha[1][1] | alpha[1][2] | .....   | alpha[1][s-2]   | alpha[1][s -1]  |
//    |     0       |      0      | alpha[2][2] | .....   | alpha[2][s-2]   | alpha[2][s -1]  |
//    |     0       |      0      |       0     |         |                 |                 |
//    |     0       |      0      |       0     |    0    | alpha[s-2][s-2] | alpha[s-2][s-1] |
//    |                                                           0         | alpha[s-1][s-1] |
//------------------------------------------------------------------------------------------------
//  |      |              |                                     |                  |
//
// K[0]   K[1]           K[2]                                 K[s-1]              y_(n+1)
//
//
// With this coefficient storage convention, to structure of the statments to advance the solution
// one timestep has the form
//
// K[0] = F(yn)
//
// For j = 1 to s-1
//
// K[j] = F(y_n + dt*( sum_(0 <= i <= j-1)  alpha[i][j-1]*K[i] )
//
// y_(n+1) = y_(n) + dt* (sum(0 <= i <= s-1) alphaCoeff[i][s-1]*K[i] )
//


class RKsteadyStateCoeff
{
    public:


RKsteadyStateCoeff(){}

~RKsteadyStateCoeff(){}

void initialize() {}


void getRKcoefficients(long stageOrder, double gamma, std::vector< std::vector<double> >& alphaCoefficients)
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

    alphaCoefficients.clear();
    alphaCoefficients.resize(degree);
    for(size_t k = 0; k < degree; k++) {alphaCoefficients[k].resize(degree,0.0);}

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

    private:


void getTchebyShiftFactors(double M, double gamma, double& delta, double& beta, RKpolynomialFunction& ChebyPoly)
{
    RKpolynomialFunction P = ChebyPoly;

    RKpolynomialFunction Pprime;
    Pprime = P.differentiate();

    PprimeCondition Pdelta;
    Pdelta.initialize(P,Pprime,M);

    double a = 0.0;
    double b;

    double tolerance = 1.0e-12;

    if(std::abs(gamma-2.0) > tolerance)
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
// using the recurrence
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




void backSolve(long N, double* c, double** A)
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

void create2Darray(double** &a, long m, long n)
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

void  delete2Darray(double** &a)
{
	double* dataPtr = a[0];
    delete [] dataPtr;
    delete [] a;
}

};

#endif
 
