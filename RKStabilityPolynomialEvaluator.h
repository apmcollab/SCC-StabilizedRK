//
// RKStabilityPolynomialEvaluator.h 
//
// Author: Chris Anderson  
// (C) UCLA 2012 
//
#include "RKsteadyStateCoeff.h"
#include <complex>


#ifndef RK_STABILITY_POLYNOMIAL_EVALUATOR_
#define RK_STABILITY_POLYNOMIAL_EVALUATOR_

class RKStabilityPolynomialEvaluator
{
public :

RKStabilityPolynomialEvaluator()
{
    this->stageOrder = 0;
    this->gamma      = 0;
    FYk              = nullptr;
    ZYk              = nullptr;
    alphaCoeff.clear();
}

RKStabilityPolynomialEvaluator(long stageOrder, double gamma)
{
    this->stageOrder = stageOrder;
    this->gamma      = gamma;
    FYk              = new double[stageOrder];
    ZYk              = new std::complex< double >[stageOrder];


    RKsteadyStateCoeff rkSteadyStateCoeff;
    rkSteadyStateCoeff.getRKcoefficients(stageOrder,gamma,alphaCoeff);
}

void initialize()
{
    this->stageOrder = 0;
    this->gamma      = 0;

    if(FYk != nullptr) delete [] FYk;
    if(ZYk != nullptr) delete [] ZYk;
    FYk  = nullptr;
    ZYk  = nullptr;
    alphaCoeff.clear();
}

void initialize(long stageOrder, double gamma)
{
    this->stageOrder = stageOrder;
    this->gamma      = gamma;

    if(FYk != nullptr) delete [] FYk;
    if(ZYk != nullptr) delete [] ZYk;

    FYk              = new double[stageOrder];
    ZYk              = new std::complex<double>[stageOrder];


    RKsteadyStateCoeff rkSteadyStateCoeff;
    rkSteadyStateCoeff.getRKcoefficients(stageOrder,gamma,alphaCoeff);
}

~RKStabilityPolynomialEvaluator()
{
	if(FYk != nullptr) delete [] FYk;
	if(ZYk != nullptr) delete [] ZYk;
}

//
// Returns the norm of the stability polynomial when evaluated
// at a point z = (realZ,imagZ)
//

double operator()(double realZ,double imagZ)
{
    std::complex<double> z(realZ,imagZ);
    std::complex<double> rho;
    std::complex<double> Ytmp;

    long i; long k;


    ZYk[0] = z;
    for(k = 1; k < stageOrder; k++)
    {
    Ytmp = 1.0;
    for(i = 0; i < k; i++)
    {
     Ytmp += alphaCoeff[i][k-1]*ZYk[i];
    }
    ZYk[k] = z*Ytmp;
    }

    rho = 1.0;
    for(k = 0; k < stageOrder; k++)
    {
    rho += alphaCoeff[k][stageOrder-1]*ZYk[k];
    }

    return abs(rho);
}

double operator()(double realZ)
{
    double rho;
    double Ytmp;
    long i; long k;

    double z= realZ;

    FYk[0] = z;
    for(k = 1; k < stageOrder; k++)
    {
    Ytmp = 1.0;
    for(i = 0; i < k; i++)
    {
     Ytmp += alphaCoeff[i][k-1]*FYk[i];
    }
    FYk[k] = z*Ytmp;
    }

    rho = 1.0;
    for(k = 0; k < stageOrder; k++)
    {
    rho += alphaCoeff[k][stageOrder-1]*FYk[k];
    }
    return rho;
}

    std::vector<std::vector<double> > alphaCoeff;

    double*     FYk;
    long stageOrder;
    double    gamma;
    
    std::complex < double >* ZYk;

};
#endif


