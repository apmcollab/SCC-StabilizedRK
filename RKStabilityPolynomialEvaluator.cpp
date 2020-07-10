//
// RKStabilityPolynomialEvaluator.cpp 
//
// Author: Chris Anderson  
// (C) UCLA 2012 
//
#include "RKStabilityPolynomialEvaluator.h"

RKStabilityPolynomialEvaluator::RKStabilityPolynomialEvaluator()
{
    this->stageOrder = 0;
    this->gamma      = 0;
    FYk              = nullptr;
    ZYk              = nullptr;
    alphaCoeff.clear();
}

RKStabilityPolynomialEvaluator::RKStabilityPolynomialEvaluator(long stageOrder, double gamma)
{
    this->stageOrder = stageOrder;
    this->gamma      = gamma;
    FYk              = new double[stageOrder];
    ZYk              = new std::complex< double >[stageOrder];


    RKsteadyStateCoeff rkSteadyStateCoeff;
    rkSteadyStateCoeff.getRKcoefficients(stageOrder,gamma,alphaCoeff);
}

void RKStabilityPolynomialEvaluator::initialize()
{
    this->stageOrder = 0;
    this->gamma      = 0;

    if(FYk != nullptr) delete [] FYk;
    if(ZYk != nullptr) delete [] ZYk;
    FYk  = nullptr;
    ZYk  = nullptr;
    alphaCoeff.clear();
}

void RKStabilityPolynomialEvaluator::initialize(long stageOrder, double gamma)
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

RKStabilityPolynomialEvaluator::~RKStabilityPolynomialEvaluator()
{
	if(FYk != nullptr) delete [] FYk;
	if(ZYk != nullptr) delete [] ZYk;
}

//
// Returns the norm of the stability polynomial when evaluated
// at a point z = (realZ,imagZ)
//

double RKStabilityPolynomialEvaluator::operator()(double realZ,double imagZ)
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

double RKStabilityPolynomialEvaluator::operator()(double realZ)
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


