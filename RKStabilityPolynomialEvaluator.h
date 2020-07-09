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

    RKStabilityPolynomialEvaluator();
    RKStabilityPolynomialEvaluator(long stageOrder, double gamma);
    void initialize();
    void initialize(long stageOrder, double gamma);

    ~RKStabilityPolynomialEvaluator();
//
// Returns the norm of the stability polynomial when evaluated
// at a point z = (realZ,imagZ)
//

    double operator()(double realZ,double imagZ);
    double operator()(double realZ);

    double** alphaCoeff;
    double* FYk;
    long stageOrder;
    double    gamma;
    
    std::complex < double >* ZYk;
};
#endif


