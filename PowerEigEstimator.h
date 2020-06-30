//
// PowerEigEstimator.h 
//
// Author: Chris Anderson  
// (C) UCLA 2012 
//
#include <iostream>
#include <iomanip>
using namespace std;

#ifndef __PowerEigEstimator__
#define __PowerEigEstimator__
//
// Class PowerEigEstimator
//
// June 22, 2004
//
template <class RKvector, class RKoperator > class
PowerEigEstimator
{
public :

PowerEigEstimator() : Vtmp()
{
    growthFactorMax = 5.0;
}

~PowerEigEstimator(){}

//
// dtInitial : initial guess for an Euler timestep that will yield accurate
// solutions.
//
// Vstar          : The solution value at which the maximal eigenvalue of the 
//                  Jacobian is to be computed. 
// F              : ODE function; we're solving dv/dt = F(v)
// relTol         : Relative accuracy tolerance for the estimate of |lambda_max|
//                  (I've been using .1)
// powMaxIter     : maximal number of power iterations
// eulerStepCount : number of Euler steps used in the construction of the initial
//                  guess. (I've used 1 or 2.)
//
// returns an estimate of |lambda_max| where |lambda_max|is the largest eigenvalue
// of the Jacobian of F at Vstar.
//
double getEigenvalueEstimate(double dtInitial, RKvector& Vstar, RKoperator& F,
double relTol,long powerMaxIter, long eulerStepCount)
{ 
    Vtmp.initialize(Vstar);
    FVstar.initialize(Vstar);
    Veig.initialize(Vstar);

    double dtEuler = determineEulerStep(dtInitial, Vstar, F);

    long k;
//
//  Take a few Euler steps to relax the transients, and then
//  determine the eigenvector approximation by taking the
//  difference between the evolved state and the initial state
//
    F.apply(Vstar,FVstar);

    Vtmp  = FVstar;
    Vtmp *= dtEuler; 

    Veig  = Vstar;
    Veig += Vtmp;
//
//  Additional steps
//
    for(k = 1; k < eulerStepCount; k++)
    {
    F.apply(Veig,Vtmp);
    Vtmp *= dtEuler;
    Veig += Vtmp;
    }
    Veig -= Vstar;
//
//  Power iteration using Veig as an initial guess
//
    double eigEps = dtEuler;
    double eigEst = 0.0;
    double eigEstOld;
    double Vnorm;

    Vnorm   = Veig.nrm2();
    Veig   *= dtEuler/Vnorm;

    long powerCount = 1;
    int doneFlag    = 0;
    
    while((powerCount <= powerMaxIter)&&(doneFlag == 0))
    {
    Veig += Vstar;
    F.apply(Veig, Vtmp);
    Vtmp     -= FVstar;
    Vnorm     = Vtmp.nrm2();
    eigEstOld = eigEst;
    eigEst    = Vnorm/eigEps;
    cout << setw(3) << powerCount << " " << setw(10) << eigEst << "  " << endl;
    Veig  = Vtmp;
    Veig *= eigEps/Vnorm;
    if(powerCount > 1) {if(fabs(eigEst/eigEstOld -1.0) < relTol) doneFlag = 1;}
    powerCount++;
    }
    
    return eigEst;
}

//
// This routine uses a bisection like procedure to determine
// a timestep for which forward Euler will be acceptable. 
// 
double getEulerStepNormRatio(double dt, RKvector& Vstar,  RKoperator& F)
{
    double v0norm;
    double v1norm;
    double normRatio;
    int zeroInitialFlag;

    zeroInitialFlag = 0;
    FV.initialize(Vstar);
    v0norm = Vstar.nrm2();
    //
    // Determine if the initial state 
    // is identically zero
    //
    if(fabs(v0norm) < 1.0e-12)
    {
    zeroInitialFlag = 1; Vkp2.initialize(Vstar);
    }
 
    // Initial Euler Steps

    Vkp1.initialize(Vstar);
    F.apply(Vstar,FV);
    Vkp1.axpy(dt,FV);
    v1norm = Vkp1.nrm2();

    if(zeroInitialFlag)
    {
    Vkp2.initialize(Vkp1);
    F.apply(Vkp1,FV);
    Vkp2.axpy(dt,FV);
    v0norm = v1norm;
    v1norm = Vkp2.nrm2();
    }
   
    normRatio = v1norm/v0norm;
    return normRatio;
}

//
// This routine uses a bisection like procedure to determine
// a timestep for which forward Euler will be acceptable. 
// 
double determineEulerStep(double dtInitial, RKvector& Vstar,  RKoperator& F)
{
    double normRatio;
    double growthFactor;
    double dtStar;


    double growthFactorMax = 4.0;
    double timestepFactor  = 0.25;
  
    int growingModeCount;
    int doneFlag   =  0;
    long maxIter   =  100;
    long iterCount =  0;

    dtStar           = dtInitial;
    normRatio        = getEulerStepNormRatio(dtStar, Vstar,  F);
    growingModeCount = 0;

    while((doneFlag == 0)&&(normRatio > 1.0)&&(iterCount < maxIter))
    {
    if(normRatio < 1.0){doneFlag = 1;}
    else
    {
        growthFactor = (1.0/dtStar)*(normRatio - 1.0);
        if(growthFactor > growthFactorMax)
        {
            dtStar    = timestepFactor*(dtStar/(normRatio + 1));
            normRatio = getEulerStepNormRatio(dtStar, Vstar,  F);
            iterCount++;
        }
        else
        {
            growingModeCount++;
            if(growingModeCount == 1) dtStar = timestepFactor/growthFactorMax;
            if(growingModeCount == 2) dtStar *= 0.5;
            if(growingModeCount == 3) doneFlag = 1;
            if(growingModeCount != 3) 
            {
            normRatio = getEulerStepNormRatio(dtStar, Vstar,  F);
            iterCount++;
            }
        }
    }
    }
    return dtStar;
}

 RKvector   Vtmp;
 RKvector   Vkp1;
 RKvector   Vkp2;
 RKvector   FV;

 double growthFactorMax;

 RKvector   Veig;
 RKvector FVstar;
};
#endif

 
