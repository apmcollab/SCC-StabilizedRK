/*
 * RKtimestepEstimator.h
 *
 *  Created on: Jul 10, 2020
 *      Author: anderson
 */

#ifndef RK_EIG_STEP_ESTIMATOR_
#define RK_EIG_STEP_ESTIMATOR_

#include <vector>

#include "StabilizedRK.h"
#include "RKsteadyStateCoeff.h"
#include "RKEigEstimator.h"
#include "SRKtimestepEstimator.h"


template <class RKvector, class RKoperator> class RKeigStepEstimator
{

public:

void setVerbose(bool val = true)
{
	verboseFlag = val;
}

void clearVerbose()
{
	verboseFlag = false;
}

void initialize(StabilizedRK <RKvector, RKoperator >& stabilizedRKmethod)
{
    verboseFlag = false;
	stabilizedRKmethodPtr = &stabilizedRKmethod;
}
//
//#################################################################################
//  Added for eigenvalue estimation and a-postiori timestep size estimation
//#################################################################################

void setEigEstOutputFlag(bool flagVal = true)
{eigEstFlag = flagVal;}

void initializeEigRoutines()
{
//
//  Set up RKEigEstimator
//
    std::vector<std::vector<double>> alphaCoeff;

    rkSteadyStateCoeff.getRKcoefficients(stageOrder, RKgamma,alphaCoeff);

    std::vector< std::vector<double> > RKcoefficients;

    RKcoefficients.resize((long)(stageOrder-1));
    for(size_t i = 0; i < stageOrder-1; ++i){RKcoefficients[i].resize(stageOrder-1,0.0);}

    for(size_t i = 0; i < stageOrder-1; i++)
    {
    for(size_t j = 0; j < stageOrder-1; j++)
    {
        RKcoefficients[i][j]    = alphaCoeff[i][j];
    }}

    rkEigEstimator.initialize(stageOrder,RKcoefficients);

    if(srkTimestepEstimator != nullptr) delete  srkTimestepEstimator;
    srkTimestepEstimator = new SRKtimestepEstimator(stageOrder,RKgamma);
}

void estimateEigenSystem(double dt)
{
    long   i;
    int info;

    stageArrayPointers.clear();
    for(i = 0; i < stageOrder; i++)
    {
    stageArrayPointers.push_back(&(stabilizedRKmethodPtr->FYk[i]));
    }
    stageScaling = 1.0/dt;
    rkEigEstimator.estimateEigenvalues(stageArrayPointers,stageScaling,Wreal,Wimag,eigCount,dt);
}

double getEstimatedTimestep(double dt,double maximalDtBound)
{
    double dtNew;
    estimateEigenSystem(dt);
    if(eigEstFlag == 1){srkTimestepEstimator->setVerboseFlag(1);}
    //
    // Note setting theta reduction factor = 0.0
    //
    thetaReductionFactor = 0.0;
    dtNew = srkTimestepEstimator->evaluateMaximalTimestep(thetaReductionFactor,&Wreal[0], &Wimag[0], eigCount,maximalDtBound);
    return dtNew;
}


double computeInitialTimestep(double dtInitial)
{
    double dt;
    double dtNew;
//
//  Evolve the system one timestep to get an initial timestep estimate. If
//  the timestep estimate is smaller than dtInitial, then roll back the solution
//  and return the reduced timestep.
//
    dt = dtInitial;
    if(verboseFlag != 0)
    {
        printf("XXXX     Initial Timestep Determination     XXXX\n");
    }
    stabilizedRKmethodPtr->Ynsave  = stabilizedRKmethodPtr->Yn;
    stabilizedRKmethodPtr->FYnsave = stabilizedRKmethodPtr->FYn;
    stabilizedRKmethodPtr->advance(stabilizedRKmethodPtr->Yn,stabilizedRKmethodPtr->FYn,dt,
                                   stabilizedRKmethodPtr->Yn,stabilizedRKmethodPtr->FYn);
    dtNew   = getEstimatedTimestep(dt,2.0*dtInitial);
    if(dtNew < dt)  // roll back
    {
    stabilizedRKmethodPtr->Yn  = stabilizedRKmethodPtr->Ynsave;
    stabilizedRKmethodPtr->FYn = stabilizedRKmethodPtr->FYnsave;
    dt  = dtNew;
    }
    else
    {
    dt         = dtInitial;
    stabilizedRKmethodPtr->totalTime += dt;
    }
    if(verboseFlag != 0)
    {printf("\nXXXX     Starting dt found  %4.4e     XXXX\n\n",dt);}

    if(stabilizedRKmethodPtr->errorCheckType == StabilizedRK <RKvector, RKoperator >::INFNORM)
    {stabilizedRKmethodPtr->residualNorm   = stabilizedRKmethodPtr->getResidualNormMaxAbs();}
    else
    {stabilizedRKmethodPtr->residualNorm   = stabilizedRKmethodPtr->getResidualNorm2();}

    return dt;
}

    bool verboseFlag;

    StabilizedRK <RKvector, RKoperator >* stabilizedRKmethodPtr;

	double RKgamma;
	int    stageOrder;

    RKEigEstimator < RKvector >  rkEigEstimator;
    std::vector<double> Wreal;
    std::vector<double> Wimag;

    std::vector < RKvector* > stageArrayPointers;

    SRKtimestepEstimator* srkTimestepEstimator;

    double   stageScaling;
    double   dtMaxReductionFactor;
    double   thetaReductionFactor;
    bool     eigEstFlag;
    long     eigCount;

    RKsteadyStateCoeff rkSteadyStateCoeff;
};
#endif /* SCC_STABILIZEDRK_RKINITIALSTEPESTIMATOR_H_ */
