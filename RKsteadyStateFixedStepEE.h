//
// RKsteadyStateFixedStepEE.h 
//
// Author: Chris Anderson  
// (C) UCLA 2012 
//
//
// ############################################################
// RKsteadyStateFixedStepEE.h 
//
// May 18, 2009
//
// Convergence to steady state using stabilized RK methods
// with a fixed timestep and eigenvalues estimated using the
// procedure of Ekeland et.al. 
//
//
//###############################################################
//
#include <cmath>
#include <vector>
using std::min;
using std::max;

#include "ArrayStructure1D.h"
#include "RKsteadyStateCoeff.h"
#include "StabilizedRK.h"
#include "ClassicRK.h"
#include "RKEigEstimator.h"
#include "SRKtimestepEstimator.h"

#ifndef __RKsteadyStateFixedStepEE__
#define __RKsteadyStateFixedStepEE__ 
/**
   A templated class for stabilized Runge-Kutta methods used to
   compute steady state solutions.

   The minimal functionality required of the classes 
   that are used in this template

   RKvector  
   ---------
   A vector class with the following member functions:

   RKvector()                                            (null constructor)
   RKvector(const RKvector&)                             (copy constructor)

   initialize(const RKvector&)                           (copy initializer)
  
   operator =                                            (duplicate assignemnt)
   operator +=                                           (incremental addition)
   operator *(double alpha)                              (scalar multiplication)

   axpy(double alpha, RKvector& x)   EE                    (*this =  alpha*x + *this)
   nrm2()                                                (2-norm of vector)
   amax()                                                (maxAbs-norm of vector)
   
   ############################################################################
   
   RKoperator 
   ----------

   An opearator class with the following two member functions:

   void apply(RKvector& Vin, RKvector & Vout)
   void output(long step, double time, double residual, RKvector& V)

   ############################################################################
*/

template <class RKvector, class RKoperator>
class RKsteadyStateFixedStepEE
{
public : 
             
RKsteadyStateFixedStepEE() : stabilizedRKmethod()
{
    evaluationCount = 0;
    stepCount       = 0;
    totalTime       = 0.0;
    
    verboseFlag       = 0; 
    outputFlag        = 0;
    stepMax              = 10000;
    currentDt            = 0.0;
    eigEstFlag           = 0;
    srkTimestepEstimator = 0;
    
    errorCheckType     = INFNORM;
}

~RKsteadyStateFixedStepEE()
{
	if(srkTimestepEstimator != 0) delete  srkTimestepEstimator;
}

void initialize(long stageOrder, double gamma, RKvector& y0, RKoperator& F)
{
    this->stageOrder = stageOrder;
    this->gamma      = gamma;

    stabilizedRKmethod.initialize(stageOrder,gamma,y0,F);
    
    ODEoperator = &F;     
        
    Yn.initialize(y0);            // Initialize local solution values 
    FYn.initialize(y0);      
    Ynm1.initialize(y0);          // Initialize local solution values     

    Dn.initialize(y0);            // Initialize difference vectors
    Dnm1.initialize(y0);          

    ODEoperator->apply(Yn,FYn);   // Compute initial residual   
    evaluationCount        = 1;

    srkTimestepEstimator = new SRKtimestepEstimator(stageOrder,gamma);
    
}



long getEvaluationCount()
{return evaluationCount + stabilizedRKmethod.getEvaluationCount();};

void resetEvaluationCount()
{stabilizedRKmethod.resetEvaluationCount(); evaluationCount = 0;}

long getStepCount()
{return stepCount;};

void resetStepCount()
{stepCount = 0;}

double getTotalTime()
{return totalTime;}

void resetTotalTime()
{totalTime = 0.0;}

void setOutputFlag(int flag)
{outputFlag = flag;};

double getCurrentTimestep()
{return currentDt;}

void setEigEstOutputFlag(int flagVal)
{eigEstFlag = flagVal;}

void setVerboseFlag(int flagVal)
{verboseFlag = flagVal;}

//
// Returns the currently computed solution 
//
RKvector getSolution()
{
	return Yn;
}
//
// Returns the 2-norm of the computed solution. 
//
double getSolutionNorm2()
{
	return Yn.nrm2();
}
//
// Returns the inf-norm of the computed solution. 
//
double getSolutionNormMaxAbs()
{
	return Yn.amax();
}
//
// Returns the 2-norm of the residual.
//
double getResidualNorm2()
{
	return FYn.nrm2();
}
//
// Returns the inf-norm of the residual.
//
double getResidualNormMaxAbs()
{
	return FYn.amax();
}

//
// Returns the currently evaluated residual in 
// the norm specified by the errorCheckType 
//

double getResidualNorm()
{
	return residualNorm;
}

//
// Sets the maximial iterations used in the adaptive time-stepping 
// method. 
//
void setMaximalIterations(long maxIter)
{
	stepMax = maxIter;
}


void initializeRKeigEstimator(long rkStageOrder, double rkGammaFactor)
{
// 
//  Set up RKEigEstimator 
// 
    double** alphaPtr = RKsteadyStateCoeff::getRKcoefficientsPtr(rkStageOrder, rkGammaFactor);
	long i; long j; 
	
	DoubleArrayStructure2D RKcoefficients(rkStageOrder-1,rkStageOrder-1);
	for(i = 0; i < rkStageOrder-1; i++)
	{
	for(j = 0; j < rkStageOrder-1; j++)
	{
		RKcoefficients(i,j)    = alphaPtr[i][j];
	}}
	
	rkEigEstimator.initialize(rkStageOrder,RKcoefficients);
}


void estimateEigenSystem(double dt)
{
	long   i;
	int info;
	
	stageArrayPointers.clear();
	for(i = 0; i < stageOrder; i++) 
    {
    stageArrayPointers.push_back(stabilizedRKmethod.FYk[i]);
    }
    stageScaling = 1.0/dt;
    info = rkEigEstimator.estimateEigenvalues(stageArrayPointers,stageScaling,Wreal,Wimag,eigCount,dt);
}

double getTimestep(double dt,double maximalDtBound)
{
	double dtTmp = dt;
	double dtNew;
	long i;
	estimateEigenSystem(dt);
	
    if(eigEstFlag == 1){srkTimestepEstimator->setVerboseFlag(1);}
	//
	// Note setting theta reduction factor = 0.0
	//
	thetaReductionFactor = 0.0;
    dtNew = srkTimestepEstimator->evaluateMaximalTimestep(thetaReductionFactor,Wreal.getDataPointer(),
    Wimag.getDataPointer(), eigCount,maximalDtBound);
	return dtNew;
}

int computeSteadyStateSolution(double dtInitial, double tol, int errorType)
{
    double maximalDtBound = dtInitial;
    double dtTemp;
    double residualSave;
    double growthMax  = 100000000.0;
    double growthRate;
    double residualNorm2;
       
    errorCheckType    = errorType;
	
    initializeRKeigEstimator(stageOrder, gamma);
	
    if(errorCheckType == INFNORM)
    {residualNorm   = getResidualNormMaxAbs();}
    else
    {residualNorm   = getResidualNorm2();}
    
    residualSave = residualNorm;
 
	
    totalTime = 0.0;
    stepCount = 0;
    
    while((stepCount < stepMax)&&(residualNorm > tol))
    {
    stabilizedRKmethod.advance(Yn,FYn,dtInitial,Yn,FYn);
    totalTime += dtInitial;
    stepCount++;
        
    if(errorCheckType == INFNORM)
    {residualNorm   = getResidualNormMaxAbs();}
    else
    {residualNorm   = getResidualNorm2();}
        
    if(verboseFlag != 0)
    {
    residualNorm2 = getResidualNorm2();
    printf("%-3d  %4.4e  %4.4e  %4.4e \n",stepCount, dtInitial,residualNorm,residualNorm2);
    }
    //
    // Reestimate timestep
    //
    currentDt = getTimestep(dtInitial,maximalDtBound);
    cout << " Estimated Timestep   : " << currentDt << endl;
    //
    // Detected only positive eigenvalues, so reset to 
    // previous timestep.
    //
    if(outputFlag != 0)
    {
        ODEoperator->output(getStepCount(),totalTime,residualNorm,Yn);
    }
    residualSave = residualNorm;
    }

    if(verboseFlag != 0)
    {
        printf("Final Residual             : %10.5e \n",residualNorm);
        printf("Total Function Evaluations : %d     \n",stabilizedRKmethod.getEvaluationCount());
        printf("Total Evolution Time       : %4.4f \n",getTotalTime());
    }
    
    return 0;
}
 
//
//  ODE and solution variables
// 
	RKoperator* ODEoperator;     // ODE


    RKvector             Yn;     // Solution
    RKvector            Ynm1;

    RKvector             Dn;     // Solution
    RKvector            Dnm1;

    RKvector            FYn;     // Residual
    
    
    double residualNorm;         // current value of the residual in 
                                 // the norm specified by errorCheckType
//
//  Stabilized RK method instance and parameters
//
    long                                       stageOrder;          
    double                                          gamma;      
    StabilizedRK<RKvector, RKoperator >      stabilizedRKmethod; 
    ClassicRK <RKvector, RKoperator >        classicRK;
//
//  Adaptive method variables
//
    long        evaluationCount; // ODE apply operator count 
    long        stepCount;
    double      totalTime;   

    enum {INFNORM, L2NORM};   

	double                tol; // stopping tolerance
    long              stepMax; // upper limit on number of steps taken
    int           verboseFlag;  // verbose output flag
    int            outputFlag;  // flag indicating the invocation of state output 
    
    
    RKEigEstimator < RKvector >  rkEigEstimator;
    DoubleArrayStructure1D Wreal;
	DoubleArrayStructure1D Wimag;
	long eigCount;
	
    std::vector < RKvector* > stageArrayPointers;
    double stageScaling;
    
    double              currentDt;
    long                eigEstFlag;
    double    dtMaxReductionFactor; 
    double    thetaReductionFactor;
    SRKtimestepEstimator* srkTimestepEstimator;
    
    int             errorCheckType;
    
    
    static void zeroTimestepError()
    {
	cout << "RKsteadyStateFixedStepEE : Initial steady state evolution failed to detect any " << endl;
	cout << "                  negative eigenvalues of linearized operator.        " << endl;
	cout << "                                                                      " << endl;
	cout << "     Try restarting with a smaller initial timestep.                  " << endl;
	cout << "            XXX  Program Halted XXX "                                   << endl;
    }
};
#endif

 
