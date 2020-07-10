//
// RKsteadyStateEE.h 
//
// Author: Chris Anderson  
// (C) UCLA 2012 
//
//
#include <cmath>
#include <vector>

#include "RKsteadyStateCoeff.h"
#include "StabilizedRK.h"
#include "ClassicRK.h"
#include "RKEigEstimator.h"
#include "SRKtimestepEstimator.h"

#ifndef RK_STEADYSTATE_EE_
#define RK_STEADYSTATE_EE_
/**
   A templated class for stabilized Runge-Kutta methods used to
   compute steady state solutions with a timestep chosen using
   Ekeland et.al. maximal eigenvalue detection.

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

   norm2()                                               (2-norm of vector)
   normInf()                                             (Inf-norm of vector)
   
   ############################################################################
   
   RKoperator 
   ----------

   An opearator class with the following two member functions:

   void apply(RKvector& Vin, RKvector & Vout)
   void output(long step, double time, double residual, RKvector& V)

   ############################################################################
*/

template <class RKvector, class RKoperator>
class RKsteadyStateEE
{
public : 
             
RKsteadyStateEE() : stabilizedRKmethod()
{
    evaluationCount = 0;
    stepCount       = 0;
    totalTime       = 0.0;
    
    verboseFlag       = 0; 
    outputFlag        = 0;
    stepMax           = 10000;
    currentDt         = 0.0;
    eigEstFlag        = 0;
    srkTimestepEstimator = 0;
    dtMaxReductionFactor = 0.95;
    thetaReductionFactor = 0.0;
    
    errorCheckType     = INFNORM;
}

~RKsteadyStateEE()
{
	if(srkTimestepEstimator != 0) delete  srkTimestepEstimator;
}

void initialize()
{
    evaluationCount = 0;
    stepCount       = 0;
    totalTime       = 0.0;
    
    verboseFlag       = 0; 
    outputFlag        = 0;
    stepMax           = 10000;
    currentDt         = 0.0;
    initialDt         = 0.0;
    maximalDtBound    = 0.0;
    eigEstFlag        = 0;
    srkTimestepEstimator = 0;
    dtMaxReductionFactor = 0.95;
    thetaReductionFactor = 0.0;
    
    errorCheckType     = INFNORM;
}
void initialize(long stageOrder, double gamma, RKvector& y0, RKoperator& F)
{
	initialize();
	
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

    Ynsave.initialize(Yn);       // Initialize rollback variables 
    FYnsave.initialize(FYn);
    
    srkTimestepEstimator = new SRKtimestepEstimator(stageOrder,gamma);
    
}

//
// Restores the computed solution to the previous timestep. The class only
// saves one previous step; so this member function can only be called once
// after the advance(...) method is called. 
//
void rollBack()       
{
	Yn =  Ynsave;
   FYn = FYnsave;
}

void setDtMaxReductionFactor(double val)
{
	dtMaxReductionFactor = val;
}

void setThetaReductionFactor(double val)
{
	thetaReductionFactor = val;
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
	return Yn.norm2();
}
//
// Returns the inf-norm of the computed solution. 
//
double getSolutionNormMaxAbs()
{
	return Yn.normInf();
}
//
// Returns the 2-norm of the residual.
//
double getResidualNorm2()
{
	return FYn.norm2();
}
//
// Returns the inf-norm of the residual.
//
double getResidualNormMaxAbs()
{
	return FYn.normInf();
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
  double** alphaPtr = rkSteadyStateCoeff.getRKcoefficientsPtr(rkStageOrder, rkGammaFactor);
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

double getTimestep(double dt)
{
	double dtNew;

	estimateEigenSystem(dt);
    if(eigEstFlag == 1){srkTimestepEstimator->setVerboseFlag(1);}
    dtNew = srkTimestepEstimator->evaluateMaximalTimestep(thetaReductionFactor, &Wreal[0],&Wimag[0], eigCount,maximalDtBound);
    dtNew *= dtMaxReductionFactor;

    //
    // Don't allow arbitrary decrease, only a fraction of original timestep
    /*
    if(dtNew < dt*0.75) 
    {
        dtNew = 0.75*dt;
    	totalTime -= dt;
    	rollBack();
    	
        if(errorCheckType == INFNORM)
        {residualNorm   = getResidualNormMaxAbs();}
        else
        {residualNorm   = getResidualNorm2();}
    
    	if(verboseFlag != 0)
        {
          printf("XXXXXXXXXXXXXXXXXXXXX Rollback  \n");
        }
    	return dtNew;
    }
    */
    //
    // Always use minimal timestep, because in steady state computations, the
    // stiff modes vanish. 
    //
    
    // if(dtNew > dtTmp) dtNew = dtTmp;

	return dtNew;
}


double computeInitialTimestep()
{
	double dt;
	double dtNew;
//
//  Evolve the system one timestep to get an initial timestep estimate. If
//  the timestep estimate is small than dtInitial, then roll back the solution
//  and return the reduced timestep. 
//
    dt = initialDt;
	if(verboseFlag != 0)
	{
		printf("XXXX     Initial Timestep Determination     XXXX\n");
	}
	Ynsave  = Yn;
	FYnsave = FYn;
    stabilizedRKmethod.advance(Yn,FYn,dt,Yn,FYn);
    dtNew   = getTimestep(dt);
    if(dtNew < dt)  // roll back
    {
    Yn  = Ynsave;
    FYn = FYnsave;
    dt  = dtNew;
    }
    else
    {
    dt         = initialDt;
    totalTime += dt;
    }
    if(verboseFlag != 0)
	{printf("\nXXXX     Starting dt found  %4.4e     XXXX\n\n",dt);}	
	 
	if(errorCheckType == INFNORM)
    {residualNorm   = getResidualNormMaxAbs();}
    else
	{residualNorm   = getResidualNorm2();}

	return dt;
}

int computeSteadyStateSolution(double dtInitial, double dtMax, double tol, int errorType)
{
    this->initialDt       = dtInitial;
    this->maximalDtBound  = dtInitial;

    double dtNew;
    double reduceFactor = .75;


    int    rollBackFlag      = 0;
    long   stepTaken         = 0;
    double currentDtSave     = 0.0;

    errorCheckType    = errorType;
    initializeRKeigEstimator(stageOrder, gamma);
	
    if(errorCheckType == INFNORM)
    {residualNorm   = getResidualNormMaxAbs();}
    else
    {residualNorm   = getResidualNorm2();}

    totalTime = 0.0;
    stepCount = 0;
    
    currentDt =  computeInitialTimestep();
	
    if(fabs(currentDt) == 0.0) {zeroTimestepError();}

    while((stepCount < stepMax)&&(residualNorm > tol))
    {
    Ynsave = Yn;
    FYnsave = FYn;
    stabilizedRKmethod.advance(Yn,FYn,currentDt,Yn,FYn);
    rollBackFlag      = 0;
    currentDtSave      = currentDt;
    //
    // Reestimate timestep
    
    totalTime += currentDt;
    stepCount++;
    stepTaken++;
       
    if(errorCheckType == INFNORM)
    {residualNorm   = getResidualNormMaxAbs();}
    else
    {residualNorm   = getResidualNorm2();}

    //
    // New (6/25/10) logic behind SteadyStateEE.
    //
    // If the new timestep estimate is much less than the previous timestep,
    // this indicates that the solution has ventured into a stiff region
    // and it is likely that the previous timestep value was computed erroneously. Therefore
    // the solution is rolled back to the previous timestep.
    //
    // The new timestep is used only if it is smaller than the current timestep, e.g.
    // the timestep can only decrease.

    dtNew     = getTimestep(currentDt);
    if(fabs(dtNew) > 1.0e-12)
    {
    if(dtNew/currentDt < 0.8)
    {
    totalTime -= currentDt;
    rollBack();
    rollBackFlag = 1;
    stepTaken --;
    currentDt = dtNew;
    }
    else
    {
    if(currentDt > dtNew) currentDt = dtNew;
    }
    }

    //
    // Only print out diagnostics of the current step has been accepted
    //
    if((verboseFlag != 0)&&(rollBackFlag == 0))
    {
    printf("%-3ld  %15.10e  %15.10e  \n",stepTaken,currentDtSave,residualNorm);
    }

    if((outputFlag != 0)&&(rollBackFlag == 0))
    {
        ODEoperator->output(stepTaken,totalTime,residualNorm,Yn);
    }
    
    }

    if(verboseFlag != 0)
    {
        printf("Final Residual             : %10.5e \n",residualNorm);
        printf("Total Function Evaluations : %ld     \n",stabilizedRKmethod.getEvaluationCount());
        printf("Total Evolution Time       : %4.4f \n",getTotalTime());
    }
    
    return 0;
}

/*
int computeSteadyStateSolution(double dtInitial, long initialCount, 
double dtMax, double tol, int errorCheckType)
{
	double dtTemp;
	
	initializeRKeigEstimator(stageOrder, gamma);
	
    if(errorCheckType == INFNORM)
    {residualNorm   = getResidualNormMaxAbs();}
    else
	{residualNorm   = getResidualNorm2();}
	
	totalTime = 0.0;
    stepCount = 0;
    
	currentDt = estimateInitialTimestep(dtInitial,initialCount, errorCheckType);
	
	if(fabs(currentDt) == 0.0)
	{zeroTimestepError();}
	
	if(currentDt > dtMax) currentDt = dtMax;

	
	while((stepCount < stepMax)&&(residualNorm > tol))
	{
		stabilizedRKmethod.advance(Yn,FYn,currentDt,Yn,FYn);
        totalTime += currentDt;
        stepCount++;
        //
        // Reestimate timestep
        //
        dtTemp    = currentDt;
        currentDt = getTimestep(currentDt);
        if(currentDt > dtMax) currentDt = dtMax;
        
        //
        // Detected only positive eigenvalues, so reset to 
        // previous timestep.
        //
        if(fabs(currentDt) == 0.0) currentDt = dtTemp;
        //
        // Change timestep if required
        //
        if(errorCheckType == INFNORM)
        {residualNorm   = getResidualNormMaxAbs();}
        else
	    {residualNorm   = getResidualNorm2();}
        
        if(verboseFlag != 0)
        {
        printf("%-3ld  %4.4e  %4.4e  \n",stepCount, currentDt,residualNorm);
        }
        
        if(outputFlag != 0)
        {
        ODEoperator->output(getStepCount(),totalTime,residualNorm,Yn);
        }
    }

    if(verboseFlag != 0)
    {
        printf("Final Residual             : %10.5e \n",residualNorm);
        printf("Total Function Evaluations : %ld     \n",stabilizedRKmethod.getEvaluationCount());
        printf("Total Evolution Time       : %4.4f \n",getTotalTime());
    }
    
    return 0;
}
*/
//
//  ODE and solution variables
// 
	RKoperator* ODEoperator;     // ODE


    RKvector             Yn;     // Solution
    RKvector            Ynm1;

    RKvector             Dn;     // Solution
    RKvector            Dnm1;

    RKvector            FYn;     // Residual

    RKvector           FYnsave;  // Roll-back temporaries
    RKvector           Ynsave;   // 
    
    
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
    
    

    RKsteadyStateCoeff rkSteadyStateCoeff;

    RKEigEstimator < RKvector >  rkEigEstimator;
    std::vector<double> Wreal;
    std::vector<double> Wimag;
	long eigCount;
	
    std::vector < RKvector* > stageArrayPointers;
    double stageScaling;
    
    double              currentDt;
    double              initialDt;
    double         maximalDtBound;
    long                eigEstFlag;
    double    dtMaxReductionFactor; 
    double    thetaReductionFactor;
    SRKtimestepEstimator* srkTimestepEstimator;
    
    int             errorCheckType;
    
    
    static void zeroTimestepError()
    {
	cout << "RKsteadyStateEE : Initial steady state evolution failed to detect any " << endl;
	cout << "                  negative eigenvalues of linearized operator.        " << endl;
	cout << "                                                                      " << endl;
	cout << "     Try restarting with a smaller initial timestep.                  " << endl;
	cout << "            XXX  Program Halted XXX "                                   << endl;
    }
};
#endif

 
