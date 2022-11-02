//
// RKsteadyState.h 
//
// Author: Chris Anderson  
// (C) UCLA 2012 
//
//
// ############################################################
// RKsteadyState.h 
//
// Version : Oct. 14,2004 
//
//  Uses StabilizedRK instances to provide 
//  the stabilized RK methods. 
//
// Chris Anderson UCLA 
//
// 10/14/04 Removed experimental code
// 06/04/04  Fixed Rollback problem 
// 02/24/06 Added code to improve convergence for oscillatory systems
//###############################################################
//

#include "StabilizedRKconstants.h"
#include "RKsteadyStateCoeff.h"
#include "StabilizedRK.h"
#include "ClassicRK.h"

#include <cmath>

#ifndef RK_STEADY_STATE_
#define RK_STEADY_STATE_
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

   norm2()                                               (2-norm of vector)
   normInf()                                             (maxAbs-norm of vector)
   
   ############################################################################
   
   RKoperator 
   ----------

   An opearator class with the following two member functions:

   void apply(RKvector& Vin, RKvector & Vout)
   void output(long step, double time, double residual, RKvector& V)

   ############################################################################
*/

template <class RKvector, class RKoperator>
class RKsteadyState
{


public : 
             
RKsteadyState() : stabilizedRKmethod()
{
    evaluationCount = 0;
    stepCount       = 0;
    totalTime       = 0.0;
    initializeAdaptiveVariables();
    errorCheckType  = RKnormType::INFNORM;
}

~RKsteadyState()
{
   // if(srkTimestepEstimator != nullptr) delete  srkTimestepEstimator;
}

void initialize(long stageOrder, double gamma, 
RKvector& y0, RKoperator& F)
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

    Ynsave.initialize(Yn);       // Initialize rollback variables 
    FYnsave.initialize(FYn);
    
    errorCheckType = RKnormType::INFNORM;
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

void setOutputFlag(bool flag = true)
{
    outputFlag = flag;
};

void clearOutputFlag()
{
    outputFlag = false;
};

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
    return Yn.getSolutionNorm2();
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


void setVerboseFlag(bool flagVal = true)
{
    verboseFlag = flagVal;
}

void clearVerboseFlag()
{
    verboseFlag = false;
}
//
// Sets the maximial iterations used in the adaptive time-stepping 
// method. 
//
void setMaximalIterations(long maxIter)
{
    stepMax = maxIter;
}
//
// Set the maximial number of rollbacks that will occur with any
// given timestep. If the number of rollbaCOMPONENTS_STABILIZEDRK_STABILIZEDRKCONSTANTS_H_cks exceedes this
// maximal count; then the algorithm just continues with the current
// step. 
// 
//
void setRollBackMax(long rollBackMaxCount)
{
    rollBackMax = rollBackMaxCount;
}
//
// Set the number of valid timesteps that are taken before an
// attempt is made to increase the timestep size. 
//
void setTimestepConstantCount(long constCount)
{
    constCountMax = constCount;
}
//
// Set the factor that is used to increase the timestep, e.g.
// dtNew = increFactor*dtOld 
//
void setTimestepIncrementFactor(double increFactor)
{
    alpha1 = increFactor;
}
//
// Set the factor that is used to decrease the timestep, e.g.
// dtNew = decreFactor*dtOld 
//
void setTimestepDecrementFactor(double decreFactor)
{
    beta1  = decreFactor;
}

void initializeAdaptiveVariables()
{
    stepMax           = 4000; // upper limit on number of steps taken
    constCountMax     = 1;    // number of timesteps to take before
                              // attempting to increase stepsize

    rollBackMax       = 2;    // maximal number of rollbacks >= 3

    alpha1            = 1.5;   // default timestep increment factor
    beta1             = .75 ;  // default timestep decrement factor

    verboseFlag       = 0; 
    outputFlag        = 0;
}
// computeSteadyStateSolution uses an adaptive timestepping routine to
// evolve the systems to steady state.
//
// The adaptive timestepping method is one that attempts to chose
// as large a timestep as possible while preserving a decrease in the
// 2-norm of the residual. 
//
// The timestep adaptation is done using the 2-norm, while the
// stopping condition can be either the 2-norm or the inf-norm.
//
void computeSteadyStateSolution(double initialTimestep, double tol,  RKnormType errorCheckType)
{
    double dt;
    double residNorm1; 
    double residNorm2;

    dt        = initialTimestep;
    initFlag  = 1;
    totalTime = 0.0;
    this->errorCheckType = errorCheckType;

    residNorm1   = getResidualNorm2();

    if(errorCheckType == RKnormType::INFNORM)
    {residualNorm   = getResidualNormMaxAbs();}
    else
    {residualNorm   = residNorm1;}
    
    long constCount    = 0;
    long rollBackCount = 0;
    double dtStar;
    long invalidStep = 0;

    resetStepCount();
    resetTotalTime();
    while((stepCount < stepMax)&&(residualNorm > tol))
    {
        stabilizedRKmethod.advance(Yn,FYn,dt,Yn,FYn);
        residNorm2   = getResidualNorm2();

        if(residNorm2 <= residNorm1) 
        {
            if(constCount < constCountMax)
            {constCount++;}
            else 
            {dt *= alpha1; constCount = 0;}
            stepCount++;
            residNorm1 = residNorm2;
            totalTime += dt;
            Ynsave     = Yn;
            FYnsave    = FYn;
        }
        else
        {
        if(initFlag == 1)            // Allow stepMax rollbacks at the initial timestep
        {rollBackCount = -stepMax;}  // so that we start off with a timestep that decreases
        else                         // the norm.
        {rollBackCount = 0;}

        initFlag = 0;
        dtStar  = dt;
        while((residNorm2 > residNorm1)&&(rollBackCount < rollBackMax))
        {
            rollBack();
            rollBackCount++;
            invalidStep++;
            dtStar *= beta1;
            stabilizedRKmethod.advance(Yn,FYn,dtStar,Yn,FYn);
            residNorm2   = getResidualNorm2();
            if(verboseFlag != 0)
            {
            printf("Rollback %-3ld   %4.4e   %4.4e\n",getStepCount(), residNorm2,dtStar);
            }
        }
        dt = dtStar;
        residNorm1 = residNorm2;
        totalTime += dt;
        stepCount++;
        Ynsave     = Yn;
        FYnsave    = FYn;
        }
        //
        // Update residual
        //

        if(errorCheckType == RKnormType::INFNORM)
        {residualNorm   = getResidualNormMaxAbs();}
        else
        {residualNorm   = residNorm1;}

        if(verboseFlag != 0)
        {
        printf("%-3ld   %4.4e   %4.4e\n",getStepCount(), residualNorm,dt);
        }
      
        if(outputFlag != 0)
        {
        ODEoperator->output(getStepCount(),totalTime,residualNorm,Yn);
        }
      

        }
        if(verboseFlag != 0)
        {
        printf("Total Function Evaluations : %ld \n",getEvaluationCount());
        printf("Total Steps : %ld  Invalid Steps : %ld \n",getStepCount(),invalidStep);
        printf("Total Evolution Time       : %4.4f \n",getTotalTime());
        }
}

void computeSteadyStateSolutionFixedStep(double dt, double tol, RKnormType errorCheckType)
{
    this->errorCheckType = errorCheckType;
    if(errorCheckType == RKnormType::INFNORM)
    {residualNorm   = getResidualNormMaxAbs();}
    else
    {residualNorm   = getResidualNorm2();}

    totalTime = 0.0;
    stepCount = 0;

    while((stepCount < stepMax)&&(residualNorm > tol))
    {
        stabilizedRKmethod.advance(Yn,FYn,dt,Yn,FYn);
        totalTime += dt;
        stepCount++;

        if(errorCheckType == RKnormType::INFNORM)
        {residualNorm   = getResidualNormMaxAbs();}
        else
        {residualNorm   = getResidualNorm2();}
        
        if(verboseFlag != 0)
        {
        printf("%-3ld   %4.4e  \n",stepCount, residualNorm);
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
}


int getMonotoneTimestep(double& currentDt, double& currentResidual, 
    long& reductionIncrement, double reduceFactor)
{
    static double d1          = 0.0;
    static double d2          = 0.0;
    static double d3          = 0.0;
    static long step          = 0;
    static double residualOld = 0.0;

    reductionIncrement = 0;

    if(step == 0) {residualOld  = currentResidual;}
    if(step == 1) {d1 = currentResidual - residualOld;}
    if(step == 2) {d2 = currentResidual - residualOld;}
    if(step >= 3) {d3 = currentResidual - residualOld;}
     
    if(step >= 3) // check for growth or oscillations
    {
        if(d1 > 0.0)
        {
        if(d2 > 0.0) {currentDt *= reduceFactor; reductionIncrement = 1;}  // residual growing 
        else 
        {if(d3 > 0.0){currentDt *= reduceFactor; reductionIncrement = 1;}} // residual oscillating 
        }
    }
    
     
    if(step >= 3) // check for constant iterations
    {
    if(reductionIncrement == 0)
    {
    if(std::abs(d3/currentResidual) < 1.0e-03){currentDt *= reduceFactor; reductionIncrement = 1;}
    }
    }
     


    if(step >= 3)
    {
    d1 = d2;
    d2 = d3;
    }

    step++;
    residualOld = currentResidual;
    return 0;
}

int computeSteadyStateSolutionFixedStepMonotone(double dt, double tol, double reduceFactor,
long maxReductions,RKnormType errorCheckType, double& finalTimestep)
{
    long reductionCount;
    long reductionIncrement;
    
    this->errorCheckType = errorCheckType;
    if(errorCheckType == RKnormType::INFNORM)
    {residualNorm   = getResidualNormMaxAbs();}
    else
    {residualNorm   = getResidualNorm2();}

    totalTime      = 0.0;
    stepCount      = 0;
    reductionCount = 0;

    while((stepCount < stepMax)&&(residualNorm > tol)&&(reductionCount < maxReductions))
    {
        getMonotoneTimestep(dt,residualNorm,reductionIncrement,reduceFactor);
        reductionCount += reductionIncrement;
        stabilizedRKmethod.advance(Yn,FYn,dt,Yn,FYn);
        //if(eigEstFlag == 1)
        //{
        //printf(" Estimated Timestep   : %4.4e \n",getEstimatedTimestep(dt,maximalDtBound));
        //}
        totalTime += dt;
        stepCount++;

        if(errorCheckType == RKnormType::INFNORM)
        {residualNorm   = getResidualNormMaxAbs();}
        else
        {residualNorm   = getResidualNorm2();}
        
        if(verboseFlag != 0)
        {
        printf("%-3ld  %4.4e  %4.4e  \n",stepCount, dt, residualNorm);
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

    finalTimestep    = dt;
    if(reductionCount == maxReductions) {return 1;}

    return 0;
}


int getClassicRKtimestep(double& currentDt, double& currentResidual, 
    long& reductionIncrement, double reduceFactor)
{
    static double d1          = 0.0;
    static double d2          = 0.0;
    static double d3          = 0.0;
    static long step          = 0;
    static double residualOld = 0.0;

    reductionIncrement = 0;

    if(step == 0) {residualOld  = currentResidual;}
    if(step == 1) {d1 = currentResidual - residualOld;}
    if(step == 2) {d2 = currentResidual - residualOld;}
    if(step >= 3) {d3 = currentResidual - residualOld;}
     
    if(step >= 3) // check for growth or oscillations
    {
        if(d1 > 0.0)
        {
        if(d2 > 0.0) {currentDt *= reduceFactor; reductionIncrement = 1;}  // residual growing 
        else 
        {if(d3 > 0.0){currentDt *= reduceFactor; reductionIncrement = 1;}} // residual oscillating 
        }
    }

    if(step >= 3) // check for constant iterations -- increase stepsize to induce damping
    {
    if(std::abs((d3/currentResidual)*(1.0/currentDt)) < 1.0e-02){currentDt *= 1.0/reduceFactor;}
    }

    if(step >= 3)
    {
    d1 = d2;
    d2 = d3;
    }

    step++;
    residualOld = currentResidual;
    return 0;
}
//
// After a specified number of time-step reductions due to an oscillatory
// residual, this routine switches to a classic 4th order Runge-Kutta 
// scheme.
//
int computeSteadyStateMonotoneOSC(double dt, double tol, double reduceFactor,
long maxReductions,long reduceCountSwitch, RKnormType errorCheckType, double& finalTimestep)
{
    double residualNormL2;
    long reductionCount;
    long reductionIncrement;

    long classicRKevaluations = 0;

    residualNormL2   = getResidualNorm2();

    this->errorCheckType = errorCheckType;
    if(errorCheckType == RKnormType::INFNORM)
    {residualNorm   = getResidualNormMaxAbs();}
    else
    {residualNorm   = residualNormL2;}

    totalTime      = 0.0;
    stepCount      = 0;
    reductionCount = 0;

    while((stepCount < stepMax)&&(residualNorm > tol)&&(reductionCount < reduceCountSwitch))
    {
        getMonotoneTimestep(dt,residualNormL2,reductionIncrement,reduceFactor);
        reductionCount += reductionIncrement;
        stabilizedRKmethod.advance(Yn,FYn,dt,Yn,FYn);
        totalTime += dt;
        stepCount++;

        residualNormL2   = getResidualNorm2();
        if(errorCheckType == RKnormType::INFNORM)
        {residualNorm   = getResidualNormMaxAbs();}
        else
        {residualNorm   = residualNormL2;}
        
        if(verboseFlag != 0)
        {
        printf("%-3ld  %4.4e  %4.4e  %4.4e\n",stepCount, dt, residualNorm, residualNormL2);
        }
        
        if(outputFlag != 0)
        {
        ODEoperator->output(getStepCount(),totalTime,residualNorm,Yn);
        }
    }

    if(reductionCount == reduceCountSwitch) // exited due to oscillations 
    {
    classicRK.initialize(*ODEoperator, Yn, 4);

    if(verboseFlag != 0)
    {
        printf("###################### \n");
        printf("Switched To Classic RK \n");
        printf("###################### \n");
    }


    while((stepCount < stepMax)&&(residualNorm > tol)&&(reductionCount < maxReductions))
    {
        getClassicRKtimestep(dt,residualNormL2,reductionIncrement,reduceFactor);
        reductionCount += reductionIncrement;
        classicRK.advance(Yn,FYn,dt,Yn,FYn);
        classicRKevaluations  += 4;
        totalTime += dt;
        stepCount++;

        residualNormL2   = getResidualNorm2();
        if(errorCheckType == RKnormType::INFNORM)
        {residualNorm   = getResidualNormMaxAbs();}
        else
        {residualNorm   = residualNormL2;}
        
        if(verboseFlag != 0)
        {
        printf("%-3ld  %4.4e  %4.4e  %4.4e\n",stepCount, dt, residualNorm, residualNormL2);
        }
    }

    }

    if(verboseFlag != 0)
    {
        printf("Final Residual             : %10.5e \n",residualNorm);
        printf("Total Function Evaluations : %ld     \n",
        classicRKevaluations + stabilizedRKmethod.getEvaluationCount());
        printf("Total Evolution Time       : %4.4f \n",getTotalTime());
    }

    finalTimestep    = dt;
    if(reductionCount == maxReductions) {return 1;}

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

    RKvector           FYnsave;  // Roll-back temporaries
    RKvector           Ynsave;   // 
    
    
    double residualNorm;         // current value of the residual in 
                                 // the norm specified by errorCheckType
//
//  Stabilized RK method instance and parameters
//
    long                                     stageOrder;
    double                                   gamma;
    StabilizedRK<RKvector, RKoperator >      stabilizedRKmethod; 
    ClassicRK <RKvector, RKoperator >        classicRK;

    RKsteadyStateCoeff rkSteadyStateCoeff;
//
//  Adaptive method variables
//
    RKnormType  errorCheckType;
    long        evaluationCount; // ODE apply operator count 
    long        stepCount;
    double      totalTime;

    double                tol; // stopping tolerance
    long              stepMax; // upper limit on number of steps taken
    long        constCountMax; // number of timesteps to take before
                               // attempting to increase stepsize

    long         rollBackMax;  // maximal number of rollbacks to take

    int              initFlag; // Flag inidicating initial startup
    double             alpha1; // initial timestep increment factor
    double             alpha2; // default timestep increment factor
    double              beta1; // default timestep decrement factor

    bool           verboseFlag;  // verbose output flag
    bool           outputFlag;  // flag indicating the invocation of state output
};
#endif

 
