//
// StabilizedRK.h 
//
// Author: Chris Anderson  
// (C) UCLA 2012 
//
//
//###############################################################
// StabilizedRK.h 
//
// Version : May 26, 2004                 Chris Anderson UCLA 
//###############################################################
//
#include "RKsteadyStateCoeff.h"

#ifndef __StabilizedRK__
#define __StabilizedRK__
/**
   A templated class for stabilized Runge-Kutta methods.

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

   axpy(double alpha, RKvector& x)                       (*this =  alpha*x + *this)
   norm2()                                                (2-norm of vector)
   normInf()                                                (maxAbs-norm of vector)
   
   ############################################################################
   
   RKoperator 
   ----------

   An opearator class with the following member function:

   void apply(RKvector& Vin, RKvector & Vout)

   ############################################################################

April 05, 2004                                            Chris Anderson  
*/

template <class RKvector, class RKoperator>
class StabilizedRK
{
public : 

//
// Instantiation of the stabilized RungeKutta method. 
//
// stageOrder = the number of stages used in the Runga-Kutta method.
//              If the problem is not very stiff, then low order 
//              (1-3) methods work best, while high order methods
//              (4-10) work better for mildly stiff, or stiff equations.
//
// gamma      = The factor that determinines the size of the stability 
//              region and the magnitude of the damping associated with of the
//              stabilized method. (0 < gamma < 2.0) The stability interval is 
//              [-gamma*stageOrder*stageOrder, 0]. 
//
//              Larger gamma gives methods with larger stability regions, but
//              smaller damping; while smaller gamma gives smaller stability 
//              regions with more damping. 
//
//              Good values for gamma are 1.5 <= gamma <= 1.75. 
//   
//              
StabilizedRK()
{
    this->stageOrder = 0;
    this->gamma      = 0; 

    alphaCoeff       = 0;
    evaluationCount  = 0;
    FYkArraySize      = 0;
    FYk               = 0;
}

~StabilizedRK()
{ 
  long i; 
  if(FYk != 0) 
  {
  for(i = 0; i < FYkArraySize; i++) {delete FYk[i];}
  delete [] FYk;
  }
  
}

void initialize(long stageOrder, double gamma, 
RKvector& y0, RKoperator& F)
{
    this->stageOrder = stageOrder;
    this->gamma      = gamma;

    evaluationCount  = 0;
 
	ODEoperator = &F;            //  set operator

	Yn.initialize(y0);           //  Solution 
    FYn.initialize(y0);          //  Residual

    Ytmp.initialize(y0);
	ODEoperator->apply(Yn,FYn);  // Initial residual

    Ynsave.initialize(Yn);       // Rollback temporaries
    FYnsave.initialize(FYn);
}

void createStageTemporaries(long stageOrder)
{
    if(FYkArraySize >= stageOrder) return;

    long k;

    RKvector** FYkPtr  = new RKvector*[stageOrder];
   
    //
    // copy existing stage temporary pointers
    //

    for(k = 0; k < FYkArraySize; k++)
    {
    FYkPtr[k] = FYk[k];
    }

    //
    // create new stage temporaries
    //

    for(k = FYkArraySize; k < stageOrder; k++)
    {
    FYkPtr[k] = new RKvector(Yn);
    }
    //
    // remove old array of pointers 
    //
    if(FYk != 0) delete [] FYk;

    // assign new array pointer

    FYk = FYkPtr;
    FYkArraySize = stageOrder;
}

void setInitialCondition(RKvector &Y0)
{
  Yn      = Y0;
  applyOp(Yn, FYn);
  Ynsave  = Yn; 
  FYnsave = FYn;
}

void advance(double dt)
{
    Ynsave  = Yn;     
    FYnsave = FYn;
    advance(Ynsave, FYnsave, stageOrder, gamma, dt, Yn, FYn);
}

void advance(RKvector &Yin,  RKvector& FYin, double dt, 
RKvector &Yout, RKvector& FYout)
{
    advance(Yin, FYin, stageOrder, gamma, dt, Yout, FYout);
}

double advance(RKvector& Yin, double dt, RKvector &Yout, RKvector& FYout)
{
    applyOp(Yin, FYn);
    advance(Yin, FYn, stageOrder, gamma, dt, Yout, FYout);
    return FYout.norm2();
}
//
//  Main RK advance routine 
//
void advance(RKvector &Yin, RKvector& FYin, long sOrder, double sFactor, 
double dt, RKvector &Yout, RKvector& FYout)
{
   long i; long k;

   createStageTemporaries(sOrder);
   alphaCoeff  = RKsteadyStateCoeff::getRKcoefficientsPtr(sOrder, sFactor);

   *FYk[0] = FYin;
   *FYk[0] *= dt;

   for(k = 1; k < sOrder; k++)
   {
   Ytmp = Yin;
   for(i = 0; i < k; i++)
   {
     Ytmp.axpy(alphaCoeff[i][k-1],*FYk[i]);
   }
   applyOp(Ytmp,*FYk[k]);
   *FYk[k] *= dt;
   }

   Ytmp = Yin;

   for(k = 0; k < sOrder; k++)
   {
   Ytmp.axpy(alphaCoeff[k][sOrder-1],*FYk[k]);
   }

   Yout = Ytmp;
//
//  Compute new residual
//
    applyOp(Yout,FYout);
}

    long getEvaluationCount()
    {return evaluationCount;};

    void resetEvaluationCount()
    {evaluationCount = 0;}
//
//  applyOp(...) 
//  applies the ODE operator and keeps track of how many times it's been
//  called.
//
void applyOp(RKvector& Vin, RKvector& Vout)
{
          ODEoperator->apply(Vin,Vout);
          evaluationCount++;
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
//  Class variables for RK evolution 
//
    long  stageOrder;            // Stage order (number of stages ) 
    double     gamma;            // Stability region from [-gamma*stageOrder*stageOrder,0]. 
                                 // Need gamma <= 2
    
    double** alphaCoeff;         // RK stage coefficients 

	RKoperator* ODEoperator;     // ODE
    long        evaluationCount; // ODE apply operator count 


    RKvector             Yn;     // Solution
    RKvector            FYn;     // Residual

    RKvector            Ytmp;     // Temporary
    RKvector**            FYk;     // Array for stage components
    long         FYkArraySize;   

    
    RKvector           FYnsave; // Roll-back temporaries
    RKvector           Ynsave;  // 
};
#endif


/*
    long i; long j;
    for(j = 0; j < stageOrder; j++)
    {
    for(i = 0; i < stageOrder; i++)
    {
    printf("%10.7e  ",alphaCoeff[i][j]);
    }
    printf("\n \n");
    }
*/ 

 
