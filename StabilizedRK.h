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

#ifndef STABILIZED_RK_
#define STABILIZED_RK_
/**
   A templated class for stabilized Runge-Kutta methods.

   The minimal functionality required of the classes 
   that are used in this template

   RKvector  
   ---------
   A vector class with the following member functions:

   RKvector()                                            (null constructor)
   RKvector(const RKvector&)                             (copy constructor)

   initialize()                                          (null initializer)
   initialize(const RKvector&)                           (copy initializer)
  
   operator =                                            (duplicate assignemnt)
   operator +=                                           (incremental addition)
   operator *(double alpha)                              (scalar multiplication)

   axpy(double alpha, RKvector& x)                       (*this =  alpha*x + *this)
   norm2()                                               (2-norm of vector)
   normInf()                                             (infinity-norm of vector)
   
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


StabilizedRK()
{
    initialize();
}

// Copy constructor - creates duplicate
// (including pointer to ODE operator)

StabilizedRK(const StabilizedRK& R)
{
	stageOrder  = R.stageOrder;
	RKgamma     = R.RKgamma;
	alphaCoeff  = R.alphaCoeff;

	ODEoperator = R.ODEoperator;

	evaluationCount = R.evaluationCount;
    FYk             = R.FYk;

	Yn.initialize(R.Yn);
    FYn.initialize(R.FYn);
    Ytmp.initialize(R.Ytmp);
    FYnsave.initialize(R.FYnsave);
    Ynsave.initialize(R.Ynsave);
}

StabilizedRK(long stageOrder, double RKgamma, RKvector& y0, RKoperator& F)
{
	initialize(stageOrder,RKgamma,y0,F);
}


~StabilizedRK()
{}

void initialize()
{
    stageOrder = 0;
    RKgamma    = 0;
    alphaCoeff.clear();

    rkSteadyStateCoeff.initialize();

	ODEoperator = nullptr;

    evaluationCount = 0;

    FYk.clear();

	Yn.initialize();
    FYn.initialize();
    Ytmp.initialize();
    FYnsave.initialize();
    Ynsave.initialize();
}
void initialize(long stageOrder, double RKgamma, RKvector& y0, RKoperator& F)
{
    this->stageOrder = stageOrder;
    this->RKgamma    = RKgamma;

    evaluationCount  = 0;
 
	ODEoperator      = &F;       //  set operator

	Yn.initialize(y0);           //  Solution 
    FYn.initialize(y0);          //  Residual

    Ytmp.initialize(y0);
	ODEoperator->apply(Yn,FYn);  // Initial residual

    Ynsave.initialize(Yn);       // Rollback temporaries
    FYnsave.initialize(FYn);

   createStageTemporaries(stageOrder);
   rkSteadyStateCoeff.getRKcoefficients(stageOrder, RKgamma, alphaCoeff);
}

void createStageTemporaries(long stageOrder)
{
    if((long)FYk.size() >= stageOrder) {return;}
    FYk.resize(stageOrder,Yn);
}

void setInitialCondition(RKvector &Y0)
{
  Yn      = Y0;
  applyOp(Yn, FYn);

  Ynsave  = Yn; 
  FYnsave = FYn;
}

void setOperator(RKoperator& F)
{
	ODEoperator = &F;
}

void advance(double dt)
{
    Ynsave  = Yn;     
    FYnsave = FYn;
    advance(Ynsave, FYnsave, stageOrder, RKgamma, dt, Yn, FYn);
}

void advance(RKvector &Yin,  RKvector& FYin, double dt, RKvector &Yout, RKvector& FYout)
{
    advance(Yin, FYin, stageOrder, RKgamma, dt, Yout, FYout);
}

double advance(RKvector& Yin, double dt, RKvector &Yout, RKvector& FYout)
{
    applyOp(Yin, FYn);
    advance(Yin, FYn, stageOrder, RKgamma, dt, Yout, FYout);
    return FYout.norm2();
}
//
//  Main RK advance routine 
//
void advance(RKvector &Yin, RKvector& FYin, long sOrder, double sFactor, 
double dt, RKvector &Yout, RKvector& FYout)
{
   long i; long k;

   if((sOrder != stageOrder)||(sFactor != RKgamma))
   {
   createStageTemporaries(sOrder);
   rkSteadyStateCoeff.getRKcoefficients(sOrder, sFactor, alphaCoeff);
   stageOrder = sOrder;
   RKgamma      = sFactor;
   }

   FYk[0] = FYin;
   FYk[0] *= dt;

   for(k = 1; k < sOrder; k++)
   {
   Ytmp = Yin;
   for(i = 0; i < k; i++)
   {
     Ytmp.axpy(alphaCoeff[i][k-1],FYk[i]);
   }
   applyOp(Ytmp,FYk[k]);
   FYk[k] *= dt;
   }

   Ytmp = Yin;

   for(k = 0; k < sOrder; k++)
   {
   Ytmp.axpy(alphaCoeff[k][sOrder-1],FYk[k]);
   }

   Yout = Ytmp;
//
//  Compute new residual
//
    applyOp(Yout,FYout);
}

long getEvaluationCount()
{
	 return evaluationCount;
};

void resetEvaluationCount()
{
	evaluationCount = 0;
}
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
    long  stageOrder;             // Stage order (number of stages )
    double     RKgamma;             // Stability region from [-gamma*stageOrder*stageOrder,0].
                                  // Need RKgamma <= 2
    
    std::vector<std::vector<double>> alphaCoeff; // RK stage coefficients

	RKoperator* ODEoperator;       // ODE
    long        evaluationCount;   // ODE apply operator count

    RKvector                  Yn;  // Solution
    RKvector                 FYn;  // Residual

    RKvector                Ytmp;  // Temporary
    std::vector<RKvector>    FYk;  // Array for stage components

    RKvector             FYnsave; // Roll-back temporaries
    RKvector              Ynsave; //

    RKsteadyStateCoeff rkSteadyStateCoeff;
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

 
