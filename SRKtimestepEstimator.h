//
// SRKtimestepEstimator.h 
//
// Author: Chris Anderson  
// (C) UCLA 2012 
//
#include <iostream>
#include <cmath>
#include <vector>
#include <cmath>

#include "RKStabilityPolynomialEvaluator.h"
#include "rkf45.h"


#ifndef SRK_TIMESTEP_ESTIMATOR_
#define SRK_TIMESTEP_ESTIMATOR_


class SRKstabilityContourODE
{
	public :
	
	SRKstabilityContourODE(long stageOrder, double gamma)
	{
	this->stageOrder   = stageOrder;
    this->gamma        = gamma;
    stabilityPoly.initialize(stageOrder,gamma);
    eps = 1.0e-08;
	}
	
	void operator()(double* t, double* y, double* yp)
	{
	double zReal = y[0];
	double zImag = y[1];
	
	double dphiX = (stabilityPoly(zReal+eps,zImag) - stabilityPoly(zReal-eps,zImag))/(2.0*eps);
	double dphiY = (stabilityPoly(zReal,zImag+eps) - stabilityPoly(zReal,zImag-eps))/(2.0*eps);
	double gradPhiNorm = sqrt(dphiX*dphiX + dphiY*dphiY);
	yp[0] = -dphiY/gradPhiNorm;
	yp[1] =  dphiX/gradPhiNorm;
	};
	
	double eps;
    long stageOrder;
    double gamma;
    RKStabilityPolynomialEvaluator stabilityPoly;
};


class SRKtimestepEstimator
{
	public :
	
	SRKtimestepEstimator(long stageOrder, double gamma) :
	stabilityContourODE(stageOrder,gamma)
	{
	this->stageOrder   = stageOrder;
    this->gamma        = gamma;	
    
    //
    // Reset gamma if required so that gamma < 2.0
    //
    if(this->gamma >= 2.0) {this->gamma = 1.98;}
    
    stabilityPoly.initialize(this->stageOrder,this->gamma);
    
    createStabilityBoundsTables();
    verboseFlag = 0;
	}
	
	void setVerboseFlag(int val)
	{
		verboseFlag = val;
	}
	
	//
	// This routine determines the maximal timestep by finding the minimum of the
	// largest timestep allowed for the each of the input eigenvalues. 
	//
	// This routine only estimates the timestep for eigenvalues with negative real parts.
	// If there are no eigenvalues with negative real parts, the routine returns a timestep
	// of size 0.0.
	//
	double evaluateMaximalTimestep(double* lambdaReal, double* lambdaImag, long valueCount,
	double maximalDtBound)
	{
	double dtSize;
	double dtMax;
	long iStart = 0;
	
    long i; 
    
	if(verboseFlag == 1)
	{
	printf(" Estimated Dominant Eigenvalues and Maximal Timestep Allowed \n");
	}
	
	dtMax  = 0.0;
	for(i = 0; i < valueCount; i++)
	{
    	if(lambdaReal[i] < 0.0)
    	{
    	dtMax  = getMaximalTimestep(lambdaReal[i], lambdaImag[i],maximalDtBound);
        
        if(verboseFlag == 1)
		{
    	printf("%2ld Re: %12.5g Im: %12.5g Dt_Max: %12.5g \n",i+1,lambdaReal[i],lambdaImag[i],dtMax);
    	}
   	 	
    	iStart = i+1;
    	break;
    	}
	}
	for(i = iStart; i < valueCount; i++)
	{
		if(lambdaReal[i] < 0.0)
		{
		dtSize = getMaximalTimestep(lambdaReal[i], lambdaImag[i],maximalDtBound);
		
		if(verboseFlag == 1)
		{
    	printf("%2ld Re: %12.5g Im: %12.5g Dt_Max: %12.5g \n",i+1,lambdaReal[i],lambdaImag[i],dtSize);
    	}
   	 	
		dtMax  = (dtSize < dtMax) ? dtSize : dtMax;
		}
	}
	return dtMax;	
	}
	
	//
	// This routine scales the maximal timstep determined by the 
	// value 
	//
	// (1 - thetaReductionFactor*sin(theta))
	//
	// where theta = arg(lamdaReal + i*lambdaImag)
	//
	double evaluateMaximalTimestep(double thetaReductionFactor,double* lambdaReal, double* lambdaImag,
	long valueCount,double maximalDtBound)
	{
	double dtSize;
	double dtMax;
	long i; 
    long iStart = 0;
	
	if(verboseFlag == 1)
	{
	printf(" Estimated Dominant Eigenvalues and Maximal Timestep Allowed \n");
	}
	
	dtMax = maximalDtBound;
	
	for(i = 0; i < valueCount; i++)
	{
    	if(lambdaReal[i] < 0.0)
    	{
    	dtMax  = getMaximalTimestep(lambdaReal[i], lambdaImag[i],maximalDtBound);
    	dtMax *= 1.0 - thetaReductionFactor*(lambdaImag[i]/(lambdaReal[i]*lambdaReal[i]+ lambdaImag[i]*lambdaImag[i]));
        
        if(verboseFlag == 1)
		{
    	printf("%2ld Re: %12.5g Im: %12.5g Dt_Max: %12.5g \n",i+1,lambdaReal[i],lambdaImag[i],dtMax);
    	}
    	
    	iStart = i+1;
    	break;
    	}
	}
	for(i = iStart; i < valueCount; i++)
	{
		if(lambdaReal[i] < 0.0)
		{
		dtSize = getMaximalTimestep(lambdaReal[i], lambdaImag[i],maximalDtBound);
		dtSize *= 1.0 - thetaReductionFactor*(lambdaImag[i]/(lambdaReal[i]*lambdaReal[i]+ lambdaImag[i]*lambdaImag[i]));
		
		if(verboseFlag == 1)
		{
    	printf("%2ld Re: %12.5g Im: %12.5g Dt_Max: %12.5g \n",i+1,lambdaReal[i],lambdaImag[i],dtSize);
    	}
    	
		dtMax  = (dtSize < dtMax) ? dtSize : dtMax;

		}
	}
	return dtMax;	
	}
	
	double evaluateStabilityFactor(double lambdaReal,double lambdaImag)
	{
	return stabilityPoly(lambdaReal,lambdaImag);
	}
	
	double getMaximalTimestep(double lambdaReal, double lambdaImag,double maximalDtBound)
	{
	double zeroVal = 1.0e-12; // size we're calling a zero value
	long   n;                 // size of the stability bound table
	double timestep;
	double angleVar;
	double lambdaNorm;
	double s;
	//
	// Check for zero lambda -- return maximalDtBound
	//
	
	if((std::abs(lambdaReal) < zeroVal)&&(std::abs(lambdaImag) < zeroVal))
    {
    return maximalDtBound;
    }
    
	//
	// Check for Re(lambda)  > 0 return a sensible value to resolve the
	// solution behavior for a first order method. 
	//
    if(lambdaReal > 0.0) 
    {
    return 0.25/lambdaReal;
    }
//
//  Real case - use maximal element of stabilityBounds 
//
	n = angleVariable.size();
	
	if(std::abs(lambdaImag) < zeroVal)
	{
	return stabilityBound[n-1]/std::abs(lambdaReal);
	}
//
//  flip input value to quadrant II if required.
//
    if(lambdaImag < 0.0) lambdaImag *= -1.0;
    
	angleVar            = std::atan2(lambdaImag,lambdaReal)/(3.141592653589793238);
	lambdaNorm          = std::sqrt(lambdaReal*lambdaReal + lambdaImag*lambdaImag);
//
//  Use bisection to find table entries bracketing the angle value
//
    long a =  0;
    long b =  n-1;
    long c = (a+b)/2;
    int istop = 0;
    
    if(std::abs(angleVar - angleVariable[a]) < zeroVal) istop = 2;
    while(istop == 0)
    {
    if(angleVar < angleVariable[c])
    {
    	b = c;
    	c = (a+b)/2;
    	if(c == a) istop = 1;
    }
    else
    {
       a = c;
       c = (a+b)/2;
       if(c == a) istop = 1;
       if(std::abs(angleVar - angleVariable[a]) < zeroVal) istop = 2;
    }
    }
    
    if(istop == 2)
    {
    return stabilityBound[a]/lambdaNorm;
    }
    
    s = (angleVar - angleVariable[a])/(angleVariable[b]-angleVariable[a]);
    timestep = (stabilityBound[a] + s*(stabilityBound[b]-stabilityBound[a]))/lambdaNorm;

	return timestep; 
    }

    void getStabilityBoundaryData(std::vector<double>& xBdry, std::vector<double>& yBdry)
    {
    xBdry = zReal;
    yBdry = zReal;
    }
    
	void createStabilityBoundsTables()
	{
	long dimension      = 2;
    double dt;
    double dt0            = 0.1;
    double dt1            = (stageOrder*stageOrder*gamma)/200.0;
	double errorTolerance = 0.001;
    
	rkf45 < SRKstabilityContourODE > rkfMethod(stabilityContourODE,dimension);
	rkfMethod.setErrorTolerance(errorTolerance);
	
	
	zReal.clear();
    zImag.clear();
    
    angleVariable.clear();
    stabilityBound.clear();
//
//  Allocate data for ODE solver
//
	long i;
	
	double t;
    double* y = new double[2];
    
    t    = 0.0;
    y[0] = 0.0;
    y[1] = 0.0;
    
    int retVal    = 0;
    
    dt = dt0;
       
    zReal.push_back(y[0]);
    zImag.push_back(y[1]);
//
//  Call rkf45 method
//
    int iStop = 0;
    while(iStop != 1)
    {
    retVal = rkfMethod.advance(dt,t,y);
    if(retVal != 0)
    {
    std::cout  << " RKF Method Error Encountered In SRKtimestepEstimator" << std::endl;
    std::cout  << " Error Code : " << rkfMethod.getErrorCode() << std::endl;
    exit(0);
    }
    if(y[1] < 0.0) 
    {
    	iStop = 1;
    }
    
    zReal.push_back(y[0]);
    zImag.push_back(y[1]);
    
    if(-y[0] > 4.0) dt = dt1;
    }
//
//  Fix up the last computed point using linear interpolation
//
	long n     = zReal.size();
	
    double s   = -zImag[n-2]/(zImag[n-1]-zImag[n-2]);
    
    zReal[n-1] = zReal[n-2] + s*(zReal[n-1]-zReal[n-2]);
    zImag[n-1] = 0.0;
//
//  Create table for computation of timestep
//
    double angleVar;
    double angleSave;
    double zNorm;
    
    angleVar  = 0.5;
    angleSave = 0.5;
    angleVariable.push_back(angleVar);
    stabilityBound.push_back(0.0);
        	
    for(i = 1; i < n; i++)
    {
    	angleVar = atan2(zImag[i],zReal[i])/(3.141592653589793238);
    	if(angleVar > angleSave)
    	{
    	angleVariable.push_back(angleVar);
    	zNorm = sqrt(zReal[i]*zReal[i] + zImag[i]*zImag[i]);
    	stabilityBound.push_back(zNorm);
    	angleSave = angleVar;
    	}
    }
	//
	// clean up 
	//
	delete [] y;
	}
	
	int verboseFlag;
	
    long stageOrder;
    double gamma;
    
    std::vector<double> zReal;
    std::vector<double> zImag;
    
    std::vector<double>  angleVariable;
    std::vector<double> stabilityBound;
    
    SRKstabilityContourODE stabilityContourODE;
    RKStabilityPolynomialEvaluator stabilityPoly;
};

#endif

