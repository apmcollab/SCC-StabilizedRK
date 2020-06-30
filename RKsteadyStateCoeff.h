//
// RKsteadyStateCoeff.h 
//
// Author: Chris Anderson  
// (C) UCLA 2012 
//
//
// ############################################################
// RKsteadyStateCoeff.h 
//
// Version : May 26, 2004 
//
// Chris Anderson UCLA 
//###############################################################
//

#include <iostream>
#include <iomanip>
using namespace std;


#ifndef __RKsteadyStateCoeff__
#define __RKsteadyStateCoeff__

class RKpolynomialFunction 
{
    public :
//
//  Constructors
//
    RKpolynomialFunction();
    RKpolynomialFunction(const RKpolynomialFunction& P);
    RKpolynomialFunction(int n);
    RKpolynomialFunction(int n, double coefficients[]);
//
//  Initialize/Resize
//
    void initialize();
    void initialize(const RKpolynomialFunction& P);
    void initialize(int n);
    void initialize(int n, double coefficients[]);
    void resizeTo(int n);

    virtual ~RKpolynomialFunction();

    RKpolynomialFunction   differentiate();

    RKpolynomialFunction operator=(const RKpolynomialFunction& P);
    RKpolynomialFunction operator+(const RKpolynomialFunction& P);
    RKpolynomialFunction operator-(const RKpolynomialFunction& P);
    RKpolynomialFunction operator*(const RKpolynomialFunction& P);
    RKpolynomialFunction operator-() const;

    RKpolynomialFunction operator*(double alpha);
    friend RKpolynomialFunction operator*(double alpha, const RKpolynomialFunction& P);

    RKpolynomialFunction operator/(double alpha);

    virtual double operator()(double x);

    const double& operator[](long i) const;
    double&  operator[](long i);

    RKpolynomialFunction shift(double p);
    RKpolynomialFunction scale(double alpha);
    int getDegree() const;
    double* getCoefficients() const;
//
//  Input/Output
//
    friend ostream&  
    operator <<(ostream& outStream, const RKpolynomialFunction& A);
//
//  Class data
//
    public  :

    double*           a;           // ptr to coefficient array
    int          degree;
    static const double zero;
};


class PprimeCondition
{
    public : 

    PprimeCondition()
    {}

    PprimeCondition(RKpolynomialFunction& P, RKpolynomialFunction& Pprime, double M)
    {
    this->P.initialize(P);
    this->Pprime.initialize(Pprime);
    this->M    = M;
    }

    void initialize(RKpolynomialFunction& P, RKpolynomialFunction& Pprime, double M)
    {
    this->P.initialize(P);
    this->Pprime.initialize(Pprime);
    this->M    = M;
    }

    double operator()(double delta) 
    {
    Pdelta = P.scale(2.0/(M-delta));
    Pdelta = Pdelta.shift(-((2.0*M)/(M-delta) - 1)*((M-delta)/2.0));
    beta   = Pdelta(0.0);
    

    PprimeDelta = Pprime.scale(2.0/(M-delta));
    PprimeDelta = PprimeDelta.shift(-((2.0*M)/(M-delta) - 1)*((M-delta)/2.0));
    PprimeDelta = PprimeDelta*(2.0/(M-delta));
    PprimeDelta = PprimeDelta/beta;

    return PprimeDelta(0.0) - 1.0;
    }

    double getRoot(double* ax, double* bx, double *tol);

    double beta;
    RKpolynomialFunction P;
    RKpolynomialFunction Pdelta;
    RKpolynomialFunction Pprime;
    RKpolynomialFunction PprimeDelta;
    
    double M;

    private : 

    double d1machForZeroin_(long *i);
};


class RKcoefficients
{
    public :

    RKcoefficients();
    RKcoefficients(long stageOrder, double gamma);

    ~RKcoefficients();

    long     stageOrder;
    double        gamma;
    double** alphaCoeff;
};

class RKsteadyStateCoeff
{
    public:

    static void getRKcoefficients(long stageOrder, double gamma, 
    double** alphaCoeff);

    static double** getRKcoefficientsPtr(long stageOrder, 
    double gamma);

    private:

    static void getTchebyShiftFactors(double M, double gamma, double& delta, 
    double& beta,RKpolynomialFunction& Cheby);

    static void backSolve(long N, double* c, double** A);
    static void create2Darray(double**& a, long m, long n);
    static void delete2Darray(double**& a); 

    static RKcoefficients** coeffCache;
    static long cacheSize;
    static long cacheStorageSize;
    static long cacheStorageIncrement;
    static void expandCache(long storageIncrement);
};


#endif
 
