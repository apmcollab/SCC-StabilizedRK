/*
 * PprimeCondition.h
 *
 *  Created on: Jul 13, 2020
 *      Author: anderson
 */

#include "RKpolynomialFunction.h"

#include <cfloat>
#include <cmath>

#ifndef P_PRIME_CONDITION_
#define P_PRIME_CONDITION_

class PprimeCondition
{
    public :

    PprimeCondition()
    {initialize();}

    PprimeCondition(RKpolynomialFunction& P, RKpolynomialFunction& Pprime, double M)
    {
    this->P.initialize(P);
    this->Pprime.initialize(Pprime);
    this->M    = M;
    }

    void initialize()
    {
    beta = 0.0;
    P.initialize();
    Pdelta.initialize();
    Pprime.initialize();
    PprimeDelta.initialize();
    M = 0.0;
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

    double beta;
    RKpolynomialFunction P;
    RKpolynomialFunction Pdelta;
    RKpolynomialFunction Pprime;
    RKpolynomialFunction PprimeDelta;

    double M;

//
//  getRoot is a localized version of Zeroin.f --
//  zeroin.f translated by f2c (version 19980403)and then
//  updated to C++ syntax.
//
//  Original Comments :
//
//  zeroin : a zero of the function  f(x) is computed in the interval ax,bx.
//
//  Input..
//
//  ax     left endpoint of initial interval
//  bx     right endpoint of initial interval
//  f      function subprogram which evaluates f(x) for any x in
//         the interval  ax,bx
//  tol    desired length of the interval of uncertainty of the
//         final result (.ge.0.)
//
//  Output..
//
//  zeroin abscissa approximating a zero of  f  in the interval ax,bx
//
//  It is assumed  that   f(ax)   and   f(bx)   have  opposite  signs
//  this is checked, and an error message is printed if this is not
//  satisfied.   zeroin  returns a zero  x  in the given interval
//  ax,bx  to within a tolerance  4*macheps*abs(x)+tol, where macheps  is
//  the  relative machine precision defined as the smallest representable
//  number such that  1.+macheps .gt. 1.
//
//  This function subprogram is a slightly  modified  translation  of
//  the algol 60 procedure  zero  given in  richard brent, algorithms for
//  minimization without derivatives, prentice-hall, inc. (1973).
//
double getRoot(double* ax, double* bx, double *tol)
{

    long c__4 = 4;
    double ret_val, d__1, d__2;
    double a, b, c__, d__, e, p, q, r__, s;
    double fa, fb, fc, xm, eps, tol1;

    eps = d1machForZeroin_(&c__4);
    tol1 = eps + 1.;

    a = *ax;
    b = *bx;
    fa = this->operator()(a);
    fb = this->operator()(b);

    if(std::abs(fa) < *tol) {return a;}
    if(std::abs(fb) < *tol) {return b;}
/*     check that f(ax) and f(bx) have different signs */
    if (fa == 0. || fb == 0.) {
	goto L20;
    }
    if (fa * (fb / std::abs(fb)) <= 0.) {
	goto L20;
    }

    printf("##########  Root Finder Error  ############# \n");
    printf("The function values at the initial interval \n");
    printf("endpoionts are not of opposite sign. \n");
    printf("The root value returned is inaccurate.\n\n");

    ret_val = 0;
    return ret_val;
L20:
    c__ = a;
    fc = fa;
    d__ = b - a;
    e = d__;
L30:
    if (std::abs(fc) >= std::abs(fb)) {
	goto L40;
    }
    a = b;
    b = c__;
    c__ = a;
    fa = fb;
    fb = fc;
    fc = fa;
L40:
    tol1 = eps * 2. * std::abs(b) + *tol * .5;
    xm = (c__ - b) * .5;
    if (std::abs(xm) <= tol1 || fb == 0.) {
	goto L150;
    }

/* see if a bisection is forced */

    if (std::abs(e) >= tol1 && std::abs(fa) > std::abs(fb)) {
	goto L50;
    }
    d__ = xm;
    e = d__;
    goto L110;
L50:
    s = fb / fa;
    if (a != c__) {
	goto L60;
    }

/* linear interpolation */

    p = xm * 2. * s;
    q = 1. - s;
    goto L70;

/* inverse quadratic interpolation */

L60:
    q = fa / fc;
    r__ = fb / fc;
    p = s * (xm * 2. * q * (q - r__) - (b - a) * (r__ - 1.));
    q = (q - 1.) * (r__ - 1.) * (s - 1.);
L70:
    if (p <= 0.) {
	goto L80;
    }
    q = -q;
    goto L90;
L80:
    p = -p;
L90:
    s = e;
    e = d__;
    if (p * 2. >= xm * 3. * q - (d__1 = tol1 * q, std::abs(d__1)) || p >= (d__2 =
	    s * .5 * q, std::abs(d__2))) {
	goto L100;
    }
    d__ = p / q;
    goto L110;
L100:
    d__ = xm;
    e = d__;
L110:
    a = b;
    fa = fb;
    if (std::abs(d__) <= tol1) {
	goto L120;
    }
    b += d__;
    goto L140;
L120:
    if (xm <= 0.) {
	goto L130;
    }
    b += tol1;
    goto L140;
L130:
    b -= tol1;
L140:
    fb = this->operator()(b);
    if (fb * (fc / std::abs(fc)) > 0.) {
	goto L20;
    }
    goto L30;
L150:
    ret_val = b;
    return ret_val;
} /* zeroin_ */

/*
   d1machForZeroin is a localized version of d1mach -- d1mach
   translated by f2c (version 19980403) and then updated to
   C++ syntax.
*/
double d1machForZeroin_(long *i)
{
	switch(*i){
 	  case 1: return DBL_MIN;
 	  case 2: return DBL_MAX;
	  case 3: return DBL_EPSILON/FLT_RADIX;
 	  case 4: return DBL_EPSILON;
 	  case 5: return std::log10(double(FLT_RADIX));
	  }
 	fprintf(stderr, "invalid argument: d1mach(%ld)\n", *i);
 	exit(1); return 0; /*/+ for compilers that complain of missing return values +/ */
}
};

#endif
