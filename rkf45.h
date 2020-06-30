//
// rkf45.h 
//
// Author: Chris Anderson  
// (C) UCLA 2012 
//

#include <cmath>

#include "stdlib.h"
#include "stdio.h"


#ifndef _rkf45_
#define _rkf45_

/* Transation of rkf45.f  by f2c and then refactored into a
   C++ templated class
   
*/

#undef abs
#undef min
#undef max
#undef TRUE_  
#undef FALSE_  

#define abs(x) ((x) >= 0 ? (x) : -(x))
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define TRUE_ (1)
#define FALSE_ (0)

template < class ODEoperator >  class rkf45
{
public:
	
rkf45(ODEoperator& f, long dimension)
{
	this->fun             = &f;
	this->dimension       = dimension;
	this->errorTolerance  = 1.0e-03;
	this->errorCode       = 0;
	this->evaluationCount = 0;
	
    work  = new double[4 + 6*dimension];
    iwork = new long[10];
	
	c_b17 = .2; // global constant used by translated routine
	c__2   = 2;	// global constant used by translated routine
}
	

~rkf45()
{
	if(work != 0)  delete [] work;
	if(iwork != 0) delete [] iwork;
}

void setErrorTolerance(double tolerance)
{
    errorTolerance  = tolerance ;
}

double getErrorTolerance() const
{
    return errorTolerance;
}

int getErrorCode() const
{
    return errorCode;
}

long getEvaluationCount() const
{
	return evaluationCount;
}

/**
@arg dt     : timestep
@arg t      : current time
@arg y      : double array containing solution at time t
\returns     
y is updated with the solution value at t + dt. The return value is 0 if advance is successful; 
otherwise the return value is the error code iflag from the rkf45 method. The value of 
iflag can also be obtained using the getErrorCode() member function. 
The following description is taken from the rkf45.f source code file. <p>

<pre>
       iflag = 3 -- integration was not completed because relative error
                    tolerance was too small. relerr has been increased
                    appropriately for continuing.
             = 4 -- integration was not completed because more than
                    3000 derivative evaluations were needed. this
                    is approximately 500 steps.
             = 5 -- integration was not completed because solution
                    vanished making a pure relative error test
                    impossible. must use non-zero abserr to continue.
                    using the one-step integration mode for one step
                    is a good way to proceed.
             = 6 -- integration was not completed because requested
                    accuracy could not be achieved using smallest
                    allowable stepsize. user must increase the error
                    tolerance before continued integration can be
                   attempted.
             = 7 -- it is likely that rkf45 is inefficient for solving
                    this problem. too much output is restricting the
                    natural stepsize choice. use the one-step integrator
                    mode.
             = 8 -- invalid input parameters
                    this indicator occurs if any of the following is
                    satisfied -   neqn .le. 0
                                  t=tout  and  iflag .ne. +1 or -1
                                  relerr or abserr .lt. 0.
                                  iflag .eq. 0  or  .lt. -2  or  .gt. 8
</pre>

The errorTolerance specified is used to set both the absolute and 
relative error tolerence values. As dicussed in the comments in
original rkf45.f, <a href="../rkf45Comments.f">rkf45Comments.f</a>,
the advance() routine should not be used with relative error control 
smaller than about 1.e-8.
*/
	
int advance(double dt, double&t, double* y)
{
//
//  Localize parameters
//
    double tout   = t +  dt; 
    double relerr = errorTolerance;
    double abserr = errorTolerance;
    long   iflag  = 1;
    errorCode     = 0;
    
    
    rkf45_(&dimension,y,&t,&tout,&relerr,&abserr,&iflag,work,iwork);
    
    evaluationCount = iwork[0];
    
    if(iflag != 2)
    {
    errorCode = iflag;
    return 1;
    }

    return 0;
}

	ODEoperator* fun;	
	long   dimension;
	
	double errorTolerance;
	long   errorCode;
	long   evaluationCount;
	

    double* work;
    long*  iwork;
//
// global constants used by translated routine 
//
	double     c_b17;
	long        c__2;	

//
// Required f2c supporting f2c library routines translated 
// into C++ equivalents
//

void s_stop(const char *s, int n)
{
	int i;
	if(n > 0)
	{
	fprintf(stderr, "STOP ");
	for(i = 0; i<n ; ++i)
	putc(*s++, stderr);
	fprintf(stderr, " statement executed\n");
	}
	exit(0);
}


double pow_dd(double *ap, double *bp)
{
	return(pow(*ap, *bp));
}

long i_sign(long *a, long *b)
{
	long x;
	x = (*a >= 0 ? *a : - *a);
	return( *b >= 0 ? x : -x);
}

double d_sign(double *a, double *b)
{
	double x;
	x = (*a >= 0 ? *a : - *a);
	return( *b >= 0 ? x : -x);
}


int rkf45_(long *neqn, double *y, double*t, double *tout, double *relerr, double *abserr, 
long*iflag, double *work, long *iwork)
{
    static long k1, k2, k3, k4, k5, k6, k1m;

/*     fehlberg fourth-fifth order runge-kutta method */

/*     written by h.a.watts and l.f.shampine */
/*                   sandia laboratories */
/*                  albuquerque,new mexico */

/*    rkf45 is primarily designed to solve non-stiff and mildly stiff */
/*    differential equations when derivative evaluations are inexpensive. */
/*    rkf45 should generally not be used when the user is demanding */
/*    high accuracy. */

/* abstract */

/*    subroutine  rkf45  integrates a system of neqn first order */
/*    ordinary differential equations of the form */
/*             dy(i)/dt = f(t,y(1),y(2),...,y(neqn)) */
/*              where the y(i) are given at t . */
/*    typically the subroutine is used to integrate from t to tout but it */
/*    can be used as a one-step integrator to advance the solution a */
/*    single step in the direction of tout.  on return the parameters in */
/*    the call list are set for continuing the integration. the user has */
/*    only to call rkf45 again (and perhaps define a new value for tout). */
/*    actually, rkf45 is an interfacing routine which calls subroutine */
/*    rkfs for the solution.  rkfs in turn calls subroutine  fehl which */
/*    computes an approximate solution over one step. */

/*    rkf45  uses the runge-kutta-fehlberg (4,5)  method described */
/*    in the reference */
/*    e.fehlberg , low-order classical runge-kutta formulas with stepsize */
/*                 control , nasa tr r-315 */

/*    the performance of rkf45 is illustrated in the reference */
/*    l.f.shampine,h.a.watts,s.davenport, solving non-stiff ordinary */
/*                 differential equations-the state of the art , */
/*                 sandia laboratories report sand75-0182 , */
/*                 to appear in siam review. */


/*    the parameters represent- */
/*      f -- subroutine f(t,y,yp) to evaluate derivatives yp(i)=dy(i)/dt */
/*      neqn -- number of equations to be integrated */
/*      y(*) -- solution vector at t */
/*      t -- independent variable */
/*      tout -- output point at which solution is desired */
/*      relerr,abserr -- relative and absolute error tolerances for local */
/*            error test. at each step the code requires that */
/*                 abs(local error) .le. relerr*abs(y) + abserr */
/*            for each component of the local error and solution vectors */
/*      iflag -- indicator for status of integration */
/*      work(*) -- array to hold information internal to rkf45 which is */
/*            necessary for subsequent calls. must be dimensioned */
/*            at least  3+6*neqn */
/*      iwork(*) -- long array used to hold information internal to */
/*            rkf45 which is necessary for subsequent calls. must be */
/*            dimensioned at least  5 */


/*  first call to rkf45 */

/*    the user must provide storage in his calling program for the arrays */
/*    in the call list  -      y(neqn) , work(3+6*neqn) , iwork(5)  , */
/*    declare f in an external statement, supply subroutine f(t,y,yp) and */
/*    initialize the following parameters- */

/*      neqn -- number of equations to be integrated.  (neqn .ge. 1) */
/*      y(*) -- vector of initial conditions */
/*      t -- starting point of integration , must be a variable */
/*      tout -- output point at which solution is desired. */
/*            t=tout is allowed on the first call only, in which case */
/*            rkf45 returns with iflag=2 if continuation is possible. */
/*      relerr,abserr -- relative and absolute local error tolerances */
/*            which must be non-negative. relerr must be a variable while */
/*            abserr may be a constant. the code should normally not be */
/*            used with relative error control smaller than about 1.e-8 . */
/*            to avoid limiting precision difficulties the code requires */
/*            relerr to be larger than an internally computed relative */
/*            error parameter which is machine dependent. in particular, */
/*            pure absolute error is not permitted. if a smaller than */
/*            allowable value of relerr is attempted, rkf45 increases */
/*            relerr appropriately and returns control to the user before */
/*            continuing the integration. */
/*      iflag -- +1,-1  indicator to initialize the code for each new */
/*            problem. normal input is +1. the user should set iflag=-1 */
/*            only when one-step integrator control is essential. in this */
/*            case, rkf45 attempts to advance the solution a single step */
/*            in the direction of tout each time it is called. since this */
/*            mode of operation results in extra computing overhead, it */
/*            should be avoided unless needed. */


/*  output from rkf45 */

/*      y(*) -- solution at t */
/*      t -- last point reached in integration. */
/*      iflag = 2 -- integration reached tout. indicates successful retur */
/*                   and is the normal mode for continuing integration. */
/*            =-2 -- a single successful step in the direction of tout */
/*                   has been taken. normal mode for continuing */
/*                   integration one step at a time. */
/*            = 3 -- integration was not completed because relative error */
/*                   tolerance was too small. relerr has been increased */
/*                   appropriately for continuing. */
/*            = 4 -- integration was not completed because more than */
/*                   3000 derivative evaluations were needed. this */
/*                   is approximately 500 steps. */
/*            = 5 -- integration was not completed because solution */
/*                   vanished making a pure relative error test */
/*                   impossible. must use non-zero abserr to continue. */
/*                   using the one-step integration mode for one step */
/*                   is a good way to proceed. */
/*            = 6 -- integration was not completed because requested */
/*                   accuracy could not be achieved using smallest */
/*                   allowable stepsize. user must increase the error */
/*                   tolerance before continued integration can be */
/*                   attempted. */
/*            = 7 -- it is likely that rkf45 is inefficient for solving */
/*                   this problem. too much output is restricting the */
/*                   natural stepsize choice. use the one-step integrator */
/*                   mode. */
/*            = 8 -- invalid input parameters */
/*                   this indicator occurs if any of the following is */
/*                   satisfied -   neqn .le. 0 */
/*                                 t=tout  and  iflag .ne. +1 or -1 */
/*                                 relerr or abserr .lt. 0. */
/*                                 iflag .eq. 0  or  .lt. -2  or  .gt. 8 */
/*      work(*),iwork(*) -- information which is usually of no interest */
/*                   to the user but necessary for subsequent calls. */
/*                   work(1),...,work(neqn) contain the first derivatives */
/*                   of the solution vector y at t. work(neqn+1) contains */
/*                   the stepsize h to be attempted on the next step. */
/*                   iwork(1) contains the derivative evaluation counter. */


/*  subsequent calls to rkf45 */

/*    subroutine rkf45 returns with all information needed to continue */
/*    the integration. if the integration reached tout, the user need onl */
/*    define a new tout and call rkf45 again. in the one-step integrator */
/*    mode (iflag=-2) the user must keep in mind that each step taken is */
/*    in the direction of the current tout. upon reaching tout (indicated */
/*    by changing iflag to 2),the user must then define a new tout and */
/*    reset iflag to -2 to continue in the one-step integrator mode. */

/*    if the integration was not completed but the user still wants to */
/*    continue (iflag=3,4 cases), he just calls rkf45 again. with iflag=3 */
/*    the relerr parameter has been adjusted appropriately for continuing */
/*    the integration. in the case of iflag=4 the function counter will */
/*    be reset to 0 and another 3000 function evaluations are allowed. */

/*    however,in the case iflag=5, the user must first alter the error */
/*    criterion to use a positive value of abserr before integration can */
/*    proceed. if he does not,execution is terminated. */

/*    also,in the case iflag=6, it is necessary for the user to reset */
/*    iflag to 2 (or -2 when the one-step integration mode is being used) */
/*    as well as increasing either abserr,relerr or both before the */
/*    integration can be continued. if this is not done, execution will */
/*    be terminated. the occurrence of iflag=6 indicates a trouble spot */
/*    (solution is changing rapidly,singularity may be present) and it */
/*    often is inadvisable to continue. */

/*    if iflag=7 is encountered, the user should use the one-step */
/*    integration mode with the stepsize determined by the code or */
/*    consider switching to the adams codes de/step,intrp. if the user */
/*    insists upon continuing the integration with rkf45, he must reset */
/*    iflag to 2 before calling rkf45 again. otherwise,execution will be */
/*    terminated. */

/*    if iflag=8 is obtained, integration can not be continued unless */
/*    the invalid input parameters are corrected. */

/*    it should be noted that the arrays work,iwork contain information */
/*    required for subsequent integration. accordingly, work and iwork */
/*    should not be altered. */


/*     if compiler checks subscripts, change work(1) to work(3+6*neqn) */

/*     compute indices for the splitting of the work array */

    /* Parameter adjustments */
    --y;
    --work;
    --iwork;

    /* Function Body */
    k1m = *neqn + 1;
    k1 = k1m + 1;
    k2 = k1 + *neqn;
    k3 = k2 + *neqn;
    k4 = k3 + *neqn;
    k5 = k4 + *neqn;
    k6 = k5 + *neqn;

/*     this interfacing routine merely relieves the user of a long */
/*     calling list via the splitting apart of two working storage */
/*     arrays. if this is not compatible with the users compiler, */
/*     he must use rkfs directly. */

    rkfs_( neqn, &y[1], t, tout, relerr, abserr, iflag, &work[1], &
	    work[k1m], &work[k1], &work[k2], &work[k3], &work[k4], &work[k5], 
	    &work[k6], &work[k6 + 1], &iwork[1], &iwork[2], &iwork[3], &iwork[
	    4], &iwork[5]);

    return 0;
} /* rkf45_ */

int rkfs_(long *neqn, double *y, double *
	t, double *tout, double *relerr, double *abserr, long *
	iflag, double *yp, double *h__, double *f1, double *
	f2, double *f3, double *f4, double *f5, double *savre,
	 double *savae, long *nfe, long *kop, long *init, 
	long *jflag, long *kflag)
{
    /* Initialized data */

    static double remin = 1e-12;
    static long maxnfe = 3000;

    /* System generated locals */
    long i__1;
    double d__1, d__2, d__3, d__4;

    static double hmin, toln, epsp1, a;
    static long k;
    static double s, scale;
    static long mflag;
    static double eeoet, ae, ee;
    static int hfaild;
    static double dt, et, u26, esttol;
    static int output;
    static double eps, rer, tol, ypk;


/*     fehlberg fourth-fifth order runge-kutta method */


/*     rkfs integrates a system of first order ordinary differential */
/*     equations as described in the comments for rkf45 . */
/*     the arrays yp,f1,f2,f3,f4,and f5 (of dimension at least neqn) and */
/*     the variables h,savre,savae,nfe,kop,init,jflag,and kflag are used */
/*     internally by the code and appear in the call list to eliminate */
/*     local retention of variables between calls. accordingly, they */
/*     should not be altered. items of possible interest are */
/*         yp - derivative of solution vector at t */
/*         h  - an appropriate stepsize to be used for the next step */
/*         nfe- counter on the number of derivative function evaluations */
/*     remin is the minimum acceptable value of relerr.  attempts */
/*     to obtain higher accuracy with this subroutine are usually */
/*     very expensive and often unsuccessful. */

    /* Parameter adjustments */
    --f5;
    --f4;
    --f3;
    --f2;
    --f1;
    --yp;
    --y;

    /* Function Body */


/*     the expense is controlled by restricting the number */
/*     of function evaluations to be approximately maxnfe. */
/*     as set, this corresponds to about 500 steps. */



/*     check input parameters */


    if (*neqn < 1) {
	goto L10;
    }
    if (*relerr < 0. || *abserr < 0.) {
	goto L10;
    }
    mflag = abs(*iflag);
    if (mflag == 0 || mflag > 8) {
	goto L10;
    }
    if (mflag != 1) {
	goto L20;
    }

/*     first call, compute machine epsilon */

    eps = 1.;
L5:
    eps /= 2.;
    epsp1 = eps + 1.;
    if (epsp1 > 1.) {
	goto L5;
    }
    u26 = eps * 26.;
    goto L50;

/*     invalid input */
L10:
    *iflag = 8;
    return 0;

/*     check continuation possibilities */

L20:
    if (*t == *tout && *kflag != 3) {
	goto L10;
    }
    if (mflag != 2) {
	goto L25;
    }

/*     iflag = +2 or -2 */
    if (*kflag == 3 || *init == 0) {
	goto L45;
    }
    if (*kflag == 4) {
	goto L40;
    }
    if (*kflag == 5 && *abserr == 0.) {
	goto L30;
    }
    if (*kflag == 6 && *relerr <= *savre && *abserr <= *savae) {
	goto L30;
    }
    goto L50;

/*     iflag = 3,4,5,6,7 or 8 */
L25:
    if (*iflag == 3) {
	goto L45;
    }
    if (*iflag == 4) {
	goto L40;
    }
    if (*iflag == 5 && *abserr > 0.) {
	goto L45;
    }

/*     integration cannot be continued since user did not respond to */
/*     the instructions pertaining to iflag=5,6,7 or 8 */
L30:
    s_stop("", (int)0);

/*     reset function evaluation counter */
L40:
    *nfe = 0;
    if (mflag == 2) {
	goto L50;
    }

/*     reset flag value from previous call */
L45:
    *iflag = *jflag;
    if (*kflag == 3) {
	mflag = abs(*iflag);
    }

/*     save input iflag and set continuation flag value for subsequent */
/*     input checking */
L50:
    *jflag = *iflag;
    *kflag = 0;

/*     save relerr and abserr for checking input on subsequent calls */
    *savre = *relerr;
    *savae = *abserr;

/*     restrict relative error tolerance to be at least as large as */
/*     2*eps+remin to avoid limiting precision difficulties arising */
/*     from impossible accuracy requests */

    rer = eps * 2. + remin;
    if (*relerr >= rer) {
	goto L55;
    }

/*     relative error tolerance too small */
    *relerr = rer;
    *iflag = 3;
    *kflag = 3;
    return 0;

L55:
    dt = *tout - *t;

    if (mflag == 1) {
	goto L60;
    }
    if (*init == 0) {
	goto L65;
    }
    goto L80;

/*     initialization -- */
/*                       set initialization completion indicator,init */
/*                       set indicator for too many output points,kop */
/*                       evaluate initial derivatives */
/*                       set counter for function evaluations,nfe */
/*                       estimate starting stepsize */

L60:
    *init = 0;
    *kop = 0;

    a = *t;
    fun->operator()(&a, &y[1], &yp[1]);
    *nfe = 1;
    if (*t != *tout) {
	goto L65;
    }
    *iflag = 2;
    return 0;


L65:
    *init = 1;
    *h__ = abs(dt);
    toln = (float)0.;
    i__1 = *neqn;
    for (k = 1; k <= i__1; ++k) {
	tol = *relerr * (d__1 = y[k], abs(d__1)) + *abserr;
	if (tol <= (float)0.) {
	    goto L70;
	}
	toln = tol;
	ypk = (d__1 = yp[k], abs(d__1));
/* Computing 5th power */
	d__1 = *h__, d__2 = d__1, d__1 *= d__1;
	if (ypk * (d__2 * (d__1 * d__1)) > tol) {
	    d__3 = tol / ypk;
	    *h__ = pow_dd(&d__3, &c_b17);
	}
L70:
	;
    }
    if (toln <= 0.) {
	*h__ = 0.;
    }
/* Computing MAX */
/* Computing MAX */
    d__3 = abs(*t), d__4 = abs(dt);
    d__1 = *h__, d__2 = u26 *max(d__3,d__4);
    *h__ = max(d__1,d__2);
    *jflag = i_sign(&c__2, iflag);


/*     set stepsize for integration in the direction from t to tout */

L80:
    *h__ = d_sign(h__, &dt);

/*     test to see if rkf45 is being severely impacted by too many */
/*     output points */

    if (abs(*h__) >= abs(dt) * 2.) {
	++(*kop);
    }
    if (*kop != 100) {
	goto L85;
    }

/*     unnecessary frequency of output */
    *kop = 0;
    *iflag = 7;
    return 0;

L85:
    if (abs(dt) > u26 * abs(*t)) {
	goto L95;
    }

/*     if too close to output point,extrapolate and return */

    i__1 = *neqn;
    for (k = 1; k <= i__1; ++k) {
/* L90: */
	y[k] += dt * yp[k];
    }
    a = *tout;
    fun->operator()(&a, &y[1], &yp[1]);
    ++(*nfe);
    goto L300;


/*     initialize output point indicator */

L95:
    output = FALSE_;

/*     to avoid premature underflow in the error tolerance function, */
/*     scale the error tolerances */

    scale = 2. / *relerr;
    ae = scale * *abserr;


/*     step by step integration */

L100:
    hfaild = FALSE_;

/*     set smallest allowable stepsize */

    hmin = u26 * abs(*t);

/*     adjust stepsize if necessary to hit the output point. */
/*     look ahead two steps to avoid drastic changes in the stepsize and */
/*     thus lessen the impact of output points on the code. */

    dt = *tout - *t;
    if (abs(dt) >= abs(*h__) * 2.) {
	goto L200;
    }
    if (abs(dt) > abs(*h__)) {
	goto L150;
    }

/*     the next successful step will complete the integration to the */
/*     output point */

    output = TRUE_;
    *h__ = dt;
    goto L200;

L150:
    *h__ = dt * .5;



/*     core integrator for taking a single step */

/*     the tolerances have been scaled to avoid premature underflow in */
/*     computing the error tolerance function et. */
/*     to avoid problems with zero crossings,relative error is measured */
/*     using the average of the magnitudes of the solution at the */
/*     beginning and end of a step. */
/*     the error estimate formula has been grouped to control loss of */
/*     significance. */
/*     to distinguish the various arguments, h is not permitted */
/*     to become smaller than 26 units of roundoff in t. */
/*     practical limits on the change in the stepsize are enforced to */
/*     smooth the stepsize selection process and to avoid excessive */
/*     chattering on problems having discontinuities. */
/*     to prevent unnecessary failures, the code uses 9/10 the stepsize */
/*     it estimates will succeed. */
/*     after a step failure, the stepsize is not allowed to increase for */
/*     the next attempted step. this makes the code more efficient on */
/*     problems having discontinuities and more effective in general */
/*     since local extrapolation is being used and extra caution seems */
/*     warranted. */


/*     test number of derivative function evaluations. */
/*     if okay,try to advance the integration from t to t+h */

L200:
    if (*nfe <= maxnfe) {
	goto L220;
    }

/*     too much work */
    *iflag = 4;
    *kflag = 4;
    return 0;

/*     advance an approximate solution over one step of length h */

L220:
    fehl_(neqn, &y[1], t, h__, &yp[1], &f1[1], &f2[1], &f3[1], &f4[1]
	    , &f5[1], &f1[1]);
    *nfe += 5;

/*     compute and test allowable tolerances versus local error estimates */
/*     and remove scaling of tolerances. note that relative error is */
/*     measured with respect to the average of the magnitudes of the */
/*     solution at the beginning and end of the step. */

    eeoet = 0.;
    i__1 = *neqn;
    for (k = 1; k <= i__1; ++k) {
	et = (d__1 = y[k], abs(d__1)) + (d__2 = f1[k], abs(d__2)) + ae;
	if (et > 0.) {
	    goto L240;
	}

/*       inappropriate error tolerance */
	*iflag = 5;
	return 0;

L240:
	ee = (d__1 = yp[k] * -2090. + (f3[k] * 21970. - f4[k] * 15048.) + (f2[
		k] * 22528. - f5[k] * 27360.), abs(d__1));
/* L250: */
/* Computing MAX */
	d__1 = eeoet, d__2 = ee / et;
	eeoet = max(d__1,d__2);
    }

    esttol = abs(*h__) * eeoet * scale / 752400.;

    if (esttol <= 1.) {
	goto L260;
    }


/*     unsuccessful step */
/*                       reduce the stepsize , try again */
/*                       the decrease is limited to a factor of 1/10 */

    hfaild = TRUE_;
    output = FALSE_;
    s = .1;
    if (esttol < 59049.) {
	s = .9 / pow_dd(&esttol, &c_b17);
    }
    *h__ = s * *h__;
    if (abs(*h__) > hmin) {
	goto L200;
    }

/*     requested error unattainable at smallest allowable stepsize */
    *iflag = 6;
    *kflag = 6;
    return 0;


/*     successful step */
/*                        store solution at t+h */
/*                        and evaluate derivatives there */

L260:
    *t += *h__;
    i__1 = *neqn;
    for (k = 1; k <= i__1; ++k) {
/* L270: */
	y[k] = f1[k];
    }
    a = *t;
    fun->operator()(&a, &y[1], &yp[1]);
    ++(*nfe);


/*                       choose next stepsize */
/*                       the increase is limited to a factor of 5 */
/*                       if step failure has just occurred, next */
/*                          stepsize is not allowed to increase */

    s = 5.;
    if (esttol > 1.889568e-4) {
	s = .9 / pow_dd(&esttol, &c_b17);
    }
    if (hfaild) {
	s = min(s,1.);
    }
/* Computing MAX */
    d__2 = s * abs(*h__);
    d__1 = max(d__2,hmin);
    *h__ = d_sign(&d__1, h__);

/*     end of core integrator */


/*     should we take another step */

    if (output) {
	goto L300;
    }
    if (*iflag > 0) {
	goto L100;
    }


/*     integration successfully completed */

/*     one-step mode */
    *iflag = -2;
    return 0;

/*     interval mode */
L300:
    *t = *tout;
    *iflag = 2;
    return 0;

} /* rkfs_ */

int fehl_(long *neqn, double *y, double *t, double *h__, double *yp, double *f1, double *f2, 
double *f3, double *f4, double *f5, double *s)
{
    /* System generated locals */
    long i__1;
    double d__1;

    /* Local variables */
    static long k;
    static double ch;


/*     fehlberg fourth-fifth order runge-kutta method */

/*    fehl integrates a system of neqn first order */
/*    ordinary differential equations of the form */
/*             dy(i)/dt=f(t,y(1),---,y(neqn)) */
/*    where the initial values y(i) and the initial derivatives */
/*    yp(i) are specified at the starting point t. fehl advances */
/*    the solution over the fixed step h and returns */
/*    the fifth order (sixth order accurate locally) solution */
/*    approximation at t+h in array s(i). */
/*    f1,---,f5 are arrays of dimension neqn which are needed */
/*    for internal storage. */
/*    the formulas have been grouped to control loss of significance. */
/*    fehl should be called with an h not smaller than 13 units of */
/*    roundoff in t so that the various independent arguments can be */
/*    distinguished. */




    /* Parameter adjustments */
    --s;
    --f5;
    --f4;
    --f3;
    --f2;
    --f1;
    --yp;
    --y;

    /* Function Body */
    ch = *h__ / 4.;
    i__1 = *neqn;
    for (k = 1; k <= i__1; ++k) {
/* L221: */
	f5[k] = y[k] + ch * yp[k];
    }
    d__1 = *t + ch;
    fun->operator()(&d__1, &f5[1], &f1[1]);

    ch = *h__ * 3. / 32.;
    i__1 = *neqn;
    for (k = 1; k <= i__1; ++k) {
/* L222: */
	f5[k] = y[k] + ch * (yp[k] + f1[k] * 3.);
    }
    d__1 = *t + *h__ * 3. / 8.;
    fun->operator()(&d__1, &f5[1], &f2[1]);

    ch = *h__ / 2197.;
    i__1 = *neqn;
    for (k = 1; k <= i__1; ++k) {
/* L223: */
	f5[k] = y[k] + ch * (yp[k] * 1932. + (f2[k] * 7296. - f1[k] * 7200.));
    }
    d__1 = *t + *h__ * 12. / 13.;
    fun->operator()(&d__1, &f5[1], &f3[1]);

    ch = *h__ / 4104.;
    i__1 = *neqn;
    for (k = 1; k <= i__1; ++k) {
/* L224: */
	f5[k] = y[k] + ch * (yp[k] * 8341. - f3[k] * 845. + (f2[k] * 29440. - 
		f1[k] * 32832.));
    }
    d__1 = *t + *h__;
    fun->operator()(&d__1, &f5[1], &f4[1]);

    ch = *h__ / 20520.;
    i__1 = *neqn;
    for (k = 1; k <= i__1; ++k) {
/* L225: */
	f1[k] = y[k] + ch * (yp[k] * -6080. + (f3[k] * 9295. - f4[k] * 5643.) 
		+ (f1[k] * 41040. - f2[k] * 28352.));
    }
    d__1 = *t + *h__ / 2.;
    fun->operator()(&d__1, &f1[1], &f5[1]);

/*     compute approximate solution at t+h */

    ch = *h__ / 7618050.;
    i__1 = *neqn;
    for (k = 1; k <= i__1; ++k) {
/* L230: */
	s[k] = y[k] + ch * (yp[k] * 902880. + (f3[k] * 3855735. - f4[k] * 
		1371249.) + (f2[k] * 3953664. + f5[k] * 277020.));
    }

    return 0;
} /* fehl_ */


};


#undef abs
#undef min
#undef max
#undef TRUE_  
#undef FALSE_  

#endif
