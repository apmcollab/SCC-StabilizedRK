//
// SVDutility.cpp 
//
// Author: Chris Anderson  
// (C) UCLA 2012 
//


#include <iostream>
using namespace std;
#include "SVDutility.h"

/*  DGELSS computes the minimum norm solution to a real linear least */
/*  squares problem: */

/*  Minimize 2-norm(| b - A*x |). */

/*  using the singular value decomposition (SVD) of A. A is an M-by-N */
/*  matrix which may be rank-deficient. */

/*  Several right hand side vectors b and solution vectors x can be */
/*  handled in a single call; they are stored as the columns of the */
/*  M-by-NRHS right hand side matrix B and the N-by-NRHS solution matrix 
*/
/*  X. */

/*  The effective rank of A is determined by treating as zero those */
/*  singular values which are less than RCOND times the largest singular 
*/
/*  value. */

/*  Arguments */
/*  ========= */

/*  M       (input) INTEGER */
/*          The number of rows of the matrix A. M >= 0. */

/*  N       (input) INTEGER */
/*          The number of columns of the matrix A. N >= 0. */

/*  NRHS    (input) INTEGER */
/*          The number of right hand sides, i.e., the number of columns */
/*          of the matrices B and X. NRHS >= 0. */

/*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*          On entry, the M-by-N matrix A. */
/*          On exit, the first min(m,n) rows of A are overwritten with */
/*          its right singular vectors, stored rowwise. */

/*  LDA     (input) INTEGER */
/*          The leading dimension of the array A.  LDA >= max(1,M). */

/*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS) */
/*          On entry, the M-by-NRHS right hand side matrix B. */
/*          On exit, B is overwritten by the N-by-NRHS solution */
/*          matrix X.  If m >= n and RANK = n, the residual */
/*          sum-of-squares for the solution in the i-th column is given */
/*          by the sum of squares of elements n+1:m in that column. */

/*  LDB     (input) INTEGER */
/*          The leading dimension of the array B. LDB >= max(1,MAX(M,N)). 
*/

/*  S       (output) DOUBLE PRECISION array, dimension (min(M,N)) */
/*          The singular values of A in decreasing order. */
/*          The condition number of A in the 2-norm = S(1)/S(min(m,n)). */

/*  RCOND   (input) DOUBLE PRECISION */
/*          RCOND is used to determine the effective rank of A. */
/*          Singular values S(i) <= RCOND*S(1) are treated as zero. */
/*          If RCOND $<$ 0, machine precision is used instead. */

/*  RANK    (output) INTEGER */
/*          The effective rank of A, i.e., the number of singular values 
*/
/*          which are greater than RCOND*S(1). */

/*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK) 
*/
/*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */

/*  LWORK   (input) INTEGER */
/*          The dimension of the array WORK. LWORK >= 1, and also: */
/*          LWORK >= 3*N+MAX(2*N,NRHS,M) if M >= N, */
/*          LWORK >= 3*M+MAX(2*M,NRHS,N) if M < N. */
/*          For good performance, LWORK should generally be larger. */

/*  INFO    (output) INTEGER */
/*          = 0:  successful exit */
/*          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/*          > 0:  the algorithm for computing the SVD failed to converge; 
*/
/*                if INFO = i, i off-diagonal elements of an intermediate 
*/
/*                bidiagonal form did not converge to zero. */

/*  ===================================================================== 
*/

extern "C"  int dgelss_(long* m,  long* n, long* nrhs,
double* a, long* lda, double* b, long* ldb, double* s, double* rcond, long* rank,
double* work, long* lwork, long* info);
//
//  Computes the minimum norm solution to a real linear least
//  squares problem:
//
//  Minimize 2-norm(| B - A*X |).
//
//  using the singular value decomposition (SVD) of A. A is an m-by-n
//  matrix which may be rank-deficient.
//
//  tol is used to determine the effective rank of A.
//  Singular values S(i) <= RCOND*tol are treated as zero.
//  If tol $<$ 0, machine precision is used instead.
//
//  The rank and singular values are  returned. 
//
// 
// A is an m by n matrix whose data is stored by rows in contiguous memory
// locations. 
//
// B is an array of size m
//
// X is an array of size n
//
// singularValues is an array of size > min(m,n)
// 
// The min(m,n) singular values are returned in singularValues. 
//
void SVDutility::svdSolve(double* A, long m, long n, double* B, double tol, 
long& rank, double* singularValues, double* X)
{
    long i; long j;
	//
	// Make a copy of Q and transpose the entries of A at the same time
	//
	double* Q = new double[m*n];
	
	for(i = 0; i < m; i++)
	{
	for(j = 0; j < n; j++)
	{
	Q[i + j*m] = A[j + i*n];
	}}


    long minDim           = m;
    if(n < minDim) minDim = n;
    long maxDim           = m;
    if(n > maxDim) maxDim = n;

    long nrhs    = 1;      // just want singular values
    long lda     = m;

    double* Bptr = new double[maxDim];
    long ldb     = maxDim;
   
    // create duplicate of B
    
    for(i = 0; i < m; i++)
    {
    Bptr[i] = B[i];
    }

    double* SPtr =  singularValues;
    double rCond = tol;
    rank         = 0;

    long lwork   = 5*maxDim;
    double* workPtr = new double[lwork];

    long info    = 0;
    dgelss_(&m,  &n, &nrhs,Q, &lda,Bptr, &ldb, SPtr, &rCond, &rank,
    workPtr, &lwork, &info);

    if(info != 0)
    cerr << "Singular Value Computation Failed : DGELSS Info = " << info << endl;
    
    // capture the solution 
    
    for(i = 0; i < n; i++)
    {
    X[i] = B[i];
    }
    //
    //  clean up
    //
	delete [] workPtr;
    delete [] Q;
    delete [] Bptr;
}
 
//
//  Computes the minimum norm solution to a real linear least
//  squares problem:
//
//  Minimize 2-norm(| B - A*X |).
//
//  using the singular value decomposition (SVD) of A. A is an m-by-n
//  matrix which may be rank-deficient. B is overwritten with the solution X.
//
//  tol is used to determine the effective rank of A.
//  Singular values S(i) <= RCOND*tol are treated as zero.
//  If tol $<$ 0, machine precision is used instead.
//
//  The rank and singular values are  returned. 
//
// 
// A is an m by n matrix whose data is stored by rows in contiguous memory
// locations. 
//
// B (input/output) is an array of size >= max(m,n)
//
//
void SVDutility::svdSolve(double* A, long m, long n, double* B, double tol)
{
    long i; long j;
	//
	// Make a copy of Q and transpose the entries of A at the same time
	//
	double* Q = new double[m*n];
	
	for(i = 0; i < m; i++)
	{
	for(j = 0; j < n; j++)
	{
	Q[i + j*m] = A[j + i*n];
	}}

    long minDim           = m;
    if(n < minDim) minDim = n;
    long maxDim           = m;
    if(n > maxDim) maxDim = n;

    long nrhs    = 1;      // just want singular values
    long lda     = m;

    double* Bptr = B;
    long ldb     = maxDim;
   
    double* SPtr =  new double[minDim];
    double rCond = tol;
    long rank     = 0;

    long lwork   = 5*maxDim;
    double* workPtr = new double[lwork];

    long info    = 0;
    dgelss_(&m,  &n, &nrhs,Q, &lda,Bptr, &ldb, SPtr, &rCond, &rank,
    workPtr, &lwork, &info);

    if(info != 0)
    cerr << "Singular Value Computation Failed : DGELSS Info = " << info << endl;
    //
    //  clean up
    //
	delete [] workPtr;
	delete [] SPtr;
    delete [] Q;
}

//
//  Computes the minimum norm solution to a real linear least
//  squares problem:
//
//  Minimize 2-norm(| B - A*X |).
//
//  using the singular value decomposition (SVD) of A. A is an m-by-n
//  matrix which may be rank-deficient. 
//
//  tol is used to determine the effective rank of A.
//  Singular values S(i) <= RCOND*tol are treated as zero.
//  If tol $<$ 0, machine precision is used instead.
//
//  The rank and singular values are  returned. 
//
// 
// A is an m by n matrix whose data is stored by rows in contiguous memory
// locations. 
//
// B an array of size m x nrhs
//
// X an array of size n x nrhs 
//
//
void SVDutility::svdSolve(double* A, long m, long n, double* B, long nrhs, 
double* X, double tol)
{
    long i; long j;
	//
	// Make a copy of Q and transpose the entries of A at the same time
	//
	double* Q = new double[m*n];
	
	for(i = 0; i < m; i++)
	{
	for(j = 0; j < n; j++)
	{
	Q[i + j*m] = A[j + i*n];
	}}
	

    long minDim           = m;
    if(n < minDim) minDim = n;
    long maxDim           = m;
    if(n > maxDim) maxDim = n;

    long lda     = m;

    double* Bptr = new double[maxDim*nrhs];
    long ldb     = maxDim;
    
    //
	// Make a copy of B and transpose entries at the same time
	//
	for(i = 0; i < m; i++)
	{
	for(j = 0; j < nrhs; j++)
	{
	Bptr[i + j*m] = B[j + i*nrhs];
	}}
	
   
    double* SPtr =  new double[minDim];
    double rCond = tol;
    long rank     = 0;

    long lwork   = 5*maxDim;
    double* workPtr = new double[lwork];

    long info    = 0;
    dgelss_(&m,  &n, &nrhs,Q, &lda,Bptr, &ldb, SPtr, &rCond, &rank,
    workPtr, &lwork, &info);

    if(info != 0)
    cerr << "Singular Value Computation Failed : DGELSS Info = " << info << endl;
    //
    // Capture output
    //
	for(i = 0; i < n; i++)
	{
	for(j = 0; j < nrhs; j++)
	{
	X[j + i*nrhs] = Bptr[i + j*n];
	}}
	
    //
    //  clean up
    //
	delete [] workPtr;
	delete [] SPtr;
	delete [] Bptr;
    delete [] Q;
}

//
// Computes the singular values of A and returns them in singularValues
//
// The singularValues array must be pre-allocated and have a 
// size > min(m,n)
// 
// The min(m,n) singular values are returned in singularValues. 
//
void SVDutility::singularValues(double* A, long m, long n,double* singularValues)
{
    long i; long j;
	//
	// Make a copy of Q and transpose the entries of A at the same time
	//
	double* Q = new double[m*n];
	
	for(i = 0; i < m; i++)
	{
	for(j = 0; j < n; j++)
	{
	Q[i + j*m] = A[j + i*n];
	}}

    long minDim           = m;
    if(n < minDim) minDim = n;
    long maxDim           = m;
    if(n > maxDim) maxDim = n;

    long nrhs    = 0;      // just want singular values
    long lda     = m;

    double* Bptr = new double[maxDim];
    long ldb     = maxDim;

    double* SPtr = singularValues;

    double rCond = -1.0;
    long   rank  = 0;

	long lwork   = 5*maxDim;
    double* workPtr = new double[lwork];

	long info    = 0;
    dgelss_(&m,  &n, &nrhs,Q, &lda,Bptr, &ldb, SPtr, &rCond, &rank,
    workPtr, &lwork, &info);

    if(info != 0)
    cerr << "Singular Value Computation Failed " << endl;
//
//  clean up
//
	delete [] workPtr;
    delete [] Bptr;
    delete [] Q;
}

