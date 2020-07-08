//
// GenEigUtility.cpp 
//
// Author: Chris Anderson  
// (C) UCLA 2012 
//
//
// Wrapper class for general matrix eigenvalue/eigenvector routines.
// 
// This class provides a C++ interface to routines implemented in the library
// libcamgeneig.a with dependencies upon camblas and camlapack. The required library 
// linking statement has the form 
//
// -lcamgeneig  -lcamlapack -lcamblas
//
#include <iostream>
using namespace std;

#include "GenEigUtility.h"

/*  Purpose */
/*  ======= */ 

/*  DGEEV computes for an N-by-N real nonsymmetric matrix A, the */
/*  eigenvalues and, optionally, the left and/or right eigenvectors. */

/*  The right eigenvector v(j) of A satisfies */
/*                   A * v(j) = lambda(j) * v(j) */
/*  where lambda(j) is its eigenvalue. */
/*  The left eigenvector u(j) of A satisfies */
/*                u(j)**H * A = lambda(j) * u(j)**H */
/*  where u(j)**H denotes the conjugate transpose of u(j). */

/*  The computed eigenvectors are normalized to have Euclidean norm */
/*  equal to 1 and largest component real. */

/*  Arguments */
/*  ========= */

/*  JOBVL   (input) CHARACTER*1 */
/*          = 'N': left eigenvectors of A are not computed; */
/*          = 'V': left eigenvectors of A are computed. */

/*  JOBVR   (input) CHARACTER*1 */
/*          = 'N': right eigenvectors of A are not computed; */
/*          = 'V': right eigenvectors of A are computed. */

/*  N       (input) long */
/*          The order of the matrix A. N >= 0. */

/*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*          On entry, the N-by-N matrix A. */
/*          On exit, A has been overwritten. */

/*  LDA     (input) long */
/*          The leading dimension of the array A.  LDA >= max(1,N). */

/*  WR      (output) DOUBLE PRECISION array, dimension (N) */
/*  WI      (output) DOUBLE PRECISION array, dimension (N) */
/*          WR and WI contain the real and imaginary parts, */
/*          respectively, of the computed eigenvalues.  Complex */
/*          conjugate pairs of eigenvalues appear consecutively */
/*          with the eigenvalue having the positive imaginary part */
/*          first. */

/*  VL      (output) DOUBLE PRECISION array, dimension (LDVL,N) */
/*          If JOBVL = 'V', the left eigenvectors u(j) are stored one */
/*          after another in the columns of VL, in the same order */
/*          as their eigenvalues. */
/*          If JOBVL = 'N', VL is not referenced. */
/*          If the j-th eigenvalue is real, then u(j) = VL(:,j), */
/*          the j-th column of VL. */
/*          If the j-th and (j+1)-st eigenvalues form a complex */
/*          conjugate pair, then u(j) = VL(:,j) + i*VL(:,j+1) and */
/*          u(j+1) = VL(:,j) - i*VL(:,j+1). */

/*  LDVL    (input) long */
/*          The leading dimension of the array VL.  LDVL >= 1; if */
/*          JOBVL = 'V', LDVL >= N. */

/*  VR      (output) DOUBLE PRECISION array, dimension (LDVR,N) */
/*          If JOBVR = 'V', the right eigenvectors v(j) are stored one */
/*          after another in the columns of VR, in the same order */
/*          as their eigenvalues. */
/*          If JOBVR = 'N', VR is not referenced. */
/*          If the j-th eigenvalue is real, then v(j) = VR(:,j), */
/*          the j-th column of VR. */
/*          If the j-th and (j+1)-st eigenvalues form a complex */
/*          conjugate pair, then v(j) = VR(:,j) + i*VR(:,j+1) and */
/*          v(j+1) = VR(:,j) - i*VR(:,j+1). */

/*  LDVR    (input) long */
/*          The leading dimension of the array VR.  LDVR >= 1; if */
/*          JOBVR = 'V', LDVR >= N. */

/*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK)) */
/*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */

/*  LWORK   (input) long */
/*          The dimension of the array WORK.  LWORK >= max(1,3*N), and */
/*          if JOBVL = 'V' or JOBVR = 'V', LWORK >= 4*N.  For good */
/*          performance, LWORK must generally be larger. */

/*          If LWORK = -1, then a workspace query is assumed; the routine */
/*          only calculates the optimal size of the WORK array, returns */
/*          this value as the first entry of the WORK array, and no error */
/*          message related to LWORK is issued by XERBLA. */

/*  INFO    (output) long */
/*          = 0:  successful exit */
/*          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/*          > 0:  if INFO = i, the QR algorithm failed to compute all the */
/*                eigenvalues, and no eigenvectors have been computed; */
/*                elements i+1:N of WR and WI contain eigenvalues which */
/*                have converged. */


/* Subroutine */ extern "C" int dgeev_(char *jobvl, char *jobvr, long *n, double *
	a, long *lda, double *wr, double *wi, double *vl, 
	long *ldvl, double *vr, long *ldvr, double *work, 
	long *lwork, long *info, int jobvl_len, int jobvr_len);

//
// The routine computeEigenvalues(...) computes the eigenvalues of a general 
// double matrix. 
//
// --Input-- 
//
// Adata  : array of doubles of size Asize*Asize specifying the matrix entries. Since
//         only the eigenvalues are computed, there is no need to transpose the data
//         to conform to a Fortran storage convention.
//
// Asize  : the dimension of the matrix
//
// Wreal  : array of doubles of size Asize 
// Wimag  : array of doubles of size Asize 
//
// Output  : Wreal and Wimag contain the real and imaginary parts of the computed eigenvalues
//         
// Returns : info, an long indicating return status (0 = ok) (!=0 see documentataion)
//
//
// ---> No bounds checking is performed on the input data <---
//
// 
// CRA April 10, 2009
//
int GenEigUtility::computeEigenvalues(double* Adata, long Asize, double* wReal, double* wImag)
{
	char* jobvl        = (char*)"N";
	char* jobvr        = (char*)"N";
	
	long            n  = Asize;
	double*         a  = new double[n*n];
	
	//
	//  The matrix data is overwritten by the eigenvalue routine
	//  so one must make a copy of it calling the eigensystem 
	//  routine. No need to worry about transposing because we aren't 
	//  computing the eigenvectors.
	//
	
	for(long i = 0; i < n*n; i++) {a[i] = Adata[i];}
	
	long lda        = n;
	double*    wr  = wReal;
	double*    wi  = wImag;
	
	double*    vl  = 0;
	double*    vr  = 0;
	long     ldvl   = 1;
  long     ldvr   = 1;
	
	long lwork      = 5*n;
	double* work    = new double[lwork];
	long info       = 0;
	int jobvl_len   = 1;
	int jobvr_len   = 1;
	
	dgeev_(jobvl, jobvr, &n,a, &lda, wr, wi, vl, &ldvl, vr, &ldvr, work, 
	&lwork, &info, jobvl_len, jobvr_len);
	
	delete [] a;
	delete [] work;
	
	if(info != 0)
	{
	 cerr << " Eigenvalue routine dgeev (general matrix) failed " << endl;
	 cerr << " Info Value: " << info << endl;
	 return info;
	}
	return info;
}


