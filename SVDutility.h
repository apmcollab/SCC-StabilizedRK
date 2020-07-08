//
// SVDutility.h 
//
// Author: Chris Anderson  
// (C) UCLA 2012 
//

//
// Wrapper class for linpack SVD routines
// 
// This class provides a C++ interface to routines implemented in the library
// libcamgeneig.a with dependencies upon camblas and camlapack. The required library 
// linking statement has the form 
//
// -lcameig  -lcamlapack -lcamblas
//  
//  The input/output data types are standard C array structures. Input arrays are
//  passed by pointer and it is assumed that the arrays are allocated in contiguous
//  memory and are stored by rows (C convention). 
//  
//  Routines with CAMdoubleMatrix and CAMdoubleVector data types are provided by the
//  CAMqrSvdUtility class. 
//
#ifndef _SVDutility_
#define _SVDutility_

class SVDutility
{
public :
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
//  A is an m by n matrix whose data is stored by rows in contiguous memory
//  locations. 
//
//  B is an array of size m
//
//  X is an array of size n
//
//  singularValues is an array of size > min(m,n)
// 
//  The min(m,n) singular values are returned in singularValues. 
//
    static void svdSolve(double* A, long m, long n, double* B, 
    double tol, long& rank, double* singularValues, double* X);
//
//  Computes the minimum norm solution to a real linear least
//  squares problem:
//
//  Minimize 2-norm(| B - A*X |) 
//
//  using the singular value decomposition (SVD) of A. A is an m-by-n
//  matrix which may be rank-deficient. The input array B is over-written
//  with the solution X. 
//
//  tol is used to determine the effective rank of A.
//  Singular values S(i) <= RCOND*tol are treated as zero.
//  If tol $<$ 0, machine precision is used instead.
//
//  A is an m by n matrix whose data is stored by rows in contiguous memory
//  locations. 
//
//  B input/output is an array of size > max(m,n)
//
//
   static void svdSolve(double* A, long m, long n, double* B, double tol);
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
//  A is an m by n matrix whose data is stored by rows in contiguous memory
//  locations. 
//
//  B is an array of size m x nrhs
//
//  X is an array of size n x nrhs
//
   static void svdSolve(double* A, long m, long n, double* B, long nrhs, 
   double* X, double tol);
//
// Computes the singular values of A and returns them in singularValues
//
// The singularValues array must be pre-allocated and have a 
// size > min(m,n)
// 
// The min(m,n) singular values are returned in singularValues. 
//
    static void singularValues(double* A, long m, long n,double* singularValues);
};
#endif  
