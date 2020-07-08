//
// GenEigUtility.h 
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
// -lcameig  -lcamlapack -lcamblas
//

#include <iostream>
using namespace std;

#ifndef __GenEigUtility__
#define __GenEigUtility__


class GenEigUtility
{
public :
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
// Returns : info, an integer indicating return status (0 = ok) (!=0 see documentataion)
//
//
// ---> No bounds checking is performed on the input data <---
//
// 
// CRA April 10, 2009
//
static int computeEigenvalues(double* Adata, long Asize, double* wReal, double* wImag);
};

#endif


