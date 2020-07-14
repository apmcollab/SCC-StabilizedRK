//
// RKEigEstimator.h 
//
// Author: Chris Anderson  
// (C) UCLA 2012 
//
#include <cmath>
#include <vector>


#ifndef SCC_LAPACK_HEADERS_

extern "C" int dgeev_(char* JOBVL, char* JOBVR, long* N,double* A, long* LDA,double* WR, double* WI, double* VL,
                     long* LDVL, double* VR, long* LDVR, double* WORK, long* LWORK, long* INFO);

extern "C" void dgesvd_(char* JOBU,char* JOBVT, long* M, long* N, double* APtr, long* LDA, double* SPtr, double* UPtr, long* LDU, double* VTPtr, long* LDVT,
                       double* WORKtmp, long* LWORK, long* INFO);
#endif

#ifndef RK_EIG_ESTIMATOR_
#define RK_EIG_ESTIMATOR_

//
/* 
	This class creates an estimate of the dominant eigenvalues of the
	Jacobian associated with an ODE using procedure based upon that
	described in the paper by Ekeland, Owren and Oines

	"Stiffness Detection and Estimation of Dominant Spectra with Explicit Runge-Kutta Methods"
    ACM Transactios on Mathematical Software, Vol. 24, No. 4, Dec. 1998,  Pg. 368-382
    
	The procedure implemented here differs from that described in the paper
    in that a modified Gram-Schmidt is used to carry out the orthogonalization 
    and the accumulation of the upper Hessenberg matrix is accomplished without
    explicitly constructing the inverse of the Runge-Kutta coefficient matrix. 

    
    CRA April 14, 2009 : Original creation
    CRA July  12, 2020 : Refactoring to remove external dependencies
*/



template <class Vector> class RKEigEstimator  
{
public:
//###############################################################
//                      initialize(...) 
//###############################################################
//
// Input: 
// 
// double RKstageOrder : 
//
// The number of stages in the Runge-Kutta method
//
// DoubleArrayStructure2D RKcoefficients :  
//
// A double array of size (RKstageOrder-1) X (RKstageOrder-1) containing a
// the Runge-Kutta coefficients used to construct the stages.
//
// If the s = RKstageOrder is the stage order of the method then
// the coefficients of the RK method are assumed to be in an s-1 X s-1
// upper triangular array containing the coefficients used to
// determine the stages K[1] through K[s-1]. (K[0] is always F(y_n))
//
//
//    | alpha[0][0] | alpha[0][1] | alpha[0][2] | .....   | alpha[0][s-2]   |
//    |     0       | alpha[1][1] | alpha[1][2] | .....   | alpha[1][s-2]   |
//    |     0       |      0      | alpha[2][2] | .....   | alpha[2][s-2]   |
//    |     0       |      0      |       0     |         |                 |
//    |     0       |      0      |       0     |    0    | alpha[s-2][s-2] |
//----------------------------------------------------------------------------
//  |      |              |                                     |
//
// K[0]   K[1]           K[2]                                 K[s-1]
//
// These coefficients correspond to a loop for constructing the RK stages of
// the form
//
// K[0] = F(y_n)
//
// For j = 1 to s-1
//
// K[j] = F(y_n + dt*( sum_(0 <= i <= j-1)  alpha[i][j-1]*K[i] )
//
// In the estimation procedure the coefficients of the RK method
// used to construct y_(n+1) from the stages are not used.
//
// 
void initialize(long RKstageOrder, std::vector< std::vector<double> >& RKcoefficients)
{
	this->RKstageOrder     = RKstageOrder; // Stage order of the RK method

    this->dependencyTol    =  1.0e-04;
    this->pseudoInverseTol =  1.0e-04;


	A1.initialize(RKcoefficients);
	r.initialize(RKstageOrder,RKstageOrder);
}

//
// The input stages are defined without a dt multiplicative factor, e.g.
// k[i] = f(y_n + dt*(sum_(j< i) a_(i,j)*k[j] )

void estimateEigenvalues(std::vector < Vector* > k,
std::vector<double>& Wreal, std::vector<double>& Wimag,
long& eigCount, double dt)
{
	int stageScalingFlag = 0;
	double stageScaling  = 0.0;

	estimateEigenvalues(k,stageScalingFlag, stageScaling, Wreal, Wimag, eigCount,dt);
}

void estimateEigenvalues(std::vector < Vector* > k, double stageScaling,
std::vector<double>& Wreal, std::vector<double>& Wimag, long& eigCount, double dt)
{
	int stageScalingFlag = 1;
	estimateEigenvalues(k,stageScalingFlag, stageScaling, Wreal, Wimag, eigCount,dt);
}

void estimateEigenvalues(std::vector < Vector* > k, int stageScalingFlag, double stageScaling,
std::vector<double>& Wreal, std::vector<double>& Wimag, long& eigCount, double dt)
{
    long   qSize = 0;
    
    long i; long j; long ii; long jj; long p;
 
    r.setToValue(0.0);

    // Create a QR factorization of [k_0, k_1, ....] until j = RKstageOrder or the 
    // the jth vector computed is linearly dependent on the previous vectors. 

    p = RKstageOrder;
    
    for(j = 0; j < p; j++) 
    {
    	
    if(j >= (long)q.size())
    {
    	if(j == 0)
    	{
    	q.clear();
    	q.resize(1,*k[0]);
    	}
    	else
    	{
    	q.push_back(*k[0]);
    	}
    }
    
    q[j] = *k[j];
    
    if(stageScalingFlag == 1) 
    {
    q[j] *= stageScaling;
    }
    
    for(i = 0; i < j; i++) 
    {
		r(i,j) =  q[i].dot(q[j]);
		q[j].axpy(-r(i,j),q[i]);
	}
    r(j,j) = std::sqrt(std::abs(q[j].dot(q[j])));
    
    //
    // Estimate dependence by forming singular values of r
    //
    rCopy.initialize(j+1,j+1);
    rCopy.setToValue(0.0);

    for(ii = 0; ii <= j; ii++)
    {
    for(jj = 0; jj <= j; jj++)
    {
    rCopy(ii,jj) = r(ii,jj);
    }}
    
    singularValues(rCopy.getDataPointer(),j+1,j+1,rSingularValues);

    if(rSingularValues[j]/rSingularValues[0] <  dependencyTol)
    {
    qSize = j;
    break;
    }
    q[j]  *= 1.0/r(j,j);
    }
    
    if(qSize == 0) {qSize = p;}
    //
    // qSize is the size of the orthonormal basis 
    // computed
    // 

    long hSize  = qSize;
    if(hSize == RKstageOrder) 
    {
    	hSize -= 1;
    }
    
    /*
    for(ii = 0; ii < j; ii++)
    {
    cout << rSingularValues(ii) << " ";
    }
    cout << endl;
    */
    
//  Ktilde = [Q^T]*([k_1, k_3, .... k_hSize-1] - [k_0, k_0, .... k_0])
//
//  The inner products in this matrix need not be computed, as the
//  required inner products are available from the R component of the
//  QR factorization of [k_0, k_1, ....] obtained above. 
// 
	Ktilde.initialize(hSize,hSize);
	Ktilde.setToValue(0.0);
	
	for(j = 1; j < hSize; j++)
	{
	for(i = 0; i <= j;    i++)
	{
	Ktilde(i,j-1) = r(i,j);
	}}
	
	for(i = 0; i < hSize; i++)
	{
	Ktilde(i,hSize-1) = r(i,hSize);
	}
	
	for(j = 0; j < hSize; j++)
	{
	Ktilde(0,j) -= r(0,0);
	}

//
// Construct the transpose of 
//
// [H]  = [Ktilde] * [A1]^-1 * [R]^-1
// 
    Z.initialize(hSize,hSize);
    Z.setToValue(0.0);
    Ht.initialize(hSize,hSize);
    Ht.setToValue(0.0);
//
// The method of construction can be accomplished without the explicit 
// construction of the inverse of A1 and R because
// 
// [H]  = [Ktilde] * [A1]^-1 * [R]^-1
//
// => 
//
// [H^T] = [R^T]^-1 * [A1^T]^-1 * Ktilde^T 
// 
// A1 and R are upper triangular so the application of the 
// inverse on the right can be accomplished by back-solving 
// with lower triangular matrices. 

    for(j = 0; j < hSize; j++)
    {
    for(i = 0; i < hSize; i++)
    {
    Z(i,j) = Ktilde(j,i);
    for(p = 0; p < i; p++)
    {
    Z(i,j) -= A1(p,i)*Z(p,j);
    }
    Z(i,j) *= 1.0/A1(i,i);
    }}
    

    for(j = 0; j < hSize; j++)
    {
    for(i = 0; i < hSize; i++)
    {
    Ht(i,j) = Z(i,j);
    for(p = 0; p < i; p++)
    {
    Ht(i,j) -= r(p,i)*Ht(p,j);
    }
    Ht(i,j) *= 1.0/r(i,i);
    }}

//
//  Compute the eigenvalues of Ht = eigenvalues of H 
// 
//  For now, just call the standard dense eigenvalue solver. This
//  could be replaced by passing the transpose of Ht to an eigensystem 
//  routine that assumes the input matrix is in upper-Hessenberg form.
// 
	eigCount = hSize;

	computeEigenvalues(Ht.getDataPointer(),hSize, Wreal, Wimag);
	
	for(i = 0; i < eigCount; i++)
	{
	Wreal[i]*= 1.0/dt;
    Wimag[i]*= 1.0/dt;
	}
}

    class DoubleArrayStructure2D
    {
    public:

    DoubleArrayStructure2D() {rows = 0; cols = 0; data.clear();}

    DoubleArrayStructure2D(long rows, long cols)
    {
    initialize(rows,cols);
    }

    void initialize()
    {
    rows = 0;
    cols = 0;
    data.clear();
    }

    void initialize(long rowSize, long colSize)
    {
    if(rowSize*colSize == 0){initialize(); return;}

    rows = rowSize;
    cols = colSize;
    data.resize(rows*cols);
    std::fill(data.begin(),data.end(),0.0);
    }

    void initialize(const std::vector< std::vector<double> >& A)
    {
    if(A.size() == 0) {initialize(); return;}
    rows = A.size();
    cols = A[0].size();

    initialize(rows,cols);

    for(long i = 0; i < rows; i++)
    {
    	for(long j = 0; j < cols; j++)
    	{
    		operator()(i,j) = A[i][j];
    	}
    }
    }

    double& operator()(long i, long j)
    {

    if((i > rows-1)||(j > cols -1))
    {
    std::cout << "Access error " << " i " << i << " j " << j << std::endl;
    }
    return data[i + j*rows];
    }

    const double& operator()(long i, long j) const
    {

    if((i > rows-1)||(j > cols -1))
    {
    std::cout << "Access error " << " i " << i << " j " << j << std::endl;
    }

    return data[i + j*rows];
    }

    void setToValue(double val)
    {
    std::fill(data.begin(),data.end(),val);
    }

    double* getDataPointer() {return &data[0];}

    long rows;
    long cols;
    std::vector<double> data;
    };


void singularValues(double* A, long m, long n,std::vector<double>& singularValues)
{
    // Duplicate input, since svd routine destroys matrix

	double* Q = new double[m*n];

    for(size_t k = 0; k < (size_t)(m*n); k++)
    {
        Q[k] = A[k];
    }

    char JOBU  = 'N';
    char JOBVT = 'N';
    long M     = m;
    long N     = n;
    long LDA   = m;

    singularValues.clear();
    singularValues.resize(std::min(M,N),0.0);

    double U  = 1.0;
    long LDU  = 1;

    double VT = 1.0;
	long LDVT = 1;

    long LWORK = 10*std::min(M,N) + 1;

    double* WORKtmp = new double[LWORK];

    long INFO = 0;

    dgesvd_(&JOBU, &JOBVT, &M, &N, Q, &LDA, &singularValues[0], &U, &LDU, &VT, &LDVT,
                       WORKtmp, &LWORK, &INFO);
    if(INFO != 0)
    {
    std::cerr << "Singular Value Computation Failed " << std::endl;
    exit(0);
    }

//  clean up

	delete [] WORKtmp;
    delete [] Q;
}

void computeEigenvalues(double* Adata, long Asize, std::vector<double>& wReal, std::vector<double>& wImag)
{
	char JOBVL        = 'N';
	char JOBVR        = 'N';

	long            N  = Asize;
	double* Acopy     = new double[N*N];
	for(long i = 0; i < N*N; i++) {Acopy[i] = Adata[i];}

	long LDA        = N;

	wReal.clear();
	wReal.resize(N,0.0);

	wImag.clear();
	wImag.resize(N,0.0);


	double   VL   = 0.0;
	double   VR   = 0.0;
	long     LDVL = 1;
    long     LDVR = 1;

	long     LWORK  = 5*N;
	double*  WORK   = new double[LWORK];
	long     INFO   = 0;

	dgeev_(&JOBVL, &JOBVR, &N, Acopy, &LDA, &wReal[0], &wImag[0], &VL, &LDVL, &VR, &LDVR, WORK,&LWORK, &INFO);


	if(INFO != 0)
	{
	 std::cerr << " Eigenvalue routine dgeev (general matrix) failed " << std::endl;
	 std::cerr << " Info Value: " << INFO << std::endl;
     exit(0);
	}

	delete [] Acopy;
	delete [] WORK;
}


//
//  Class data members
//
	long          RKstageOrder;
	double       dependencyTol;
    double    pseudoInverseTol;
	
	std::vector< Vector > q;

    DoubleArrayStructure2D Ht;
    DoubleArrayStructure2D r;
    DoubleArrayStructure2D A1;
    DoubleArrayStructure2D Ktilde;
    DoubleArrayStructure2D Z;
    
    DoubleArrayStructure2D rCopy;

    std::vector<double>   rSingularValues;

	long eigCount;
};

#endif

