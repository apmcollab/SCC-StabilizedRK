//
// RKEigEstimator.h 
//
// Author: Chris Anderson  
// (C) UCLA 2012 
//
#include <cmath>
#include <vector>

#ifndef RK_EIG_ESTIMATOR_
#define RK_EIG_ESTIMATOR_


extern "C"  int dgelss_(long* m,  long* n, long* nrhs,double* a, long* lda, double* b, long* ldb, double* s,
double* rcond, long* rank, double* work, long* lwork, long* info);



extern "C" int dgeev_(char *jobvl, char *jobvr, long *n, double *a, long *lda, double *wr, double *wi, double *vl,
long *ldvl, double *vr, long *ldvr, double *work, long *lwork, long *info);

//
/* 
	This class creates an estimate of the dominant eigenvalues of the
	Jacobian associated with an ODE using procedure based on that
	described in the paper by Ekeland, Owren and Oines
	"Stiffness Detection and Estimation of Dominant Spectra with Explicit Runge-Kutta Methods"
    
    ACM Transactios on Mathematical Software, Vol. 24, No. 4, Dec. 1998,  Pg. 368-382
    
	The procedure implemented here differs from that described in the paper
    in that a modified Gram-Schmidt is used to carry out the orthogonalization 
    and the accumulation of the upper Hessenberg matrix is accomplished without
    explicitly constructing the inverse of the Runge-Kutta coefficient matrix. 
    
    CRA April, 14, 2009
*/


//
// In this procedure it is assumed that the RK method has the form
// 
// k[0] = f(yn)
// 
// k[i] = f(y_n + dt*(sum_(0 <= j< i) a(i,j)*k[j] ) for i = 1 to RKstageOrder-1
//
// y_(n+1) = y_(n) + dt*(sum(0 <= i <= RKstageOrder-1) b(i)*k[i] ) 
// 
//
// Here a(i,j) are the Runge-Kutta stage coefficents. The final accumulation 
// coefficients b(j) are not used.
//
// Note that the stage std::vectors k[i] are not multiplied by dt in this
// representation of the Runge-Kutta method. 
// 


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
// A double array containing the Runge-Kutta stage
// coefficients. The entries of this array have the structure
//
// a(2,1) a(3,1)  ...    a(s,1)
//    0   a(3,2)  ...      *
//    *      *             *
//    *            *       *
//    0      *     0    a(s,s-1)
//
// where a(i,j) are the coefficients of the RK method when written in the
// form described in above. 
// 
// The ith column of this matrix corresponds to the coefficients used to 
// accumulate the ith stage vector.
// 
void initialize(long RKstageOrder, std::vector< std::vector<double> >& RKcoefficients)
{
	this->RKstageOrder     = RKstageOrder; // Stage order of the RK method

    this->dependencyTol    =  1.0e-04;
    this->pseudoInverseTol =  1.0e-04;

	// Coefficients packed in format specified by Ekeland's paper
	

    std::cout << "xxxxxxxxxxxxxxxxxxxxxxxxxx" << std::endl;

    for(long j = 0; j < RKcoefficients.size(); j++)
    {
    for(long i = 0; i < RKcoefficients.size(); i++)
    {
    std::cout << RKcoefficients[i][j] << " ";
    }
    std::cout <<"\n";
    }

    std::cout << "xxxxxxxxxxxxxxxxxxxxxxxx" << std::endl;



	A1.initialize(RKcoefficients);

	std::cout << "ggggggggggggggggggggggg" << std::endl;

    for(long j = 0; j < RKcoefficients.size(); j++)
    {
    for(long i = 0; i < RKcoefficients.size(); i++)
    {
    std::cout << A1(i,j) << " ";
    }
    std::cout <<"\n";
    }

    std::cout << "ggggggggggggggggg" << std::endl;



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
    //
    // Create a QR factorization of [k_0, k_1, ....] until j = RKstageOrder or the 
    // the jth vector computed is linearly dependent on the previous vectors. 
    //

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
    
    //kNorm  = sqrt(fabs(q[j].dot(q[j])));
    for(i = 0; i < j; i++) 
    {
		r(i,j) =  q[i].dot(q[j]);
		q[j].axpy(-r(i,j),q[i]);
	}
    r(j,j) = sqrt(fabs(q[j].dot(q[j])));
    
    //
    // Estimate dependence by forming singular values of r
    //
    rCopy.initialize(j+1,j+1);
    rCopy.setToValue(0.0);

    rSingularValues.resize(j+1);
    std::fill(rSingularValues.begin(), rSingularValues.end(), 0.0);

    for(ii = 0; ii <= j; ii++)
    {
    for(jj = 0; jj <= j; jj++)
    {
    rCopy(ii,jj) = r(ii,jj);
    }}
    
    singularValues(rCopy.getDataPointer(),j+1,j+1,&rSingularValues[0]);
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
// construcion of the inverse of A1 and R because
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
/*
   
    DoubleArrayStructure2D           rTrans;
    DoubleArrayStructure2D          A1Trans;
    DoubleArrayStructure2D      KtildeTrans;
    
    rTrans.initialize(hSize,hSize);
    A1Trans.initialize(hSize,hSize);
    KtildeTrans.initialize(hSize,hSize);
    
    for(i = 0; i < hSize; i++)
    {
    for(j = 0; j < hSize; j++)
    {
     rTrans(j,i)     = r(i,j);
    A1Trans(j,i)     = A1(i,j);
    KtildeTrans(j,i) = Ktilde(i,j);
    }} 
    
    rSingularValues.initialize(hSize);
    SVDutility::singularValues(A1Trans.getDataPointer(),hSize,hSize,rSingularValues.getDataPointer());
    cout << "A1 singular values" << endl;
    for(i = 0; i < hSize; i++)
    {
    cout << rSingularValues(i)/rSingularValues(0) << " ";  
    }
    cout << endl;
    
    rSingularValues.initialize(hSize);
    SVDutility::singularValues(rTrans.getDataPointer(),hSize,hSize,rSingularValues.getDataPointer());
    cout << "r singular values" << endl;
    for(i = 0; i < hSize; i++)
    {
    cout << rSingularValues(i)/rSingularValues(0) << " ";  
    }
    cout << endl;
//
//  Evaluate [H^T]  using the pseudo inverse 
//
    SVDutility::svdSolve(A1Trans.getDataPointer(), hSize, hSize, 
    KtildeTrans.getDataPointer(), hSize,Z.getDataPointer(),pseudoInverseTol);
    
    SVDutility::svdSolve(rTrans.getDataPointer(), hSize, hSize, 
    Z.getDataPointer(), hSize,Ht.getDataPointer(),pseudoInverseTol);
    
*/
 


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


   
    /*
    cout << " Z " << endl;
    for(j = 0; j < hSize; j++)
    {
    for(i = 0; i < hSize; i++)
    {
    cout << Z(j,i) << " ";
    }
    cout << endl;
    }
    cout << endl;
    
     
    cout << " Ht " << endl;
    for(j = 0; j < hSize; j++)
    {
    for(i = 0; i < hSize; i++)
    {
    cout << Ht(j,i) << " ";
    }
    cout << endl;
    }
    cout << endl;
    exit(0);

    cout << " R " << endl;
    for(j = 0; j < hSize; j++)
    {
    for(i = 0; i < hSize; i++)
    {
    cout << r(j,i) << " ";
    }
    cout << endl;
    }
    cout << endl;
    */
    

//
//  Compute the eigenvalues of Ht = eigenvalues of H 
// 
//  For now, just call the standard dense eigenvalue solver. This
//  could be replaced by passing the transpose of Ht to an eigensystem 
//  routine that assumes the input matrix is in upper-Hessenberg form.
// 
	eigCount = hSize;

	Wreal.clear();
	Wreal.resize(hSize);

	Wimag.clear();
	Wimag.resize(hSize);

    long info;
    
    /*
    cout << " HT " << endl;
    cout << endl;
    
    for(i = 0; i < hSize; i++)
    {
    for(j = 0; j < hSize; j++)
    {
        cout <<Ht(i,j) << " ";
        }
        cout << endl;
        }
    */ 
    
	computeEigenvalues(Ht.getDataPointer(),hSize, &Wreal[0], &Wimag[0]);
	
	
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


void singularValues(double* A, long m, long n,double* singularValues)
{
    long i; long j;


	double* Q = new double[m*n];

	for(i = 0; i < m; i++)
	{
	for(j = 0; j < n; j++)
	{
	Q[i + j*m] = A[i + j*m];
	}}

	//for(size_t k = 0; k < (size_t)(m*n); ++k)
	//{
	//Q[k] = A[k];
	//}

    long maxDim = std::max(m,n);


    long nrhs    = 0;  // just want singular values
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
    std::cerr << "Singular Value Computation Failed " << std::endl;
//
//  clean up
//
	delete [] workPtr;
    delete [] Bptr;
    delete [] Q;
}

void computeEigenvalues(double* Adata, long Asize, double* wReal, double* wImag)
{
	char jobvl        = 'N';
	char jobvr        = 'N';

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
	std::vector<double> wr(n,0.0);
	std::vector<double> wi(n,0.0);

	double    vl  = 0.0;
	double    vr  = 0.0;
	long     ldvl   = 1;
    long     ldvr   = 1;

	long lwork      = 5*n;
	double* work    = new double[lwork];
	long info       = 0;

	dgeev_(&jobvl, &jobvr, &n, a, &lda, &wr[0], &wi[0], &vl, &ldvl, &vr, &ldvr, work,&lwork, &info);

    for(long i = 0; i < n; i++)
    {
    wReal[i] = wr[i];
    wImag[i] = wi[i];
    }

	delete [] a;
	delete [] work;

	if(info != 0)
	{
	 std::cerr << " Eigenvalue routine dgeev (general matrix) failed " << std::endl;
	 std::cerr << " Info Value: " << info << std::endl;
     exit(0);
	}
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

