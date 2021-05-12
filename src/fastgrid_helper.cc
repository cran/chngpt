#ifndef fastgrid_H
#define fastgrid_H


#include "fastgrid_helper.h"
#include "matrix.h"
#include "distributions.h"
#include "stat.h"
#include "la.h"
#include "ide.h"
#include "smath.h"

#include <R.h>           // needed to use Rprintf()
#include <R_ext/Utils.h> // needed to allow user interrupts

#include <Rdefines.h>
#include <Rinternals.h>

#include <float.h> //DBL_EPSILON
#include <R_ext/Lapack.h>
#include <Rmath.h>

#define RUNIF runif
#define PRINTF Rprintf
#define MAX(A,B)    ((A) > (B) ? (A) : (B))
#define MIN(A,B)    ((A) < (B) ? (A) : (B))


using namespace scythe;
using namespace std;



// modified from ide.h, but because there was lots of error when SCYTHE_LAPACK is defined, it cannot be used
// modifications include removing of lapack
Matrix<> qr_solve (const Matrix<>& A, const Matrix<>& b)
{
    SCYTHE_DEBUG_MSG("Using lapack/blas for qr_solve");
    SCYTHE_CHECK_10(A.isNull(), scythe_null_error, "A is NULL")
    SCYTHE_CHECK_10(A.rows() != b.rows(), scythe_conformation_error,
        "A and b do not conform");

    /* Do decomposition */
   
    // Set up working variables
    Matrix<> QR = A;
    double* QRarray = QR.getArray(); // input/output array pointer
    int rows = (int) QR.rows();
    int cols = (int) QR.cols();
//    Matrix<unsigned int> pivot(cols, 1); // pivot vector
//    int* parray = (int*) pivot.getArray(); // pivot vector array pointer
    Matrix<> tau = Matrix<>(rows < cols ? rows : cols, 1);
    double* tarray = tau.getArray(); // tau output array pointer
    double tmp, *work; // workspace vars
    int lwork, info;   // workspace size var and error info var

//    for (unsigned int i=0; i<pivot.size(); i++) PRINTF("%u ", pivot(i)); PRINTF("\n");       

    // the pivot code from ide.h does not work for some reason - pivot seem to get the wrong value
//    // Get workspace size
//    lwork = -1;
//    dgeqp3_(&rows, &cols, QRarray, &rows, parray, tarray, &tmp, &lwork, &info);
//    SCYTHE_CHECK_10(info != 0, scythe_lapack_internal_error, "Internal error in LAPACK routine dgeqp3");
//    lwork = (int) tmp;
//    work = new double[lwork];
//    // run the routine for real
//    dgeqp3_(&rows, &cols, QRarray, &rows, parray, tarray, work, &lwork, &info);
//    SCYTHE_CHECK_10(info != 0, scythe_lapack_internal_error, "Internal error in LAPACK routine dgeqp3");
//    delete[] work;    

    // first get workspace size
    lwork = -1;
    dgeqrf_(&rows, &cols, QRarray, &rows, tarray, &tmp, &lwork, &info);
    SCYTHE_CHECK_10(info != 0, scythe_lapack_internal_error, "Internal error in LAPACK routine dgeqp3");    
    lwork = (int) tmp;
    work = new double[lwork];    
    // now run the routine for real
    dgeqrf_(&rows, &cols, QRarray, &rows, tarray, work, &lwork, &info);
    SCYTHE_CHECK_10(info != 0, scythe_lapack_internal_error, "Internal error in LAPACK routine dgeqp3");    
    delete[] work;    

//    //pivot-=1;
//    for (int idx=0; idx<cols; idx++) pivot(idx) = pivot(idx) - 1;
//    for (unsigned int i=0; i<pivot.size(); i++) PRINTF("%u ", pivot(i)); PRINTF("\n");       

    /* Now solve the system */
    
    // working vars
    int nrhs = (int) b.cols();
    Matrix<> bb = b;
    double* barray = bb.getArray();
    int taudim = (int) tau.size();

    // Get workspace size
    lwork = -1;
    dormqr_("L", "T", &rows, &nrhs, &taudim, QRarray, &rows,
                    tarray, barray, &rows, &tmp, &lwork, &info);

    SCYTHE_CHECK_10(info != 0, scythe_lapack_internal_error,
        "Internal error in LAPACK routine dormqr");

    // And now for real
    lwork = (int) tmp;
    work = new double[lwork];
    dormqr_("L", "T", &rows, &nrhs, &taudim, QRarray, &rows,
                    tarray, barray, &rows, work, &lwork, &info);

    SCYTHE_CHECK_10(info != 0, scythe_lapack_internal_error,
        "Internal error in LAPACK routine dormqr");

    dtrtrs_("U", "N", "N", &taudim, &nrhs, QRarray, &rows, barray,
                    &rows, &info);

    SCYTHE_CHECK_10(info > 0, scythe_type_error, "Matrix is singular");
    SCYTHE_CHECK_10(info < 0, scythe_lapack_internal_error,
        "Internal error in LAPACK routine dtrtrs");

    delete[] work;

    Matrix<> result(A.cols(), b.cols(), false);
    for (int i = 0; i < cols; ++i)
      result(i, _) = bb(i, _);

//    for (unsigned int i=0; i<pivot.size(); i++) PRINTF("%u ", pivot(i)); PRINTF("\n");       
//    for (int i=0; i<rows; i++) PRINTF("%f ", bb(i)); PRINTF("\n");       


    return result;
}





// taken from R's src/main/random.c
void SampleReplace(int k, int n, int *y)
{
    int i;
#ifndef SCYTHE_COMPILE_DIRECT    
    GetRNGstate();    
#endif
    for (i = 0; i < k; i++) y[i] = n * unif_rand() + 1;
#ifndef SCYTHE_COMPILE_DIRECT    
    PutRNGstate();    
#endif
}

// 'x' is preallocated of length n
// output: 'y' of length 'k' preallocated
void SampleNoReplace(int k, int n, int *y, int *x)
{
    int i, j;
#ifndef SCYTHE_COMPILE_DIRECT    
    GetRNGstate();    
#endif
    for (i = 0; i < n; i++) x[i] = i;
    for (i = 0; i < k; i++) {
        j = (int)((double)n * RUNIF(0.0,1.0));
        y[i] = x[j] + 1;
        x[j] = x[--n];
    }
#ifndef SCYTHE_COMPILE_DIRECT    
    PutRNGstate();    
#endif
}


void make_symmetric(double* matrix, int rows)
{
      for (int i = 1; i < rows; ++i)
        for (int j = 0; j < i; ++j)
          matrix[i * rows + j] = matrix[j * rows + i];
}
  // it is not clear to me whether crossprod is calling lapack or not. crossprod1 is the way I make sure it is
  // if a row-major matrix is passed as A, it will be transposed automatically
Matrix<> crossprod1(const Matrix<>& A)
{
    SCYTHE_DEBUG_MSG("Using lapack/blas for crossprod");
    // Set up some constants
    const double zero = 0.0;
    const double one = 1.0;

    // Set up return value and arrays
    Matrix<> res(A.cols(), A.cols(), false);
    double* Apnt = A.getArray();
    double* respnt = res.getArray();
    int rows = (int) A.rows();
    int cols = (int) A.cols();
    //for (int i=0; i<rows*cols; i++) PRINTF("%f ", Apnt[i]); PRINTF("\n");       

    dsyrk_("L", "T", &cols, &rows, &one, Apnt, &rows, &zero, respnt,
                   &cols);
    make_symmetric(respnt, cols); 

    return res;
}


Matrix<> tcrossprod1(const Matrix<>& A)
{
    SCYTHE_DEBUG_MSG("Using lapack/blas for crossprod");
    // Set up some constants
    const double zero = 0.0;
    const double one = 1.0;

    // Set up return value and arrays
    Matrix<> res(A.rows(), A.rows(), false); // this is changed from crossprod1
    double* Apnt = A.getArray();
    double* respnt = res.getArray();
    int rows = (int) A.rows();
    int cols = (int) A.cols();
    //for (int i=0; i<rows*cols; i++) PRINTF("%f ", Apnt[i]); PRINTF("\n");       

    // several arguments are changed compared to crossprod1
    dsyrk_("L", "N", &rows, &cols, &one, Apnt, &rows, &zero, respnt, &rows); 
    make_symmetric(respnt, rows); 

    return res;
}




Matrix<> myqr_getQ (const Matrix<>& A)
{
    // Set up working variables
    double* QRarray = A.getArray(); // input/output array pointer
    int rows = (int) A.rows();
    int cols = (int) A.cols();
    Matrix<> tau = Matrix<>(rows < cols ? rows : cols, 1);
    double* tarray = tau.getArray(); // tau output array pointer
    double tmp, *work; // workspace vars
    int lwork, info;   // workspace size var and error info var
    
    // QR decomposition
    // dgeqrf is faster than dgeqr2
    // note that dgeqpf solves a different equation with column pivoting of the input matrix
    
    // first get workspace size
    lwork = -1;
    dgeqrf_(&rows, &cols, QRarray, &rows, tarray, &tmp, &lwork, &info);
    SCYTHE_CHECK_10(info != 0, scythe_lapack_internal_error, "Internal error in LAPACK routine dgeqp3");    
    lwork = (int) tmp;
    work = new double[lwork];    
    // now run the routine for real
    dgeqrf_(&rows, &cols, QRarray, &rows, tarray, work, &lwork, &info);
    SCYTHE_CHECK_10(info != 0, scythe_lapack_internal_error, "Internal error in LAPACK routine dgeqp3");    
    delete[] work;    
//    dgeqp3_(&rows, &cols, QRarray, &rows, parray, tarray, &tmp, &lwork, &info);
//    dgeqf3_(&rows, &cols, QRarray, &rows, parray, tarray, work, &lwork, &info);

    //pivot -= 1;
    
    // Create orthogonal matrix Q (in QRarray)
    // Get workspace size
    lwork = -1;
    dorgqr_(&rows, &cols, &cols, QRarray, &rows, tarray, &tmp, &lwork, &info);
    lwork = (int) tmp;
    work = new double[lwork];
    // run the routine for real
    dorgqr_(&rows, &cols, &cols, QRarray, &rows, tarray, work, &lwork, &info);
    delete[] work;
    //for (int i=0; i<rows; i++) {for (int j=0; j<cols; j++)  PRINTF("%f ", A(i,j)); PRINTF("\n");}            

    Matrix<> result(A);
//    for(int i =0; i < rows; ++i) {
//        memcpy(_Q+row*rank, _A+row*_n, sizeof(double)*(rank));
//    }
        
    return result;
}


// replace Z and Y with B and r
// B is Q_1, the first p columns of Q
// H = Q_1 * Q_1' (https://stats.stackexchange.com/questions/139969/speeding-up-hat-matrices-like-xxx-1x-projection-matrices-and-other-as)
void _preprocess(Matrix<double, Row>& Z, Matrix<double, Row>& Y) {
     
    // compute Q
    Matrix<> Q = myqr_getQ (Z);       
    // assign Q to Z
    int p=Z.cols(), n=Z.rows(); for (int i=0; i<n; i++) for (int j=0; j<p; j++) Z(i,j) = Q(i,j);
    // compute resid
    Y = Y - tcrossprod1(Q) * Y; 
//    PRINTF("Z r\n"); for (int i=0; i<n; i++) {for (int j=0; j<p; j++)  PRINTF("%f ", Q(i,j));   PRINTF("%f ", Y(i,0));   PRINTF("\n");}                
     
     // the old implementation based on inverting X'X
//    int p=Z.cols();
//    Matrix<> A = invpd(crossprod(Z));
//    Eigen Aeig = eigen(A);
//    Matrix <double,Row,Concrete> A_eig (p, p, true, 0);
//    for (int j=0; j<p; j++) A_eig(j,j)=sqrt(Aeig.values(j));
//    Y = Y - Z * A * (t(Z) * Y);                       // save r in Y
//    Z = Z * (Aeig.vectors * A_eig * t(Aeig.vectors)); // save B in Z    
////    PRINTF("Z r\n"); for (int i=0; i<Z.rows(); i++) {for (int j=0; j<p; j++) PRINTF("%f ", Z(i,j));   PRINTF("%f ", Y(i,0)); PRINTF("\n");}            

}




#endif /* fastgrid_H */


//struct Eigen {
//   Matrix<> values;
//   Matrix<> vectors;
//};
//
//Eigen eigen (const Matrix<>& A, bool vectors=true)
//{
//    SCYTHE_DEBUG_MSG("Using lapack/blas for eigen");
//    SCYTHE_CHECK_10(! A.isSquare(), scythe_dimension_error,
//        "Matrix not square");
//    SCYTHE_CHECK_10(A.isNull(), scythe_null_error,
//        "Matrix is NULL");
//    // Should be symmetric but rounding errors make checking for this
//    // difficult.
//
//    // Make a copy of A
//    Matrix<> AA = A;
//
//    // Get a point to the internal array and set up some vars
//    double* Aarray = AA.getArray(); // internal array points
//    int order = (int) AA.rows();    // input matrix is order x order
//    double dignored = 0;            // we do not use this option
//    int iignored = 0;               // or this one
//    double abstol = 0.0;            // tolerance (default)
//    int m;                          // output value
//    Matrix<> result;                // result matrix
//    char getvecs[1];                // are we getting eigenvectors?
//    if (vectors) {
//      getvecs[0] = 'V';
//      result = Matrix<>(order, order + 1, false);
//    } else {
//      result = Matrix<>(order, 1, false);
//      getvecs[0] = 'N';
//    }
//    double* eigenvalues = result.getArray(); // pointer to result array
//    int* isuppz = new int[2 * order];        // indices of nonzero eigvecs
//    double tmp;   // inital temporary value for getting work-space info
//    int lwork, liwork, *iwork, itmp; // stuff for workspace
//    double *work; // and more stuff for workspace
//    int info = 0;  // error code holder
//
//    // get optimal size for work arrays
//    lwork = -1;
//    liwork = -1;
//    dsyevr_(getvecs, "A", "L", &order, Aarray, &order, &dignored,
//        &dignored, &iignored, &iignored, &abstol, &m, eigenvalues, 
//        eigenvalues + order, &order, isuppz, &tmp, &lwork, &itmp,
//        &liwork, &info);
//    SCYTHE_CHECK_10(info != 0, scythe_lapack_internal_error,
//        "Internal error in LAPACK routine dsyevr");
//    lwork = (int) tmp;
//    liwork = itmp;
//    work = new double[lwork];
//    iwork = new int[liwork];
//
//    // do the actual operation
//    dsyevr_(getvecs, "A", "L", &order, Aarray, &order, &dignored,
//        &dignored, &iignored, &iignored, &abstol, &m, eigenvalues, 
//        eigenvalues + order, &order, isuppz, work, &lwork, iwork,
//        &liwork, &info);
//    SCYTHE_CHECK_10(info != 0, scythe_lapack_internal_error,
//        "Internal error in LAPACK routine dsyevr");
//
//    delete[] isuppz;
//    delete[] work;
//    delete[] iwork;
//    
//    Eigen resobj;
//    if (vectors) {
//      resobj.values = result(_, 0);
//      resobj.vectors = result(0, 1, result.rows() -1, result.cols() - 1);
//    } else {
//      resobj.values = result;
//    }
//
//    return resobj;
//}


