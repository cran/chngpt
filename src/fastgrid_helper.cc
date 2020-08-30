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



struct Eigen {
   Matrix<> values;
   Matrix<> vectors;
};

Eigen eigen (const Matrix<>& A, bool vectors=true)
{
    SCYTHE_DEBUG_MSG("Using lapack/blas for eigen");
    SCYTHE_CHECK_10(! A.isSquare(), scythe_dimension_error,
        "Matrix not square");
    SCYTHE_CHECK_10(A.isNull(), scythe_null_error,
        "Matrix is NULL");
    // Should be symmetric but rounding errors make checking for this
    // difficult.

    // Make a copy of A
    Matrix<> AA = A;

    // Get a point to the internal array and set up some vars
    double* Aarray = AA.getArray(); // internal array points
    int order = (int) AA.rows();    // input matrix is order x order
    double dignored = 0;            // we do not use this option
    int iignored = 0;               // or this one
    double abstol = 0.0;            // tolerance (default)
    int m;                          // output value
    Matrix<> result;                // result matrix
    char getvecs[1];                // are we getting eigenvectors?
    if (vectors) {
      getvecs[0] = 'V';
      result = Matrix<>(order, order + 1, false);
    } else {
      result = Matrix<>(order, 1, false);
      getvecs[0] = 'N';
    }
    double* eigenvalues = result.getArray(); // pointer to result array
    int* isuppz = new int[2 * order];        // indices of nonzero eigvecs
    double tmp;   // inital temporary value for getting work-space info
    int lwork, liwork, *iwork, itmp; // stuff for workspace
    double *work; // and more stuff for workspace
    int info = 0;  // error code holder

    // get optimal size for work arrays
    lwork = -1;
    liwork = -1;
    dsyevr_(getvecs, "A", "L", &order, Aarray, &order, &dignored,
        &dignored, &iignored, &iignored, &abstol, &m, eigenvalues, 
        eigenvalues + order, &order, isuppz, &tmp, &lwork, &itmp,
        &liwork, &info);
    SCYTHE_CHECK_10(info != 0, scythe_lapack_internal_error,
        "Internal error in LAPACK routine dsyevr");
    lwork = (int) tmp;
    liwork = itmp;
    work = new double[lwork];
    iwork = new int[liwork];

    // do the actual operation
    dsyevr_(getvecs, "A", "L", &order, Aarray, &order, &dignored,
        &dignored, &iignored, &iignored, &abstol, &m, eigenvalues, 
        eigenvalues + order, &order, isuppz, work, &lwork, iwork,
        &liwork, &info);
    SCYTHE_CHECK_10(info != 0, scythe_lapack_internal_error,
        "Internal error in LAPACK routine dsyevr");

    delete[] isuppz;
    delete[] work;
    delete[] iwork;
    
    Eigen resobj;
    if (vectors) {
      resobj.values = result(_, 0);
      resobj.vectors = result(0, 1, result.rows() -1, result.cols() - 1);
    } else {
      resobj.values = result;
    }

    return resobj;
}

// get B and r and save in Z and Y
void _preprocess(Matrix<double, Row>& Z, Matrix<double, Row>& Y) {
    int p=Z.cols();
    Matrix<> A = invpd(crossprod(Z));
    Eigen Aeig = eigen(A);
    Matrix <double,Row,Concrete> A_eig (p, p, true, 0);
    for (int j=0; j<p; j++) A_eig(j,j)=sqrt(Aeig.values(j));
    Y = Y - Z * A * (t(Z) * Y);                       // save r in Y
    Z = Z * (Aeig.vectors * A_eig * t(Aeig.vectors)); // save B in Z
    //PRINTF("A\n"); for (i=0; i<p; i++) {for (j=0; j<p; j++)  PRINTF("%f ", A(i,j)); PRINTF("\n");}            
    //PRINTF("eigen values\n"); for (j=0; j<p; j++) PRINTF("%f ", Aeig.values(j)); PRINTF("\n");
    //PRINTF("B\n"); for (i=0; i<n; i++) {for (j=0; j<p; j++) PRINTF("%f ", B(i,j)); PRINTF("\n");}            
}
//            // R code:    
//            A <- solve(t(Z.sorted) %*% Z.sorted)
//            H <- Z.sorted %*% A %*% t(Z.sorted)
//            r <- as.numeric((diag(n)-H) %*% y.sorted)
//            a.eig <- eigen(A)   
//            A.sqrt <- a.eig$vectors %*% diag(sqrt(a.eig$values)) %*% solve(a.eig$vectors)
//            B = Z.sorted %*% A.sqrt  
         





#endif /* fastgrid_H */
