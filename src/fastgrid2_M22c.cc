
//////////////////////////////////////////////////////////////////////////
// 
// This software is distributed under the terms of the GNU GENERAL
// PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
// file for more information.
//
//////////////////////////////////////////////////////////////////////////

#ifndef CHNGPTBOOT_CC
#define CHNGPTBOOT_CC


// the following leads to many problem
//#ifndef SCYTHE_LAPACK 
//#define SCYTHE_LAPACK


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


using namespace std;
using namespace scythe;

// copied from random.c, sample with replacement
static void SampleReplace(int k, int n, int *y)
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


inline void make_symmetric(double* matrix, int rows)
{
      for (int i = 1; i < rows; ++i)
        for (int j = 0; j < i; ++j)
          matrix[i * rows + j] = matrix[j * rows + i];
}
  // it is not clear to me whether crossprod is calling lapack or not. crossprod1 is the way I make sure it is
  // if a row-major matrix is passed as A, it will be transposed automatically
inline Matrix<> crossprod1(const Matrix<>& A)
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


// copied from ide.h

struct Eigen {
   Matrix<> values;
   Matrix<> vectors;
};

inline Eigen
eigen (const Matrix<>& A, bool vectors=true)
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

// first p-1 col of X and Y are changed
inline void _preprocess(Matrix<double,Row>& X, Matrix<double,Row>& Y) {
    int i,j; 
    int p=X.cols();    
    int n=X.rows();    
    
//// R code:    
//    A <- solve(t(Z.sorted) %*% Z.sorted)
//    H <- Z.sorted %*% A %*% t(Z.sorted)
//    r <- as.numeric((diag(n)-H) %*% y.sorted)
//    a.eig <- eigen(A)   
//    A.sqrt <- a.eig$vectors %*% diag(sqrt(a.eig$values)) %*% solve(a.eig$vectors)
//    B = Z.sorted %*% A.sqrt
                
    Matrix<> Z=X(0,0,X.rows()-1,p-2);
    Matrix<> A = invpd(crossprod(Z));
    
    // eigen decomposition to get B
    Eigen Aeig = eigen(A);
    Matrix <double,Row,Concrete> A_eig (p-1, p-1, true, 0);
    for (j=0; j<p-1; j++) A_eig(j,j)=sqrt(Aeig.values(j));
    Matrix<> B = Z * (Aeig.vectors * A_eig * t(Aeig.vectors));
    
    // replace the first p-1 column of X with B
    for (j=0; j<p-1; j++) X(_,j)=B(_,j);
    //PRINTF("A\n"); for (i=0; i<p-1; i++) {for (j=0; j<p-1; j++)  PRINTF("%f ", A(i,j)); PRINTF("\n");}            
    //PRINTF("eigen values\n"); for (j=0; j<p-1; j++) PRINTF("%f ", Aeig.values(j)); PRINTF("\n");
    //PRINTF("B\n"); for (i=0; i<n; i++) {for (j=0; j<p-1; j++) PRINTF("%f ", B(i,j)); PRINTF("\n");}            
    
    Matrix <> r = Y - Z * A * (t(Z) * Y);
    // replace Y with r
    for (i=0; i<n; i++) Y(i)=r(i);
}
         


extern "C" {
       
double _fastgrid2M22c_search(Matrix<double,Row>& X, Matrix<double,Row>& Y, vector<double>& w, bool wAllOne, 
    int * thresholdIdx, bool skipping, vector<double>& thresholds, 
    int n, int p, int nThresholds, 
    Matrix<double,Row>& Xcusum, vector<double>& Ycusum, vector<double>& Wcusum, vector<double>& xcusum, vector<double>& x2cusum, vector<double>& x3cusum, vector<double>& xrcusum, Matrix<double,Row>& xBcusum, 
    Matrix<double,Row>& Xcusum2, vector<double>& Ycusum2, vector<double>& Wcusum2, vector<double>& xcusum2, vector<double>& x2cusum2, vector<double>& x3cusum2, vector<double>& xrcusum2, Matrix<double,Row>& xBcusum2, 
//    Matrix<double,Row>& Xcpy, Matrix<double,Row>& Ycpy, // debug use
    double* logliks)
{

    // loop index. 
    int i,j,tt;    
    int k; // k has a meaning in the algorithm 
	int chosen=0;//initialize to get rid of compilation warning
	
    // covariate to be thresholed
    vector<double> x_cpy(n);
    for(i=0; i<n; i++) x_cpy[i] = X(i,p-1);        
    //for (i=0; i<n; i++) PRINTF("%f ", x_cpy[i]); PRINTF("\n");
        
    // compute cusum of X and Y in the forward order (upper hinge)
    Xcusum(0,_) = X(0,_)* w[0];    
    Ycusum[0] = Y(0)* w[0];        
    Wcusum[0] = w[0];              
    xcusum[0]    = x_cpy[0] * w[0]; // xcusum is the last columnn of Xcusum, but perhaps better to be its own variable
    x2cusum[0]   = pow(x_cpy[0],2) * w[0];        
    x3cusum[0]   = pow(x_cpy[0],3) * w[0];        
    xrcusum[0]   = x_cpy[0] * Y(0)* w[0];        
    xBcusum(0,_) = x_cpy[0] * X(0,_)* w[0];    
    for (i=1; i<n; i++) {
        Xcusum(i,_)  = Xcusum(i-1,_)    + X(i,_)* w[i]; 
        Ycusum[i]    = Ycusum[i-1]      + Y(i)* w[i]; 
        Wcusum[i]    = Wcusum[i-1]      + w[i]; 
        xcusum[i]    = xcusum[i-1]      + x_cpy[i]* w[i]; 
        x2cusum[i]   = x2cusum[i-1]     + pow(x_cpy[i],2)* w[i]; 
        x3cusum[i]   = x3cusum[i-1]     + pow(x_cpy[i],3)* w[i]; 
        xrcusum[i]   = xrcusum[i-1]     + x_cpy[i]* Y(i)* w[i]; 
        xBcusum(i,_) = xBcusum(i-1,_)   + x_cpy[i]* X(i,_)* w[i]; 
    }
    
    // compute cusum of X and Y in the reverse order (hinge)
    for (i=n-1; i>=1; i--) {
        Xcusum2(i,_)  = Xcusum(n-1,_)    - Xcusum(i-1,_);
        Ycusum2[i]    = Ycusum[n-1]      - Ycusum[i-1]; 
        Wcusum2[i]    = Wcusum[n-1]      - Wcusum[i-1]; 
        xcusum2[i]    = xcusum[n-1]      - xcusum[i-1]; 
        x2cusum2[i]   = x2cusum[n-1]     - x2cusum[i-1]; 
        x3cusum2[i]   = x3cusum[n-1]     - x3cusum[i-1]; 
        xrcusum2[i]   = xrcusum[n-1]     - xrcusum[i-1]; 
        xBcusum2(i,_) = xBcusum(n-1,_)   - xBcusum(i-1,_); 
    }
    Xcusum2(0,_)  = Xcusum(n-1,_);
    Ycusum2[0]    = Ycusum[n-1]; 
    Wcusum2[0]    = Wcusum[n-1]; 
    xcusum2[0]    = xcusum[n-1]; 
    x2cusum2[0]   = x2cusum[n-1]; 
    x3cusum2[0]   = x3cusum[n-1]; 
    xrcusum2[0]   = xrcusum[n-1]; 
    xBcusum2(0,_) = xBcusum(n-1,_); 


    //for (i=0; i<n; i++) {for (j=0; j<p; j++)  PRINTF("%f ", Xcusum(i,j)); PRINTF("\n");}        
    //for (i=0; i<n; i++) PRINTF("%f ", Ycusum[i]); PRINTF("\n");
    //for (i=0; i<n; i++) {for (j=0; j<p; j++)  PRINTF("%f ", Xcusum2(i,j)); PRINTF("\n");}        
    //for (i=0; i<n; i++) PRINTF("%f ", Ycusum2[i]); PRINTF("\n");
        
    // Step 2: Initialize VV, Vr and VB     
    Matrix <double,Row,Concrete> VV(3, 3), Vr(3, 1), VB(3, p-1);    
    VV(1,2)=VV(2,1)=0;
    // x-e
    for(i=0; i<n; i++) X(i,p-1)=x_cpy[i]-thresholds[0]; 
    if(wAllOne) {
        VV(0,0)=0; for (i=0; i<n; i++) VV(0,0) += pow(X(i,p-1),2); 
        
        Vr(0,0)=0; for (i=0; i<n; i++) Vr(0,0) += X(i,p-1)*Y(i); 
        
        for (j=0; j<p-1; j++) {
            VB(0,j)=0; for (i=0; i<n; i++) VB(0,j) += X(i,p-1) * X(i,j); // remember the first p-1 column of X is B after _preprocess()
        }
    } else {
    }
    // threshold x (stored in the last column of X) at the first threshold value
    // Upper hinge (x-e)-
    for(i=0; i<thresholdIdx[0]; i++) X(i,p-1)=x_cpy[i]-thresholds[0];
    for(i=thresholdIdx[0]; i<n; i++) X(i,p-1)=0; 
    if(wAllOne) {
        VV(1,1)=0; for (i=0; i<n; i++) VV(1,1) += pow(X(i,p-1),4); 
        VV(0,1)=0; for (i=0; i<n; i++) VV(0,1) += pow(X(i,p-1),3);
        VV(1,0)=VV(0,1);
        
        Vr(1,0)=0; for (i=0; i<n; i++) Vr(1,0) += pow(X(i,p-1),2)*Y(i); 
        
        for (j=0; j<p-1; j++) {
            VB(1,j)=0; for (i=0; i<n; i++) VB(1,j) += pow(X(i,p-1),2) * X(i,j); // remember the first p-1 column of X is B after _preprocess()
        }
    } else {
    }
    // threshold x (stored in the last column of X) at the first threshold value
    // Hinge (x-e)+
    for(i=0; i<thresholdIdx[0]; i++) X(i,p-1)=0;  
    for(i=thresholdIdx[0]; i<n; i++) X(i,p-1)=x_cpy[i]-thresholds[0]; 
    if(wAllOne) {
        VV(2,2)=0; for (i=0; i<n; i++) VV(2,2) += pow(X(i,p-1),4); 
        VV(0,2)=0; for (i=0; i<n; i++) VV(0,2) += pow(X(i,p-1),3);
        VV(2,0)=VV(0,2);
        
        Vr(2,0)=0; for (i=0; i<n; i++) Vr(2,0) += pow(X(i,p-1),2)*Y(i); 
        
        for (j=0; j<p-1; j++) {
            VB(2,j)=0; for (i=0; i<n; i++) VB(2,j) += pow(X(i,p-1),2) * X(i,j); // remember the first p-1 column of X is B after _preprocess()
        }
    } else {
    }
    //for (i=0; i<3; i++) {for (j=0; j<3; j++)   PRINTF("%f ", VV(i,j)); PRINTF("\n");}        
    //for (i=0; i<3; i++) {for (j=0; j<1; j++)   PRINTF("%f ", Vr(i,j)); PRINTF("\n");}        
    //for (i=0; i<3; i++) {for (j=0; j<p-1; j++) PRINTF("%f ", VB(i,j)); PRINTF("\n");}        
    
    // Step 4 (Step 3 = Step 4b)
    double delta, s, crit, crit_max=R_NegInf, e;// s is defined as difference in squared threshold
    for(tt=0; tt<nThresholds; tt++) { 
            
        // Step 4a: update VV, Vr and VB     
        if(tt>0) {
            e=thresholds[tt-1];
            delta= thresholds[tt]-e; //X(tt,p-1)-X(tt-1,p-1); //PRINTF("%f ", delta); PRINTF("\n"); 
            s=pow(thresholds[tt],2)-pow(e, 2);
            
            // x-e part
            VV(0,0) += -2*delta*xcusum[n-1] + 2*n*delta*e + n * pow(delta,2); 
            Vr[0] += -delta * Ycusum[n-1];
            for (j=0; j<p-1; j++) {
                VB(0,j) +=   -delta * Xcusum(n-1,j); 
            }
            
            // upper hinge part
            // e.g, when tt=1, we are at the second threshold, we want to update vv, vr and vB for e_2
            //            tt=1 and tt+1=2 
            k  = thresholdIdx[tt-1]; // k is 1-based index of x, e_t = x_k
            VV(1,1) += -4*delta*x3cusum[k-1] + 2*(2*pow(delta,2)+4*delta*e+s)*x2cusum[k-1] - (4*delta*s+4*delta*pow(e,2)+4*e*s)*xcusum[k-1] + k*s*s + 2*k*s*pow(e,2); 
            VV(0,1) += -3*delta*x2cusum[k-1] + (2*pow(delta,2)+4*delta*e+s) * xcusum[k-1] - k*delta*s - k*s*e - k*delta*pow(e,2); 
            VV(1,0) = VV(0,1);
            Vr[1] += -2*delta * xrcusum[k-1] + s * Ycusum[k-1] ;
            for (j=0; j<p-1; j++) {
                VB(1,j) += -2*delta * xBcusum(k-1,j) + s * Xcusum(k-1,j); 
            }

            // hinge part
            k  = n-k; // the current threshold is the kth largest 
            VV(2,2) += -4*delta*x3cusum2[n-k] + 2*(2*pow(delta,2)+4*delta*e+s)*x2cusum2[n-k] - (4*delta*s+4*delta*pow(e,2)+4*e*s)*xcusum2[n-k] + k*s*s + 2*k*s*pow(e,2); 
            VV(0,2) += -3*delta*x2cusum2[n-k] +   (2*pow(delta,2)+4*delta*e+s) *xcusum2[n-k] -  k*delta*s - k*s*e - k*delta*pow(e,2); 
            VV(2,0) = VV(0,2);
            Vr[2] += -2*delta * xrcusum2[n-k] + s * Ycusum2[n-k];
            for (j=0; j<p-1; j++) {
                VB(2,j) += -2*delta * xBcusum2(n-k,j) + s * Xcusum2(n-k,j); 
            }                        
        }
        
        if(tt==1) {
            //for (i=0; i<3; i++) {for (j=0; j<3; j++)   PRINTF("%f ", VV(i,j)); PRINTF("\n");}        
            //for (i=0; i<3; i++) {for (j=0; j<1; j++)   PRINTF("%f ", Vr(i,j)); PRINTF("\n");}        
            //for (i=0; i<3; i++) {for (j=0; j<p-1; j++) PRINTF("%f ", VB(i,j)); PRINTF("\n");}        
        }
            
        // Step 4b: compute Y'H_eY - Y'HY
        crit = (t(Vr) * invpd(VV - VB * t(VB)) * Vr) (0,0);
        //for debugging, index in the if condition can be changed from 0 to 1, ...
//        if(tt==1) {
////            Matrix<double,Row,Concrete> Xe(n,p+1);
////            for (i=0; i<n; i++) for (j=0; j<p-1; j++) Xe(i,j) = Xcpy(i,j); 
////            for (i=0; i<n; i++) Xe(i,p-1) = X(i,p-1); 
////            for (i=0; i<n; i++) Xe(i,p) = pow(X(i,p-1),2); 
////            //for (i=0; i<n; i++) {for (j=0; j<p+1; j++)  PRINTF("%f ", Xe(i,j)); PRINTF("\n");}        
////            Matrix<> Z= Xe(0,0,Xe.rows()-1,p-2);
////            Matrix<> Ve=Xe(0,p-1,Xe.rows()-1,p);
//            //for (i=0; i<n; i++) {for (j=0; j<2; j++)  PRINTF("%f ", Ve(i,j)); PRINTF("\n");}                    
//            
//            //Matrix<> tmpM = t(Ve) * Z * invpd(t(Z) * Z) * t(Z) * Ve;
//            //PRINTF("%f %f %f %f", tmpM(0,0), tmpM(0,1), tmpM(1,0), tmpM(1,1)); PRINTF("\n");
//            Matrix<> tmpM1=VB * t(VB);     
//            PRINTF("%f %f %f %f", tmpM1(0,0), tmpM1(0,1), tmpM1(1,0), tmpM1(1,1)); PRINTF("\n"); 
//            Matrix<> tmpM2=VV - VB * t(VB);     
//            PRINTF("%f %f %f %f", tmpM2(0,0), tmpM2(0,1), tmpM2(1,0), tmpM2(1,1)); PRINTF("\n"); 
//            
//            //double tmp = (t(Ycpy) * Xe * invpd(t(Xe) * Xe) * t(Xe) * Ycpy) (0,0);
//            //PRINTF("Y H_e Y %f ", tmp); PRINTF("\n");        
//        }
            
        logliks[tt] = crit;
        if(crit>=crit_max) {
            chosen = tt;
            crit_max=crit;
        }             
    }        
    //PRINTF("logliks: \n");  for(i=0; i<nThresholds; i++) PRINTF("%f ", logliks[i]); PRINTF("\n");  

    return thresholds[chosen]; 
    
}
       

  // assume X and Y are sorted in chngptvar from small to large
  // assume last col of X is chngptvar, which will be updated as move through the grid
  // thresholdIdx are 1-based index, which define the grid of thresholds
  // For fastgrid2, the meaning of the first two variables are B and r instead 
SEXP fastgrid2M22c_gaussian(
     SEXP u_X, SEXP u_Y, 
     SEXP u_W, SEXP u_wAllOne, 
     SEXP u_thresholdIdx, SEXP u_skipping,
     SEXP u_nBoot, SEXP u_nSub)
{
    double* X_dat = REAL(u_X);
    double* Y_dat=REAL(u_Y); 
    double* W_dat=REAL(u_W);    
    bool wAllOne=asLogical(u_wAllOne)==1;
    int *thresholdIdx=INTEGER(u_thresholdIdx);
    bool skipping=asLogical(u_skipping)==1;
    int nBoot = asInteger(u_nBoot);
    
    const int n = nrows(u_X);
    const int p = ncols(u_X); // number of predictors, including the thresholded variable, for quad, it is actually 1-less than the predictors
    int nThresholds=length(u_thresholdIdx);
    //PRINTF("thresholdIdx: "); for (int i=0; i<nThresholds; i++) PRINTF("%i ", thresholdIdx[i]); PRINTF("\n");
        
    // The rows and colns are organized in a way now that they can be directly casted and there is no need to do things as in the JSS paper on sycthe or MCMCpack MCMCmetrop1R.cc
    Matrix<double,Col,Concrete> Xcol (n, p, X_dat); //column major     
    Matrix<double,Row,Concrete> X(Xcol); // convert to row major so that creating bootstrap datasets can be faster and to pass to _grid_search
    Matrix<double,Row,Concrete> Y(n, 1, Y_dat); // define as a matrix instead of vector b/c will be used in matrix operation
    vector<double> W(W_dat, W_dat + n);
//    Matrix<double,Row,Concrete> Xcpy(X), Ycpy(Y); // debug use
    
	// these variables are reused within each bootstrap replicate
    Matrix <double,Row,Concrete> Xcusum (n, p), xBcusum (n, p);
	vector<double> Ycusum(n), Wcusum(n), xcusum (n), x2cusum (n), x3cusum (n), xrcusum (n); 
    Matrix <double,Row,Concrete> Xcusum2 (n, p), xBcusum2 (n, p);
	vector<double> Ycusum2(n), Wcusum2(n), xcusum2 (n), x2cusum2 (n), x3cusum2 (n), xrcusum2 (n); 
	vector<double> thresholds(nThresholds); 
	
	int i,j;
 
	
    if (nBoot<0.1) {
    // a single search
    
        SEXP _logliks=PROTECT(allocVector(REALSXP, nThresholds));
        double *logliks=REAL(_logliks);       
        
        double yhy=0;
        if (p>1) {
        // if p is 0, we cannot do what is in this condition
            // compute Y'HY. this is not needed to find e_hat, but good to have for comparison with other estimation methods
            Matrix<> Z=X(0,0,X.rows()-1,p-2);
            yhy = ((t(Y) * Z) * invpd(crossprod(Z)) * (t(Z) * Y))(0);
            //PRINTF("yhy %f ", yhy); PRINTF("\n"); 
             
            // Compute B and r
            _preprocess(X, Y);
            //PRINTF("B\n"); for (int i=0; i<n; i++) {for (int j=0; j<p; j++) PRINTF("%f ", X(i,j)); PRINTF("\n");}
            //PRINTF("Y\n"); for (int i=0; i<n; i++) PRINTF("%f ", Y(i)); PRINTF("\n"); 
        }
            
        for(i=0; i<nThresholds; i++) thresholds[i]=X(thresholdIdx[i]-1,p-1);

        _fastgrid2M22c_search(X, Y, W, wAllOne, thresholdIdx, skipping, thresholds, 
                             n, p, nThresholds,
                             Xcusum, Ycusum, Wcusum, xcusum, x2cusum, x3cusum, xrcusum, xBcusum, 
                             Xcusum2, Ycusum2, Wcusum2, xcusum2, x2cusum2, x3cusum2, xrcusum2, xBcusum2, 
//                             Xcpy, Ycpy, // debug use
                             logliks); 
        for(i=0; i<nThresholds; i++) logliks[i]=logliks[i]+yhy; 
        //PRINTF("logliks : \n");  for(i=0; i<nThresholds; i++) PRINTF("%f ", logliks[i]); PRINTF("\n");  
            
        UNPROTECT(1);
        return _logliks;

    } else {
    // bootstrap
    
      //  //output variables: logliks will not be returned to R, estimates from each bootstrap copy will be stored in coef and returned
        double * logliks = (double *) malloc((nThresholds) * sizeof(double));
        SEXP _coef=PROTECT(allocVector(REALSXP, nBoot*(p+2+1)));// p+2 slopes, 1 threshold
        double *coef=REAL(_coef);    
        
    	// these variables are reused within each bootstrap replicate
    	vector<int> index(n);
        Matrix <double,Row,Concrete> Xb(n,p), Yb(n,1), Xbreg(n,p+2); // Xbreg is used in regression, it contains the quadratic term
    	vector<double> Wb(n); 
    	double e_hat;
    	
        for (int b=0; b<nBoot; b++) {        
            // create bootstrap dataset, note that index is 1-based
            SampleReplace(n, n, &(index[0]));
            // Step 1: sort
            sort (index.begin(), index.end());
            
            for (i=0; i<n; i++) { //note that index need to -1 to become 0-based
                Xb(i,_)=X(index[i]-1,_); 
                Yb(i)  =Y(index[i]-1); 
                Wb[i]  =W[index[i]-1]; 
            } 
//            for (i=0; i<n; i++) { Xb(i,_)=X(i,_); Yb(i)=Y(i); Wb[i]=W[i]; } // debug use, can be used to compare with non-boot
            //for (i=0; i<n; i++) {for (j=0; j<p; j++)  PRINTF("%f ", Xb(i,j)); PRINTF("\n");} 
            
            for(i=0; i<nThresholds; i++) thresholds[i]=Xb(thresholdIdx[i]-1,p-1);
            
            // Step 2: Compute B and r
            if (p>1) {
            // if p is 0, we cannot do what is in this condition
                _preprocess(Xb, Yb);            
            }
            e_hat = _fastgrid2M22c_search(Xb, Yb, Wb, wAllOne, thresholdIdx, skipping, thresholds, 
                                 n, p, nThresholds,
                                 Xcusum, Ycusum, Wcusum, xcusum, x2cusum, x3cusum, xrcusum, xBcusum, 
                                 Xcusum2, Ycusum2, Wcusum2, xcusum2, x2cusum2, x3cusum2, xrcusum2, xBcusum2, 
//                             Xcpy, Ycpy, // debug use
                                 logliks); 

            //****** fit model at the selected threshold and save results in coef            
            // since after _preprocess, Xb and Yb have changed to B and r, we need to put X and Y back
            // We put it in Xbreg because we need quadratic term
            for (i=0; i<n; i++) { //note that index need to -1 to become 0-based
                Yb(i)=Y(index[i]-1); 
                for (j=0; j<p; j++) Xbreg(i,j)=X(index[i]-1,j);
                // x-e
                Xbreg(i,p-1) -= e_hat; 
                // (x-e)-^2
                if(X(index[i]-1,p-1)<e_hat) Xbreg(i,p)  =pow(Xbreg(i,p-1),2); else Xbreg(i,p)=0;
                // (x-e)+^2
                if(X(index[i]-1,p-1)>e_hat) Xbreg(i,p+1)=pow(Xbreg(i,p-1),2); else Xbreg(i,p+1)=0;
            }     
            Matrix <> beta_hat = invpd(crossprod(Xbreg)) * (t(Xbreg) * Yb);
            
            for (j=0; j<p+2; j++) coef[b*(p+3)+j]=beta_hat[j];             
            coef[b*(p+3)+p+2] = e_hat;       
                 
            //PRINTF("e_hat %f\n", e_hat); 
            //PRINTF("beta_hat\n"); for (i=0; i<p+1+2; i++) PRINTF("%f ", beta_hat(i)); PRINTF("\n");
            //PRINTF("Xbreg_e\n"); for (i=0; i<n; i++) {for (j=0; j<p+1+2; j++)  PRINTF("%f ", Xbreg(i,j)); PRINTF("\n");} 
            //PRINTF("Yb\n"); for (i=0; i<n; i++) PRINTF("%f ", Yb(i)); PRINTF("\n");
            
        } 
         
        UNPROTECT(1);
        free(logliks);
        return _coef;
    }
    
}

}  // end extern C

#endif
//#endif
