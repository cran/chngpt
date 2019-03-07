
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
    int order = (int) AA.rows();    // input matrix is order x_order
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


extern "C" {

//variables such as Xcusum are defined outside this function since when bootstrapping, we do not want to allocate the memory over and over again
int _twoD_search(
    Matrix<double,Row>& B, vector<double>& r, vector<double>& x,
    vector<double>& w, bool wAllOne, 
    vector<int>& stratified_by,
    int n, int p, 
    int n0, int n1,
    int nThresholds1, int nThresholds2, 
    vector<int>& thresholdIdx1, vector<int>& thresholdIdx2, 
    vector<double>& thresholds1, vector<double>& thresholds2, 
    Matrix<double,Row>& Bcusum, vector<double>& rcusum, vector<double>& Wcusum, vector<double>& xcusum,
    bool isUpperHinge,
    double* logliks)
{

    const int TH=2; // assume there are two thresholds
	
    // loop index
    int i,j,t1,t2,iB,k;    
       
    //PRINTF("thresholds1\n"); for (int i=0; i<nThresholds1; i++) PRINTF("%f ", thresholds1[i]); PRINTF("\n");
    //PRINTF("thresholds2\n"); for (int i=0; i<nThresholds2; i++) PRINTF("%f ", thresholds2[i]); PRINTF("\n"); 
        
    // save a copy of x for potential later use, e.g. if we do threshold thinning b/c x gets updated in Step 4
    vector<double> x_cpy(n);
    for(i=0; i<n; i++) x_cpy[i] = x[i];
    
    // for hinge or segmented, compute cusum of B and r in the reverse order
    // for upperhinge, compute cusum of B and r in the forward order
    if (!isUpperHinge) {
        // first group from 1 to n0: u==0
        Bcusum(n0-1,_) = B(n0-1,_)*w[n0-1]; 
        rcusum[n0-1] = r[n0-1]*w[n0-1];     
        xcusum[n0-1] = x[n0-1]*w[n0-1];     
        Wcusum[n0-1] = w[n0-1];             
        //          
        for (i=n0-2; i>=0; i--) {
            Bcusum(i,_) = Bcusum(i+1,_) + B(i,_)* w[i]; 
            rcusum[i]   = rcusum[i+1]   + r[i]* w[i]; 
            xcusum[i]   = xcusum[i+1]   + x[i]* w[i]; 
            Wcusum[i]   = Wcusum[i+1]   + w[i]; 
        }
        // second group from n0+1 to n: u==1
        Bcusum(n-1,_) = B(n-1,_)*w[n -1]; 
        rcusum[n-1] = r[n-1]*w[n -1];     
        xcusum[n-1] = x[n-1]*w[n -1];     
        Wcusum[n-1] = w[n-1];   
        //          
        for (i=n-2; i>=n0; i--) {
            Bcusum(i,_) = Bcusum(i+1,_) + B(i,_)* w[i]; 
            rcusum[i]   = rcusum[i+1]   + r[i]* w[i]; 
            xcusum[i]   = xcusum[i+1]   + x[i]* w[i]; 
            Wcusum[i]   = Wcusum[i+1]   + w[i]; 
        }
    } else {
    // upperhinge model
//        Bcusum(0,_) = B(0,_)* w[0];    
//        rcusum[0] = r(0)* w[0];        
//        Wcusum[0] = w[0];              
//        x_r_cusum[0]   = x[0] * r(0)* w[0];        
//        x_B_cusum(0,_) = x[0] * B(0,_)* w[0];    
//        for (i=1; i<n; i++) {
//            Bcusum(i,_)    = Bcusum(i-1,_)    + B(i,_)* w[i]; 
//            rcusum[i]      = rcusum[i-1]      + r[i]* w[i]; 
//            Wcusum[i]      = Wcusum[i-1]      + w[i]; 
//            x_r_cusum[i]   = x_r_cusum[i-1]   + x[i]* r[i]* w[i]; 
//            x_B_cusum(i,_) = x_B_cusum(i-1,_) + x[i]* B(i,_)* w[i]; 
//        }
    }
    //for (i=0; i<n; i++) {for (j=0; j<p; j++)  PRINTF("%f ", Bcusum(i,j)); PRINTF("\n");}        
    
    Matrix <double,Row,Concrete> vB (TH, p, true, 0), vv (TH, TH, true, 0), vr (TH, 1, true, 0);
    Matrix <double,Row,Concrete> vB2(TH, p, true, 0), vv2(TH, TH, true, 0), vr2(TH, 1, true, 0); // a second copy for updating across the second threshold
    double delta;
    
    // Step 3: Initialize vv, vr and vB 
    if (!isUpperHinge) {
        // (x-e)+
        for(i=0;  i<thresholdIdx1[0]; i++) x[i]=0;  
        for(i=thresholdIdx1[0]; i<n0; i++) x[i]=x[i]-thresholds1[0]; 
        for(i=n0; i<thresholdIdx2[0]; i++) x[i]=0;  
        for(i=thresholdIdx2[0]; i<n ; i++) x[i]=x[i]-thresholds2[0]; 
    } else {
//        // (x-e)-
//        for(i=0; i<thresholdIdx[0]; i++) x[i]=x[i]-thresholds[0];  
//        for(i=thresholdIdx[0]; i<n; i++) x[i]=0; 
    }
    //for (i=0; i<n; i++) {PRINTF("%f ", x_cpy[i]); PRINTF("%f ", xcusum[i]); PRINTF("%f ", x[i]); PRINTF("%f ", r[i]); PRINTF("\n");}
    
    if(wAllOne) {
        // we may be able to make this more efficient by changing summation index since we know some of the x's are 0
        // vv is v'v
        for (i=0 ; i<n0; i++) vv(0,0) += pow(x[i],2); 
        for (i=n0; i<n ; i++) vv(1,1) += pow(x[i],2); 
        for (i=0 ; i<n0; i++) vr(0,0) += x[i]*r[i]; 
        for (i=n0; i<n ; i++) vr(1,0) += x[i]*r[i]; 
        for (j=0; j<p; j++) {
            for (i=0; i<n0; i++) vB(0,j) += x[i] * B(i,j);
            for (i=n0; i<n; i++) vB(1,j) += x[i] * B(i,j);
        }
    } else {
    }
    // save a copy for updating across the first row
    vv2=vv; vr2=vr; vB2=vB; 
    //for (i=0; i<TH; i++) {for (j=0; j<TH; j++) PRINTF("%f ", vv(i,j)); PRINTF("\n");}
    //for (i=0; i<TH; i++) PRINTF("%f ", vr[i]); PRINTF("\n");
    //for (i=0; i<TH; i++) {for (j=0; j<p; j++) PRINTF("%f ", vB(i,j)); PRINTF("\n");}
    
    int which=0, chosen;
    double crit, crit_max;
    
    // Step 4: Compute criterion function
    crit = (t(vr) * invpd(vv - vB * t(vB)) * vr)[0];
    logliks[which++] = crit;
    chosen=which-1; crit_max=crit;
      
    // Step 5:
    for(t2=0; t2<nThresholds2; t2++) { 
        
        // update vv2/vr2/vB2 for the first row
        if(t2>0) {
            // Step 5a for updating across the first row: update vv2, vr2 and vB2
            delta= thresholds2[t2]-thresholds2[t2-1]; 
            if (!isUpperHinge) {
                iB = thresholdIdx2[t2]-1; 
                k  = n-thresholdIdx2[t2]+1; // the current threshold is the kth largest                     
                vv2(1,1) -= 2 * delta * (xcusum[iB] - k*thresholds2[t2-1]) -  k * pow(delta,2); 
                vr2(1,0) -= delta * rcusum[iB];
                for (j=0; j<p; j++) vB2(1,j) -= delta * Bcusum(iB,j);                 
            } else {
                // upperhinge model
            }
            // copy to vv/vr/vB for updating columns
            vv=vv2; vr=vr2; vB=vB2; 

            // Step 5b: Compute criterion function
            crit = (t(vr2) * invpd(vv2 - vB2 * t(vB2)) * vr2)[0];
            logliks[which++] = crit;
            if(crit>=crit_max) {chosen = which-1; crit_max=crit;}
        }
        
        // update vv/vr/vB for columns
        for(t1=1; t1<nThresholds1; t1++) {         
            // Step 5a for updating columns: update vv, vr and vB
            delta= thresholds1[t1]-thresholds1[t1-1]; 
            if (!isUpperHinge) {
                iB = thresholdIdx1[t1]-1; 
                k  = n0-thresholdIdx1[t1]+1; // the current threshold is the kth largest                     
                vv(0,0) -= 2 * delta * (xcusum[iB] - k*thresholds1[t1-1]) -  k * pow(delta,2); 
                vr(0,0) -= delta * rcusum[iB];
                for (j=0; j<p; j++) vB(0,j) -= delta * Bcusum(iB,j);                 
            } else {
//                // upperhinge model                                   
//                // e.g, when t=1, we are at the second threshold, we want to update vv, vr and vB for e_2
//                //            t=1 and t+1=2 
//                k  = thresholdIdx[t-1]; // k is 1-based index of x, e_t = x_k
//                vv -= 2*delta*(Bcusum(k-1,p-1) - k*thresholds[t-1]) - k * pow(delta,2);            
//                vr -= delta * rcusum[k-1];
//                for (j=0; j<p-1; j++) vB[j] -= delta * Bcusum(k-1,j);     
            }
            
            // Step 5b: Compute criterion function
            crit = (t(vr) * invpd(vv - vB * t(vB)) * vr)[0];
            logliks[which++] = crit;
            if(crit>=crit_max) {chosen = which-1; crit_max=crit;}                     
        } // end for t1
        
    } // end for t2

    return chosen; 
    
}

// X gets written over when the function is done
void _preprocess2(Matrix<double,Row>& X, Matrix<double,Row>& Y, vector<double>& r) {
    int i,j; 
    int p=X.cols();    
    int n=X.rows();    
    
//// R code:    for reference only
//    A <- solve(t(Z.sorted) %*% Z.sorted)
//    H <- Z.sorted %*% A %*% t(Z.sorted)
//    r <- as.numeric((diag(n)-H) %*% y.sorted)
//    a.eig <- eigen(A)   
//    A.sqrt <- a.eig$vectors %*% diag(sqrt(a.eig$values)) %*% solve(a.eig$vectors)
//    B = Z.sorted %*% A.sqrt
                
    Matrix<> A = invpd(crossprod(X));
    
    Matrix <> resid = Y - X * A * (t(X) * Y);
    // copy to r
    for (i=0; i<n; i++) r[i]=resid[i];
    
    // eigen decomposition to get B
    Eigen Aeig = eigen(A);
    Matrix <double,Row,Concrete> A_eig (p, p, true, 0);
    for (j=0; j<p; j++) A_eig(j,j)=sqrt(Aeig.values(j));
    Matrix<> B = X * (Aeig.vectors * A_eig * t(Aeig.vectors));    
    // Do this last: replace the first p column of X with B
    for (j=0; j<p; j++) X(_,j)=B(_,j);
    
    //PRINTF("A\n"); for (i=0; i<p; i++) {for (j=0; j<p; j++)  PRINTF("%f ", A(i,j)); PRINTF("\n");}            
    //PRINTF("eigen values\n"); for (j=0; j<p; j++) PRINTF("%f ", Aeig.values(j)); PRINTF("\n");
    //PRINTF("B\n"); for (i=0; i<n; i++) {for (j=0; j<p; j++) PRINTF("%f ", B(i,j)); PRINTF("\n");}            
}
         


  // For fastgrid,
  // assume X and Y are sorted 
  // thresholdIdx are 1-based index, which define the grid of thresholds
  // For fastgrid2, the meaning of the first two variables are B and r instead 
SEXP twoD_gaussian(
     SEXP u_X, 
     SEXP u_x, 
     SEXP u_Y, 
     SEXP u_W, SEXP u_wAllOne, 
     SEXP u_stratified_by, 
     //SEXP u_thresholdIdx1, SEXP u_thresholdIdx2, 
     SEXP u_lb, SEXP u_ub, 
     SEXP u_nBoot, 
     SEXP u_isUpperHinge) {
     
    double* X_dat = REAL(u_X);
    double* x_dat = REAL(u_x);
    double* Y_dat=REAL(u_Y); 
    double* W_dat=REAL(u_W);    
    bool wAllOne=asLogical(u_wAllOne)==1;
    int* stratified_by_dat = INTEGER(u_stratified_by);
//    int *thresholdIdx1=INTEGER(u_thresholdIdx1);
//    int *thresholdIdx2=INTEGER(u_thresholdIdx2);
//    int n0 = asInteger(u_n0);
//    int n1 = asInteger(u_n1);
    double lb = asReal(u_lb);
    double ub = asReal(u_ub);
    int nBoot = asInteger(u_nBoot);
    bool isUpperHinge=asLogical(u_isUpperHinge)==1;
    
    const int n = nrows(u_X);
    const int p = ncols(u_X); // number of predictors, not including the thresholed variable since u_added is separated out
//    int nThresholds1=length(u_thresholdIdx1);
//    int nThresholds2=length(u_thresholdIdx2);
    //PRINTF("thresholdIdx: "); for (int i=0; i<nThresholds; i++) PRINTF("%i ", thresholdIdx[i]); PRINTF("\n");
        
    Matrix<double,Col,Concrete> Xcol (n, p, X_dat); //The rows and colns are organized in a way now that they can be directly casted into column major, and there is no need to do things as in the JSS paper on sycthe or MCMCpack MCMCmetrop1R.cc
    Matrix<double,Row,Concrete> X(Xcol); // convert to row major so that creating bootstrap datasets can be faster and to pass to _grid_search
    Matrix<double,Row,Concrete> Y(n, 1, Y_dat); // define as a matrix instead of vector b/c will be used in matrix operation
    vector<double> x(x_dat, x_dat + n); // this does not need to be converted into vectors, but better to do so for consistency with bootstrap code
    vector<int> stratified_by(stratified_by_dat, stratified_by_dat + n); // this does not need to be converted into vectors, but better to do so for consistency with bootstrap code
    vector<double> W(W_dat, W_dat + n); // this does not need to be converted into vectors, but better to do so for consistency with bootstrap code
    vector<double> r(n); // residual vector
    
	// these variables are reused within each bootstrap replicate
    Matrix <double,Row,Concrete> Xcusum(n, p);
    vector<double> rcusum(n), Wcusum(n), xcusum(n); 
	
	int i,j;

	if (nBoot<0.1) {
    // a single search
    
        // define thresholds and thresholdIdx
        int n1=0; for(i=0; i<n; i++) n1+=stratified_by[i];
        int n0=n-n1;
        int nLower1=(int)round(n0*lb)+1, nUpper1=(int)round(n0*ub);
        int nLower2=(int)round(n1*lb)+1, nUpper2=(int)round(n1*ub);
        int nThresholds1=nUpper1-nLower1+1;
        int nThresholds2=nUpper2-nLower2+1;
        vector<int>    thresholdIdx1(nThresholds1), thresholdIdx2(nThresholds2);
    	vector<double> thresholds1  (nThresholds1), thresholds2  (nThresholds2); 
        for(i=0; i<nThresholds1; i++) {
            thresholdIdx1[i]=nLower1+i;
            thresholds1[i]=x[thresholdIdx1[i]-1];
        }
        for(i=0; i<nThresholds2; i++) {
            thresholdIdx2[i]=nLower2+i+n0;
            thresholds2[i]=x[thresholdIdx2[i]-1];
        }		
        //PRINTF("n0 %d n1 %d\n", n0, n1); 
        //PRINTF("nThresholds1 %d nThresholds2 %d\n", nThresholds1, nThresholds2); 
        //for (i=0; i<nThresholds1; i++) PRINTF("%f ", thresholds1[i]); PRINTF("\n"); 
        //for (i=0; i<nThresholds2; i++) PRINTF("%f ", thresholds2[i]); PRINTF("\n"); 
    
        SEXP _logliks=PROTECT(allocVector(REALSXP, nThresholds1*nThresholds2));
        double *logliks=REAL(_logliks);       
        
        // compute Y'HY. this is not needed to find e_hat, but good to have for comparison with other estimation methods
        // this needs to be done before step 2 since X gets changed by _preprocess2
        double yhy = ((t(Y) * X) * invpd(crossprod(X)) * (t(X) * Y))(0);
        
        // Step 1. Sort
        // Data input should be already sorted.
         
        // Step 2. Compute B and r. X is replaced by B when the function is done.
        _preprocess2(X, Y, r);
        //PRINTF("B\n"); for (int i=0; i<n; i++) {for (int j=0; j<p; j++) PRINTF("%f ", X(i,j)); PRINTF("\n");}
        //PRINTF("Y\n"); for (int i=0; i<n; i++) PRINTF("%f ", Y(i)); PRINTF("\n"); 
        
        _twoD_search(X, r, x,
                        W, wAllOne, 
                        stratified_by,
                        n, p, 
                        n0, n1,
                        nThresholds1, nThresholds2,
                        thresholdIdx1, thresholdIdx2, 
                        thresholds1, thresholds2, 
                        Xcusum, rcusum, Wcusum, xcusum, 
                        isUpperHinge,
                        logliks); 
        for(i=0; i<nThresholds1*nThresholds2; i++) logliks[i]=logliks[i]+yhy; 
        //PRINTF("yhy: %.4f\n", yhy);  
        //for(int t1=0; t1<nThresholds1; t1++) { {for(int t2=0; t2<nThresholds2; t2++) PRINTF("%.4f ", logliks[t1+n1*t2]);} PRINTF("\n");}
            
        UNPROTECT(1);
        return _logliks;

    } else {
    // bootstrap
    // difference from single search is that we are not return liklihoods but parameter estimates from bootstrap datasets
    
        // output variables
        // logliks will not be returned to R, estimates from each bootstrap copy will be stored in coef and returned
        // to avoid allocating space for logliks for each bootstrap replicate, we allocate a big enough space just once
        double * logliks = (double *) malloc((n*n) * sizeof(double));
        SEXP _coef=PROTECT(allocVector(REALSXP, nBoot*(p+2+2)));// p slopes, 2 thresholds and 2 threshold-related slopes
        double *coef=REAL(_coef);    
        
    	// these variables are reused within each bootstrap replicate
    	vector<int> index(n);
        Matrix <double,Row,Concrete> Xb(n,p), Yb(n,1), Xeb(n,p+2);
    	vector<double> xb(n), Wb(n), rb(n); 
    	vector<int> stratified_by_b(n);
    	double e_hat, f_hat;
    	int chosen;
    	
        for (int b=0; b<nBoot; b++) {        
            // create bootstrap dataset, note that index is 1-based
            SampleReplace(n, n, &(index[0]));
            
            // Step 1: sort
            sort (index.begin(), index.end());
            
            for (i=0; i<n; i++) { //note that index need to -1 to become 0-based
                Xb(i,_)=X(index[i]-1,_); 
                Yb[i]  =Y[index[i]-1]; 
                xb[i]  =x[index[i]-1]; 
                Wb[i]  =W[index[i]-1]; 
                stratified_by_b[i]=stratified_by[index[i]-1]; 
            } 
            //for (i=0; i<n; i++) PRINTF("%f ", xb[i]); PRINTF("\n");
            //for (i=0; i<n; i++) PRINTF("%d ", stratified_by_b[i]); PRINTF("\n");
            //for (i=0; i<n; i++) {for (j=0; j<p; j++)  PRINTF("%f ", Xb(i,j)); PRINTF("\n");} 
            //for (i=0; i<n; i++) { Xb(i,_)=X(i,_); Yb(i)=Y(i); Wb[i]=W[i]; } // debug use, can be used to compare with non-boot
            
            // define thresholds and thresholdIdx
            // these need to be redefined for every bootstrap replicate b/c n0 and n1 could change in bootstrap
            int n1=0; for(i=0; i<n; i++) n1+=stratified_by_b[i];
            int n0=n-n1;
            int nLower1=(int)round(n0*lb)+1, nUpper1=(int)round(n0*ub);
            int nLower2=(int)round(n1*lb)+1, nUpper2=(int)round(n1*ub);
            int nThresholds1=nUpper1-nLower1+1;
            int nThresholds2=nUpper2-nLower2+1;
            vector<int>    thresholdIdx1(nThresholds1), thresholdIdx2(nThresholds2);
        	vector<double> thresholds1  (nThresholds1), thresholds2  (nThresholds2); 
            for(i=0; i<nThresholds1; i++) {
                thresholdIdx1[i]=nLower1+i;
                thresholds1[i]=xb[thresholdIdx1[i]-1];
            }
            for(i=0; i<nThresholds2; i++) {
                thresholdIdx2[i]=nLower2+i+n0;
                thresholds2[i]=xb[thresholdIdx2[i]-1];
            }		
            //PRINTF("n0 %d n1 %d\n", n0, n1); 
            //PRINTF("nThresholds1 %d nThresholds2 %d\n", nThresholds1, nThresholds2); 
            //for (i=0; i<nThresholds1; i++) PRINTF("%f ", thresholds1[i]); PRINTF("\n"); 
            //for (i=0; i<nThresholds2; i++) PRINTF("%f ", thresholds2[i]); PRINTF("\n"); 
    
            // Step 2:
            _preprocess2(Xb, Yb, rb);            
            
            // Step 3-5
            chosen = _twoD_search(Xb, rb, xb,
                        Wb, wAllOne, 
                        stratified_by_b,
                        n, p, 
                        n0, n1, 
                        nThresholds1, nThresholds2, 
                        thresholdIdx1, thresholdIdx2, 
                        thresholds1, thresholds2, 
                        Xcusum, rcusum, Wcusum, xcusum, 
                        isUpperHinge,
                        logliks); 
            e_hat=thresholds1[chosen % nThresholds1];
            f_hat=thresholds2[chosen / nThresholds1];
            //PRINTF("chosen %d\n", chosen); 
            //PRINTF("e_hat %f\n", e_hat); PRINTF("f_hat %f\n", f_hat); 

            // fit model at the selected threshold and save results in coef
            
            // after _preprocess2, Xb have changed, thus we need to copy from X again
            for (i=0; i<n; i++) { 
                for (j=0; j<p; j++)
                    Xeb(i,j)=X(index[i]-1,j); 
            }
            // after search, xb have changed, thus we need to copy from x again
            if (!isUpperHinge) {
                // (x-e)+
                for(i=0; i<n0; i++) {
                    Xeb(i,p) = x[index[i]-1]<e_hat?0:x[index[i]-1]-e_hat;  
                    Xeb(i,p+1) = 0; 
                }
                for(i=n0; i<n; i++) {
                    Xeb(i,p) = 0; 
                    Xeb(i,p+1) = x[index[i]-1]<f_hat?0:x[index[i]-1]-f_hat;                  
                }
            } else {
//                // (x-e)-
//                for(i=0; i<n; i++) 
//                    if(Xb(i,p-1)<e_hat) Xb(i,p-1) -= e_hat; else Xb(i,p-1) = 0;  
            }
            //PRINTF("Xeb\n"); for (i=0; i<n; i++) {for (j=0; j<p+2; j++)  PRINTF("%f ", Xeb(i,j)); PRINTF("\n");} 
            //PRINTF("Yb\n"); for (i=0; i<n; i++) PRINTF("%f ", Yb(i)); PRINTF("\n");
            //PRINTF("Xb\n"); for (i=0; i<n; i++) {for (j=0; j<p; j++)  PRINTF("%f ", Xb(i,j)); PRINTF("\n");} 
            //PRINTF("Yb\n"); for (i=0; i<n; i++) PRINTF("%f ", Yb(i)); PRINTF("\n");
//            for (i=0; i<n; i++) { Xb(i,_)=X(i,_); Yb(i)=Y(i); } // debug use
            
                        
            Matrix <> beta_hat = invpd(crossprod(Xeb)) * (t(Xeb) * Yb);
            //PRINTF("beta_jat\n"); for (i=0; i<p+2; i++) PRINTF("%f ", beta_hat[i]); PRINTF("\n");
            
            for (j=0; j<p+2; j++) coef[b*(p+2)+j]=beta_hat[j];            
            coef[b*(p+2)+p+2] = e_hat;
            coef[b*(p+2)+p+3] = f_hat;

        } 
         
        UNPROTECT(1);
        free(logliks);
        return _coef;
    }
    
}

}  // end extern C

#endif
//#endif
