
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

// taken from R's src/main/random.c
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
       

// X is actually cbind(B,x) and Y is actually r
  //variables such as Xcusum are defined outside this function since when bootstrapping, we do not want to allocate the memory over and over again
double _fastgrid2step_search(Matrix<double,Row>& X, Matrix<double,Row>& Y, vector<double>& w, bool wAllOne, 
    int * thresholdIdx, bool skipping, vector<double>& thresholds, 
    int n, int p, int nThresholds,
    Matrix<double,Row>& Xcusum, vector<double>& Ycusum, vector<double>& Wcusum, 
    double* logliks)
{

    // loop index. 
    int i,j,t;    
    int iX; 
	int chosen=0;//initialize to get rid of compilation warning
	
    vector<double> x_cpy(n);
   	vector<double> x_r_cusum (n);
    // store additional Cusum
    Matrix <double,Row,Concrete> x_B_cusum (n, p); // x^2 is the last column
	if(skipping) { // guess if the thresholds are thinned. we could have change the .Call interface, but it is too much work
        // save a copy of x for later use
        for(i=0; i<n; i++) x_cpy[i] = X(i,p-1);        
    }
        
    Xcusum(n-1,_) = X(n-1,_)*w[n-1]; 
    Ycusum[n-1] = Y(n-1)*w[n-1];     
    Wcusum[n-1] = w[n-1];             
    for (i=n-2; i>=0; i--) {
        Xcusum(i,_)    = Xcusum(i+1,_)    + X(i,_)* w[i]; 
        Ycusum[i]      = Ycusum[i+1]      + Y(i)* w[i]; 
        Wcusum[i]      = Wcusum[i+1]      + w[i]; 
	}
    if(skipping) {
        x_r_cusum[n-1]   = x_cpy[n-1]*Y(n-1)   * w[n-1]; 
        x_B_cusum(n-1,_) = x_cpy[n-1]*X(n-1,_) * w[n-1]; 
        for (i=n-2; i>=0; i--) {
            x_r_cusum[i]   = x_r_cusum[i+1]   + x_cpy[i]*Y(i)* w[i]; 
            x_B_cusum(i,_) = x_B_cusum(i+1,_) + x_cpy[i]*X(i,_)* w[i]; 
        }
    }    
        

    
    //for (i=0; i<n; i++) {for (j=0; j<p; j++)  PRINTF("%f ", Xcusum(i,j)); PRINTF("\n");}        
    //for (i=0; i<n; i++) PRINTF("%f ", Ycusum[i]); PRINTF("\n");
    
    vector<double> vB(p-1); 
    double vv, vr;
    double delta, crit, crit_max=R_NegInf;
    
    // Step 2: Initialize vv, vr and vB 
    
    // threshold x, which is the last column of X, at the first threshold value
    // (x-e)+
    for(i=0; i<thresholdIdx[0]; i++) X(i,p-1)=0;  
    for(i=thresholdIdx[0]; i<n; i++) X(i,p-1)=1; 
    //for (j=0; j<n; j++) PRINTF("%f ", X(j,p-1)); PRINTF("\n");
    
    vv=0; vr=0; 
    if(wAllOne) {
        for (i=0; i<n; i++) vv += pow(X(i,p-1),2); 
        for (i=0; i<n; i++) vr += X(i,p-1)*Y(i); 
        for (j=0; j<p-1; j++) {
            vB[j]=0; 
            for (i=0; i<n; i++) 
                vB[j] += X(i,p-1) * X(i,j);
        }
    } else {
    }
    
    // Step 4 (Step 3 = Step 4b)
    for(t=0; t<nThresholds; t++) { 
            
        // Step 4a: update vv, vr and vB
        if(t>0) {
            delta= thresholds[t]-thresholds[t-1]; //X(t,p-1)-X(t-1,p-1); //PRINTF("%f ", delta); PRINTF("\n"); 
            // only update when delta>0
            if(delta>0) {
                iX = thresholdIdx[t]-1; 
                //k  = n-thresholdIdx[t]+1; // the current threshold is the kth largest             
                vv -= 1; 
                vr -= Y[iX];
                for (j=0; j<p-1; j++) vB[j] -= X(iX,j);             
            }
        }
        
        // Step 4b: compute Y'H_eY - Y'HY
        // first compute vB'vB, abusing the notation a little bit
        crit=0; for (j=0; j<p-1; j++) crit += pow(vB[j],2); 
        // now compute crit
        crit = vr*vr / (vv - crit);
        logliks[t] = crit;
        if(crit>=crit_max) {
            chosen = t;
            crit_max=crit;
        }             
    }        
    //PRINTF("logliks: \n");  for(i=0; i<nThresholds; i++) PRINTF("%f ", logliks[i]); PRINTF("\n");  

    return thresholds[chosen]; 
    
}


  // For fastgrid,
  // assume X and Y are sorted in chngptvar from small to large
  // assume last col of X is chngptvar, which will be updated as move through the grid
  // thresholdIdx are 1-based index, which define the grid of thresholds
  // For fastgrid2, the meaning of the first two variables are B and r instead 
SEXP fastgrid2step_gaussian(
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
    int nSub = asInteger(u_nSub);
    
    const int n = nrows(u_X);
    const int p = ncols(u_X); // number of predictors, including the thresholed variable
    int nThresholds=length(u_thresholdIdx);
    //PRINTF("thresholdIdx: "); for (int i=0; i<nThresholds; i++) PRINTF("%i ", thresholdIdx[i]); PRINTF("\n");
        
    // The rows and colns are organized in a way now that they can be directly casted and there is no need to do things as in the JSS paper on sycthe or MCMCpack MCMCmetrop1R.cc
    Matrix<double,Col,Concrete> Xcol (n, p, X_dat); //column major     
    Matrix<double,Row,Concrete> X(Xcol); // convert to row major so that creating bootstrap datasets can be faster and to pass to _grid_search
    Matrix<double,Row,Concrete> Y(n, 1, Y_dat); // define as a matrix instead of vector b/c will be used in matrix operation
    vector<double> W(W_dat, W_dat + n);
    
	
	int i,j;

	if (nBoot<0.1) {
    // a single search
    
    	vector<double> thresholds(nThresholds); 
        Matrix <double,Row,Concrete> Xcusum (n, p);
    	vector<double> Ycusum(n), Wcusum(n); 
    
        SEXP _logliks=PROTECT(allocVector(REALSXP, nThresholds));
        double *logliks=REAL(_logliks);       
        
        // compute Y'HY. this is not needed to find e_hat, but good to have for comparison with other estimation methods
        double yhy=0;
        if (p>1) {
        // if p is 0, we cannot do what is in this condition         
           Matrix<> Z=X(0,0,X.rows()-1,p-2);
           yhy = ((t(Y) * Z) * invpd(crossprod(Z)) * (t(Z) * Y))(0);
            
            // Step 2. Compute B and r
            _preprocess(X, Y);
            //PRINTF("B\n"); for (int i=0; i<n; i++) {for (int j=0; j<p; j++) PRINTF("%f ", X(i,j)); PRINTF("\n");}
            //PRINTF("Y\n"); for (int i=0; i<n; i++) PRINTF("%f ", Y(i)); PRINTF("\n"); 
        }
        
        for(i=0; i<nThresholds; i++) thresholds[i]=X(thresholdIdx[i]-1,p-1);

        _fastgrid2step_search(X, Y, W, wAllOne, thresholdIdx, skipping, thresholds, 
                             n, p, nThresholds,
                             Xcusum, Ycusum, Wcusum, 
                             logliks); 
        for(i=0; i<nThresholds; i++) logliks[i]=logliks[i]+yhy; 
        //PRINTF("logliks : \n");  for(i=0; i<nThresholds; i++) PRINTF("%f ", logliks[i]); PRINTF("\n");  
            
        UNPROTECT(1);
        return _logliks;

    } else {
    // bootstrap
    
      //  //output variables: logliks will not be returned to R, estimates from each bootstrap copy will be stored in coef and returned
        SEXP _coef=PROTECT(allocVector(REALSXP, nBoot*(p+1)));// p slopes, 1 threshold
        double *coef=REAL(_coef);    
        
    	if (nSub==0) nSub=n;
    	
    	// these variables are reused within each bootstrap replicate
        Matrix <double,Row,Concrete> Xcusum (nSub, p);
    	vector<double> Ycusum(nSub), Wcusum(nSub); 
    
    	vector<int> index(nSub), index2(n);
        Matrix <double,Row,Concrete> Xb(nSub,p), Yb(nSub,1);
    	vector<double> Wb(nSub); 
    	double e_hat;
    	int minThresholdIdx=thresholdIdx[0], maxThresholdIdx=thresholdIdx[nThresholds-1], nThresholdsb;
    	
        for (int b=0; b<nBoot; b++) {        
            // create bootstrap dataset, note that index is 1-based
            if (nSub==n) {
                SampleReplace(n, n, &(index[0]));
            } else {
                SampleNoReplace(nSub, n, &(index[0]), &(index2[0]));
            }
            //for (i=0; i<nSub; i++) { PRINTF("%i ", index[i]); } PRINTF("\n"); 
            
            // Step 1: sort
            sort (index.begin(), index.end());
            
            for (i=0; i<nSub; i++) { //note that index need to -1 to become 0-based
                Xb(i,_)=X(index[i]-1,_); 
                Yb(i)  =Y(index[i]-1); 
                Wb[i]  =W[index[i]-1]; 
            } 
            // debug use, can be used to compare with non-boot
            //for (i=0; i<n; i++) { Xb(i,_)=X(i,_); Yb(i)=Y(i); Wb[i]=W[i]; } 
            //for (i=0; i<n; i++) {for (j=0; j<p; j++)  PRINTF("%f ", Xb(i,j)); PRINTF("%f \n", Yb(i)); } 
            
            // recompute thresholdIdx, nThresholds, thresholds
            
//            // this only works when there is no m-out-of-n subsampling       
//            nThresholdsb=nThresholds;
//            vector<double> thresholdsb;
//            vector<int> thresholdIdxb;
//            for(i=0; i<nThresholdsb; i++) {
//                thresholdIdxb.push_back(thresholdIdx[i]); //thresholdIdx and thresholdIdxb are the same
//                thresholdsb.push_back(Xb(thresholdIdx[i]-1,p-1)); // make new thresholds vector that match the resampled dataset
//            }
            
            // this works there is subsampling
            vector<double> thresholdsb;
            vector<int> thresholdIdxb;
            for (i=0; i<nSub; i++) { 
                if (index[i]>=minThresholdIdx && index[i]<=maxThresholdIdx) {
                    thresholdIdxb.push_back (i+1); // threshold id is 1-based
                    thresholdsb.push_back   (Xb(i,p-1));// make new thresholds vector
                }
            }
            nThresholdsb=thresholdsb.size();
            
            //for(i=0; i<nThresholdsb; i++) PRINTF("%i ", thresholdIdxb[i]); PRINTF("\n"); 

            //PRINTF("%i \n", nThresholdsb);
            
            double * logliks = (double *) malloc((nThresholdsb) * sizeof(double));

            if (p>1) {
            // if p is 0, we cannot do what is in this condition
                // Compute B and r
                _preprocess(Xb, Yb);            
            }
            
            e_hat = _fastgrid2step_search(Xb, Yb, Wb, wAllOne, thresholdIdxb.data(), skipping, thresholdsb, 
                                 nSub, p, nThresholdsb,
                                 Xcusum, Ycusum, Wcusum, 
                                 logliks); 
            //PRINTF("e_hat %f\n", e_hat); 

            // fit model at the selected threshold and save results in coef
            
            // since after _preprocess, Xb and Yb have changed to B and r, we need to put X and Y back
            for (i=0; i<nSub; i++) { //note that index need to -1 to become 0-based
                Xb(i,_)=X(index[i]-1,_); 
                Yb(i)  =Y(index[i]-1); 
            }
            // debug use, can be used to compare with non-boot
            //for (i=0; i<n; i++) { Xb(i,_)=X(i,_); Yb(i)=Y(i); Wb[i]=W[i]; } 
            //PRINTF("Xb\n"); for (i=0; i<n; i++) {for (j=0; j<p; j++)  PRINTF("%f ", Xb(i,j)); PRINTF("\n");} 
            //PRINTF("Yb\n"); for (i=0; i<n; i++) PRINTF("%f ", Yb(i)); PRINTF("\n");
            // create x_e at e_hat
            // (x-e)+
            for(i=0; i<nSub; i++) 
                if(Xb(i,p-1)>e_hat) Xb(i,p-1) = 1; else Xb(i,p-1) = 0;  
            //PRINTF("Xb_e\n"); for (i=0; i<n; i++) {for (j=0; j<p; j++)  PRINTF("%f ", Xb(i,j)); PRINTF("\n");} 
            //PRINTF("Yb\n"); for (i=0; i<n; i++) PRINTF("%f ", Yb(i)); PRINTF("\n");
            
                        
            Matrix <> beta_hat = invpd(crossprod(Xb)) * (t(Xb) * Yb);
            for (j=0; j<p; j++) coef[b*(p+1)+j]=beta_hat[j];             
            
            coef[b*(p+1)+p] = e_hat;
            
            free(logliks);
        } 
         
        UNPROTECT(1);
        return _coef;
    }
    
}

}  // end extern C

#endif
//#endif
