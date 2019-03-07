// grid_search and boot_grid_search were removed because it is done as a comparator and did not implemented weights, thinned thresholds etc.


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


extern "C" {



         
  //vectors such as thresholds are defined outside this function since when bootstrapping, we don't want to allocate the memory over and over again
void _fastgrid_search(Matrix<double,Row>& X, vector<double>& Y, vector<double>& w, bool wAllOne, int * thresholdIdx,
    int n, int p, int nThresholds, bool isUpperHinge,
    Matrix<double,Row>& Xcusum, vector<double>& Ycusum, vector<double>& Wcusum, vector<double>& thresholds, vector<double>& Cps, 
    double* ans1, double* ans2)
{

    // loop index. iX is the index of the x vector, i is the index of threshold vector
    int i,iX,j;
    int k;
    
	int chosen=0;//initialize to get rid of compilation warning
    vector<double> C(p); 

    // for hinge or segmented, compute cusum of X and Y in the reverse order
    // for upperhinge, compute cusum of X and Y in the forward order
    if (!isUpperHinge) {
        Xcusum(n-1,_) = X(n-1,_)* w[n-1]; for (i=n-2; i>=0; i--) Xcusum(i,_) = Xcusum(i+1,_) + X(i,_)* w[i]; 
        Ycusum[n-1] = Y[n-1]* w[n-1];     for (i=n-2; i>=0; i--) Ycusum[i] = Ycusum[i+1] + Y[i]* w[i]; 
        Wcusum[n-1] = w[n-1];             for (i=n-2; i>=0; i--) Wcusum[i] = Wcusum[i+1] + w[i]; 
    } else {
        Xcusum(0,_) = X(0,_)* w[0]; for (i=1; i<n; i++) Xcusum(i,_) = Xcusum(i-1,_) + X(i,_)* w[i]; 
        Ycusum[0] = Y[0]* w[0];     for (i=1; i<n; i++) Ycusum[i] = Ycusum[i-1] + Y[i]* w[i]; 
        Wcusum[0] = w[0];           for (i=1; i<n; i++) Wcusum[i] = Wcusum[i-1] + w[i]; 
    }
    //for (i=0; i<n; i++) {for (j=0; j<p; j++)  PRINTF("%f ", Xcusum(i,j)); PRINTF("\n");}        
    //for (i=0; i<n; i++) PRINTF("%f ", Ycusum[i]); PRINTF("\n");
    
    // grid of e
    for(i=0; i<nThresholds; i++) thresholds[i]=X(thresholdIdx[i]-1,p-1);

    Matrix <double> A, Ainv, Ainv_save;
    double delta, rss, rss_max=R_NegInf;
    
    // Step 2: Initialize A:=X'X and C:=X'Y at the first point in the grid
    // compute X_e, which is the last column of X
    if (!isUpperHinge) {
        // (x-e)+
        for(iX=0; iX<thresholdIdx[0]; iX++) X(iX,p-1)=0;  
        for(iX=thresholdIdx[0]; iX<n; iX++) X(iX,p-1)=X(iX,p-1)-thresholds[0]; 
    } else {
        // (x-e)-
        for(iX=0; iX<thresholdIdx[0]; iX++) X(iX,p-1)=X(iX,p-1)-thresholds[0];  
        for(iX=thresholdIdx[0]; iX<n; iX++) X(iX,p-1)=0; 
    }
    //for (j=0; j<n; j++) PRINTF("%f ", X(j,p-1)); PRINTF("\n");
    
    if(wAllOne) {
        for (j=0; j<p; j++) {C[j]=0; for (int cntr=0; cntr<n; cntr++) C[j] += X(cntr,j) * Y[cntr];}
        A = crossprod1(X); 
    } else {
        for (j=0; j<p; j++) {
            C[j]=0; 
            for (int cntr=0; cntr<n; cntr++) 
                C[j] += X(cntr,j) * Y[cntr] * w[cntr];
        }
        // mutiply X with sqrt(w) before crossprod. Note that this changes X! Don't use X after this
        for (iX=0; iX<n; iX++) X(iX,_)=X(iX,_)*sqrt(w[iX]);
        A = crossprod1(X);
    }
    Cps[0]=C[p-1];// save C[p-1] to be used for computing parameter estimate when done
            
    // Step 4:
    for(i=0; i<nThresholds; i++) {            
           
        // Step 4a: update A:=X'X and C:=X'Y 
        // wAllOne handled by replacing n-i with Wcusum(i)
        if(i>0) {
            delta= thresholds[i]-thresholds[i-1]; //X(i,p-1)-X(i-1,p-1); //PRINTF("%f ", delta); PRINTF("\n"); 
            if (!isUpperHinge) {
                iX = thresholdIdx[i]-1; 
                // update A
                A(p-1,p-1) += pow(delta,2) * Wcusum[iX]; // + Delta' * W * Delta
                for (j=0; j<p-1; j++) A(p-1,j) -= delta * Xcusum(iX,j); // - Delta' * W * X
                for (j=0; j<p-1; j++) A(j,p-1) = A(p-1,j);  // - X' * W * Delta
                A(p-1,p-1) -= 2 * delta * (Xcusum(iX,p-1) - thresholds[i-1]*Wcusum[iX]); // the last element of both Delta' * W * X and X' * W * Delta
                // update C
                C[p-1] -= delta * Ycusum[iX];
            } else {
                k  = thresholdIdx[i-1]; // k is 1-based index of x, e_t = x_k
                // update A
                A(p-1,p-1) += pow(delta,2) * Wcusum[k-1]; // + Delta' * W * Delta
                for (j=0; j<p-1; j++) A(p-1,j) -= delta * Xcusum(k-1,j); // - Delta' * W * X
                for (j=0; j<p-1; j++) A(j,p-1) = A(p-1,j);  // - X' * W * Delta
                A(p-1,p-1) -= 2 * delta * (Xcusum(k-1,p-1) - thresholds[i-1]*Wcusum[k-1]); // the last element of both Delta' * W * X and X' * W * Delta
                // update C
                C[p-1] -= delta * Ycusum[k-1];                
            }
            Cps[i]=C[p-1];// save C[p-1]
        //for (int cntr=0; cntr<p; cntr++) for (j=0; j<p; j++)  PRINTF("%f ", A(cntr,j)); PRINTF("\n");        
        }
        
        // Step 4b (or step 3 when i=0): compute Y'HY                
        Ainv = invpd(A);
        rss=0; for (j=0; j<p; j++) for (int cntr=0; cntr<p; cntr++) rss += C[j] * C[cntr] * Ainv(j,cntr); // Y'HY
        ans1[i] = rss;
        if(rss>=rss_max) {
            chosen = i;
            Ainv_save=Ainv;
            rss_max=rss;
        }             
    }        

    // save results: estimated coefficients and threshold
    ans2[p]=thresholds[chosen];
    C[p-1]=Cps[chosen];
    for (int jj=0; jj<p; jj++) { ans2[jj]=0; for (j=0; j<p; j++) ans2[jj]+=Ainv_save(jj,j)*C[j]; }            
    
}

  // For fastgrid,
  // assume X and Y are sorted in chngptvar from small to large
  // assume last col of X is chngptvar, which will be updated as move through the grid
  // thresholdIdx are 1-based index, which define the grid of thresholds
  // For fastgrid2, the meaning of the first two variables are B and r instead 
SEXP fastgrid_gaussian(SEXP u_X, SEXP u_Y, SEXP u_W, SEXP u_wAllOne, SEXP u_thresholdIdx, SEXP u_nBoot, SEXP u_isUpperHinge)
{
    double* uX_dat = REAL(u_X);
    double* Y_dat=REAL(u_Y);
    double* W_dat=REAL(u_W);    
    bool wAllOne=asLogical(u_wAllOne)==1;
    int *thresholdIdx=INTEGER(u_thresholdIdx);
    int nBoot = asInteger(u_nBoot);
    bool isUpperHinge=asLogical(u_isUpperHinge)==1;
    //for (int i=0; i<nThresholds; i++) PRINTF("%i ", thresholdIdx[i]); PRINTF("\n");
    
    const int n = nrows(u_X);
    const int p = ncols(u_X);
    int nThresholds=length(u_thresholdIdx);
        
    // The rows and colns are organized in a way now that they can be directly casted and there is no need to do things as in the JSS paper on sycthe or MCMCpack MCMCmetrop1R.cc
    Matrix<double,Col,Concrete> Xcol (n, p, uX_dat); //column major     
    Matrix<double,Row,Concrete> X(Xcol); // convert to row major so that creating bootstrap datasets can be faster and to pass to _grid_search
    //for (i=0; i<n; i++) {for (j=0; j<p; j++)  PRINTF("%f ", Xb(i,j)); PRINTF("\n");}        
    
	// these variables are reused within each bootstrap replicate
	vector<int> index(n);
    Matrix <double,Row,Concrete> Xb (n, p), Xcusum (n, p, true, 0);
	vector<double> thresholds(nThresholds), Cps(nThresholds), Ycusum(n), Yb(n), Wb(n); 
	vector<double> Wcusum(n); //double rsses[nThresholds], Wcusum[n];
    	
	if (nBoot<0.1) {
    // a single search
        //output variables: logliks will be stored in ans and returned, stats won't be returned
        SEXP _ans=PROTECT(allocVector(REALSXP, nThresholds));
        double *ans=REAL(_ans);        
        double * stats = (double *) malloc((p+1) * sizeof(double));
        
        for (int i=0; i<n; i++) Yb[i]=Y_dat[i]; // cannot pass Y_dat directly b/c type mismatch
        for (int i=0; i<n; i++) Wb[i]=W_dat[i]; 
        _fastgrid_search(X, Yb, Wb, wAllOne, thresholdIdx, n, p, nThresholds, isUpperHinge,
                            Xcusum, Ycusum, Wcusum, thresholds, Cps, ans, stats); 
        //PRINTF("logliks : \n");  for(int i=0; i<nThresholds; i++) PRINTF("%f ", ans[i]); PRINTF("\n");  
            
        UNPROTECT(1);
        free(stats);
        return _ans;

    } else {
    // bootstrap
        //output variables: rsses won't be returned to R, estimates from each bootstrap copy will be stored in ans and returned
        double * rsses = (double *) malloc((nThresholds) * sizeof(double));
        SEXP _ans=PROTECT(allocVector(REALSXP, nBoot*(p+1)));// p slopes, 1 threshold, 1 goodness of fit stat
        double *ans=REAL(_ans);    
        
        for (int b=0; b<nBoot; b++) {        
            // create bootstrap dataset, note that index is 1-based
            SampleReplace(n, n, &(index[0]));
            // Step 1: sort
            sort (index.begin(), index.end());
            for (int i=0; i<n; i++) { //note that index need to -1 to become 0-based
                Xb(i,_)=X(index[i]-1,_); 
                Yb[i]=Y_dat[index[i]-1]; 
                Wb[i]=W_dat[index[i]-1]; 
            } 
            //PRINTF("Xb\n"); for (int i=0; i<n; i++) {for (int j=0; j<p; j++)  PRINTF("%f ", Xb(i,j)); PRINTF("\n");} 
            //PRINTF("Yb\n"); for (int i=0; i<n; i++) PRINTF("%f ", Yb[i]); PRINTF("\n");
            _fastgrid_search(Xb, Yb, Wb, wAllOne, thresholdIdx, n, p, nThresholds, isUpperHinge,
                                 Xcusum, Ycusum, Wcusum, thresholds, Cps, rsses, ans+b*(p+1)); // rsses not used here
        } 
         
        UNPROTECT(1);
        free(rsses);
        return _ans;
    }
    
}
  

    
}

#endif
//#endif

