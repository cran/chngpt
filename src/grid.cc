
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


using namespace std;
using namespace scythe;



extern "C" {

double _grid_search(Matrix<double,Row>& X, Matrix<double,Row>& Y, vector<double>& w, 
    int * thresholdIdx, vector<double>& thresholds, 
    int n, int p, int nThresholds,
    double* logliks)
{

    // loop index. 
    int i,j;    
    int k; // k has a meaning in the algorithm 
    int chosen=0;//initialize to get rid of compilation warning

    double delta, crit, crit_max=R_NegInf;
    
    vector<double> x(n), x_cpy(n); // x_cpy is needed for upperHinge b/c as we progress from left to right, the original vals of x's are lost
    for(i=0; i<n; i++) x[i] = X(i,p-1);
//    if (isUpperHinge) for(i=0; i<n; i++) x_cpy[i] = x[i];
    
    // multiply X and Y with sqrt(w), but first save x
    for(i=0; i<n; i++) {
        X(i,_)=X(i,_)*sqrt(w[i]);
        Y(i) *= sqrt(w[i]);
    }

    // compute X_e, which is the last column of X
    //PRINTF("X[,p-1]\n"); for (int ii=0; ii<n; ii++) PRINTF("%f ", X(ii,p-1)); PRINTF("\n");
//    if (!isUpperHinge) {
//        // (x-e)+
//        for(i=0; i<thresholdIdx[0]; i++) x[i]=0;  
//        for(i=thresholdIdx[0]; i<n; i++) x[i]=x[i]-thresholds[0]; 
//    } else {
        // (x-e)-
        for(i=0; i<thresholdIdx[0]; i++) x[i]=x[i]-thresholds[0];  
        for(i=thresholdIdx[0]; i<n; i++) x[i]=0; 
//    }
    //PRINTF("X[,p-1]\n"); for (int ii=0; ii<n; ii++) PRINTF("%f ", X(ii,p-1)); PRINTF("\n");
    
    //PRINTF("nThresholds %i\n", nThresholds); 
    for(i=0; i<nThresholds; i++) {             
        // update X_e
        if(i>0) {
            delta=thresholds[i]-thresholds[i-1]; 
            //PRINTF("delta %f\n", delta); 
//            if (!isUpperHinge) {
//                // (x-e)+
//                for (k=thresholdIdx[i-1]; k<thresholdIdx[i]-1; k++) x[k]=0; //needed for thinned threshold
//                for (k=thresholdIdx[i]-1; k<n; k++) x[k]=x[k]-delta;
//            } else {
                // (x-e)-
                for (k=0; k<thresholdIdx[i-1]; k++) x[k]=x[k]-delta;
                for (k=thresholdIdx[i-1]; k<thresholdIdx[i]-1; k++) x[k]=x_cpy[k]-thresholds[i]; // needed for thinned threshold
//            }
            //PRINTF("x\n"); for (int ii=0; ii<n; ii++) PRINTF("%f ", x[ii]); PRINTF("\n");
        }
        
        // replace the last column of X with weighted x
        for(j=0; j<n; j++) {
            X(j,p-1)=x[j]*sqrt(w[j]);
        }
        //PRINTF("X[,p-1]\n"); for (int ii=0; ii<n; ii++) PRINTF("%f ", X(ii,p-1)); PRINTF("\n");

        // compute criterion function         
        Matrix<> H = X * invpd(crossprod1(X)) * t(X); 
        crit = (t(Y) * H * Y)(0);
        logliks[i] = crit;
        if(crit>=crit_max) {
            chosen = i;
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
// For grid, the meaning of the first two variables are B and r instead 
SEXP gridC_gaussian(
     SEXP u_X, SEXP u_Y, 
     SEXP u_W,
     SEXP u_thresholdIdx, SEXP u_skipping,
     SEXP u_nBoot, SEXP u_nSub)
{
    double* X_dat = REAL(u_X);
    double* Y_dat=REAL(u_Y); 
    double* W_dat=REAL(u_W);    

    int *thresholdIdx=INTEGER(u_thresholdIdx);
    int nBoot = Rf_asInteger(u_nBoot);
//    bool isUpperHinge=asLogical(u_isUpperHinge)==1;
    
    int i,j;

    const int n = Rf_nrows(u_X);
    const int p = Rf_ncols(u_X); // number of predictors, including the thresholed variable
    int nThresholds=Rf_length(u_thresholdIdx);
    //for (i=0; i<nThresholds; i++) PRINTF("%i ", thresholdIdx[i]); PRINTF("\n");
        
    // The rows and colns are organized in a way now that they can be directly casted and there is no need to do things as in the JSS paper on sycthe or MCMCpack MCMCmetrop1R.cc
    Matrix<double,Col,Concrete> Xcol (n, p, X_dat); //column major     
    Matrix<double,Row,Concrete> X(Xcol); // convert to row major so that creating bootstrap datasets can be faster and to pass to _grid_search
    Matrix<double,Row,Concrete> Y(n, 1, Y_dat); // define as a matrix instead of vector b/c will be used in matrix operation
    vector<double> W(W_dat, W_dat + n);
    
    vector<double> thresholds(nThresholds); 
    
    if (nBoot<0.1) {
    // a single search
    
        SEXP _logliks=PROTECT(Rf_allocVector(REALSXP, nThresholds));
        double *logliks=REAL(_logliks);        
        
        for(i=0; i<nThresholds; i++) thresholds[i]=X(thresholdIdx[i]-1,p-1);
        //PRINTF("thresholds: "); for(i=0; i<nThresholds; i++) PRINTF("%f ", thresholds[i]); PRINTF("\n");

        _grid_search(X, Y, W, 
                             thresholdIdx, thresholds, 
                             n, p, nThresholds,
                             logliks); 
        //PRINTF("logliks : \n");  for(i=0; i<nThresholds; i++) PRINTF("%f ", logliks[i]); PRINTF("\n");  
            
        UNPROTECT(1);
        return _logliks;

    } else {
    // bootstrap
    
      //  //output variables: logliks will not be returned to R, estimates from each bootstrap copy will be stored in coef and returned
        double * logliks = (double *) malloc((nThresholds) * sizeof(double));
        SEXP _coef=PROTECT(Rf_allocVector(REALSXP, nBoot*(p+1)));// p slopes, 1 threshold
        double *coef=REAL(_coef);    
        
        // these variables are reused within each bootstrap replicate
        vector<int> index(n);
        Matrix <double,Row,Concrete> Xb(n,p), Yb(n,1);
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
//            PRINTF("index\n"); for (i=0; i<n; i++) PRINTF("%i ", index[i]); PRINTF("\n");
//            PRINTF("Xb\n"); for (i=0; i<n; i++) PRINTF("%f ", Xb(i,p-1)); PRINTF("\n");
//            PRINTF("Wb\n"); for (i=0; i<n; i++) PRINTF("%f ", Wb[i]); PRINTF("\n");
            
            for(i=0; i<nThresholds; i++) thresholds[i]=Xb(thresholdIdx[i]-1,p-1);
            
            e_hat = _grid_search(Xb, Yb, Wb, 
                                 thresholdIdx, thresholds, 
                                 n, p, nThresholds,
                                 logliks); 
            //PRINTF("e_hat %f\n", e_hat); 

            //////////////////////////////////////////////////
            // fit model at the selected threshold and save results in coef
            
            // after _grid_search, the last col of Xb changed
            for (i=0; i<n; i++) { //note that index need to -1 to become 0-based
                Xb(i,p-1)=X(index[i]-1,p-1); 
            }
            //PRINTF("Xb\n"); for (i=0; i<n; i++) {for (j=0; j<p; j++)  PRINTF("%f ", Xb(i,j)); PRINTF("\n");} 
            //PRINTF("Yb\n"); for (i=0; i<n; i++) PRINTF("%f ", Yb(i)); PRINTF("\n");
//            for (i=0; i<n; i++) { Xb(i,_)=X(i,_); Yb(i)=Y(i); } // debug use

            // create x_e at e_hat
//            if (!isUpperHinge) {
//                // (x-e)+
//                for(i=0; i<n; i++) 
//                    if(Xb(i,p-1)<e_hat) Xb(i,p-1) = 0; else Xb(i,p-1) -= e_hat;  
//            } else {
                // (x-e)-
                for(i=0; i<n; i++) 
                    if(Xb(i,p-1)<e_hat) Xb(i,p-1) -= e_hat; else Xb(i,p-1) = 0;  
//            }
            //PRINTF("Xb_e\n"); for (i=0; i<n; i++) {for (j=0; j<p; j++)  PRINTF("%f ", Xb(i,j)); PRINTF("\n");} 
            //PRINTF("Yb\n"); for (i=0; i<n; i++) PRINTF("%f ", Yb(i)); PRINTF("\n");
            
            // other columns of X and Y have already been weighted in the call _fast...
            for(i=0; i<n; i++) Xb(i,p-1)=Xb(i,p-1)*sqrt(Wb[i]);
                        
            Matrix <> beta_hat = invpd(crossprod(Xb)) * (t(Xb) * Yb);
            for (j=0; j<p; j++) coef[b*(p+1)+j]=beta_hat[j];             
            
            coef[b*(p+1)+p] = e_hat;
        } 
         
        UNPROTECT(1);
        free(logliks);
        return _coef;
    }
    
}

}  // end extern C

#endif
//#endif
