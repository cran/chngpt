
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

  // assume X and Y are sorted in chngptvar from small to large
  // assume last col of X is chngptvar
  // thresholdIdx are 1-based index, which define the grid of thresholds
SEXP fastgrid2_gaussian(
     SEXP u_model,
     SEXP u_X, 
     SEXP u_Y, 
     SEXP u_W,
     SEXP u_thresholdIdx, 
     SEXP u_verbose,
     // bootstrap arguments
     SEXP u_nBoot, 
     SEXP u_nSub,
     SEXP u_withReplacement,
     SEXP u_sieveY
){

    int model = asInteger(u_model);
    int ncut = model==111?2:1; 
    
    double* X_dat = REAL(u_X);
    double* Y_dat=REAL(u_Y); 
    double* W_dat=REAL(u_W);  
    bool sieveBoot = !isNull(u_sieveY);
    double* sieveY_dat; 
    if (sieveBoot) sieveY_dat=REAL(u_sieveY); else sieveY_dat=Y_dat; // the latter is just to get rid of no-initialization warning
    
    int *thresholdIdx=INTEGER(u_thresholdIdx);
    int verbose=asInteger(u_verbose); 
    int nBoot = asInteger(u_nBoot);
    
    const int n = nrows(u_X);
    const int p = ncols(u_X); 
    int nThresholds=length(u_thresholdIdx);
    int nSub = asInteger(u_nSub);
    if (nSub==0) nSub=n;
    int withReplacement=asInteger(u_withReplacement); 
    //PRINTF("thresholdIdx: "); for (int i=0; i<nThresholds; i++) PRINTF("%i ", thresholdIdx[i]); PRINTF("\n");
        
    // The rows and colns are organized in a way now that they can be directly casted and there is no need to do things as in the JSS paper on sycthe or MCMCpack MCMCmetrop1R.cc
//PRINTF("DEBUG n %i p %i \n", n, p); 

// works okay under volta
//R version 3.5.3 (2019-03-11)
//Platform: x86_64-pc-linux-gnu (64-bit)
//Running under: Ubuntu 14.04.5 LTS
//
//Matrix products: default
//BLAS/LAPACK: /app/easybuild/software/OpenBLAS/0.2.18-GCC-5.4.0-2.26-LAPACK-3.6.1/lib/libopenblas_prescottp-r0.2.18.so

//but not under rhino
//R version 3.6.2 (2019-12-12)
//Platform: x86_64-pc-linux-gnu (64-bit)
//Running under: Ubuntu 18.04.4 LTS
//
//Matrix products: default
//BLAS/LAPACK: /app/software/OpenBLAS/0.3.7-GCC-8.3.0/lib/libopenblas_haswellp-r0.3.7.so

    Matrix<double,Col,Concrete> Xcol (n, p, X_dat); //column major     
//PRINTF("DEBUG 0\n"); 
    Matrix<double,Row,Concrete> X(Xcol); // convert to row major so that creating bootstrap datasets can be faster and to pass to _grid_search
    Matrix<double,Row,Concrete> Y(n, 1, Y_dat); // define as a matrix instead of vector b/c will be used in matrix operation    
    vector<double> w(W_dat, W_dat + n);

    vector<double> thresholds(nThresholds); 
    
    // importantly, both X and Y are weighted but not x because x will be compared to e
    vector<double> x(n);
    int i,j;
    for (i=0; i<n; i++) {
        w[i]=sqrt(w[i]);
        x[i] = X(i,p-1);    
        X(i,_)=X(i,_)*w[i];
        Y(i,0)=Y(i,0)*w[i];
    }
    
    Matrix<double, Row> Z=X(0,0,X.rows()-1,p-2); // in the R code we make sure p>1
    
    // define these here to avoid the need to re-allocate memory when doing bootstrap
    Matrix <double,Row,Concrete> Bcusum (n, p-1);
    vector<double> rcusum(n), Wcusum(n), xcusum(n); 
    // for higher order models we need to define more variables. 
    // to make it more efficient for the basic segmented models, we change n to nsmall
    int nsmall = model>10?n:1;
    Matrix <double,Row,Concrete> xBcusum (nsmall, p-1), x2Bcusum (nsmall, p-1), BcusumR (nsmall, p-1), xBcusumR (nsmall, p-1), x2BcusumR (nsmall, p-1);
    vector<double> x2cusum (nsmall), x3cusum (nsmall), xrcusum (nsmall), x4cusum (nsmall), x5cusum (nsmall), x2rcusum (nsmall), rcusumR(nsmall);
    vector<double> WcusumR(nsmall), xcusumR(nsmall), x2cusumR(nsmall), x3cusumR(nsmall), x4cusumR(nsmall), x5cusumR(nsmall), xrcusumR(nsmall), x2rcusumR(nsmall); 

    if (nBoot<0.1) {
    // a single search        
        if(model==111) {
            SEXP _ehatvec=PROTECT(allocVector(REALSXP, 2));
            double *ehatvec=REAL(_ehatvec);
            if (p>1) _preprocess(Z,Y);
            for(i=0; i<nThresholds; i++) thresholds[i]=x[thresholdIdx[i]-1];
            M111_search (Z, Y, x, w, 
                         thresholdIdx, thresholds, nThresholds, 
                         Bcusum, rcusum, Wcusum, xcusum,
                         x2cusum,
                         ehatvec, verbose); 
            UNPROTECT(1);
            return _ehatvec;        
            
        } else {
            SEXP _logliks=PROTECT(allocVector(REALSXP, nThresholds));
            double *logliks=REAL(_logliks);       
    
            // compute Y'HY. this is not needed to find e_hat, but good to have for comparison with other estimation methods
            double yhy = 0;
            
            // Step 2. Compute B and r, which are saved in Z and Y
            if (p>1) {
                yhy = ((t(Y) * Z) * invpd(crossprod(Z)) * (t(Z) * Y))(0);
                _preprocess(Z,Y);
            }
                
            for(i=0; i<nThresholds; i++) thresholds[i]=x[thresholdIdx[i]-1];
            if(model==5) { 
                   Mstep_search(Z, Y, x, w, 
                                 thresholdIdx, thresholds, nThresholds,
                                 Bcusum, rcusum, Wcusum,
                                 logliks);                              
            } else if(model==10) {
                   M10_search(Z, Y, x, w, 
                                 thresholdIdx, thresholds, nThresholds,
                                 Bcusum, rcusum, Wcusum, xcusum,
                                 logliks);                              
            } else if(model==20) {
                   M20_search (Z, Y, x, w, 
                                 thresholdIdx, thresholds, nThresholds, 
                                 Bcusum, rcusum, Wcusum, xcusum, x2cusum, x3cusum, xrcusum, xBcusum, 
                                 logliks); 
            } else if(model==204) {
                   M20c_search (Z, Y, x, w, 
                                 thresholdIdx, thresholds, nThresholds, 
                                 Bcusum, rcusum, Wcusum, xcusum, x2cusum, x3cusum, xrcusum, xBcusum, 
                                 logliks); 
            } else if(model==22) {
                   M22_search (Z, Y, x, w, 
                                 thresholdIdx, thresholds, nThresholds, 
                                 Bcusum, rcusum, Wcusum, xcusum, x2cusum, x3cusum, xrcusum, xBcusum, 
                                 BcusumR, rcusumR, WcusumR, xcusumR, x2cusumR, x3cusumR, xrcusumR, xBcusumR, 
                                 logliks); 
            } else if(model==224) {
                   M22c_search(Z, Y, x, w, 
                                 thresholdIdx, thresholds, nThresholds,
                                 Bcusum, rcusum, Wcusum, xcusum, x2cusum, x3cusum, xrcusum, xBcusum, 
                                 BcusumR, rcusumR, WcusumR, xcusumR, x2cusumR, x3cusumR, xrcusumR, xBcusumR, 
                                 logliks); 
            } else if(model==30) {
                   M30_search (Z, Y, x, w, 
                                 thresholdIdx, thresholds, nThresholds, 
                                 Bcusum, rcusum, Wcusum, xcusum, x2cusum, x3cusum, x4cusum, x5cusum, xrcusum, x2rcusum,  xBcusum, x2Bcusum, 
                                 logliks); 
            } else if(model==334) {
                   M33c_search (Z, Y, x, w, 
                                 thresholdIdx, thresholds, nThresholds, 
                                 Bcusum, rcusum, Wcusum, xcusum, x2cusum, x3cusum, x4cusum, x5cusum, xrcusum, x2rcusum,  xBcusum, x2Bcusum, 
                                 BcusumR, rcusumR, WcusumR, xcusumR, x2cusumR, x3cusumR, x4cusumR, x5cusumR, xrcusumR, x2rcusumR, xBcusumR, x2BcusumR, 
                                 logliks); 
            } else PRINTF("wrong \n");
            
            
            for(i=0; i<nThresholds; i++) logliks[i]=logliks[i]+yhy; 
            if(verbose>0) {PRINTF("search logliks: ");  for(i=0; i<nThresholds; i++) PRINTF("%f ", logliks[i]); PRINTF("\n");}            
            UNPROTECT(1);
            return _logliks;
        }

    } else if (nSub==n) {
    // Efron bootstrap

        //output variables: logliks will not be returned to R, estimates from each bootstrap copy will be stored in coef and returned
        // for getting results from M111
        double * ehatvec = (double *) malloc(ncut * sizeof(double));//stores the index of the largest likelihood
        // for getting results from other models
        double * logliks = (double *) malloc((nThresholds) * sizeof(double));
        
        // p is the number of other covariates + 1
        int p_coef=0;
        if        (model==5)  { p_coef=p+1; 
        } else if (model==10) { p_coef=p+1;
        } else if (model==20) { p_coef=p+2;
        } else if (model==204){ p_coef=p+1;
        } else if (model==22) { p_coef=p+4;
        } else if (model==224){ p_coef=p+3;
        } else if (model==30) { p_coef=p+3;
        } else if (model==334){ p_coef=p+4;
        } else if (model==111){ p_coef=p+3;
        } else PRINTF("wrong \n");;
        
        SEXP _coef=PROTECT(allocVector(REALSXP, nBoot*(p_coef)));
        double *coef=REAL(_coef);    

        Matrix <double,Row,Concrete> Zb(n,p-1), Xbreg(n,p_coef-ncut), Yb(n,1);
        vector<double> wb(n), xb(n);
        vector<int> index(n);         
        double e_hat=0;
        
        for (int b=0; b<nBoot; b++) {   
            
            if (!sieveBoot) {
                // create bootstrap dataset, note that index is 1-based
                SampleReplace(n, n, &(index[0]));
                // Step 1: sort
                sort (index.begin(), index.end());                                
                for (i=0; i<n; i++) { //note that index need to minus 1 to become 0-based
                    Zb(i,_)=Z(index[i]-1,_); 
                    Yb(i,0)=Y(index[i]-1,0); 
                    xb[i]  =x[index[i]-1];    
                    wb[i]  =w[index[i]-1]; 
                } 
            } else {   
               // use sieveY_dat for Y, keep x the same
                for (i=0; i<n; i++) { 
                    index[i]=i+1;// needed because index is used later as well
                    Zb(i,_)=Z(i,_); 
                    Yb(i,0)=sieveY_dat[b*n+i]; 
                    xb[i]  =x[i];    
                    wb[i]  =w[i]; 
                }             
            }
            //for (i=0; i<n; i++) {for (j=0; j<p-1; j++)  PRINTF("%f ", Zb(i,j)); PRINTF("%f %f %f ", Yb(i,0), xb[i], wb[i]); PRINTF("\n");} 
            //if(b==757 | b==0) {for (i=0; i<n; i++) PRINTF("%d,", index[i]); PRINTF("\n");}                
            
            for(i=0; i<nThresholds; i++) thresholds[i]=xb[thresholdIdx[i]-1];
            if (verbose>=2) { PRINTF("thresholds: "); for(i=0; i<nThresholds; i++) PRINTF("%f ", thresholds[i]); PRINTF("\n");}
            
            // compute Y'HY. This is not needed to find e_hat, but good to have for comparison with other estimation methods
            double yhyb=0;
            if(verbose>0 && p>1) yhyb = ((t(Yb) * Zb) * invpd(crossprod(Zb)) * (t(Zb) * Yb))(0);
    
            // Compute B and r, which are saved in Z and Y
            if (p>1) _preprocess(Zb,Yb);

            if(model==5) { 
                   e_hat=Mstep_search(Zb, Yb, xb, wb, 
                                 thresholdIdx, thresholds, nThresholds,
                                 Bcusum, rcusum, Wcusum,
                                 logliks);                              
            } else if(model==10) {
                   e_hat=M10_search(Zb, Yb, xb, wb, 
                                 thresholdIdx, thresholds, nThresholds,
                                 Bcusum, rcusum, Wcusum, xcusum,
                                 logliks);                              
            } else if(model==20) {
                   e_hat=M20_search(Zb, Yb, xb, wb, 
                                 thresholdIdx, thresholds, nThresholds, 
                                 Bcusum, rcusum, Wcusum, xcusum, x2cusum, x3cusum, xrcusum, xBcusum, 
                                 logliks); 
            } else if(model==204) {
                   e_hat=M20c_search(Zb, Yb, xb, wb, 
                                 thresholdIdx, thresholds, nThresholds, 
                                 Bcusum, rcusum, Wcusum, xcusum, x2cusum, x3cusum, xrcusum, xBcusum, 
                                 logliks); 
            } else if(model==22) {
                   e_hat=M22_search(Zb, Yb, xb, wb, 
                                 thresholdIdx, thresholds, nThresholds, 
                                 Bcusum, rcusum, Wcusum, xcusum, x2cusum, x3cusum, xrcusum, xBcusum, 
                                 BcusumR, rcusumR, WcusumR, xcusumR, x2cusumR, x3cusumR, xrcusumR, xBcusumR, 
                                 logliks); 
            } else if(model==224) {
                   e_hat=M22c_search(Zb, Yb, xb, wb, 
                                 thresholdIdx, thresholds, nThresholds,
                                 Bcusum, rcusum, Wcusum, xcusum, x2cusum, x3cusum, xrcusum, xBcusum, 
                                 BcusumR, rcusumR, WcusumR, xcusumR, x2cusumR, x3cusumR, xrcusumR, xBcusumR, 
                                 logliks); 
            } else if(model==30) {
                   e_hat=M30_search(Zb, Yb, xb, wb, 
                                 thresholdIdx, thresholds, nThresholds, 
                                 Bcusum, rcusum, Wcusum, xcusum, x2cusum, x3cusum, x4cusum, x5cusum, xrcusum, x2rcusum,  xBcusum, x2Bcusum, 
                                 logliks); 
            } else if(model==334) {
                   e_hat=M33c_search(Zb, Yb, xb, wb, 
                                 thresholdIdx, thresholds, nThresholds, 
                                 Bcusum, rcusum, Wcusum, xcusum, x2cusum, x3cusum, x4cusum, x5cusum, xrcusum, x2rcusum,  xBcusum, x2Bcusum, 
                                 BcusumR, rcusumR, WcusumR, xcusumR, x2cusumR, x3cusumR, x4cusumR, x5cusumR, xrcusumR, x2rcusumR, xBcusumR, x2BcusumR, 
                                 logliks); 
            } else if(model==111) {
                   M111_search (Zb, Yb, xb, wb, 
                                 thresholdIdx, thresholds, nThresholds, 
                                 Bcusum, rcusum, Wcusum, xcusum,
                                 x2cusum,
                                 ehatvec, verbose); 
            } else PRINTF("wrong \n");
            
            if(verbose>0) {
                if(model==111) {  
                    PRINTF("est thresholds: "); for(i=0; i<2; i++) PRINTF("%f ", ehatvec[i]); PRINTF("\n");
                } else {
                    for(i=0; i<nThresholds; i++) logliks[i]=logliks[i]+yhyb; 
                    PRINTF("boot logliks: "); for(i=0; i<nThresholds; i++) PRINTF("%f ", logliks[i]); PRINTF("\n");
                }
            }
                    
            // fit model at the selected threshold and save results in coef 
            // since Xb and Yb are changed during search, we need to copy from X 
            for (i=0; i<n; i++) { 
                
                if (sieveBoot) {
                    Yb(i,0)=sieveY_dat[b*n+i]; 
                } else {
                    Yb(i,0)=Y(index[i]-1,0); 
                }                
                for (j=0; j<p-1; j++) Xbreg(i,j)=X(index[i]-1,j);
                
                // thresholded covariates
                if (model==5) {     // step model               
                    if(x[index[i]-1]>e_hat) Xbreg(i,p-1) =1; else Xbreg(i,p-1) = 0; 
                
                } else if(model==10) {                    
                    if(x[index[i]-1]<e_hat) Xbreg(i,p-1) =x[index[i]-1]- e_hat; else Xbreg(i,p-1) = 0; // (x-e)- 
                
                } else if (model==20) {                    
                    if(x[index[i]-1]<e_hat) Xbreg(i,p-1) =x[index[i]-1]- e_hat; else Xbreg(i,p-1) = 0;// (x-e)-                    
                    Xbreg(i,p)=pow(Xbreg(i,p-1),2);// (x-e)-^2
                
                } else if (model==204) {                    
                    if(x[index[i]-1]<e_hat) Xbreg(i,p-2) =x[index[i]-1]- e_hat; else Xbreg(i,p-2) = 0;// (x-e)-                    
                    Xbreg(i,p-1)=pow(Xbreg(i,p-2),2);// (x-e)-^2
                    Xbreg(i,p-2) =x[index[i]-1]; // x
                
                } else if (model==22) {                    
                    if(x[index[i]-1]<e_hat) Xbreg(i,p-1) =x[index[i]-1]- e_hat; else Xbreg(i,p-1) = 0;// (x-e)-                    
                    Xbreg(i,p)=pow(Xbreg(i,p-1),2);// (x-e)-^2
                    // populate with weighted x
                    Xbreg(i,p+1)=X(index[i]-1,p-1);                    
                    if(x[index[i]-1]>e_hat) Xbreg(i,p+1) =x[index[i]-1]- e_hat; else Xbreg(i,p+1) = 0;// (x-e)+                    
                    Xbreg(i,p+2)=pow(Xbreg(i,p+1),2);// (x-e)+^2
                
                } else if (model==224) {                    
                    Xbreg(i,p-1) =x[index[i]-1]- e_hat; // x-e                    
                    if(x[index[i]-1]<e_hat) Xbreg(i,p)  =pow(Xbreg(i,p-1),2); else Xbreg(i,p)=0;// (x-e)-^2                    
                    if(x[index[i]-1]>e_hat) Xbreg(i,p+1)=pow(Xbreg(i,p-1),2); else Xbreg(i,p+1)=0;// (x-e)+^2                       
                    Xbreg(i,p-1) =x[index[i]-1]; // x
                       
                } else if (model==30) {                       
                    if(x[index[i]-1]<e_hat) Xbreg(i,p-1) =x[index[i]-1]- e_hat; else Xbreg(i,p-1) = 0;//(x-e)-
                    Xbreg(i,p)=pow(Xbreg(i,p-1),2);//(x-e)-^2
                    Xbreg(i,p+1)=pow(Xbreg(i,p-1),3);//(x-e)-^3
    
                } else if (model==334) {                    
                    Xbreg(i,p-1) =x[index[i]-1]- e_hat; // x-e
                    Xbreg(i,p)=pow(Xbreg(i,p-1),2);// (x-e)^2                    
                    if(x[index[i]-1]<e_hat) Xbreg(i,p+1)=pow(Xbreg(i,p-1),3); else Xbreg(i,p+1)=0;// (x-e)-^3                    
                    if(x[index[i]-1]>e_hat) Xbreg(i,p+2)=pow(Xbreg(i,p-1),3); else Xbreg(i,p+2)=0;// (x-e)+^3
                    Xbreg(i,p-1) =x[index[i]-1]; // x
                    
                } else if(model==111) {                    
                    if(x[index[i]-1]<ehatvec[0]) Xbreg(i,p-1) =x[index[i]-1]- ehatvec[0]; else Xbreg(i,p-1) = 0; // (x-e)- 
                    if(x[index[i]-1]<ehatvec[1]) Xbreg(i,p)   =x[index[i]-1]- ehatvec[1]; else Xbreg(i,p) = 0; // (x-f)- 
                
                } else PRINTF("wrong \n");               
                 
                // add weight to thresholded covariates
                for (j=p-1; j<p_coef-ncut; j++) Xbreg(i,j)=Xbreg(i,j)*wb[i];
            }            
            
            Matrix <> beta_hat = invpd(crossprod(Xbreg)) * (t(Xbreg) * Yb);
            if(model==111) {
                for (int j=0; j<p_coef-ncut; j++) coef[b*p_coef+j]=beta_hat[j];                         
                coef[(b+1)*p_coef-2] = *ehatvec;
                coef[(b+1)*p_coef-1] = *(ehatvec+1);
            } else {
                for (int j=0; j<p_coef-1; j++) coef[b*p_coef+j]=beta_hat[j];                         
                coef[(b+1)*p_coef-1] = e_hat;
            }
            
//            PRINTF("Xbreg\n"); for (i=0; i<n; i++) {for (j=0; j<p+ncut-1; j++)  PRINTF("%f ", Xbreg(i,j)); PRINTF("\n");} 
//            PRINTF("Yb\n"); for (i=0; i<n; i++) PRINTF("%f ", Yb(i)); PRINTF("\n");
//            for (int j=0; j<p_coef-ncut; j++) PRINTF("%f ", beta_hat[j]); PRINTF("\n");
//            //for (i=0; i<n; i++) { Xb(i,_)=X(i,_); Yb(i)=Y(i); } // debug use
        }
        
        UNPROTECT(1);
        free(logliks);
        free(ehatvec);
        return _coef;
        
    } else {
        // subsampling bootstrap
              
        int p_coef=0;
        if(model==5) { 
            p_coef=p+1; 
        } else PRINTF("wrong model in m-out-of-n bootstrap\n");;
                
        SEXP _coef=PROTECT(allocVector(REALSXP, nBoot*(p_coef)));
        double *coef=REAL(_coef);    

        // redefine these variables with nSub rows
        Matrix <double,Row,Concrete> Bcusum (nSub, p-1);
        vector<double> rcusum(nSub), Wcusum(nSub); 
    
        Matrix <double,Row,Concrete> Zb(nSub,p-1), Xbreg(nSub,p_coef-1), Yb(nSub,1);
        vector<double> wb(nSub), xb(nSub); 
        vector<int> index(nSub), index2(n);
        double e_hat=0;
        
        int nThresholdsb, minThresholdIdx=thresholdIdx[0], maxThresholdIdx=thresholdIdx[nThresholds-1];

        for (int b=0; b<nBoot; b++) {        
            if(withReplacement>0) SampleReplace (nSub, n, &(index[0])); else SampleNoReplace (nSub, n, &(index[0]), &(index2[0]));
            // Step 1: sort
            sort (index.begin(), index.end());
            
            for (i=0; i<nSub; i++) { //note that index need to -1 to become 0-based
                Zb(i,_)=Z(index[i]-1,_); 
                Yb(i,0)=Y(index[i]-1,0); 
                xb[i]  =x[index[i]-1];    
                wb[i]  =w[index[i]-1]; 
            } 
            //for (i=0; i<n; i++) {for (j=0; j<p-1; j++)  PRINTF("%f ", Zb(i,j)); PRINTF("%f %f %f ", Yb(i,0), xb[i], wb[i]); PRINTF("\n");} 
            
//            // this only works when there is no m-out-of-n subsampling       
//            nThresholdsb=nThresholds;
//            vector<double> thresholdsb;
//            vector<int> thresholdIdxb;
//            for(i=0; i<nThresholdsb; i++) {
//                thresholdIdxb.push_back(thresholdIdx[i]); //thresholdIdx and thresholdIdxb are the same
//                thresholdsb.push_back(Xb(thresholdIdx[i]-1,p-1)); // make new thresholds vector that match the resampled dataset
//            }
            
            // needs to do this for subsampling
            vector<double> thresholdsb;
            vector<int> thresholdIdxb;
            for (i=0; i<nSub; i++) { 
                if (index[i]>=minThresholdIdx && index[i]<=maxThresholdIdx) {
                    thresholdIdxb.push_back (i+1); // threshold id is 1-based
                    thresholdsb.push_back   (xb[i]);// make new thresholds vector
                }
            }
            nThresholdsb=thresholdsb.size();
            
            double * logliks = (double *) malloc((nThresholdsb) * sizeof(double));

            // Compute B and r, which are saved in Z and Y
            if (p>1) _preprocess(Zb,Yb);

            if(model==5) { 
                   e_hat=Mstep_search(Zb, Yb, xb, wb, 
                                 thresholdIdxb.data(), thresholdsb, nThresholdsb,
                                 Bcusum, rcusum, Wcusum,
                                 logliks);                              
            } else PRINTF("wrong \n");
                
            // fit model at the selected threshold and save results in coef 
            // since Xb and Yb are changed during search, we need to copy from X 
            for (i=0; i<nSub; i++) { 
                
                Yb(i,0)=Y(index[i]-1,0); 
                for (j=0; j<p-1; j++) Xbreg(i,j)=X(index[i]-1,j);
                
                // thresholded covariates
                if (model==5) {     // step model               
                    if(x[index[i]-1]>e_hat) Xbreg(i,p-1) =1; else Xbreg(i,p-1) = 0; // (x-e)-                 
                } else PRINTF("wrong \n");               
                 
                // add weight to thresholded covariates
                for (j=p-1; j<p_coef-1; j++) Xbreg(i,j)=Xbreg(i,j)*wb[i];
            }
            
            
            Matrix <> beta_hat = invpd(crossprod(Xbreg)) * (t(Xbreg) * Yb);
            for (int j=0; j<p_coef-1; j++) coef[b*p_coef+j]=beta_hat[j];                         
            coef[(b+1)*p_coef-1] = e_hat;
            
            free(logliks);
        
            //PRINTF("Xbreg\n"); for (i=0; i<n; i++) {for (j=0; j<p; j++)  PRINTF("%f ", Xbreg(i,j)); PRINTF("\n");} 
            //PRINTF("Yb\n"); for (i=0; i<n; i++) PRINTF("%f ", Yb(i)); PRINTF("\n");
            //for (int j=0; j<p_coef-1; j++) PRINTF("%f ", beta_hat[j]); PRINTF("\n");
            //for (i=0; i<n; i++) { Xb(i,_)=X(i,_); Yb(i)=Y(i); } // debug use
        }
        
        UNPROTECT(1);
        return _coef;
        
    }
    
} // end function


}  // end extern C

#endif


