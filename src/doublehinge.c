
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


#include <float.h> //DBL_EPSILON
#include <R.h>           // needed to use Rprintf()
#include <R_ext/Utils.h> // needed to allow user interrupts
#include <R_ext/Lapack.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>

#define RUNIF runif
#define PRINTF Rprintf



/*
    mse = sapply(chngpts.1, function(e.1) {
          sapply(chngpts.2, function(e.2) {
            y.hat=ifelse(x<=e.1, lower.y, ifelse(x>=e.2, upper.y, lower.y+(upper.y-lower.y)/(e.2-e.1)*(x-e.1) ))
            sum((y-y.hat)**2)
    })
    })
*/



SEXP double_hinge_fit(
     SEXP u_X, SEXP u_Y, 
     SEXP u_chngpts_1, SEXP u_chngpts_2, 
     SEXP u_lower_y, SEXP u_upper_y,
     SEXP u_nBoot)
{
     
    double* X = REAL(u_X);
    double* Y=REAL(u_Y); 
    double* chngpts_1=REAL(u_chngpts_1);    
    double* chngpts_2=REAL(u_chngpts_2);    
    double lower_y = asReal(u_lower_y);
    double upper_y = asReal(u_upper_y);
    int nBoot = asInteger(u_nBoot);
    
    int n_1=length(u_chngpts_1);
    int n_2=length(u_chngpts_2);
    int n  =length(u_X);
    
    int i, j, k;
    double x, y_hat, mse, e_1, e_2;    
    
//	if (nBoot<0.1) {
//    // a single search
      
    // allocate memory for return variables and other variables
    SEXP _coef=PROTECT(allocVector(REALSXP, 3*(nBoot+1)));// p slopes, 1 threshold
    double *coef=REAL(_coef);    
    double * mses = (double *) malloc(n_1 * n_2 * sizeof(double));
    
    
    // fill mses
    int ind=0;
    int which_min=0; 
    double mse_min=INFINITY;
    for(i=0; i<n_1; i++) {
    for(j=0; j<n_2; j++) {
        e_1=chngpts_1[i];
        e_2=chngpts_2[j];
        mse=0;
        for(k=0; k<n; k++) {
            x=X[k];
            //work out y_hat
            y_hat=lower_y;
            if (x>e_1) {
                y_hat=upper_y;
                if (x<e_2) {
                    y_hat=lower_y+(upper_y-lower_y)/(e_2-e_1)*(x-e_1);
                }
            }
            mse+=pow(Y[k]-y_hat, 2);
        }
        if(mse<mse_min) {
            mse_min=mse;
            which_min=ind;
        }
        mses[ind++]=mse; 
    }
    }
    
    
    i=which_min / n_2;
    j=which_min % n_2;
    //PRINTF("%i %i %i ", which_min, i, j); PRINTF("\n");
    
    e_1=chngpts_1[i];
    e_2=chngpts_2[j];
    coef[0]=e_1; 
    coef[1]=e_2; 
    coef[2]=(upper_y-lower_y)/(e_2-e_1); 
    
    

/*
    	// these variables are reused within each bootstrap replicate
    	vector<int> index(n);
        Matrix <double,Row,Concrete> Xb(n,p), Yb(n,1), Xbreg(n,p+1); // Xbreg is used in regression, it contains the quadratic term
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
*/            
    

    UNPROTECT(1);
    free(mses);
    return _coef;
    
}


#endif
