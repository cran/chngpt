//////////////////////////////////////////////////////////////////////////
// 
// This software is distributed under the terms of the GNU GENERAL
// PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
// file for more information.
//
//////////////////////////////////////////////////////////////////////////

#ifndef CHNGPTBOOT_CC
#define CHNGPTBOOT_CC


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



extern "C" {
  
  // assume X and Y are sorted
  // nLower and nUpper are 1-based index
  SEXP boot_grid_search_fast(SEXP _X, SEXP _Y, SEXP _nLower, SEXP _nUpper, SEXP _B)
  {
    int i,j; //for loop index
    // put _X and _Y into Matrixes Xori and Y
    // note that the rows and colns are organized in a way now that they can be directly casted and there is no need to do things as in MCMCpack MCMCmetrop1R.cc
    double* X_dat = REAL(_X);
    const int n = nrows(_X);
    const int p = ncols(_X);
    Matrix<double,Col,Concrete> Xcol (n, p, X_dat); //column major 
    // convert to row major so that creating bootstrap datasets can be faster   
    Matrix<double,Row,Concrete> Xori(Xcol); // row major, name this Xori instead of X so that there is less a mistake of using X when Xb should be used below
    //for (i=0; i<n; i++) {for (j=0; j<p; j++)  PRINTF("%f ", Xb(i,j)); PRINTF("\n");}        
    double *Y=REAL(_Y);
    // bootstrap replicate
    double B = asReal(_B);
    // bounds
    int nLower=asInteger (_nLower);
    int nUpper=asInteger (_nUpper);
    
    // output
    SEXP _ans=PROTECT(allocVector(REALSXP, B*(p+2)));// p slopes, 1 threshold, 1 goodness of fit stat
    double *ans=REAL(_ans);    
    
	// these variables are reused within each bootstrap replicate
	vector<int> index(n);
    Matrix <double,Row,Concrete> Xb (n, p), Xcusum (n, p, true, 0);
    Matrix <double> H, J, A, Ainv, Ainv_save;
	vector<double> thresholds(nUpper-nLower+1), zps(nUpper-nLower+1); 
	vector<double> z(p), Ycusum(n), Yb(n); 
	int chosen=0;//initialize to get rid of compilation warning
	double efinal;
	Matrix <> estimatedSlope(p, 1);
	
	// loop through bootstrap replicates
    for (int b=0; b<B; b++) {
        
        // create bootstrap dataset
        SampleReplace(n, n, &(index[0]));// note that index is 1-based
        sort (index.begin(), index.end());
        for (i=0; i<n; i++) { Xb(i,_)=Xori(index[i]-1,_); Yb[i]=Y[index[i]-1]; } //-1 to become 0-based
        //for (i=0; i<n; i++) {for (j=0; j<p; j++)  PRINTF("%f ", Xb(i,j)); PRINTF("\n");}        
        
        // save a copy of thresholds because Xb gets updated as as loop through candidate thresholds
        // this is only a subset of X
        for(i=nLower-1; i<nUpper; i++) thresholds[i-nLower+1]=Xb(i,p-1);
        
        // compute cusum of Xb and Y in the reverse order
        Xcusum(n-1,_) = Xb(n-1,_);
        for (i=n-2; i>=0; i--) Xcusum(i,_) = Xcusum(i+1,_) + Xb(i,_); 
        //for (i=0; i<n; i++) {for (j=0; j<p; j++)  PRINTF("%f ", Xcusum(i,j)); PRINTF("\n");}        
        Ycusum[n-1] = Yb[n-1];
        for (i=n-2; i>=0; i--) Ycusum[i] = Ycusum[i+1] + Yb[i]; 
        //for (i=0; i<n; i++) PRINTF("%f ", Ycusum[i]); PRINTF("\n");
        
        // loop through candidate thresholds, Xb gets updated in the process
        double delta, rss, rss_min=INFINITY;
        // first set the x lower than nLower to 0
        for(i=0; i<nLower-1; i++) Xb(i,p-1)=0;                
        for(i=nLower-1; i<nUpper; i++) {    
            // Update Xb, (p-1)th column only
            delta=Xb(i,p-1); //delta is e in the first iteration and delta of success in the latter iterations
            for (j=i; j<n; j++) Xb(j,p-1)=Xb(j,p-1)-delta; //for (j=0; j<n; j++) PRINTF("%f ", Xb(j,p-1)); PRINTF("\n");
            //PRINTF("%f ", delta); PRINTF("\n");    
            
            // update A
            if(i==nLower-1) {
                A = crossprod(Xb); 
            } else {
                A(p-1,p-1) += (n-i)*pow(delta,2);
                for (j=0; j<p-1; j++) A(p-1,j) -= delta * Xcusum(i,j);
                for (j=0; j<p-1; j++) A(j,p-1) = A(p-1,j);  
                A(p-1,p-1) -= 2 * delta * (Xcusum(i,p-1) - (n-i)*thresholds[i-nLower]);
                //for (int k=0; k<p; k++) for (j=0; j<p; j++)  PRINTF("%f ", A(k,j)); PRINTF("\n");        
            }
            
//            // compute rss, this code won't run anymore because we have changed the type of Yb from Matrix to vector
//            H = Xb * invpd(A) * t(Xb); 
//            //the next line offer a faster way than: -(t(Yb) * H * Yb)(0);
//            rss=0; for (j=0; j<n; j++) for (int k=0; k<n; k++) rss -= Yb[j] * Yb[k] * H(j,k);
//            //PRINTF("%f ", rss); PRINTF("\n"); 
            
            // a faster way to compute rss
            // update z
            if(i==nLower-1) {
                //z = t(Xb) * Yb; 
                for (j=0; j<p; j++) {
                    z[j]=0;
                    for (int k=0; k<n; k++) z[j] += Xb(k,j) * Yb[k];
                }
            } else {
                z[p-1] -= delta * Ycusum[i];
            }
            // save z[p-1]
            zps[i-nLower+1]=z[p-1];
            Ainv = invpd(A);
            rss=0; for (j=0; j<p; j++) for (int k=0; k<p; k++) rss -= z[j] * z[k] * Ainv(j,k);
            //PRINTF("%f ", rss); PRINTF("\n");
            
            if(rss<=rss_min) {
                chosen = i-(nLower-1);
                Ainv_save=Ainv;
                rss_min=rss;
            }             
        }        
        
        // save the estimated coefficients and threshold
        efinal=thresholds[chosen];
        z[p-1]=zps[chosen];
        for (i=0; i<p; i++) {
            ans[b*(p+2)+i]=0;
            for (j=0; j<p; j++) ans[b*(p+2)+i]+=Ainv_save(i,j)*z[j];
        }            
        ans[b*(p+2)+p]=efinal;
        ans[b*(p+2)+p+1]=-rss_min;
        
    } // end bootstrap loop
     
    UNPROTECT(1);
    return _ans;
  }

  
  
// an optimized C/C++ implementation that serves as a reference for more advanced algebraic optimization
  // assume X and Y are sorted
  // nLower and nUpper are 1-based index
  SEXP boot_grid_search(SEXP _X, SEXP _Y, SEXP _nLower, SEXP _nUpper, SEXP _B)
  {
    int i,j; //for loop index
    // put _X and _Y into Matrixes Xori and Y
    // note that the rows and colns are organized in a way now that they can be directly casted and there is no need to do things as in MCMCpack MCMCmetrop1R.cc
    double* X_dat = REAL(_X);
    const int n = nrows(_X);
    const int p = ncols(_X);
    Matrix<double,Col,Concrete> Xcol (n, p, X_dat); //column major 
    // convert to row major so that creating bootstrap datasets can be faster   
    Matrix<double,Row,Concrete> Xori(Xcol); // row major        
    //for (i=0; i<n; i++) {for (j=0; j<p; j++)  PRINTF("%f ", Xb(i,j)); PRINTF("\n");}        
    double *Y_dat=REAL(_Y);
    Matrix <> Y (n, 1, Y_dat);
    // bootstrap replicate
    double B = asReal(_B);
    // bounds
    int nLower=asInteger (_nLower);
    int nUpper=asInteger (_nUpper);
    
    // output
    SEXP _ans=PROTECT(allocVector(REALSXP, B*(p+2)));// p slopes, 1 threshold, 1 goodness of fit stat
    double *ans=REAL(_ans);    
    
	// these variables are reused within each bootstrap replicate
	vector<int> index(n);
    Matrix <double,Row,Concrete> Xb (n, p);
    Matrix <> Yb (n, 1);
    Matrix <double> H, J;
	vector<double> logliks(nUpper-nLower+1); 
	vector<double> thresholds(nUpper-nLower+1); 
	int chosen;
	double efinal;
	Matrix <> estimatedSlope(p, 1);
	
	// loop through bootstrap replicates
    for (int b=0; b<B; b++) {
        
        // fill index, note that it is 1-based
        SampleReplace(n, n, &(index[0]));
        // sort index
        sort (index.begin(), index.end());
        // create bootstrap dataset
        for (i=0; i<n; i++) { 
            Xb(i,_)=Xori(index[i]-1,_); //-1 to become 0-based
            Yb(i)=Y(index[i]-1); 
        }
        
        // save a copy of x as thresholds because they get changed below
        for(i=nLower-1; i<nUpper; i++) thresholds[i-nLower+1]=Xb(i,p-1);
        // set the x lower than nLower to 0
        for(i=0; i<nLower-1; i++) Xb(i,p-1)=0;        
                
        // loop through candidate thresholds
        double delta;
        double rss;
        for(i=nLower-1; i<nUpper; i++) {    
            // Update the change point variable in Xb. 
            delta=Xb(i,p-1); //delta is e in the first iteration
            for (j=i; j<n; j++) Xb(j,p-1)=Xb(j,p-1)-delta; //for (j=0; j<n; j++) PRINTF("%f ", Xb(j,p-1)); PRINTF("\n");
            H = Xb * invpd(crossprod(Xb)) * t(Xb); // hat matrix
            //the next line offer a faster way than: -(t(Yb) * H * Yb)(0);
            rss=0; for (j=0; j<n; j++) for (int k=0; k<n; k++) rss -= Yb[j] * Yb[k] * H(j,k);
            logliks[i-(nLower-1)] = -rss; // since Yb'Yb does not depend on threshold, there is no need to compute it
        }        
        //for(i=nLower-1; i<nUpper; i++)  PRINTF("%f ", logliks[i-(nLower-1)]); PRINTF("\n");
        
        // save the estimated coefficients and threshold
        chosen = distance(logliks.begin(), max_element(logliks.begin(), logliks.end()));
        efinal=thresholds[chosen];
        // compute slope estimate. If slope estimates are computed withint the previous loop, it is slower
        for(i=nLower-1+chosen; i<n; i++) Xb(i,p-1)=Xori(index[i]-1,p-1)-efinal; //-1 to become 0-based
        //for (i=0; i<n; i++) {for (j=0; j<p; j++)  PRINTF("%f ", Xb(i,j)); PRINTF("\n");}        
        estimatedSlope = invpd(crossprod(Xb)) * t(Xb) * Yb;
        for (i=0; i<p; i++) ans[b*(p+2)+i]=estimatedSlope(i);
        ans[b*(p+2)+p]=efinal;
        ans[b*(p+2)+p+1]=logliks[chosen];
        
    }
     
    UNPROTECT(1);
    return _ans;
  }

  

  
  // the function that performs grid search and returns the best chngpt 
  // assume X is sorted in chngptvar from small to large
  // nLower and nUpper are 1-based index
  SEXP grid_search(SEXP _X, SEXP _Y, SEXP _nLower, SEXP _nUpper)
  {

    // put _X and _Y into Matrixes X and Y
    // note that the rows and colns are organized in a way now that they can be directly casted and there is no need to do things as in the JSS paper on sycthe
    int i,j;
    double* X_dat = REAL(_X);
    const int n = nrows(_X);
    const int p = ncols(_X);
    Matrix <> X (n, p, X_dat);
    double *Y_dat=REAL(_Y);
    Matrix <> Y (n, 1, Y_dat);
    int nLower=asInteger (_nLower);
    int nUpper=asInteger (_nUpper);
    // output
    SEXP _ans=PROTECT(allocVector(REALSXP, nUpper-nLower+1));
    double *ans=REAL(_ans);    
    
    // set the x lower than nLower to 0
    for(i=0; i<nLower-1; i++) X(i,p-1)=0;
    
    Matrix<double> H;
    double delta;
    for(i=nLower-1; i<nUpper; i++) {    
        // update the change point variable in X
        delta=X(i,p-1); // delta is e in the first iteration 
        for (j=i; j<n; j++) X(j,p-1)=X(j,p-1)-delta;
        //for (j=0; j<n; j++) PRINTF("%f ", X(j,p-1)); PRINTF("\n");
        H = - X * invpd(crossprod(X)) * t(X); // - hat matrix
        for (j=0; j<n; j++) H(j,j)=1+H(j,j); // I - H
        ans[i-(nLower-1)] = (t(Y) * H * Y)(0);
    }
          
    UNPROTECT(1);
    return _ans;
  }
  
  
  // unit testing for performance comparison
  SEXP performance_unit_test(SEXP _X, SEXP _Y, SEXP _B, SEXP _I)
  {

    int i,j; //for loop index
    // put _X and _Y into Matrixes X and Y
    // note that the rows and colns are organized in a way now that they can be directly casted and there is no need to do things as in MCMCpack MCMCmetrop1R.cc
    double* X_dat = REAL(_X);
    const int n = nrows(_X);
    const int p = ncols(_X);
    Matrix<double,Col,Concrete> Xcol (n, p, X_dat); //column major 
    // convert to row major so that creating bootstrap datasets can be faster   
    Matrix<double,Row,Concrete> X(Xcol); // row major        
    //for (i=0; i<n; i++) {for (j=0; j<p; j++)  PRINTF("%f ", Xb(i,j)); PRINTF("\n");}        
    double *Y_dat=REAL(_Y);
    Matrix <> Y (n, 1, Y_dat);
    // bootstrap replicate
    int B = asInteger(_B);
    int I = asInteger(_I);

    // output
    SEXP _ans=PROTECT(allocVector(REALSXP, 1));
    double *ans=REAL(_ans);    
    ans[0]=0;
    
//    // compare X and Xcol in computing H
//    Matrix<double> H;
//    if(I==1) for(i=1; i<B; i++) H = - X * invpd(crossprod(X)) * t(X); // - hat matrix  
//    if(I==2) for(i=1; i<B; i++) H = - Xcol * invpd(crossprod(Xcol)) * t(Xcol); // - hat matrix  

    // compare two ways to compute t(Y) * H * Y
    Matrix<double> H = X * t(X);
    double out=0;
    if(I==1) for(i=1; i<B; i++) out=(t(Y) * H * Y)[0]; 
    if(I==2) for(i=1; i<B; i++) {for (j=1; j<n; j++) for (int k=1; k<n; k++) out += Y[j] * Y[k] * H(j,k);} 
          
          
    UNPROTECT(1);
    return _ans;
  }
  
    
}

#endif





//
////an older copy
//  // assume X and Y are sorted
//  // nLower and nUpper are 1-based index
////  SEXP boot_chngpt(SEXP _X, SEXP _Y, SEXP _indices, SEXP _nLower, SEXP _nUpper, SEXP _B)
//  SEXP boot_grid_search_2(SEXP _X, SEXP _Y, SEXP _nLower, SEXP _nUpper, SEXP _B)
//  {
//    int i; //for loop index
//    // put _X and _Y into Matrixes X and Y
//    // note that the rows and colns are organized in a way now that they can be directly casted and there is no need to do things as in MCMCpack MCMCmetrop1R.cc
//    double* X_dat = REAL(_X);
//    const int n = nrows(_X);
//    const int p = ncols(_X);
//    Matrix<double,Col,Concrete> Xcol (n, p, X_dat); //column major 
//    // convert to row major so that creating bootstrap datasets can be faster   
//    Matrix<double,Row,Concrete> X(Xcol); // row major        
//    //for (i=0; i<n; i++) {for (int j=0; j<p; j++)  PRINTF("%f ", Xb(i,j)); PRINTF("\n");}        
//    double *Y_dat=REAL(_Y);
//    Matrix <> Y (n, 1, Y_dat);
//    // bootstrap replicate
//    double B = asReal(_B);
//    // bounds
//    int nLower=asInteger (_nLower);
//    int nUpper=asInteger (_nUpper);
//    
////    //indices
////	int i_indices=0;
////    int* indices_dat= INTEGER(_indices);
////    Matrix <> indices ((int)B, n, indices_dat);    
//    
//    // output
//    SEXP _ans=PROTECT(allocVector(REALSXP, B*(p+2)));// p slopes, 1 threshold, 1 goodness of fit stat
//    double *ans=REAL(_ans);    
//    
//	vector<int> index(n);
//	// these variables are reused within each bootstrap replicate
//    Matrix <double,Row,Concrete> Xb (n, p);
//    Matrix <> Yb (n, 1);
//	Matrix <double,Row,Concrete> estimates(nUpper-nLower+1, p);
//    Matrix<double> H, J, H1;
//	vector<double> logliks(nUpper-nLower+1); 
//	vector<double> thresholds(nUpper-nLower+1); 
//	int chosen;
//    for (int b=0; b<B; b++) {
//        
//        // fill index, note that it is 1-based
//        SampleReplace(n, n, &(index[0]));
////        for (i=0; i<n; i++) index[i]=indices_dat[i_indices++];
////        for (i=0; i<n; i++) PRINTF("%i ", index[i]); PRINTF("\n");              
//        // sort index
//        sort (index.begin(), index.end());
//        // create bootstrap dataset
//        for (i=0; i<n; i++) { 
//            Xb(i,_)=X(index[i]-1,_); //-1 to become 0-based
//            Yb(i)=Y(index[i]-1); 
//        }
//        
//        // save a copy of x as thresholds because they get changed below
//        for(i=nLower-1; i<nUpper; i++) thresholds[i-nLower+1]=Xb(i,p-1);
//        // set the x lower than nLower to 0
//        for(i=0; i<nLower-1; i++) Xb(i,p-1)=0;        
//        
//        // regression models
//        double e;
//        for(i=nLower-1; i<nUpper; i++) {    
//            // update the change point variable in Xb
//            // note that Xb(,p-1) saves the delta between successive x after the first iteration
//            e=Xb(i,p-1); 
//            for (int j=i; j<n; j++) Xb(j,p-1)=Xb(j,p-1)-e;
//            //for (int j=0; j<n; j++) PRINTF("%f ", Xb(j,p-1)); PRINTF("\n");
//            J = invpd(crossprod(Xb)) * t(Xb);
//            estimates(i-(nLower-1),_)=J * Yb;
//            //H = Xb * J; // - hat matrix
//            H1 = Xb * invpd(crossprod(Xb)) * t(Xb); 
//            //for (int j=0; j<n; j++) H(j,j)=1+H(j,j); // I - H
//            logliks[i-(nLower-1)] = (t(Yb) * H1 * Yb)(0);
//        }        
//        for(i=nLower-1; i<nUpper; i++)  PRINTF("%f ", logliks[i-(nLower-1)]); PRINTF("\n");
//        
//        // save the estimated coefficients and threshold
//        chosen = distance(logliks.begin(), max_element(logliks.begin(), logliks.end())); PRINTF("%i ", chosen); PRINTF("\n");
//        for (i=0; i<p; i++) ans[b*(p+2)+i]=estimates(chosen,i);
//        ans[b*(p+2)+p]=thresholds[chosen];
//        ans[b*(p+2)+p+1]=logliks[chosen];
//        
//    }
//     
//    UNPROTECT(1);
//    return _ans;
//  }
