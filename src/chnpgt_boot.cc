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



extern "C" {
  
  // it is not clear to me whether crossprod is calling lapack or not. crossprod1 is the way I make sure it is
  inline void make_symmetric(double* matrix, int rows)
  {
      for (int i = 1; i < rows; ++i)
        for (int j = 0; j < i; ++j)
          matrix[i * rows + j] = matrix[j * rows + i];
  }
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

  // assume X is sorted in chngptvar from small to large
  // nLower and nUpper are 1-based index
  void _fastgrid_search(Matrix<double,Row>& X, vector<double>& Y, double * w, bool wAllOne, int nLower, int nUpper, int n, int p, 
    Matrix<double,Row>& Xcusum, vector<double>& Ycusum, vector<double>& Wcusum, vector<double>& thresholds, vector<double>& Cps, 
    double* ans1, double* ans2)
  {

    // loop index
    int i,j;
    
	int chosen=0;//initialize to get rid of compilation warning
    vector<double> C(p); 

    // compute cusum of X and Y in the reverse order
    if(wAllOne) {
        Xcusum(n-1,_) = X(n-1,_);
        for (i=n-2; i>=0; i--) Xcusum(i,_) = Xcusum(i+1,_) + X(i,_); 
        Ycusum[n-1] = Y[n-1];
        for (i=n-2; i>=0; i--) Ycusum[i] = Ycusum[i+1] + Y[i]; 
    } else  {
        Xcusum(n-1,_) = X(n-1,_) * w[n-1];
        for (i=n-2; i>=0; i--) Xcusum(i,_) = Xcusum(i+1,_) + X(i,_) * w[i]; 
        Ycusum[n-1] = Y[n-1] * w[n-1];
        for (i=n-2; i>=0; i--) Ycusum[i] = Ycusum[i+1] + Y[i] * w[i]; 
    }
    Wcusum[n-1] = w[n-1];
    for (i=n-2; i>=0; i--) Wcusum[i] = Wcusum[i+1] + w[i]; 
    //for (i=0; i<n; i++) {for (j=0; j<p; j++)  PRINTF("%f ", Xcusum(i,j)); PRINTF("\n");}        
    //for (i=0; i<n; i++) PRINTF("%f ", Ycusum[i]); PRINTF("\n");
    
    // grid of e
    for(i=nLower-1; i<nUpper; i++) thresholds[i-nLower+1]=X(i,p-1);

    Matrix <double> A, Ainv, Ainv_save;
    double delta, rss, rss_min=R_PosInf;
    
    // Step 2: compuate A:=X_e'X_e and C:=X_e'Y at the first point in the grid
    // compute X_e
    for(i=0; i<nLower-1; i++) X(i,p-1)=0;  
    // change from nLower-1 to n-1, but has to be from back to front
    for(i=n-1; i>nLower-2; i--) X(i,p-1)=X(i,p-1)-X(nLower-1,p-1); 
    //for (j=0; j<n; j++) PRINTF("%f ", X(j,p-1)); PRINTF("\n");
    if(wAllOne) {
        for (j=0; j<p; j++) {C[j]=0; for (int k=0; k<n; k++) C[j] += X(k,j) * Y[k];}
        A = crossprod1(X); 
    } else {
        for (j=0; j<p; j++) {C[j]=0; for (int k=0; k<n; k++) C[j] += X(k,j) * Y[k] * w[k];}
        // mutiply X with sqrt(w) before crossprod. Note that this changes X! Don't use X after this
        for (i=0; i<n; i++) X(i,_)=X(i,_)*sqrt(w[i]);
        A = crossprod1(X);
    }
    Cps[0]=C[p-1];// save C[p-1]
            
    // loop through candidate thresholds
    for(i=nLower-1; i<nUpper; i++) {            
        // Step 4a: update A:=X'X and C:=X'Y
        // wAllOne handled by replacing n-i with Wcusum(i)
        if(i>nLower-1) {
            delta= thresholds[i-nLower+1]-thresholds[i-nLower]; //X(i,p-1)-X(i-1,p-1); //PRINTF("%f ", delta); PRINTF("\n"); 
            A(p-1,p-1) += pow(delta,2) * Wcusum[i]; // Delta' * W * Delta
            for (j=0; j<p-1; j++) A(p-1,j) -= delta * Xcusum(i,j); // Delta' * W * X
            for (j=0; j<p-1; j++) A(j,p-1) = A(p-1,j);  // X' * W * Delta
            A(p-1,p-1) -= 2 * delta * (Xcusum(i,p-1) - thresholds[i-nLower]*Wcusum[i]); // the last element of both Delta' * W * X and X' * W * Delta
            C[p-1] -= delta * Ycusum[i];
            Cps[i-nLower+1]=C[p-1];// save C[p-1]
            //for (int k=0; k<p; k++) for (j=0; j<p; j++)  PRINTF("%f ", A(k,j)); PRINTF("\n");        
        }
        
        // Step 3 and 4b: compute Y'HY                
        Ainv = invpd(A);
        rss=0; for (j=0; j<p; j++) for (int k=0; k<p; k++) rss -= C[j] * C[k] * Ainv(j,k);
        ans1[i-(nLower-1)] = rss;
        if(rss<=rss_min) {
            chosen = i-(nLower-1);
            Ainv_save=Ainv;
            rss_min=rss;
        }             
    }        

    // save results: estimated coefficients and threshold
    ans2[p]=thresholds[chosen];
    C[p-1]=Cps[chosen];
    for (i=0; i<p; i++) { ans2[i]=0; for (j=0; j<p; j++) ans2[i]+=Ainv_save(i,j)*C[j]; }            
    ans2[p+1]=-rss_min;
    
  }


  // assume X is sorted in chngptvar from small to large
  // nLower and nUpper are 1-based index
  SEXP fastgrid_search(SEXP u_X, SEXP u_Y, SEXP u_W, SEXP u_wAllOne, SEXP u_nLower, SEXP u_nUpper)
  {
    // input
    double* uX_dat = REAL(u_X);
    double *Y_dat=REAL(u_Y);
    double *W=REAL(u_W);
    bool wAllOne=asLogical(u_wAllOne)==1;
    const int n = nrows(u_X);
    const int p = ncols(u_X);
    int nLower=asInteger (u_nLower);
    int nUpper=asInteger (u_nUpper);
    
    // output
    SEXP _ans=PROTECT(allocVector(REALSXP, nUpper-nLower+1));
    double *ans=REAL(_ans);        
    
    // The rows and colns are organized in a way now that they can be directly casted and there is no need to do things as in the JSS paper on sycthe
    Matrix<double,Col,Concrete> Xcol (n, p, uX_dat); //column major 
    Matrix<double,Row,Concrete> X(Xcol); // convert to row major to pass to _fastgrid_search    
    vector<double> Y(n);
    for (int i=0; i<n; i++) Y[i]=Y_dat[i]; 

    Matrix <double,Row,Concrete> Xcusum (n, p, true, 0);
	vector<double> thresholds(nUpper-nLower+1), Cps(nUpper-nLower+1), Ycusum(n); 			
	vector<double> Wcusum(n);//double stats[p+2], Wcusum[n]; 
    double * stats = (double *) malloc((p+2) * sizeof(double));
    _fastgrid_search(X, Y, W, wAllOne, nLower, nUpper, n, p, Xcusum, Ycusum, Wcusum, thresholds, Cps, ans, stats); // stats not used here
        
    UNPROTECT(1);
    free(stats);
    return _ans;
  }
       
  
  // assume X and Y are sorted
  // nLower and nUpper are 1-based index
  SEXP boot_fastgrid_search(SEXP u_X, SEXP u_Y, SEXP u_W, SEXP u_wAllOne, SEXP u_nLower, SEXP u_nUpper, SEXP _B)
  {
    // put u_X and u_Y into Matrixes Xori and Y
    // note that the rows and colns are organized in a way now that they can be directly casted and there is no need to do things as in MCMCpack MCMCmetrop1R.cc
    double* uX_dat = REAL(u_X);
    const int n = nrows(u_X);
    const int p = ncols(u_X);
    double B = asReal(_B);
    int nLower=asInteger (u_nLower);
    int nUpper=asInteger (u_nUpper);
    
    // output
    SEXP _ans=PROTECT(allocVector(REALSXP, B*(p+2)));// p slopes, 1 threshold, 1 goodness of fit stat
    double *ans=REAL(_ans);    
    
    Matrix<double,Col,Concrete> Xcol (n, p, uX_dat); //column major 
    // convert to row major so that creating bootstrap datasets can be faster   
    Matrix<double,Row,Concrete> Xori(Xcol); // row major, name this Xori instead of X so that there is less a mistake of using X when Xb should be used below
    //for (i=0; i<n; i++) {for (j=0; j<p; j++)  PRINTF("%f ", Xb(i,j)); PRINTF("\n");}        
    double *Y=REAL(u_Y);
    double *W=REAL(u_W);
    bool wAllOne=asLogical(u_wAllOne)==1;
    
	// these variables are reused within each bootstrap replicate
	vector<int> index(n);
    Matrix <double,Row,Concrete> Xb (n, p), Xcusum (n, p, true, 0);
	vector<double> thresholds(nUpper-nLower+1), Cps(nUpper-nLower+1), Ycusum(n), Yb(n); 
	vector<double> Wcusum(n); //double rsses[nUpper-nLower+1], Wcusum[n];
    double * rsses = (double *) malloc((nUpper-nLower+1) * sizeof(double));

	
    for (int b=0; b<B; b++) {        
        // create bootstrap dataset, note that index is 1-based
        SampleReplace(n, n, &(index[0]));
        // Step 1: sort
        sort (index.begin(), index.end());
        for (int i=0; i<n; i++) { Xb(i,_)=Xori(index[i]-1,_); Yb[i]=Y[index[i]-1]; } //note that index need to -1 to become 0-based
        //for (i=0; i<n; i++) {for (j=0; j<p; j++)  PRINTF("%f ", Xb(i,j)); PRINTF("\n");} 
        _fastgrid_search(Xb, Yb, W, wAllOne, nLower, nUpper, n, p, Xcusum, Ycusum, Wcusum, thresholds, Cps, rsses, ans+b*(p+2)); // rsses not used here
    } 
     
    UNPROTECT(1);
    free(rsses);
    return _ans;
  }

  // weights not implemented
  // an optimized implementation that serves as a reference for more advanced algebraic optimization
  // assume X and Y are sorted
  // nLower and nUpper are 1-based index
  SEXP boot_grid_search(SEXP u_X, SEXP u_Y, SEXP u_W, SEXP u_wAllOne, SEXP u_nLower, SEXP u_nUpper, SEXP _B)
  {
    int i,j; //for loop index
    // put u_X and u_Y into Matrixes Xori and Y
    // note that the rows and colns are organized in a way now that they can be directly casted and there is no need to do things as in MCMCpack MCMCmetrop1R.cc
    double* uX_dat = REAL(u_X);
    const int n = nrows(u_X);
    const int p = ncols(u_X);
    Matrix<double,Col,Concrete> Xcol (n, p, uX_dat); //column major 
    // convert to row major so that creating bootstrap datasets can be faster   
    Matrix<double,Row,Concrete> Xori(Xcol); // row major        
    //for (i=0; i<n; i++) {for (j=0; j<p; j++)  PRINTF("%f ", Xb(i,j)); PRINTF("\n");}        
    double *Y_dat=REAL(u_Y);
    Matrix <> Y (n, 1, Y_dat);
    // bootstrap replicate
    double B = asReal(_B);
    // bounds
    int nLower=asInteger (u_nLower);
    int nUpper=asInteger (u_nUpper);
    //double *W=REAL(u_W);
    //bool wAllOne=asLogical(u_wAllOne)==1;
    
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
            H = Xb * invpd(crossprod1(Xb)) * t(Xb); // hat matrix
            //the next line offer a faster way than: -(t(Yb) * H * Yb)(0);
            rss=0; for (j=0; j<n; j++) for (int k=0; k<n; k++) rss -= Yb[j] * Yb[k] * H(j,k);
            logliks[i-(nLower-1)] = -rss; // since Yb'Yb does not depend on threshold, there is no need to compute it
        }        
        //for(i=nLower-1; i<nUpper; i++)  PRINTF("%f ", logliks[i-(nLower-1)]); PRINTF("\n");
        
        // save the estimated coefficients and threshold
        chosen = distance(logliks.begin(), max_element(logliks.begin(), logliks.end()));
        efinal=thresholds[chosen];
        // compute slope estimate. Alternative, slope estimates could be computed withint the previous loop, but it would be slower
        for(i=nLower-1+chosen; i<n; i++) Xb(i,p-1)=Xori(index[i]-1,p-1)-efinal; //-1 to become 0-based
        //for (i=0; i<n; i++) {for (j=0; j<p; j++)  PRINTF("%f ", Xb(i,j)); PRINTF("\n");}        
        estimatedSlope = invpd(crossprod1(Xb)) * t(Xb) * Yb;
        for (i=0; i<p; i++) ans[b*(p+2)+i]=estimatedSlope(i);
        ans[b*(p+2)+p]=efinal;
        ans[b*(p+2)+p+1]=logliks[chosen];
        
    }
     
    UNPROTECT(1);
    return _ans;
  }

  // weights not implemented
  // the function that performs grid search and returns the best chngpt 
  // assume X is sorted in chngptvar from small to large
  // nLower and nUpper are 1-based index
  SEXP grid_search(SEXP u_X, SEXP u_Y, SEXP u_W, SEXP u_wAllOne, SEXP u_nLower, SEXP u_nUpper)
  {

    // put u_X and u_Y into Matrixes X and Y
    // note that the rows and colns are organized in a way now that they can be directly casted and there is no need to do things as in the JSS paper on sycthe
    int i,j;
    double* uX_dat = REAL(u_X);
    const int n = nrows(u_X);
    const int p = ncols(u_X);
    Matrix <> X (n, p, uX_dat);
    double *Y_dat=REAL(u_Y);
    Matrix <> Y (n, 1, Y_dat);
    int nLower=asInteger (u_nLower);
    int nUpper=asInteger (u_nUpper);
    //double *W=REAL(u_W);
    //bool wAllOne=asLogical(u_wAllOne)==1;
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
        H = - X * invpd(crossprod1(X)) * t(X); // - hat matrix
        for (j=0; j<n; j++) H(j,j)=1+H(j,j); // I - H
        ans[i-(nLower-1)] = (t(Y) * H * Y)(0);
    }
          
    UNPROTECT(1);
    return _ans;
  }
  
  
  // unit testing for performance comparison
  SEXP performance_unit_test(SEXP u_X, SEXP u_Y, SEXP _B, SEXP _I)
  {

    int i,j; //for loop index
    // put u_X and u_Y into Matrixes X and Y
    // note that the rows and colns are organized in a way now that they can be directly casted and there is no need to do things as in MCMCpack MCMCmetrop1R.cc
    double* uX_dat = REAL(u_X);
    const int n = nrows(u_X);
    const int p = ncols(u_X);
    Matrix<double,Col,Concrete> Xcol (n, p, uX_dat); //column major 
    // convert to row major so that creating bootstrap datasets can be faster   
    Matrix<double,Row,Concrete> X(Xcol); // row major        
    //for (i=0; i<n; i++) {for (j=0; j<p; j++)  PRINTF("%f ", Xb(i,j)); PRINTF("\n");}        
    double *Y_dat=REAL(u_Y);
    Matrix <> Y (n, 1, Y_dat);
    // bootstrap replicate
    int B = asInteger(_B); B=B*1;
    int I = asInteger(_I); I=I*1;

    // output
    SEXP _ans=PROTECT(allocVector(REALSXP, 1));
    double *ans=REAL(_ans);    
    ans[0]=0; // to avoid compiler warning
    
//    // compare X and Xcol in computing H
//    Matrix<double> H;
//    if(I==1) for(i=1; i<B; i++) H = - X * invpd(crossprod1(X)) * t(X); // - hat matrix  
//    if(I==2) for(i=1; i<B; i++) H = - Xcol * invpd(crossprod1(Xcol)) * t(Xcol); // - hat matrix  

//    // compare two ways to co//mpute t(Y) * H * Y
//    Matrix<double> H = X * t(X);
//    double out=0;
//    if(I==1) for(i=1; i<B; i++) out=(t(Y) * H * Y)[0]; 
//    if(I==2) for(i=1; i<B; i++) {for (j=1; j<n; j++) for (int k=1; k<n; k++) out += Y[j] * Y[k] * H(j,k);} 

    // compare crossprod1 row-major and column-major input
    Matrix<double> A=crossprod1(Xcol);
    Matrix<double> A1=crossprod1(X);
    for (i=0; i<p; i++) {for (j=0; j<p; j++)  PRINTF("%f ", A(i,j)); PRINTF("\n");}        
    for (i=0; i<p; i++) {for (j=0; j<p; j++)  PRINTF("%f ", A1(i,j)); PRINTF("\n");}        
    
    double* Xpnt = X.getArray();
    double* Xcolpnt = Xcol.getArray();
    for (i=0; i<n*p; i++) PRINTF("%f ", Xcolpnt[i]); PRINTF("\n");       
    for (i=0; i<n*p; i++) PRINTF("%f ", Xpnt[i]); PRINTF("\n");       
    PRINTF("%i %i\n",(int) X.rows(), (int) X.cols());       
    PRINTF("%i %i\n",(int) Xcol.rows(),(int) Xcol.cols());       
          
          
    UNPROTECT(1);
    return _ans;
  }
  
    
}

#endif
//#endif

