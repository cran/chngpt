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


  // unit testing for performance comparison
SEXP performance_unit_test(SEXP u_X, SEXP u_Y, SEXP u_nBoot, SEXP u_I)
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
    int nBoot = asInteger(u_nBoot); nBoot=nBoot*1;
    int I = asInteger(u_I); I=I*1;

    // output
    SEXP _ans=PROTECT(allocVector(REALSXP, 1));
    double *ans=REAL(_ans);    
    ans[0]=0; // to avoid compiler warning
    
//    // compare X and Xcol in computing H
//    Matrix<double> H;
//    if(I==1) for(i=1; i<nBoot; i++) H = - X * invpd(crossprod1(X)) * t(X); // - hat matrix  
//    if(I==2) for(i=1; i<nBoot; i++) H = - Xcol * invpd(crossprod1(Xcol)) * t(Xcol); // - hat matrix  

//    // compare two ways to co//mpute t(Y) * H * Y
//    Matrix<double> H = X * t(X);
//    double out=0;
//    if(I==1) for(i=1; i<nBoot; i++) out=(t(Y) * H * Y)[0]; 
//    if(I==2) for(i=1; i<nBoot; i++) {for (j=1; j<n; j++) for (int cntr=1; cntr<n; cntr++) out += Y[j] * Y[cntr] * H(j,cntr);} 

    // compare crossprod1 row-major and column-major input
    Matrix<double> A=crossprod1(Xcol);
    Matrix<double> A1=crossprod1(X);
    for (i=0; i<p; i++) {for (j=0; j<p; j++)  PRINTF("%f ", A(i,j)); PRINTF("\n");}        
    for (i=0; i<p; i++) {for (j=0; j<p; j++)  PRINTF("%f ", A1(i,j)); PRINTF("\n");}        
    
    
    double* Xpnt = X.getArray();
    double* Xcolpnt = Xcol.getArray();
    for (i=0; i<n*p; i++) PRINTF("%f ", Xcolpnt[i]); 
    PRINTF("\n");       
    for (i=0; i<n*p; i++) PRINTF("%f ", Xpnt[i]); 
    PRINTF("\n");       
    PRINTF("%i %i\n",(int) X.rows(), (int) X.cols());       
    PRINTF("%i %i\n",(int) Xcol.rows(),(int) Xcol.cols());       
          
          
    UNPROTECT(1);
    return _ans;
    
}
  
    
}

#endif
//#endif


