
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

double M33c_search(
    Matrix<double,Row>& B, Matrix<double,Row>& r, vector<double>& x, vector<double>& w, 
    int * thresholdIdx, vector<double>& thresholds, int nThresholds,
    // forward cumulative sum
    Matrix<double,Row>& Bcusum, vector<double>& rcusum, vector<double>& Wcusum, vector<double>& xcusum, 
    vector<double>& x2cusum, vector<double>& x3cusum, vector<double>& x4cusum, vector<double>& x5cusum, 
    vector<double>& xrcusum, vector<double>& x2rcusum, Matrix<double,Row>& xBcusum, Matrix<double,Row>& x2Bcusum, 
    // reverse cumulative sum
    Matrix<double,Row>& BcusumR, vector<double>& rcusumR, vector<double>& WcusumR, vector<double>& xcusumR, 
    vector<double>& x2cusumR, vector<double>& x3cusumR, vector<double>& x4cusumR, vector<double>& x5cusumR, 
    vector<double>& xrcusumR, vector<double>& x2rcusumR, Matrix<double,Row>& xBcusumR, Matrix<double,Row>& x2BcusumR, 
    double* logliks)
{

    // loop index. 
    int i,j,tt;    
    int k; // k has a meaning in the algorithm 
    int chosen=0;//initialize to get rid of compilation warning
    
    int n=B.rows();
    int p=B.cols()+1;
    
    // need a copy of x in the calculation
    vector<double> x_cpy(n);
    for(i=0; i<n; i++) x_cpy[i] = x[i];        
    
    if(p>1) BcusumR(n-1,_) = B(n-1,_)* w[n-1];    
    rcusumR[n-1] = r(n-1)* w[n-1];        
    WcusumR[n-1] = w[n-1]* w[n-1];              
    xcusumR[n-1]    = x[n-1] * w[n-1]* w[n-1]; // xcusumR is the last columnn of BcusumR, but perhaps better to be its own variable
    x2cusumR[n-1]   = pow(x[n-1],2) * w[n-1]* w[n-1];        
    x3cusumR[n-1]   = pow(x[n-1],3) * w[n-1]* w[n-1];        
    x4cusumR[n-1]   = pow(x[n-1],4) * w[n-1]* w[n-1];        
    x5cusumR[n-1]   = pow(x[n-1],5) * w[n-1]* w[n-1];        
    xrcusumR[n-1]   = x[n-1] * r(n-1)* w[n-1];        
    x2rcusumR[n-1]  = pow(x[n-1],2) * r(n-1)* w[n-1];        
    if(p>1) xBcusumR(n-1,_) = x[n-1] * B(n-1,_)* w[n-1];    
    if(p>1) x2BcusumR(n-1,_)= pow(x[n-1],2) * B(n-1,_)* w[n-1];    
    for (i=n-2; i>=0; i--) {
        if(p>1) BcusumR(i,_)  = BcusumR(i+1,_)    + B(i,_)* w[i]; 
        rcusumR[i]    = rcusumR[i+1]      + r(i)* w[i]; 
        WcusumR[i]    = WcusumR[i+1]      + w[i]* w[i]; 
        xcusumR[i]    = xcusumR[i+1]      + x[i]* w[i]* w[i]; 
        x2cusumR[i]   = x2cusumR[i+1]     + pow(x[i],2)* w[i]* w[i]; 
        x3cusumR[i]   = x3cusumR[i+1]     + pow(x[i],3)* w[i]* w[i]; 
        x4cusumR[i]   = x4cusumR[i+1]     + pow(x[i],4)* w[i]* w[i]; 
        x5cusumR[i]   = x5cusumR[i+1]     + pow(x[i],5)* w[i]* w[i]; 
        xrcusumR[i]   = xrcusumR[i+1]     + x[i]* r(i)* w[i]; 
        x2rcusumR[i]  = x2rcusumR[i+1]    + pow(x[i],2)* r(i)* w[i]; 
        if(p>1) xBcusumR(i,_) = xBcusumR(i+1,_)   + x[i]* B(i,_)* w[i]; 
        if(p>1) x2BcusumR(i,_)= x2BcusumR(i+1,_)  + pow(x[i],2)* B(i,_)* w[i]; 
    }
    
    
    if(p>1) Bcusum(0,_) = B(0,_)* w[0];    
    rcusum[0] = r(0)* w[0];        
    Wcusum[0] = w[0]* w[0];              
    xcusum[0]    = x[0] * w[0]* w[0]; 
    x2cusum[0]   = pow(x[0],2) * w[0]* w[0];        
    x3cusum[0]   = pow(x[0],3) * w[0]* w[0];        
    x4cusum[0]   = pow(x[0],4) * w[0]* w[0];        
    x5cusum[0]   = pow(x[0],5) * w[0]* w[0];        
    xrcusum[0]   = x[0] * r(0)* w[0];        
    x2rcusum[0]  = pow(x[0],2) * r(0)* w[0];        
    if(p>1) xBcusum(0,_) = x[0] * B(0,_)* w[0];    
    if(p>1) x2Bcusum(0,_)= pow(x[0],2) * B(0,_)* w[0];    
    for (i=1; i<n; i++) {
        if(p>1) Bcusum(i,_)  = Bcusum(i-1,_)    + B(i,_)* w[i]; 
        rcusum[i]    = rcusum[i-1]      + r(i)* w[i]; 
        Wcusum[i]    = Wcusum[i-1]      + w[i]* w[i]; 
        xcusum[i]    = xcusum[i-1]      + x[i]* w[i]* w[i]; 
        x2cusum[i]   = x2cusum[i-1]     + pow(x[i],2)* w[i]* w[i]; 
        x3cusum[i]   = x3cusum[i-1]     + pow(x[i],3)* w[i]* w[i]; 
        x4cusum[i]   = x4cusum[i-1]     + pow(x[i],4)* w[i]* w[i]; 
        x5cusum[i]   = x5cusum[i-1]     + pow(x[i],5)* w[i]* w[i]; 
        xrcusum[i]   = xrcusum[i-1]     + x[i]* r(i)* w[i]; 
        x2rcusum[i]  = x2rcusum[i-1]    + pow(x[i],2)* r(i)* w[i]; 
        if(p>1) xBcusum(i,_) = xBcusum(i-1,_)   + x[i]* B(i,_)* w[i]; 
        if(p>1) x2Bcusum(i,_)= x2Bcusum(i-1,_)  + pow(x[i],2)* B(i,_)* w[i]; 
    }
    
    //for (i=0; i<n; i++) {for (j=0; j<p; j++)  PRINTF("%f ", Bcusum(i,j)); PRINTF("\n");}        
    //for (i=0; i<n; i++) PRINTF("%f ", rcusum[i]); PRINTF("\n");
    
    
    // Step 2: Initialize VV, Vr and VB         
    
    Matrix <double,Row,Concrete> VV(4, 4), Vr(4, 1), VB(4, p-1);
    double e, d, d2, d3, crit, crit_max=R_NegInf;// d2 is defined as difference in squared threshold
    
    VV(2,3)=VV(3,2)=0;
    
    // x-e, not thresholded
    for(i=0; i<n; i++) x[i]=x[i]-thresholds[0]; 
    VV(0,0)=0; for (i=0; i<n; i++) VV(0,0) += pow(x[i],2)* w[i]* w[i]; 
    VV(1,1)=0; for (i=0; i<n; i++) VV(1,1) += pow(x[i],4)* w[i]* w[i]; 
    VV(0,1)=0; for (i=0; i<n; i++) VV(0,1) += pow(x[i],3)* w[i]* w[i]; 
    VV(1,0)=VV(0,1);
    
    Vr(0,0)=0; for (i=0; i<n; i++) Vr(0,0) += x[i]*r(i)* w[i]; 
    Vr(1,0)=0; for (i=0; i<n; i++) Vr(1,0) += pow(x[i],2)*r(i)* w[i]; 
    
    for (j=0; j<p-1; j++) {
        VB(0,j)=0; for (i=0; i<n; i++) VB(0,j) += x[i] * B(i,j)* w[i]; // remember the first p-1 column of B is B after _preprocess()
        VB(1,j)=0; for (i=0; i<n; i++) VB(1,j) += pow(x[i],2) * B(i,j)* w[i]; // remember the first p-1 column of B is B after _preprocess()
    }
    // threshold x (stored in the last column of B) at the first threshold value
    // Upper hinge (x-e)-
    for(i=0; i<thresholdIdx[0]; i++) x[i]=x_cpy[i]-thresholds[0];
    for(i=thresholdIdx[0]; i<n; i++) x[i]=0; 
    VV(2,2)=0; for (i=0; i<n; i++) VV(2,2) += pow(x[i],6)* w[i]* w[i]; 
    VV(0,2)=0; for (i=0; i<n; i++) VV(0,2) += pow(x[i],4)* w[i]* w[i];
    VV(1,2)=0; for (i=0; i<n; i++) VV(1,2) += pow(x[i],5)* w[i]* w[i];
    VV(2,0)=VV(0,2);
    VV(2,1)=VV(1,2);
    
    Vr(2,0)=0; for (i=0; i<n; i++) Vr(2,0) += pow(x[i],3)*r(i)* w[i]; 
    
    for (j=0; j<p-1; j++) {
        VB(2,j)=0; for (i=0; i<n; i++) VB(2,j) += pow(x[i],3) * B(i,j)* w[i]; // remember the first p-1 column of B is B after _preprocess()
    }
    // threshold x (stored in the last column of B) at the first threshold value
    // Hinge (x-e)+
    for(i=0; i<thresholdIdx[0]; i++) x[i]=0;  
    for(i=thresholdIdx[0]; i<n; i++) x[i]=x_cpy[i]-thresholds[0]; 
    VV(3,3)=0; for (i=0; i<n; i++) VV(3,3) += pow(x[i],6)* w[i]* w[i]; 
    VV(0,3)=0; for (i=0; i<n; i++) VV(0,3) += pow(x[i],4)* w[i]* w[i];
    VV(1,3)=0; for (i=0; i<n; i++) VV(1,3) += pow(x[i],5)* w[i]* w[i];
    VV(3,0)=VV(0,3);
    VV(3,1)=VV(1,3);
    
    Vr(3,0)=0; for (i=0; i<n; i++) Vr(3,0) += pow(x[i],3)*r(i)* w[i]; 
    
    for (j=0; j<p-1; j++) {
        VB(3,j)=0; for (i=0; i<n; i++) VB(3,j) += pow(x[i],3) * B(i,j)* w[i]; // remember the first p-1 column of B is B after _preprocess()
    }
        
    // Step 4 (Step 3 = Step 4b)
    for(tt=0; tt<nThresholds; tt++) { 
            
        // Step 4a: update VV, Vr and VB     
        if(tt>0) {
            d=thresholds[tt]-thresholds[tt-1]; //B(tt,p-1)-B(tt-1,p-1); //PRINTF("%f ", d); PRINTF("\n"); 
            d2=pow(thresholds[tt],2) - pow(thresholds[tt-1],2);
            d3=pow(thresholds[tt],3) - pow(thresholds[tt-1],3);
            
            // e.g, when tt=1, we are at the second threshold, we want to update vv, vr and vB for e_2
            //            tt=1 and tt+1=2 
            k  = thresholdIdx[tt-1]; // k is 1-based index of x, e_t = x_k
            e = thresholds[tt-1];
                        
            double c7=-2*d*xcusum[n-1] + 2*Wcusum[n-1]*d*e + Wcusum[n-1]* pow(d,2); 
            double c8=-4*d*x3cusum[k-1] + (6*d*e+3*d2+3*d*d)*x2cusum[k-1] - (3*d*e*e+3*d2*e+d3+3*d*d2)*xcusum[k-1] + Wcusum[k-1]*d*e*e*e + Wcusum[k-1]*d3*e + Wcusum[k-1]*d*d3;
            double c9=-5*d*x4cusum[k-1] + (12*d*e+4*d2+6*d*d)*x3cusum[k-1] - (9*d*e*e+9*d2*e+d3+9*d*d2)*x2cusum[k-1] + (2*d*e*e*e+6*d2*e*e+2*d3*e+2*d*d3+3*d2*d2)*xcusum[k-1] - (Wcusum[k-1]*d2*e*e*e + Wcusum[k-1]*d3*e*e + Wcusum[k-1]*d2*d3); 
            double c10=-6*d*x5cusum[k-1] + (18*d*e+6*d2+9*d*d)*x4cusum[k-1] - (18*d*e*e+18*d2*e+2*d3+18*d*d2)*x3cusum[k-1] + (6*d*e*e*e+18*d2*e*e+6*d3*e+6*d*d3+9*d2*d2)*x2cusum[k-1] - (6*d2*e*e*e+6*d3*e*e+6*d2*d3)*xcusum[k-1] + (2*Wcusum[k-1]*d3*e*e*e + Wcusum[k-1]*d3*d3); 
                                        
            double c11= -3*d*x2cusum[n-1] + (4*d*e+d2+2*d*d)*xcusum[n-1] -(Wcusum[n-1]*d*e*e + Wcusum[n-1]*d2*e + Wcusum[n-1]*d*d2);
            double c12= -4*d*x3cusumR[k] + (6*d*e+3*d2+3*d*d)*x2cusumR[k] - (3*d*e*e+3*d2*e+d3+3*d*d2)*xcusumR[k] + WcusumR[k]*d*e*e*e + WcusumR[k]*d3*e + WcusumR[k]*d*d3;
            double c13= -4*d*x3cusum[n-1] + (8*d*e+2*d2+4*d*d)*x2cusum[n-1] - (4*d*e*e+4*d2*e+4*d*d2)*xcusum[n-1] + 2*Wcusum[n-1]*d2*e*e + Wcusum[n-1]*d2*d2;
            double c14= -5*d*x4cusumR[k] + (12*d*e+4*d2+6*d*d)*x3cusumR[k] - (9*d*e*e+9*d2*e+d3+9*d*d2)*x2cusumR[k] + (2*d*e*e*e+6*d2*e*e+2*d3*e+2*d*d3+3*d2*d2)*xcusumR[k] - (WcusumR[k]*d2*e*e*e + WcusumR[k]*d3*e*e + WcusumR[k]*d2*d3); 
            double c15= -6*d*x5cusumR[k] + (18*d*e+6*d2+9*d*d)*x4cusumR[k] - (18*d*e*e+18*d2*e+2*d3+18*d*d2)*x3cusumR[k] + (6*d*e*e*e+18*d2*e*e+6*d3*e+6*d*d3+9*d2*d2)*x2cusumR[k] - (6*d2*e*e*e+6*d3*e*e+6*d2*d3)*xcusumR[k] + (2*WcusumR[k]*d3*e*e*e + WcusumR[k]*d3*d3); 
                            
            VV(0,0) += c7;
            VV(1,1) += c13;
            VV(2,2) += c10;
            VV(3,3) += c15;
            VV(0,1) += c11; VV(1,0) = VV(0,1);
            VV(0,2) += c8;  VV(2,0) = VV(0,2);
            VV(0,3) += c12; VV(3,0) = VV(0,3);
            VV(1,2) += c9;  VV(2,1) = VV(1,2);
            VV(1,3) += c14; VV(3,1) = VV(1,3);

            Vr[0] +=   -d * rcusum[n-1];
            Vr[1] += -2*d * xrcusum[n-1] + d2 * rcusum[n-1] ;
            Vr[2] += -3*d * x2rcusum[k-1] +3*d2 * xrcusum[k-1] -d3 * rcusum[k-1] ;
            Vr[3] += -3*d * x2rcusumR[k]  +3*d2 * xrcusumR[k]  -d3 * rcusumR[k] ;
            
            for (j=0; j<p-1; j++) {
                VB(0,j) +=   -d * Bcusum(n-1,j); 
                VB(1,j) += -2*d * xBcusum(n-1,j) + d2 * Bcusum(n-1,j); 
                VB(2,j) += -3*d * x2Bcusum(k-1,j) +3*d2 * xBcusum(k-1,j) -d3 * Bcusum(k-1,j); 
                VB(3,j) += -3*d * x2BcusumR(k ,j) +3*d2 * xBcusumR(k ,j) -d3 * BcusumR(k ,j); 
            }
    
//            if(tt==1) {
//                for (i=0; i<4; i++) {for (j=0; j<4; j++)   PRINTF("%f ", VV(i,j)); PRINTF("\n");}        
//                for (i=0; i<4; i++) {for (j=0; j<1; j++)   PRINTF("%f ", Vr(i,j)); PRINTF("\n");}        
//                for (i=0; i<4; i++) {for (j=0; j<p-1; j++) PRINTF("%f ", VB(i,j)); PRINTF("\n");}        
//            }
        } // end if tt==0

//            if(tt==0) {
//                for (i=0; i<4; i++) {for (j=0; j<4; j++)   PRINTF("%f ", VV(i,j)); PRINTF("\n");}        
//                for (i=0; i<4; i++) {for (j=0; j<1; j++)   PRINTF("%f ", Vr(i,j)); PRINTF("\n");}        
//                for (i=0; i<4; i++) {for (j=0; j<p-1; j++) PRINTF("%f ", VB(i,j)); PRINTF("\n");}        
//            }

                   
        // Step 4b: compute r'H_eY - r'HY
        if (p>1) {
            crit = ( t(Vr) * chol_solve(VV - VB * t(VB), Vr) ) (0,0);
        } else { // VB has 0 columns
            crit = ( t(Vr) * chol_solve(VV, Vr) ) (0,0);
        }
            
        logliks[tt] = crit;
        if(crit>=crit_max) {
            chosen = tt;
            crit_max=crit;
        }             
    }        
    //PRINTF("logliks: \n");  for(i=0; i<nThresholds; i++) PRINTF("%f ", logliks[i]); PRINTF("\n");  

    return thresholds[chosen]; 
    
}


#endif
//#endif
