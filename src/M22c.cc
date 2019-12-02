
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

       
double M22c_search(
    Matrix<double,Row>& B, Matrix<double,Row>& r, vector<double>& x, vector<double>& w, 
    int * thresholdIdx, vector<double>& thresholds, int nThresholds, 
    Matrix<double,Row>& Bcusum, vector<double>& rcusum, vector<double>& Wcusum, vector<double>& xcusum, 
    vector<double>& x2cusum, vector<double>& x3cusum, vector<double>& xrcusum, Matrix<double,Row>& xBcusum, 
    Matrix<double,Row>& BcusumR, vector<double>& rcusumR, vector<double>& WcusumR, vector<double>& xcusumR, 
    vector<double>& x2cusumR, vector<double>& x3cusumR, vector<double>& xrcusumR, Matrix<double,Row>& xBcusumR, 
    double* logliks)
{

    // loop index. 
    int i,j,tt;    
    int k; // k has a meaning in the algorithm 
    int chosen=0;//initialize to get rid of compilation warning
    
    int n=B.rows();
    int p=B.cols()+1;
    
    // a copy is needed because it goes back and forth
    vector<double> x_cpy(n);
    for(i=0; i<n; i++) x_cpy[i] = x[i];        
    //for (i=0; i<n; i++) PRINTF("%f ", x_cpy[i]); PRINTF("\n");
        
    // compute cusum of B and r in the forward order (upper hinge)
    if(p>1) Bcusum(0,_) = B(0,_)* w[0];    
    rcusum[0] = r(0)* w[0];        
    Wcusum[0] = w[0]* w[0];              
    xcusum[0]    = x_cpy[0] * w[0]* w[0]; // xcusum is the last columnn of Bcusum, but perhaps better to be its own variable
    x2cusum[0]   = pow(x_cpy[0],2) * w[0]* w[0];        
    x3cusum[0]   = pow(x_cpy[0],3) * w[0]* w[0];        
    xrcusum[0]   = x_cpy[0] * r(0)* w[0];        
    if(p>1) xBcusum(0,_) = x_cpy[0] * B(0,_)* w[0];    
    for (i=1; i<n; i++) {
        if(p>1) Bcusum(i,_)  = Bcusum(i-1,_)    + B(i,_)* w[i]; 
        rcusum[i]    = rcusum[i-1]      + r(i)* w[i]; 
        Wcusum[i]    = Wcusum[i-1]      + w[i]* w[i]; 
        xcusum[i]    = xcusum[i-1]      + x_cpy[i]* w[i]* w[i]; 
        x2cusum[i]   = x2cusum[i-1]     + pow(x_cpy[i],2)* w[i]* w[i]; 
        x3cusum[i]   = x3cusum[i-1]     + pow(x_cpy[i],3)* w[i]* w[i]; 
        xrcusum[i]   = xrcusum[i-1]     + x_cpy[i]* r(i)* w[i]; 
        if(p>1) xBcusum(i,_) = xBcusum(i-1,_)   + x_cpy[i]* B(i,_)* w[i]; 
    }
    
    // compute cusum of B and r in the reverse order (hinge)
    for (i=n-1; i>=1; i--) {
        if(p>1) BcusumR(i,_)  = Bcusum(n-1,_)    - Bcusum(i-1,_);
        rcusumR[i]    = rcusum[n-1]      - rcusum[i-1]; 
        WcusumR[i]    = Wcusum[n-1]      - Wcusum[i-1]; 
        xcusumR[i]    = xcusum[n-1]      - xcusum[i-1]; 
        x2cusumR[i]   = x2cusum[n-1]     - x2cusum[i-1]; 
        x3cusumR[i]   = x3cusum[n-1]     - x3cusum[i-1]; 
        xrcusumR[i]   = xrcusum[n-1]     - xrcusum[i-1]; 
        if(p>1) xBcusumR(i,_) = xBcusum(n-1,_)   - xBcusum(i-1,_); 
    }
    if(p>1) BcusumR(0,_)  = Bcusum(n-1,_);
    rcusumR[0]    = rcusum[n-1]; 
    WcusumR[0]    = Wcusum[n-1]; 
    xcusumR[0]    = xcusum[n-1]; 
    x2cusumR[0]   = x2cusum[n-1]; 
    x3cusumR[0]   = x3cusum[n-1]; 
    xrcusumR[0]   = xrcusum[n-1]; 
    if(p>1) xBcusumR(0,_) = xBcusum(n-1,_); 


    //for (i=0; i<n; i++) {for (j=0; j<p; j++)  PRINTF("%f ", Bcusum(i,j)); PRINTF("\n");}        
    //for (i=0; i<n; i++) PRINTF("%f ", rcusum[i]); PRINTF("\n");
    //for (i=0; i<n; i++) {for (j=0; j<p; j++)  PRINTF("%f ", BcusumR(i,j)); PRINTF("\n");}        
    //for (i=0; i<n; i++) PRINTF("%f ", rcusumR[i]); PRINTF("\n");
        
    // Step 2: Initialize VV, Vr and VB     
    Matrix <double,Row,Concrete> VV(3, 3), Vr(3, 1), VB(3, p-1);    
    VV(1,2)=VV(2,1)=0;
    // x-e
    for(i=0; i<n; i++) x[i]=x_cpy[i]-thresholds[0]; 
    VV(0,0)=0; for (i=0; i<n; i++) VV(0,0) += pow(x[i],2)* w[i]* w[i]; 
    
    Vr(0,0)=0; for (i=0; i<n; i++) Vr(0,0) += x[i]*r(i)* w[i]; 
    
    for (j=0; j<p-1; j++) {
        VB(0,j)=0; for (i=0; i<n; i++) VB(0,j) += x[i] * B(i,j)* w[i]; // remember the first p-1 column of B is B after _preprocess()
    }
    // threshold x (stored in the last column of B) at the first threshold value
    // Upper hinge (x-e)-
    for(i=0; i<thresholdIdx[0]; i++) x[i]=x_cpy[i]-thresholds[0];
    for(i=thresholdIdx[0]; i<n; i++) x[i]=0; 
    VV(1,1)=0; for (i=0; i<n; i++) VV(1,1) += pow(x[i],4)* w[i]* w[i]; 
    VV(0,1)=0; for (i=0; i<n; i++) VV(0,1) += pow(x[i],3)* w[i]* w[i];
    VV(1,0)=VV(0,1);
    
    Vr(1,0)=0; for (i=0; i<n; i++) Vr(1,0) += pow(x[i],2)*r(i)* w[i]; 
    
    for (j=0; j<p-1; j++) {
        VB(1,j)=0; for (i=0; i<n; i++) VB(1,j) += pow(x[i],2) * B(i,j)* w[i]; // remember the first p-1 column of B is B after _preprocess()
    }
    // threshold x (stored in the last column of B) at the first threshold value
    // Hinge (x-e)+
    for(i=0; i<thresholdIdx[0]; i++) x[i]=0;  
    for(i=thresholdIdx[0]; i<n; i++) x[i]=x_cpy[i]-thresholds[0]; 
    VV(2,2)=0; for (i=0; i<n; i++) VV(2,2) += pow(x[i],4)* w[i]* w[i]; 
    VV(0,2)=0; for (i=0; i<n; i++) VV(0,2) += pow(x[i],3)* w[i]* w[i];
    VV(2,0)=VV(0,2);
    
    Vr(2,0)=0; for (i=0; i<n; i++) Vr(2,0) += pow(x[i],2)*r(i)* w[i]; 
    
    for (j=0; j<p-1; j++) {
        VB(2,j)=0; for (i=0; i<n; i++) VB(2,j) += pow(x[i],2) * B(i,j)* w[i]; // remember the first p-1 column of B is B after _preprocess()
    }
    //for (i=0; i<3; i++) {for (j=0; j<3; j++)   PRINTF("%f ", VV(i,j)); PRINTF("\n");}        
    //for (i=0; i<3; i++) {for (j=0; j<1; j++)   PRINTF("%f ", Vr(i,j)); PRINTF("\n");}        
    //for (i=0; i<3; i++) {for (j=0; j<p-1; j++) PRINTF("%f ", VB(i,j)); PRINTF("\n");}        
    
    // Step 4 (Step 3 = Step 4b)
    double delta, s, crit, crit_max=R_NegInf, e;// s is defined as difference in squared threshold
    for(tt=0; tt<nThresholds; tt++) { 
            
        // Step 4a: update VV, Vr and VB     
        if(tt>0) {
            e=thresholds[tt-1];
            delta= thresholds[tt]-e; //B(tt,p-1)-B(tt-1,p-1); //PRINTF("%f ", delta); PRINTF("\n"); 
            s=pow(thresholds[tt],2)-pow(e, 2);
            
            // x-e part
            VV(0,0) += -2*delta*xcusum[n-1] + 2*Wcusum[n-1]*delta*e + Wcusum[n-1] * pow(delta,2); 
            Vr[0] += -delta * rcusum[n-1];
            for (j=0; j<p-1; j++) {
                VB(0,j) +=   -delta * Bcusum(n-1,j); 
            }
            
            // upper hinge part
            // e.g, when tt=1, we are at the second threshold, we want to update vv, vr and vB for e_2
            //            tt=1 and tt+1=2 
            k  = thresholdIdx[tt-1]; // k is 1-based index of x, e_t = x_k
            VV(1,1) += -4*delta*x3cusum[k-1] + 2*(2*pow(delta,2)+4*delta*e+s)*x2cusum[k-1] - (4*delta*s+4*delta*pow(e,2)+4*e*s)*xcusum[k-1] + Wcusum[k-1]*s*s + 2*Wcusum[k-1]*s*pow(e,2); 
            VV(0,1) += -3*delta*x2cusum[k-1] + (2*pow(delta,2)+4*delta*e+s) * xcusum[k-1] - Wcusum[k-1]*delta*s - Wcusum[k-1]*s*e - Wcusum[k-1]*delta*pow(e,2); 
            VV(1,0) = VV(0,1);
            Vr[1] += -2*delta * xrcusum[k-1] + s * rcusum[k-1] ;
            for (j=0; j<p-1; j++) {
                VB(1,j) += -2*delta * xBcusum(k-1,j) + s * Bcusum(k-1,j); 
            }

            // hinge part
            k  = n-k; // the current threshold is the kth largest 
            VV(2,2) += -4*delta*x3cusumR[n-k] + 2*(2*pow(delta,2)+4*delta*e+s)*x2cusumR[n-k] - (4*delta*s+4*delta*pow(e,2)+4*e*s)*xcusumR[n-k] + WcusumR[n-k]*s*s + 2*WcusumR[n-k]*s*pow(e,2); 
            VV(0,2) += -3*delta*x2cusumR[n-k] +   (2*pow(delta,2)+4*delta*e+s) *xcusumR[n-k] -  WcusumR[n-k]*delta*s - WcusumR[n-k]*s*e - WcusumR[n-k]*delta*pow(e,2); 
            VV(2,0) = VV(0,2);
            Vr[2] += -2*delta * xrcusumR[n-k] + s * rcusumR[n-k];
            for (j=0; j<p-1; j++) {
                VB(2,j) += -2*delta * xBcusumR(n-k,j) + s * BcusumR(n-k,j); 
            }                        
        }
        
        if(tt==1) {
            //for (i=0; i<3; i++) {for (j=0; j<3; j++)   PRINTF("%f ", VV(i,j)); PRINTF("\n");}        
            //for (i=0; i<3; i++) {for (j=0; j<1; j++)   PRINTF("%f ", Vr(i,j)); PRINTF("\n");}        
            //for (i=0; i<3; i++) {for (j=0; j<p-1; j++) PRINTF("%f ", VB(i,j)); PRINTF("\n");}        
        }
            
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
