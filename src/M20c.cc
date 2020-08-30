
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


using namespace std;
using namespace scythe;


double M20c_search(
    Matrix<double,Row>& B, Matrix<double,Row>& r, vector<double>& x, vector<double>& w, 
    int * thresholdIdx, vector<double>& thresholds, int nThresholds,
    // forward cumulative sum
    Matrix<double,Row>& Bcusum, vector<double>& rcusum, vector<double>& Wcusum, vector<double>& xcusum, 
    vector<double>& x2cusum, vector<double>& x3cusum, vector<double>& xrcusum, Matrix<double,Row>& xBcusum, 
    double* logliks)
{

    // loop index. 
    int i,j,t;    
    int k; // k has a meaning in the algorithm 
    int chosen=0;//initialize to get rid of compilation warning
    
    int n=B.rows();
    int p=B.cols()+1;
    
    if(p>1) Bcusum(0,_) = B(0,_)* w[0];    
    rcusum[0] = r(0)* w[0];        
    Wcusum[0] = w[0]* w[0];              
    xcusum[0]    = x[0] * w[0]* w[0]; // xcusum is the last columnn of Bcusum, but perhaps better to be its own variable
    x2cusum[0]   = pow(x[0],2) * w[0]* w[0];        
    x3cusum[0]   = pow(x[0],3) * w[0]* w[0];        
    xrcusum[0]   = x[0] * r(0)* w[0];        
    if(p>1) xBcusum(0,_) = x[0] * B(0,_)* w[0];    
    for (i=1; i<n; i++) {
        if(p>1) Bcusum(i,_)  = Bcusum(i-1,_)    + B(i,_)* w[i]; 
        rcusum[i]    = rcusum[i-1]      + r(i)* w[i]; 
        Wcusum[i]    = Wcusum[i-1]      + w[i]* w[i]; 
        xcusum[i]    = xcusum[i-1]      + x[i]* w[i]* w[i]; 
        x2cusum[i]   = x2cusum[i-1]     + pow(x[i],2)* w[i]* w[i]; 
        x3cusum[i]   = x3cusum[i-1]     + pow(x[i],3)* w[i]* w[i]; 
        xrcusum[i]   = xrcusum[i-1]     + x[i]* r(i)* w[i]; 
        if(p>1) xBcusum(i,_) = xBcusum(i-1,_)   + x[i]* B(i,_)* w[i]; 
    }    
    //for (i=0; i<n; i++) {for (j=0; j<p; j++)  PRINTF("%f ", Bcusum(i,j)); PRINTF("\n");}        
    //for (i=0; i<n; i++) PRINTF("%f ", rcusum[i]); PRINTF("\n");
        
        
    // Step 2: Initialize VV, Vr and VB     
    
    // threshold x (stored in the last column of B) at the first threshold value
    // (x-e)-
    for(i=0; i<thresholdIdx[0]; i++) x[i]=x[i]-thresholds[0];  
    for(i=thresholdIdx[0]; i<n; i++) x[i]=0; 
//    PRINTF("thresholds[0] %f ", thresholds[0]); PRINTF("\n");  
//    for (i=0; i<n; i++) PRINTF("%f ", r(i)); PRINTF("\n");    
//    for (j=0; j<n; j++) PRINTF("%f ", B(j,p-1)); PRINTF("\n");    
    
    vector<double> vB(p-1); 
    double vv, vr;
    double delta, s, crit, crit_max=R_NegInf;

    vv=0; vr=0; 
    for (i=0; i<n; i++) vv += pow(x[i],4)*w[i]* w[i]; 
    for (i=0; i<n; i++) vr += pow(x[i],2)*r(i)* w[i]; 
    for (j=0; j<p-1; j++) {
        vB[j]=0; 
        for (i=0; i<n; i++) vB[j] += pow(x[i],2) * B(i,j)* w[i];
    }
    
    // Step 4 (Step 3 = Step 4b)
    for(t=0; t<nThresholds; t++) {             
        // Step 4a: update vv, vr and vB
        if(t>0) {
            delta= thresholds[t]-thresholds[t-1]; //B(t,p-1)-B(t-1,p-1); //PRINTF("%f ", delta); PRINTF("\n"); 
            s=pow(thresholds[t],2)-pow(thresholds[t-1], 2);
            
            // e.g, when t=1, we are at the second threshold, we want to update vv, vr and vB for e_2
            //            t=1 and t+1=2 
            k  = thresholdIdx[t-1]; // k is 1-based index of x, e_t = x_k
            
            vv += -4*delta*x3cusum[k-1] + 2*(2*pow(delta,2)+4*delta*thresholds[t-1]+s)*x2cusum[k-1] - (4*delta*s+4*delta*pow(thresholds[t-1],2)+4*thresholds[t-1]*s)*xcusum[k-1] + Wcusum[k-1]*s*s + 2*Wcusum[k-1]*s*pow(thresholds[t-1],2); 
                      
            vr += -2*delta * xrcusum[k-1] + s * rcusum[k-1] ;

            for (j=0; j<p-1; j++) vB[j] += -2*delta * xBcusum(k-1,j) + s * Bcusum(k-1,j); 
        }
        
        // Step 4b: compute r'H_eY - r'HY
        // first compute vB'vB, abusing the notation a little bit
        crit=0; for (j=0; j<p-1; j++) crit += pow(vB[j],2); 
        // now compute crit
        crit = vr*vr / (vv - crit);

//        if(t==0) {
//            PRINTF("vv %f \n", vv); 
//            PRINTF("vr %f \n", vr); 
//            PRINTF("vB\n"); for (j=0; j<p-1; j++)  PRINTF("%f ", vB[j]); PRINTF("\n");
//        }

        logliks[t] = crit;
        if(crit>=crit_max) {
            chosen = t;
            crit_max=crit;
        }             
    }        
    
    return thresholds[chosen]; 
    
}



#endif
//#endif
