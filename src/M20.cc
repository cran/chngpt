
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


double M20_search(
    Matrix<double,Row>& B, Matrix<double,Row>& r, vector<double>& x, vector<double>& w, 
    int * thresholdIdx, vector<double>& thresholds, int nThresholds,
    // forward cumulative sum
    Matrix<double,Row>& Bcusum, vector<double>& rcusum, vector<double>& Wcusum, vector<double>& xcusum, 
    vector<double>& x2cusum, vector<double>& x3cusum, vector<double>& xrcusum, Matrix<double,Row>& xBcusum, 
    double* logliks)
{

    // loop index. 
    int i,j,tt;    
    int k; // k has a meaning in the algorithm 
    int chosen=0;//initialize to get rid of compilation warning
    
    int n=B.rows();
    int p=B.cols()+1;
    
    // for hinge or segmented, compute cusum of B and r in the reverse order
    // for upperhinge, compute cusum of B and r in the forward order
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
    //for (j=0; j<n; j++) PRINTF("%f ", B(j,p-1)); PRINTF("\n");    
    
    Matrix <double,Row,Concrete> VV(2, 2), Vr(2, 1), VB(2, p-1);
    double delta, s, crit, crit_max=R_NegInf;// s is defined as difference in squared threshold
    
    VV(0,0)=0; for (i=0; i<n; i++) VV(0,0) += pow(x[i],2)*w[i]*w[i]; 
    VV(1,1)=0; for (i=0; i<n; i++) VV(1,1) += pow(x[i],4)*w[i]*w[i]; 
    VV(0,1)=0; for (i=0; i<n; i++) VV(0,1) += pow(x[i],3)*w[i]*w[i];
    VV(1,0)=VV(0,1);
    
    Vr(0,0)=0; for (i=0; i<n; i++) Vr(0,0) += x[i]*r(i)*w[i]; 
    Vr(1,0)=0; for (i=0; i<n; i++) Vr(1,0) += pow(x[i],2)*r(i)*w[i]; 
    
    for (j=0; j<p-1; j++) {
        VB(0,j)=0; for (i=0; i<n; i++) VB(0,j) += x[i] * B(i,j)*w[i]; // remember the first p-1 column of B is B after _preprocess()
        VB(1,j)=0; for (i=0; i<n; i++) VB(1,j) += pow(x[i],2) * B(i,j)*w[i]; // remember the first p-1 column of B is B after _preprocess()
    }
    
    // Step 4 (Step 3 = Step 4b)
    for(tt=0; tt<nThresholds; tt++) { 
            
        // Step 4a: update VV, Vr and VB     
        if(tt>0) {
            delta= thresholds[tt]-thresholds[tt-1]; //B(tt,p-1)-B(tt-1,p-1); //PRINTF("%f ", delta); PRINTF("\n"); 
            s=pow(thresholds[tt],2)-pow(thresholds[tt-1], 2);

            // e.g, when tt=1, we are at the second threshold, we want to update vv, vr and vB for e_2
            //            tt=1 and tt+1=2 
            k  = thresholdIdx[tt-1]; // k is 1-based index of x, e_t = x_k
            VV(0,0) += -2*delta*xcusum[k-1] + 2*Wcusum[k-1]*delta*thresholds[tt-1] + Wcusum[k-1] * pow(delta,2); 
            VV(1,1) += -4*delta*x3cusum[k-1] + 2*(2*pow(delta,2)+4*delta*thresholds[tt-1]+s)*x2cusum[k-1] - (4*delta*s+4*delta*pow(thresholds[tt-1],2)+4*thresholds[tt-1]*s)*xcusum[k-1] + Wcusum[k-1]*s*s + 2*Wcusum[k-1]*s*pow(thresholds[tt-1],2); 
            VV(0,1) += -3*delta*x2cusum[k-1] + (2*pow(delta,2)+4*delta*thresholds[tt-1]+s) * xcusum[k-1] - Wcusum[k-1]*delta*s - Wcusum[k-1]*s*thresholds[tt-1] - Wcusum[k-1]*delta*pow(thresholds[tt-1],2); 
            VV(1,0) = VV(0,1);
            Vr[0] += -delta * rcusum[k-1];
            Vr[1] += -2*delta * xrcusum[k-1] + s * rcusum[k-1] ;
            for (j=0; j<p-1; j++) {
                VB(0,j) +=   -delta * Bcusum(k-1,j); 
                VB(1,j) += -2*delta * xBcusum(k-1,j) + s * Bcusum(k-1,j); 
            }
        }        
        
        // Step 4b: compute r'H_eY - r'HY
        if (p>1) {
            crit = ( t(Vr) * chol_solve(VV - VB * t(VB), Vr) ) (0,0);
        } else { // VB has 0 columns
            crit = ( t(Vr) * chol_solve(VV, Vr) ) (0,0);
        }
        //crit = ( t(Vr) * invpd(VV - VB * t(VB)) * Vr ) (0,0);

//        if(tt==0) {
//            PRINTF("VV\n"); 
//            for (i=0; i<2; i++) PRINTF("%f %f %f \n", VV(i,0), VV(i,1)); 
//            PRINTF("Vr %f %f \n", Vr(0,0), Vr(1,0)); 
//            PRINTF("VB\n"); for (i=0; i<2; i++) {for (j=0; j<p-1; j++) PRINTF("%f ", VB(i,j)); PRINTF("\n");}
//            Matrix<> tmp=VV - VB * t(VB);
//            PRINTF("VV - VB * t(VB)\n"); for (i=0; i<2; i++) PRINTF("%f %f \n", tmp(i,0), tmp(i,1)); 
//            Matrix<> tmp2=invpd(VV - VB * t(VB))* Vr;
//            PRINTF("invpd(VV - VB * t(VB))* Vr\n"); PRINTF("%f %f \n", tmp2(0,0), tmp2(1,0)); 
//            Matrix<> tmp3=chol_solve(VV - VB * t(VB), Vr);
//            PRINTF("chol_solve(VV - VB * t(VB), Vr)\n"); PRINTF("%f %f \n", tmp3(0,0), tmp3(1,0)); 
//        }

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
