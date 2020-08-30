
//////////////////////////////////////////////////////////////////////////
// 
// This software is distributed under the terms of the GNU GENERAL
// PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
// file for more information.
//
//////////////////////////////////////////////////////////////////////////

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


double M10_search(
    //input
    Matrix<double, Row>& B, Matrix<double, Row>& r, vector<double>& x, vector<double>& w, 
    int * thresholdIdx, vector<double>& thresholds, int nThresholds,
    //output
    Matrix<double,Row>& Bcusum, vector<double>& rcusum, vector<double>& Wcusum, vector<double>& xcusum, 
    double* logliks)
{

    // loop index. 
    int i,j,t;    
    int k;
    //int iX, skipped; // k has a meaning in the algorithm 
    int chosen=0;//initialize to get rid of compilation warning
    
    int n=B.rows();
    int p=B.cols()+1;
    
    xcusum[0] = x[0]* w[0]* w[0];    
    rcusum[0] = r(0)* w[0];        
    Wcusum[0] = w[0]* w[0];    
    if(p>1) Bcusum(0,_) = B(0,_)* w[0];    // if we take out the if condition, there will be memory issue when p is 1
    for (i=1; i<n; i++) {
        xcusum[i]           = xcusum[i-1]      + x[i]* w[i]* w[i];     // x is not prescaled, which is why w is multiplied twice
        rcusum[i]           = rcusum[i-1]      + r(i)* w[i];           // r is prescaled
        Wcusum[i]           = Wcusum[i-1]      + w[i]* w[i];           // w is sqrt of true weights
        if(p>1) Bcusum(i,_) = Bcusum(i-1,_)    + B(i,_)* w[i];         // B is prescaled
    }
    //for (i=0; i<n; i++) {for (j=0; j<p; j++)  PRINTF("%f ", Bcusum(i,j)); PRINTF("\n");}        
    //for (i=0; i<n; i++) PRINTF("%f ", rcusum[i]); PRINTF("\n");
    
    vector<double> vB(p-1); 
    double vv, vr;
    double delta, crit, crit_max=R_NegInf;
    
    // Step 2: Initialize vv, vr and vB 
    
    // (x-e)-
    for(i=0; i<thresholdIdx[0]; i++) x[i]=x[i]-thresholds[0];  
    for(i=thresholdIdx[0]; i<n; i++) x[i]=0; 
    //for (j=0; j<n; j++) PRINTF("%f ", x[j]); PRINTF("\n");
    
    vv=0; vr=0; 
    for (i=0; i<n; i++) vv += pow(x[i],2)*w[i]* w[i]; 
    for (i=0; i<n; i++) vr += x[i]*r(i)*w[i]; 
    for (j=0; j<p-1; j++) {
        vB[j]=0; 
        for (i=0; i<n; i++) vB[j] += x[i] * B(i,j)*w[i];
    }
    
    // Step 4 (Step 3 = Step 4b)
    for(t=0; t<nThresholds; t++) { 
            
        // Step 4a: update vv, vr and vB
        if(t>0) {
            delta= thresholds[t]-thresholds[t-1]; //B(t,p-1)-B(t-1,p-1); //PRINTF("%f ", delta); PRINTF("\n"); 
            // e.g, when t=1, we are at the second threshold, we want to update vv, vr and vB for e_2
            //            t=1 and t+1=2 
            k  = thresholdIdx[t-1]; // k is 1-based index of x, e_t = x_k
            vv -= 2*delta*(xcusum[k-1] - Wcusum[k-1]*thresholds[t-1]) - Wcusum[k-1] * pow(delta,2);            
            vr -= delta * rcusum[k-1];
            for (j=0; j<p-1; j++) vB[j] -= delta * Bcusum(k-1,j); 
        }
        
        // Step 4b: compute r'H_eY - r'HY
        // first compute vB'vB, abusing the notation a little bit
        crit=0; for (j=0; j<p-1; j++) crit += pow(vB[j],2); 
        // now compute crit
        crit = vr*vr / (vv - crit);
        logliks[t] = crit;
        if(crit>=crit_max) {
            chosen = t;
            crit_max=crit;
        }             
    }        
    //PRINTF("logliks: \n");  for(i=0; i<nThresholds; i++) PRINTF("%f ", logliks[i]); PRINTF("\n");  

    return thresholds[chosen]; 
    
}



