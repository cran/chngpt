
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

double Mstep_search(
    Matrix<double,Row>& B, 
    Matrix<double,Row>& r, 
    vector<double>& x,
    vector<double>& w, 
    int * thresholdIdx, vector<double>& thresholds, int nThresholds,
    Matrix<double,Row>& Bcusum, vector<double>& rcusum, vector<double>& Wcusum, 
    double* logliks)
{

    // loop index. 
    int i,j,t;    
    int iX; 
    int chosen=0;//initialize to get rid of compilation warning    
    int n=B.rows();
    int p=B.cols()+1;
    
    if(p>1) Bcusum(n-1,_) = B(n-1,_)*w[n-1]; 
    rcusum[n-1] = r(n-1)*w[n-1];     
    Wcusum[n-1] = w[n-1]*w[n-1];             
    for (i=n-2; i>=0; i--) {
        if(p>1) Bcusum(i,_)    = Bcusum(i+1,_)    + B(i,_)* w[i]; 
        rcusum[i]      = rcusum[i+1]      + r(i)* w[i]; 
        Wcusum[i]      = Wcusum[i+1]      + w[i]* w[i]; 
    }
    //for (i=0; i<n; i++) {for (j=0; j<p; j++)  PRINTF("%f ", Bcusum(i,j)); PRINTF("\n");}        
    //for (i=0; i<n; i++) PRINTF("%f ", rcusum[i]); PRINTF("\n");
    
    vector<double> vB(p-1); 
    double vv, vr;
    double crit, crit_max=R_NegInf;
    
    // Step 2: Initialize vv, vr and vB 
    
    // threshold x, which is the last column of B, at the first threshold value
    // (x-e)+
    for(i=0; i<thresholdIdx[0]; i++) x[i]=0;  
    for(i=thresholdIdx[0]; i<n; i++) x[i]=1; 
    //for (j=0; j<n; j++) PRINTF("%f ", B(j,p-1)); PRINTF("\n");
    
    vv=0; vr=0; 
    for (i=0; i<n; i++) vv += pow(x[i],2)* w[i]* w[i]; 
    for (i=0; i<n; i++) vr += x[i]*r(i)* w[i]; 
    for (j=0; j<p-1; j++) {
        vB[j]=0; 
        for (i=0; i<n; i++) 
            vB[j] += x[i] * B(i,j)* w[i];
    }
    
    // Step 4 (Step 3 = Step 4b)
    for(t=0; t<nThresholds; t++) { 
            
        // Step 4a: update vv, vr and vB
        if(t>0) {
            iX = thresholdIdx[t]-1; 
            //k  = n-thresholdIdx[t]+1; // the current threshold is the kth largest             
            vv -= w[iX]*w[iX]; 
            vr -= r[iX]*w[iX];
            for (j=0; j<p-1; j++) vB[j] -= B(iX,j)*w[iX];             
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

#endif
//#endif
