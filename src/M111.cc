
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



double get_crit(Matrix<double,Row>& VV, Matrix<double,Row>& Vr, Matrix<double,Row>& VB, int p, double crit_max, int verbose) {
       
    Matrix<> tmp;// tmp is symmetric
    if (p>1) {
        tmp=VV - VB * t(VB);
    } else {
        tmp=VV;
    }

    double det=tmp(0,0)*tmp(1,1)-tmp(0,1)*tmp(0,1);
    /* there are other, perhaps numerically more stable ways to compute determinant, e.g.  */
    //det_mean = (tmp(0,0) + tmp(1,1) + tmp(0,1) + tmp(1,0))/4;
    //det=(tmp(0,0)-det_mean)*(tmp(1,1)-det_mean) - (tmp(0,1)-det_mean)*(tmp(0,1)-det_mean) + det_mean * (tmp(0,0) + tmp(1,1) - 2*tmp(0,1));
    /* but for numerical stability, this does not seem the key step. The key seems to be in how det is used in calculating crit */
    
    double crit = tmp(1,1)*Vr(0,0)*Vr(0,0)/det + tmp(0,0)*Vr(1,0)*Vr(1,0)/det - 2*tmp(0,1)*Vr(0,0)*Vr(1,0)/det;
    // by distributing det into the paranthesis, numerical stability is much better when det is very small
    // the usual way to compute crit is not numerically stable and may get strange results when det is very small
    //crit = (tmp(1,1)*Vr(0,0)*Vr(0,0) + tmp(0,0)*Vr(1,0)*Vr(1,0) - 2*tmp(0,1)*Vr(0,0)*Vr(1,0))/det;
    // log and exp do not help either
    //crit = exp(log(tmp(1,1)*Vr(0,0)*Vr(0,0) + tmp(0,0)*Vr(1,0)*Vr(1,0) - 2*tmp(0,1)*Vr(0,0)*Vr(1,0)) - log(det));

//    //when the two threshold values are close to each other, the design matrix is highly collinear and the matrix inversion is not stable. 
//    //This hack throws away a grid point if the criterion function computed at that point is 20% larger than the previous point on the column. 
//    //But the 1.2 cutoff is arbitrary, it fails when dataset is small
//    if (crit>crit_max*1.2) {
//        crit = R_NegInf;
//        PRINTF("\n######## Warning from M111_search: potential numerical precision problem, usually won't affect the result ########\n\n");                                      
//        if(verbose) {
//            PRINTF("det %.20f\n", det);                                         
//        }
//    } 
    
    return (crit);
}


double M111_search(
    Matrix<double,Row>& B, Matrix<double,Row>& r, vector<double>& x, vector<double>& w, 
    int * thresholdIdx, vector<double>& thresholds, int nThresholds,
    // forward cumulative sum
    Matrix<double,Row>& Bcusum, vector<double>& rcusum, vector<double>& Wcusum, vector<double>& xcusum, 
    vector<double>& x2cusum, //vector<double>& x2crosscusum, vector<double>& xrcusum, Matrix<double,Row>& xBcusum, 
    double* ehatvec, int verbose)
{    
    int i,j,tt,s; // loop index, t cannot be used as index because t is also transpose
    int k,g; // k has a meaning in the algorithm 
    double x_e=0; // a place to hold repeatedly used values such as x_cpy[k]-thresholds[tt]
    
    double rel_tol=0.000001; // use this to decide whether two thresholds are the same, useful when doing bootstrap
    
    int chosen1, chosen2;
    
    int n=B.rows();
    int p=B.cols()+1;
    
    if(p>1) Bcusum(0,_) = B(0,_)* w[0];    
    rcusum[0] = r(0)* w[0];        
    Wcusum[0] = w[0]* w[0];              
    xcusum[0]    = x[0] * w[0]* w[0]; 
    x2cusum[0]   = pow(x[0],2) * w[0]* w[0];        
//    x3cusum[0]   = pow(x[0],3) * w[0]* w[0];        
//    xrcusum[0]   = x[0] * r(0)* w[0];        
//    if(p>1) xBcusum(0,_) = x[0] * B(0,_)* w[0];    
    for (i=1; i<n; i++) {
        if(p>1) Bcusum(i,_)  = Bcusum(i-1,_)    + B(i,_)* w[i]; 
        rcusum[i]    = rcusum[i-1]      + r(i)* w[i]; 
        Wcusum[i]    = Wcusum[i-1]      + w[i]* w[i]; 
        xcusum[i]    = xcusum[i-1]      + x[i]* w[i]* w[i]; 
        x2cusum[i]   = x2cusum[i-1]     + pow(x[i],2)* w[i]* w[i]; 
//        x3cusum[i]   = x3cusum[i-1]     + pow(x[i],3)* w[i]* w[i]; 
//        xrcusum[i]   = xrcusum[i-1]     + x[i]* r(i)* w[i]; 
//        if(p>1) xBcusum(i,_) = xBcusum(i-1,_)   + x[i]* B(i,_)* w[i]; 
    }    
    //for (i=0; i<n; i++) {for (j=0; j<p; j++)  PRINTF("%f ", Bcusum(i,j)); PRINTF("\n");}        
    //for (i=0; i<n; i++) PRINTF("%f ", rcusum[i]); PRINTF("\n");
        

    // Step 2: Initialize VV, Vr and VB     
    
    vector<double> x2=x; // needed to initialize VV2 
    vector<double> x_cpy=x; // needed during update
    
    // compute (x-e)- and (x2-e)- at the first and second threshold, respectively
    for(i=0; i<thresholdIdx[0]; i++) x[i]=x[i]-thresholds[0];  
    for(i=thresholdIdx[0]; i<n; i++) x[i]=0; 
    for(i=0; i<thresholdIdx[1]; i++) x2[i]=x2[i]-thresholds[1];  
    for(i=thresholdIdx[1]; i<n; i++) x[i]=0; 
    
    Matrix <double,Row,Concrete> VV(2, 2), Vr(2, 1), VB(2, p-1); // for updating columns
    Matrix <double,Row,Concrete> VV2(2, 2), Vr2(2, 1), VB2(2, p-1); // for updating 1st row
    double delta, crit, crit_max=R_NegInf;// s is defined as difference in squared threshold
    
    VV2(0,0)=0; for (i=0; i<thresholdIdx[0]; i++) VV2(0,0) += x[i]*x[i]*w[i]*w[i]; 
    //VV2(1,1)=0; for (i=0; i<thresholdIdx[1]; i++) VV2(1,1) += x2[i]*x2[i]*w[i]*w[i]; PRINTF("%f\n", VV2(1,1));
    VV2(1,1) = x2cusum[thresholdIdx[1]-1] + thresholdIdx[1]*pow(thresholds[1],2) - 2*thresholds[1]*xcusum[thresholdIdx[1]-1]; //PRINTF("%f\n", VV2(1,1));
    VV2(0,1)=0; for (i=0; i<thresholdIdx[0]; i++) VV2(0,1) += x[i]*x2[i]*w[i]*w[i]; 
    VV2(1,0)=VV2(0,1);
    
    Vr2(0,0)=0; for (i=0; i<thresholdIdx[0]; i++) Vr2(0,0) += x[i] *r(i)*w[i]; 
    Vr2(1,0)=0; for (i=0; i<thresholdIdx[1]; i++) Vr2(1,0) += x2[i]*r(i)*w[i]; 
    
    for (j=0; j<p-1; j++) {
        VB2(0,j)=0; for (i=0; i<thresholdIdx[0]; i++) VB2(0,j) += x[i] * B(i,j)*w[i]; 
        VB2(1,j)=0; for (i=0; i<thresholdIdx[1]; i++) VB2(1,j) += x2[i]* B(i,j)*w[i]; 
    }
                
//    PRINTF("VV2:\n"); for (i=0; i<2; i++) {for (j=0; j<2; j++) PRINTF("%f ", VV2(i,j)); PRINTF("\n");}
//    PRINTF("Vr2:\n"); for (i=0; i<2; i++) PRINTF("%f ", Vr2[i]); PRINTF("\n");
//    PRINTF("VB2:\n"); for (i=0; i<2; i++) {for (j=0; j<p-1; j++) PRINTF("%f ", VB2(i,j)); PRINTF("\n");}
//    Matrix<> tmp=VV2 - VB2 * t(VB2);
//    PRINTF("VV - VB * t(VB)\n"); for (i=0; i<2; i++) PRINTF("%f %f \n", tmp(i,0), tmp(i,1)); 
//    Matrix<> tmp2=invpd(VV2 - VB2 * t(VB2))* Vr2;
//    PRINTF("invpd(VV - VB * t(VB))* Vr\n"); PRINTF("%f %f \n", tmp2(0,0), tmp2(1,0)); 
//    Matrix<> tmp3=chol_solve(VV2 - VB2 * t(VB2), Vr2);
//    PRINTF("chol_solve(VV - VB * t(VB), Vr)\n"); PRINTF("%f %f \n", tmp3(0,0), tmp3(1,0)); 


    // Step 4: Compute criterion function
    // We do not save the matrix of logLik because it takes O(n^2) memory
    // if two thresholds are identical in bootstrap, skip
    if (fabs((thresholds[1]-thresholds[0])/thresholds[0])<rel_tol) {
        crit=R_NegInf; 
    } else {
        crit = get_crit(VV2, Vr2, VB2, p, R_NegInf, verbose);
        if (verbose>=2) PRINTF("0 1 %.10f %.10f %f\n", thresholds[1], thresholds[0], crit);                                         
    }
    crit_max=crit;
    chosen1=0; chosen2=1; 
    //if(verbose) PRINTF("crit %.0f", crit);         
            
    // Step 5:
    if (nThresholds<=2) return 0; //nothing to do if there are only two thresholds since we will have one for e and one for f
    
    // work on the upper triangle of the thresholds matrix, going from left to right
    for(s=2; s<nThresholds; s++) { 

        // work on the first row
        
        delta= thresholds[s]-thresholds[s-1]; 
        k = thresholdIdx[0]; // k is 1-based index of x, e_t = x_k
        g = thresholdIdx[s-1];
        x_e = x_cpy[g]-thresholds[s];
                    
        VV2(0,1) += -delta * (xcusum[k-1] - Wcusum[k-1]*thresholds[0]);
        VV2(1,0) = VV2(0,1);        
        // two ways to compute this        
        //VV2(1,1) += Wcusum[g-1]*delta*delta + 2*delta*(Wcusum[g-1]*thresholds[s-1] - xcusum[g-1]) + x_e*x_e;        
        VV2(1,1) = x2cusum[g-1] + Wcusum[g-1]*pow(thresholds[s],2) - 2*thresholds[s]*xcusum[g-1];
        
        Vr2(1,0) += x_e*r[g] - delta*rcusum[g-1];
        for (j=0; j<p-1; j++) VB2(1,j) += x_e*B(g,j) - delta*Bcusum(g-1,j); 
        
        if (fabs((thresholds[s]-thresholds[0])/thresholds[0])<rel_tol) {
            crit=R_NegInf; // if two thresholds are identical in bootstrap
        } else {
            crit = get_crit(VV2, Vr2, VB2, p, crit_max, verbose);
            if (verbose>=2) PRINTF("0 %d %.10f %.10f %f\n", s, thresholds[s], thresholds[0], crit);                                         
        }    
        if(crit>=crit_max) {
            chosen1=0; chosen2 = s; 
            crit_max=crit;
        }    
                
        // work down the column        
        VV=VV2; Vr=Vr2; VB=VB2; // operate on this copy and leave VV2 for row updates
        for(tt=1; tt<s; tt++) { 
            
            if(verbose>=3) PRINTF("%d %d %f %f\n", tt, s, thresholds[tt], thresholds[s]);  
            delta= thresholds[tt]-thresholds[tt-1]; 
            k  = thresholdIdx[tt-1]; // k is 1-based index of x, e_t = x_k
            x_e = x_cpy[k]-thresholds[tt];
            
            VV(0,1) += x_e*(x_cpy[k]-thresholds[s]) - delta*(xcusum[k-1] - Wcusum[k-1]*thresholds[s]);
            VV(1,0) = VV(0,1);
            
            //gets same results between windows and linux, but results not good
            //VV(0,0) += Wcusum[k-1]*delta*delta + 2*delta*(Wcusum[k-1]*thresholds[tt-1] - xcusum[k-1]) + x_e*x_e; 
            //gets close but diff results between windows and linux, but good results, perhaps VV(1,1) and VV(0,0) need to be computed differently for det to be good
            //VV(0,0) += Wcusum[k-1]*thresholds[tt]*thresholds[tt]- Wcusum[k-1]*thresholds[tt-1]*thresholds[tt-1] - 2*delta*xcusum[k-1] + x_e*x_e; 
            // use x2cusum
            VV(0,0) = x2cusum[k-1] + Wcusum[k-1]*pow(thresholds[tt],2) - 2*thresholds[tt]*xcusum[k-1];
            
            Vr(0,0) += x_e*r[k] - delta*rcusum[k-1];
            for (j=0; j<p-1; j++) VB(0,j) += x_e*B(k,j) - delta*Bcusum(k-1,j);
                                    
            if (fabs((thresholds[s]-thresholds[tt])/thresholds[tt])<rel_tol) {
                crit=R_NegInf; // if two thresholds are identical in bootstrap
            } else {
                crit = get_crit(VV, Vr, VB, p, crit_max, verbose);
                if (verbose>=2) PRINTF("%d %d %.10f %.10f %f\n", tt, s, thresholds[tt], thresholds[s], crit);                                         
            }    
            if(crit>=crit_max) {            
                chosen1 = tt; chosen2 = s; 
                crit_max=crit;
            } 
                
        } // end loop for column update
        //if(verbose) PRINTF("\n"); // print end of line after the crit        

    }
    //PRINTF("logliks: \n");  for(i=0; i<nThresholds; i++) PRINTF("%f ", logliks[i]); PRINTF("\n");  

    *ehatvec=thresholds[chosen1];
    *(ehatvec+1)=thresholds[chosen2];

    return 0; 
    
}



#endif
//#endif
