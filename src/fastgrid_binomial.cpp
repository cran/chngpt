//////////////////////////////////////////////////////////////////////////
// 
// This software is distributed under the terms of the GNU GENERAL
// PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
// file for more information.
//
// Author: Tao Yang
//
//////////////////////////////////////////////////////////////////////////


// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;                   // use the Armadillo library for matrix computations
using namespace Rcpp;

//namespace myf{                                
          // define IRLS for solving beta(e)
  vec beta_e(double e, vec y, mat Z, vec x){
    int n = y.size();
    vec mu = zeros(n);
    vec res = zeros(n);
    double err = 1;
    double tol = 1e-6;
    int k = 0;
    int p = Z.n_cols+1;
    mat H(p,p);
    mat Eps = zeros(p,p);
    Eps.diag() = tol*ones(p);
    mat X(n,p);
    
    X.cols(0,(p-2)) = Z;
    X.col(p-1) = (x-e)%(x>e);
    vec beta0 = zeros(p);
    vec beta1 = zeros(p);
    mat D(n,p);
    vec out = NA_REAL*ones(p);
    while (err>tol){
      k = k+1;
      if (k>200){
        return out;
      }
//    cout << k << endl;
      mu = exp(X*beta0)/(1+exp(X*beta0));
      res = y-mu;
      for (int i=0;i<p;i++){
        D.col(i) = mu%(1-mu);
      }
      H = trans(X)*(D%X);
      if (det(H)<tol){
        H = H+Eps;
      }
      beta1 = inv(H)*trans(X)*res+beta0;
      err = sum(abs(beta1-beta0));
      beta0 = beta1;
    }
    return beta0;
  }
//}

// [[Rcpp::export]]

vec _fastgrid_binomial(vec y, mat Z, vec x, vec e_seq){       // fast grid search for model type hinge/segmented
  int n = y.size();
  vec res(n);
  int n_grid = e_seq.size();
  double e = e_seq[0];
  int p = Z.n_cols+2;
  vec beta0(p-1);
  vec out(p);
  mat X(n,p-1);
  X.cols(0,p-3) = Z;
  X.col(p-2) = (x-e)%(x>e);
  mat D(n,p-1);
  mat H(p-1,p-1);

  mat Beta (n_grid,p-1);
  vec l(n_grid);                               // vector of likelihood
    
  beta0 = beta_e(e,y,Z,x);

  Beta.row(0) = trans(beta0);
  vec mu = exp(X*beta0)/(1+exp(X*beta0));
  l[0] = sum(log((y==1)%mu+(y==0)%(1-mu)));
  double tol=1e-5;
  mat Eps = zeros(p-1,p-1);
  Eps.diag() = tol*ones(p-1);
    
  for (int t=1;t<n_grid;t++){               // use (5) to approximate beta[t]
//  cout << t << endl;
    e = e_seq[t];
    X.col(p-2) = (x-e)%(x>e);
    mu = exp(X*beta0)/(1+exp(X*beta0));
    for (int j=0;j<(p-1);j++){
      D.col(j) = mu%(1-mu);
    }
    res = y-mu;
    H = trans(X)*(D%X);
    if (det(H)<tol){
      H = H+Eps;
    }
    beta0 = inv(H)*trans(X)*res+beta0;
    mu = exp(X*beta0)/(1+exp(X*beta0));
    l[t] = sum(log((y==1)%mu+(y==0)%(1-mu)));
    Beta.row(t) = trans(beta0);
  }
  return l;
    
//  int ind = l.index_max();
//  e = e_seq[ind];                   // identify ehat which achieves maximum likelihood
//  beta0 = myf::beta_e(e,y,Z,x);     // calculate exact betahat
//  out.subvec(0,p-2) = beta0;
//  out[p-1] = e;
//  return out;
}






