#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;                   

#include <R.h>           // needed to use unif_rand(), Rprintf()
vec beta_e(double e, vec y, mat Z, vec x);
vec _fastgrid_binomial(vec y, mat Z, vec x, vec e_seq);

// copied from random.c, sample with replacement
static void SampleReplace(int k, int n, int *y) {for (int i = 0; i < k; i++) y[i] = n * unif_rand() + 1;}

// for sorting
static int compare_function(const void *a,const void *b) {return *((int *) a) - *((int *) b);}


RcppExport SEXP fastgrid_binomial(
           SEXP u_X, SEXP u_Y, SEXP u_W, 
           SEXP u_wAllOne, SEXP u_thresholdIdx, SEXP u_nBoot, 
           SEXP u_isUpperHinge) {
BEGIN_RCPP
    
    /*    
    double* W_dat=REAL(u_W);    
    bool wAllOne=asLogical(u_wAllOne)==1;
    bool isUpperHinge=asLogical(u_isUpperHinge)==1;    
    */
    
    int nBoot = Rcpp::as<int >(u_nBoot);
    
    // get Z and x and y
    Rcpp::NumericMatrix uX_dat(u_X); // creates Rcpp matrix from SEXP
    int n = uX_dat.nrow(), p = uX_dat.ncol();
    arma::mat Z(uX_dat.begin(), n, p-1, false); // reuses memory and avoids extra copy
    arma::colvec x(uX_dat.begin()+n*(p-1), n, false);
    Rcpp::NumericVector Y_dat(u_Y); // creates Rcpp vector from SEXP
    arma::colvec y(Y_dat.begin(), Y_dat.size(), false);
    
    // define e_seq 
    Rcpp::NumericVector thresholdIdx(u_thresholdIdx); // creates Rcpp vector from SEXP
    //Rcout << thresholdIdx << std::endl; //all are printed on the same line
    int nThresholds = thresholdIdx.size();
    arma::colvec e_seq(nThresholds, fill::zeros);
    
    
	if (nBoot<0.1) {
    // a single search
        Rcpp::RObject logLiks;
        for (int i=0; i<nThresholds; i++) e_seq[i] = x[thresholdIdx[i]-1];
        //Rcout << e_seq << std::endl;// each element is printed on a new line
        logLiks = Rcpp::wrap( _fastgrid_binomial(y, Z, x, e_seq) );
        return logLiks;

    } else {
    // bootstrap
        arma::mat coefs(p+1, nBoot, fill::zeros);
        arma::colvec yb(n, fill::zeros);
        arma::colvec xb(n, fill::zeros);
        arma::mat Zb(n, p-1, fill::zeros);
       	int* index=new int [n];
       	arma::colvec logLiks;
       	double e;
        vec beta0(p);
    
        for (int b=0; b<nBoot; b++) {                    
            // create bootstrap dataset, note that index is 1-based
            SampleReplace(n, n, &(index[0]));
            // Step 1: sort
            qsort (index, n, sizeof(*index), compare_function);
            
            for (int i=0; i<n; i++) { //note that index need to -1 to become 0-based
                Zb.row(i)=Z.row(index[i]-1);
                yb[i]=y[index[i]-1]; 
                xb[i]=x[index[i]-1]; 
                //Wb[i]=W_dat[index[i]-1]; 
            } 
            for (int i=0; i<nThresholds; i++) e_seq[i] = xb[thresholdIdx[i]-1];
            
            logLiks = _fastgrid_binomial(yb, Zb, xb, e_seq);
            int ind = logLiks.index_max();
            e = e_seq[ind];                   // identify ehat which achieves maximum likelihood
            beta0 = beta_e(e,yb,Zb,xb);     // calculate exact betahat
            coefs.submat(0,b,p-1,b) = beta0;
            coefs(p,b) = e;        
        } 
        delete index;
        return Rcpp::wrap(coefs);
    }
END_RCPP
}

