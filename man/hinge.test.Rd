\name{hinge.test}
\alias{hinge.test}
\title{
A non-nested hypothesis testing problem for threshold regression models
}
\description{
Test a hinge effect against a linear effect
}
\usage{
hinge.test(formula, cov.interest, family = c("binomial", "gaussian"), data, thres = NA,
    lb.quantile = 0.1, ub.quantile = 0.9, chngpts.cnt = 10, method = c("FDB", "B", "DB"),
    boot.B = 10000, B2 = NA, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
%%     ~~Describe \code{formula} here~~
}
  \item{cov.interest}{
%%     ~~Describe \code{cov.interest} here~~
}
  \item{family}{
%%     ~~Describe \code{family} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{thres}{
If supplied, this will be the threshold value to use in the hinge model.
}
  \item{lb.quantile}{
    lower bound of threshold candidates in quantile
}
  \item{ub.quantile}{
    upper bound of threshold candidates in quantile
}
  \item{chngpts.cnt}{
    number of candidate thresholds 
}
  \item{method}{
    type of test. FDB: false double bootstrap, B: parametric bootstrap, DB: double bootstrap.
}
  \item{boot.B}{
    number of parametric bootstrap replicates for B and FDB
}
  \item{B2}{
    number of inner bootstrap replicates for DB
  }
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
}
\value{
A list of class htest
\item{p.value}{P-value}
\item{chngpts}{Vector of change points evaluated}
\item{TT}{Standardized absolute score statistics}
\item{V.S.hat}{Estimated variance-covariance matrix of the score statistics}
}
\references{
He, Fong, Fouda, Permar. A non-nested hypothesis testing problem for threshold regression model, under review}
\author{
Zonglin He
}
\examples{

dat=sim.hinge(threshold.type = 'NA',family = 'binomial',thres='NA',X.ditr = 'norm',mu.X = c(0,0,0),
    coef.X = c(0,.5,.5,.4),cov.X = diag(3),eps.sd = 1,seed = 1,n=100)
test=hinge.test(Y~X1+X2, "x", family="binomial", data=dat,'method'='FDB',boot.B=10)
test

}
