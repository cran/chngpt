\name{chngpt.logistic}
\alias{chngpt.logistic}
\alias{chngpt.logistic.2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate change point logistic model
}
\description{
Estimate change point logistic model
}

\usage{
chngpt.logistic   (formula.null, formula.chngpt, data, tol = 1e-04, 
    maxit = 1000, verbose=FALSE)
chngpt.logistic.2 (formula.null, formula.chngpt, data, tol = 1e-04, 
    maxit = 1000, verbose=FALSE)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula.null}{
formula. Null model formula.}
  \item{formula.chngpt}{
formula. Change point variable formula.
}
  \item{data}{
data frame.
}
  \item{tol}{
numeric. Stopping criterion on the coefficient estimate.
}
  \item{maxit}{
integer. Maximum number of iterations in the outer loop of optimization.
}
  \item{verbose}{Boolean.}
}
\details{

chngpt.logistic.2 searches over all distinct change points. It only works for main effect for now.

}
\value{
A an object of type  with the following components
  \item{converged}{Boolean}
  \item{coefficients}{vector. Estimated coefficients. The last element, named ".chngpt", is the estimated change point}
  \item{test}{htest. Max score test results}
  \item{iter}{integer. Number of iterations}
}
\examples{

library(RUnit)

seed=1
label="sigmoid4"; x.distr="norm"
beta=log(1)
e.=3.4
alpha=sim.alphas[[paste(label,"_",x.distr,sep="")]][paste(e.,"", sep=""),paste(beta,"",sep="")]
if(is.null(alpha)) stop("alpha not found")
data=sim.sigmoid(label, n=250, seed=seed, alpha, beta, x.distr=x.distr, e.=e., b.=-Inf)

formula.null=y~z; formula.chngpt=~x*z
fit.1 = chngpt.logistic   (formula.null, formula.chngpt, data, tol=1e-4, maxit=1e3) 

formula.null=y~z; formula.chngpt=~x
fit.2 = chngpt.logistic.2 (formula.null, formula.chngpt, data, tol=1e-4, maxit=1e3) 

checkException(
    chngpt.logistic.2 (formula.null=y~z, formula.chngpt=~x*z, data, tol=1e-4, maxit=1e3)
)


}
