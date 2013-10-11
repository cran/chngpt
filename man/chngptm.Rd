\name{chngptm}
\alias{chngptm}
\alias{coef.chngptm}
\alias{vcov.chngptm}
\alias{print.chngptm}
\alias{summary.chngptm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate change point logistic model
}
\description{
Estimate change point logistic model
}

\usage{
chngptm   (formula.null, formula.chngpt, data, search.all.thresholds=NULL, tol = 1e-04,
    maxit = 100, chngpt.init=NULL, verbose=FALSE)

\method{coef}{chngptm} (object, ...)
\method{vcov}{chngptm} (object, ...)
\method{print}{chngptm} (x, ...)
\method{summary}{chngptm} (object, ...)

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
  \item{search.all.thresholds}{
Boolean. If NULL, then when sample size is less than 500, it is set to TRUE. If TRUE, search all thresholds. 
}
  \item{tol}{
numeric. Stopping criterion on the coefficient estimate.
}
  \item{maxit}{
integer. Maximum number of iterations in the outer loop of optimization.
}
  \item{chngpt.init}{
numeric. Initial value for the change point.
}
  \item{verbose}{Boolean.}

  \item{x}{chngptm fit object.}
  \item{object}{chngptm fit object.}
  \item{\dots}{arguments passed to or from methods}


}
\details{
Please note that the performance of the estimated variance of the coefficient estimates have not been studied thoroughly.
}
\value{
A an object of type chngptm with the following components
  \item{converged}{Boolean}
  \item{coefficients}{vector. Estimated coefficients. The last element, named ".chngpt", is the estimated change point}
  \item{test}{htest. Max score test results}
  \item{iter}{integer. Number of iterations}
}
\examples{

data=sim.sigmoid("sigmoid4", n=250, seed=1, alpha=sim.alphas[["sigmoid4_norm"]]["3.4","0"], 
    beta=0, x.distr="norm", e.=3.4, b.=-Inf)

fit.1 = chngptm (formula.null=y~z, formula.chngpt=~x, data, tol=1e-4, maxit=1e3) 
print(fit.1)
summary(fit.1)

\dontrun{
# not run because otherwise the examples take >5s and that is a problem for R CMD check
# has interaction
fit.2 = chngptm (formula.null=y~z, formula.chngpt=~x*z, data, tol=1e-4, maxit=1e3
    , search.all.thresholds=TRUE) 
print(fit.2)
summary(fit.2)
}

}
