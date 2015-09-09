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
chngptm (formula.1, formula.2, data, family,
  type=c("step","hinge","segmented","stegmented"),
  est.method=NULL, 
  lb.quantile=.1, ub.quantile=.9, grid.size=500,
  tol=1e-4, maxit=1e2, chngpt.init=NULL, 
  weights=NULL,
  verbose=FALSE,
  ...) 

\method{coef}{chngptm} (object, ...)
\method{vcov}{chngptm} (object, ...)
\method{print}{chngptm} (x, ...)
\method{summary}{chngptm} (object, ...)

}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula.1}{
The part of formula that is free of terms involving thresholded variables}
  \item{formula.2}{
The part of formula that is only composed of thresholded variables}
  \item{data}{
data frame.
}
  \item{type}{step: flat before and after change point; hinge: flat before and slope after change point; segmented: slope before and after change point}
  \item{est.method}{string. Better leave it at NULL.}
  \item{grid.size}{integer. }
  \item{lb.quantile}{lower bound of the search range for change point estimate}
  \item{ub.quantile}{upper bound of the search range for change point estimate}
  \item{tol}{
numeric. Stopping criterion on the coefficient estimate.
}
  \item{maxit}{
integer. Maximum number of iterations in the outer loop of optimization.
}
  \item{chngpt.init}{
numeric. Initial value for the change point.
}

  \item{family}{passed to glm}
  \item{weights}{passed to glm}
  \item{verbose}{Boolean.}

  \item{x}{chngptm fit object.}
  \item{object}{chngptm fit object.}
  \item{\dots}{arguments passed to glm or coxph}


}
\details{
Without lb.quantile and ub.quantile, finite sample performance of estimator drops considerably!
When est.method is sigmoidapprox, Newton-Raphson is done with initial values chosen by change point hypothesis testing. The testing procedure may be less subjective to finite sample volatility. 

}
\value{
A an object of type chngptm with the following components
  \item{converged}{Boolean}
  \item{coefficients}{vector. Estimated coefficients. The last element, named ".chngpt", is the estimated change point}
  \item{test}{htest. Max score test results}
  \item{iter}{integer. Number of iterations}
}
\examples{

dat=sim.chngpt("sigmoid4", type="step", n=200, seed=1, 
    beta=0, x.distr="norm", e.=3.4, b.=-Inf)

fit.1 = chngptm (formula.1=y~z, formula.2=~x, family="binomial", dat, type="step") 
print(fit.1)
summary(fit.1)

\dontrun{
# not run because otherwise the examples take >5s and that is a problem for R CMD check
# has interaction
fit.2 = chngptm (formula.1=y~z, formula.2=~x*z, family="binomial", dat, tol=1e-4, maxit=1e3) 
print(fit.2)
summary(fit.2)

fit.1 = chngptm (formula.1=y~z, formula.2=~x, family="binomial", dat, type="hinge") 
print(fit.1)
summary(fit.1)

library(survival)
test1 <- data.frame(time=c(4,3,1,1,2,2,3), 
              status=c(1,1,1,0,1,1,0), 
              x=c(0,2,1,1,1,0,0), 
              sex=c(0,0,0,0,1,1,1)) 
fit=chngptm(formula.1=Surv(time, status)~1, formula.2=~x, data=test1, family="coxph", type="step")

}

}
