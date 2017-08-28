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

chngptm (formula.1, formula.2, family, data,
  type=c("step","hinge","segmented","segmented2","stegmented"), 
  est.method=c("default","smoothapprox","grid"), est.method.boot=c("grid","smoothapprox"),
  var.type=c("none","robust","model","smooth","robusttruth","bootstrap","all"),
  aux.fit=NULL,
  test.inv.ci=TRUE, boot.test.inv.ci=TRUE,
  lb.quantile=.1, ub.quantile=.9, grid.search.max=500, ci.bootstrap.size=500, alpha=0.05, 
  save.boot=FALSE, m.out.of.n=FALSE,
  b.transition=Inf,
  tol=1e-4, maxit=1e2, chngpt.init=NULL, search.bound=10,
  weights=NULL, verbose=FALSE,
  ...) 


\method{coef}{chngptm} (object, ...)
\method{vcov}{chngptm} (object, ...)
\method{print}{chngptm} (x, ...)
\method{summary}{chngptm} (object, var.type=NULL, verbose=FALSE, ...)

}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula.1}{The part of formula that is free of terms involving thresholded variables}
  \item{formula.2}{The part of formula that is only composed of thresholded variables}
  \item{family}{string. coxph or any valid argument that can be passed to glm. But variance estimate is only available for binomial and gaussian (only model-based for latter)}
  \item{data}{data frame.}
  \item{type}{types of threshold effects. segmented2 differs from segmented in parameterization. }
  \item{b.transition}{Numeric. Controls whether threshold model or smooth transition model. Default to Inf, which correponds to threshold model}
  \item{est.method}{string. Better leave it at NULL. grid: grid search; smoothapprox: smooth approximation}
  \item{var.type}{string. Different methods for estimating covariance matrix and constructing confidence intervals}
  \item{aux.fit}{a model fit object that is needed for model-robust estimation of covariance matrix}
  \item{grid.search.max}{integer. }
  \item{test.inv.ci}{Boolean, whether or not to find test-inversion confidence interval for threshold }
  \item{ci.bootstrap.size}{integer, number of bootstrap }
  \item{alpha}{double, norminal type I error rate}
  \item{save.boot}{Boolean, whether to save bootstrap samples}
  \item{lb.quantile}{lower bound of the search range for change point estimate}
  \item{ub.quantile}{upper bound of the search range for change point estimate}
  \item{tol}{Numeric. Stopping criterion on the coefficient estimate.}
  \item{maxit}{integer. Maximum number of iterations in the outer loop of optimization.}
  \item{chngpt.init}{numeric. Initial value for the change point.}
  \item{weights}{passed to glm}
  \item{verbose}{Boolean.}

  \item{x}{chngptm fit object.}
  \item{object}{chngptm fit object.}
  \item{\dots}{arguments passed to glm or coxph}
  \item{m.out.of.n}{whether to perform m out of n bootstrap}

  \item{est.method.boot}{estimation search method for bootstrap}
  \item{boot.test.inv.ci}{whether to get test inversion CI under bootstrap}
  \item{search.bound}{bounds for search for sloping parameters}

}
\details{
Without lb.quantile and ub.quantile, finite sample performance of estimator drops considerably!
When est.method is smoothapprox, Newton-Raphson is done with initial values chosen by change point hypothesis testing. The testing procedure may be less subjective to finite sample volatility. 

}
\value{
A an object of type chngptm with the following components
  \item{converged}{Boolean}
  \item{coefficients}{vector. Estimated coefficients. The last element, named ".chngpt", is the estimated change point}
  \item{test}{htest. Max score test results}
  \item{iter}{integer. Number of iterations}
}
\references{

Fong, Y., Huang, Y., Gilbert, P., Permar S. (2017) Estimation and hypothesis testing under four types of change point models using chngpt, BMC Bioinformatics, under review.

Fong, Y., Di, C., Huang, Y., Gilbert, P. (2017) Model-robust inference for continuous threshold regression models, Biometrics, 73(2):452-462.

Pastor-Barriuso, R. and Guallar, E. and Coresh, J. (2003) Transition models for change-point estimation in logistic regression. Statistics in Medicine. 22:13141

}
\examples{

dat=sim.chngpt("sigmoid4", type="step", n=200, seed=1, beta=1, alpha=-1, x.distr="norm", e.=4)

fit.1=chngptm(formula.1=y~z, formula.2=~x, family="binomial", dat, type="step", est.method="grid")
print(fit.1)
summary(fit.1)

\dontrun{
# not run because otherwise the examples take >5s and that is a problem for R CMD check
# has interaction

library(survival)
test1 <- data.frame(time=c(4,3,1,1,2,2,3), 
              status=c(1,1,1,0,1,1,0), 
              x=c(0,2,1,1,1,0,0), 
              sex=c(0,0,0,0,1,1,1)) 
fit=chngptm(formula.1=Surv(time, status)~1, formula.2=~x, data=test1, family="coxph", type="step")

}

}
