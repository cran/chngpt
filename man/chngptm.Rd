\name{chngptm}
\alias{chngptm}
\alias{chngptm.xy}
\alias{coef.chngptm}
\alias{residuals.chngptm}
\alias{vcov.chngptm}
\alias{plot.chngptm}
\alias{print.chngptm}
\alias{summary.chngptm}
\alias{lincomb}
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
  # segmented2 is the model studied in Cheng 2008
  weights=NULL, # this arg is kept here due to historical reasons
  est.method=c("default","smoothapprox","grid","fastgrid"), 
  var.type=c("none","robust","model","smooth","robusttruth","bootstrap","all"), 
  aux.fit=NULL, 
  lb.quantile=.1, ub.quantile=.9, grid.search.max=5000, 
  test.inv.ci=TRUE, boot.test.inv.ci=FALSE, 
  # test.inv.ci is passed to local functions, boot.test.inv.ci is global in this function
  ci.bootstrap.size=1000, alpha=0.05, save.boot=FALSE, m.out.of.n=FALSE, 
  # grid.search.max is the maximum number of grid points used in grid search
  b.transition=Inf,# controls whether threshold model or smooth transition model
  tol=1e-4, maxit=1e2, chngpt.init=NULL, search.bound=10,
  keep.best.fit=TRUE, # best.fit is needed for making prediction and plotting
  verbose=FALSE, ...) 

chngptm.xy(x, y, type=c("step","hinge","segmented","segmented2","stegmented"),
    ...)

\method{coef}{chngptm} (object, ...)
\method{residuals}{chngptm} (object, ...)
\method{vcov}{chngptm} (object, var.type=NULL, ...)
\method{print}{chngptm} (x, ...)
\method{plot}{chngptm} (x, which = NULL, xlim = NULL, lwd = 2, lcol = "red",
 add = FALSE, add.points = TRUE, add.ci = TRUE, breaks
 = 20, ...) 
\method{summary}{chngptm} (object, var.type=NULL, expo=FALSE, verbose=FALSE, ...)

lincomb(object, comb, alpha=0.05)

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
  \item{add.points}{Boolean.}
  \item{add.ci}{Boolean.}
  \item{add}{Boolean.}
  \item{breaks}{integer.}
  \item{keep.best.fit}{Boolean.}
  \item{y}{outcome}

  \item{x}{chngptm fit object.}
  \item{object}{chngptm fit object.}
  \item{\dots}{arguments passed to glm or coxph}
  \item{m.out.of.n}{whether to perform m out of n bootstrap}

  \item{boot.test.inv.ci}{whether to get test inversion CI under bootstrap}
  \item{search.bound}{bounds for search for sloping parameters}
  \item{which}{an integer}
  \item{xlim}{xlim}
  \item{lwd}{lwd}
  \item{lcol}{line col}
  
  \item{comb}{a vector of combination coefficients that will be used to form an inner product with the estimated slope}
  \item{expo}{If family is binomial and expo is TRUE, coefficients summary will be shown on the scale of odds ratio instead of slopes}

}
\details{
Without lb.quantile and ub.quantile, finite sample performance of estimator drops considerably!
When est.method is smoothapprox, Newton-Raphson is done with initial values chosen by change point hypothesis testing. The testing procedure may be less subjective to finite sample volatility. 

If var.method is bootstrap, summary of fitted model contains p values for each estimated slope. These p values are approximate p-values, obtained assuming that the bootstrap distributions are normal. 

}
\value{
A an object of type chngptm with the following components
  \item{converged}{Boolean}
  \item{coefficients}{vector. Estimated coefficients. The last element, named ".chngpt", is the estimated change point}
  \item{test}{htest. Max score test results}
  \item{iter}{integer. Number of iterations}
}
\references{

Fong, Y. (2018) Fast bootstrap confidence intervals for continuous threshold linear regression. Journal of Computational and Graphical Statistics, in press.

Fong, Y., Di, C., Huang, Y., Gilbert, P. (2017) Model-robust inference for continuous threshold regression models, Biometrics, 73(2):452-462.

Pastor-Barriuso, R. and Guallar, E. and Coresh, J. (2003) Transition models for change-point estimation in logistic regression. Statistics in Medicine. 22:13141

}
\examples{

    
# threshold linear regression
dat.1=sim.chngpt("thresholded", "segmented", n=200, seed=1, beta=1, alpha=-1, x.distr="norm", e.=4,
    family="gaussian")
# default search method is fastgrid
fit.1=chngptm(formula.1=y~z, formula.2=~x, family="gaussian", dat.1, type="segmented", 
    est.method="fastgrid", var.type="bootstrap", ci.bootstrap.size=10)
summary(fit.1)
head(resid(fit.1))
plot(fit.1)
# with weights
fit.1.a=chngptm(formula.1=y~z, formula.2=~x, family="gaussian", dat.1, type="segmented", 
    est.method="fastgrid", var.type="bootstrap", weights=ifelse(dat.1$x<3.5,100,1)
    , ci.bootstrap.size=10)
summary(fit.1.a)
plot(fit.1.a)


# threshold logistic regression
dat.2=sim.chngpt("thresholded", "step", n=200, seed=1, beta=1, alpha=-1, x.distr="norm", e.=4, 
    family="binomial")

fit.2=chngptm(formula.1=y~z, formula.2=~x, family="binomial", dat.2, type="step",est.method="grid")
summary(fit.2) 
# no variance estimates available for discontinuous threshold models such as step
# vcov(fit.2$best.fit) gives the variance estimates for the best model conditional on threshold est

# also supports cbind() formula on left hand side
set.seed(1)
dat.2$success=rbinom(nrow(dat.2), 10, 1/(1 + exp(-dat.2$eta)))
dat.2$failure=10-dat.2$success
fit.2a=chngptm(formula.1=cbind(success,failure)~z, formula.2=~x, family="binomial", dat.2, 
    type="step")


# threshold Cox regression, no variance estimates available
library(survival)
dat.3 <- data.frame(time=c(4,3,1,1,2,2,3), status=c(1,1,1,0,1,1,0), x=c(0,2,1,1,1,0,0), 
    sex=c(0,0,0,0,1,1,1)) 
fit.3=chngptm(formula.1=Surv(time, status)~1, formula.2=~x, data=dat.3, family="coxph",type="step")
summary(fit.3)


# Poisson
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- as.integer(gl(3,1,9))
treatment <- gl(3,3)
print(d.AD <- data.frame(treatment, outcome, counts))
fit.4=chngptm(formula.1=counts ~ treatment, formula.2=~outcome, data=d.AD, family="poisson", 
    type="segmented", var.type="bootstrap", verbose=1, ci.bootstrap.size=10)
summary(fit.4)




\dontrun{
# not run because otherwise the examples take >5s and that is a problem for R CMD check

# an interaction model. est.method will be set to grid inside the function. fastgrid not
# available for models with interaction terms yet
fit=chngptm(formula.1=mpg~cyl + disp + hp, formula.2=~hp*drat, mtcars, type="segmented", 
    family="gaussian", var.type="bootstrap", ci.bootstrap.size=100)
summary(fit)


}

}
