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

chngptm (formula.1, formula.2, family, data, type = c("hinge",
         "quadhinge", "cubichinge", "x4hinge", "upperhinge",
         "quadupperhinge", "cubicupperhinge", "x4upperhinge",
         "M21", "M12", "M22", "M22c", "M31", "M13", "M33c",
         "segmented", "segmented2", "step", "stegmented"),
         formula.strat = NULL, weights = NULL, offset = NULL,
         est.method = c("default", "fastgrid2", "fastgrid", "grid",
         "smoothapprox"), var.type =
         c("default", "none", "robust", "model", "bootstrap",
         "all"), aux.fit = NULL, lb.quantile = 0.1, ub.quantile
         = 0.9, grid.search.max = Inf, test.inv.ci = TRUE,
         boot.test.inv.ci = FALSE, ci.bootstrap.size = 1000,
         alpha = 0.05, save.boot = TRUE, m.out.of.n = 0,
         b.transition = Inf, tol = 1e-04, maxit = 100,
         chngpt.init = NULL, search.bound = 10, keep.best.fit =
         TRUE, ncpus = 1, verbose = FALSE, ...)
         
chngptm.xy(x, y, type=c("step","hinge","segmented","segmented2","stegmented"),
    ...)

\method{coef}{chngptm} (object, ...)
\method{residuals}{chngptm} (object, ...)
\method{vcov}{chngptm} (object, var.type=NULL, ...)
\method{print}{chngptm} (x, ...)
\method{plot}{chngptm} (x, which = NULL, xlim = NULL, lwd = 2, lcol = "red",
                 lty = 1, add = FALSE, add.points = TRUE, add.ci =
                 TRUE, breaks = 20, mark.chngpt = FALSE, xlab = NULL,
                 ylab = NULL, ...) 
\method{summary}{chngptm} (object, var.type = NULL, expo = FALSE,
         show.slope.post.threshold = FALSE, verbose = FALSE,
         ...)

lincomb(object, comb, alpha=0.05)

}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula.1}{The part of formula that is free of terms involving thresholded variables}
  \item{formula.2}{The part of formula that is only composed of thresholded variables}
  \item{formula.strat}{stratification formula}
  \item{family}{string. coxph or any valid argument that can be passed to glm. But variance estimate is only available for binomial and gaussian (only model-based for latter)}
  \item{data}{data frame.}
  \item{type}{types of threshold effects. segmented2 differs from segmented in parameterization. }
  \item{b.transition}{Numeric. Controls whether threshold model or smooth transition model. Default to Inf, which correponds to threshold model}
  \item{est.method}{default: estimation algorithm will be chosen optimally; fastgrid2: a super fast grid search algorithm, limited to linear regression; grid: plain grid search, works for almost all models; smoothapprox: approximates the likelihood function using a smooth function, only works for some models. fastgrid = fastgrid2, kept for backward compatibility}
  \item{var.type}{string. Different methods for estimating covariance matrix and constructing confidence intervals}
  \item{aux.fit}{a model fit object that is needed for model-robust estimation of covariance matrix}
  \item{grid.search.max}{The maximum number of grid points used in grid search. When doing fast grid search, grid.search.max is set to Inf internally because it does not take more time to examine all potential thresholds.}
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
  \item{ncpus}{Number of cores to use if the OS is not Windows.}
  \item{keep.best.fit}{Boolean.}
  \item{y}{outcome}
  \item{show.slope.post.threshold}{boolean}

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
  \item{mark.chngpt}{mark.chngpt}
  \item{xlab}{xlab}
  \item{ylab}{ylab}
  \item{offset}{offset}
  \item{lty}{lty}
  
  \item{comb}{a vector of combination coefficients that will be used to form an inner product with the estimated slope}
  \item{expo}{If family is binomial and expo is TRUE, coefficients summary will be shown on the scale of odds ratio instead of slopes}

}
\details{
Without lb.quantile and ub.quantile, finite sample performance of estimator drops considerably!
When est.method is smoothapprox, Newton-Raphson is done with initial values chosen by change point hypothesis testing. The testing procedure may be less subjective to finite sample volatility. 

If var.method is bootstrap, summary of fitted model contains p values for each estimated slope. These p values are approximate p-values, obtained assuming that the bootstrap distributions are normal. 

When var.method is bootstrap and the OS is not Windows, the boot package we use under the hood takes advantage of \code{ncpus} cores through parallel::mclapply. 

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
types=c("hinge", "segmented", "quadhinge", "cubichinge")
for (type in types) {
    fit=chngptm(formula.1=logratio~1, formula.2=~range, lidar, type=type, family="gaussian", 
        var.type="bootstrap")
    summary(fit)
    plot(fit)
}


# with weights
dat.1=sim.chngpt("thresholded", "segmented", n=200, seed=1, beta=1, alpha=-1, x.distr="norm", e.=4,
    family="gaussian")
fit.1.a=chngptm(formula.1=y~z, formula.2=~x, family="gaussian", dat.1, type="segmented", 
    est.method="fastgrid", var.type="bootstrap", weights=ifelse(dat.1$x<3.5,100,1)
    , ci.bootstrap.size=10)
summary(fit.1.a)
plot(fit.1.a)


# threshold logistic regression
dat.2=sim.chngpt("thresholded", "step", n=200, seed=1, beta=1, alpha=-1, x.distr="norm", e.=4, 
    family="binomial")

fit.2=chngptm(formula.1=y~z, formula.2=~x, family="binomial", dat.2, type="step", est.method="grid")
summary(fit.2) 
# no variance estimates available for discontinuous threshold models such as step
# vcov(fit.2$best.fit) gives the variance estimates for the best model conditional on threshold est

# also supports cbind() formula on left hand side
set.seed(1)
dat.2$success=rbinom(nrow(dat.2), 10, 1/(1 + exp(-dat.2$eta)))
dat.2$failure=10-dat.2$success
fit.2a=chngptm(formula.1=cbind(success,failure)~z, formula.2=~x, family="binomial", dat.2, 
    type="step")


# Poisson example
counts <- c(18,17,15,20,10,20,25,13,12,33,35)
x <- 1:length(counts)
print(d.AD <- data.frame(x, counts))
fit.4=chngptm(formula.1=counts ~ 1, formula.2=~x, data=d.AD, family="poisson", 
    type="segmented", var.type="bootstrap", verbose=1, ci.bootstrap.size=1)
summary(fit.4)




\dontrun{
# Not run because otherwise the examples take >5s and that is a problem for R CMD check

# coxph example
library(survival)
fit=chngptm(formula.1=Surv(time, status) ~ ph.ecog, formula.2=~age, data=lung, family="coxph",
    type="segmented", var.type="bootstrap", ci.bootstrap.size=10)
summary(fit)


# one interaction term (mtcars is part of R default installation)
# est.method will be grid as fastgrid not available for models with interaction terms yet
fit=chngptm(formula.1=mpg ~ hp, formula.2=~hp*drat, mtcars, type="segmented", 
    family="gaussian", var.type="bootstrap", ci.bootstrap.size=10)
summary(fit)

# interaction, upperhinge model, bootstrap
fit=chngptm(formula.1=mpg ~ hp, formula.2=~hp*drat, mtcars, type="upperhinge", 
    family="gaussian", var.type="bootstrap", ci.bootstrap.size=10)
summary(fit)

# more than one interaction term
# subsampling bootstrap confidence interval for step model
fit=chngptm(formula.1=mpg~hp+wt, formula.2=~hp*drat+wt*drat, mtcars, type="step",
    family="gaussian", var.type="bootstrap", ci.bootstrap.size=10)
summary(fit)

# step model, subsampling bootstrap confidence intervals
fit=chngptm(formula.1=mpg~hp, formula.2=~drat, mtcars, type="step",
    family="gaussian", var.type="bootstrap", ci.bootstrap.size=10, verbose=TRUE)
summary(fit)

# higher order threshold models
dat=sim.chngpt(mean.model="thresholded", threshold.type="M22", n=500, seed=1, 
    beta=c(32,2,10, 10), x.distr="norm", e.=6, b.transition=Inf, family="gaussian", 
    alpha=0, sd=0, coef.z=0)
fit.0=chngptm(formula.1=y~z, formula.2=~x, dat, type="M22", family="gaussian", 
    est.method="fastgrid2"); plot(fit.0)

dat=sim.chngpt(mean.model="thresholded", threshold.type="M22c", n=500, seed=1, 
    beta=c(32,2,32, 10), x.distr="norm", e.=6, b.transition=Inf, family="gaussian", 
    alpha=0, sd=0, coef.z=0)
fit.0=chngptm(formula.1=y~z, formula.2=~x, dat, type="M22c", family="gaussian", 
    est.method="fastgrid2"); plot(fit.0)
    


}

}
