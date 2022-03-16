\name{chngptm}
\alias{chngptm}
\alias{chngptm.xy}
\alias{chngptm.xy}
\alias{chngptm.xy}
\alias{coef.chngptm}
\alias{residuals.chngptm}
\alias{vcov.chngptm}
\alias{plot.chngptm}
\alias{predict.chngptm}
\alias{print.chngptm}
\alias{summary.chngptm}
\alias{logLik.chngptm}
\alias{AIC.chngptm}
\alias{lincomb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Threshold Models Estimation
}
\description{
Estimate threshold generalized linear models, Cox proportional hazards models, and linear mixed models. Supports 14 types of two-phase (one threshold) models and 1 type of three-phase (two thresholds) model.
}

\usage{

chngptm (formula.1, formula.2, family, data, type = c("hinge",
 "M01", "M02", "M03", "M04", "upperhinge", "M10",
 "M20", "M30", "M40", "M21", "M12", "M21c", "M12c",
 "M22", "M22c", "M31", "M13", "M33c", "segmented",
 "M11", "segmented2", "M111", "step", "stegmented"),
 formula.strat = NULL, weights = NULL, offset = NULL,
 REML = TRUE, re.choose.by.loglik = FALSE, est.method =
 c("default", "fastgrid2", "fastgrid", "grid",
 "smoothapprox"), var.type = c("default", "none",
 "robust", "model", "bootstrap", "all"), aux.fit =
 NULL, lb.quantile = 0.05, ub.quantile = 0.95,
 grid.search.max = Inf, test.inv.ci = TRUE,
 boot.test.inv.ci = FALSE, bootstrap.type =
 c("nonparametric", "wild", "sieve", "wildsieve",
 "awb"), m.out.of.n = 0, subsampling = 0, order.max =
 10, ci.bootstrap.size = 1000, alpha = 0.05, save.boot
 = TRUE, b.transition = Inf, tol = 1e-04, maxit = 100,
 chngpt.init = NULL, search.bound = 10, keep.best.fit =
 TRUE, ncpus = 1, verbose = FALSE, ...)
         
chngptm.xy(x, y, type=c("step","hinge","segmented","segmented2","stegmented"),
    ...)

\method{coef}{chngptm} (object, ...)
\method{residuals}{chngptm} (object, ...)
\method{vcov}{chngptm} (object, var.type=NULL, ...)
\method{print}{chngptm} (x, ...)
\method{predict}{chngptm} (object, newdata = NULL, 
 type = c("link", "response", "terms"), ...)
\method{plot}{chngptm} (x, which = NULL, xlim = NULL, ylim = NULL, lwd = 2,
         lcol = "red", lty = 1, add = FALSE, add.points = TRUE,
         add.ci = TRUE, breaks = 20, mark.chngpt = TRUE, xlab =
         NULL, ylab = NULL, plot.individual.line = FALSE, main
         = "", y.adj = NULL, auto.adj.y = FALSE, transform =
         NULL, ...)
\method{summary}{chngptm} (object, var.type = NULL, expo = FALSE,
 show.slope.post.threshold = FALSE, verbose = FALSE,
 boot.type = "perc", ...)
\method{logLik}{chngptm} (object, ...)
\method{AIC}{chngptm} (object, ...)
 

lincomb(object, comb, alpha = 0.05, boot.type = "perc")

}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula.1}{The part of formula that is free of terms involving thresholded variables}
  \item{formula.2}{The part of formula that is only composed of thresholded variables}
  \item{formula.strat}{stratification formula}
  \item{family}{string. coxph or any valid argument that can be passed to glm. But variance estimate is only available for binomial and gaussian (only model-based for latter)}
  \item{data}{data frame.}
  \item{type}{type}
  \item{transform}{transform}
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
  \item{newdata}{newdata}
  \item{object}{chngptm fit object.}
  \item{\dots}{arguments passed to glm or coxph}
  \item{m.out.of.n}{sample size for m-out-of-n bootstrap, default 0 for not doing this type of bootstrap}
  \item{subsampling}{sample size for subsampling bootstrap, default 0 for not doing this type of bootstrap}

  \item{boot.test.inv.ci}{whether to get test inversion CI under bootstrap}
  \item{search.bound}{bounds for search for sloping parameters}
  \item{which}{an integer}
  \item{y.adj}{y.adj}
  \item{auto.adj.y}{auto.adj.y}
  \item{xlim}{xlim}
  \item{ylim}{ylim}
  \item{lwd}{lwd}
  \item{lcol}{line col}
  \item{mark.chngpt}{mark.chngpt}
  \item{xlab}{xlab}
  \item{ylab}{ylab}
  \item{offset}{offset}
  \item{lty}{lty}
  \item{boot.type}{lty}
  \item{bootstrap.type}{nonparametric: the default, classical Efron bootstrap, works for homoscedastic and heteroscedastic indepdendent errors; 
                        sieve: works for homoscedastic autocorrelated errors; 
                        wild: works for heteroscedastic independent errors; 
                        wildsieve: works for heteroscedastic autocorrelated errors; 
                        awb: autoregressive wild bootstrap, also works for heteroscedastic autocorrelated errors, but performance may not be as good as wildsieve}
  \item{order.max}{order of autocorrelation for autocorrelated errors in sieve and wildsieve bootstrap}
  
  \item{comb}{a vector of combination coefficients that will be used to form an inner product with the estimated slope}
  \item{expo}{If family is binomial and expo is TRUE, coefficients summary will be shown on the scale of odds ratio instead of slopes}

  \item{REML}{mixed model fitting - should the estimates be chosen to optimize the REML criterion for a fixed threshold}
  \item{re.choose.by.loglik}{mixed model fitting - should the estimates be chosen to optimize likelihood (REML nor not) or goodness of fit}
  \item{plot.individual.line}{boolean}
  \item{main}{character string}
}
\details{
Without lb.quantile and ub.quantile, finite sample performance of estimator drops considerably!\cr
When est.method is smoothapprox, Newton-Raphson is done with initial values chosen by change point hypothesis testing. The testing procedure may be less subjective to finite sample volatility. \cr
\cr
If var.method is bootstrap, summary of fitted model contains p values for each estimated slope. These p values are approximate p-values, obtained assuming that the bootstrap distributions are normal. \cr
\cr
When var.method is bootstrap and the OS is not Windows, the boot package we use under the hood takes advantage of \code{ncpus} cores through parallel::mclapply. \cr
\cr
lincomb can be used to get the estimate and CI for a linear combination of slopes. 
}
\value{
A an object of type chngptm with the following components
  \item{converged}{Boolean}
  \item{coefficients}{vector. Estimated coefficients. The last element, named ".chngpt", is the estimated change point}
  \item{test}{htest. Max score test results}
  \item{iter}{integer. Number of iterations}
}
\references{

Son, H, Fong, Y. (2020) Fast Grid Search and Bootstrap-based Inference for Continuous Two-phase Polynomial Regression Models, Environmetrics, in press.

Elder, A., Fong, Y. (2020) Estimation and Inference for Upper Hinge Regression Models, Environmental and Ecological Statistics, 26(4):287-302.

Fong, Y. (2019) Fast bootstrap confidence intervals for continuous threshold linear regression, Journal of Computational and Graphical Statistics, 28(2):466-470.

Fong, Y., Huang, Y., Gilbert, P., Permar S. (2017) chngpt: threshold regression model estimation and inference, BMC Bioinformatics, 18(1):454.

Fong, Y., Di, C., Huang, Y., Gilbert, P. (2017) Model-robust inference for continuous threshold regression models, Biometrics, 73(2):452-462.

Pastor-Barriuso, R. and Guallar, E. and Coresh, J. (2003) Transition models for change-point estimation in logistic regression. Statistics in Medicine. 22:13141


}
\examples{

# also see the vignette for examples
    
# threshold linear regression
# for actual use, set ci.bootstrap.size to default or higher
par(mfrow=c(2,2))
types=c("hinge", "segmented", "M02", "M03")
for (type in types) {
    fit=chngptm(formula.1=logratio~1, formula.2=~range, lidar, type=type, family="gaussian", 
        var.type="bootstrap", ci.bootstrap.size=100)
    print(summary(fit))
    for (i in 1:3) plot(fit, which=i)
    out=predict(fit)
    plot(lidar$range, out, main=type)
}


# with weights
dat.1=sim.chngpt("thresholded", "segmented", n=200, seed=1, beta=1, alpha=-1, x.distr="norm", e.=4,
    family="gaussian")
fit.1.a=chngptm(formula.1=y~z, formula.2=~x, family="gaussian", dat.1, type="segmented", 
    est.method="fastgrid", var.type="bootstrap", weights=ifelse(dat.1$x<3.5,100,1)
    , ci.bootstrap.size=10)
summary(fit.1.a)
plot(fit.1.a)
# fit.1.a$vcov$boot.samples

\dontrun{
# likelihood test, combination of slopes
dat=sim.chngpt("thresholded", "segmented", n=200, seed=1, beta=1, alpha=-1, x.distr="norm", e.=4,
    family="gaussian")
fit=chngptm(y~z, ~x, family="gaussian", dat, type="segmented", ci.bootstrap.size=100)
fit.0=lm(y~1,dat)
# likelihood ratio test using lmtest::lrtest
library(lmtest)
lrtest(fit, fit.0)
# estimate the slope after threshold using lincomb function in the chngpt package
lincomb(fit, c(0,0,1,1))
}


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

fit.4a=chngptm(formula.1=counts ~ 1, formula.2=~x, data=d.AD, family="quasipoisson", 
    type="segmented", var.type="bootstrap", verbose=1, ci.bootstrap.size=1)



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
fit=chngptm(formula.1=mpg ~ hp, formula.2=~hp*drat, mtcars, type="M10", 
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
    

# examples of aux.fit
fit.0=glm(yy~zz+ns(xx,df=3), data, family="binomial")
fit = chngptm (formula.1=yy~zz, formula.2=~xx, family="binomial", data, type="hinge", 
    est.method="smoothapprox", var.type="all", verbose=verbose, aux.fit=fit.0, 
    lb.quantile=0.1, ub.quantile=0.9, tol=1e-4, maxit=1e3)



}

# example of random intercept
dat=sim.twophase.ran.inte(threshold.type="segmented", n=50, seed=1)
fit = chngptm (formula.1=y~z+(1|id), formula.2=~x, family="gaussian", dat, 
    type="segmented", est.method="grid", var.type="bootstrap", ci.bootstrap.size=1)
plot(fit)
out=predict(fit, re.form=NA)
plot(dat$x, out)
out.1=predict(fit, type="response", re.form=NULL)# includes re
plot(dat$x, out.1, type="p", xlab="x")






}
