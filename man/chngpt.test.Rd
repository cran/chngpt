\name{chngpt.test}
\alias{antoch.test}
\alias{chngpt.test}
\alias{plot.chngpt.test}
\title{
Change Point Tests
}
\description{
Hypothesis testing for change point covariate in logistic regression and linear regression.
}
\usage{

chngpt.test (formula.null, formula.chngpt, family=c("binomial","gaussian"), data, 
    type=c("step","hinge","segmented","stegmented"),
    test.statistic=c("lr","score"), # support for score is gradually descreasing
    chngpts=NULL, lb.quantile=.1, ub.quantile=.9, 
    chngpts.cnt=50, #this is set to 25 if int is weighted.two.sided or weighted.one.sided
    prec.weights=NULL,
    p.val.method=c("MC","param.boot"), 
    mc.n=5e4, # 1e3 won't cut it, the p values estimated could be smaller than nominal
    boot.B=1e4,
    robust=FALSE,
    keep.fits=FALSE, verbose=FALSE
) 


antoch.test (formula, data, chngpt.var, plot.=FALSE)

\method{plot}{chngpt.test}(x, by.percentile=TRUE, both=FALSE, main=NULL, ...) 

}
\arguments{
  \item{formula.null}{
formula for the null model. The outcome needs to be a 0/1 variable.
}
  \item{formula.chngpt}{
formula for the change point model. 
For example, suppose formula.null=y~z and we want to test whether I(x>cutff) is a significant predictor, formula.chngpt=~x
If instead we are interested in testing the null that neither I(x>cutff) nor z*I(x>cutff) is a significant predictor, formula.chngpt=~x*z
}
  \item{data}{data frame.}
  \item{family}{Currently only linear and logistic regression are supported.}
  \item{type}{step: flat before and after change point; hinge: flat before and slope after change point; segmented: slope before and after change point}
  \item{test.statistic}{method for testing main effects of some threshold model.}
  \item{chngpts}{
A grid of potential change points to maximize over. If not supplied, they will be 
set to a vector of length chngpt.cnt equally spaced between lb.quantile and ub.quantile.
}
  \item{robust}{
Boolean.
}
  \item{lb.quantile}{
number. The lower bound in the search for change point in the unit of quantile.
}
  \item{ub.quantile}{
number. The upper bound in the search for change point in the unit of quantile.
}
  \item{chngpts.cnt}{
integer. Number of potential change points to maximize over.
}
  \item{mc.n}{
integer. Number of multivariate normal samples to generate in the Monte Carlo procedure to evaluate p-value.
}

  \item{verbose}{
Boolean. 
}
  \item{chngpt.var}{
string. Name of the predictor to detect change point
}

  \item{plot.}{
Boolean. Whether to make a plot.
}
  \item{formula}{
formula.
}
  \item{x}{
An object of type chngpt.test. 
}
  \item{\dots}{arguments passed to or from methods}
  \item{by.percentile}{
%%     ~~Describe \code{expr} here~~
}
  \item{both}{
%%     ~~Describe \code{expr} here~~
}
  \item{main}{
%%     ~~Describe \code{expr} here~~
}
  \item{prec.weights}{
%%     ~~Describe \code{expr} here~~
}
  \item{p.val.method}{
%%     ~~Describe \code{expr} here~~
}
  \item{boot.B}{
%%     ~~Describe \code{expr} here~~
}
  \item{keep.fits}{
%%     ~~Describe \code{expr} here~~
}

}
\details{

The model under the altnerative is the model under the null plus terms involving the threshold. For example, when the type is segmented and formula.null=~z, formula.chngpt=~x, the model under the null is ~z+x and the model under the alternative is ~z+x+(x-e)_+.

If there are missing values in the chngpt formula, those rows will be removed from the whole dataset, including null model and chngpt model.

antoch.test is only implemented for main effect only and is based on Antoch et al. (2004). Also see Fong et al. (2014).

}
\value{
A list of class htest and chngpt.test
\item{p.value}{P-value}
\item{family}{Family from input}
\item{method}{Method from input}
}
\references{

Fong Y, Di C, and Permar S. (2015) Change-Point Testing in Logistic Regression Models with Interaction Term. Statistics in Medicine. 34:1483--1494

Pastor-Barriuso, R. and Guallar, E. and Coresh, J. (2003) Transition models for change-point estimation in logistic regression. Statistics in Medicine. 22:13141

Antoch, J. and Gregoire, G. and Jaruskova, D. (2004) Detection of structural changes in generalized linear models. Statistics and probability letters. 69:315

}
\examples{


dat=sim.chngpt("thresholded", "step", n=200, seed=1, beta=1, alpha=-1, x.distr="norm", e.=4, 
    family="binomial")
test=chngpt.test(formula.null=y~z, formula.chngpt=~x, dat, type="step", family="binomial",
    mc.n=10)
test
plot(test)

dat=sim.chngpt("thresholded", "segmented", n=200, seed=1, beta=1, alpha=-1, x.distr="norm", e.=4,
    family="binomial")
test=chngpt.test(formula.null=y~z, formula.chngpt=~x, dat, type="segmented", family="binomial",
    mc.n=10)
test
plot(test)

\dontrun{
# not run because otherwise the examples take >5s and that is a problem for R CMD check

# has interaction
test = chngpt.test(formula.null=y~z, formula.chngpt=~x*z, dat, type="step", family="binomial")
test
plot(test)
}


}
