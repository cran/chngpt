\name{chngpt.test}
\alias{antoch.test}
\alias{chngpt.test}
\alias{plot.chngpt.test}
\title{
Change Point Tests
}
\description{
Hypothesis testing for change point covariate in logistic regression.
}
\usage{
chngpt.test (formula.null, formula.chngpt, data, 
    type=c("step","hinge","segmented","stegmented"),
    main.method=c("lr","score"),
    chngpts=NULL, lb.quantile=.1, ub.quantile=.9, 
    chngpts.cnt=50, # this is set to 25 if int is weighted.two.sided or weighted.one.sided
    single.weight=1,
    mc.n=5e4, 
    prob.weights=NULL,
    verbose=FALSE 
) 

antoch.test (formula, data, chngpt.var, plot.=FALSE)

\S3method{plot}{chngpt.test}(x, \dots)

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
  \item{data}{
data frame.
}
  \item{type}{step: flat before and after change point; hinge: flat before and slope after change point; segmented: slope before and after change point}
  \item{main.method}{method for testing main effects of some threshold model.}
  \item{chngpts}{
A grid of potential change points to maximize over. If not supplied, they will be 
set to a vector of length chngpt.cnt equally spaced between lb.quantile and ub.quantile.
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
  \item{single.weight}{
numeric. Only used when interaction.method="weighted.single.arg".
}
  \item{prob.weights}{
numeric. Not yet fully supported.
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

}
\details{

If there are missing values in the chngpt formula, those rows will be removed from the whole dataset, including null model and chngpt model.

The test results, e.g. p-value, should not change if the same function is called twice, because internally seed to random number generator is set to 1 before MC and set back to the R state afterwards

antoch.test is only implemented for main effect only and is based on Antoch et al. (2004). Also see Fong et al. (2014).
}
\value{
A list of class htest and chngpt.test
\item{p.value}{P-value}
\item{chngpts}{Vector of change points evaluated}
\item{TT}{Standardized absolute score statistics}
\item{V.S.hat}{Estimated variance-covariance matrix of the score statistics}
}
\references{

Fong Y, Di C, and Permar S. (2014) Change-Point Testing in Logistic Regression Models with Interaction Term. Statistics in Medicine.

Pastor-Barriuso, R. and Guallar, E. and Coresh, J. (2003) Transition models for change-point estimation in logistic regression. Statistics in Medicine. 22:13141

Antoch, J. and Gregoire, G. and Jaruskova, D. (2004) Detection of structural changes in generalized linear models. Statistics and probability letters. 69:315

}
\examples{

seed=2
beta=log(.6)
data=sim.chngpt("sigmoid2", type="step", n=200, seed=seed, beta, x.distr="norm", e.=4.5, b.=-Inf)

test = chngpt.test(formula.null=y~z, formula.chngpt=~x, data, type="step")
test
plot(test)


\dontrun{
# not run because otherwise the examples take >5s and that is a problem for R CMD check
# has interaction
test = chngpt.test(formula.null=y~z, formula.chngpt=~x*z, data)
test
plot(test)
}

#antoch.test(y~1, dat, "x.star")


}
