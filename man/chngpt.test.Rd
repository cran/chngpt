\name{chngpt.test}
\alias{antoch.test}
\alias{chngpt.test}
\alias{plot.chngpt.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Change Point Tests
}
\description{
Performs the change point tests.
}
\usage{
chngpt.test (formula.null, formula.chngpt, data, 
    interaction.method=c(
        "lr.mc", "lr.pastor", 
        "weighted.two.sided", "weighted.one.sided", "weighted.single.arg", 
        "main.itxn", "main.only", "itxn.only"), 
    chngpts=NULL, lb.quantile=.1, ub.quantile=.9, 
    chngpts.cnt=50, 
    b.=-30, single.weight=1,
    p.val.method=c("max"), ret.p.val=TRUE, mc.n=5e4, 
    prob.weights=NULL,
    verbose=FALSE 
)

antoch.test (formula, data, chngpt.var, plot.=FALSE)

\S3method{plot}{chngpt.test}(x, \dots)

}
\arguments{
  \item{formula.null}{
formula for the null model.
}
  \item{formula.chngpt}{
formula for the change point model. If interaction, one of the covariates must be in the null model
}
  \item{data}{
data frame.
}
  \item{chngpt.var}{
string. Name of the predictor to detect change point
}
  \item{chngpts}{
number or vector of numbers. Potential change points to evaluate. Optional. If no supplied, 
set to be a vector of length chngpt.cnt equally spaced between lb.quantile and ub.quantile.
}
  \item{lb.quantile}{
number. The lower bound in the search for change point in the unit of quantile.
}
  \item{ub.quantile}{
number. The upper bound in the search for change point in the unit of quantile.
}
  \item{chngpts.cnt}{
integer. Number of potential change points to evalute between the lower bound and the upper bound.
}
  \item{mc.n}{
integer. Number of multivariate normal samples to generate to evaluate p-value.
}
  \item{ret.p.val}{
Boolean. 
}
  \item{p.val.method}{
string. 
}
  \item{b.}{
numeric. Slope of logistic function that approximates the step function 
}
  \item{interaction.method}{
string. See details.
}
  \item{single.weight}{
numeric. 
}
  \item{prob.weights}{
numeric. 
}
  \item{verbose}{
Boolean. 
}
  \item{plot.}{
Boolean. 
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
    If there are missing values in the chngpt formula, those rows will be removed from the whole dataset, including null model and chngpt model

Different values of interaction.method correspond to maximum of the following.
\preformatted{
    lr.mc: likelihood ratio statistics
    weighted.two.sided: a series of weighted transformation of scores that cover both directions of interaction
    weighted.one.sided: a series of weighted transformation of scores that assume the direction of interaction coefficient is known
    main.only: 
}
}
\value{
A list of class chngpt.test
\item{p.value}{P-value}
\item{chngpts}{Vector of change points evaluated}
\item{TT}{Standardized absolute score statistics}
\item{V.S.hat}{Estimated variance-covariance matrix of the score statistics}
}
\references{
Fong Y, Di C, and Permar S. (2013) Change-Point Testing in Logistic Regression Models with Interaction Term. Submitted.
}
\examples{

seed=2
beta=log(.6)
alpha=sim.alphas[["sigmoid2_norm"]]["4.5", paste(log(.6),"",sep="")]
if(is.null(alpha)) stop("alpha not found")
data=sim.sigmoid("sigmoid2", n=250, seed=seed, alpha, beta, x.distr="norm", e.=4.5, b.=-Inf)


test = chngpt.test(formula.null=y~z, formula.chngpt=~x, data)
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
