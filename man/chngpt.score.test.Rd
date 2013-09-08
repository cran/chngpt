\name{chngpt.score.test}
\alias{antoka.test}
\alias{chngpt.score.test}
\alias{plot.chngpt.score.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Change Point Tests
}
\description{
Performs the change point tests.
}
\usage{
chngpt.score.test (formula.null, formula.chngpt, data, interaction.only=TRUE, 
    chngpts=NULL, lb.quantile=.1, ub.quantile=.9, chngpts.cnt=50, b.=-30, mc.n=5e4,
    ret.p.val=TRUE, verbose=TRUE)

antoka.test (formula, data, chngpt.var, plot.=FALSE)

\S3method{plot}{chngpt.score.test}(x, \dots)

}
\arguments{
  \item{formula.null}{
formula for the null model.
}
  \item{formula.chngpt}{
formula for the change point model.
}
\item{interaction.only}{
    boolean. 
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
  \item{b.}{
number. Slope of the sigmoid function in the derivation of the score statistics.
}
  \item{mc.n}{
integer. Number of multivariate normal samples to generate to evaluate p-value.
}
  \item{ret.p.val}{
Boolean. 
}
  \item{verbose}{
Boolean. 
}
  \item{formula}{
formula.
}
  \item{plot.}{
Boolean. 
}
  \item{x}{
An object of type chngpt.score.test. 
}
  \item{\dots}{arguments passed to or from methods}

}
\details{
    If there are missing values in the chngpt formula, those rows will be removed from the whole dataset, including null model and chngpt model
}
\value{
A list of class chngpt.score.test
\item{p.value}{P-value}
\item{chngpts}{Vector of change points evaluated}
\item{TT}{Standardized absolute score statistics}
\item{V.S.hat}{Estimated variance-covariance matrix of the score statistics}
}

\examples{

seed=2
beta=log(.6)
alpha=sim.alphas[["sigmoid2_norm"]]["4.5", paste(log(.6),"",sep="")]
if(is.null(alpha)) stop("alpha not found")
data=sim.sigmoid("sigmoid2", n=250, seed=seed, alpha, beta, x.distr="norm", e.=4.5, b.=-Inf)


test = chngpt.score.test(formula.null=y~z, formula.chngpt=~x, data)
test
plot(test)


test = chngpt.score.test(formula.null=y~z, formula.chngpt=~x*z, data)
test
plot(test)


#antoka.test(y~1, dat, "x.star")


}
