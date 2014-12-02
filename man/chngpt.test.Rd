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
    interaction.method=c(
        "lr.mc", "lr.pastor", 
        "weighted.two.sided", "weighted.one.sided", "weighted.single.arg", 
        "main.itxn", "main.only", "itxn.only"), 
    chngpts=NULL, lb.quantile=.1, ub.quantile=.9, chngpts.cnt=50, b.=-30, 
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
  \item{interaction.method}{
Only applicable when formula.chngpt contains interaction. See details.
}
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
  \item{b.}{
numeric. Slope of logistic function that approximates the step function 
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

When formula.chngpt does not contain interaction term, a test based on the maximum of score statistics is performed. 
When formula.chngpt contains interaction term, one of the following tests is performed depending on the value of interaction.method. Unless otherwise specified, method for determining p-values is based on Monte Carlo (Fong et al. 2014).
\preformatted{
    lr.mc: test statistic is maximum of likelihood ratio statistics. 
    lr.pastor: test statistic is maximum of likelihood ratio statistics. Method for determining p-values is based on inequality (Pastor-Barriuso et al. 2003)
    weighted.two.sided: test statistic is maximum of weighted combination of scores that covers both directions of interaction.
    weighted.one.sided: test statistic is maximum of weighted combination of scores that assumes power is enhanced when beta1 and beta2 have the same sign. 
    weighted.single.arg: test statistic is a weighted combination of scores using weight provided in single.weight. 
    main.itxn: test statistic is maximum of two vectors of scores corresponding to the score for the main effect and the score for the interaction effect.
    main.only: test statistic is maximum of one vectors of scores corresponding to the score for the main effect.
    itxn.only: test statistic is maximum of one vectors of scores corresponding to the score for the interaction effect.
}

If there are missing values in the chngpt formula, those rows will be removed from the whole dataset, including null model and chngpt model

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
