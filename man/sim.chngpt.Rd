\name{sim.chngpt}
\alias{sim.chngpt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulation Function
}
\description{
Generate simulation datasets for change point Monte Carlo studies.
}
\usage{

sim.chngpt (
    mean.model=c("thresholded","thresholdedItxn","quadratic","quadratic2b","cubic2b",
        "exp","flatHyperbolic","z2","z2hinge","z2segmented","z2linear"), 
    threshold.type=c("NA","step","hinge","segmented","segmented2","stegmented"),
#segmented2 differs from segmented in parameterization, it is studied in Cheng 2008
    b.transition=Inf,
    family=c("binomial","gaussian"), 
    x.distr=c("norm","norm3","norm6","imb","lin","mix","gam","zbinary","gam1","gam2",
        "fixnorm", "fixnorm3", "fixnorm6"), # gam1 is a hack to allow e. be different
    e.=NULL, mu.x=4.7, sd.x=NULL, sd=0.3, 
    alpha=NULL, alpha.candidate=NULL, coef.z=log(1.4), beta=NULL, beta.itxn=NULL, 
    n, seed, 
    weighted=FALSE, # sampling weights
    heteroscedastic=FALSE,
    verbose=FALSE) 

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{threshold.type}{
string. Types of threshold effect to simulate, only applicable when label does not start with sigmoid.
}
  \item{family}{
string. Glm family.
}
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{seed}{
%%     ~~Describe \code{seed} here~~
}
  \item{weighted}{
%%     ~~Describe \code{beta} here~~
}
  \item{beta}{
%%     ~~Describe \code{beta} here~~
}
  \item{coef.z}{
numeric. Coefficient for z.
}
  \item{beta.itxn}{
numeric. Coefficient for z.
}
  \item{alpha}{
numeric, intercept.
}
  \item{mu.x}{
numeric
}
  \item{sd.x}{
numeric
}
  \item{mean.model}{
numeric
}
  \item{x.distr}{
string. Possible values: norm (normal distribution), gam (gamma distribution). gam1 is a hack to allow e. be different
}
  \item{e.}{
%%     ~~Describe \code{e.} here~~
}
  \item{verbose}{
  Boolean
}
  \item{b.transition}{
%%     ~~Describe \code{b.} here~~
}
  \item{sd}{
%%     ~~Describe \code{b.} here~~
}
  \item{alpha.candidate}{
Candidate values of alpha, used in code to determine alpha values
}

  \item{heteroscedastic}{Boolean.}

}

\details{
mean.model, threshold.type and b.transition all affect mean models. 
}
\value{
A data frame with following columns: 
    \item{y}{0/1 outcome}
    \item{x}{observed covariate that we are interested in}
    \item{x.star}{unobserved covariate that underlies x}
    \item{z}{additional covariate}
In addition, columns starting with 'w' are covariates that we also adjust in the model; columns starting with 'x' are covariates derived from x.
}
\examples{

seed=2
par(mfrow=c(2,2))
dat=sim.chngpt(mean.model="thresholded", threshold.type="hinge", family="gaussian", beta=0, n=200, 
    seed=seed, alpha=-1, x.distr="norm", e.=4, heteroscedastic=FALSE)
plot(y~z, dat)
dat=sim.chngpt(mean.model="thresholded", threshold.type="hinge", family="gaussian", beta=0, n=200, 
    seed=seed, alpha=-1, x.distr="norm", e.=4, heteroscedastic=TRUE)
plot(y~z, dat)
dat=sim.chngpt(mean.model="z2", threshold.type="hinge", family="gaussian", beta=1, n=200, 
    seed=seed, alpha=1, x.distr="norm", e.=4, heteroscedastic=FALSE)
plot(y~z, dat)
dat=sim.chngpt(mean.model="z2", threshold.type="hinge", family="gaussian", beta=1, n=200, 
    seed=seed, alpha=1, x.distr="norm", e.=4, heteroscedastic=TRUE)
plot(y~z, dat)

}
