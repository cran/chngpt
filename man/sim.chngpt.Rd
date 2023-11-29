\name{sim.chngpt}
\alias{sim.twophase.ran.inte}
\alias{sim.threephase}
\alias{sim.chngpt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulation Function
}
\description{
Generate simulation datasets for change point Monte Carlo studies.
}
\usage{

sim.chngpt (mean.model = c("thresholded", "thresholdedItxn",
 "quadratic", "quadratic2b", "cubic2b", "exp",
 "flatHyperbolic", "z2", "z2hinge", "z2segmented",
 "z2linear", "logistic"), threshold.type = c("NA",
 "M01", "M02", "M03", "M10", "M20", "M30", "M11",
 "M21", "M12", "M22", "M22c", "M31", "M13", "M33c",
 "hinge", "segmented", "upperhinge", "segmented2",
 "step", "stegmented"), b.transition = Inf, family =
 c("binomial", "gaussian"), x.distr = c("norm",
 "norm3", "norm6", "imb", "lin", "mix", "gam",
 "zbinary", "gam1", "gam2", "fixnorm", "unif"), e. =
 NULL, mu.x = 4.7, sd.x = NULL, sd = 0.3, mu.z = 0,
 alpha = NULL, alpha.candidate = NULL, coef.z =
 log(1.4), beta = NULL, beta.itxn = NULL,
 logistic.slope = 15, n, seed, weighted = FALSE,
 heteroscedastic = FALSE, ar = FALSE, verbose = FALSE)

sim.twophase.ran.inte(threshold.type, n, seed)

sim.threephase(n, seed, gamma = 1, e = 3, beta_e = 5, f = 7, beta_f = 2, coef.z = 1)



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
n
}
  \item{mu.z}{
n
}
  \item{seed}{
seed
}
  \item{weighted}{
beta
}
  \item{beta}{
beta
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
e.
}
  \item{verbose}{
  Boolean
}
  \item{b.transition}{
b.
}
  \item{sd}{
b.
}
  \item{ar}{
autocorrelation
}
  \item{alpha.candidate}{
Candidate values of alpha, used in code to determine alpha values
}
  \item{e}{
e
}
  \item{beta_e}{
beta_e
}
  \item{f}{
f
}
  \item{beta_f}{
beta_f
}
  \item{logistic.slope}{
beta_f
}
  \item{gamma}{
beta_f
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
