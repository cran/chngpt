\name{sim.sigmoid}
\alias{sim.sigmoid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulation Function
}
\description{
Simulation Function
}
\usage{
sim.sigmoid (label, n, seed, alpha, beta, coef.z=log(1.4), x.distr="norm", e.=NULL, 
    b.=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{label}{
string. Simulate scenario, see details.
}
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{seed}{
%%     ~~Describe \code{seed} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{beta}{
%%     ~~Describe \code{beta} here~~
}
  \item{coef.z}{
numeric. Coefficient for z.
}
  \item{x.distr}{
string. Possible values: norm (normal distribution), gam (gamma distribution)
}
  \item{e.}{
%%     ~~Describe \code{e.} here~~
}
  \item{b.}{
%%     ~~Describe \code{b.} here~~
}
}
\details{
When label is "sigmoid1", an intercept only model is the data generative model. When label is "sigmoid2", a binary covariate z is also part of the data generative model. 
}
\value{
A data frame with following columns: 
    \item{y}{0/1 outcome}
    \item{x}{observed covariate that we are interested in}
    \item{x.star}{unobserved covariate that underlies x}
    \item{z}{additional covariate}
In addition, columns starting with 'w' are covariates that we also adjust in the model; columns starting with 'x' are covariates derived from x.
}
