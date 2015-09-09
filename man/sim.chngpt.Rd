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

sim.chngpt  (label, 
    type=c("step","hinge","segmented","stegmented"),
    n, seed, 
    beta, coef.z=log(1.4), 
    x.distr="norm", 
    e., b.=-Inf,
    alpha.candidate=NULL, verbose=FALSE) 

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{label}{
string. Simulate scenario, see details.
}
  \item{type}{
string. Types of threshold effect.
}
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{seed}{
%%     ~~Describe \code{seed} here~~
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
  \item{verbose}{Boolean}
  \item{b.}{
%%     ~~Describe \code{b.} here~~
}
  \item{alpha.candidate}{candidate values of alpha, used in code to determine alpha values}
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
