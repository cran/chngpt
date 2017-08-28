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
    label=c("sigmoid2","sigmoid3","sigmoid4","sigmoid5","sigmoid6","quadratic",
            "quadratic2b","cubic2b","exp","flatHyperbolic"), 
    n, seed, 
    type=c("NA","step","hinge","segmented","segmented2","stegmented"),
    family=c("binomial","gaussian"),
    x.distr=c("norm","norm3","norm6","imb","lin","mix","gam","zbinary","gam1","gam2"),
    e.=NULL, b.transition=Inf,
    beta=NULL, coef.z=log(1.4), alpha=NULL,
    sd=0.3, mu=4.7, sd.x=NULL,
    alpha.candidate=NULL, verbose=FALSE) 

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{label}{
string. Simulate scenario, see details.
}
  \item{type}{
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
  \item{beta}{
%%     ~~Describe \code{beta} here~~
}
  \item{coef.z}{
numeric. Coefficient for z.
}
  \item{alpha}{
numeric, intercept.
}
  \item{mu}{
numeric
}
  \item{sd.x}{
numeric
}
  \item{x.distr}{
string. Possible values: norm (normal distribution), gam (gamma distribution). gam1 is a hack to allow e. be different
}
  \item{e.}{
%%     ~~Describe \code{e.} here~~
}
  \item{verbose}{Boolean}
  \item{b.transition}{
%%     ~~Describe \code{b.} here~~
}
  \item{sd}{
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
