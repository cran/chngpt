\name{convert.coef predictx threshold.func}
\alias{convert.coef}
\alias{predictx}
\alias{threshold.func}
\title{
Helper functions 
}
\description{
Some helper functions. predictx returns confidence bands for predictions as functions of the change point variable. threshold.func returns thresholded covariates.
}
\usage{

convert.coef(coef.0, threshold.type)

predictx(fit, boot.ci.type = c("perc", "basic", "symm"), alpha
 = 0.05, xx = NULL, verbose = FALSE, return.boot =
 FALSE, include.intercept = FALSE, get.simultaneous =
 TRUE)
         
threshold.func(threshold.type, coef, xx, x.name, include.intercept=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{include.intercept}{
%%     ~~Describe \code{coef.0} here~~
}
  \item{coef.0}{
%%     ~~Describe \code{coef.0} here~~
}
  \item{threshold.type}{
%%     ~~Describe \code{threshold.type} here~~
}
  \item{get.simultaneous}{
%%     ~~Describe \code{threshold.type} here~~
}
  \item{return.boot}{
%%     ~~Describe \code{threshold.type} here~~
}

  \item{fit}{
%%     ~~Describe \code{fit} here~~
}
  \item{boot.ci.type}{
%%     ~~Describe \code{boot.ci.type} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}

  \item{coef}{
%%     ~~Describe \code{coef} here~~
}
  \item{xx}{
%%     ~~Describe \code{xx} here~~
}
  \item{x.name}{
%%     ~~Describe \code{x.name} here~~
}
}
