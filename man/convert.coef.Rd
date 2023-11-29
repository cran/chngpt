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
coef.0
}
  \item{coef.0}{
coef.0
}
  \item{threshold.type}{
threshold.type
}
  \item{get.simultaneous}{
threshold.type
}
  \item{return.boot}{
threshold.type
}

  \item{fit}{
fit
}
  \item{boot.ci.type}{
boot.ci.type
}
  \item{alpha}{
alpha
}
  \item{verbose}{
verbose
}

  \item{coef}{
coef
}
  \item{xx}{
xx
}
  \item{x.name}{
x.name
}
}
