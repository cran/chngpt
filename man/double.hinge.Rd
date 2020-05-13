\name{double.hinge}
\alias{double.hinge}
\alias{plot.double.hinge}
\alias{fitted.double.hinge}
\alias{residuals.double.hinge}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit Double Hinge Models
}
\description{
Fit double hinge models.
}
\usage{
double.hinge(x, y, lower.y = NULL, upper.y = NULL,
         var.type = c("none", "bootstrap"), ci.bootstrap.size =
         1000, alpha = 0.05, save.boot = TRUE, ncpus = 1)

\method{plot}{double.hinge} (x, which = NULL, xlim = NULL, 
 lwd = 2, lcol = "red",
 lty = 1, add.points = TRUE, add.ci = TRUE, breaks =
 20, mark.chngpt = FALSE, xlab = NULL, ylab = NULL,
 ...) 
\method{fitted}{double.hinge} (object, ...) 
\method{residuals}{double.hinge} (object, ...) 
}
\arguments{
  \item{object}{
%%     ~~Describe \code{x} here~~
}
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{lower.y}{
%%     ~~Describe \code{lower.y} here~~
}
  \item{upper.y}{
%%     ~~Describe \code{upper.y} here~~
}
  \item{var.type}{
%%     ~~Describe \code{var.type} here~~
}
  \item{ci.bootstrap.size}{
%%     ~~Describe \code{ci.bootstrap.size} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{save.boot}{
%%     ~~Describe \code{save.boot} here~~
}
  \item{ncpus}{
%%     ~~Describe \code{ncpus} here~~
}
  \item{lcol}{
%%     ~~Describe \code{ncpus} here~~
}
  \item{lwd}{
%%     ~~Describe \code{ncpus} here~~
}

 'which' 'xlim' 'lty' 'add.points' 'add.ci' 'breaks' 'mark.chngpt'
 'xlab' 'ylab'


  \item{which}{
%%     ~~Describe \code{x} here~~
}
  \item{xlim}{
%%     ~~Describe \code{x} here~~
}
  \item{lty}{
%%     ~~Describe \code{x} here~~
}
  \item{add.points}{
%%     ~~Describe \code{x} here~~
}
  \item{add.ci}{
%%     ~~Describe \code{x} here~~
}
  \item{breaks}{
%%     ~~Describe \code{x} here~~
}
  \item{mark.chngpt}{
%%     ~~Describe \code{x} here~~
}
  \item{xlab}{
%%     ~~Describe \code{x} here~~
}
  \item{ylab}{
%%     ~~Describe \code{x} here~~
}


  \item{\dots}{arguments passed along}

}

\details{
If lower.y and upper.y are not supplied, min(y) is taken as the function value when x is less than or equal to the first threshold, and max(y) is taken as the function value when x is greater than or equal to the second threshold.

If the function is expected to be decreasing between the two thresholds, lower.y and upper.y should be supplied to ensure the correct fit.

mse is residual sum of squares

}
