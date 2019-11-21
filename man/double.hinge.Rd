\name{double.hinge}
\alias{double.hinge}
\alias{plot.double.hinge}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit Double Hinge Models
}
\description{
Fit double hinge models.
}
\usage{
double.hinge(x, y, mid.x = NULL, lb.quantile = 0.05, ub.quantile = 0.95, 
    lower.y = 0, upper.y = 100, var.type = c("none", "bootstrap"), 
    ci.bootstrap.size = 1000, alpha = 0.05, save.boot = TRUE, ncpus = 1)

\method{plot}{double.hinge} (x, lcol=1, lwd=2, ...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{mid.x}{
%%     ~~Describe \code{mid.x} here~~
}
  \item{lb.quantile}{
%%     ~~Describe \code{lb.quantile} here~~
}
  \item{ub.quantile}{
%%     ~~Describe \code{ub.quantile} here~~
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
  \item{\dots}{arguments passed along}

}
