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
         1000, alpha = 0.05, save.boot = TRUE, ncpus = 1, 
         boot.ci.type=c("percentile","symmetric"))

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
x
}
  \item{x}{
x
}
  \item{y}{
y
}
  \item{lower.y}{
lower.y
}
  \item{upper.y}{
upper.y
}
  \item{var.type}{
var.type
}
  \item{boot.ci.type}{
var.type
}
  \item{ci.bootstrap.size}{
ci.bootstrap.size
}
  \item{alpha}{
alpha
}
  \item{save.boot}{
save.boot
}
  \item{ncpus}{
ncpus
}
  \item{lcol}{
ncpus
}
  \item{lwd}{
ncpus
}

  \item{which}{
x
}
  \item{xlim}{
x
}
  \item{lty}{
x
}
  \item{add.points}{
x
}
  \item{add.ci}{
x
}
  \item{breaks}{
x
}
  \item{mark.chngpt}{
x
}
  \item{xlab}{
x
}
  \item{ylab}{
x
}


  \item{\dots}{arguments passed along}

}

\details{
If lower.y and upper.y are not supplied, min(y) is taken as the function value when x is less than or equal to the first threshold, and max(y) is taken as the function value when x is greater than or equal to the second threshold.

If the function is expected to be decreasing between the two thresholds, lower.y and upper.y should be supplied to ensure the correct fit.

mse is residual sum of squares

}
