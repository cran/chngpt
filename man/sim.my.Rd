\name{sim.my}
\alias{sim.my}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate data
}
\description{
Simulate data
}
\usage{
sim.my(n, seed, label, alpha, beta, e. = NULL, b. = NULL, tr. = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
Sample size}
  \item{seed}{
Seed for random number generator
}
  \item{label}{
A character string which specifies the simulation scenario. sigmoid4, sigmoidgam4, elbow4
}
  \item{alpha}{
regression parameter
}
  \item{beta}{
regression parameter
}
  \item{e.}{
inflection point for the logistic transformation (the log scale)
}
  \item{b.}{
slope for the logistic transformation
}
  \item{tr.}{
threshold point
}
}
\details{
When the label starts with elbow, the transformation on x.star is elbow shaped.
When the label starts with sigmoid, the transformation on x.star is sigmoid shaped. 
Data simulated from logit(Pr(Y==1))=alpha + beta*(transformed x.star).

}
\value{
A data frame with columns: y, x.star, x.star.expit (if label starts with sigmoid), x.star.tr (if label starts with elbow), x.bin.med (x.star dichotomized at median), x.tri (x.star trichotomized at tertiles). 
}
\examples{

alpha=-1; beta=log(0.2)
e.=5; b.=-30; t.=1
dat=sim.my(n=250, seed=1, label="sigmoid4", alpha, beta, e.=e., b.=b.)

}
