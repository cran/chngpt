\name{sim.pastor}
\alias{sim.pastor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate data according to one of the scenarios considered in Pastor-Barriuso et al 2003
}
\description{
Simulate data according to one of the scenarios considered in Pastor-Barriuso et al 2003
}
\usage{
sim.pastor(seed)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seed}{
Seed for the random number generator.
}
}
\value{
A data frame with columns: y, x.star, x.star.expit, and x.bin.med (x.star dichotomized at median). 
}
\examples{

dat=sim.pastor(seed=1)

}
