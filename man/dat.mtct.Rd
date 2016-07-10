\name{dat.mtct}
\alias{dat.mtct}
\docType{data}
\title{
An Example Dataset
}
\description{
A dataset from the immune correlates study of Maternal To Child Transmission of HIV-1
}
\usage{data("dat.mtct")}
\format{
  A data frame with 236 observations on the following 3 variables.
  \describe{
    \item{\code{y}}{a numeric vector}
    \item{\code{birth}}{a factor with levels \code{C-section} \code{Vaginal}}
    \item{\code{NAb_SF162LS}}{a numeric vector}
  }
}
\references{
Permar, S. R., Fong, Y., Nathan Vandergrift, Genevieve G. Fouda, Peter Gilbert, Georgia D. Tomaras, Feng Gao and Barton F. Haynes et al.
(2015) Maternal HIV-1 Envelope variable loop 3-specific IgG responses and reduced risk of perinatal
transmission. Journal of Clinical Investigation, 125(7):2702:2706.

Fong, Y., Huang, Y., Gilbert, P. (2015) Estimation and hypothesis testing under four types of change
point models using chngpt, under review.

}
\examples{
data(dat.mtct)
## maybe str(dat.mtct) ; plot(dat.mtct) ...
}
\keyword{datasets}
