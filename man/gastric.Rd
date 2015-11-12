\name{gastric}
\alias{gastric}

\docType{data}

\title{
  Gastric Cancer Survival Data
}

\description{ Contains the survival times of 45 gastrointestinal tumor
  patients who were treated with both chemotherapy and radiotherapy. It
  has both exact and right-censored observations. }

\usage{gastric}

\format{ A data frame with 30 observations and 3 variables:
  
  L: left-end points of the interval-censored survival times;

  R: right-end points of the interval-censored survival times.}

\source{ Klein and Moeschberger (2003), page 224. }

\references{ Klein, J. P. and Moeschberger, M. L. (2003).  \emph{Survival
  Analysis: Techniques for Censored and Truncated Data (2nd ed.)}.
  Springer. }

\seealso{ \code{\link{npsurv}}. }

\examples{
data(gastric)
r = Uhaz(gastric, deg=2)
plot(r$h, fn="h")   # hazard function
plot(r$h, fn="d")   # density function
plot(r$h, fn="s")   # survival function
}

\keyword{datasets}

