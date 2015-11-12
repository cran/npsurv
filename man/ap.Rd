\name{ap}
\alias{ap}

\docType{data}

\title{
  Angina Pectoris Survival Data
}

\description{ Contains the survival times in years from the time of
  diagnosis for 2418 male patients with angina pectoris. Some patients
  are lost to follow-up, hence giving right-censored observations. Each
  integer-valued survival time is treated as being censored within a
  one-year interval. }

\usage{ap}

\format{

  A data frame with 30 observations and 3 variables:
  
  \code{L}: left-end point of an interval-censored retraction time;
  
  \code{R}: right-end point of an interval-censored retraction time;
  
  \code{count}: number of patients in the interval.

}

\source{ Lee and Wang (2003), page 92. }

\references{ Lee, E. T. and Wang, J. W. (2003). \emph{Statistical Methods for
Survival Data Analysis}. Wiley. }

\seealso{ \code{\link{npsurv}}. }

\examples{
data(ap)
r = Uhaz(ap, deg=2)           # smooth U-shaped hazard
plot(r$h, fn="h")             # hazard
plot(r$h, fn="d")             # density

# NPMLE and shape-restricted estimation
plot(npsurv(ap), fn="s")      # survival under no shape restriction
plot(r$h, fn="s", add=TRUE)   # survival with smooth U-shaped hazard
}

\keyword{datasets}

