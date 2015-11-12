\name{acfail}
\alias{acfail}

\docType{data}

\title{
  Air Conditioner Failure Data
}

\description{ Contains the number of operating hours between successive
  failure times of the air conditioning systems in Boeing airplanes }

\usage{acfail}

\format{ A numeric vector storing the failure times.}
  
\source{ Proschan (1963) }

\references{ Proschan, F. (1963). Theoretical explanation of observed
  decreasing failure rate. \emph{Technometrics}, \bold{5}, 375-383. }

\seealso{ \code{\link{Uhaz}}. }

\examples{
data(acfail)
r = Uhaz(acfail, deg=2)
plot(r$h, fn="h")
plot(r$h, fn="d")
}

\keyword{datasets}

