\name{logLikuh}
\alias{logLikuh}

\title{Computes the Log-likelihood Value of a U-shaped Hazard Function}

\description{
  
  \code{logLikuh} returns the log-likelihood value of a U-shaped hazard
  function, given a data set.
  
}

\usage{

logLikuh(h, data)

}

\arguments{
  
  \item{h}{an object of class \code{uh}. }
  
  \item{data}{numeric vector or matrix for exact or interval-censored
    observations, or an object of class \code{icendata}.}
  
}

\value{

  Log-likelihood value evaluated at \code{h}, given \code{data}.
  
}

\author{ Yong Wang <yongwang@auckland.ac.nz>}

\references{

  Wang, Y. and Fani, S. (2015). Nonparametric maximum likelihood
  computation of a U-shaped hazard function. (In preparation)

}

\seealso{ \code{\link{Uhaz}}, \code{\link{icendata}}, \code{\link{plot.uh}} }

\examples{
data(ap)
(h0 = uh(.2, NULL, NULL, NULL, NULL, 15, 1))   # Uniform hazard
plot(h0, ylim=c(0,.3))
logLikuh(h0, ap)

r = Uhaz(ap, deg=2)
r$ll
logLikuh(r$h, ap)
plot(r$h, add=TRUE, col="red3")
}

\keyword{ function }   % at least one, from doc/KEYWORDS
