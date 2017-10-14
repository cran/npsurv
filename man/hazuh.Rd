\name{hazuh}
\alias{hazuh}
\alias{chazuh}
\alias{survuh}
\alias{denuh}

\title{Distributional Functions given a U-shaped Hazard Function}

\description{
  
  Given an object of class \code{uh}: 
  
  \code{hazuh} computes the hazard values;
  
  \code{chazuh} computes the cumulative hazard values;
  
  \code{survuh} computes the survival function values;
  
  \code{denuh} computes the density function values.
  
}

\usage{

hazuh(t, h)
chazuh(t, h)
survuh(t, h)
denuh(t, h)

}

\arguments{
  
  \item{t}{time points at which the function is to be evaluated. }
  
  \item{h}{an object of class \code{uh}.}
  
}

\value{

  A numeric vector of the function values. 
  
}

\author{ Yong Wang <yongwang@auckland.ac.nz>}

\references{

  Wang, Y. and Fani, S. (2017). Nonparametric maximum likelihood
  computation of a U-shaped hazard function. \emph{Statistics and
  Computing}, (in print).

}

\seealso{ \code{\link{Uhaz}}, \code{\link{icendata}}, \code{\link{plot.uh}} }

\examples{
data(ap)
h = Uhaz(icendata(ap), deg=2)$h
hazuh(0:15, h)     # hazard
chazuh(0:15, h)    # cumulative hazard
survuh(0:15, h)    # survival probability
denuh(0:15, h)     # density
}

\keyword{ function }   % at least one, from doc/KEYWORDS
