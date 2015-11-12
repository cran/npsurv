\name{km}

\alias{km}

\title{Kaplan-Meier Estimation}
\description{
  
  \code{km} computes the nonparametric maximum likelihood esimate
  (NPMLE) of a survival function for right-censored data.
  
}

\usage{

km(data, w=1)

}

\arguments{
  
  \item{data}{vector or matrix, or an object of class \code{icendata}.}
  
  \item{w}{weights/multiplicities of observations.}

}

\details{

  For details about the arguments, see \code{icendata}.

}

\value{

  A list consisting of:
  
  \item{f}{NPMLE, an object of class \code{idf}.}
  
  \item{ll}{log-likelihood value of the NPMLE \code{f}.}
  
}

\author{ Yong Wang <yongwang@auckland.ac.nz>}

\references{

  Kaplan, E. L. and Meier, P. (1958). Nonparametric estimation from
  incomplete observations. \emph{Journal of the American Statistical
  Association}, \bold{53}, 457-481.

}

\seealso{ \code{\link{icendata}}, \code{\link{npsurv}}, \code{\link{idf}}. }

\examples{
x = cbind(1:5, c(1,Inf,3,4,Inf))
(f = km(x)$f)
plot(f)

data(leukemia)
i = leukemia[,"group"] == "Placebo"
plot(km(leukemia[i,1:2])$f, xlim=c(0,40), col="green3") # placebo
plot(km(leukemia[!i,1:2])$f, add=TRUE)                  # 6-MP
}
\keyword{ function }   % at least one, from doc/KEYWORDS
