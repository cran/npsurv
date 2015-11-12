\name{idf}
\alias{idf}
\alias{idf.object}
\alias{print.idf}

%- Also NEED an '\alias' for EACH other topic documented here.

\title{Interval Distribution Function}
\description{
  
  \code{idf} creates an object of class \code{idf}. An \code{idf} object
  stores a distribution function defined on a set of intervals.

}

\usage{

idf(left, right, p)
\method{print}{idf}(x, ...)

}

\arguments{
  
  \item{left, right}{left and right endpoints of intervals on which the
    distribution function is defined.}
  
  \item{p}{probabilities allocated to the intervals. Probability values
    will be normalized inside the function. }
  
  \item{x}{an object of class \code{idf}. }
  
  \item{...}{other arguments for printing. }
  
}

\details{

  When left and right endpoints are equal, the intervals reduce to
  points.
  
  \code{print.idf} prints an object of class \code{idf} as a three-coumn
  matrix.
  
}

\value{

  \item{left, right}{left and right endpoints of intervals on which the
  distribution function is defined.}
  
  \item{p}{probabilities allocated to the intervals.}

}

\author{ Yong Wang <yongwang@auckland.ac.nz>}

\seealso{ \code{\link{icendata}}, \code{\link{Deltamatrix}}, \code{\link{npsurv}}. }

\examples{
idf(1:5, 1:5*3-2, c(1,1,2,2,4))
npsurv(cbind(1:5, 1:5*3-2))$f    # NPMLE 
}

\keyword{ function }   % at least one, from doc/KEYWORDS
