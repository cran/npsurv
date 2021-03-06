% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/npsurv.R
\name{idf}
\alias{idf}
\alias{idf.object}
\alias{print.idf}
\title{Interval Distribution Function}
\usage{
idf(left, right, p)
\method{print}{idf}(x, ...)
}
\arguments{
\item{left, right}{left and right endpoints of intervals on which the
distribution function is defined.}

\item{p}{probabilities allocated to the intervals. Probability values will
be normalized inside the function.}

\item{x}{an object of class \code{idf}.}

\item{...}{other arguments for printing.}
}
\value{
\item{left, right}{left and right endpoints of intervals on which the
distribution function is defined.}

\item{p}{probabilities allocated to the intervals.}
}
\description{
Class \code{idf} can be used to store a distribution function
defined on a set of intervals. There are several functions
associated with the class.
}
\details{
\code{idf} creates an object of class \code{idf}. An \code{idf} object
stores a distribution function defined on a set of intervals.


When left and right endpoints are identical, the intervals just
represent exact points.

\code{print.idf} prints an object of class \code{idf} as a three-coumn
matrix.
}
\examples{

idf(1:5, 1:5*3-2, c(1,1,2,2,4))
npsurv(cbind(1:5, 1:5*3-2))$f    # NPMLE 

}
\seealso{
\code{\link{icendata}}, \code{\link{Deltamatrix}},
\code{\link{npsurv}}.
}
\author{
Yong Wang <yongwang@auckland.ac.nz>
}
\keyword{function}
