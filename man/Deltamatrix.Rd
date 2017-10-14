\name{Deltamatrix}
\alias{Deltamatrix}

\title{Delta matrix}

\description{
  
  \code{Deltamatrix} computes the Delta matrix, along with maximal
  intersection intervals, for a set of intervals.

}

\usage{

Deltamatrix(LR) 

}

\arguments{
  
  \item{LR}{two-column matrix, each row of which stores an censoring
  interval of the form \eqn{(L_i, R_i]}{(Li, Ri]}.  If \eqn{L_i =
  R_i}{Li = Ri}, it is an exact observation. }
  
}

\details{

  An intersection interval is a nonempty intersection of any combination
  of the given intervals, and a maximal intersection interval is an
  intersection interval that contains no other intersection interval.

  The Delta matrix is a matrix of indicators (\code{TRUE} or
  \code{FALSE}). The rows correspond to the given interval-censored
  observations, and the columns the maximal intersection intervals. A
  \code{TRUE} value of the (i,j)-th element means that the i-th
  observation covers the j-th maximal intersection interval, and a
  \code{FALSE} value means the opposite.

}

\value{

  A list with components:
  
  \item{left}{left endpoints of the maximal intersection intervals.}
  
  \item{right}{right endpoints of the maximal intersection intervals.}
  
  \item{Delta}{logical matrix, for the Delta matrix.}

}

\author{ Yong Wang <yongwang@auckland.ac.nz>}

\references{

  Wang, Y. (2008). Dimension-reduced nonparametric maximum likelihood
  computation for interval-censored data.  \emph{Computational
  Statistics & Data Analysis}, \bold{52}, 2388-2402.

}

\seealso{ \code{\link{icendata}}, \code{\link{idf}}. }

\examples{
(x = cbind(1:5,1:5*3-2))
Deltamatrix(x)
}
\keyword{ function }   % at least one, from doc/KEYWORDS
