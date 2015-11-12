\name{icendata}
\alias{icendata}
\alias{is.icendata}
\alias{icendata.object}

\title{Class of Interval-censored Data}

\description{
  
  Function \code{icendata} creates an object of class 'icendata', which
  can save interval-censored and exact observations.

  Function \code{is.icendata} simply checks if an object is of class
  'icendata'.

}

\usage{

 icendata(x, w=1)
 is.icendata(x)

}

\arguments{
  
  \item{x}{vector or matrix.}
  
  \item{w}{weights or multiplicities of the observations.}
  
}

\details{

  If \code{x} is a vector, it contains only exact observations, with
  weights given in \code{w}.

  If \code{x} is a two-column matrix, it contains interval-censored
  observations and stores their left and right end-points in the first
  and second column, respectively. If the left and right end-points are
  equal, then the observation is exact. Weights are provided by
  \code{w}.

  If \code{x} is a three-column matrix, it contains interval-censored
  observations and stores their left and right end-points in the first
  and second column, respectively. The weight of each observation is the
  third-column value multiplied by the corresponding weight value in
  \code{w}.

  It is useful to turn interval-censored (and exact) observations into
  the format imposed by \code{icendata} so that they can be processed in
  a standardized format. Also, exact and interval-censored observations
  are stored separately in this format and can hence be dealt with more
  easily. Most functions in the package \code{npsurv} first ensure that
  the data has this format before processing.

  Observations of zero weights are removed. Identical observations are
  aggregated.

}

\value{
  
  \item{t}{numeric vector, storing exact observations.}
  
  \item{wt}{numeric vector, storing the weights of exact observations.}
  
  \item{o}{two-column numeric matrix, storing interval-censored observations.}
  
  \item{wo}{numeric vector, storing the weights of interval-censored
  observations.}
  
  \item{i1}{logical vector, indicating whether exact observations are
    less than \code{upper}.}
  
  \item{upper}{the largest finite value of \code{t} and \code{o}.}

  \item{u}{numeric vector, containing 0 and all unique finite values in
    \code{t} and \code{o}.}
  
}

\author{ Yong Wang <yongwang@auckland.ac.nz> }

\references{

  Wang, Y. (2008). Dimension-reduced nonparametric maximum likelihood
  computation for interval-censored data. Computational Statistics &
  Data Analysis, 52, 2388-2402.

  Wang, Y. and S. Fani (2015). \emph{Nonparametric maximum likelihood
  computation of a U-shaped hazard function}. (In Preparation).
  
}

\seealso{ \code{\link{npsurv}}, \code{\link{Uhaz}}. }

\examples{
data(ap)
(x = icendata(ap))
is.icendata(x)

data(gastric)
icendata(gastric)

data(leukemia)
i = leukemia[,"group"] == "6-MP"
icendata(leukemia[i,1:2])
}

\keyword{ function }
