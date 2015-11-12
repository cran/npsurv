\name{npsurv}
\alias{npsurv}
\alias{npsurv.object}

\title{Nonparametric Survival Function Estimation}
\description{
  
  \code{npsurv} computes the nonparametric maximum likelihood esimate
  (NPMLE) of a survival function for general interval-censored data.

}

\usage{

npsurv(data, w=1, maxit=100, tol=1e-6, verb=0)

}

\arguments{
  
  \item{data}{vector or matrix, or an object of class \code{icendata}.}
  
  \item{w}{weights or multiplicities of the observations.}

  \item{maxit}{maximum number of iterations.}

  \item{tol}{ tolerance level for stopping the algorithm. It is used as
    the threshold on the increase of the log-likelihood after each
    iteration. }

  \item{verb}{ verbosity level for printing intermediate results at each
    iteration. }

}

\details{

  If \code{data} is a vector, it contains only exact observations, with
  weights given in \code{w}.

  If \code{data} is a matrix with two columns, it contains
  interval-censored observations, with the two columns storing their
  left and right end-points, respectively. If the left and right
  end-points are equal, then the observation is exact. Weights are
  provided by \code{w}.

  If \code{data} is a matrix with three columns, it contains
  interval-censored observations, with the first two columns storing
  their left and right end-points, respectively. The weight of each
  observation is the third-column value multiplied by the corresponding
  weight value in \code{w}.

  The algorithm used for the computing the NPMLE is the constrained
  Newton method (CNM) that is proposed by Wang (2008).

  Inside the function, it detects if data has only right censoring, and
  if so, the Kaplan-Meier estimate is computed directly by function
  \code{km}.

}

\value{

  An object of class \code{npsurv}, which is a list consisting of:
  
  \item{f}{NPMLE, an object of class \code{idf}.}
  
  \item{upper}{largest finite value in the data.}
  
  \item{convergence}{= \code{TRUE}, converged successfully;

    = \code{FALSE}, maximum number of iterations reached.}
  
  \item{ll}{log-likelihood value of the NPMLE \code{f}.}
  
  \item{maxgrad}{maximum gradient value of the NPMLE \code{f}.}
  
  \item{numiter}{number of iterations used.}
  
}

\author{ Yong Wang <yongwang@auckland.ac.nz>}

\references{

  Wang, Y. (2008). Dimension-reduced nonparametric maximum likelihood
  computation for interval-censored data. \emph{Computational Statistics &
  Data Analysis}, \bold{52}, 2388-2402.

}

\seealso{ \code{\link{icendata}}, \code{\link{Deltamatrix}},
  \code{\link{idf}}, \code{\link{km}}. }

\examples{
## all exact observations
data(acfail)
plot(npsurv(acfail))

## right-censored (and exact) observations
data(gastric)
plot(npsurv(gastric))

data(leukemia)
i = leukemia[,"group"] == "Placebo"
plot(npsurv(leukemia[i,1:2]), xlim=c(0,40), col="blue") # placebo
plot(npsurv(leukemia[!i,1:2]), add=TRUE, col="red")     # 6-MP

## purely interval-censored data
data(ap)
plot(npsurv(ap))

data(cancer)
cancerRT = with(cancer, cancer[group=="RT",1:2])
plot(npsurv(cancerRT), xlim=c(0,60))                  # survival of RT 
cancerRCT = with(cancer, cancer[group=="RCT",1:2])
plot(npsurv(cancerRCT), add=TRUE, col="green")        # survival of RCT 
}
\keyword{ function }
