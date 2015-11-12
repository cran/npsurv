\name{Uhaz}
\alias{Uhaz}
\alias{Uhaz.object}

\title{U-shaped Hazard Function Estimation}
\description{
  
  \code{Uhaz} computes the nonparametric maximum likelihood esimate
  (NPMLE) of a U-shaped hazard function from exact or interval-censored
  data, or a mix of the two types of data.

}

\usage{

Uhaz(data, w=1, deg=1, maxit=100, tol=1e-6, verb=0)

}

\arguments{
  
  \item{data}{vector or matrix, or an object of class \code{icendata}.}
  
  \item{w}{weights or multiplicities of the observations.}

  \item{deg}{nonnegative real number for spline degree (i.e., p in the
  formula below).}

  \item{maxit}{maximum number of iterations.}

  \item{tol}{ tolerance level for stopping the algorithm. It is used as
    the threshold on the increase of the log-likelihood after each
    iteration. }

  \item{verb}{ verbosity level for printing intermediate results in each
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

  The algorithm used for the computing the NPMLE of a hazard function
  under the U-shape restriction is is proposed by Wang and Fani
  (2015). Such a hazard function is given by

  \deqn{h(t) = \alpha + \sum_{j = 1}^k \nu_j(\tau_j - t)_+^p + \sum_{j
      = 1}^{m} \mu_j (t-\eta_j)_+^p,}{     h(t) = alpha + sum_{j=1}^k nu_j (tau_j - t)_+^p 
                  + sum_{j=1}^m mu_j (t - eta_j)_+^p,}

  where \eqn{\alpha,\nu_j,\mu_j \ge 0}{alpha, nu_j, mu_j \ge 0},
  \eqn{\tau_1 < \cdots < \tau_k \le \eta_1 < \cdots < \eta_m,}{tau_1 < ... < tau_k <= eta_1 < ... < eta_m,} and \eqn{p \ge 0}{p >= 0} is the
  the spline degree which determines 
  the smoothness of the U-shaped hazard. As p increases, the
  family of hazard functions becomes increasingly smoother, but at the
  time, smaller. When p = 0, the hazard function is U-shaped, as studied
  by Bray et al. (1967). When p = 1, the hazard function is convex, as
  studied by Jankowski and Wellner (2009a,b).

  Note that \code{deg} (i.e., p in the above mathematical display) can
  take on any nonnegative real value.
  
}

\value{

  An object of class \code{Uhaz}, which is a list consisting of:
  
  \item{convergence}{= \code{TRUE}, converged successfully;

    = \code{FALSE}, maximum number of iterations reached.}
  
  \item{grad}{gradient values at the knots.}
  
  \item{numiter}{number of iterations used.}
  
  \item{ll}{log-likelihood value of the NPMLE \code{h}.}
  
  \item{h}{NPMLE of the U-shaped hazard function, an object of class
    \code{uh}.}
  
}

\author{ Yong Wang <yongwang@auckland.ac.nz>}

\references{

  Bray, T. A., Crawford, G. B., and Proschan, F. (1967).  \emph{Maximum
  Likelihood Estimation of a U-shaped Failure Rate Function}. Defense
  Technical Information Center.
  
  Jankowski, H. K. and Wellner, J. A. (2009a). Computation of
  nonparametric convex hazard estimators via profile
  methods. \emph{Journal of Nonparametric Statistics}, \bold{21},
  505-518.

  Jankowski, H. K. and Wellner, J. A. (2009b). Nonparametric estimation
  of a convex bathtub-shaped hazard function. \emph{Bernoulli},
  \bold{15}, 1010-1035.
  
  Wang, Y. and Fani, S. (2015). Nonparametric maximum likelihood
  computation of a U-shaped hazard function. (In preparation)

}

\seealso{ \code{\link{icendata}}, \code{\link{nzmort}}. }

\examples{
## Interval-censored observations
data(ap)
(r = Uhaz(ap, deg=0))
plot(r, ylim=c(0,.3))
for(i in 1:6) plot(Uhaz(ap, deg=i/2), add=TRUE, col=i+1)

## Exact observations
data(nzmort)
x = with(nzmort, nzmort[ethnic=="maori",])[,1:2]   # Maori mortality
(h0 = Uhaz(x[,1]+0.5, x[,2], deg=0)$h)     # U-shaped hazard
(h1 = Uhaz(x[,1]+0.5, x[,2], deg=1)$h)     # convex hazard
(h2 <- Uhaz(x[,1]+0.5, x[,2], deg=2)$h)    # smooth U-shaped hazard

plot(h0, pch=2)              # plot hazard functions
plot(h1, add=TRUE, col="green3", pch=1)
plot(h2, add=TRUE, col="red3", pch=19)

age = 0:max(x[,1])                         # plot densities
count = integer(length(age))
count[x[,"age"]+1] = x[,"deaths"]
barplot(count/sum(count), space=0, col="lightgrey")
axis(1, pos=NA, at=0:10*10)
plot(h0, fn="d", add=TRUE, pch=2)
plot(h1, fn="d", add=TRUE, col="green3", pch=1)
plot(h2, fn="d", add=TRUE, col="red3", pch=19)

plot(h0, fn="s", pch=2)      # plot survivals
plot(h1, fn="s", add=TRUE, col="green3", pch=1)
plot(h2, fn="s", add=TRUE, col="red3", pch=19)

## Exact and right-censored observations
data(gastric)
plot(h0<-Uhaz(gastric, deg=0)$h)           # plot hazard functions
plot(h1<-Uhaz(gastric, deg=1)$h, add=TRUE, col="green3")
plot(h2<-Uhaz(gastric, deg=2)$h, add=TRUE, col="red3")

plot(npsurv(gastric), fn="s", col="grey")  # plot survival functions
plot(h0, fn="s", add=TRUE)
plot(h1, fn="s", add=TRUE, col="green3")
plot(h2, fn="s", add=TRUE, col="red3")
}

\keyword{ function }   % at least one, from doc/KEYWORDS
