\name{uh}
\alias{uh}
\alias{uh.object}
\alias{print.uh}

\title{U-shaped Hazard Function}
\description{
  
  \code{uh} creates an object of class \code{uh}, which stores a
  U-shaped hazard function. 

  \code{print.uh} prints an object of class \code{uh}.
  
}

\usage{

uh(alpha, tau, nu, eta, mu, upper=Inf, deg=1, collapse=TRUE)
\method{print}{uh}(x, ...)

}

\arguments{
  
  \item{alpha}{a nonnegative value, for the constant coefficient.}
  
  \item{tau}{vector of nonnegative real values, for left knots. }
  
  \item{nu}{vector of nonnegative values, for masses associated with the
    left knots. }
  
  \item{eta}{vector of nonnegative real values, for right knots. }
  
  \item{mu}{vector of nonnegative real values, for masses associated with the
    right knots. }

  \item{upper}{a positive value, at which point the hazard starts to
    become infinite. }

  \item{deg}{nonnegative real number for spline degree (i.e., p in the formula below).}

  \item{collapse}{logical, indicating if identical knots should be collapsed.}

  \item{x}{an object of class \code{uh}.}
  
  \item{...}{other auguments for printing.}
  
}

\details{

  A U-shape hazard function, as generalized by Wang and Fani (2017), is
  given by

  \deqn{h(t) = \alpha + \sum_{j = 1}^k \nu_j(\tau_j - t)_+^p + \sum_{j
      = 1}^{m} \mu_j (t-\eta_j)_+^p,}{     h(t) = alpha + sum_{j=1}^k nu_j (tau_j - t)_+^p 
                  + sum_{j=1}^m mu_j (t - eta_j)_+^p,}

  where \eqn{\alpha,\nu_j,\mu_j \ge 0}{alpha, nu_j, mu_j \ge 0},
  \eqn{\tau_1 < \cdots < \tau_k \le \eta_1 < \cdots < \eta_m,}{tau_1 < ... < tau_k <= eta_1 < ... < eta_m,} and \eqn{p \ge 0}{p >= 0} is the
  the spline degree which determines 
  the smoothness of the U-shaped hazard. As p increases, the
  family of hazard functions becomes increasingly smoother, but at the
  same time, smaller. When \eqn{p = 0}{p = 0}, the hazard function is
  U-shaped, as studied by Bray et al. (1967). When \eqn{p = 1}{p = 1},
  the hazard function is convex, as studied by Jankowski and Wellner
  (2009a,b).

  \code{print.uh} prints an object of class \code{uh}. While
  \code{alpha}, \code{upper} and \code{deg} are printed as they are,
  \code{tau} and \code{nu} are printed as a two-column matrix, and so
  are \code{eta} and \code{mu}.
  
}

\value{

  \code{uh} returns an object of class \code{uh}. It is a list with
  components \code{alpha}, \code{tau}, \code{nu}, \code{eta}, \code{mu},
  \code{upper} and \code{deg}, which store their corresponding values as
  described above.
  
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
  
  Wang, Y. and Fani, S. (2017). Nonparametric maximum likelihood
  computation of a U-shaped hazard function. \emph{Statistics and
  Computing}, (in print).

}

\seealso{ \code{\link{Uhaz}}, \code{\link{icendata}}, \code{\link{plot.uh}} }

\examples{
(h0 = uh(3, 2, 3, 4, 5, 7, deg=0))              # deg = 0
plot(h0, ylim=c(0,20))
(h1 = uh(4, 2, 3, 5, 6, 7, deg=1))              # deg = 1
plot(h1, add=TRUE, col="green3")
(h2 = uh(1, 1:2, 3:4, 5:6, 7:8, 9, deg=2))      # deg = 2
plot(h2, add=TRUE, col="red3")
}

\keyword{ function }   % at least one, from doc/KEYWORDS
