\name{marijuana}
\alias{marijuana}

\docType{data}

\title{
  Angina Pectoris Survival Data
}

\description{ Contains the answers of 191 California high school
  students to the question: "When did you first use marijuana?". An
  answer can be an exact age, or "I have never used it", which gives
  rise to a right-censored observation, or "I have used it but cannot
  recall just when the first time was", which gives rise to a
  left-censored observation. }

\usage{marijuana}

\format{ A data frame with 21 observations and 3 variables:
  
  L: left-end point of an interval-censored time;

  R: right-end point of an interval-censored time;

 count: number of students in the interval.}

\source{ Turnbull and Weiss (1978). See also Klein and Moeschberger
  (1997), page 17. }

\references{ Turnbull and Weiss (1978). A likelihood ratio statistic
  fortesting goodness of fit with randomly censored data. \emph{Biometrics},
  \bold{34}, 367-375.
  
  Klein and Moeschberger (2003). \emph{Survival Analysis: Techniques for
  Censored and Truncated Data} (2nd ed.). Springer
  
}

\seealso{ \code{\link{npsurv}}. }

\examples{
data(marijuana)
r = Uhaz(marijuana, deg=2)
plot(r$h, fn="h")
plot(r$h, fn="s")
}

\keyword{datasets}

