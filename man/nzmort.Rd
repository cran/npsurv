\name{nzmort}
\alias{nzmort}

\docType{data}

\title{
  New Zealand Mortality in 2000
}

\description{ Contains the number of deaths of Maori and Non-Maori
  people at each age in New Zealand in 2000. }

\usage{nzmort}

\format{ A data frame with 210 observations and 3 variables:
  
  age: at which age the deaths occurred;

  deaths: number of people died at the age;

  ethnic: either Maori or Non-Maori.}

\details{Data contains no age with zero death.}

\source{ \url{http://www.mortality.org/} }

\seealso{ \code{\link{Uhaz}}. }

\examples{
data(nzmort)
x = with(nzmort, nzmort[ethnic=="maori",])[,1:2]      # Maori mortality
# x = with(nzmort, nzmort[ethnic!="maori",])[,1:2]    # Non-Maori mortality

## As exact observations
# Plot hazard functions
h0 = Uhaz(x[,1]+0.5, x[,2], deg=0)$h    # U-shaped hazard
plot(h0, fn="h", col="green3", pch=2)
h1 = Uhaz(x[,1]+0.5, x[,2], deg=1)$h    # convex hazard
plot(h1, fn="h", add=TRUE, pch=1)
h2 = Uhaz(x[,1]+0.5, x[,2], deg=2)$h    # smooth U-shaped hazard
plot(h2, fn="h", add=TRUE, col="red3")

# Plot densities
age = 0:max(x[,1])
count = integer(length(age))
count[x[,"age"]+1] = x[,"deaths"]
barplot(count/sum(count), space=0, col="lightgrey", ylab="Density")
axis(1, pos=NA, at=0:10*10)
plot(h0, fn="d", add=TRUE, col="green3", pch=2)
plot(h1, fn="d", add=TRUE, col="blue3", pch=1)
plot(h2, fn="d", add=TRUE, col="red3", pch=19)

## As interval-censored observations
# Plot hazard functions
x2 = cbind(x[,1], x[,1]+1, x[,2])
h0 = Uhaz(x2, deg=0)$h      # U-shaped hazard
plot(h0, fn="h", col="green3", pch=2)
h1 = Uhaz(x2, deg=1)$h      # convex hazard
plot(h1, fn="h", add=TRUE, pch=1)
h2 = Uhaz(x2, deg=2)$h      # smooth U-shaped hazard
plot(h2, fn="h", add=TRUE, col="red3", pch=1)

# Plot densities
barplot(count/sum(count), space=0, col="lightgrey")
axis(1, pos=NA, at=0:10*10)
plot(h0, fn="d", add=TRUE, col="green3", pch=2)
plot(h1, fn="d", add=TRUE, col="blue3", pch=1)
plot(h2, fn="d", add=TRUE, col="red3", pch=19)
}

\keyword{datasets}

