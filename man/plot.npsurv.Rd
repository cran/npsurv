\name{plot.npsurv}
\alias{plot.npsurv}
\alias{plot.idf}
\alias{plotsurvidf}
\alias{plotgradidf}

\title{Plot Functions for Nonparametric Survival Estimation}

\description{
  
  \code{plot.npsurv} and \code{plot.idf} are wrapper functions that call
  either \code{plotsurvidf} or \code{plotgradidf}.

  \code{plotsurvidf} plots the survival function of the nonparametric
  maximum likelihood estimate (NPMLE). 

  \code{plotgradidf} plots the gradient function of the NPMLE.

}

\usage{

\method{plot}{npsurv}(x, ...)
\method{plot}{idf}(x, data, fn=c("surv","grad"), ...)
plotsurvidf(f, style=c("box","uniform","left","right","midpoint"),
            xlab="Time", ylab="Survival Probability", col="blue3", fill=NA,  
            add=FALSE, lty=1, lty.inf=2, xlim, ...)
plotgradidf(f, data, w=1, col1="red3", col2="blue3", 
            xlab="Survival Time", ylab="Gradient", xlim, ...)

}
	    
\arguments{
  
  \item{x}{an object of class \code{npsurv} (i.e., an output of function
    \code{npsurv}) or an object of class \code{idf}.}

  \item{fn}{either "surv" or "grad", to indicate plotting either the
    survival or the gradient function.}
  
  \item{f}{an object of class \code{idf}.}

  \item{style}{for how to plot the survival function on a "maximal
    intersection interval":

    = \code{box}, plot a rectangle, which shows the uncertainty of
    probability allocation within the interval;

    = \code{uniform}, treat it as a uniform distribution and hence the
    diagonal line of the rectangle is plotted;

    = \code{left}, plot only the left side of the rectangle;

    = \code{right}, plot only the right side of the rectangle;

    = \code{midpoint}, plot a vertical line at the midpoint of the interval.

  }
  
  \item{xlab, ylab}{x- or y-axis label. }
  
  \item{add}{= \code{TRUE}, adds the curve to the existing plot;

    = \code{FALSE}, plots the curve in a new one.}

  \item{col}{color for all line segments, including box/rectangle borders.}

  \item{fill}{color for filling a box/rectangle.}

  \item{lty}{line type}

  \item{lty.inf}{line type for the rectangle that may extend to infinity.}
  
  \item{data}{vector or matrix that stores observations, or an object of
  class \code{icendata}.}

  \item{w}{additional weights/multiplicities of the observations stored in
    \code{x}. }
  
  \item{col1}{color for drawing maximal intersection intervals allocated
    with positive probabilities.}

  \item{col2}{color for drawing all gradients and the maximal
    intersection intervals allocated with zero probabilities.}

  \item{xlim}{x-coordinate limit points.}

  \item{...}{arguments for other graphical parameters (see \code{par}).}
  
}

\details{

  \code{plotsurvidf} by default chooses a less saturated color for
  \code{fill} than \code{col}.
  
  \code{plotgradidf} plots gradient values as vertical lines located as
  the left endpoints of the maximal intersection intervals. Each maximal
  intersection interval is plotted as a wider line on the horizontal
  zero-gradient line, with a circle to represent the open left endpoint
  of the interval and a solid point the closed right endpoint of the
  interval. The maximal intersection intervals allocated with positive
  probabilities have zero gradients, and hence no vertical lines are
  drawn for them.
  
}

\author{ Yong Wang <yongwang@auckland.ac.nz>}

\references{

  Wang, Y. (2008). Dimension-reduced nonparametric maximum likelihood
  computation for interval-censored data. \emph{Computational Statistics &
  Data Analysis}, \bold{52}, 2388-2402.

}

\seealso{ \code{\link{icendata}}, \code{\link{idf}}, \code{\link{npsurv}}. }

\examples{
data(ap)
plot(r<-npsurv(ap))              # survival function
plot(r$f, ap, fn="g")            # all gradients virtually zeros.

data(cancer)
cancerRT = with(cancer, cancer[group=="RT",1:2])
plot(rt<-npsurv(cancerRT), xlim=c(0,60))                  # survival of RT 
cancerRCT = with(cancer, cancer[group=="RCT",1:2])
plot(rct<-npsurv(cancerRCT), add=TRUE, col="green3") # survival of RCT 
## as uniform dististrbutions.
plot(rt, add=TRUE, style="uniform", col="blue3")
plot(rct, add=TRUE, style="uniform", col="green3")

## plot gradients; must supply data
plot(rt, cancerRT, fn="g")        # for group RT
plotgradidf(rct$f, cancerRCT)   # or, for group RCT
}

\keyword{ function }   % at least one, from doc/KEYWORDS
