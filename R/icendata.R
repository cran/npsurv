######################## #
# Interval-censored data #
######################## #

# Also allows for exact observations included.

##'Class of Interval-censored Data
##'
##'
##' Class \code{icendata} can be used to store general
##' interval-censored data, which may possibly contain exact
##' observations.There are several functions associated with the
##' class.
##' 
##'Function \code{icendata} creates an object of class 'icendata', which can be
##'used to save both interval-censored and exact observations.
##'
##'Function \code{is.icendata} simply checks if an object is of class
##''icendata'.
##'
##'
##'If \code{x} is a vector, it contains only exact observations, with weights
##'given in \code{w}.
##'
##'If \code{x} is a two-column matrix, it contains interval-censored
##'observations and stores their left and right endpoints in the first and
##'second column, respectively. If the left and right endpoints are equal, then
##'the observation is exact. Weights are provided by \code{w}.
##'
##'If \code{x} is a three-column matrix, it contains interval-censored
##'observations and stores their left and right endpoints in the first and
##'second column, respectively. The weight of each observation is the
##'third-column value multiplied by the corresponding weight value in \code{w}.
##'
##'It is useful to turn interval-censored (and exact) observations into the
##'format imposed by \code{icendata} so that they can be processed in a
##'standardized format by other functions. Also, exact and interval-censored
##'observations are stored separately in this format and can hence be dealt
##'with more easily. Most functions in the package \code{npsurv} first ensure
##'that the data has this format before processing.
##'
##'Observations of zero weights are removed. Identical observations are
##'aggregated.
##'
##'An interval-valued observation is either \eqn{(L_i, R_i]}{(Li, Ri]} if
##'\eqn{L_i < R_i}{Li < Ri}, or \eqn{[L_i, R_i]}{[Li, Ri]} if \eqn{L_i =
##'R_i}{Li = Ri}.
##'
##'@aliases icendata is.icendata icendata.object
##'@param x vector or matrix.
##'@param w weights or multiplicities of the observations.
##'@return
##'
##'\item{t}{numeric vector, storing exact observations.}
##'
##'\item{wt}{numeric vector, storing the weights of exact observations.}
##'
##'\item{o}{two-column numeric matrix, storing interval-censored observations.}
##'
##'\item{wo}{numeric vector, storing the weights of interval-censored
##'observations.}
##'
##'\item{i1}{logical vector, indicating whether exact observations are less
##'than \code{upper}.}
##'
##'\item{upper}{the largest finite value of \code{t} and \code{o}.}
##'
##'\item{u}{numeric vector, containing 0 and all unique finite values in
##'\code{t} and \code{o}.}
##'@author Yong Wang <yongwang@@auckland.ac.nz>
##'@seealso \code{\link{npsurv}}, \code{\link{Uhaz}}.
##'@references
##'
##'Wang, Y. (2008). Dimension-reduced nonparametric maximum likelihood
##'computation for interval-censored data. Computational Statistics & Data
##'Analysis, 52, 2388-2402.
##'
##'Wang, Y. and Fani, S. (2017). Nonparametric maximum likelihood computation
##'of a U-shaped hazard function. \emph{Statistics and Computing}, (in print).
##'@keywords class function
##'@examples
##'
##'data(ap)
##'(x = icendata(ap))
##'is.icendata(x)
##'
##'data(gastric)
##'icendata(gastric)
##'
##'data(leukemia)
##'i = leukemia[,"group"] == "6-MP"
##'icendata(leukemia[i,1:2])
##'
##'@usage
##'icendata(x, w=1)
##'is.icendata(x)
##' 
##'@export icendata
##'@export is.icendata

icendata = function(x, w=1) {
  if(is.null(x)) return(NULL)
  if(is.icendata(x)) {
    if(all(w == 1)) return(x)
    w = rep(w, length = length(x$t) + nrow(x$o))
    if(length(x$t) > 0) x$wt = x$wt * w[1:length(x$wt)]
    if(nrow(x$o) > 0) x$wo = x$wo * w[length(x$wt)+1:nrow(x$o)]
    return(x)
  }
  z = vector("list", 7)
  names(z) = c("t", "wt", "o", "wo", "i1", "upper", "u")
  if(is.vector(x)) x = cbind(x, x)
  if(!is.matrix(x)) x = as.matrix(x)
  if(ncol(x) == 3) {w = w * x[,3]; x = x[,1:2]}
  if(length(w) != nrow(x)) w = rep(w, len=nrow(x))
  iw = w > 0
  w = w[iw]
  x = x[iw,,drop=FALSE]
  o = order(x[,1], x[,2])
  x = x[o,]
  w = w[o]
  id = c(TRUE, diff(x[,1]) > 0 | diff(x[,2]) > 0)
  id[is.na(id)] = FALSE            # for Inf's
  w = aggregate(w, by=list(group=cumsum(id)), sum)[,2]
  x = x[id,]
  i = x[,1] == x[,2]
  z$t = x[i,1]
  names(z$t) = NULL
  z$wt = w[i]
  z$o = x[!i,1:2,drop=FALSE]
  dimnames(z$o) = list(NULL, c("L","R"))
  z$wo = w[!i]
  z$upper = max(x[,1])
  z$i1 = z$t != z$upper
  z$u = sort(unique(c(0, pmin(c(x[,1], x[,2]), z$upper))))
  class(z) = "icendata"
  z
}

is.icendata = function(x) "icendata" %in% class(x)

# is.rightcensored.icendata = function(x) all(x$o[,2] == Inf)

expand.icendata = function(x) {
  if(!is.icendata(x)) x = icendata(x)
  z = vector("list", 7)
  names(z) = c("t", "wt", "o", "wo", "i1", "upper", "u")
  z$upper = x$upper
  if(length(x$t) > 0) {
    z$t = rep(x$t, x$wt)
    z$wt = rep(1, length(z$t))
    z$i1 = z$t != z$upper
  }
  else z$t = z$wt = numeric(0)
  if(nrow(x$o) > 0) {
    z$o = cbind(rep(x$o[,1], x$wo), rep(x$o[,2], x$wo))
    z$wo = rep(1, nrow(z$o))
    colnames(z$o) = c("L","R")
  }
  else {z$o = matrix(nrow=0, ncol=2); z$wo = numeric(0)}
  z$u = x$u
  class(z) = "icendata"
  z
}

