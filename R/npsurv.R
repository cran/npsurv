# ----------------------------------------------------------------------- #
# Nonparametric maximum likelihood estimation from interval-censored data #
# ----------------------------------------------------------------------- #



##'Nonparametric Survival Function Estimation
##'
##'
##'\code{npsurv} computes the nonparametric maximum likelihood esimate (NPMLE)
##'of a survival function for general interval-censored data.
##'
##'
##'If \code{data} is a vector, it contains only exact observations, with
##'weights given in \code{w}.
##'
##'If \code{data} is a matrix with two columns, it contains interval-censored
##'observations, with the two columns storing their left and right end-points,
##'respectively. If the left and right end-points are equal, then the
##'observation is exact. Weights are provided by \code{w}.
##'
##'If \code{data} is a matrix with three columns, it contains interval-censored
##'observations, with the first two columns storing their left and right
##'end-points, respectively. The weight of each observation is the third-column
##'value multiplied by the corresponding weight value in \code{w}.
##'
##'The algorithm used for computing the NPMLE is either the constrained Newton
##'method (CNM) (Wang, 2008), or the hierachical constrained Newton method
##'(HCNM) (Wang and Taylor, 2013) when there are a large number of maximal
##'intersection intervals.
##'
##'Inside the function, it examines if data has only right censoring, and if
##'so, the Kaplan-Meier estimate is computed directly by function \code{km}.
##'
##'An interval-valued observation is either \eqn{(L_i, R_i]}{(Li, Ri]} if
##'\eqn{L_i < R_i}{Li < Ri}, or \eqn{[L_i, R_i]}{[Li, Ri]} if \eqn{L_i =
##'R_i}{Li = Ri}.
##'
##'@aliases npsurv npsurv.object
##'@param data vector or matrix, or an object of class \code{icendata}.
##'@param w weights or multiplicities of the observations.
##'@param maxit maximum number of iterations.
##'@param tol tolerance level for stopping the algorithm. It is used as the
##'threshold on the increase of the log-likelihood after each iteration.
##'@param verb verbosity level for printing intermediate results at each
##'iteration.
##'@return
##'
##'An object of class \code{npsurv}, which is a list with components:
##'
##'\item{f}{NPMLE, an object of class \code{idf}.}
##'
##'\item{upper}{largest finite value in the data.}
##'
##'\item{convergence}{= \code{TRUE}, converged successfully;
##'
##'= \code{FALSE}, maximum number of iterations reached.}
##'
##'\item{method}{method used internally, either \code{cnm} or \code{hcnm}.}
##'
##'\item{ll}{log-likelihood value of the NPMLE \code{f}.}
##'
##'\item{maxgrad}{maximum gradient value of the NPMLE \code{f}.}
##'
##'\item{numiter}{number of iterations used.}
##'@author Yong Wang <yongwang@@auckland.ac.nz>
##'@seealso \code{\link{icendata}}, \code{\link{Deltamatrix}},
##'\code{\link{idf}}, \code{\link{km}}.
##'@references
##'
##'Wang, Y. (2008). Dimension-reduced nonparametric maximum likelihood
##'computation for interval-censored data. \emph{Computational Statistics &
##'Data Analysis}, \bold{52}, 2388-2402.
##'
##'Wang, Y. and Taylor, S. M. (2013). Efficient computation of nonparametric
##'survival functions via a hierarchical mixture formulation. \emph{Statistics
##'and Computing}, \bold{23}, 713-725.
##'@keywords function
##'@examples
##'
##'## all exact observations
##'data(acfail)
##'plot(npsurv(acfail))
##'
##'## right-censored (and exact) observations
##'data(gastric)
##'plot(npsurv(gastric))
##'
##'data(leukemia)
##'i = leukemia[,"group"] == "Placebo"
##'plot(npsurv(leukemia[i,1:2]), xlim=c(0,40), col="blue") # placebo
##'plot(npsurv(leukemia[!i,1:2]), add=TRUE, col="red")     # 6-MP
##'
##'## purely interval-censored data
##'data(ap)
##'plot(npsurv(ap))
##'
##'data(cancer)
##'cancerRT = with(cancer, cancer[group=="RT",1:2])
##'plot(npsurv(cancerRT), xlim=c(0,60))                  # survival of RT 
##'cancerRCT = with(cancer, cancer[group=="RCT",1:2])
##'plot(npsurv(cancerRCT), add=TRUE, col="green")        # survival of RCT 
##'
##'@export npsurv
npsurv = function(data, w=1, maxit=100, tol=1e-6, verb=0) {
  x2 = icendata(data, w)
  if(nrow(x2$o) == 0 || all(x2$o[,2] == Inf)) { # exact or right-censored only
    r0 = km(x2)
    r = list(f=r0$f, upper=max(x2$t, x2$o[,1]), convergence=TRUE, ll=r0$ll,
        maxgrad=0, numiter=1)
    return(structure(r, class="npsurv"))
  }
  x = rbind(cbind(x2$t, x2$t), x2$o)
  nx = nrow(x)
  w = c(x2$wt, x2$wo)
  wr = sqrt(w)
  n = sum(w)
  upper = x2$upper
  dmat = Deltamatrix(x)
  left = dmat$left
  right = dmat$right
  D = dmat$Delta
  m = length(left)
  p = double(m)
  i = rowSums(D) != 1
  j = colSums(D[!i,,drop=FALSE]) > 0
  j[c(1,m)] = TRUE
  repeat {                 # Initial p must ensure P > 0
    jm = which.max(colSums(D[i,,drop=FALSE]))
    j[jm] = TRUE
    i[D[,jm]] = FALSE
    if( sum(i) == 0 ) break
  }
  p = colSums(w * D) * j
  p = p / sum(p)
  if(m >= 200) {                     ## Turn to HCNM
    r = hcnm(w=w, D=D, p0=p, maxit=maxit, tol=tol, verb=verb)
    j = r$pf > 0
    f = idf(left[j], right[j], r$pf[j]) 
    r = list(f=f, upper=upper, convergence=r$convergence, method="hcnm", ll=r$ll,
             maxgrad=r$maxgrad, numiter=r$numiter)
    return(structure(r, class="npsurv"))
  }
  
  P = drop(D %*% p)
  ll = sum( w * log(P) )
  converge = FALSE
  for(i in 1:maxit) {
    p.old = p
    ll.old = ll
    S = D / pmax(P, 1e-100)
    d = colSums(w * S)
    dmax = max(d) - n
    if(verb > 0) {
      cat("##### Iteration", i, "#####\n")
      cat("Log-likelihood: ", signif(ll, 6), "\n")
    }
    if(verb > 1) cat("Maximum gradient: ", signif(dmax, 6), "\n")
    if(verb > 2) {cat("Probability vector:\n"); print(p)} 
    j[which(j)-1 + aggregate(d, by=list(group=cumsum(j)), which.max)[,2]] = TRUE
    pj = pnnls(wr * S[,j,drop=FALSE], wr * 2, sum=1)$x
    p[j] = pj / sum(pj)
    alpha = 1                # line search
    pd = p - p.old
    lld = sum(d * pd)
    p.alpha = p
    repeat {
      P.alpha = drop(D %*% p.alpha)
      ll.alpha = sum(w * log(P.alpha))
      if(ll.alpha >= ll + alpha * lld * .33)
        { p = p.alpha; P = P.alpha; ll = ll.alpha; break }
      if((alpha <- alpha * .5) < 1e-10) break
      p.alpha = p.old + alpha * pd
    }
    j = p > 0
    if( ll <= ll.old + tol ) {converge=TRUE; break}
  }
  f = idf(left[j], right[j], p[j])
  r = list(f=f, upper=upper, convergence=converge, method="cnm", ll=ll,
      maxgrad=max(crossprod(w/P, D))-n, numiter=i)
  structure(r, class="npsurv")
}

# LR    matrix of intervals

# An interval is either (Li, Ri] if Li < Ri, or [Li, Ri] if Li = Ri. 



##'Delta matrix
##'
##'
##'\code{Deltamatrix} computes the Delta matrix, along with maximal
##'intersection intervals, for a set of intervals.
##'
##'
##'An intersection interval is a nonempty intersection of any combination of
##'the given intervals, and a maximal intersection interval is an intersection
##'interval that contains no other intersection interval.
##'
##'The Delta matrix is a matrix of indicators (\code{TRUE} or \code{FALSE}).
##'The rows correspond to the given interval-censored observations, and the
##'columns the maximal intersection intervals. A \code{TRUE} value of the
##'(i,j)-th element means that the i-th observation covers the j-th maximal
##'intersection interval, and a \code{FALSE} value means the opposite.
##'
##'@param LR two-column matrix, each row of which stores an censoring interval
##'of the form \eqn{(L_i, R_i]}{(Li, Ri]}.  If \eqn{L_i = }{Li = Ri}\eqn{
##'R_i}{Li = Ri}, it is an exact observation.
##'@return
##'
##'A list with components:
##'
##'\item{left}{left endpoints of the maximal intersection intervals.}
##'
##'\item{right}{right endpoints of the maximal intersection intervals.}
##'
##'\item{Delta}{logical matrix, for the Delta matrix.}
##'@author Yong Wang <yongwang@@auckland.ac.nz>
##'@seealso \code{\link{icendata}}, \code{\link{idf}}.
##'@references
##'
##'Wang, Y. (2008). Dimension-reduced nonparametric maximum likelihood
##'computation for interval-censored data.  \emph{Computational Statistics &
##'Data Analysis}, \bold{52}, 2388-2402.
##'@keywords function
##'@examples
##'
##'(x = cbind(1:5,1:5*3-2))
##'Deltamatrix(x)
##'
##'@export Deltamatrix
Deltamatrix = function(LR) {
  L = LR[,1]
  R = LR[,2]
  ic = L != R             # inverval-censored
  nc = sum(ic)
  # tol = max(R[R!=Inf]) * 1e-8
  if(nc > 0) {
    L1 = L[ic] + max(R[R!=Inf]) * 1e-8       # open left endpoints
    LRc = cbind(c(L1, R[ic]), c(rep(0,nc), rep(1,nc)), rep(1:nc, 2))
    LRc.o = LRc[order(LRc[,1]),]
    j = which(diff(LRc.o[,2]) == 1)
    left = L[ic][LRc.o[j,3]]
    right = R[ic][LRc.o[j+1,3]]
  }
  else left = right = numeric(0)
  if(nrow(LR) - nc > 0) {
    ut = unique(L[!ic])
    jin = colSums(outer(ut, left, ">") & outer(ut, right, "<=")) > 0
    left = c(ut, left[!jin])     # remove those that contain exact obs.
    right = c(ut, right[!jin])
    o = order(left, right)
    left = left[o]
    right = right[o]
  }
  ## D = outer(L, left, "<=") & outer(R, right, ">=") 
  D = outer(L, left, "<=") & outer(R, right, ">=") &
    (outer(L, right, "<") | outer(R, left, "=="))  

  dimnames(D) = names(left) = names(right) = NULL
  list(left=left, right=right, Delta=D)
}

# interval distribution function, i.e., a distribution function defined on
# a set of intervals.

# left      Left endpoints of the intervals
# right     Right endpoints of the intervals
# p         Probability masses allocated to the intervals



##'Interval Distribution Function
##'
##'
##' Class \code{idf} can be used to store a distribution function
##' defined on a set of intervals. There are several functions
##' associated with the class.
##'
##' \code{idf} creates an object of class \code{idf}. An \code{idf} object
##'stores a distribution function defined on a set of intervals.
##'
##'
##'When left and right endpoints are identical, the intervals just
##' represent exact points.
##'
##'\code{print.idf} prints an object of class \code{idf} as a three-coumn
##'matrix.
##'
##'@aliases idf idf.object print.idf
##'@param left,right left and right endpoints of intervals on which the
##'distribution function is defined.
##'@param p probabilities allocated to the intervals. Probability values will
##'be normalized inside the function.
##'@param x an object of class \code{idf}.
##'@param ... other arguments for printing.
##'@return
##'
##'\item{left, right}{left and right endpoints of intervals on which the
##'distribution function is defined.}
##'
##'\item{p}{probabilities allocated to the intervals.}
##'@author Yong Wang <yongwang@@auckland.ac.nz>
##'@seealso \code{\link{icendata}}, \code{\link{Deltamatrix}},
##'\code{\link{npsurv}}.
##'@keywords function
##'@examples
##'
##'idf(1:5, 1:5*3-2, c(1,1,2,2,4))
##'npsurv(cbind(1:5, 1:5*3-2))$f    # NPMLE 
##'
##'@usage
##'idf(left, right, p)
##'\method{print}{idf}(x, ...)
##' 
##'@export idf
##'@export print.idf
idf = function(left, right, p) {
  if(length(left) != length(right)) stop("length(left) != length(right)")
  names(left) = names(right) = names(p) = NULL
  p = rep(p, length=length(left))
  f = list(left=left, right=right, p=p/sum(p))
  structure(f, class="idf")
}

print.idf = function(x, ...) {
  print(cbind(left=x$left, right=x$right, p=x$p), ...)
}

# Kaplan-Meier estimate of the survival function for right-censored data



##'Kaplan-Meier Estimation
##'
##'
##'\code{km} computes the nonparametric maximum likelihood esimate (NPMLE) of a
##'survival function for right-censored data.
##'
##'
##'For details about the arguments, see \code{icendata}.
##'
##'@param data vector or matrix, or an object of class \code{icendata}.
##'@param w weights/multiplicities of observations.
##'@return
##'
##'A list with components:
##'
##'\item{f}{NPMLE, an object of class \code{idf}.}
##'
##'\item{ll}{log-likelihood value of the NPMLE \code{f}.}
##'@author Yong Wang <yongwang@@auckland.ac.nz>
##'@seealso \code{\link{icendata}}, \code{\link{npsurv}}, \code{\link{idf}}.
##'@references
##'
##'Kaplan, E. L. and Meier, P. (1958). Nonparametric estimation from incomplete
##'observations. \emph{Journal of the American Statistical Association},
##'\bold{53}, 457-481.
##'@keywords function
##'@examples
##'
##'x = cbind(1:5, c(1,Inf,3,4,Inf))
##'(f = km(x)$f)
##'plot(f)
##'
##'data(leukemia)
##'i = leukemia[,"group"] == "Placebo"
##'plot(km(leukemia[i,1:2])$f, xlim=c(0,40), col="green3") # placebo
##'plot(km(leukemia[!i,1:2])$f, add=TRUE)                  # 6-MP
##'
##'@export km
km = function(data, w=1) {
  x = icendata(data, w)
  if(any(x$o[,2] != Inf))
    stop("Not all observations are exact or right-censored")
  if(nrow(x$o) == 0) {              # no right-censored observations
    f = idf(x$t, x$t, x$wt)
    ll = sum(x$wt * log(f$p))
    return(list(f=f, ll=ll))
  }
  c = colSums(x$wo * outer(x$o[,1], x$t, "<"))
  n = sum(x$wt, x$wo)                            # number of observations
  r = n - c - c(0,cumsum(x$wt))[1:length(x$t)]   # no. at risk
  S = cumprod(1 - x$wt/r)                        # survival prob.
  # tab = cbind(x$t, x$wt, c, r, S)
  p = rev(diff(rev(c(1,S,0))))
  dc = x$wt + c
  if(max(x$t) > max(x$o[,1])) {
    f = idf(x$t, x$t, p[-length(p)])
    ll = sum( x$wt * log(f$p) )
  }
  else {
    f = idf(c(x$t,max(x$o[,1])), c(x$t,Inf), p)
    ll = sum(c(x$wt, n - sum(x$wt)) * log(f$p))
  }
  list(f=f, ll=ll)
}

####  Plot functions



##'Plot Functions for Nonparametric Survival Estimation
##'
##'Functions for plotting nonparametric survival functions and related ones.
##' 
##'\code{plot.npsurv} and \code{plot.idf} are wrapper functions that call
##'either \code{plotsurvidf} or \code{plotgradidf}.
##'
##'\code{plotsurvidf} plots the survival function of the nonparametric maximum
##'likelihood estimate (NPMLE).
##'
##'\code{plotgradidf} plots the gradient function of the NPMLE.
##'
##'
##'\code{plotsurvidf} by default chooses a less saturated color for \code{fill}
##'than \code{col}.
##'
##'\code{plotgradidf} plots gradient values as vertical lines located as the
##'left endpoints of the maximal intersection intervals. Each maximal
##'intersection interval is plotted as a wider line on the horizontal
##'zero-gradient line, with a circle to represent the open left endpoint of the
##'interval and a solid point the closed right endpoint of the interval. The
##'maximal intersection intervals allocated with positive probabilities have
##'zero gradients, and hence no vertical lines are drawn for them.
##'
##'@aliases plot.npsurv plot.idf plotsurvidf plotgradidf
##'@param x an object of class \code{npsurv} (i.e., an output of function
##'\code{npsurv}) or an object of class \code{idf}.
##'@param fn either "surv" or "grad", to indicate plotting either the survival
##'or the gradient function.
##'@param f an object of class \code{idf}.
##'@param style for how to plot the survival function on a "maximal
##'intersection interval":
##'
##'= \code{box}, plot a rectangle, which shows the uncertainty of probability
##'allocation within the interval;
##'
##'= \code{uniform}, treat it as a uniform distribution and hence the diagonal
##'line of the rectangle is plotted;
##'
##'= \code{left}, plot only the left side of the rectangle;
##'
##'= \code{right}, plot only the right side of the rectangle;
##'
##'= \code{midpoint}, plot a vertical line at the midpoint of the interval.
##'
##'@param xlab,ylab x- or y-axis label.
##'@param add = \code{TRUE}, adds the curve to the existing plot;
##'
##'= \code{FALSE}, plots the curve in a new one.
##'@param col color for all line segments, including box/rectangle borders.
##'@param fill color for filling a box/rectangle. By default, a lighter
##'semi-transparent color is used.
##'@param lty line type
##'@param lty.inf line type for the rectangle that may extend to infinity.
##'@param data vector or matrix that stores observations, or an object of class
##'\code{icendata}.
##'@param w additional weights/multiplicities of the observations stored in
##'\code{x}.
##'@param col1 color for drawing maximal intersection intervals allocated with
##'positive probabilities.
##'@param col2 color for drawing all gradients and the maximal intersection
##'intervals allocated with zero probabilities.
##'@param xlim x-coordinate limit points.
##'@param ... arguments for other graphical parameters (see \code{par}).
##'@author Yong Wang <yongwang@@auckland.ac.nz>
##'@seealso \code{\link{icendata}}, \code{\link{idf}}, \code{\link{npsurv}}.
##'@references
##'
##'Wang, Y. (2008). Dimension-reduced nonparametric maximum likelihood
##'computation for interval-censored data. \emph{Computational Statistics &
##'Data Analysis}, \bold{52}, 2388-2402.
##'@keywords function
##'@examples
##'
##'data(ap)
##'plot(r<-npsurv(ap))              # survival function
##'plot(r$f, ap, fn="g")            # all gradients virtually zeros.
##'
##'data(cancer)
##'cancerRT = with(cancer, cancer[group=="RT",1:2])
##'plot(rt<-npsurv(cancerRT), xlim=c(0,60))                  # survival of RT 
##'cancerRCT = with(cancer, cancer[group=="RCT",1:2])
##'plot(rct<-npsurv(cancerRCT), add=TRUE, col="green3") # survival of RCT 
##'## as uniform dististrbutions.
##'plot(rt, add=TRUE, style="uniform", col="blue3")
##'plot(rct, add=TRUE, style="uniform", col="green3")
##'
##'## plot gradients; must supply data
##'plot(rt, cancerRT, fn="g")        # for group RT
##'plotgradidf(rct$f, cancerRCT)   # or, for group RCT
##'
##'@usage
##'\method{plot}{npsurv}(x, ...)
##'\method{plot}{idf}(x, data, fn=c("surv","grad"), ...)
##'plotsurvidf(f, style=c("box","uniform","left","right","midpoint"),
##'            xlab="Time", ylab="Survival Probability", col="blue3", fill=0,  
##'            add=FALSE, lty=1, lty.inf=2, xlim, ...)
##'plotgradidf(f, data, w=1, col1="red3", col2="blue3", 
##'            xlab="Survival Time", ylab="Gradient", xlim, ...)
##'
##'@export plot.npsurv
##'@export plot.idf
##'@export plotsurvidf
##'@export plotgradidf

plot.npsurv = function(x, ...) plot(x$f, ...)

plot.idf = function(x, data, fn=c("surv","grad"), ...) {
  fn = match.arg(fn)
  fnR = getFunction(paste("plot",fn,"idf",sep=""))
  switch(fn, "surv" = fnR(x, ...), "grad" = fnR(x, data, ...)  )
}

plotgradidf = function(f, data, w=1, col1="red3", col2="blue3", 
    xlab="Survival Time", ylab="Gradient", xlim, ...) {
  x2 = icendata(data, w)
  x = rbind(cbind(x2$t, x2$t), x2$o)
  w = c(x2$wt, x2$wo)
  dmat = Deltamatrix(x)
  D = dmat$Delta
  if(missing(xlim)) {
    upper = max(dmat$left, dmat$right[f$right<Inf])
    xlim = range(0, upper * 1.05)
  }
  m = length(dmat$left)
  p = double(m)
  p[dmat$left %in% f$left & dmat$right %in% f$right] = f$p
  # g = colSums(w * D / (D %*% p)[,1]) - sum(w)
  P = (D %*% p)[,1]
  g = crossprod(w/P, D)[1,] - sum(w)
  plot(dmat$left, g, type="h", col=col2, xlab=xlab, ylab=ylab, xlim=xlim, ...)
  lines(xlim, c(0,0), lty=1)
  j = p > 0
  ms = sum(j)
  points(dmat$left[!j], rep(0,m-ms), pch=1, col=col2, cex=1)
  points(dmat$right[!j], rep(0, m-ms), pch=20, col=col2, cex=.8)
  segments(dmat$left[!j], rep(0, m-ms),
           pmin(dmat$right[!j], xlim[2]), rep(0, m-ms),
           col=col2, lwd=3)
  points(dmat$left[j], rep(0,ms), pch=1, col=col1, cex=1)
  points(dmat$right[j], rep(0, ms), pch=20, col=col1, cex=.8)
  segments(dmat$left[j], rep(0, ms), pmin(dmat$right[j], xlim[2]), rep(0, ms),
           col=col1, lwd=3)
} 

plotsurvidf = function(f, style=c("box","uniform","left","right","midpoint"),
    xlab="Time", ylab="Survival Probability", col="blue3", fill=0,  
    add=FALSE, lty=1, lty.inf=2, xlim, ...) {
  style = match.arg(style)
  k = length(f$left)
  S = 1 - cumsum(f$p)
  upper = max(f$left, f$right[f$right != Inf])
  if(max(f$right) == Inf) point.inf = upper * 1.2
  else point.inf = upper
  if( missing(xlim) ) xlim = c(0, point.inf)
  m = length(f$p)
  if(!is.na(fill) && fill==0) {
    fill.hsv = drop(rgb2hsv(col2rgb(col))) * c(1, .3, 1)
    fill = hsv(fill.hsv[1], fill.hsv[2], fill.hsv[3], .3)
  }
  switch(style,
         box = {
           d = c(f$left[1], rep(f$right, rep(2,k)), f$right[k]) # right
           s = rep(c(1,S), rep(2,k+1))
           if(f$right[k] == Inf) d[2*k] = upper
           else d[2*k+2] = upper
           if( !add ) plot(d, s, type="n", col=col, xlim=xlim,
                           xlab=xlab, ylab=ylab, lty=lty, ...)
           if(style == "box") {
             Sc = c(1, S)
             j = which(f$right > f$left)
             rect(f$left[j], Sc[j+1], f$right[j], Sc[j], border=col,
                  col=fill)
           }
           lines(d, s, col=col, lty=lty, ...)
           lines(c(upper, point.inf), c(S[k-1],S[k-1]), col=col,
                 lty=lty.inf)
           if(f$right[k] != Inf) {       # left
             d = rep(c(f$left,f$right[k]), rep(2,k+1))
             s = c(1,rep(S, rep(2,k)),0)
           }
           else {
             d = rep(f$left, c(rep(2,k-1), 1))
             s = c(1,rep(S[-k], rep(2,k-1)))
           }
           add = TRUE
         }, 
         left = { d = rep(c(f$left,f$right[k]), rep(2,k+1))
                  s = c(1,rep(S, rep(2,k)),0)
                  d[2*k+2] = upper
                },
         right = { d = c(f$left[1], rep(f$right, rep(2,k)), f$right[k])
                   s = rep(c(1,S), rep(2,k+1))
                   if(f$right[k] == Inf) d[2*k] = upper
                   else d[2*k+2] = upper
                 },
         midpoint = { d1 = (f$left + f$right) / 2
                      d = c(f$left[1], rep(d1, rep(2,k)), f$right[k])
                      if(f$right[k] == Inf) d[2*k] = upper
                      else d[2*k+2] = upper
                      s = rep(c(1,S), rep(2,k+1))
                    },
         uniform = { d = c(rbind(f$left,f$right), rep(f$right[k],2))
                     if(f$right[k] == Inf) d[2*k] = upper
                     else d[2*k+2] = upper
                     s = c(1,rep(S, rep(2,k)),S[k])
                   }     )
  if( add ) lines(d, s, col=col, lty=lty,  ...)
  else plot(d, s, type="l", col=col, xlim=xlim, xlab=xlab, ylab=ylab,
            lty=lty, ...)
  abline(h=0, col="black")
  lines(c(0,f$left[1]), c(1,1), col=col)
  if(f$right[k] < Inf)
    lines(c(upper, point.inf), rep(0,2), col=col, lty=lty)
  else points(upper, S[k-1], col=col, pch=20)
}

## ==========================================================================
##
## Hierarchical CNM: a variant of the Constrained Newton Method for finding
## the NPMLE survival function of a data set containing interval censoring.
## This is a new method to build on those in the Icens and MLEcens
## packages.  It uses the idea of block subsets of the S matrix to move
## probability mass among blocks of candidate support intervals.
##
## Usage (parameters and return value) is similar to the methods in package
## Icens, although note the transposed clique matrix.
##
## Arguments:
##   data: Data
##   w:  Weights
##   D: Clique matrix, n*m (note, transposed c.f. Icens::EMICM,
##      MLEcens::reduc).  The clique matrix may contain conditional
##      probabilities rather than just membership flags, for use in HCNM
##      recursively calling itself.
##   p0: Vector (length m) of initial estimates for the probabilities of
##      the support intervals.
##   maxit: Maximum number of iterations to perform
##   tol: Tolerance for the stopping condition (in log-likelihood value)
##   blockpar:
##     NA or NULL  means choose a value based on the data (using n and r)
##     ==0  means same as cnm (don't do blocks)
##      <1  means nblocks is this power of sj, e.g. 0.5 for sqrt
##      >1  means exactly this block size (e.g. 40)
##   recurs.maxit: For internal use only: maximum number of iterations in
##      recursive calls
##   depth: For internal use only: depth of recursion
##   verb: For internal use only: depth of recursion
##
## Author: Yong Wang and Stephen S. Taylor
##
## Reference: Wang, Y. and Taylor, S. M. (2013). Efficient computation of
## nonparametric survival functions via a hierarchical mixture
## formulation. Statistics and Computing, 23, 713-725.
##
## ==========================================================================

hcnm = function(data, w=1, D=NULL, p0=NULL, maxit=100, tol=1e-6,
                blockpar=NULL, recurs.maxit=2, depth=1, verb=0) {
  if(missing(D)) {
    x2 = icendata(data, w)
    if(nrow(x2$o) == 0 || all(x2$o[,2] == Inf)) { # exact or right-censored only
      r0 = km(x2)
      r = list(f=r0$f, convergence=TRUE, ll=r0$ll, maxgrad=0, numiter=1)
      class(r) = "npsurv"
      return(r)
    }
    x = rbind(cbind(x2$t, x2$t), x2$o)
    nx = nrow(x)
    w = c(x2$wt, x2$wo)
    dmat = Deltamatrix(x)
    left = dmat$left
    right = dmat$right
    intervals = cbind(left, right)
    D = dmat$Delta
  }
  else {
    if (missing(p0)) stop("Must provide 'p0' with D.")
    if (!missing(data)) warning("D and data both provided.  LR ignored!")
    nx = nrow(D)
    w = rep(w, length=nx)
    intervals = NULL
  }
  n = sum(w)
  wr = sqrt(w)
  converge = FALSE
  m = ncol(D)
  m1 = 1:m
  nblocks = 1
  # maxdepth = depth
  i = rowSums(D) == 1
  r = mean(i)         # Proportion of exact observations
  if(is.null(p0)) {
    ## Derive an initial p vector.
    j = colSums(D[i,,drop=FALSE]) > 0
    while(any(c(FALSE,(i <- rowSums(D[,j,drop=FALSE])==0)))) {
      j[which.max(colSums(D[i,,drop=FALSE]))] = TRUE
    }
    p = colSums(w * D) * j
  }
  else { if(length(p <- p0) != m) stop("Argument 'p0' is the wrong length.") }
  p = p / sum(p)
  P = drop(D %*% p)
  ll = sum(w * log(P))
  evenstep = FALSE
  
  for(iter in 1:maxit) {
    p.old = p
    ll.old = ll
    S = D / P
    g = colSums(w * S)
    dmax = max(g) - n
    if(verb > 0) {
      cat("##### Iteration", i, "#####\n")
      cat("Log-likelihood: ", signif(ll, 6), "\n")
    }
    if(verb > 1) cat("Maximum gradient: ", signif(dmax, 6), "\n")
    if(verb > 2) {cat("Probability vector:\n"); print(p)} 
    j = p > 0
    if(depth==1) {
      s = unique(c(1,m1[j],m))
      if (length(s) > 1) for (l in 2:length(s)) {
        j[s[l-1] + which.max(g[s[l-1]:s[l]]) - 1] = TRUE
      }
    }
    sj = sum(j)
    ## BW: matrix of block weights: sj rows, nblocks columns
    if(is.null(blockpar) || is.na(blockpar))
      ## Default blockpar based on log(sj)
      iter.blockpar = ifelse(sj < 30, 0,
                             1 - log(max(20,10*log(sj/100)))/log(sj))
    else iter.blockpar = blockpar
    if(iter.blockpar==0 | sj < 30) {
      nblocks = 1
      BW = matrix(1, nrow=sj, ncol=1)
    }
    else {
      nblocks = max(1, if(iter.blockpar>1) round(sj/iter.blockpar)
                       else floor(min(sj/2, sj^iter.blockpar)))
      i = seq(0, nblocks, length=sj+1)[-1]
      if(evenstep) {
        nblocks = nblocks + 1
        BW = outer(round(i)+1, 1:nblocks, "==")
      }
      else BW = outer(ceiling(i), 1:nblocks, "==")
      storage.mode(BW) = "numeric"
    }

    for(block in 1:nblocks) {
      jj = logical(m)
      jj[j] = BW[,block] > 0
      sjj = sum(jj)
      if (sjj > 1 && (delta <- sum(p.old[jj])) > 0) {
        Sj = S[,jj]
        res = pnnls(wr * Sj, wr * drop(Sj %*% p.old[jj]) + wr, sum=delta)
        if (res$mode > 1) warning("Problem in pnnls(a,b)")
        p[jj] = p[jj] +  BW[jj[j],block] *
          (res$x * (delta / sum(res$x)) - p.old[jj])
      }
    }
    
    ## Maximise likelihood along the line between p and p.old
    p.gap = p - p.old              # vector from old to new estimate
    ## extrapolated rise in ll, based on gradient at old estimate
    ll.rise.gap = sum(g * p.gap) 
    alpha = 1
    p.alpha = p
    ll.rise.alpha = ll.rise.gap
    repeat {
      P = drop(D %*% p.alpha)
      ll = sum(w * log(P))
      if(ll >= ll.old && ll + ll.rise.alpha <= ll.old) {
        p = p.alpha               # flat land reached
        converge = TRUE
        break
      }
      if(ll > ll.old && ll >= ll.old + ll.rise.alpha * .33) {
        p = p.alpha               # Normal situation:  new ll is higher
        break
      }
      if((alpha <- alpha * 0.5) < 1e-10) {
        p = p.old
        P = drop(D %*% p)
        ll = ll.old
        converge = TRUE
        break
      }
      p.alpha = p.old + alpha * p.gap
      ll.rise.alpha = alpha * ll.rise.gap
    }
    if(converge) break

    if (nblocks > 1) {
      ## Now jiggle p around among the blocks
      Q = sweep(BW,1,p[j],"*")  # Matrix of weighted probabilities: [sj,nblocks]
      q = colSums(Q)            # its column sums (total in each block)
      ## Now Q is n*nblocks Matrix of probabilities for mixture components
      Q = sweep(D[,j] %*% Q, 2, q, "/")  
      if (any(q == 0)) {
        warning("A block has zero probability!")
      }
      else {
        ## Recursively call HCNM to allocate probability among the blocks 
        res = hcnm(w=w, D=Q, p0=q, blockpar=iter.blockpar,
                   maxit=recurs.maxit, recurs.maxit=recurs.maxit,
                   depth=depth+1)
        # maxdepth = max(maxdepth, res$maxdepth)
        if (res$ll > ll) {
          p[j] = p[j] * (BW %*% (res$pf / q))
          P = drop(D %*% p)
          ll = sum(w * log(P))  # should match res$lval
        }
      }
    }
    if(iter > 2) if( ll <= ll.old + tol ) {converge=TRUE; break}
    evenstep = !evenstep
  }
  list(pf=p, intervals=intervals, convergence=converge, method="hcnm", ll=ll,
       maxgrad=max(crossprod(w/P, D))-n, numiter=iter)
}
