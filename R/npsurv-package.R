

##'Air Conditioner Failure Data
##'
##'Contains the number of operating hours between successive failure times of
##'the air conditioning systems in Boeing airplanes
##'
##'
##'@name acfail
##'@docType data
##'@format A numeric vector storing the failure times.
##'@seealso \code{\link{Uhaz}}.
##'@references Proschan, F. (1963). Theoretical explanation of observed
##'decreasing failure rate. \emph{Technometrics}, \bold{5}, 375-383.
##'@source Proschan (1963)
##'@keywords datasets
##'@examples
##'
##'data(acfail)
##'r = Uhaz(acfail, deg=2)
##'plot(r$h, fn="h")
##'plot(r$h, fn="d")
##'
NULL





##'Angina Pectoris Survival Data
##'
##'Contains the survival times in years from the time of diagnosis for 2418
##'male patients with angina pectoris. Some patients are lost to follow-up,
##'hence giving right-censored observations. Each integer-valued survival time
##'is treated as being censored within a one-year interval.
##'
##'
##'@name ap
##'@docType data
##'@format
##'
##'A data frame with 30 observations and 3 variables:
##'
##'\code{L}: left-end point of an interval-censored retraction time;
##'
##'\code{R}: right-end point of an interval-censored retraction time;
##'
##'\code{count}: number of patients in the interval.
##'@seealso \code{\link{npsurv}}.
##'@references Lee, E. T. and Wang, J. W. (2003). \emph{Statistical Methods for
##'Survival Data Analysis}. Wiley.
##'@source Lee and Wang (2003), page 92.
##'@keywords datasets
##'@examples
##'
##'data(ap)
##'r = Uhaz(ap, deg=2)           # smooth U-shaped hazard
##'plot(r$h, fn="h")             # hazard
##'plot(r$h, fn="d")             # density
##'
##'# NPMLE and shape-restricted estimation
##'plot(npsurv(ap), fn="s")      # survival under no shape restriction
##'plot(r$h, fn="s", add=TRUE)   # survival with smooth U-shaped hazard
##'
NULL





##'Breast Retraction Times after Beast Cancer Treatments.
##'
##'Contains the breast retraction times in months for 94 breast cancer patients
##'who received either radiation therapy or radiation therapy plus adjuvant
##'chemotherapy.
##'
##'
##'@name cancer
##'@docType data
##'@format A data frame with 94 observations and 3 variables:
##'
##'L: left-end points of the interval-censored retraction times;
##'
##'R: right-end points of the interval-censored retraction times;
##'
##'group: either \code{RT} (radiation therapy) or \code{RCT} (radiation therapy
##'plus adjuvant chemotherapy).
##'@seealso \code{\link{npsurv}}.
##'@references Finkelstein, D. M. and R. A. Wolfe (1985). A semiparametric
##'model for regression analysis of interval-censored failure time data.
##'\emph{Biometrics}, \bold{41}, pp.933-945.
##'@source Finkelstein and Wolfe (1985).
##'@keywords datasets
##'@examples
##'
##'data(cancer)
##'i = cancer$group == "RT"
##'plot(npsurv(cancer[i,1:2]), xlim=c(0,60))
##'plot(npsurv(cancer[!i,1:2]), add=TRUE, col="green3")
##'
NULL





##'Gastric Cancer Survival Data
##'
##'Contains the survival times of 45 gastrointestinal tumor patients who were
##'treated with both chemotherapy and radiotherapy. It has both exact and
##'right-censored observations.
##'
##'
##'@name gastric
##'@docType data
##'@format A data frame with 30 observations and 3 variables:
##'
##'L: left-end points of the interval-censored survival times;
##'
##'R: right-end points of the interval-censored survival times.
##'@seealso \code{\link{npsurv}}, \code{\link{Uhaz}}.
##'@references Klein, J. P. and Moeschberger, M. L. (2003).  \emph{Survival
##'Analysis: Techniques for Censored and Truncated Data (2nd ed.)}.  Springer.
##'@source Klein and Moeschberger (2003), page 224.
##'@keywords datasets
##'@examples
##'
##'data(gastric)
##'plot(npsurv(gastric), col="grey")      # survival function
##'plot(h0<-Uhaz(gastric, deg=0)$h, fn="s", add=TRUE, col="green3")
##'plot(h1<-Uhaz(gastric, deg=1)$h, fn="s", add=TRUE)
##'plot(h2<-Uhaz(gastric, deg=2)$h, fn="s", add=TRUE, col="red3")
##'
##'plot(h0, fn="h", col="green3")         # hazard function
##'plot(h1, fn="h", add=TRUE)
##'plot(h2, fn="h", add=TRUE, col="red3")
##'
##'plot(h0, fn="d", col="green3")         # density function
##'plot(h1, fn="d", add=TRUE)
##'plot(h2, fn="d", add=TRUE, col="red3") 
##'
##'
NULL





##'Remission Times for Acute Leukemia Patients
##'
##'Contains remission times in weeks of 42 acute leukemia patients, who
##'received either the treatment of drug 6-mercaptopurine or the placebo
##'treatment. Each remission time is either exactly observed or right-censored.
##'
##'
##'@name leukemia
##'@docType data
##'@format A data frame with 42 observations and 3 variables:
##'
##'L: left-end points of the interval-censored remission times in weeks;
##'
##'R: right-end points of the interval-censored remission times;
##'
##'group: either 6-MP (6-mercaptopurine) or Placebo.
##'@seealso \code{\link{npsurv}}.
##'@references Freireich, E. O. et al. (1963). The effect of 6-mercaptopmine on
##'the duration of steroid induced remission in acute leukemia. \emph{Blood},
##'\bold{21}, 699-716.
##'@source Freireich et al. (1963).
##'@keywords datasets
##'@examples
##'
##'data(leukemia)
##'i = leukemia[,"group"] == "Placebo"
##'plot(npsurv(leukemia[i,1:2]), xlim=c(0,40), col="green3") # placebo
##'plot(npsurv(leukemia[!i,1:2]), add=TRUE)                  # 6-MP
##'
##'## Treat each remission time as interval-censored:
##'x = leukemia
##'ii = x[,1] == x[,2]
##'x[ii,2] = x[ii,1] + 1
##'plot(npsurv(x[i,1:2]), xlim=c(0,40), col="green3")        # placebo
##'plot(npsurv(x[!i,1:2]), add=TRUE)                         # 6-MP
##'
NULL





##'Angina Pectoris Survival Data
##'
##'Contains the answers of 191 California high school students to the question:
##'"When did you first use marijuana?". An answer can be an exact age, or "I
##'have never used it", which gives rise to a right-censored observation, or "I
##'have used it but cannot recall just when the first time was", which gives
##'rise to a left-censored observation.
##'
##'
##'@name marijuana
##'@docType data
##'@format A data frame with 21 observations and 3 variables:
##'
##'L: left-end point of an interval-censored time;
##'
##'R: right-end point of an interval-censored time;
##'
##'count: number of students in the interval.
##'@seealso \code{\link{npsurv}}.
##'@references Turnbull and Weiss (1978). A likelihood ratio statistic
##'fortesting goodness of fit with randomly censored data. \emph{Biometrics},
##'\bold{34}, 367-375.
##'
##'Klein and Moeschberger (2003). \emph{Survival Analysis: Techniques for
##'Censored and Truncated Data} (2nd ed.). Springer
##'@source Turnbull and Weiss (1978). See also Klein and Moeschberger (1997),
##'page 17.
##'@keywords datasets
##'@examples
##'
##'data(marijuana)
##'r = Uhaz(marijuana, deg=2)
##'plot(r$h, fn="h")
##'plot(r$h, fn="s")
##'
NULL





##'New Zealand Mortality in 2000
##'
##'Contains the number of deaths of Maori and Non-Maori people at each age in
##'New Zealand in 2000.
##'
##'Data contains no age with zero death.
##'
##'@name nzmort
##'@docType data
##'@format A data frame with 210 observations and 3 variables:
##'
##'age: at which age the deaths occurred;
##'
##'deaths: number of people died at the age;
##'
##'ethnic: either Maori or Non-Maori.
##'@seealso \code{\link{Uhaz}}.
##'@source \url{https://www.mortality.org/}
##'@keywords datasets
##'@examples
##'
##'data(nzmort)
##'x = with(nzmort, nzmort[ethnic=="maori",])[,1:2]      # Maori mortality
##'# x = with(nzmort, nzmort[ethnic!="maori",])[,1:2]    # Non-Maori mortality
##'
##'## As exact observations
##'# Plot hazard functions
##'h0 = Uhaz(x[,1]+0.5, x[,2], deg=0)$h    # U-shaped hazard
##'plot(h0, fn="h", col="green3", pch=2)
##'h1 = Uhaz(x[,1]+0.5, x[,2], deg=1)$h    # convex hazard
##'plot(h1, fn="h", add=TRUE, pch=1)
##'h2 = Uhaz(x[,1]+0.5, x[,2], deg=2)$h    # smooth U-shaped hazard
##'plot(h2, fn="h", add=TRUE, col="red3")
##'
##'# Plot densities
##'age = 0:max(x[,1])
##'count = integer(length(age))
##'count[x[,"age"]+1] = x[,"deaths"]
##'barplot(count/sum(count), space=0, col="lightgrey", ylab="Density")
##'axis(1, pos=NA, at=0:10*10)
##'plot(h0, fn="d", add=TRUE, col="green3", pch=2)
##'plot(h1, fn="d", add=TRUE, col="blue3", pch=1)
##'plot(h2, fn="d", add=TRUE, col="red3", pch=19)
##'
##'## As interval-censored observations
##'# Plot hazard functions
##'x2 = cbind(x[,1], x[,1]+1, x[,2])
##'h0 = Uhaz(x2, deg=0)$h      # U-shaped hazard
##'plot(h0, fn="h", col="green3", pch=2)
##'h1 = Uhaz(x2, deg=1)$h      # convex hazard
##'plot(h1, fn="h", add=TRUE, pch=1)
##'h2 = Uhaz(x2, deg=2)$h      # smooth U-shaped hazard
##'plot(h2, fn="h", add=TRUE, col="red3", pch=1)
##'
##'# Plot densities
##'barplot(count/sum(count), space=0, col="lightgrey")
##'axis(1, pos=NA, at=0:10*10)
##'plot(h0, fn="d", add=TRUE, col="green3", pch=2)
##'plot(h1, fn="d", add=TRUE, col="blue3", pch=1)
##'plot(h2, fn="d", add=TRUE, col="red3", pch=19)
##'
NULL



