#' coffee
#' @description While individual calibrated radiocarbon dates can span several centuries, combining multiple dates together with any chronological constraints can make a chronology much more robust and precise. This package uses Bayesian methods to enforce the chronological ordering of radiocarbon and other dates, for example for trees with multiple radiocarbon dates spaced at exactly known intervals (e.g., every 10 annual rings). Another example is sites where the relative chronological position of the dates is taken into account - the ages of dates further down a site must be older than those of dates further up.
#'
#' examples include:
#' 
#' rings()
#' 
#' strat()
#' @docType package
#' @author Maarten Blaauw <maarten.blaauw@qub.ac.uk> J. Andres Christen <jac@cimat.mx>
#' @importFrom grDevices dev.cur dev.off pdf dev.copy2pdf grey rgb dev.list extendrange
#' @importFrom graphics abline box curve hist image layout legend lines par plot points polygon segments rect axis mtext
#' @importFrom stats approx dbeta density dgamma dnorm dunif lm quantile rnorm weighted.mean coef median cor runif var as.ts
#' @importFrom utils read.csv read.table write.table packageName txtProgressBar setTxtProgressBar
#' @import IntCal
#' @name coffee
NULL

# ensure that functions such as pMC.age etc. are available immediately:
library(IntCal)

# chronological order forcing of fossils and environmental events, some in blocks, some after and before options

# also include (in future updates?): strat with gaps, sections w mixed dates, add uniform dist correction

# some ideas/blurbs for future updates:

# columns within the .csv files: lab ID, age, error, position (see below), calibration curve, model, outlier, ...
# position: 0, 1, 2, etc. for ordered dates. If some are within a Phase but unordered within the phase, they carry all the same number, e.g., 2, 2, 2. If they are ordered within a phase, perhaps do 2.0, 2.1, 2.2, etc.
# model: could be date, after, before, exponential, uniform, normal, (all requiring 1 or 2 parameters; date & error)

# no need for 'after' or 'before', because this can already be inferred from the position column? Add an entry with no information other than the order? For non-dated levels...

# strat has a modified prior to solve the problem mentioned by Steier and Rom 2000 (https://doi.org/10.1017/S0033822200058999), as well as Nicholls & Jones 2001. Is this OK? Check with Andres

# careful with referencing all relevant papers, incl. bwigg, OxCal, ...

# should optional columns for delta.R/STD and t.a/t.b be given, for tree and stat?

# add MCMC checks a la baconvergence for strat?

# also provide 95% ranges (or hpds?) hpds sometimes are >100%???

# an animation of strat would be fun, with 0 burnin and a short run. Draw coloured dots of the simulated ages and how they move through the iterations. Leave shadows of previous dots?

# References

#' Buck CE, Kenworthy JB, Litton CD, Smith AFM, 1991. Combining archaeological and radiocarbon information: a Bayesian approach to calibration. Antiquity 65, 808-821.

#' Christen JA, 2003. Bwigg: An Internet facility for Bayesian radiocarbon wiggle-matching. Internet Archaeology 13. \doi{10.11141/ia.13.2}

#' Christen JA, Fox C 2010. A General Purpose Sampling Algorithm for Continuous Distributions (the t-walk). Bayesian Analysis 5, 263-282.

### some generic functions

