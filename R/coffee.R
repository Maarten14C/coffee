#' coffee
#' @description While individual calibrated radiocarbon dates can span several centuries, combining multiple dates together with any chronological constraints can make a chronology much more robust and precise. This package uses Bayesian methods to enforce the chronological ordering of radiocarbon and other dates, for example for trees with multiple radiocarbon dates spaced at exactly known intervals (e.g., every 10 annual rings; Christen 2003). Another example is sites where the relative chronological position of the dates is taken into account - the ages of dates further down a site must be older than those of dates further up (e.g., Buck et al. 1991). MCMC runs are done using the t-walk (Christen and Fox 2010).
#'
#' examples include:
#' 
#' rings()
#' 
#' strat()
#'
#' Note that several R packages exist to run and/or extract results from OxCal (see OxcAAR and https://github.com/gavinsimpson/roxcal).
#' However, coffee is meant as a stand-alone package that doesn't rely on having OxCal installed.
#' Additionally, the data files of coffee are meant to be more easily readable and writable by humans.
#'
#' References
#' Buck CE, Kenworthy JB, Litton CD, Smith AFM, 1991. Combining archaeological and radiocarbon information: a Bayesian approach to calibration. Antiquity 65, 808-821.
#'
#' Christen JA, 2003. Bwigg: An Internet facility for Bayesian radiocarbon wiggle-matching. Internet Archaeology 13. \doi{10.11141/ia.13.2}
#'
#' Christen JA, Fox C 2010. A General Purpose Sampling Algorithm for Continuous Distributions (the t-walk). Bayesian Analysis 5, 263-282
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

# some ideas/blurbs for future updates:

# strats with uncertain sizes of gaps, sections w multiple dates in a position (that should date to the exact same year), starts and ends

# columns within the .csv files: lab ID, age, error, position (see below), calibration curve, model, outlier, ...
# position: 0, 1, 2, etc. for ordered dates. If some are within a Phase but unordered within the phase, they carry all the same number, e.g., 2, 2, 2. If they are ordered within a phase, perhaps do 2.0, 2.1, 2.2, etc.
# model: could be date, after, before, exponential, uniform, normal, (all requiring 1 or 2 parameters; date & error)

# no need for 'after' or 'before', because this can already be inferred from the position column? Add an entry with no information other than the order? For non-dated levels...

# add MCMC checks a la baconvergence for strat?

# also provide 95% ranges (or hpds?) hpds sometimes are >100%?

# an animation of strat would be fun, with 0 burnin and a short run. Draw coloured dots of the simulated ages and how they move through the iterations. Leave shadows of previous dots?



### some internal functions to write and read things

# function to load results in global environment - copied from rbacon package
# parameter position defaults to 1, which equals an assignment to the global environment
assign_to_global <- function(key, val, pos=1) {
  assign(key, val, envir=as.environment(pos) )
}


assign_dir <- function(umbrella, name, option.name, ask=TRUE, talk=TRUE) {
  # first check if the umbrella folder exists already. If it doesn't, ask for permission to create it
  if(!dir.exists(umbrella)) {
    if(ask) {
      wdir <- FALSE
      ans <- readline(paste0("OK to create the folder ", umbrella, " and place files of this and future runs there? (Y/n) "))
      ans <- tolower(substr(ans,1,1))[1]
      if(ans=="y" || ans=="")
        wdir <- dir.create(umbrella, FALSE) else
          stop("No problem. Please provide an alternative folder location for ", option.name, call.=FALSE)
      if(!wdir)
        stop("cannot write into the current directory.\nPlease set ", option.name, " to somewhere where you have writing access, e.g. Desktop or ~.", call.=FALSE)
    } else
        dir.create(umbrella)
  }

  # now write the folder for the run.
  # We've already been given permission from the user to make folders and write files within the folder provided.
  rundir <- file.path(umbrella, name)
  if(talk)
    message("The run's files will be put in this folder: ", rundir)
  if(!dir.exists(rundir))
    wdir <- dir.create(rundir)
  return(rundir)
}

