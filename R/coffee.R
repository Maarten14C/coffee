# todo:
# check strat, are blocks dealt with correctly?
# plot: gaps are not plotted on the right depth if there's a block above them

# dates close to 0 14C BP (e.g. 260 +- 20) go wrong - too close to end cc? Add option youngest.age???

# colouring of dates in blocks doesn't go well

# add optional circle symbol below each age-modelled mode, coloured according to the info$overlap value

# add a strat option for multiple cores which share events?

# add MCMC checks a la baconvergence for strat? It is using IAT already to propose thinning values

# an animation of strat would be fun, with 0 burnin and a short run. Draw coloured dots of the simulated ages and how they move through the iterations. Leave shadows of previous dots?

# preparing the rings function to facilitate animations would also be helpful

# ensure that functions such as pMC.age etc. are available immediately:
library(rintcal)

### some internal functions to write and read things

# function to load results in global environment - copied from rbacon package
# parameter position defaults to 1, which equals an assignment to the global environment
assign_to_global <- function(key, val, pos=1) {
  assign(key, val, envir=as.environment(pos) )
}


assign_dir <- function(umbrella, name, option.name, ask=FALSE, talk=TRUE) {
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

