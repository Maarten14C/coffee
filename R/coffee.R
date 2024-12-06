# todo:

# make some cases to test before submitting any update to CRAN (e.g. with cc=0, postbomb, blocks, undated levels, ...)

# make it possible to write in the cc column gap_gamma instead of 13, etc. g, n, e, u

# plot: gaps are not plotted on the correct depth if there's a block above them

# error if the bottom has a block. Can add a gap and an undated level below it and then take it away.... not nice. Adding just an undated level causes error in if (x[below.block[i]] < max(x[thisblock]))  

# add a strat option for multiple cores which share events?

# add MCMC checks a la baconvergence for strat? It is using IAT already to propose thinning values

# an animation of strat would be fun, with 0 burnin and a short run. Draw coloured dots of the simulated ages and how they move through the iterations. Leave shadows of previous dots?

# preparing the rings function to facilitate animations would also be helpful

# ensure that functions such as pMCtoC14 etc. are available immediately:
# not loading rice now because it's only needed for hpd, which has a bug for distributions from short iterations
# library(rice)

### some internal functions to write and read things

# function to load results in global environment - copied from rbacon package
# parameter position defaults to 1, which equals an assignment to the global environment
assign_to_global <- function(key, val, pos=1) 
  assign(key, val, envir=as.environment(pos) )



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

