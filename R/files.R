
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



