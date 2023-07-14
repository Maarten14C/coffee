


# internal functions to speed up reading and writing files, using the data.table R package if present
fastread <- function(fnam, ...)
  if("data.table" %in% (.packages()))
    as.data.frame(fread(fnam, ...)) else
      read.table(fnam, ...)



fastwrite <- function(out, fnam, ...)
  if("data.table" %in% (.packages()))
    fwrite(as.data.frame(out), fnam, ...) else
      write.table(out, fnam, ...)



# read the dets file of a strat site
read.strat <- function(name="mystrat", strat.dir="strats", sep=",", normal=TRUE, delta.R=0, delta.STD=0, t.a=3, t.b=4, cc=1) {
  stratdir <- assign_dir(strat.dir, name, "strat.dir", ask=FALSE)
  if(name %in% c("block_example", "undated_example"))
	file.copy(system.file(file.path("extdata", 
      paste0(name, ".csv")), package="coffee"),
        stratdir, recursive=TRUE, overwrite=FALSE)  
  dat <- read.table(file.path(stratdir, paste0(name, ".csv")), header=TRUE, sep=sep)

  # some initial checks of the data file
  if(ncol(dat) < 4)
    stop("Need at least 4 columns in strat file. Please check", call.=FALSE)
  if(ncol(dat) == 4)
    dat <- cbind(dat, cc) # then add a column noting the calibration curve
  # we need increasing stratigraphical positions - larger numbers are further down
  if(min(diff(dat[,4])) < 0)
    stop("Positions (column 4) should be in chronological order, increasing downward", call.=FALSE)
  if(min(dat[,5]) < 0)
    stop("Unexpected calibration curve entry (column 5) in strat file. Please check", call.=FALSE)

  # age offsets and t.a/t.b can also be provided
  if(ncol(dat) == 5) { # no columns for age offsets or t.a/t.b in .csv file...
    dat[,2] <- dat[,2] - delta.R # but they can be provided as options
    dat[,3] <- sqrt(dat[,3]^2 + delta.STD^2)
  } else
    if(ncol(dat) > 5)
      if(tolower(colnames(dat))[6] %in% c("delta.r", "delta", "offset", "reservoir")) {
        dat[,2] <- dat[,2] - dat[,6] # delta.mean
        dat[,3] <- sqrt(dat[,3]^2 + dat[,7]^2) # delta.STD
      }
    if(ncol(dat) == 9) { # then we expect t.a and t.b in columns 8 resp. 9
      t.a <- dat[,8]
      t.b <- dat[,9]
    }
  return(list(dets=dat, normal=normal, t.a=t.a, t.b=t.b, delta.R=delta.R, delta.STD=delta.STD))
}




# set initial values for the MCMC run
ballpark.strat <- function(method=1, dat, gaps, postbomb=postbomb, normal=normal, t.a=t.a, t.b=t.b, cc.dir=cc.dir) {
  ballpark.x <- rep(NA, nrow(dat))
  for(i in 1:nrow(dat))
    if(dat[i,5] %in% 0:4) { # then it's a date
      tmp <- caldist(dat[i,2], dat[i,3], cc = dat[i,5], postbomb = postbomb, normal = normal, t.a = t.a, t.b = t.b, cc.dir=cc.dir)
      ballpark.x[i] <- weighted.mean(tmp[,1], tmp[,2], na.rm=TRUE)
    } 
  for(i in which(dat[,5] == 10)) {
    surround <- c(ballpark.x[i+1],ballpark.x[i-1])
    ballpark.x[i] <- runif(1, min(surround), max(surround))	
  }

  # place the ages in order
  dated.pos <- which(dat[,5] %in% c(0:4,10))
  if(method == 1) {
    ballpark.lm <- lm(ballpark.x ~ dat[dated.pos,4], weights=dat[dated.pos,3]^2, na.action=na.omit)
    ballpark.x <- ballpark.lm$fitted.values
  } else
      if(method == 2) 
        ballpark.x <- sort(ballpark.x)

  x0 <- sort(jitter(ballpark.x, factor=0.5)) # initial ball-park age estimates
  xp0 <- sort(jitter(ballpark.x, factor=0.5))  # twalk needs two sets of estimates

  # here set any exact spans
  if(length(gaps) > 0)
    for(i in 1:nrow(gaps)) {
      younger <- max(1, which(dat[,4] < gaps[i,4]))
      x0[younger+1] <- x0[younger] + gaps[i,2]
      xp0[younger+1] <- xp0[younger] + gaps[i,2]
    }

  return(rbind(x0,xp0))
}

