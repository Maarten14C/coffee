


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
  if(name %in% c("block_example", "undated_example", "gaps_example"))
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



# From the dets file, find the relevant information on positions, blocks and gaps
structure <- function(dat) {
  # prepare for later calculations of ages, dates and gaps
  is.age <- dat[,5] %in% c(0:4,10)
  is.date <- dat[,5] %in% 0:4
  is.undated <- dat[,5] == 10
  xlength <- length(which(is.age)) # correct? Takes into account possible gaps
  pos.ages <- which(is.age)
  pos.dates <- which(is.date)
  pos.undated <- which(is.undated)
  dets <- dat[pos.ages,] # so no gaps

  # find any dates in blocks
  is.block <- which(duplicated(dat[,4])) # doesn't include the first instance of a duplicated entry...
  blocks <- c()
  if(length(is.block) > 0) {
    is.block <- dat[,4] %in% dat[is.block,4] # ... so here we find /all/ duplicate entries
    has.blocks <- TRUE
    blocks <- unique(dat[is.block,4])
  } else
    has.blocks <- FALSE
  pos.blocks <- c()
  if(length(is.block) > 0)
    pos.blocks <- which(is.block)

  # now deal with any gaps
  is.gap <- dat[,5] > 10
  pos.gaps <- which(is.gap)
  gaps <- c()
  m <- rep(2, xlength)
  if(length(pos.gaps) > 0)
    gaps <- dat[pos.gaps,]
  m[1] <- xlength # topmost date has M of xlength

  # fill vars r (row), p (position), d (dated), u (undated), b (block, where 0 = outside block, 1 = first block, etc.), g (gap)
  r <- 1:xlength
  p <- 1 # deals with position entries that are not integers
  k <- 0
  for(i in 2:xlength) {
	step <- dat[i,4] - dat[i-1,4]
    if(step > 0) # we've moved down a position
      p[i] <- p[i-1] +1
    if(step < 0)
      stop("positions in .csv file should be in chronological order")
	if(step == 0)
	  p[i] <- p[i-1]
  }

  # dates and undated levels which need an age
  d <- rep(0, xlength)
  d[which(dets[,5]>0)] <- 1
  u <- rep(0, xlength)
  u[which(dets[,5] == 10)] <- 1
  
  within.block <- c()
  # dates within blocks

  b <- rep(0, xlength)
  if(length(is.block) > 0) {
    above.block <- rep(NA, length(blocks))
    below.block <- above.block
    within.block <- list(c())
    for(i in 1:length(blocks)) {
      these <- which(p %in% blocks[i])
      b[these] <- i
      above.block[i] <- max(which(dets[,4] < blocks[i]))
      below.block[i] <- min(which(dets[,4] > blocks[i]))
      within.block[[i]] <- (above.block[i]+1) : (below.block[i]-1)
    }
  } else {
    above.block <- c()
    below.block <- c()
  }

  from <- (1:xlength)-1
  from[1] <- 1
  to <- from+1

  # for calculating l.x, any entries containing gaps should be removed
  if(length(gaps) > 0) {
    from <- from * (d+u)
    from <- from[from>0]
    to <- to * (d+u)
    to <- to[to>0]
  }
  to[1] <- xlength # topmost date has uniform time span covering the entire range

  # gaps
  g <- rep(0, xlength)
  has.gaps <- FALSE
  if(max(dat[,5]) > 10) { # then we have gaps
    has.gaps <- TRUE
    g[which(dat[,5] == 11)] <- 11 # exact
    g[which(dat[,5] == 12)] <- 12 # normal
    g[which(dat[,5] == 13)] <- 13 # gamma
    g[which(dat[,5] == 14)] <- 14 # uniform
  }

  gaps.younger <- c(); is.exact <- c(); is.normal <- c(); is.gamma <- c(); is.uniform <- c();
  has.exact <- c(); has.normal <- c(); has.gamma <- c(); has.uniform <- c();
  above.exact <- c(); above.normal <- c(); above.gamma <- c(); above.uniform <- c();
  exact.val1 <- c(); normal.val1 <- c(); normal.val2 <- c();
  gamma.val1 <- c(); gamma.val2 <-c();
  uniform.val1 <- c(); uniform.val2 <- c()
  if(length(gaps) > 0) {
    gaps.younger <- c()
    above.exact <- which(dat[,5] == 11) - 1
    above.normal <- which(dat[,5] == 12) - 1
    above.gamma <- which(dat[,5] == 13) - 1
    above.uniform <- which(dat[,5] == 14) - 1

    # now find the corresponding dets entries, based on the positions
    above.exact <- which(dets[,4] %in% dat[above.exact,4])
    above.normal <- which(dets[,4] %in% dat[above.normal,4])
    above.gamma <- which(dets[,4] %in% dat[above.gamma,4])
    above.uniform <- which(dets[,4] %in% dat[above.uniform,4])

    # positions in the gap variable
    is.exact <- which(gaps[,5] == 11)
    is.normal <- which(gaps[,5] == 12)
    is.gamma <- which(gaps[,5] == 13)
    is.uniform <- which(gaps[,5] == 14)

    exact.val1 <- gaps[is.exact,2]
    normal.val1 <- gaps[is.normal,2]
    normal.val2 <- gaps[is.normal,3]
    gamma.val1 <- gaps[is.gamma,2]
    gamma.val2 <- gaps[is.gamma,3]
    uniform.val1 <- gaps[is.uniform,2]
    uniform.val2 <- gaps[is.uniform,3]

    has.exact <- ifelse(length(above.exact) > 0, TRUE, FALSE)
    has.normal <- ifelse(length(above.normal) > 0, TRUE, FALSE)
    has.gamma <- ifelse(length(above.gamma) > 0, TRUE, FALSE)
    has.uniform <- ifelse(length(above.uniform) > 0, TRUE, FALSE)
    for(i in 1:nrow(gaps))
      gaps.younger[i] <- max(1, which(dat[,4] < gaps[i,4]))
  }

  struc <- list(dets=dets, r=r, p=p, d=d, u=u, b=b, g=g, m=m, xlength=xlength, from=from, to=to,
    is.age=is.age, is.date=is.date, is.undated=is.undated,
    pos.ages=pos.ages, pos.dates=pos.dates, pos.undated=pos.undated,
    has.blocks=has.blocks, blocks=blocks, is.block=is.block, pos.blocks=pos.blocks,
    above.block=above.block, below.block=below.block, within.block=within.block,
    has.gaps=has.gaps, gaps=gaps, is.gap=is.gap, pos.gaps=pos.gaps,
    gaps.younger=gaps.younger,
    has.exact=has.exact, is.exact=is.exact, above.exact=above.exact, exact.val1=exact.val1,
    has.normal=has.normal, is.normal=is.normal, above.normal=above.normal, normal.val1=normal.val1, normal.val2=normal.val2,
    has.gamma=has.gamma, is.gamma=is.gamma, above.gamma=above.gamma, gamma.val1=gamma.val1, gamma.val2=gamma.val2,
    has.uniform=has.uniform, is.uniform=is.uniform, above.uniform=above.uniform, uniform.val1=uniform.val1, uniform.val2=uniform.val2
    )

  return(struc)
}



# set initial values for the MCMC run
ballpark.strat <- function(method=1, dat, gaps, postbomb=postbomb, normal=normal, t.a=t.a, t.b=t.b, cc.dir=cc.dir) {
  is.exact <- c()
  ballpark.x <- rep(NA, nrow(dat))
  for(i in 1:nrow(dat))
    if(dat[i,5] %in% 0:4) { # then it's a date
      tmp <- caldist(dat[i,2], dat[i,3], cc = dat[i,5], postbomb = postbomb, normal = normal, t.a = t.a, t.b = t.b, cc.dir=cc.dir)
      ballpark.x[i] <- weighted.mean(tmp[,1], tmp[,2], na.rm=TRUE)
     # if(dat[i,3] == 0)
     #   is.exact <- c(is.exact, i)
     # if(is.null(dat[i,3]))
     #   is.exact <- c(is.exact, i)
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
  } 
  if(method == 2) 
    ballpark.x <- sort(ballpark.x)

  x0 <- sort(jitter(ballpark.x, factor=0.5)) # initial ball-park age estimates
  xp0 <- sort(jitter(ballpark.x, factor=0.5)) # twalk needs two sets of estimates

  # deal with exact ages (they have 0 or empty entries in error column)
  if(length(is.exact) > 0) {
    x0[is.exact] <- dat[is.exact,2]
    xp0[is.exact] <- dat[is.exact,2]
  }

  # here set any exact spans
  if(length(gaps) > 0)
    for(i in 1:nrow(gaps)) {
      younger <- max(1, which(dat[,4] < gaps[i,4]))
      x0[younger+1] <- x0[younger] + gaps[i,2]
      xp0[younger+1] <- xp0[younger] + gaps[i,2]
    }

  if(min(diff(x0)) < 0 || min(diff(xp0)) < 0)
    warning("reversals present in initial values")

  return(rbind(x0,xp0))
}



# scissors and thinner functions adapted from those of the rbacon R package

#' @name scissors
#' @title Remove the first n iterations.
#' @description Removes iterations of the MCMC time series, and then updates the output files.
#' @details The strat function will perform thousands millions of MCMC iterations, although usually only a fraction of these will be stored. 
#' The remaining MCMC iterations should be well mixed (the upper panel
#' of the fit of the iterations should show no undesirable features such as trends or sudden systematic drops or rises).
#' If the run has a visible remaining burn-in, scissors can be used to remove them.
#' To remove, e.g., the first 300 iterations, type \code{scissors(300)}. To remove the last 300 iterations, type \code{scissors(-300)}. To remove iterations 300 to 600, type \code{scissors(300:600)}.
#'
#' @param burnin Number of iterations to remove of the iterative time series. If this value is higher than the amount of remaining iterations,
#' a warning is given and the iterations are not removed. If the provided number is negative, the iterations will be removed from the end of the run, not from the start. If a range is given, this range of iterations is removed.
#' @param set Detailed information of the current run, stored within this session's memory as variable \code{info}.
#' @param write Whether or not to write the changes to the output file. Defaults to TRUE.
#' @param save.info Whether or not to store a variable `info' in the session which contains the run input, output and settings. Defaults to \code{save.info=TRUE}.
#' @return NA
#'
#' @export
scissors <- function(burnin, set=get('info'), write=TRUE, save.info=TRUE) {
  output <- fastread(paste0(set$strat.dir, ".out"))
  energy <- fastread(paste0(set$strat.dir, "_energy.out"))[,1]

  if(length(burnin) > 1) {
    if(length(burnin) >= nrow(output))
      stop("cannot remove that many iterations, there would be none left!", call.=FALSE)
    output <- output[-burnin,]
    energy <- energy[-burnin]
  } else {
      if(abs(burnin) >= nrow(output))
        stop("cannot remove that many iterations, there would be none left!", call.=FALSE)
      if(burnin > 0) {
        output <- output[-(1:burnin),]
		energy <- energy[-(1:burnin)]
      } else {
          output <- output[-((nrow(output)-abs(burnin)):nrow(output)),]
          energy <- energy[-((length(energy)-abs(burnin)):length(energy))]
        }
    }

  if(write) {
    fastwrite(output, paste0(set$strat.dir, ".out"), col.names=FALSE, row.names=FALSE)
	fastwrite(energy, paste0(set$strat.dir, "_energy.out"), col.names=FALSE, row.names=FALSE)
  }
  set$output <- output
  set$Tr <- nrow(output)
  set$Us <- energy 
  set$Ups <- energy
  if(save.info)
    assign_to_global("info", set)
  invisible(set)
}



#' @name thinner
#' @title Thin iterations.
#' @description Randomly thin iterations by a given proportion, for example if autocorrelation is visible within the MCMC series.
#' @details From all iterations, a proportion is removed with to-be-removed iterations sampled randomly among all iterations.
#' @param proportion Proportion of iterations to remove. Should be between 0 and 1. Default \code{proportion=0.1}.
#' @param set Detailed information of the current run, stored within this session's memory as variable \code{info}.
#' @param write Whether or not to write the changes to the output file. Defaults to TRUE.
#' @param save.info Whether or not to store a variable `info' in the session which contains the run input, output and settings. Defaults to \code{save.info=TRUE}.
#' @return NA
#'
#' @export
thinner <- function(proportion=0.1, set=get('info'), write=TRUE, save.info=TRUE) {
  output <- fastread(paste0(set$strat.dir, ".out"))
  energy <- fastread(paste0(set$strat.dir, "_energy.out"))[,1]
	
  if(proportion >= 1)
    stop("cannot remove that many iterations, there would be none left!", call.=FALSE)
  proportion <- sample(nrow(output), proportion*nrow(output))
  output <- output[-proportion,]
  energy <- energy[-proportion]

  if(write) {
    fastwrite(output, paste0(set$strat.dir, ".out"), col.names=FALSE, row.names=FALSE)
	fastwrite(energy, paste0(set$strat.dir, "_energy.out"), col.names=FALSE, row.names=FALSE)
  }
  set$output <- output
  set$Tr <- nrow(output)
  set$Us <- energy 
  set$Ups <- energy
  if(save.info)
    assign_to_global("info", set)
  invisible(set)
}
