
#' @name strat
#' @title Model chronologically ordered dates
#' @description Model radiocarbon dates (or dates that are already on the cal BP scale) of a deposit that is known to have accumulated over time, and for which therefore the dated positions can be safely assumed to are in chronological order.
#' @details Dates further down the sequence should have older ages than dates further up, even though owing to scatter, the dates themselves might not be in exact chronological order. The amount of scatter, the laboratory error and an offset can also be modelled.
#' The function reads in a .csv file of a specific format. The first column contains the names of the dates/information, the second column has the age(s) (uncalibrated for radiocarbon dates, as they will be calibrated during the modelling), the third column their errors, the fourth column their position (see below), and the fifth column cc, the calibration curve information. Additional columns for the reservoir effect (delta.R and delta.STD) and the student-t model (t.a and t.b) can be added, much like rbacon .csv files.
#' The positions of the dates (column 4) should be entered with the topmost, youngest levels first, and then working downward toward the oldest levels. The topmost position gets the lowest number (e.g., 0), and each subsequent entry should have a higher position number to ensure that the levels are ordered in time. Dates in 'blocks' where there is no known age ordering between the dates in the block (but where that block is known to be older than the level above it and younger than the level below it) should all get the same position in column 4.  
#' The function does not only deal with dates (radiocarbon or otherwise), but can also model undated levels and a range of gaps between the dated levels. This is done mostly through column 5 in the .csv files, where a 0 is for dates on the cal BP scale, 1 for radiocarbon dates that require calibration with IntCal20, 2 with Marine20, 3 with SHCal20, 4 a custom calibration curve; additional information can be provided by adding entries for undated levels (cc=10), gaps of exactly known length (cc=11), normally distributed gap lengths (cc=12), or gamma distributed gap lengths (cc=13).
#' The age estimates are obtained through a t-walk MCMC run (Christen and Fox 2010). In this process, initial ball-park point estimates for the ages of each dated depth are given, checked for their chronological ordering (and for the sizes of any gaps) and then modified through many iterations. For each iteration, a random dated depth is chosen and its age changed by just a little nudge, a check is performed to ensure that all age estimates remain in chronological order (and that gap sizes remain obeyed), and the 'energy' or likelihood of the age estimates is calculated (iterations where all ages fit well within the calibrated distributions receive a higher energy; see \code{l.calib}). 
#' Then this iteration with the updated group of age estimates is either accepted or rejected. The acceptance probability depends on the iteration's energy; if its energy is higher than that of the previous iteration it is accepted, but if it is lower, it is accepted with a probability proportional to its relative energy. Therefore, over many iterations the process will 'learn' from the data and find high-energy combinations of parameter values that fit with the prior constraints that the ages should be ordered chronologically.
#' Because the iterations are based on a process of modifying values of one parameter each iteration, and because some iterations will not be accepted, the MCMC output will often have a large degree of dependence between neighbouring iterations. Therefore, some thinning will have to be done, by storing only one every few iterations (default 20). Also, since the initial ball-park estimates could be quite wrong, the first 100 or so iterations should also be discarded (burnin). 
#' It is thus important to check the time-series of the energy after the run. We don't want to see a remaining burn-in at the start, and we don't want to see a noticeable 'structure' where iterations remain in approximately or entirely the same spot for a long time. Instead, an ideal run will look like white noise.
#' @param name Name of the stratigraphy dataset. Defaults to \code{"mystrat"}.
#' @param strat.dir The directory where the folders of the individual stratigraphies live. Defaults to \code{treedir="strats"}.
#' @param its Amount of iterations to be run. Setting this to low numbers (e.g., 1000) will result in fast but less stable and less reliable runs. Higher values will take longer but result in more stable and robust runs. Defaults to \code{50000}. Aim to set this to such values that at least 3000 iterations remain after removing the burnin and thinning.
#' @param burnin Amount of iterations to remove at the start of the run. Defaults to \code{100}.
#' @param thinning After running all iterations, only some will be stored. For example, if thinning is set at the default \code{50}, only every 50th MCMC iteration will be stored, and the others will be discarded. This is to remove the dependence between neighbouring MCMC iterations. Defaults to a value calculated from the MCMC run itself.
#' @param internal.thinning Does internal thinning during the MCMC process, storing only every 'internal.thinning' MCMC iterations.
#' @param min.its Minimum number of (remaining) iterations, below which a warning is given. Defaults for now to 3,000 iterations.
#' @param write.MCMC Especially longer runs or sites with many dates can take up lots of memory. For such cases, MCMC iterations are stored in temporary files rather than in memory. Defaults to TRUE. 
#' @param init.ages By default, the ballpark age estimates to feed the MCMC are calculated automatically, however they can also be provided manually, as two rows of values (for x0 and xp0) which have to obey the assumptions (e.g., chronological ordering).
#' @param ballpark.method Initial, ballpark values for the initial ages (if not provided by the option `init.ages`). Can be 1 (based on a linear model) or 2 (based on sorted random draws).
#' @param y.scale The scale of the vertical axis of the main plot. This can be the positions of the dated levels (`positions`) or their position order (`dates`). 
#' @param showrun Whether or not to show how the MCMC process is progressing during the run. Defaults to \code{FALSE}.
#' @param sep Separator for the fields in the .csv file. Defaults to a comma.
#' @param normal Calculations can be done assuming that the measurements are normally distributed. By default this is set to FALSE and a student-t distribution is used (Christen and Perez 2009)
#' @param delta.R The ages can be modelled to have an offset. The mean is 0 by default.
#' @param delta.STD The error of the offset. Set to 0 by default.
#' @param t.a First parameter for the student-t distribution (defaults to 3; higher numbers make the distribution approximate the normal distribution more).
#' @param t.b Second parameter for the student-t distribution (defaults to 4; higher numbers make the distribution 
#' @param cc Calibration curve to be used. Could be 1 (IntCal20; default), 2 (Marine20), 3 (SHCal20) or 4 (custom curve).
#' @param cc.dir Directory of calibration curve. Keep empty for the default value.
#' @param postbomb Negative C-14 ages should be calibrated using a postbomb curve. This could be 1 (northern-hemisphere region 1), 2 (NH region 2), 3 (NH region 3), 4 (southern hemisphere regions 1-2), or 5 (SH region 3).
#' @param BCAD The calendar scale of graphs and age output-files is in \code{cal BP} by default, but can be changed to BC/AD using \code{BCAD=TRUE}. Needs more work probably.
#' @param ask Whether or not to ask if a folder should be made (if required).
#' @param talk Whether or not to provide feedback on folders written into and on what is happening.
#' @param talk Whether or not show the progress of the MCMC run.
#' @param clean.garbage Whether or not to clean up the memory 'garbage collection' after a run. Recommendable if you have many dates or long runs.
#' @param age.span Expected age span. Defaults to run from the current year in AD to 55e3 which is the current cal BP limit for C-14 dates. If older, non-14C dates are present, age.span is set to the larger of the radiocarbon limit or twice the age of the oldest non-radiocarbon age.
#' @param ... Options for the plot. See \code{plot.strat}.
#' @return a variable 'info' which contains the dating and modelling information to produce a plot. Also calls the function \code{draw.strat} to produce a plot of the results.
#' @examples
#' \dontrun{
#' tmp <- tempdir()
#' strat(, strat.dir=tmp, its=1000, thinning=1, internal.thinning=1)
#' }
#' @references
#' Bronk Ramsey C, 1995. Radiocarbon calibration and analysis of stratigraphy: The OxCal program. Radiocarbon 37, 425 – 430.
#'
#' Buck CE, Kenworthy JB, Litton CD, Smith AFM, 1991. Combining archaeological and radiocarbon information: a Bayesian approach to calibration. Antiquity 65, 808-821.
#'
#' Buck et al. 1999. BCal: an on-line Bayesian radiocarbon calibration tool. Internet Archaeology 7. 
#'
#' Christen JA, Fox C 2010. A general purpose sampling algorithm for continuous distributions (the t-walk). Bayesian Analysis 5, 263-282. 
#'
#' Nicholls G, Jones M 2001. Radiocarbon dating with temporal order constraints. Journal of the Royal Statistical Society: Series C (Applied Statistics) 50, 503-521.
#' @export
strat <- function(name="mystrat", strat.dir="strats", its=5e4, burnin=100, thinning=c(), internal.thinning=c(), min.its=3e3, write.MCMC=TRUE, init.ages=c(), ballpark.method=2, y.scale="dates", showrun=FALSE, sep=",", normal=FALSE, delta.R=0, delta.STD=0, t.a=3, t.b=4, cc=1, cc.dir=c(), postbomb=FALSE, BCAD=FALSE, ask=FALSE, talk=TRUE, show.progress=TRUE, clean.garbage=TRUE, age.span=c(), ...) {
  start.time <- as.numeric(format(Sys.time(), "%s"))
  info <- read.strat(name, strat.dir, sep, normal, delta.R, delta.STD, t.a, t.b, cc)
  dat <- info$dets
  struc <- structure(dat) # find the structure of the data frame. Will be saved into 'info' later on
  dets <- dat[which(struc$is.age),]
  gaps <- dat[which(struc$is.gap),]

  # sanity checks
  if(talk)
    if(its > 5e6)
      message("Very many iterations, have some coffee")
  if(length(thinning) > 0)
    if(thinning < 1 || (round(thinning) != thinning) )
      stop("Thinning has to be at least 1 (and an integer)", call.=FALSE)
  if(length(internal.thinning) == 0)
     internal.thinning <- nrow(dat) else
        if(internal.thinning < 1 || (round(internal.thinning) != internal.thinning))
          stop("Internal.thinning has to be at least 1 (and an integer)", call.=FALSE)
  if(burnin < 0 || (round(burnin) != burnin))
    stop("Burnin has to be at least 0 (and an integer)", call.=FALSE)

  # set reasonable maximum age based on the dates
  if(length(age.span) == 0) {
    age.span <- max(2*(dets[,2]+dets[,3]), 55e3) # C14 or cal BP ages
    age.span <- age.span+(as.numeric(format(Sys.time(), "%Y"))-1950)
  }

  # load any calibration curve(s) we'll need
  ccs <- unique(dets[,5])
  cc.1 <- c(); cc.2 <- c(); cc.3 <- c(); cc.4 <-c()
  rc <- dets[which(dets[,5] %in% 1:4), 2] ### ok to not use dat?
  if(length(rc) > 0)
    if(min(rc) < 100) { # has postbomb dates (or close to being postbomb), which are terrestrial only. Check how rintcal's calibrate covers this
      if(1 %in% ccs) # NH
        cc.1 <- rintcal::glue.ccurves(1, postbomb, cc.dir)
      if(3 %in% ccs) # SH
        cc.3 <- rintcal::glue.ccurves(3, postbomb, cc.dir)
    } else {
        if(1 %in% ccs) # NH
          cc.1 <- rintcal::ccurve(1, FALSE, cc.dir)
        if(3 %in% ccs) # SH
          cc.3 <- rintcal::ccurve(3, FALSE, cc.dir)
      }
  if(2 %in% ccs)
    cc.2 <- rintcal::ccurve(2, FALSE, cc.dir)
  if(4 %in% ccs)
    cc.4 <- rintcal::ccurve(4, FALSE, cc.dir)

  # prepopulate non-varying values for the energy function
  dets.cc0 <- which(dets[,5] == 0) 
  has.cc0 <- ifelse(length(dets.cc0) > 0, TRUE, FALSE)
  dets.cc1 <- which(dets[,5] == 1)
  dets.cc2 <- which(dets[,5] == 2)  
  dets.cc3 <- which(dets[,5] == 3)
  dets.cc4 <- which(dets[,5] == 4)

  # prepopulate non-varying values for the support function
  xpos <- unique(dets[,4])
  is.simple <- min(diff(dets[,4])) > 0 # no blocks
  thispos <- list(); prevpos <- c(); nextpos <- c()
  if(!is.simple)
    for(i in 2:(length(xpos)-1)) {
      thispos[[i]] <- which(dets[,4] == xpos[i]) # can have multiple depths in a block
      prevpos[i] <- max(which(dets[,4] == xpos[i-1]))
      nextpos[i] <- min(which(dets[,4] == xpos[i+1]))
    }

  # sample initial age estimates to seed the MCMC (if init.ages not provided)
  if(length(init.ages) == 0)
    init.ages <- ballpark.strat(method=ballpark.method, dat=dets, gaps=struc$gaps, postbomb=postbomb, normal=normal, t.a=t.a, t.b=t.b, cc.dir=cc.dir)
  x0 <- init.ages[1,]; xp0 <- init.ages[2,]

  dets <<- dets; x0 <<- x0; info <<- info; struc <<- struc # tmp

  # calculates how well the proposed ages x fit with the dates and any constraints on ordering and gaps:
  energy <- function(x, Dets=dets, Normal=normal, ta=t.a, tb=t.b, M=struc$m,
    curves=ccs, cc1=cc.1, cc2=cc.2, cc3=cc.3, cc4=cc.4, Age.span=age.span,
    Dets.cc0=dets.cc0, Has.cc0=has.cc0, Dets.cc1=dets.cc1, Dets.cc2=dets.cc2, Dets.cc3=dets.cc3, Dets.cc4=dets.cc4,
    from=struc$from, to=struc$to, has.blocks=struc$has.blocks, blocks=struc$blocks,
    within.block=struc$within.block, above.block=struc$above.block, below.block=struc$below.block,
    has.gaps=struc$has.gaps, has.exact=struc$has.exact, has.normal=struc$has.normal, has.gamma=struc$has.gamma, has.uniform=struc$has.uniform,
    above.exact=struc$above.exact, above.normal=struc$above.normal, above.gamma=struc$above.gamma, above.uniform=struc$above.uniform,
    exact.val1=struc$exact.val1, normal.val1=struc$normal.val1, normal.val2=struc$normal.val2, gamma.val1=struc$gamma.val1,
    gamma.val2=struc$gamma.val2, uniform.val1=struc$uniform.val1, uniform.val2=struc$uniform.val2) {

    # calculcate the likelihoods of the dates
    x.cc0 <- 0
    if(Has.cc0)
      x.cc0 <- calib(Dets[Dets.cc0,2], Dets[Dets.cc0,3], x[Dets.cc0], 0, Normal, ta, tb)
    x.cc1 <- cc.energy(1, Dets, Dets.cc1, cc1, x, curves, Normal, ta, tb) # NH
    x.cc2 <- cc.energy(2, Dets, Dets.cc2, cc2, x, curves, Normal, ta, tb) # Marine
    x.cc3 <- cc.energy(3, Dets, Dets.cc3, cc3, x, curves, Normal, ta, tb) # SH
    x.cc4 <- cc.energy(4, Dets, Dets.cc4, cc4, x, curves, Normal, ta, tb) # custom ccurve

    if(has.blocks) { # then check for reversals, and adapt from and to to calculate correct l.span
	  reversal <- FALSE  
      for(i in 1:length(blocks)) {
        thisblock <- within.block[[i]]
        if (x[above.block] > min(x[thisblock]) || x[below.block] < max(x[thisblock]) ) {
          reversal <- TRUE
		  break
		} else {
		  o <- order(x[thisblock])
          to[above.block] <- thisblock[o[1]] # position of youngest age 
          from[below.block] <- thisblock[o[length(o)]] # ... and of oldest age
        } 
	  }
	  if(reversal)
        l.x <- -Inf else
          l.x <- mapply(l.span, x[from], x[to], Age.span, M)
    } else
      l.x <- mapply(l.span, x[from], x[to], Age.span, M)

    # now recalculate the energy between dates with any provided gap information
    if(has.gaps) {
      if(has.exact)
        l.x[above.exact+1] <- ifelse((x[above.exact+1] - x[above.exact]) == exact.val1, 0, -Inf)
      if(has.normal)
        l.x[above.normal+1] <- dnorm(x[above.normal+1] - x[above.normal],
          normal.val1, normal.val2, log=TRUE)
      if(has.gamma)
        l.x[above.gamma+1] <- dnorm(x[above.gamma+1] - x[above.gamma],
          gamma.val1, gamma.val2, log=TRUE)
      if(has.uniform)
        l.x[above.uniform+1] <- dunif(x[above.uniform+1] - x[above.uniform],
          uniform.val1, uniform.val2, log=TRUE)
    }

    return(-1 * sum(l.x, x.cc0, x.cc1, x.cc2, x.cc3, x.cc4))
  }

  # support function, ascertains ordering of all ages
  support <- function(ages)
    return(!min(diff(ages)) < 0)

  # write output to files if required
  out.fl <- c(); energy.fl <- c()
  if(write.MCMC) {
    out.fl <- file.path(strat.dir, name, paste0(name, "_tmp.out"))
    energy.fl <- file.path(strat.dir, name, paste0(name, "_tmp_energy.out"))
  }

  # run the twalk, but storing only output and energy (Us)
  tw <- Runtwalk(Tr=its, Obj=energy, Supp=support, dim=nrow(dets), x0=x0, xp0=xp0,
    thinning=internal.thinning, out.fl=out.fl, energy.fl=energy.fl, show.progress=show.progress)
  output <- tw$output; Us <- tw$Us;

  # post-run, remove any burn-in iterations
  if(burnin > 0) {
    output <- output[-(1:burnin),]
    Us <- Us[-(1:burnin)]
  }
  message("Removed a burn-in of ", burnin)

  # additional thinning after the run
  if(length(thinning) == 0) # by default, take the thinning value as suggested by rtwalk's IAT
    thinning <- round(max(1, coffee::IAT(tw, 0, burnin, nrow(output))))
  message("Thinning the MCMC by storing every ", thinning, " iterations")

  thinning <- seq(1, nrow(output), by=thinning)
  output <- output[thinning,]
  Us <- Us[thinning]
  if(nrow(output) < min.its)
    message(paste("Fewer than", min.its, "MCMC iterations remaining, please consider a longer run by increasing 'its'"))

  rm(tw)
  if(clean.garbage)
    gc()

  # save output and variables for future use
  write.table(output, file.path(strat.dir, name, paste0(name, ".out")), row.names=FALSE, quote=FALSE, col.names=FALSE)
  write.table(Us, file.path(strat.dir, name, paste0(name, "_energy.out")), row.names=FALSE, quote=FALSE, col.names=FALSE)
  info$dets <- dets
  info$dat <- dat
  info$struc <- struc
  info$output <- output
  info$Us <- Us
  assign_to_global("info", info)

  # draw the dates and relative information
  if(length(dat[,5] == 10) > 0)
    y.scale <- "positions"  
  dates <- draw.strat(name, info, struc, BCAD=BCAD, strat.dir=strat.dir, y.scale=y.scale, cc.dir=cc.dir, postbomb=postbomb, ybottom.lab=y.scale, ...)

  if(length(struc$pos.gaps) > 0)
    for(i in rev(struc$pos.gaps))
      struc$pos.dates[struc$pos.dates>i] <- struc$pos.dates[struc$pos.dates>i] - 1 # why this?
  overlapping <- overlap(dates$ages, dates$probs, output[,struc$pos.dates])
  info$overlap <- overlapping
  info$name <- name
  info$strat.dir <- file.path(strat.dir, name, name)
  assign_to_global("info", info) # to ensure new information is available after running the strat function

  if(talk) { 
    o <- order(overlapping)
    pos.dates <- struc$pos.dates[o]
    message("Mean overlap between calibrated and modelled dates ",
      round(mean(overlapping),2), "%, ranging from ", 
      round(min(overlapping),2), "% (date ", pos.dates[1], ") to ",
      round(max(overlapping),2), "% (date ", pos.dates[length(pos.dates)], ")")

    took <- as.numeric(format(Sys.time(), "%s")) - start.time
    if(took < 60) 
      message("Run took ", round(took, 2), " seconds") else
      if(took < 3600)
        message("Run took ", round(took/60, 2), " minutes") else
        if(took < 86400)
          message("Run took ", round(took/3600, 2), " hours") else
            message("Run took ", round(took/86400, 2), " days")
  } 
}



# internal function to find the calibrated likelihoods for specific calibration curves and values x
cc.energy <- function(j, Dets, these, cc, x, curves, Normal, ta, tb)
  if(j %in% curves) { # has dates requiring this calibration curve
    cc.y <- approx(cc[,1], cc[,2], x[these], rule=2)$y # rule=2 OK? once outside range, will float freely
    cc.er <- approx(cc[,1], cc[,3], x[these], rule=2)$y

    this.energy <- calib(Dets[these,2], Dets[these,3], cc.y, cc.er, Normal, ta, tb)
    return(sum(this.energy))
  } else # no dates which need this curve
    return(0)



# function to find the calibrated probability of x, given a calibration curve
# required for function cc.energy and also inside the energy function
calib <- function(y, er, cc.y, cc.er, normal, ta, tb)
  if(normal)
    dnorm(y, cc.y, sqrt(er^2 + cc.er^2), log=TRUE) else
      log((tb + ((y-cc.y)^2) / (2*(sqrt(cc.er^2 + er^2)^2))) ^ (-1*(ta+0.5)))



# internal function to calculate span likelihoods.
# According to Nicholls & Jones 2002 (Radiocarbon 44),
# a non-biased function for the span prior phi, given the data theta, is:
# f(phi,theta) = ( 1/ (R - d(phi)) ) ( (1/ (d_phi)^(M-1) ),
#   where R = A - P (P<=A) is the total possible span
# (by default set from 55e3 cal BP to the current year in AD),
# M is the number of events (e.g., dates, boundaries), and
# d(phi) = phi_0 - phi_m is the (modelled) span between the boundary ages phi_m and phi_0.
# f(phi,theta) = ( 1/ (R - d(phi)) ) ( (1/ (d_phi)^(M-1) )
# as log: -log((R - dphi) - ((M-1) * log(dphi)))
l.span <- function(max.x, min.x, R, M) {
  dphi <- abs(max.x - min.x)
  return(-log((R - dphi) - ((M-1) * log(dphi))))
}



#' @name ages.undated
#' @title Model ages between two dated levels
#' @description Model ages of undated levels, by for each MCMC iteration finding the age of the layer above and of the layer below, and sampling a random age from a uniform distribution between the age estimates of the two ages. 
#' @param position Position of the to-be-estimated undated layer. Should be larger than the layer above but smaller than the layer below it.
#' @param set This option reads the 'info' variable, which contains the data and the model output.
#' @param draw Whether or not to draw the age distribution.
#' @export
ages.undated <- function(position, set=get('info'), draw=TRUE) {
  positions <- which(set$dets[,5] %in% 0:4) # dated levels
  above <- max(which(set$dets[positions,4] < position))
  below <- min(which(set$dets[positions,4] > position))
  ages.above <- set$output[,above]
  ages.below <- set$output[,below]
  ages <- runif(1:nrow(set$output), ages.above, ages.below)
  if(draw)
    plot(density(ages), main="", xlab="cal BP", ylab="")
  return(ages)
}



# internal function to calculate how well the modelled distributions fit with the calibrated ones
overlap <- function(calibs, probs, modelled, n=1e3) {
  same <- c()
  for(i in 1:ncol(modelled)) {
    calib <- sample(calibs[,i], n, replace=TRUE, prob=probs[,i])
    model <- sample(modelled[,i], n, replace=TRUE)
    smaller <- length(which(calib-model < 0))
    larger <- length(which(calib-model > 0))
    same[i] <- 1-abs(smaller-larger)/n
  }
  return(100*same) # as percentage
}
