
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
#' @param clean.garbage Whether or not to clean up the memory 'garbage collection' after a run. Recommendable if you have many dates or long runs.
#' @param oldest.age Oldest expected ages. Defaults to 55e3 which is the current cal BP limit for C-14 dates. If older, non-14C dates are present, oldest.age is set to the larger of the radiocarbon limit or twice the age of the oldest non-radiocarbon age. 
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
strat <- function(name="mystrat", strat.dir="strats", its=5e4, burnin=100, thinning=c(), internal.thinning=c(), write.MCMC=TRUE, init.ages=c(), ballpark.method=2, y.scale="dates", showrun=FALSE, sep=",", normal=FALSE, delta.R=0, delta.STD=0, t.a=3, t.b=4, cc=1, cc.dir=c(), postbomb=FALSE, BCAD=FALSE, ask=FALSE, talk=TRUE, clean.garbage=TRUE, oldest.age=c(), ...) {	
  start.time <- as.numeric(format(Sys.time(), "%s"))
  info <- read.strat(name, strat.dir, sep, normal, delta.R, delta.STD, t.a, t.b, cc)
  assign_to_global("info", info)
  dat <- info$dets

  # set reasonable maximum age based on the dates
  if(length(oldest.age) == 0) {
    these <- which(dat[,5] %in% 0:4) # only dates
    oldest.age <- max(2*(dat[these,2]+dat[these,3]), 55e3) # C14 or cal BP ages
  }

  # sanity checks
  if(talk) 
    if(its > 5e6)
      message("Very many iterations, I suggest coffee")
  if(length(thinning) > 0)
    if(thinning < 1)
      stop("Thinning has to be at least 1 (and an integer)", call.=FALSE)
  if(length(internal.thinning) == 0)
	  internal.thinning <- nrow(dat) else
        if(internal.thinning < 1)
          stop("Internal.thinning has to be at least 1 (and an integer)", call.=FALSE)
  if(burnin < 0)
    stop("Burnin has to be at least 0 (and an integer)", call.=FALSE)
  
  # prepare for later calculations of ages, dates and gaps
  pos.ages <- which(dat[,5] %in% c(0:4,10)) # dates and undated levels for which age estimates are required
  pos.dates <- which(dat[,5] %in% 0:4) # only dates
  pos.gaps <- which(dat[,5] > 10) # exact gaps, or normal, gamma or uniform
  dets <- dat[pos.ages,]
  gaps <- c()
  if(length(pos.gaps) > 0)
    gaps <- dat[pos.gaps,]

  # sample initial age estimates to seed the MCMC (if init.ages not provided)
  if(length(init.ages) == 0)
    init.ages <- ballpark.strat(method=ballpark.method, dat=dets, gaps=gaps, postbomb=postbomb, normal=normal, t.a=t.a, t.b=t.b, cc.dir=cc.dir)
  x0 <- init.ages[1,]; xp0 <- init.ages[2,]

  # load any calibration curve(s) we'll need
  ccs <- unique(dets[,5])
  cc.1 <- c(); cc.2 <- c(); cc.3 <- c(); cc.4 <-c()
  rc <- dets[which(dets[,5] %in% 1:4), 2]

  if(length(rc) > 0)
    if(min(rc) < 0) { # has postbomb dates, which are terrestrial only
      if(1 %in% ccs) # NH
        cc.1 <- rintcal::glue.ccurves(1, postbomb, cc.dir)
      if(3 %in% ccs) # SH
        cc.3 <- rintcal::glue.ccurves(3, postbomb, cc.dir)
    } else {
        if(1 %in% ccs) # NH
          cc.1 <- rintcal::ccurve(1, FALSE, cc.dir)
        if(3 %in% ccs) # SH; Marine doesn't have postbomb curves
          cc.3 <- rintcal::ccurve(3, FALSE, cc.dir)
      }

  if(2 %in% ccs)
    cc.2 <- rintcal::ccurve(2, FALSE, cc.dir)
  if(4 %in% ccs)
    cc.4 <- rintcal::ccurve(4, FALSE, cc.dir)

  # find the calibrated likelihoods
  calib <- function(y, er, cc.y, cc.er, normal, ta, tb)
    if(normal)
      dnorm(y, cc.y, sqrt(er^2 + cc.er^2)) else
        (tb + ((y-cc.y)^2) / (2*(sqrt(cc.er^2 + er^2)^2))) ^ (-1*(ta+0.5))

  # According to Nicholls & Jones 2002 (Radiocarbon 44), a non-biased function
  # for the span prior phi, given the data theta, is:
  # f(phi,theta) = ( 1/ (R - d(phi)) ) ( (1/ (d_phi)^(M-1) ),
  #   where R = A - P (P<=A) is the total possible span (by default set at oldest.age=55e3 cal BP to AD 2023,
  #   the current limits of C-14 calibration), M is the number of events (e.g., dates, boundaries), and
  #   d(phi) = phi_0 - phi_m is the (modelled) span between the boundary ages phi_m and phi_0.
  # f(phi,theta) = ( 1/ (R - d(phi)) ) ( (1/ (d_phi)^(M-1) )
  # as log: -log((R - dphi) - ((M-1) * log(dphi))) # (M-1) was M
  l.span <- function(max.x, min.x, R, M) {
    if(length(R) == 0)
      R <- 55e3+(as.numeric(format(Sys.time(), "%Y"))-1950) # C14 'realm'
    dphi <- abs(max.x - min.x)
    return(-log((R - dphi) - ((M-1) * log(dphi)))) # M-1 was M
  }

  # calculates how well the proposed ages x fit with the dates and any constraints on ordering and gaps:
  energy <- function(x, Dets=dets, curves=ccs, cc1=cc.1, cc2=cc.2, cc3=cc.3, cc4=cc.4, Normal=normal, ta=t.a, tb=t.b) {
    # calculcate the likelihoods of the dates
    if(0 %in% curves) { # no calibration necessary
      these <- which(Dets[,5] == 0)
      x.cc0 <- calib(Dets[these,2], Dets[these,3], x[these], 0, Normal, ta, tb)
      x.cc0 <- log(x.cc0)
    } else
        x.cc0 <- 0

    cc.energy <- function(j, cc, x)
      if(j %in% curves) { # has dates requiring this calibration curve
        these <- which(Dets[,5] == j)
        cc.y <- approx(cc[,1], cc[,2], x[these], rule=2)$y # rule=2 OK? once outside range, will float freely
        cc.er <- approx(cc[,1], cc[,3], x[these], rule=2)$y
        this.energy <- calib(Dets[these,2], Dets[these,3], cc.y, cc.er, Normal, ta, tb)
        return(sum(log(this.energy)))
      } else # no dates which need this curve
        return(0)

    x.cc1 <- cc.energy(1, cc1, x) # IntCal20
    x.cc2 <- cc.energy(2, cc2, x) # Marine20
    x.cc3 <- cc.energy(3, cc3, x) # SHCal20
    x.cc4 <- cc.energy(4, cc4, x) # custom ccurve

    # now the contextual likelihoods: order, undated levels, time spans (exact, gamma, normal)
    # the span of the uppermost, youngest age should be unrestricted within the span prior
    # It has to be able to move more than the other dates, which will all be expressed as gaps from their younger neighbour
    l.x <- rep(0, nrow(Dets))
    l.x[1] <- l.span(max(x, na.rm=TRUE), min(x, na.rm=TRUE), oldest.age, length(x))

    # now the rest, relative to the first one
    for(i in 2:nrow(Dets)) {
      younger <- max(1, which(Dets[,4] < Dets[i,4])) # find the oldest date that is younger than x[i]. Deals with dates in 'blocks'
      older <- min(nrow(Dets), which(Dets[,4] > Dets[i,4])) # youngest age older than x[i]

    #### l.span should consider only 2 or 3 ages here (x_i and boundaries x_min, x_max), not all x
      if(Dets[i,5] %in% 0:4)
        l.x[i] <- l.span(x[i], x[younger], oldest.age, 2) else
          l.x[i] <- l.span(x[younger], x[older], oldest.age, 3) # undated
    }

    # now fill any gaps, i.e. time spans
    if(length(pos.gaps) > 0) {
      for(i in 1:nrow(gaps)) {
        younger <- max(1, which(Dets[,4] < gaps[i,4]))

        # exactly known gap (e.g., n tree rings); value in column 2
        if(gaps[i,5] == 11) {
          gap <- (x[younger+1] - x[younger]) - gaps[i,2]
          l.x[younger+1] <- log(ifelse(gap == 0, 1, 0))
        }

        # normally distributed time span; mean in column 2, sdev in column 3
        if(gaps[i,5] == 12)
          l.x[younger+1] <- log(dnorm(x[younger+1]-x[younger], gaps[i,2], gaps[i,3]))

        # gamma; mean in column 2, shape in column 3
        if(gaps[i,5] == 13)
          l.x[younger+1] <- log(dgamma(x[younger+1]-x[younger], gaps[i,3], gaps[i,3]/gaps[i,2]))

        # uniform; minimum in column 2, maximum in column 3. Perhaps leave out as no known use
        if(gaps[i,5] == 14)
          l.x[younger+1] <- log(dunif(x[younger+1]-x[younger], gaps[i,2], gaps[i,3]))
      }
    }
  return(-1 * sum(l.x, x.cc0, x.cc1, x.cc2, x.cc3, x.cc4))
  }

  # check that the ages x comply with chronological ordering
    support <- function(ages, position=dets[,4]) {
    OK <- TRUE
    ages <- ages[!is.na(ages)]

    if(min(diff(position)) > 0) { # then all in order, conventional case
      if(min(diff(ages)) < 0)
        OK <- FALSE
    } else { # more complicated case
        xpos <- unique(position)
        for(i in 2:(length(xpos)-1)) {
          thispos <- which(position == xpos[i])
          prevpos <- which(position == xpos[i-1])
          nextpos <- which(position == xpos[i+1])
          if(min(ages[thispos]) < max(ages[prevpos]) ||
            max(ages[thispos]) > min(ages[nextpos])) {
              OK <- FALSE
              break # then no need to continue this loop
          }
        }
      }
    return(OK)
  }

  # run the twalk with the above energy and support, two sets of initial estimates, and its iterations
  out.fl <- c(); energy.fl <- c()
  if(write.MCMC) {
	  out.fl <- file.path(strat.dir, name, paste0(name, "_tmp.out"))
	  energy.fl <- file.path(strat.dir, name, paste0(name, "_tmp_energy.out"))
  }
  info <- Runtwalk(Tr=its, Obj=energy, Supp=support, dim=nrow(dets), x0=x0, xp0=xp0, thinning=internal.thinning, out.fl=out.fl, energy.fl=energy.fl)
  if(burnin > 0) {
    info$output <- info$output[-(1:burnin),]
    info$Us <- info$Us[-(1:burnin)]
  }
  message("Removed a burn-in of ", burnin)

  # additional thinning after the run
  if(length(thinning) == 0) # by default, take the thinning value as suggested by rtwalk's IAT
    thinning <- round(max(1, coffee::IAT(info, 0, burnin, info$k - burnin)))
  message("Thinning the MCMC by storing every ", thinning, " iterations")

  thinning <- seq(1, nrow(info$output), by=thinning)
  info$output <- info$output[thinning,]
  info$Us <- info$Us[thinning]
  if(nrow(info$output) < 3e3)
    message("Fewer than 3k MCMC iterations remaining, please consider a longer run by increasing 'its'")

  if(clean.garbage)
    gc()

  # save output for future use
  write.table(info$output, file.path(strat.dir, name, paste0(name, ".out")), row.names=FALSE, quote=FALSE, col.names=FALSE)
  write.table(info$Us, file.path(strat.dir, name, paste0(name, "_energy.out")), row.names=FALSE, quote=FALSE, col.names=FALSE)

  info$dets <- dat
  assign_to_global("info", info)

  if(length(dat[,5] == 10) > 0)
    y.scale <- "positions"  

  dates <- draw.strat(name, info, BCAD=BCAD, strat.dir=strat.dir, y.scale=y.scale, cc.dir=cc.dir, postbomb=postbomb, ybottom.lab=y.scale, ...)

  if(length(pos.gaps) > 0)
    for(i in rev(pos.gaps))
	  pos.dates[pos.dates>i] <- pos.dates[pos.dates>i] - 1  
  overlapping <- overlap(dates$ages, dates$probs, info$output[,pos.dates])
  info$overlap <- overlapping
  assign_to_global("info", info)

  if(talk) { 
	o <- order(overlapping)
	pos.dates <- pos.dates[o]  
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
