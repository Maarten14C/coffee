


#' @name strat
#' @title Model chronologically ordered dates
#' @description Model radiocarbon dates (or dates that are already on the cal BP scale) of a deposit that is known to have accumulated over time, and for which therefore the dated depths can be safely assumed to are in chronological order.
#' @details The calculations are made in a Bayesian framework. For each iteration, it is checked that all modelled age estimates are in chronological order - if not, the iteration is not accepted. See also Buck et al. 1991 and Nicholls & Jones 2001. Other software that enables Bayesian chronological ordering includes Bcal (Buck et al. 1999) and OxCal (Sequence model; Bronk Ramsey 1995).
#' @description Model the radiocarbon or cal BP ages of a deposit that is known to have accumulated over time, and for which therefore the dated depths can be safely assumed to be in chronological order.
#' @details Dates further down the sequence should have older ages than dates further up, even though owing to scatter, the dates themselves might not be in exact chronological order. The amount of scatter, the laboratory error and an offset can also be modelled.
#' The age estimates are obtained through a t-walk MCMC run (Christen and Fox 2010). In this process, initial ball-park point estimates for the ages of each dated depth are given, and then modified through many iterations. For each iteration, a random dated depth is chosen and its age changed by just a little nudge, a check is performed to ensure that all age estimates remain in chronological order, and the 'energy' or likelihood of the age estimates is calculated (iterations where all ages fit well within the calibrated distributions receive a higher energy; see \code{l.calib}). 
#' Then this iteration with the updated group of age estimates is either accepted or rejected. The acceptance probability depends on the iteration's energy; if its energy is higher than that of the previous iteration it is accepted, but if it is lower, it is accepted with a probability proportional to its relative energy. Therefore, over many iterations the process will 'learn' from the data and find high-energy combinations of parameter values that fit with the prior constraints that the ages should be ordered chronologically.
#' Because the iterations are based on a process of modifying values of one parameter each iteration, and because some iterations will not be accepted, the MCMC output will often have a large degree of dependence between neighbouring iterations. Therefore, some thinning will have to be done, by storing only one every few iterations (default 20). Also, since the initial ball-park estimates could be quite wrong, the first 100 or so iterations should also be discarded (burnin). 
#' It is thus important to check the time-series of the energy after the run. We don't want to see a remaining burn-in at the start, and we don't want to see a noticeable 'structure' where iterations remain in approximately or entirely the same spot for a long time. Instead, an ideal run will look like white noise.
#' @param name Name of the stratigraphy dataset. Defaults to \code{"mystrat"}.
#' @param strat.dir The directory where the folders of the individual stratigraphies live. Defaults to \code{treedir="strats"}.
#' @param its Amount of iterations to be run. Setting this to low numbers (e.g., 1000) will result in fast but less stable and less reliable runs. Higher values will take longer but result in more stable and robust runs. Defaults to \code{50000}. Aim to set this to such values that at least 3000 iterations remain after removing the burnin and thinning.
#' @param burnin Amount of iterations to remove at the start of the run. Defaults to \code{100}.
#' @param thinning After running all iterations, only some will be stored. For example, if thinning is set at the default \code{50}, only every 50th MCMC iteration will be stored, and the others will be discarded. This is to remove the dependence between neighbouring MCMC iterations. Defaults to a value calculated from the MCMC run itself.
#' @param init.ages By default, the ballpark age estimates to feed the MCMC are calculated automatically, however they can also be provided manually.
#' @param span Extent by which the uniform prior should expand beyond the youngest and oldest initial age estimates. Defaults to \code{span=5000}. 
#' @param showrun Whether or not to show how the MCMC process is progressing during the run. Defaults to \code{FALSE}.
#' @param sep Separator for the fields in the .csv file. Defaults to a comma.
#' @param normal Calculations can be done assuming that the measurements are normally distributed. By default this is set to FALSE and a student-t distribution is used (Christen and Perez 2009)
#' @param delta.R The ages can be modelled to have an offset. The mean is 0 by default.
#' @param delta.STD The error of the offset. Set to 0 by default.
#' @param t.a First parameter for the student-t distribution (defaults to 3; higher numbers make the distribution approximate the normal distribution more).
#' @param t.b Second parameter for the student-t distribution (defaults to 4; higher numbers make the distribution 
#' @param cc Calibration curve to be used. Could be 1 (IntCal20; default), 2 (Marine20), 3 (SHCal20) or 4 (custom curve).
#' @param postbomb Negative C-14 ages should be calibrated using a postbomb curve. This could be 1 (northern-hemisphere region 1), 2 (NH region 2), 3 (NH region 3), 4 (southern hemisphere regions 1-2), or 5 (SH region 3).
#' @param BCAD The calendar scale of graphs and age output-files is in \code{cal BP} by default, but can be changed to BC/AD using \code{BCAD=TRUE}.
#' @param ask Whether or not to ask if a folder should be made (if required).
#' @param talk Whether or not to provide feedback on folders written into.
#' @param draw Whether or not to draw plots.
#' @param ... Options for the plot. See \code{plot.strat}.
#' @return a variable 'info' which contains the dating and modelling information to produce a plot. Also calls the function \code{draw.strat} to produce a plot of the results.
#' @examples
#' \dontrun{
#' tmp <- tempdir()
#' strat(, strat.dir=tmp)
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
strat <- function(name="mystrat", strat.dir="strats", its=5e4, burnin=100, thinning=c(), init.ages=c(), span=5e3, showrun=FALSE, sep=",", normal=TRUE, delta.R=0, delta.STD=0, t.a=3, t.b=4, cc=1, postbomb=FALSE, BCAD=FALSE, ask=TRUE, talk=TRUE, draw=TRUE, ...) {
  stratdir <- assign_dir(strat.dir, name, "strat.dir", ask, talk)
  dat <- read.table(file.path(stratdir, paste0(name, ".csv")), header=TRUE, sep=sep)

  # file sanity checks
  OK <- TRUE
  if(ncol(dat) > 5 || ncol(dat) < 4) OK <- FALSE # only 4 or 5 columns
  if(min(diff(dat[,4])) < 0) OK <- FALSE # we need increasing stratigraphical positions - larger numbers are further down
  if(min(dat[,5]) < 0 || max(dat[,5]) > 4) OK <- FALSE # can use C14 or cal BP dates
  if(!OK)
    stop("Unexpected values in strat file. Please check", call.=FALSE)
  if(length(thinning) > 0)
    if(thinning < 1)
      stop("Thinning has to be at least 1 (and an integer)", call.=FALSE)
  if(burnin < 0)
    stop("Burnin has to be at least 0 (and an integer)", call.=FALSE)

  # age offsets and t.a/t.b can also be provided
  if(ncol(dat) == 5) {
    dat[,2] <- dat[,2] - delta.R
    dat[,3] <- sqrt(dat[,3]^2 + delta.STD^2)
  } else
    if(ncol(dat) > 5)
      if(tolower(colnames(dat))[6] %in% c("delta.r", "delta", "offset", "reservoir")) {
        dat[,2] <- dat[,2] - dat[,6]
        dat[,3] <- sqrt(dat[,3]^2 + dat[,7]^2)
      }
    if(ncol(dat) == 9) {
      t.a <- dat[,8]
      t.b <- dat[,9]
    }


  # find initial age estimates to seed the MCMC. Assuming NH ccurve but shouldn't matter because these are initial points only
  calcurve <- ccurve(cc, postbomb)
  min.y <- which(dat[,2] == min(dat[,2]))[1]
  max.y <- which(dat[,2] == max(dat[,2]))[1]
  if(dat[min.y,5] == 0)
    min.y <- dat[min.y,2] else 
      min.y <- calcurve[max(1,which(calcurve[,2] <= dat[min.y,2])),1]
  if(dat[max.y,5] == 0)
    max.y <- dat[max.y,2] else
      max.y <- calcurve[min(nrow(calcurve), which(calcurve[,2] >= dat[max.y,2])),1]
  if(length(init.ages) == 0)
    init.ages <- sort(runif(nrow(dat), min.y, max.y))
  x0 <- sort(jitter(init.ages)) # initial ball-park age estimates
  xp0 <- sort(jitter(init.ages))  # twalk needs two sets of estimates
  
  # load any calibration curve(s) we'll need
  ccs <- unique(dat[,5])
  cc.1 <- c(); cc.2 <- c(); cc.3 <- c(); cc.4 <-c()
  if(min(dat[which(dat[,5] > 0), 2]) < 0) { # has postbomb dates, which are terrestrial only
    if(1 %in% ccs) # NH
      cc.1 <- IntCal::glue.ccurves(1, postbomb)
    if(3 %in% ccs) # SH
      cc.3 <- IntCal::glue.ccurves(3, postbomb)
  } else {
      if(1 %in% ccs) # NH
        cc.1 <- IntCal::ccurve(1, FALSE)
      if(3 %in% ccs) # SH
        cc.3 <- IntCal::ccurve(3, FALSE)
    }
  if(2 %in% ccs)
    cc.2 <- IntCal::ccurve(2, FALSE)
  
  # find the calibrated likelihoods
  calib <- function(y, er, cc.y, cc.er, normal, ta, tb)
    if(normal)
      dnorm(y, cc.y, sqrt(er^2 + cc.er^2)) else
        (tb + ((y-cc.y)^2) / (2*(sqrt(cc.er^2 + er^2)^2))) ^ (-1*(ta+0.5))

  # (re)define the functions relevant for Runtwalk inline

  # uniform prior for the total span
  # to deal with prior problem reported by Steier & Rom 2000 and Nicholls & Jones 2001. OK?
  # According to Nicholls & Jones 2002 (Radiocarbon 44), a non-biased function for the span prior phi, given the data theta, is:
  # f(phi,theta) = ( 1/ (R - d(phi)) ) ( 1/ (d_phi)^M-1 )
  # where R = A - P (P<=A) is the possible span (set at 55e3, the current limits of C-14 calibration), M is the number of dates/events, and d(phi) = phi_0 - phi_m is the (modelled) span between the boundary dates phi_m and phi_0.
  # This becomes:
  # span <- -log( (55e3 - dphi))  - (length(x)-2) * log( dphi)

  # span <- log((1/((max(init.ages) - min(init.ages)) + span))^(length(init.ages)-2))

  energy <- function(x, dets=dat, curves=ccs, cc1=cc.1, cc2=cc.2, cc3=cc.3, cc4=ccurve(4), Normal=normal, ta=t.a, tb=t.b) {
    dphi <- max(x) - min(x)
    span <- -log( (55e3 - dphi))  - (length(x)-2) * log( dphi) # adapted per advice by Andres Christen and Marco Aquino

    if(0 %in% curves) {
      dat <- dets[which(dets[,5] == 0),]
      x0 <- log(calib(dat[,2], dat[,3], x, 0, Normal, ta, tb))
    } else x0 <- 0

    cc.energy <- function(i, cc, x)
      if(i %in% curves) {
        these <- which(dets[,5] == i)
        cc.y <- approx(cc[,1], cc[,2], x[these])$y
        cc.er <- approx(cc[,1], cc[,3], x[these])$y
        return(log(calib(dets[these,2], dets[these,3], cc.y, cc.er, Normal, ta, tb)))
      } else
          return(0)

    x1 <- cc.energy(1, cc1, x)
    x2 <- cc.energy(2, cc2, x)
    x3 <- cc.energy(3, cc3, x)
    x4 <- cc.energy(4, cc4, x)

    return(-1 * sum(span, x0, x1, x2, x3, x4))
  }

  support <- function(ages)
    if(min(diff(ages)) > 0)
      return(TRUE) else
        return(FALSE)

  info <- Runtwalk(Tr=its, Obj=energy, Supp=support, x0=x0, xp0=xp0, PlotLogPost=ifelse(showrun, TRUE, FALSE) )

  if(burnin > 0) {
    info$output <- info$output[-(1:burnin),]
    info$Us <- info$Us[-(1:burnin)]
    #info$Tr <- info$Tr - burnin
  }
  message("\nremoved a burn-in of ", burnin)
  if(length(thinning) == 0) # by default, take the thinning value as suggested by rtwalk's IAT
    thinning <- round(max(1, coffee::IAT(info, 0, burnin, info$Tr - burnin)))
  message("thinning the MCMC by storing every by ", thinning, " iterations")
  thinning <- seq(1, its-burnin, by=thinning)
  info$output <- info$output[thinning,]
  info$Us <- info$Us[thinning]
  if(length(thinning) < 3e3)
    message("Fewer than 3k MCMC iterations remaining, please consider a longer run by increasing its")
  # save output for future use
  write.table(info$output, file.path(strat.dir, name, paste0(name, ".out")), row.names=FALSE, quote=FALSE, col.names=FALSE)
  write.table(info$Us, file.path(strat.dir, name, paste0(name, "_energy.out")), row.names=FALSE, quote=FALSE, col.names=FALSE)

  info$dets <- dat
  
  assign_to_global <- function(x, value, pos=1)
    assign(x, value, envir=as.environment(pos) )
  assign_to_global("info", info)

  if(draw)
    draw.strat(name, info, BCAD=BCAD, strat.dir=stratdir, ...)
}



#' @name draw.strat
#' @title plot the dates and model of chronologically ordered dated depths 
#' @description A plot with two panels. The top panel shows the MCMC output. The bottom panel shows the individually calibrated dates (in downward light gray) as well as the modelled ages constrained by chronological ordering (upward dark-grey) and lines with the hpd ranges (black).
#' @param name Name of the stratigraphy dataset. Defaults to \code{"mystrat"}.
#' @param set This option reads the 'info' variable, which contains the data and the model output. 
#' @param strat.dir The directory where the folders of the individual stratigraphies live. Defaults to \code{start.dir="strats"}.
#' @param sep Separator for the fields in the .csv file. Defaults to a comma.
#' @param calibrated.ex Exaggeration of the heights of the calibrated distributions. Defaults to 0.5. Note that more precise dates peak higher than dates with lower precision.
#' @param calibrated.mirror Whether or not the individually calibrated (but not the modelled) distributions should be drawn both up and down, quite a bit like fish or swans. Defaults to FALSE.
#' @param calibrated.up Whether the calibrated distributions should be drawn upward or downward (the default, resembling the reflections of islands in the sea, or swimming animals if you wish)
#' @param modelled.ex Exaggaration of the heights of the age-modelled distributions. Defaults to 0.5. Note that more precise ages peak higher than ages with lower precision.
#' @param modelled.mirror Whether or not the age-modelled distributions should be drawn both up and down, quite a bit like fish or swans. Defaults to FALSE.
#' @param modelled.up Whether the age-modelled distributions should be drawn downward or upward (the default, resembling islands in the sea)
#' @param BCAD The calendar scale of graphs and age output-files is in \code{cal BP} by default, but can be changed to BC/AD using \code{BCAD=TRUE}.
#' @param threshold Value below which probabilities should not be drawn any more. Defaults to 0.001 of the distribution's peak.
#' @param xtop.lab The label for the x-axis of the top panel showing the MCMC run. Defaults to \code{"iterations"}.
#' @param ytop.lab The label for the y-axis of the top panel showing the MCMC run. Defaults to \code{"energy"}.
#' @param xbottom.lab The label for the x-axis of the bottom panel showing the age-model output. Defaults to \code{"cal BP"} or \code{"BC/AD"}.
#' @param ybottom.lab The label for the y-axis of the bottom panel showing the age-model output. Defaults to \code{"position"}.
#' @param calibrated.col Colour of the inside of the unmodelled, calibrated ages. Defaults to semi-transparent light grey, \code{rgb(0,0,0,0.5}.
#' @param calibrated.border Colour of the border of the unmodelled, calibrated ages. Defaults to nothing, NA.
#' @param modelled.col Colour of the inside of the modelled ages. Defaults to semi-transparent dark grey, \code{rgb(0,0,0,0.5}.
#' @param modelled.border Colour of the border of the modelled ages. Defaults to semi-transparent dark grey, \code{rgb(0,0,0,0.5}.
#' @param range.col Colour of the hpd ranges. Defaults to \code{"black"}.
#' @param simulation Whether or not the data are part of a simulated stratigraphy.
#' @param simulation.col If the data are part of a simulated stratigraphy, the 'true' ages can also be plotted.
#' @param mgp Axis text margins (where should titles, labels and tick marks be plotted). Defaults to \code{mgp=c(1.7, .7, .0)}.
#' @param mar.top Margins around the top panel. Defaults to \code{mar.top=c(3,3,1,1)}.
#' @param mar.bottom Margins around the bottom panel. Defaults to \code{mar.bottom=c(3,3,0.5,1)}.
#' @param heights Relative heights of the two panels in the plot. Defaults to 0.3 for the top and 0.7 for the bottom panel.
#' @param iterations.warning Whether or not to plot a warning if there are < 3000 iterations, too few for a reliable MCMC run.
#' @param warning.loc Location of the warning
#' @param warning.col Colour of the warning - defaults to red.
#' @return A plot with two panels showing the MCMC run and the calibrated and modelled ages.
#' @author Maarten Blaauw
#' @export
draw.strat <- function(name="mystrat", set=get('info'), strat.dir="strats", sep=",", calibrated.ex=.5, calibrated.mirror=FALSE, calibrated.up=TRUE, modelled.ex=0.5, modelled.mirror=FALSE, modelled.up=FALSE, BCAD=FALSE, threshold=0.001, xtop.lab=c(), ytop.lab=c(), xbottom.lab=c(), ybottom.lab="position", calibrated.col=rgb(0, 0, 0, 0.2), calibrated.border=NA, modelled.col=rgb(0,0,0,0.5), modelled.border=rgb(0,0,0,0.5), range.col="black", simulation=FALSE, simulation.col=grey(0.5), mgp=c(2, 0.7, 0), mar.top=c(3,3,1,1), mar.bottom=c(3,3,0.5,1), heights=c(0.3, 0.7), iterations.warning=TRUE, warning.loc="bottomleft", warning.col="red") {
  layout(matrix(1:2, ncol=1), heights=heights)
  oldpar <- par(no.readonly=TRUE)
  on.exit(par(oldpar))
  par(mar=mar.top, mgp=mgp)

  # read the data and model output
  if(length(set) == 0) {
    ages <- read.table(file.path(strat.dir, name, paste0(name, ".out")), header=FALSE)
    energy <- read.table(file.path(strat.dir, name, paste0(name, "_energy.out")), header=FALSE)
    dets <- read.table(file.path(strat.dir, name, paste0(name, ".csv")), header=TRUE, sep=",")
  } else {
      energy <- set$Us
      ages <- set$output
      dets <- set$dets
    }

  if(length(xtop.lab) == 0)
    xtop.lab <- "iterations"
  if(length(ytop.lab) == 0)
    ytop.lab <- "energy"
  plot(-1*energy, type="l", xlab=xtop.lab, ylab=ytop.lab)
  if(iterations.warning)
    if(nrow(set$output) < 3e3) # then MCMC run likely not long enough to be reliable
      legend(warning.loc, legend="short MCMC run - needs more iterations", bty="n", cex=.8, text.col=warning.col)

  par(mar=mar.bottom)
  draw.dates(dets[,2], dets[,3], dets[,4], dets[,5], ex=calibrated.ex, mirror=calibrated.mirror, up=calibrated.up, col=calibrated.col, border=calibrated.border, BCAD=BCAD, draw.hpd=FALSE, threshold=threshold, normalise=TRUE, cal.lab=xbottom.lab, y.lab=ybottom.lab)

  maxdens <- 0
  for(i in 1:ncol(ages))
    maxdens <- max(maxdens, density(ages[,i])$y)

  if(BCAD)
    ages <- 1950 - ages

  hpds <- list()
  for(i in 1:ncol(ages)) {
    age <- density(ages[,i])
    age$y <- age$y / maxdens
    hpds[[i]] <- IntCal::hpd(cbind(age$x, age$y))

    if(modelled.mirror)
      pol <- cbind(c(age$x, rev(age$x)), dets[i,4]-modelled.ex*c(age$y, -rev(age$y))) else
        if(modelled.up)
          pol <- cbind(c(min(age$x), age$x, max(age$x)), dets[i,4]+modelled.ex*c(0, age$y, 0)) else
            pol <- cbind(c(min(age$x), age$x, max(age$x)), dets[i,4]-modelled.ex*c(0, age$y, 0))
    polygon(pol, col=modelled.col, border=modelled.border)
    if(simulation)
      segments(ifelse(BCAD, 1950-dets[i,1], dets[i,1]), i, ifelse(BCAD, 1950-dets[i,1], dets[i,1]), ncol(ages)+10, lty=2, col=simulation.col) # then don't be afraid to show the truth!
    if(!is.na(range.col))
      segments(hpds[[i]][,1], i, hpds[[i]][,2], i, col=range.col, lwd=2)
  }
   set$hpds <- hpds
   
   assign_to_global("info", set)
}


