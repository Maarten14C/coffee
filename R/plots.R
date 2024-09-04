#' @name draw.strat
#' @title plot the dates and model of chronologically ordered dated depths
#' @description A plot with two panels. The top panel shows the MCMC output. The bottom panel shows the individually calibrated dates (in downward light gray) as well as the modelled ages constrained by chronological ordering (upward dark-grey) and lines with the hpd ranges (black). Any similarity with swimming elephants or island chains is coincidental.
#' @param name Name of the stratigraphy dataset. Defaults to \code{"mystrat"}.
#' @param set This option reads the 'info' variable, which contains the data and the model output.
#' @param structure Information about the structure (e.g., blocks, gaps, dates, undated levels) of the dataset. Filled automatically.
#' @param y.scale The scale of the vertical axis of the main plot. This can be the positions of the dated levels (`positions`) or their position order (`dates`). 
#' @param strat.dir The directory where the folders of the individual stratigraphies live. Defaults to \code{start.dir="strats"}.
#' @param cc.dir Directory of calibration curve. Keep empty for the default value.
#' @param sep Separator for the fields in the .csv file. Defaults to a comma.
#' @param postbomb Negative C-14 ages should be calibrated using a postbomb curve. This could be 1 (northern-hemisphere region 1), 2 (NH region 2), 3 (NH region 3), 4 (southern hemisphere regions 1-2), or 5 (SH region 3).
#' @param calibrated.ex Exaggeration of the heights of the calibrated distributions. Calculated automatically by default. Note that more precise dates peak higher than dates with lower precision.
#' @param calibrated.mirror Whether or not the individually calibrated (but not the modelled) distributions should be drawn both up and down, quite a bit like fish or swans. Defaults to FALSE.
#' @param calibrated.up Whether the calibrated distributions should be drawn upward or downward (the default, resembling the reflections of islands in the sea, or swimming animals if you wish)
#' @param modelled.ex Exaggeration of the heights of the age-modelled distributions. Calculated automatically by default. Note that more precise ages peak higher than ages with lower precision.
#' @param modelled.mirror Whether or not the age-modelled distributions should be drawn both up and down, quite a bit like fish or swans. Defaults to FALSE.
#' @param modelled.up Whether the age-modelled distributions should be drawn downward or upward (the default, resembling islands in the sea)
#' @param BCAD The calendar scale of graphs and age output-files is in \code{cal BP} by default, but can be changed to BC/AD using \code{BCAD=TRUE}.
#' @param threshold Value below which probabilities should not be drawn any more. Defaults to 0.001 of the distribution's peak.
#' @param xtop.lab The label for the x-axis of the top panel showing the MCMC run. Defaults to \code{"iterations"}.
#' @param ytop.lab The label for the y-axis of the top panel showing the MCMC run. Defaults to \code{"energy"}.
#' @param xbottom.lab The label for the x-axis of the bottom panel showing the age-model output. Defaults to \code{"cal BP"} or \code{"BC/AD"}.
#' @param ybottom.lab The label for the y-axis of the bottom panel showing the age-model output. Defaults to \code{"position"}.
#' @param calibrated.col Colour of the inside of the unmodelled, calibrated ages. Defaults to semi-transparent light grey, \code{rgb(0,0,0,0.5}.
#' @param calibrated.border Colour of the border of the unmodelled non-14C ages. Defaults to nothing, NA.
#' @param calBP.col Colour of the inside of the unmodelled non-14C ages. Defaults to semi-transparent light grey, \code{rgb(0,0,0,0.5}.
#' @param calBP.border Colour of the border of the unmodelled, calibrated ages. Defaults to nothing, NA.
#' @param modelled.col Colour of the inside of the modelled ages. Defaults to semi-transparent dark grey, \code{rgb(0,0,0,0.5}.
#' @param modelled.border Colour of the border of the modelled ages. Defaults to semi-transparent dark grey, \code{rgb(0,0,0,0.5}.
#' @param range.col Colour of the hpd ranges. Defaults to \code{"black"}.
#' @param block.col Colour of the field indicating unordered dates within a 'block'. Defaults to light-blue, \code{rgb(0,0,1,.05)}.
#' @param gap.col Colour of the diagonal lines and parameters of any gaps.
#' @param gap.pos Plotting position of the gap information.
#' @param simulation Whether or not the data are part of a simulated stratigraphy.
#' @param simulation.col If the data are part of a simulated stratigraphy, the 'true' ages can also be plotted.
#' @param pos.lim Limit of the main y-axis.
#' @param age.lim Limit of the main x-axis.
#' @param mgp Axis text margins (where should titles, labels and tick marks be plotted). Defaults to \code{mgp=c(1.7, .7, .0)}.
#' @param mar.top Margins around the top panel. Defaults to \code{mar.top=c(3,3,1,1)}.
#' @param mar.bottom Margins around the bottom panel. Defaults to \code{mar.bottom=c(3,3,0.5,1)}.
#' @param heights Relative heights of the two panels in the plot. Defaults to 0.3 for the top and 0.7 for the bottom panel.
#' @param iterations.warning Whether or not to plot a warning if there are < 3000 iterations, too few for a reliable MCMC run.
#' @param min.its The minimum amount of iterations, below which a warning is printed (if \code{iterations.warning=TRUE}). Defaults to 1,000. 
#' @param warning.loc Location of the warning
#' @param warning.col Colour of the warning - defaults to red.
#' @return A plot with two panels showing the MCMC run and the calibrated and modelled ages.
#' @author Maarten Blaauw
#' @export
draw.strat <- function(name="mystrat", set=get('info'), structure=set$struc, y.scale="positions", strat.dir="strats", cc.dir=c(), sep=",", postbomb=FALSE, calibrated.ex=c(), calibrated.mirror=FALSE, calibrated.up=TRUE, modelled.ex=c(), modelled.mirror=FALSE, modelled.up=FALSE, BCAD=FALSE, threshold=0.001, xtop.lab=c(), ytop.lab=c(), xbottom.lab=c(), ybottom.lab="position", calibrated.col=rgb(0, 0, 0, 0.2), calibrated.border=NA, calBP.col=rgb(0, 0, 0, 0.2), calBP.border=NA, modelled.col=rgb(0,0,0,0.5), modelled.border=rgb(0,0,0,0.5), range.col="black", block.col=rgb(0,0,1,.05), gap.col="blue", gap.pos=1, simulation=FALSE, simulation.col=grey(0.5), pos.lim=c(), age.lim=c(), mgp=c(2, 0.7, 0), mar.top=c(3,3,1,1), mar.bottom=c(3,3,0.5,1), heights=c(0.3, 0.7), iterations.warning=TRUE, min.its=1e3, warning.loc="bottomleft", warning.col="red") {
  layout(matrix(1:2, ncol=1), heights=heights)
  oldpar <- par(no.readonly = TRUE)
  #on.exit(par(oldpar))
  newpar <- par(mar=mar.top, mgp=mgp)

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
  ## add ability to reload info?

  if(length(xtop.lab) == 0)
    xtop.lab <- "iterations"
  if(length(ytop.lab) == 0)
    ytop.lab <- "energy"
  plot(-1*energy, type="l", xlab=xtop.lab, ylab=ytop.lab)
  if(iterations.warning)
    if(nrow(set$output) < min.its) # then MCMC run likely not long enough to be reliable
      legend(warning.loc, legend="short MCMC run - needs more iterations", bty="n", cex=.8, text.col=warning.col)

  dates <- which(dets[,5] %in% 0:4)
  ages.positions <- dets[which(dets[,5] %in% c(0:4,10)),4]
  if(y.scale[1] == "dates")
   ages.positions <- 1:length(ages.positions)

  if(y.scale[1] == "positions")
    y.scale <- dets[dates,4]
  if(y.scale[1] == "dates")
    y.scale <- 1:length(dates)
  if(length(pos.lim) == 0)
    pos.lim <- extendrange(y.scale, f=0.1)

  # calculate highest prob of all dates, prior to drawing
  calib.mx <- 0
  for(i in 1:length(dates)) {
    tmp <- caldist(dets[dates[i],2], dets[dates[i],3], cc = dets[dates[i],5], cc.dir = cc.dir)
    calib.mx <- max(calib.mx, tmp[, 2])
  }
  calib.ex <- .05 / (calib.mx * nrow(dets)^1.2)
  if(length(calibrated.ex) == 0) 
    calibrated.ex <- calib.ex  
  
  mod.mx <- 0
  for(i in 1:ncol(ages))
    mod.mx <- max(mod.mx, density(ages[,i])$y)
  if(length(modelled.ex) == 0)    	
    modelled.ex <- .02 / (mod.mx * nrow(dets)^1.2)

  par(mar=mar.bottom)
  Dates <- draw.dates(dets[dates,2], dets[dates,3], y.scale, dets[dates,5], cc.dir=cc.dir, postbomb=postbomb, ex=calibrated.ex, mirror=calibrated.mirror, up=calibrated.up, col=calibrated.col, border=calibrated.border, cal.col=calBP.col, cal.border=calBP.border, BCAD=BCAD, draw.hpd=FALSE, threshold=threshold, normalise=FALSE, cal.lab=xbottom.lab, d.lab=ybottom.lab, age.lim=age.lim, d.lim=pos.lim)

  if(set$struc$has.blocks) {
    above <- set$struc$above.block
    below <- set$struc$below.block
    y.above <- c(); y.below <- c()
    for(i in 1:length(above)) {
      y.above[i] <- mean(y.scale[above[i]+c(0,1)])
      y.below[i] <- mean(y.scale[below[i]-c(0,1)])
    }
    rect(min(.1*dets[,2],-9e6), y.above, max(2*dets[,2], 9e6), y.below, col=block.col, border=NA)
  }
 
  maxdens <- 0
  for(i in 1:ncol(ages))
    maxdens <- max(maxdens, density(ages[,i])$y)
  if(BCAD)
    ages <- 1950 - ages

  hpds <- list()
  for(i in 1:ncol(ages)) {
    age <- density(ages[,i])
    age$y <- age$y / maxdens
    hpds[[i]] <- rice::hpd(cbind(age$x, age$y))

    if(modelled.mirror)
      pol <- cbind(c(age$x, rev(age$x)), ages.positions[i]-modelled.ex*c(age$y, -rev(age$y))) else
        if(modelled.up)
          pol <- cbind(c(min(age$x), age$x, max(age$x)), ages.positions[i]+modelled.ex*c(0, age$y, 0)) else
            pol <- cbind(c(min(age$x), age$x, max(age$x)), ages.positions[i]-modelled.ex*c(0, age$y, 0))
    polygon(pol, col=modelled.col, border=modelled.border)
    if(simulation)
      segments(ifelse(BCAD, 1950-dets[i,1], dets[i,1]), i, ifelse(BCAD, 1950-dets[i,1], dets[i,1]), ncol(ages)+10, lty=2, col=simulation.col) # then don't be afraid to show the truth!
    if(!is.na(range.col))
      segments(hpds[[i]][,1], ages.positions[i], hpds[[i]][,2], ages.positions[i], col=range.col, lwd=2)
  }

  if(structure$has.gaps) {
    gaps <- structure$gaps
    for(i in 1:nrow(gaps)) {
      # above <- max(which(structure$p < gaps[i,4])) # commented 29 April 2024
      above <- max(which(dets[,4] < gaps[i,4])) # p is the position index, whereas gaps[i,4] is the position "number"
      side <- ifelse(gap.pos==1, 1, 2)
      if(gap.pos==1) {
        age.above <- max(hpds[[above]][,-3])
        age.below <- max(hpds[[above+1]][,-3])
        } else {
          age.above <- min(hpds[[above]][,-3])
          age.below <- min(hpds[[above+1]][,-3])
        }
      segments(age.above, dets[above,4],
        age.above+gaps[i,2], dets[above+1,4], col=gap.col)
      adj <- ifelse(gap.pos==1, c(1,.5), c(0,.5))
      xpos <- ifelse(gap.pos==1, age.below, age.above)
      if(gaps[i,5] == 11)
        txt <- paste0("Ex(", gaps[i,2], ")") else
          if(gaps[i,5] == 12)
            txt <- paste0("N(", gaps[i,2], ",", gaps[i,3], ")") else
              if(gaps[i,5] == 13)
                txt <- paste0("Ga(", gaps[i,2], ",", gaps[i,3], ")")
      text(xpos, mean(c(dets[above,4], dets[above+1,4])), txt, cex=0.5, col=gap.col, adj=adj)
    }
  }
  
  #set$hpds <- hpds
  #set$dates <- dates 
  #assign_to_global("info", set)
  Dates$hpds <- hpds 
  invisible(Dates)
}



#' @name draw.rings
#' @title plot the dates and model of a wiggle-match dated tree
#' @description A plot with two panels. The top panel shows the calibrated distributions (in blue) and the wiggle-match age-modelled age estimates for each dated ring (grey). The bottom panel shows the fit of the uncalibrated radiocarbon dates (steelblue dots and lab error bars) to the calibration curve (green), as well as the age estimate of the oldest/starting ring (grey) and its hpd range (black).
#' @param name The name of the tree. The .csv file should be saved under a folder named exactly the same as \code{name}, and the folder should live under the \code{tree.dir} folder. Default names include Ulandryk4 and mytree.
#' @param tree.dir The directory where the folders of the individual trees live. Defaults to \code{treedir="trees"}.
#' @param sep Separator for the fields in the .csv file. Defaults to a comma.
#' @param normal Calculations can be done assuming that the measurements are normally distributed. By default this is set to FALSE and a student-t distribution is used (Christen and Perez 2009).
#' @param dat If \code{plot.rings} is called from within \code{rings()}, both dat and out are provided. If an existing run has to be plotted again, dat and out are read from the files in the folder named \code{name}.
#' @param out If \code{plot.rings} is called from within \code{rings()}, both dat and out are provided. If an existing run has to be plotted again, dat and out are read from the files in the folder named \code{name}.
#' @param cc Calibration curve to be used. Could be 1 (IntCal20; default), 2 (Marine20), 3 (SHCal20) or 4 (custom curve).
#' @param postbomb Negative C-14 ages should be calibrated using a postbomb curve. This could be 1 (northern-hemisphere region 1), 2 (NH region 2), 3 (NH region 3), 4 (southern hemisphere regions 1-2), or 5 (SH region 3).
#' @param BCAD The calendar scale of graphs and age output-files is in \code{cal BP} by default, but can be changed to BC/AD using \code{BCAD=TRUE}.
#' @param t.a First parameter for the student-t distribution (defaults to 3; higher numbers make the distribution approximate the normal distribution more).
#' @param t.b Second parameter for the student-t distribution (defaults to 4; higher numbers make the distribution approximate the normal distribution more).
#' @param x.lim Limits of the x-axes. Calculated automatically by default.
#' @param x1.axis Whether or not to plot the upper x-axis (slightly redundant since the bottom axis shows the values already). Defaults to TRUE.
#' @param x1.labels Whether or not to plot the labels of the upper x-axis. (slightly redundant since the bottom axis shows the values already). Defaults to FALSE.
#' @param x1.lab The labels for the calendar axis of the upper panel. Defaults to empty.
#' @param rev.x Whether or not to reverse the x-axis. Defaults to \code{FALSE}.
#' @param y1.lab The labels for the y-axis. Defaults to \code{"rings"}.
#' @param y1.lim Limits of the top y-axis. Calculated automatically if left empty.
#' @param y2.lim Limits of the bottom y-axis. Calculated automatically if left empty.
#' @param x2.lab The labels for the bottom calendar axis (default \code{age.lab="cal BP"} or \code{"BC/AD"} if \code{BCAD=TRUE}).
#' @param y2.lab The labels for the bottom y-axis. Defaults to 14C BP with superscript 14, so \code{expression(""^14*C~BP)}.
#' @param ex Exaggeration of the heights of the calibrated distributions. Defaults to 0.05 so there is plenty space to plot many distributions.
#' @param plot.cc Whether or not to plot a panel with the calibration curve.
#' @param plot.dists Whether or not to plot a panel with the distributions.
#' @param mar.1 Margins around the first/top plot.
#' @param mar.2 Margins around the first/bottom plot.
#' @param mgp Axis text margins (where should titles, labels and tick marks be plotted). Defaults to \code{mgp=c(1.7, .7, .0)}.
#' @param dist.res Resolution of the plot of the distribution. The default is \code{500}, which results in smooth plots.
#' @param date.col Colour of the uncalibrated dates when plotted on top of the calibration curve. Defaults to \code{"steelblue"}.
#' @param cc.col Colour of the calibration curve. Defaults to semi-transparent darkgreen, \code{cc.col=rgb(0 0.5, 0, 0.5)}.
#' @param dist.col Colour of the age-modelled distribution. Defaults to semi-transparent grey, \code{dist.col=rgb(0,0,0,0.5)}.
#' @param calib.col Colour of the calibrated distributions. Defaults to semi-transparent blue, \code{dist.col=rgb(0,0,1,0.5)}.
#' @param range.col Colour of the hpd ranges. Defaults to \code{"black"}.
#' @param set.layout By default, the layout of the two plots is set automatically (2 plots in one column).
#' @return A plot with the calibrated distributions of the individual dates and the wiggle-match distributions (top), and the dates on the calibration curve together with the age distribution for the earliest ring, 0.
#' @examples
#'   treedir <- tempdir()
#'   rings("Ulandryk", tree.dir=treedir, draw=FALSE)
#'   draw.rings("Ulandryk", tree.dir=treedir)
#' @author Maarten Blaauw, J. Andres Christen
#' @export
draw.rings <- function(name="mytree", tree.dir="trees", sep=",", normal=TRUE, dat=c(), out=c(), cc=1, postbomb=FALSE, BCAD=FALSE, t.a=3, t.b=4, x.lim=c(), x1.axis=TRUE, x1.labels=FALSE, x1.lab=c(), rev.x=FALSE, y1.lab=c(), y1.lim=c(), y2.lim=c(), x2.lab=c(), y2.lab=c(), ex=0.05, plot.cc=TRUE, plot.dists=TRUE, mar.1=c(1,3,1,1), mar.2=c(3,3,0,1), mgp=c(1.7, .7, 0), dist.res=500, date.col="steelblue", cc.col=rgb(0, 0.5, 0, 0.5), dist.col=rgb(0,0,0,0.5), calib.col=rgb(0,0,1,0.25), range.col="black", set.layout=TRUE) {

  # usually draw.rings is called from within the function rings. But runs could also be reloaded, and that would happen if dat or out are not provided. Then this is read from the file.
  if(length(dat) == 0)
    dat <- read.table(file.path(tree.dir, name, paste0(name, ".csv")), header=TRUE, sep=sep)
  if(length(out) == 0)
    out <- read.table(file.path(tree.dir, name, paste0(name, "_probs.txt")), header=TRUE, sep="\t")
  if(length(cc) == 1)
    cc <- ccurve(1, postbomb)

  if(length(y2.lim) > 0)
    C14.lim <- y2.lim else
      C14.lim <- extendrange(c(dat[,2]-dat[,3], dat[,2]+dat[,3]))
  if(length(x.lim) == 0)
    cal.lim <- range(out[,1], out[,1] - (max(dat[,4])-min(dat[,4]))) else
      ifelse(BCAD, cal.lim <- 1950 - x.lim, cal.lim <- x.lim)

  age.min <- min(out[,1])
  age.max <- age.min + (max(dat[,4]) - min(dat[,4])) # how many rings
  age.lim <- extendrange(c(age.min, age.max), f=.1) # add 10% to each end
  best.age <- out[which(out[,2] == max(out[,2])),1]
  exg <- 100*ex
  ex <- ex * (max(C14.lim)-min(C14.lim)) / max(out[,2])
  agepol <- cbind(c(min(out[,1]), out[,1], max(out[,1])),
    min(C14.lim)+c(0, 2.5*ex*out[,2], 0))

  o <- order(out[,2], decreasing=TRUE)
  rng <- cbind(out[o,1], cumsum(out[o,2]))
  rng <- range(rng[which(rng[,2] <= .95),1])
  if(BCAD)
    message("best age: ", 1950-best.age, " BC/AD, 95% range: ", 1950-rng[2], " to ", 1950-rng[1], " BC/AD, ", rng[2]-rng[1], " yr") else
      message("best age: ", best.age, " cal BP, 95% range: ", rng[2], " to ", rng[1], " cal BP, ", rng[2]-rng[1], " yr")

  ccpol <- cc#[cc.min:cc.max,]
  ccpol <- cbind(c(ccpol[,1], rev(ccpol[,1])), c(ccpol[,2]-ccpol[,3], rev(ccpol[,2]+ccpol[,3])))

  if(set.layout)
    if(plot.cc && plot.dists)
      layout(matrix(1:2, ncol=1)) else
      if(plot.cc || plot.dists)
        layout(1)

  oldpar <- par(no.readonly=TRUE)
  on.exit(par(oldpar))
  par(mar=mar.1, mgp=mgp)
    xaxt <- ifelse(BCAD, "n", "s")
    if(length(x1.lab) == 0)
      x1.lab <- ""
    if(length(y1.lab) == 0)
      y1.lab <- "rings"
    if(length(y1.lim) == 0)
      y1.lim <- rev(extendrange(dat[,4], f=c(0, .1)))
      #instead, find the last date and the height of the probdist
    if(length(x2.lab) == 0)
      x2.lab <- ifelse(BCAD, "BC/AD", "cal BP")
    if(length(y2.lab) == 0)
      y2.lab <- expression(""^14*C~BP)
    if(rev.x)
      cal.lim <- cal.lim[2:1]

    if(plot.dists) {
      plot(0, type="n", xlim=cal.lim, ylim=y1.lim, xlab=x1.lab, ylab=y1.lab, bty="l", xaxt="n")
      if(x1.axis)
        if(x1.labels) {
          if(BCAD)
            axis(1, pretty(cal.lim), 1950 - pretty(cal.lim)) else
              axis(1, pretty(cal.lim))
        } else axis(1, labels=FALSE)
      yrs <- seq(min(cal.lim), max(cal.lim), length=dist.res)
      for(i in 1:nrow(dat)) {
        prob <- l.calib(yrs, dat[i,2], dat[i,3], cc, normal=normal, t.a=t.a, t.b=t.b)
        calib <- cbind(yrs, prob/max(prob))
        calib <- calib[-calib[,2]<1e-2,]
        calib.pol <- cbind(c(min(calib[,1]), calib[,1], max(calib[,1])),
          dat[i,4]-c(0, exg*calib[,2], 0))
        polygon(calib.pol, col=calib.col, border=calib.col)
		postpol <- cbind(out[,1] - dat[i,4], dat[i,4]-exg*out[,2]/max(out[,2]))
        polygon(postpol, col=dist.col, border=dist.col)
      }
    }

    if(plot.cc) {
      if(plot.dists) # then plotting both
        par(mar=mar.2)
      plot(best.age-dat[,4], dat[,2], pch=20, col=date.col, xlab=x2.lab, ylab=y2.lab, xlim=cal.lim, ylim=C14.lim, bty="l", yaxs="i", xaxt=xaxt)
      if(BCAD)
        axis(1, pretty(cal.lim), 1950 - pretty(cal.lim))
      polygon(ccpol, col=cc.col, border=cc.col) # calibration curve
      segments(best.age-dat[,4], dat[,2]-dat[,3], best.age-dat[,4], dat[,2]+dat[,3], col=date.col)
      polygon(agepol, col=dist.col, border=dist.col)
      segments(rng[1], min(C14.lim), rng[2], min(C14.lim), lwd=5, col=range.col)
    }
}
