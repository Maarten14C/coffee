

# renamed to 'rings' because tree() is already taken by cli R

#' @name rings
#' @title wiggle-match C-14 dating of a tree
#' @description Produce a Bayesian wiggle-match date of a tree dated with multiple C-14 dates at exactly known spacing (e.g., every 10 tree-ring years).
#' @details The calculations are based on Bwigg (Christen and Litton 1995; Christen 2003). 
#' In OxCal, this is called a D_Sequence (Bronk Ramsey et al. 2001).
#' 
#' Since only one parameter has to be estimated, a MCMC approach is not necessary nor recommended, and results are calculated analytically.
#'
#' Files for tree wiggle-matching should contain the following columns: lab ID, C-14 age, error, ring, cc. Rings are counted from the outward ring (0 year old) inwards. cc should either be 1 (IntCal20; northern hemisphere #' terrestrial, 2 (Marine20, though I've never heard of marine trees), 3 (SHCal20; southern hemisphere) or 4 (custom curve). The tree files should be in plain-text and fields separated by commas.
#'
#' Default tree files include Ulandryk4 (Christen 2003) and a simulated tree (see \code{sim.tree()}).
#'
#' By default, the data are calibrated assuming a student-t distribution, which has wider tails than the normal distribution and deals well with outliers ().
#' @param name name of the tree. The .csv file should be saved under a folder named exactly the same as \code{name}, and the folder should live under the \code{treedir} folder. Default names include Ulandryk4 and mytree.
#' @param ring The ring for which the age estimate should be calculated. Defaults to the youngest ring, 0.
#' @param tree.dir The directory where the folders of the individual trees live. Defaults to \code{treedir="trees"}.
#' @param sep Separator for the fields in the .csv file. Defaults to a comma.
#' @param normal Calculations can be done assuming that the measurements are normally distributed. By default this is set to FALSE and a student-t distribution is used (Christen and Perez 2009)
#' @param delta.R The ages can be modelled to have an offset. The mean is 0 by default.
#' @param delta.STD The error of the offset. Set to 0 by default.
#' @param t.a First parameter for the student-t distribution (defaults to 3; higher numbers make the distribution approximate the normal distribution more).
#' @param t.b Second parameter for the student-t distribution (defaults to 4; higher numbers make the distribution approximate the normal distribution more).
#' @param ask Whether or not to ask if new folders should be written (if required)
#' @param age.steps Steps in years for the calculations. Defaults to 1, every year.
#' @param cutoff Value below which probabilities are no longer taken into account. Defaults to 0.000001.
#' @param cc Calibration curve to be used. Could be 1 (IntCal20; default), 2 (Marine20), 3 (SHCal20) or 4 (custom curve).
#' @param postbomb Negative C-14 ages should be calibrated using a postbomb curve. This could be 1 (northern-hemisphere region 1), 2 (NH region 2), 3 (NH region 3), 4 (southern hemisphere regions 1-2), or 5 (SH region 3).
#' @param BCAD The calendar scale of graphs and age output-files is in \code{cal BP} by default, but can be changed to BC/AD using \code{BCAD=TRUE}.
#' @param times The range of years to be calculated, as multiples of the uncertainties of the data points. E.g. if the lab error of the oldest date is 20 years, and times is set to 5, the calculation range will be extended by 20*5 years.
#' @param talk Whether or not to provide feedback on folders written into.
#' @param draw Whether or not to draw plots.
#' @param ... Options for the plot. See \code{draw.rings}.
#' @return the probabilities for the relevant calendar years.
#' @examples
#'   rings("Ulandryk4", tree.dir=tempdir())
#'   rings("mytree", tree.dir=tempdir())
#' @author Maarten Blaauw, J. Andres Christen
#' @references
#' Bronk Ramsey C, van der Plicht J, Weninger B, 2001. 'Wiggle matching' radiocarbon dates. Radiocarbon 43, 381–389.
#'
#' Christen JA, Litton CD, 1995. A Bayesian approach to wiggle-matching. Journal of Archaeological Science 22, 719-725 \doi{10.1016/0305-4403(95)90002-0}
#'
#' Christen JA, 2003. Bwigg: An Internet facility for Bayesian radiocarbon wiggle-matching. Internet Archaeology 13. \doi{10.11141/ia.13.2}
#'
#' Christen JA, Perez S, 2009. A new robust statistical model for radiocarbon data. Radiocarbon 51, 1047-1059
#' @export
rings <- function(name="mytree", ring=0, tree.dir="trees", sep=",", normal=FALSE, delta.R=0, delta.STD=0, t.a=3, t.b=4, ask=TRUE, age.steps=1, cutoff=1e-6, cc=1, postbomb=FALSE, BCAD=FALSE, times=3, talk=TRUE, draw=TRUE, ...) {

  tree.dir <- assign_dir(tree.dir, name, "tree.dir", ask, talk)
  if(name %in% c("mytree", "Ulandryk4")) {
    fileCopy <- system.file(paste0("extdata/", name, ".csv"), package="coffee")
    file.copy(fileCopy, tree.dir, recursive = TRUE, overwrite=FALSE)
  }
  dat <- read.table(file.path(tree.dir, paste0(name, ".csv")), header=TRUE, sep=",")

  # file sanity checks; ring is tree ring number
  OK <- TRUE
  if(min(diff(dat[,4])) < 0) OK <- FALSE # we need increasing positions
  if(min(dat[,5]) < 1 || max(dat[,5]) > 4) OK <- FALSE # only use C14 data
  if(ncol(dat) < 4) OK <- FALSE # only 4 or 5 columns
  if(!OK)
    stop("Unexpected values in tree file. Please check", call.=FALSE)

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

  # read the calibration curve, select the relevant section and interpolate to (by default) years
  # MCMC not required or recommendable here; all calendar years in the relevant section are calculated
  cc <- ccurve(cc, postbomb)
  cc.min <- max(1, which(cc[,2] <= (max(dat[,2]+times*dat[,3]))))
  cc.max <- min(nrow(cc), which(cc[,2] >= min(dat[,2]-times*dat[,3])))
  yrseq <- seq(cc[cc.max,1], cc[cc.min,1], by=age.steps)
  mu <- approx(cc[,1], cc[,2], yrseq)$y
  sigma <- approx(cc[,1], cc[,3], yrseq)$y
  ccc <- cbind(yrseq, mu, sigma)
  #yrseq <- cc[cc.min:cc.max,1]

  probs <- rep(0, length(yrseq))
  for(i in 1:length(yrseq)) {
    probs1 <- l.calib(yrseq[i]+dat[,4], dat[,2], dat[,3], cc=ccc, normal=normal, t.a=t.a, t.b=t.b)
    probs[i] <- prod(probs1) # this OK? Or better the sum of the logs?
  }
  probs <- probs / sum(probs) # normalise
  OK <- which(probs > cutoff)
  out <- cbind(yrseq[OK], probs[OK])
  colnames(out) <- c("yr", "prob")
  write.table(out, file.path(tree.dir, paste0(name, "_probs.txt")), sep="\t", row.names=FALSE, quote=FALSE)
  if(draw)
    draw.rings(name=name, dat=dat, out=out, cc=cc, BCAD=BCAD, normal=normal, t.a=t.a, t.b=t.b, ...)
  invisible(list(dat, out))
}



#' @name draw.rings
#' @title plot the dates and model of a wiggle-match dated tree
#' @description A plot with two panels. The top panel shows the calibrated distributions (in blue) and the wiggle-match age-modelled age estimates for each dated ring (grey). The bottom panel shows the fit of the uncalibrated radiocarbon dates (steelblue dots and lab error bars) to the calibration curve (green), as well as the age estimate of the youngest ring (grey) and its hpd range (black).
#' @param name The name of the tree. The .csv file should be saved under a folder named exactly the same as \code{name}, and the folder should live under the \code{tree.dir} folder. Default names include Ulandryk4 and mytree.
#' @param tree.dir The directory where the folders of the individual trees live. Defaults to \code{treedir="trees"}.
#' @param sep Separator for the fields in the .csv file. Defaults to a comma.
#' @param normal Calculations can be done assuming that the measurements are normally distributed. By default this is set to FALSE and a student-t distribution is used (Christen and Perez 2009)
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
#' @param x2.lab The labels for the bottom calendar axis (default \code{age.lab="cal BP"} or \code{"BC/AD"} if \code{BCAD=TRUE}).
#' @param y2.lab The labels for the bottom y-axis. Defaults to 14C BP with superscript 14, so \code{expression(""^14*C~BP)}.
#' @param ex Exaggeration of the heights of the calibrated distributions. Defaults to 0.2 so there is plenty space to plot many distributions.
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
#' @return A plot with the calibrated distributions of the individual dates and the wiggle-match distributions (top), and the dates on the calibration curve together with the age distribution for ring 0.
#' @examples
#'   treedir <- tempdir()
#'   rings("Ulandryk4", tree.dir=treedir, draw=FALSE)
#'   draw.rings("Ulandryk4", tree.dir=treedir)
#' @author Maarten Blaauw, J. Andres Christen
#' @export
draw.rings <- function(name="mytree", tree.dir="trees", sep=",", normal=TRUE, dat=c(), out=c(), cc=1, postbomb=FALSE, BCAD=FALSE, t.a=3, t.b=4, x.lim=c(), x1.axis=TRUE, x1.labels=FALSE, x1.lab=c(), rev.x=FALSE, y1.lab=c(), x2.lab=c(), y2.lab=c(), ex=0.2, plot.cc=TRUE, plot.dists=TRUE, mar.1=c(1,3,1,1), mar.2=c(3,3,0,1), mgp=c(1.7, .7, 0), dist.res=500, date.col="steelblue", cc.col=rgb(0, 0.5, 0, 0.5), dist.col=rgb(0,0,0,0.5), calib.col=rgb(0,0,1,0.25), range.col="black", set.layout=TRUE) {

  # usually draw.rings is called from within the function rings. But runs could also be reloaded, and that would happen if dat or out are not provided. Then this is read from the file.
  if(length(dat) == 0)
    dat <- read.table(file.path(tree.dir, name, paste0(name, ".csv")), header=TRUE, sep=sep)
  if(length(out) == 0)
    out <- read.table(file.path(tree.dir, name, paste0(name, "_probs.txt")), header=TRUE, sep="\t")
  if(length(cc) == 1)
    cc <- ccurve(1, postbomb)

  C14.lim <- extendrange(c(dat[,2]-dat[,3], dat[,2]+dat[,3]))
  if(length(x.lim) == 0)
    cal.lim <- range(out[,1], out[,1] + (max(dat[,4])-min(dat[,4]))) else
      ifelse(BCAD, cal.lim <- 1950 - x.lim, cal.lim <- x.lim)

  age.min <- min(out[,1])
  age.max <- age.min + (max(dat[,4]) - min(dat[,4])) # how many rings
  age.lim <- extendrange(c(age.min, age.max), f=.1) # add 10% to each end
  best.age <- out[which(out[,2] == max(out[,2])),1]
  ex <- ex * (max(C14.lim)-min(C14.lim)) / max(out[,2])
  agepol <- cbind(c(min(out[,1]), out[,1], max(out[,1])), 
    min(C14.lim)+c(0, ex*out[,2], 0))

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
    if(length(x2.lab) == 0)
      x2.lab <- ifelse(BCAD, "BC/AD", "cal BP")
    if(length(y2.lab) == 0)
      y2.lab <- expression(""^14*C~BP)
    if(rev.x)
      cal.lim <- cal.lim[2:1]

    if(plot.dists) {
      plot(0, type="n", xlim=cal.lim, ylim=range(dat[,4]), xlab=x1.lab, ylab=y1.lab, bty="l", xaxt="n")
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
          dat[i,4]+c(0, 20*calib[,2], 0))
        polygon(calib.pol, col=calib.col, border=calib.col)
        postpol <- cbind(dat[i,4]+out[,1], dat[i,4]+10*out[,2]/max(out[,2]))
        polygon(postpol, col=dist.col, border=dist.col)
      }
    }

    if(plot.cc) {
      if(plot.dists) # then plotting both
        par(mar=mar.2)
      plot(best.age+dat[,4], dat[,2], pch=20, col=date.col, xlab=x2.lab, ylab=y2.lab, xlim=cal.lim, ylim=C14.lim, bty="l", yaxs="i", xaxt=xaxt)
      if(BCAD)
        axis(1, pretty(cal.lim), 1950 - pretty(cal.lim))  
      polygon(ccpol, col=cc.col, border=cc.col) # calibration curve
      segments(best.age+dat[,4], dat[,2]-dat[,3], best.age+dat[,4], dat[,2]+dat[,3], col=date.col)
      polygon(agepol, col=dist.col, border=dist.col)
      segments(rng[1], min(C14.lim), rng[2], min(C14.lim), lwd=5, col=range.col)
    }
}


