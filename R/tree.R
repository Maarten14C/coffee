# renamed to 'rings' because tree() is already taken by cli R

#' @name rings
#' @title wiggle-match C-14 dating of a tree
#' @description Produce a Bayesian wiggle-match date of a tree dated with multiple C-14 dates at exactly known spacing (e.g., every 10 tree-ring years).
#' @details The calculations are based on Bwigg (Christen and Litton 1995; Christen 2003). 
#' In OxCal, this is called a D_Sequence (Bronk Ramsey et al. 2001).
#' 
#' Since only one parameter has to be estimated (the age of the earliest, innermost ring),
#' a MCMC approach is not necessary nor recommended, and results are calculated analytically.
#'
#' The tree files should be in plain-text and fields separated by commas, and the file's extension should be ".csv".
#' The files should start with a line contain the following headers: "lab ID", "C-14 age", "error", "ring", "cc", separated by commas. Then each row should have the corresponding values, also separated by commas.
#' Rings are counted from the inner ring (0 year old) outwards, so, forward in time. The file should start with the youngest rings, then work downward until reaching the oldest, bottommost dated rings. 
#' cc should either be 1 (IntCal20; northern hemisphere terrestrial, 2 (Marine20, though we've never heard of marine trees), 3 (SHCal20; southern hemisphere) or 4 (custom curve). 
#'
#' The default tree is called Ulandryk (Kuzman et al. 2004). As an alternative, a tree can be simulated (see \code{sim.tree()}).
#'
#' By default, the data are calibrated assuming a student-t distribution, which has wider tails than the normal distribution and deals well with scatter and outliers.
#' @param name Name of the tree. The .csv file should be saved under a folder named exactly the same as \code{name}, and the folder should live under the \code{treedir} folder. The default is Ulandryk.
#' @param tree.dir The directory where the folders of the individual trees live. Defaults to \code{tree.dir="trees"}.
#' @param sep Separator for the fields in the .csv file. Defaults to a comma.
#' @param normal Calculations can be done assuming that the measurements are normally distributed. By default this is set to FALSE and a student-t distribution is used (Christen and Perez 2009)
#' @param delta.R The ages can be modelled to have an offset. The mean is 0 by default.
#' @param delta.STD The error of the offset. Set to 0 by default.
#' @param t.a First parameter for the student-t distribution (defaults to 3; higher numbers make the distribution approximate the normal distribution more).
#' @param t.b Second parameter for the student-t distribution (defaults to 4; higher numbers make the distribution approximate the normal distribution more).
#' @param ask Whether or not to ask if new folders should be written (if required)
#' @param age.steps Steps in years for the calculations. Defaults to 1, every year.
#' @param cutoff Value below which probabilities are no longer taken into account. Defaults to 0.000001.
#' @param prob After the run, a fit of the model with the dates is calculated, as the ratio of model iterations that fit the hpd ranges of the dates. Defaults to the \code{prob=0.95} hpd ranges.
#' @param cc Calibration curve to be used, for glueing to a postbomb curve. Could be 1 (IntCal20; default), 2 (Marine20), 3 (SHCal20) or 4 (custom curve). Normally not used, except in the case where there are postbomb dates (requiring the 'gluing' of pre- and postbomb curves).
#' @param postbomb Negative C-14 ages should be calibrated using a postbomb curve. This could be 1 (northern-hemisphere region 1), 2 (NH region 2), 3 (NH region 3), 4 (southern hemisphere regions 1-2), or 5 (SH region 3).
#' @param BCAD The calendar scale of graphs and age output-files is in \code{cal BP} by default, but can be changed to BC/AD using \code{BCAD=TRUE}.
#' @param times The range of years to be calculated, as multiples of the uncertainties of the data points. E.g. if the lab error of the oldest date is 20 years, and times is set to 5, the calculation range will be extended by 20*5 years.
#' @param talk Whether or not to provide feedback on folders written into.
#' @param draw Whether or not to draw plots.
#' @param ... Options for the plot. See \code{draw.rings}.
#' @return the probabilities for the relevant calendar years.
#' @examples
#'   rings("Ulandryk", tree.dir=tempdir())
#' @author Maarten Blaauw, J. Andres Christen
#' @references
#' Bronk Ramsey C, van der Plicht J, Weninger B, 2001. 'Wiggle matching' radiocarbon dates. Radiocarbon 43, 381–389.
#'
#' Christen JA, Litton CD, 1995. A Bayesian approach to wiggle-matching. Journal of Archaeological Science 22, 719-725 \doi{10.1016/0305-4403(95)90002-0}
#'
#' Christen JA, 2003. Bwigg: An Internet facility for Bayesian radiocarbon wiggle-matching. Internet Archaeology 13. \doi{10.11141/ia.13.2}
#'
#' Christen JA, Perez S, 2009. A new robust statistical model for radiocarbon data. Radiocarbon 51, 1047-1059
#'
#' Kuzmin Y, Slusarenko I, Hajdas I, Bonani G, Christen JA. 2004. The comparison of 14C wiggle-matching results for the ‘floating’ tree-ring chronology of the Ulandryk-4 Burial Ground (Altai Mountains, Siberia). Radiocarbon 46, 943–948.
#'
#' @export
rings <- function(name="mytree", tree.dir="trees", sep=",", normal=FALSE, delta.R=0, delta.STD=0, t.a=3, t.b=4, ask=TRUE, age.steps=1, cutoff=1e-6, prob=.95, cc=1, postbomb=FALSE, BCAD=FALSE, times=3, talk=TRUE, draw=TRUE, ...) {

  tree.dir <- assign_dir(tree.dir, name, "tree.dir", ask=FALSE, talk)
  if(name %in% "Ulandryk") {
    fileCopy <- system.file("extdata/Ulandryk.csv", package="coffee")
    file.copy(fileCopy, file.path(tree.dir), recursive = TRUE, overwrite=FALSE)
  }
  dat <- read.table(file.path(tree.dir, paste0(name, ".csv")), header=TRUE, sep=",")

  # file sanity checks; ring is tree ring number
  OK <- TRUE
  if(max(diff(dat[,4])) > 0) OK <- FALSE # we need decreasing ring positions
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
  if(ncol(dat) > 4)
    cc <- unique(dat[which(dat[,5] %in% 1:4),5])[1] # only use 1 entry
  cc <- ccurve(cc, postbomb)
  
  cc.min <- max(1, which(cc[,2] <= (max(dat[,2]+times*dat[,3]))))
  cc.max <- min(nrow(cc), which(cc[,2] >= min(dat[,2]-times*dat[,3])))
  yrseq <- seq(cc[cc.max,1], cc[cc.min,1], by=age.steps)
  mu <- approx(cc[,1], cc[,2], yrseq)$y
  sigma <- approx(cc[,1], cc[,3], yrseq)$y
  ccc <- cbind(yrseq, mu, sigma)

  probs <- rep(0, length(yrseq))
  for(i in 1:length(yrseq)) {
    probs1 <- l.calib(yrseq[i]-dat[,4], dat[,2], dat[,3], cc=ccc, normal=normal, t.a=t.a, t.b=t.b)
    probs[i] <- prod(probs1)
  }
  probs <- probs / sum(probs) # normalise
  OK <- which(probs > cutoff)
  out <- cbind(yrseq[OK], probs[OK])
  colnames(out) <- c("yr", "prob")
  write.table(out, file.path(tree.dir, paste0(name, "_probs.txt")), sep="\t", row.names=FALSE, quote=FALSE)
  if(draw)
    draw.rings(name=name, dat=dat, out=out, cc=cc, BCAD=BCAD, normal=normal, t.a=t.a, t.b=t.b, ...)
  
  # calculate offset between the (uncalibrated) dates and the calibration curve
  x0 <- out[which(out[,2] == max(out[,2])),1]
  x <- x0 - dat[,4]
  mu <- approx(cc[,1], cc[,2], x)$y
  sigma <- approx(cc[,1], cc[,3], x)$y
  offset <- sqrt((dat[,2] - mu)^2) / sqrt(dat[,3]^2 + sigma^2)
  message("Offset (in standard deviations), mean ", round(mean(offset),2), ", from ",
    round(min(offset),2), " (date ", which(offset==min(offset)), ") to ", 
	round(max(offset),2), " (date ", which(offset==max(offset)), ")\n")
  
  # for each modelled age, find whether it fits any of the date's hpd ranges
  fits <- c()
  for(i in 1:nrow(dat)) {
	calib <- caldist(dat[i,2], dat[i,3], dat[i,5])
    this.hpd <- rbind(hpd(calib, prob))
    inside <- rep(0, nrow(out))
    for(j in 1:nrow(this.hpd)) { 
  	  rng <- range(this.hpd[j,1:2])
      inside <- inside + ((out[,1] - dat[i,4] > rng[1]) * (out[,1] - dat[i,4] < rng[2]))
  	}  
  	fits[i] <- 100*length(which(inside > 0)) / nrow(out) # ratio
  }
  o <- order(fits)  
  message(round(mean(fits),2),
	"% of the model ages fit within the ", 100*prob, 
	"% hpd ranges of the dates, with worst-fitting date ", 
	o[1], " (",
    round(min(fits),2), 
	"%) and best-fitting date ", 
	o[length(o)], 
	" (", round(max(fits),2), "%)")

  invisible(list(dates=dat, probs=out, offset=offset, fit=fits))
}


