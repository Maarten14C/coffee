


#' @name sim.rings
#' @title Simulate the radiocarbon dating of tree-rings
#' @description Simulate the dense radiocarbon dating of a tree or other deposit with exactly known yearly rings, and thus with gaps of exactly known age. The radiocarbon dates are assumed to have a degree of lab error and scatter. A (constant) offset can also be modelled.
#' @param name Name of the simulated tree-ring set. Defaults to \code{"mytree"}.
#' @param age.start Starting age of the simulated tree.
#' @param length Length of the sequence (if gaps are given as constant). Could be set at e.g 400 for an oak, but many trees will not live as long so shorter sequences could also make sense.
#' @param gaps How many calendar years there are between the dated rings. Can be set constant (one value, e.g. 20), or alternatively the gaps can be provided as a list of values.
#' @param offset The ages could be offset by some constant value. Defaults to 0.
#' @param scatter There is always a degree of scatter between measurements, and the amount of scatter can be modelled using, e.g., \code{scatter=2*error}. Set at 0 to model radiocarbon dates that are 100\% faithful to the calibration curve (very unlikely!).
#' @param error Laboratory error of the radiocarbon dates as percentage of the mean. Defaults to 0.02.
#' @param min.error Minimum laboratory to be reported. Defaults to 10 (C-14 year).
#' @param tree.dir The directory where the folders of the individual trees live. Defaults to \code{tree.dir="trees"}.
#' @param sep Separator for the fields in the .csv file. Defaults to a comma.
#' @param cc Calibration curve to be used. Could be 1 (IntCal20; default), 2 (Marine20), 3 (SHCal20) or 4 (custom curve).
#' @param postbomb Negative C-14 ages (younger than 0 cal BP or AD 1950) should be modelled using a postbomb curve. This could be 1 (northern-hemisphere region 1), 2 (NH region 2), 3 (NH region 3), 4 (southern hemisphere regions 1-2), or 5 (SH region 3).
#' @param ask Ask if a folder may be made and files written into it
#' @return A file containing 5 columns: the simulated calendar ages, the radiocarbon ages, their errors, the rings (starting with the youngest year and working backward in time), and the calibration curve to be used.
#' @examples
#'   treedir <- tempdir()
#'   sim.rings("manyrings", age.start=1000, length=400, gaps=10, tree.dir=treedir)
#'   rings("manyrings", tree.dir=treedir)
#' @author Maarten Blaauw, J. Andres Christen
#' @export
sim.rings <- function(name="mytree", age.start=1000, length=400, gaps=20, offset=0, scatter=1.05, error=0.03, min.error=15, tree.dir="trees", sep=",", cc=1, postbomb=FALSE, ask=TRUE) {
  if(length(gaps) == 1)
    sim.tree <- seq(age.start-length, age.start, by=gaps) else
      sim.tree <- rev(age.start - gaps)

  cc <- ccurve(cc, postbomb)
  if(age.start-length < 0)
    if(postbomb)
      cc <- glue.ccurves(cc,postbomb) else
        stop("for years younger than 0 cal BP (after AD 1950), a postbomb curve has to be defined (e.g., 1, 2, 3, 4 or 5", .call=FALSE)

  sim.y <- approx(cc[,1], cc[,2], sim.tree)$y + offset
  sim.er <- error * sim.y
  sim.er[sim.er < min.error] <- min.error
  sim.er <- round(sim.er)
  sim.y <- round(sim.y + rnorm(length(sim.tree), 0, scatter*sim.er))

  sim.tree <- cbind(sim.tree, sim.y, sim.er, max(sim.tree)-sim.tree, 1) # should be more flexible regarding calibration curve
  colnames(sim.tree) <- c("lab ID", "age", "error", "ring", "cc")
  
  treedir <- file.path(tree.dir, name)
  if(!dir.exists(treedir))
    assign_dir(tree.dir, name, "tree.dir", ask)
  write.table(sim.tree, file.path(treedir, paste0(name, ".csv")), row.names=FALSE, sep=sep, quote=FALSE)
}



#' @name sim.strat
#' @title Simulate the radiocarbon dating of random depths of a sediment which has accumulated over time.
#' @description Simulate the radiocarbon dating (or with dates that are already on the cal BP scale) of a deposit that is known to have accumulated over time, and for which therefore the dated depths can be safely assumed to be in chronological order.
#' @details Dates further down the sequence should have older ages than dates further up, even though owing to scatter, the dates themselves might not be in exact chronological order. The modelling is performed in a Bayesian framework (see \code{strat}). The amount of scatter, the laboratory error and an offset can also be modelled.
#' @param name Name of the simulated stratigraphy. Defaults to \code{"mystrat"}.
#' @param age.min Minimum age of the simulation.
#' @param length Length of the sequence.
#' @param n The amount of dated depths.
#' @param offset The ages could be offset by some constant value. Defaults to 0.
#' @param scatter There is always a degree of scatter between measurements, and the amount of scatter can be modelled using, e.g., \code{scatter=2*error}. Set at 0 to model radiocarbon dates that are 100\% faithful to the calibration curve (very unlikely!).
#' @param error Laboratory error of the radiocarbon dates as percentage of the mean. Defaults to 0.02.
#' @param min.error Minimum laboratory to be reported. Defaults to 10 (C-14 year).
#' @param rounded Rounding of the simulated calendar years. Rounds to single years (0 decimals) by default.
#' @param strat.dir The directory where the folders of the individual trees live. Defaults to \code{tree.dir="trees"}.
#' @param sep Separator for the fields in the .csv file. Defaults to a comma.
#' @param cc Calibration curve to be used. Could be 1 (IntCal20; default), 2 (Marine20), 3 (SHCal20) or 4 (custom curve).
#' @param postbomb Negative C-14 ages (younger than 0 cal BP or AD 1950) should be modelled using a postbomb curve. This could be 1 (northern-hemisphere region 1), 2 (NH region 2), 3 (NH region 3), 4 (southern hemisphere regions 1-2), or 5 (SH region 3).
#' @param ask Ask if a folder may be made and files written into it 
#' @return A file containing 5 columns: the simulated calendar ages, the radiocarbon ages, their errors, their relative positions (starting with the youngest one at top, and counting upward going down the sequence), and the calibration curve to be used.
#' @examples
#'   stratdir <- tempdir()
#'   sim.strat("ordered.mud", age.min=1000, length=5000, n=10, strat.dir=stratdir)
#' @author Maarten Blaauw, J. Andres Christen
#' @export
sim.strat <- function(name="mystrat", age.min=4321, length=800, n=5, offset=0, scatter=2*error, error=0.02, min.error=10, rounded=0, strat.dir="strats", sep=",", cc=1, postbomb=FALSE, ask=TRUE) {
  truth <- cumsum(runif(n)) # randomly increasing
  truth <- round(age.min + (length * (truth/max(truth))), digits=rounded)

  if(age.min > 0)
    ccc <- ccurve(cc, postbomb) else
      if(postbomb)
        ccc <- glue.ccurves(cc,postbomb) else
          stop("for years younger than 0 cal BP (after AD 1950), a postbomb curve has to be defined (e.g., 1, 2, 3, 4 or 5", .call=FALSE)

  # simulate the dating
  ages <- approx(ccc[,1], ccc[,2], truth)$y # the C14 ages
  ages <- round(ages + rnorm(length(ages), 0, scatter*ages)) # add scatter
  errors <- round(error * ages) # the lab errors
  errors[errors<min.error] <- min.error # very small errors are unlikely
  this.strat <- cbind(truth, ages, errors, 1:n, cc)
  colnames(this.strat) <- c("lab ID", "age", "error", "position", "cc")
  
  stratdir <- file.path(strat.dir, name)
  if(!dir.exists(stratdir))
    assign_dir(strat.dir, name, "strat.dir", ask)
  write.table(this.strat, file.path(stratdir, paste0(name, ".csv")), sep=sep, quote=FALSE, row.names=FALSE)
}


