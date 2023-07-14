####  Program (package) with the t-walk implementation in R
####  see http://www.cimat.mx/~jac/twalk/
####  Author and copyright holder J. Andres Christen

### this version adapted from the Rtwalk CRAN package version 1.8.0 (2015-09-18), which was released under the GPL-3 license. Instances of 'cat' were replaced with 'message' and txtProgressBar has been added. Some if statements have been rewritten as they were not necessary (see commented lines)

###############################################################
#### Some auxiliary functions and constants:

IntProd <- function(x) { sum(x*x) } ## square of the norm
DotProd <- function(x, y) { sum(x*y) } ## dot product



########## h1 function for the "traverse" kernel 
Simh1 <- function( dim, pphi, x, xp, beta) {
  phi <- (runif(dim) < pphi)
  rt <- NULL
  for(i in 1:dim)
    if(phi[i])
      rt <- append( rt, xp[i] + beta*(xp[i] - x[i])) else
        rt <- append( rt, x[i])
  list(rt=rt, nphi=sum(phi))
}



### Simulation of the beta parameter for kernel h1
Simfbeta <- function(at)
  if(runif(1) < (at-1)/(2*at))
    exp(1/(at + 1)*log(runif(1))) else
      exp(1/(1 - at)*log(runif(1)))



########## h function for the "walk" kernel 
Simh2 <- function( dim, pphi, aw, x, xp) {
  u <- runif(dim)
  phi <- (runif(dim) < pphi)
  z <- (aw/(1+aw))*(aw*u^2 + 2*u -1)
  z <- z*phi
  list( rt=x + (x - xp)*z, nphi=sum(phi))
}



########## h function for the "blow" kernel 
Simh3 <- function( dim, pphi, x, xp) {
  phi <- (runif(dim) < pphi)
  sigma <- max(phi*abs(xp - x))
  x + sigma*rnorm(dim)*phi # what does this do? It doesn't define anything
  list(rt=xp*phi + sigma*rnorm(dim)*phi + x*(1-phi), nphi=sum(phi), phi=phi)
}



## -log of g3, the density of h_b
G3U <- function( nphi, phi, h, x, xp) {
  sigma <- max(phi*abs(xp - x)) ## different sigma for the proposal and the reverse, but same phi
  if(nphi > 0)
    (nphi/2)*log(2*pi) + nphi*log(sigma) + 0.5*IntProd(h - xp)/(sigma^2) else
      0
}



########## h function for the "hop" kernel 
Simh4 <- function( dim, pphi, x, xp) {
  phi <- (runif(dim) < pphi)
  sigma <- max(phi*abs(xp - x))/3
  rt <- NULL
  for(i in 1:dim)
    if (phi[i])
      rt <- append( rt, x[i] + sigma*rnorm(1)) else
        rt <- append( rt, x[i])

  return(list(rt=rt, nphi=sum(phi), phi=phi))
}



log2pi <- log(2*pi); log3 <- log(3)
## -log of g4, the density of h_h.
G4U <- function( nphi, phi, h, x, xp) {
  sigma <- max(phi*abs(xp - x))/3 ## different sigma for the proposal and the reverse, but same phi
  if(nphi > 0)
    (nphi/2)*log2pi - nphi*log3 + nphi*log(sigma) + 0.5*9*IntProd((h - x))/(sigma^2) else
      0
}



############ This is the twalk implementation.  It requires the "energy" of the
############ objective function ie. U = -log f and its support.  Also the initial
############ values x0 and x0p.  Tr is the number of iterations required.
############ The dimension is dim, defaults to the global variable n
############ in the calling NAMESPACE for backward compatibility. 
############ The rest of the parameters are for dynamically plotting the trajectories
############ of the twalk for objectives of dim=2. Otherwise set
############ PlotObj=FALSE. See examples.R for details.
############ 4 kernels: traverse, walk, blow or hop
############ MB: removed plot options, now only stores every 'thinning' iterations
Runtwalk <- function(Tr, Obj, Supp, dat, dim = length(x0), x0=x0, xp0=xp0, at=6, aw=1.5, pphi=min( dim, 4)/dim, F1=0.4918, F2=F1+0.4918, F3=F2+0.0082, thinning=100, cumulative=FALSE, out.fl=c(), energy.fl=c(), ...) {

  ## Initial values
  x <- x0
  xp <- xp0

  if(Supp(x) && Supp(xp)) { # values in support
    flush(stdout())
    U <- Obj(x, ...)
    Up <- Obj(xp, ...)
    acc <- 0
    flush(stdout()) # why this again?
    #rec <- matrix(0, ncol=(2+2*dim)) ## to save U, Up, x and xp
    rec <- array(c(U, Up, x, xp))
    recacc <- array(c(0, 0))
    flush(stdout()) # third time? Probably doesn't delay things by much as only called 3 times
  } else {
      stop(paste("Initial values out of support,\n  x=", x, "\n xp=", xp))
      Tr <- 0
    }

  if(any(abs(x0 - xp0) <= 0)) {
    stop(paste("Not all entries of initial values different,\n  x=", x, "\n xp=", xp))
    Tr <- 0
  }

  every <- ceiling(Tr/thinning) # find how many sub-runs to run
  pb <- txtProgressBar(min=0, max=Tr, style = 3, char=">")

  if(length(out.fl) > 0) { # then we'll save the its in files while running
    if(file.exists(out.fl))
     file.remove(out.fl)
    if(file.exists(energy.fl))
      file.remove(energy.fl)
  }
	  
  acc <- 0
  j <- 0
  for(i in 1:Tr) {
    if(i %% 100)
      setTxtProgressBar(pb, i)
    move <- OneMove(dim=dim, Obj=Obj, Supp=Supp, x, U, xp, Up, at=at, aw=aw, pphi=pphi, F1=F1, F2=F2, F3=F3, ...)
    if(runif(1) < move$A) { # proposed iteration accepted
      tmp.recacc <- c(move$funh, move$nphi/dim)
      acc <- acc + move$nphi/dim
      U <- move$propU
      Up <- move$propUp
      x <- move$y
      xp <- move$yp
    }  else # proposed iteration rejected
         tmp.recacc <- c(move$funh, 0)

    # sub-runs to reduce the amount of iterations to be stored
    j <- j+1
    if(j == thinning) { # then store the iteration
	  if(length(out.fl) > 0) {
        fastwrite(t(c(U, Up, x, xp)), out.fl, append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)		  
        fastwrite(t(tmp.recacc), energy.fl, append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
	  } else {
        rec <- rbind(rec, c(U, Up, x, xp))
        recacc <- rbind(recacc, tmp.recacc) # add
      }
      j <- 0 # reset
    }
  }
  message("\nDone running...")

  if(length(out.fl) > 0) {
    rec <- fastread(out.fl, header=FALSE, sep=",") 
    recacc <- fastread(energy.fl, header=FALSE, sep=",")
  }

  list( dim=dim, Tr=Tr, acc=acc, Us=rec[,1], Ups=rec[,2],
    output=rec[,2+(1:dim)], outputp=rec[,2+dim+(1:dim)],
    recacc=recacc, k=every) # added k as parameter
}



# choose a random kernel to move a parameter
# 4 kernels: traverse, walk, blow or hop
# traverse if runif(1) <= F1, else walk if <= F2, else blow if <= F3, else hop
OneMove <- function( dim, Obj, Supp, x, U, xp, Up, at=6, aw=1.5, pphi=min( dim, 4)/dim, F1=0.4918, F2=0.9836, F3=0.9918, ...) {

  ker <- runif(1)  ## randomly choose a kernel

  if(ker < F1) { ## the t-walk, kernel h1: traverse
    #dir <- runif(1)
    funh <- 1

    # if ((0 <= dir) && (dir < 0.5)) {
    if(runif(1) < 0.5) { # no need to test for dir >= 0
      beta <- Simfbeta(at)
      tmp <- Simh1(dim, pphi, xp, x, beta)
      yp <- tmp$rt
      nphi <- tmp$nphi
      y  <- x
      propU <- U

      if(Supp(yp)) { # all parameters in support
        propUp <- Obj(yp, ...)
        ## The proposal is symmetric
        if(nphi == 0) ### Nothing moved
          A <- 1 else
            A <- exp((U - propU) + (Up - propUp) +  (nphi-2)*log(beta))
        } else {
            propUp <- NULL
            A <- 0  ## out of support, not accepted
          }
    } else {
        # if ((0.5 <= dir) && (dir < 1.0)) no need for this additional test
        beta <- Simfbeta(at)
        tmp <- Simh1( dim, pphi, x, xp, beta) # this differs from above when dir < 0.5
        y <- tmp$rt
        nphi <- tmp$nphi
        yp  <- xp
        propUp <- Up

        if(Supp(y)) {
          propU <- Obj(y, ...)
          ## The proposal is symmetric
          if(nphi == 0)
            A <- 1 else ### Nothing moved
              A <- exp((U - propU) + (Up - propUp) +  (nphi-2)*log(beta))
        } else {
            propU <- NULL
            A <- 0  ## out of support, not accepted
          }
      }
  } else # MB Aug 2021

      if(ker < F2) { ## the t-walk, kernel h2: walk
        # if ((F1 <= ker) && (ker < F2)) { ## the t-walk, kernel h2: walk
        #  dir <- runif(1)
        funh <- 2

#      if ((0 <= dir) && (dir < 0.5))  ## x as pivot
        if(runif(1) < 0.5) { ## x as pivot
          tmp <- Simh2( dim, pphi, aw, xp, x)
          yp <- tmp$rt
          nphi <- tmp$nphi

          y  <- x
          propU <- U

          if(Supp(yp) && all(abs(yp - y) > 0)) { # this differs from kernel h1
            propUp <- Obj(yp, ...)
            A <- exp((U - propU) + (Up - propUp))
          } else {
              propUp <- NULL
              A <- 0  ## out of support, not accepted
            }
        } else { ## xp as pivot
#			if ((0.5 <= dir) && (dir < 1.0)) { ## xp as pivot
            tmp <- Simh2( dim, pphi, aw, x, xp)
            y <- tmp$rt
            nphi <- tmp$nphi

            yp  <- xp
            propUp <- Up

            if(Supp(y) && all(abs(yp - y) > 0)) {
              propU <- Obj(y, ...)
              A <- exp((U - propU) + (Up - propUp))
            } else {
                propU <- NULL
                A <- 0  ## out of support, not accepted
              }
          }
      } else # MB Aug 2021

          if(ker < F3) {  ## the t-walk, kernel h3: blow

#		if ((F2 <= ker) && (ker < F3)) { ## the t-walk, kernel h3: blow
          #  dir <- runif(1)
            funh <- 3

#			if ((0 <= dir) && (dir < 0.5))  ## x as pivot
            if(runif(1) < 0.5) { ## x as pivot {
              tmp <- Simh3( dim, pphi, xp, x)
              yp <- tmp$rt
              nphi <- tmp$nphi
              phi <- tmp$phi

              y  <- x
              propU <- U
              if(Supp(yp) && all(yp != x)) {
                propUp <- Obj(yp, ...)
                W1 <- G3U( nphi, phi, yp, xp,  x)
                W2 <- G3U( nphi, phi, xp, yp,  x)
                A <- exp((U - propU) + (Up - propUp) +  (W1 - W2))
              } else {
                  propUp <- NULL
                  A <- 0  ## out of support, not accepted
                }
            } else { ## xp as pivot

#			if ((0.5 <= dir) && (dir < 1.0)) { ## xp as pivot
                tmp <- Simh3( dim, pphi, x, xp)
                y <- tmp$rt
                nphi <- tmp$nphi
                phi <- tmp$phi

                yp  <- xp
                propUp <- Up

                if(Supp(y) && all(y != xp)) {
                  propU <- Obj(y, ...)
                  W1 <- G3U( nphi, phi, y,  x, xp)
                  W2 <- G3U( nphi, phi, x,  y, xp)
                  A <- exp((U - propU) + (Up - propUp) +  (W1 - W2))
                } else {
                    propU <- NULL
                    A <- 0  ## out of support, not accepted
                  }
              }
          } else { ## the t-walk, kernel h4: hop
              # if (F3 <= ker) { ## the t-walk, kernel h4: hop
              # dir <- runif(1)
              funh <- 4
              if(runif(1) < 0.5) { ## x as pivot
                tmp <- Simh4( dim, pphi, xp, x)
                yp <- tmp$rt
                nphi <- tmp$nphi
                phi <- tmp$phi

                y  <- x
                propU <- U

                if(Supp(yp) && all(yp != x)) {
                  propUp <- Obj(yp, ...)
                  W1 <- G4U( nphi, phi, yp, xp,  x)
                  W2 <- G4U( nphi, phi, xp, yp,  x)
                  A <- exp((U - propU) + (Up - propUp) +  (W1 - W2))
                } else {
                    propUp <- NULL
                    A <- 0  ## out of support, not accepted
                  }
              } else { ## xp as pivot
#			if ((0.5 <= dir) && (dir < 1.0)) { ## xp as pivot
                  tmp <- Simh4( dim, pphi, x, xp)
                  y <- tmp$rt
                  nphi <- tmp$nphi
                  phi <- tmp$phi

                  yp  <- xp
                  propUp <- Up

                  if(Supp(y) && all(y != xp)) {
                    propU <- Obj(y, ...)
                    W1 <- G4U( nphi, phi, y,  x, xp)
                    W2 <- G4U( nphi, phi, x,  y, xp)
                    A <- exp((U - propU) + (Up - propUp) +  (W1 - W2))
                  } else {
                     propU <- NULL
                      A <- 0  ## out of support, not accepted
                    }
                }
          }

 if(is.nan(A) || is.na(A))  #### debugging line
    message("Rtwalk: ERROR in evaluating the objective.  Value returned by objective function:", propU)

  return(list( y=y, propU=propU, yp=yp, propUp=propUp, A=A, funh=funh, nphi=nphi))
}




######### Basic functions to analyse the output of Runtwalk, commonly saved in info

######### To plot a time series of the log of the posteriors
PlotLogObj <- function(info, from=0, to=length(info$Us))
  plot( from:(to-1), -info$Us[(from+1):to], type="l",
    ylab="Log of Objective", xlab="Iteration", main="")



######### To plot a histogram of any parameter
PlotHist <- function( info, par=1, from=0, xlab=paste("Parameter", par), main="", ...)
  if (info$dim == 1)
    hist( info$output[from:(info$Tr)], xlab=xlab, main=main, ...) else
      hist( info$output[from:(info$Tr), par], xlab=xlab, main=main, ...)



SaveOutput <- function( info, file, pars=1:(info$dim), from=1, to=info$Tr, row.names=FALSE, col.names=paste("X", pars), ...)
  if(info$dim == 1)
    write.table(info$output[from:(info$Tr)], file=file, row.names=row.names, col.names=col.names, ...) else
        write.table(info$output[from:(info$Tr), pars], file=file, row.names=row.names, col.names=col.names, ...)



############ To plot an outline of the output:
####  Calculate IAT and acceptance ratios
Ana <- function(info, from=1, to=info$Tr, par=0, file="") {
  sel <- from:to

  accrt <- rep(0, 4)
  for(h in 1:4) {
    selh <- which(info$recacc[sel,1] == h)
    accrt[h] <- sum(info$recacc[selh,2])/length(selh)
  }

  #### No plots
  Tint <- IAT( info, par=par)  ### defined below
  itmap = which(-info$Us == max(-info$Us))[1]

  message( "Ratio of moved coordinates per it=\n",
    accrt[1], accrt[2], accrt[3], accrt[4],
    "\ndim=", info$dim, "AcceptanceRatio=", info$acc/info$Tr,
    "MAPlogPost=", -info$Us[itmap], "IAT=", Tint, "IAT/dim=", Tint/info$dim,"\n\n")

  if(file != "")
    message(file=file, info$dim,
      accrt[1], accrt[2], accrt[3], accrt[4], info$acc/info$Tr,
      -info$Us[itmap], Tint/info$dim,"\n")
}



### Plot time series of parameters
TS <- function(info, pars=1:(info$dim), from=1, to=info$Tr, prime=FALSE) {
  sel <- from:to
  if(length(pars) <= 10) {
    if(info$dim == 1)
      if(!prime)
        plot(as.ts(as.matrix(info$output[sel])), main="x") else
          plot(as.ts(as.matrix(info$outputp[sel])), main="xp")
        else
          if(!prime)
              plot(as.ts(as.matrix(info$output[sel, pars])),  main="x") else
                plot(as.ts(as.matrix(info$outputp[sel, pars])), main="xp")
    } else
      message("Cannot print time series for more than 10 parameters, select a subset with arg. pars\n\n")
}



###### These functions are for calculating and plotting autcorrelations and
###### Integrated Autocorrelation Times



GetAutoCorr <- function( info, par=0, from=1, to=info$Tr, lag=30*info$dim)
  if(par>0)
    cor( info$output[from:(to-lag), par], info$output[(from+lag):to, par]) else
      cor( info$Us[from:(to-lag)], info$Us[(from+lag):to])



GetAutoCov <- function( dt, lags) {
  n <- length(dt)
  aut <- rep( 0.0, length(lags))
  mu <- mean(dt)

  for (i in 1:length(lags)) {
    lg <- lags[i]
    aut[i] <- sum( (dt[1:(n-lg)]-mu)*(dt[(lg+1):n]-mu) )/n
    #cat( i, lg, aut[i], "\n")
  }

  return(aut)
}



#' @name IAT
#' @title calculate the Integrated Autocorrelation Time
#' @description Calculate the Tntegrated Autocorrelation Time, which gives the proposed value for thinning. E.g., if the IAT is 80, it is good to thin the MCMC run by storing only every 80 iterations. This method is slower than GetAutoCov, but much better.
#' @param set This option reads the 'info' variable, which contains the data and the model output.
#' @param par The parameter to test. Defaults to 0.
#' @param from The first of the iterations of the MCMC run to check. Defaults to the first one.
#' @param to The last of the iterations of the MCMC run to check. Defaults to the last one.
#' @return The IAT value
#' @author Andres Christen
#' @export
IAT <- function(set, par=0, from=1, to) {
  ## -lag/log(GetAutoCorr( info, lag, par=par, from=from, to=to))
  ## we get the desired time series, the parameter and the from - to selection
  if(par>0) {
    if(set$dim > 1)
      dt <- set$output[from:to, par] else
        dt <- set$output[from:to]
  } else
    dt <-  set$Us[from:to]

  n <- to-from
  mu <- mean(dt)  ### with its mean and variance
  s2 <- var(dt)

  ### The maximum lag is half the sample size
  maxlag <- max( 3, floor(n/2))

  #### The gammas are sums of two consecutive autocovariances
  Ga <- rep(0,2)  ## two consecutive gammas

  lg <- 0
  Ga[1] <- s2  #sum( (dt[1:(n-lg)]-mu)*(dt[(lg+1):n]-mu) )/n
  lg <- 1
  Ga[1] <- Ga[1] + sum( (dt[1:(n-lg)]-mu)*(dt[(lg+1):n]-mu) )/n

  m <- 1
  lg <- 2*m
  Ga[2] <- sum( (dt[1:(n-lg)]-mu)*(dt[(lg+1):n]-mu) )/(n-lg)
  lg <- 2*m+1
  Ga[2] <- Ga[2] + sum( (dt[1:(n-lg)]-mu)*(dt[(lg+1):n]-mu) )/(n-lg)

  IAT <- Ga[1]/s2  ### Add the autocorrelations

  ### RULE: while Gamma stays positive and decreasing
  while((Ga[2] > 0.0) && (Ga[2] < Ga[1])) {
    m <- m+1
    if(2*m+1 > maxlag) {
      message("Not enough data, maxlag=", maxlag, "\n")
      break
    }
    Ga[1] <- Ga[2]

    lg <- 2*m
    Ga[2] <- sum( (dt[1:(n-lg)]-mu)*(dt[(lg+1):n]-mu) )/n
    lg <- 2*m+1
    Ga[2] <- Ga[2] + sum( (dt[1:(n-lg)]-mu)*(dt[(lg+1):n]-mu) )/n

    IAT <- IAT + Ga[1]/s2
  }

  IAT <- -1 + 2*IAT   ##Calculates the IAT from the gammas
  return(IAT)
}
