# definitions of distributions; x is the proposed value, y is the measured value (with error er)


# could also be used as approx to exponential, if shape < 1
l.gamma <-  function(x, y, shape) 
  return(dgamma(x, shape, shape/y))


l.uniform <- function(x, y, delta) # assuming cal BP
  return(dunif(x, min, min-delta))


l.after <- function(x, y) # assuming cal BP, == returns 0
  return(ifelse(x < y, 1, 0))


l.before <- function(x, y) # assuming cal BP, == returns 0
  return(ifelse(x > y, 1, 0))



# calibrate radiocarbon dates
l.calib <- function(x, y, er, cc=rintcal::ccurve(1,FALSE), normal=FALSE, t.a=3, t.b=4) {
  cc.x <- approx(cc[,1], cc[,2], x)$y
  cc.er <- approx(cc[,1], cc[,3], x)$y
  if(normal)
    prob <- dnorm(y, cc.x, sqrt(cc.er^2 + er^2)) else
      prob <- (t.b + ((y-cc.x)^2) / (2*(sqrt(er^2+cc.er^2)^2))) ^ (-1*(t.a+0.5))
  prob[is.na(prob)] <- 0
  return(prob)
      }



# for dates that don't need calibration
l.calBP <- function(x, y, er, normal=FALSE, t.a=3, t.b=4)
  if(normal)
    return(dnorm(x, y, er)) else
      return((t.b + ((y-x)^2) / (2*(er^2))) ^ (-1*(t.a+0.5)))


# some ideas, no ideas if they will work

# the functions l.before, l.after, l.uniform and perhaps l.gamma could be used as being e.g. younger than the date further down. How deal with that? 

# add a span function, x +- y yrs (e.g. for tree rings or varves with uncertainty) Or is that simply l.uniform + l.normal?

# flag as outlier (1) if moving y to x would involve a large jump... Only use for radiocarbon dates and normally distributed dates? Or not use at all, and only use student-t instead?
l.move <- function(x, y, er, times) 
  if(abs(x-y) > times*er)
    return(1) else 
      return(0)

