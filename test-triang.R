library(rbenchmark)
library(pestim)


dtriang2 <- function(x, xMin, xMax, xMode){

  if (any(xMin > xMax)) stop('xMax must be greater than xMin.')
  if (!(all(xMode >= xMin) & all(xMode <= xMax))) stop('xMode must be between xMin and xMax.')


  n <- length(x)
  if (length(xMin) == 1) xMin <- rep(xMin, n)
  if (length(xMax) == 1) xMax <- rep(xMax, n)
  if (length(xMode) == 1) xMode <- rep(xMode, n)

  output <- x
  output[x <= xMin | x >= xMax] <- 0
  range1 <- (x > xMin & x <= xMode)
  output[range1] <- (2 * (x[range1] - xMin[range1])) / ((xMax[range1] - xMin[range1]) * (xMode[range1] - xMin[range1]))
  range2 <- (x >= xMode & x < xMax)
  output[range2] <- (2 * (xMax[range2] - x[range2])) / ((xMax[range2] - xMin[range2]) * (xMax[range2] - xMode[range2]))
  return(output)
}

ptriang2 <- function(q, xMin, xMax, xMode){

  if (any(xMin > xMax)) stop('xMax must be greater than xMin.')
  if (!(all(xMode >= xMin) & all(xMode <= xMax))) stop('xMode must be between xMin and xMax.')

  n <- length(q)
  if (length(xMin) == 1) xMin <- rep(xMin, n)
  if (length(xMax) == 1) xMax <- rep(xMax, n)
  if (length(xMode) == 1) xMode <- rep(xMode, n)

  output <- q
  output[q <= xMin] <- 0
  output[q >= xMax] <- 1
  range1 <- (q > xMin & q <= xMode)
  output[range1] <- ((output[range1] - xMin[range1])^2) / ((xMax[range1] - xMin[range1]) * (xMode[range1] - xMin[range1]))
  range2 <- (q > xMode & q < xMax)
  output[range2] <- 1 - ((output[range2] - xMax[range2])^2) / ((xMax[range2] - xMin[range2]) * (xMax[range2] - xMode[range2]))
  return(output)
}

rtriang2 <- function(n, xMin, xMax, xMode){
  if (any(xMin > xMax)) stop('xMax must be greater than xMin.')
  if (!(all(xMode >= xMin) & all(xMode <= xMax))) stop('xMode must be between xMin and xMax.')
  u <- runif(n)
  mc <- match.call()
  mc[[1L]] <- qtriang2
  mc[['n']] <- NULL
  mc[['q']] <- u
  output <- eval(mc)
  return(output)
}


qtriang2 <- function(q, xMin, xMax, xMode){

  if (any(xMin > xMax)) stop('xMax must be greater than xMin.')
  if (!(all(xMode >= xMin) & all(xMode <= xMax))) stop('xMode must be between xMin and xMax.')

  n <- length(q)

  output <- q
  range1 <- (q < (xMode - xMin) / (xMax - xMin))
  output[range1] <- xMin + sqrt(q[range1] * (xMax - xMin) * (xMode - xMin))
  range2 <- ( q > (xMode - xMin) / (xMax - xMin))
  output[range2] <- xMax - sqrt((1 - q[range2]) * (xMax - xMin) * (xMax - xMode))
  return(output)
}

x<-function() {
  set.seed(1)
  return (pestim::rtriang(1e6, 0, 3, 1))
}

y<-function() {
  set.seed(1)
  return (rtriang2(1e6, 0, 3, 1))
}

stopifnot(identical(x(),y()))

set.seed(1)
hist(pestim::rtriang(1e10, 0, 3, 1))
#set.seed(1)
rtriang2(1e6, 0, 3, 1)

res <- benchmark(pestim::rtriang(1e6, 0, 3, 1), rtriang2(1e6, 0, 3, 1), order="relative")
res[, 1:4]


