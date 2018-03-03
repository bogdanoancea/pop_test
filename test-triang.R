library(rbenchmark)
library(pestim)


rtriang2 <- function(n, xMin, xMax, xMode){
  if (any(xMin > xMax)) stop('xMax must be greater than xMin.')
  if (!(all(xMode >= xMin) & all(xMode <= xMax))) stop('xMode must be between xMin and xMax.')
  set.seed(123)
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

output <- pestim::rtriang(1e6, 0, 3, 1)
stopifnot(identical(pestim::rtriang(1e6, 0, 3, 1), rtriang2(1e6, 0, 3, 1)))

res <- benchmark(pestim::rtriang(1e6, 0, 3, 1), rtriang2(1e6, 0, 3, 1), order="relative")
res[, 1:4]
