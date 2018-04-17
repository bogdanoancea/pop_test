library(rbenchmark)
library(pestim)


rlambda_old <- function(n, nMNO, nReg, fu, fv, flambda, relTol = 1e-6, nSim = 1e4, nStrata = c(1, 1e2), verbose = FALSE, nThreads = RcppParallel::defaultNumThreads()){

  nCells <- length(nMNO)
  if (length(nReg) != nCells) stop('nReg and nMNO must have the same length.')
  #if (nCells != 1) stop('Only one cell at a time.')

  if (nCells == 1){

    ######  Computing lambdaOpt
    lambdaOpt <- modeLambda(nMNO, nReg, fu, fv, flambda, relTol, nSim, nStrata, verbose, nThreads)

    ###### Computing rejection rate  #####################
    if (verbose) cat('Computing rejection rate...')
    f <- function(x){dlambda(x, nMNO, nReg, fu, fv, flambda, relTol, nSim, nStrata, verbose, nThreads)$probLambda}
    location <- lambdaOpt
    scale <- sqrt(flambda$shape * flambda$scale^2)
    F0 <- pcauchy(0, location = location, scale = scale)
    g <- function(x){
      #dgamma(x, shape = nMNO + 1, scale = lambdaOpt / nMNO)
      dcauchy(x, location = location, scale = scale) / (1 - F0)
      #dg(x, nMNO, nReg, fu, fv, flambda, relTol, nSim, nStrata)
    }

    fun <- function(x){g(x) / f(x)}

    optimC <- 1.05 / optimise(fun, interval = c(max(location - scale, 0), location + scale))$objective
    if (verbose) cat(' ok.\n')

    if (verbose) cat('Generating and accepting/rejecting values...\n')
    u <- runif(2 * n)
    x <- qcauchy(F0 + u * ( 1 - F0), location = location, scale = scale)
    #x <- rgamma(10 * n, shape = nMNO + 1, scale = lambdaOpt / nMNO)
    v <- runif(n)
    if (verbose) cat('   of target distribution...')
    fx <- f(x)
    if (verbose) cat(' ok.\n')
    if (verbose) cat('   of candidate distribution...')
    gx <- dcauchy(x, location = location, scale = scale)
    if (verbose) cat(' ok.\n')
    output <- x[v <= fx / (optimC * gx)]
    if (verbose) cat(paste0(length(output), ' points selected.\n'))
    while (length(output) < n) {
      u <- runif(2 * n)
      x <- qcauchy(F0 + u * ( 1 - F0), location = location, scale = scale)
      #x <- rgamma(10 * n, shape = nMNO + 1, scale = lambdaOpt / nMNO)
      v <- runif(n)
      fx <- f(x)
      gx <- g(x)
      aux <- x[v <= fx / (optimC * gx)]
      if (verbose) cat(paste0(length(output), ' points selected.\n'))
      output <- c(output, aux)
    }
    output <- output[1:n]
    if (verbose) cat(' ok.\n')
    return(output)

  } else {

    output <- lapply(seq(along = nMNO), function(i){

      rlambda(n, nMNO[i], nReg[i], fu[[i]], fv[[i]], flambda[[i]], relTol, nSim, nStrata, verbose, nThreads)

    })

    output <- Reduce(cbind, output)
    dimnames(output) <- NULL
    return(output)
  }
}



res <- benchmark(rlambda(500, nMNO = 20, nReg = 115, fu = list('unif', xMin = 0.3, xMax = 0.5),
                         fv = list('unif', xMin = 100, xMax = 120),
                         flambda = list('gamma', shape = 11, scale = 12)),

                 rlambda_old(500, nMNO = 20, nReg = 115, fu = list('unif', xMin = 0.3, xMax = 0.5),
                         fv = list('unif', xMin = 100, xMax = 120),
                         flambda = list('gamma', shape = 11, scale = 12))
                 , order="relative")
res[, 1:4]
