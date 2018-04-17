library(rbenchmark)
library(pestim)


modeLambda_old <- function(nMNO, nReg, fu, fv, flambda, relTol = 1e-6, nSim = 1e4, nStrata = c(1, 1e2), verbose = FALSE, nThreads = RcppParallel::defaultNumThreads()){

  nCells <- length(nMNO)
  if (length(nReg) != nCells) stop('nReg and nMNO must have the same length.')

  mc <- match.call()
  mc[[1L]] <- NULL

  if (nCells == 1) {
    if (verbose) cat('Searching maximum...\n')

    x0 <- (nMNO + nReg) / 2
    h <- x0 / 2
    f0lr <- dlambda(c(x0, max(x0 - h, 0), x0 + h), nMNO, nReg, fu, fv, flambda, relTol, nSim, nStrata, verbose, nThreads)$probLambda
    diff0lr <- abs(f0lr[2:3] - f0lr[1])
    stopCrit <- FALSE

    while (!stopCrit){

      index.max <- which.max(f0lr)

      if (index.max == 1) {

        h <- h / 2
        f0lr <- c(f0lr[1], dlambda(c(max(x0 - h, 0), x0 + h), nMNO, nReg, fu, fv, flambda, relTol, nSim, nStrata, verbose, nThreads)$probLambda)

      } else if (index.max == 2) {

        x0 <- max(x0 - h, 0)
        f0lr <- c(f0lr[2], dlambda(max(x0 - h, 0), nMNO, nReg, fu, fv, flambda, relTol, nSim, nStrata, verbose, nThreads)$probLambda, f0lr[1])

      } else if (index.max == 3) {

        x0 <- x0 + h
        f0lr <- c(f0lr[3], f0lr[1], dlambda(x0 + h, nMNO, nReg, fu, fv, flambda, relTol, nSim, nStrata, verbose, nThreads)$probLambda)

      }
      diff0lr <- abs(f0lr[2:3] - f0lr[1])

      stopCrit <- any(diff0lr <= f0lr[1] * relTol) | h < .Machine$double.eps

    }

    xMax <- c(x0, max(x0 - h, 0), x0 + h)[index.max]
    if (verbose) cat(' ok.\n')
    return(xMax)

  } else {

    output <- sapply(seq(along = nMNO), function(i){

      if (verbose) cat(paste0('Computing for cell ', i, '...\n'))
      locnMNO <- nMNO[i]
      locnReg <- nReg[i]
      locfu.Pars <- lapply(fu[-1], '[', i)
      locfu <- c(fu[[1L]], locfu.Pars)
      locfv.Pars <- lapply(fv[-1], '[', i)
      locfv <- c(fv[[1L]], locfv.Pars)
      locflambda.Pars <- lapply(flambda[-1], '[', i)
      locflambda <- c(flambda[[1L]], locflambda.Pars)
      locMode <- modeLambda(locnMNO, locnReg, locfu, locfv, locflambda, relTol, nSim, nStrata, verbose, nThreads)
      if (verbose) cat(' ok.\n')
      return(locMode)
    })
    return(output)

  }
}

res <- benchmark(modeLambda(nMNO = c(20, 17, 25), nReg = c(115, 123, 119),
                            fu = list(list('unif', xMin = 0.3, xMax = 0.5),
                                      list('unif', xMin = 0.35, xMax = 0.45),
                                      list('unif', xMin = 0.25, xMax = 0.43)),
                            fv = list(list('gamma', shape = 11, scale = 12),
                                      list('gamma', shape = 12, scale = 12.3),
                                      list('gamma', shape = 13, scale = 11.5)),
                            flambda = list(list('gamma', shape = 11, scale = 12),
                                           list('gamma', shape = 12, scale = 12.3),
                                           list('gamma', shape = 13, scale = 12))),
                 modeLambda_old(nMNO = c(20, 17, 25), nReg = c(115, 123, 119),
                            fu = list(list('unif', xMin = 0.3, xMax = 0.5),
                                      list('unif', xMin = 0.35, xMax = 0.45),
                                      list('unif', xMin = 0.25, xMax = 0.43)),
                            fv = list(list('gamma', shape = 11, scale = 12),
                                      list('gamma', shape = 12, scale = 12.3),
                                      list('gamma', shape = 13, scale = 11.5)),
                            flambda = list(list('gamma', shape = 11, scale = 12),
                                           list('gamma', shape = 12, scale = 12.3),
                                           list('gamma', shape = 13, scale = 12)))
                 , order="relative")
res[, 1:4]

