library(pestim)
library(data.table)
library(rbenchmark)



modelCheckInd2 <- function(nSimPar, nMNO, nReg, fu, fv, flambda, relTol = 1e-6, nSim = 1e6,
                          nStrata = c(1, 1e2, 1e2), verbose = FALSE,
                          nThreads = RcppParallel::defaultNumThreads()){
cat("I am here")
  if (length(nMNO) == 1){

    uvlambda <- ruvlambda(nSimPar, nMNO, nReg, fu, fv, flambda, relTol, nSim, nStrata, verbose, nThreads)

    nMNOrep <- lapply(1:nSimPar, function(i){

      output <- rNMNO(nSimPar,
                      lambda = uvlambda[i][['lambda']],
                      u = uvlambda[i][['u']],
                      v = uvlambda[i][['v']])
      return(output)
    })

    nSim2 = nSimPar**2
    setnames(DT, 'nMNO', 'nMNOrep')
    DT[, nMNOrep2 := nMNOrep ** 2]
    DT[, difnMNOrep1 := (nMNOrep - nMNO)]
    DT[, relDifnMNOrep1 := (nMNOrep - nMNO) / nMNO]
    DT[, difnMNOrep2 := (nMNOrep - nMNO) ** 2]
    DT[, relDifnMNOrep2 := ((nMNOrep - nMNO) / nMNO) ** 2]
    indicators <- DT[, list(B = sum(difnMNOrep1) / nSim2,
                            relB = sum(relDifnMNOrep1) / nSim2,
                            m2 = sum(nMNOrep2) / nSim2,
                            e2 = ( sum(nMNOrep) / nSim2 ) ** 2,
                            rele2 = ( sum(relDifnMNOrep1) / nSim2 ) ** 2,
                            MSE = sum(difnMNOrep2) / nSim2,
                            relMSE = sum(relDifnMNOrep2) / nSim2)]
    cat(indicators)

    indicators[, V := m2 - e2]
    indicators[, relV := relMSE - rele2]
    indicators[, nMNO := nMNO]
    indicators[, nReg := nReg]
    indicators <- indicators[, c('nMNO', 'nReg', 'B', 'relB', 'V', 'relV', 'MSE', 'relMSE'), with = FALSE]
    return(indicators[])

  } else {

    output <- lapply(seq(along = nMNO), function(i){

      modelCheckInd2(nSimPar, nMNO[i], nReg[i], fu[[i]], fv[[i]], flambda[[i]],
                    relTol, nSim, nStrata, verbose, nThreads)
    })

    output <- rbindlist(output)
    return(output[])
  }
}



a <- function() {
  set.seed(2)
  modelCheckInd2(nSimPar = 10, nMNO = c(29, 31), nReg = c(123, 119),
                 fu = list(list('unif', xMin = 0.2, xMax = 0.25),
                           list('unif', xMin = 0.21, xMax = 0.26)),
                 fv = list(list('unif', xMin = 115, xMax = 130),
                           list('unif', xMin = 114, xMax = 131)),
                 flambda = list(list('gamma', shape = 21, scale = 123 / 20),
                                list('gamma', shape = 11, scale = 124 / 10)))
}


b <- function() {
  set.seed(2)
  pestim::modelCheckInd(nSimPar = 10, nMNO = c(29, 31), nReg = c(123, 119),
                        fu = list(list('unif', xMin = 0.2, xMax = 0.25),
                                  list('unif', xMin = 0.21, xMax = 0.26)),
                        fv = list(list('unif', xMin = 115, xMax = 130),
                                  list('unif', xMin = 114, xMax = 131)),
                        flambda = list(list('gamma', shape = 21, scale = 123 / 20),
                                       list('gamma', shape = 11, scale = 124 / 10)))

}
stopifnot(identical(a(), b() ) )


res <- benchmark(a(), b(), order="relative", replications = 1)

res[, 1:4]

