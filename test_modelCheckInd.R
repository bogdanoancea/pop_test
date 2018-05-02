library(pestim)
library(data.table)
library(rbenchmark)



modelCheckInd2 <- function(nSimPar, nMNO, nReg, fu, fv, flambda, relTol = 1e-6, nSim = 1e6,
                          nStrata = c(1, 1e2, 1e2), verbose = FALSE,
                          nThreads = RcppParallel::defaultNumThreads()){
  n <- length(nMNO)
  if (n == 1){
    uvlambda <- ruvlambda(nSimPar, nMNO, nReg, fu, fv, flambda, relTol, nSim, nStrata, verbose, nThreads)
    nMNOrep <- lapply(1:nSimPar, function(i){
      output <- rNMNO(nSimPar,
                      lambda = uvlambda[i][['lambda']],
                      u = uvlambda[i][['u']],
                      v = uvlambda[i][['v']])
      return(output)
    })

    nSim2 = nSimPar**2
    DT <- rbindlist(nMNOrep)
    setnames(DT, 'nMNO', 'nMNOrep')
    DT[, nMNOrep2 := nMNOrep ** 2]
    DT[, difnMNOrep1 := (nMNOrep - nMNO)]
    DT[, relDifnMNOrep1 := difnMNOrep1 / nMNO]
    DT[, difnMNOrep2 := difnMNOrep1 ** 2]
    DT[, relDifnMNOrep2 := (relDifnMNOrep1) ** 2]
    indicators <- DT[, list(B = sum(difnMNOrep1) / nSim2,
                            relB = sum(relDifnMNOrep1) / nSim2,
                            m2 = sum(nMNOrep2) / nSim2,
                            e2 = ( sum(nMNOrep) / nSim2 ) ** 2,
                            rele2 = ( sum(relDifnMNOrep1) / nSim2 ) ** 2,
                            MSE = sum(difnMNOrep2) / nSim2,
                            relMSE = sum(relDifnMNOrep2) / nSim2)]

    indicators[, V := m2 - e2]
    indicators[, relV := relMSE - rele2]
    indicators[, nMNO := nMNO]
    indicators[, nReg := nReg]
    indicators <- indicators[, c('nMNO', 'nReg', 'B', 'relB', 'V', 'relV', 'MSE', 'relMSE'), with = FALSE]
    return(indicators[])

  } else {

    if(n < nThreads) {
      cl <- makeCluster(n)
    } else {
      cl <- makeCluster(nThreads)
    }
    registerDoParallel(cl)
      output <- foreach(i=1:n, .packages=c("pestim", "data.table"), .export=c("modelCheckInd2"),.combine = rbind, .options.snow = list(preschedule = TRUE)) %dopar% {
          outLocal <- modelCheckInd2(nSimPar, nMNO[i], nReg[i], fu[[i]], fv[[i]], flambda[[i]],
                                 relTol, nSim, nStrata, verbose, nThreads)
    }
    stopCluster(cl)
    return(output)

  }
}

fu = list(list('unif', xMin = 0.2, xMax = 0.25),
          list('unif', xMin = 0.2, xMax = 0.26),
          list('unif', xMin = 0.41, xMax = 0.66),
          list('unif', xMin = 0.31, xMax = 0.46),
          list('unif', xMin = 0.1, xMax = 0.16),
          list('unif', xMin = 0.71, xMax = 0.86),
          list('unif', xMin = 0.35, xMax = 0.56),
          list('unif', xMin = 0.29, xMax = 0.46)
          )

fv = list(list('unif', xMin = 115, xMax = 130),
          list('unif', xMin = 114, xMax = 131),
          list('unif', xMin = 11, xMax = 23),
          list('unif', xMin = 110, xMax = 121),
          list('unif', xMin = 141, xMax = 161),
          list('unif', xMin = 124, xMax = 171),
          list('unif', xMin = 114, xMax = 151),
          list('unif', xMin = 142, xMax = 172)
          )

flambda = list(list('gamma', shape = 21, scale = 123 / 20),
               list('gamma', shape = 11, scale = 131 / 10),
               list('gamma', shape = 11, scale = 23 / 10),
               list('gamma', shape = 11, scale = 121 / 10),
               list('gamma', shape = 11, scale = 161 / 10),
               list('gamma', shape = 11, scale = 171 / 10),
               list('gamma', shape = 11, scale = 151 / 10),
               list('gamma', shape = 11, scale = 172 / 10)
          )

a <- function() {
  set.seed(2)
  modelCheckInd(nSimPar = 10, nMNO = c(29, 31, 8, 25, 31, 28, 30, 21), nReg = c(123, 119, 19, 111, 124, 130, 132, 155),
                 fu,
                 fv,
                 flambda )
}


b <- function() {
  set.seed(2)
  pestim::modelCheckInd(nSimPar = 10, nMNO = c(29, 31, 8, 25, 31, 28, 30, 21), nReg = c(123, 119, 19, 111, 124, 130, 132, 155),
                        fu ,
                        fv ,
                        flambda)

}
stopifnot(identical(a(), b() ) )


res <- benchmark(a(), b(), order="relative", replications = 8)

res[, 1:4]

