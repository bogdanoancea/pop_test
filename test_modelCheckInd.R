
modelCheckInd <- function(nSimPar, nMNO, nReg, fu, fv, flambda, relTol = 1e-6, nSim = 1e6,
                          nStrata = c(1, 1e2, 1e2), verbose = FALSE,
                          nThreads = RcppParallel::defaultNumThreads()){

  if (length(nMNO) == 1){

    uvlambda <- ruvlambda(nSimPar, nMNO, nReg, fu, fv, flambda, relTol, nSim, nStrata, verbose, nThreads)

    nMNOrep <- lapply(1:nSimPar, function(i){

      output <- rNMNO(nSimPar,
                      lambda = uvlambda[i][['lambda']],
                      u = uvlambda[i][['u']],
                      v = uvlambda[i][['v']])
      return(output)
    })

    DT <- rbindlist(nMNOrep)
    setnames(DT, 'nMNO', 'nMNOrep')
    DT[, nMNOrep2 := nMNOrep ** 2]
    DT[, difnMNOrep1 := (nMNOrep - nMNO)]
    DT[, relDifnMNOrep1 := (nMNOrep - nMNO) / nMNO]
    DT[, difnMNOrep2 := (nMNOrep - nMNO) ** 2]
    DT[, relDifnMNOrep2 := ((nMNOrep - nMNO) / nMNO) ** 2]
    indicators <- DT[, list(B = sum(difnMNOrep1) / nSimPar**2,
                            relB = sum(relDifnMNOrep1) / nSimPar**2,
                            m2 = sum(nMNOrep2) / nSimPar**2,
                            e2 = ( sum(nMNOrep) / nSimPar**2 ) ** 2,
                            relm2 = sum(relDifnMNOrep2) / nSimPar**2,
                            rele2 = ( sum(relDifnMNOrep1) / nSimPar**2 ) ** 2,
                            MSE = sum(difnMNOrep2) / nSimPar ** 2,
                            relMSE = sum(relDifnMNOrep2) / nSimPar ** 2)]
    indicators[, V := m2 - e2]
    indicators[, relV := relm2 - rele2]
    indicators[, nMNO := nMNO]
    indicators[, nReg := nReg]
    indicators <- indicators[, c('nMNO', 'nReg', 'B', 'relB', 'V', 'relV', 'MSE', 'relMSE'), with = FALSE]
    return(indicators[])

  } else {

    output <- lapply(seq(along = nMNO), function(i){

      modelCheckInd(nSimPar, nMNO[i], nReg[i], fu[[i]], fv[[i]], flambda[[i]],
                    relTol, nSim, nStrata, verbose, nThreads)
    })

    output <- rbindlist(output)
    return(output[])
  }
}
