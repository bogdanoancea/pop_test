library(data.table)
nReg <- 97
nMNO <- 19
nPar <- 10
radShares <- seq(from = nMNO / nReg, to = 0.005, length.out = nPar)
radPopSizes <- round(seq(from = 0.25 * nReg, to = 1, length.out = nPar))
alpha <- 1
flambda <- list('gamma', shape = alpha + 1, scale = nReg / alpha)
results.Mean <- matrix(NA, ncol = nPar, nrow = nPar)
results.Median <- matrix(NA, ncol = nPar, nrow = nPar)
results.Mode <- matrix(NA, ncol = nPar, nrow = nPar)
for (radShare.index in seq(along = radShares)) {
  for (radPopSize.index in seq(along = radPopSizes)) {
    
    um <- nMNO / nReg - radShares[radShare.index]
    uM <- nMNO / nReg + radShares[radShare.index]
    fu <- list('unif', xMin = um, xMax = uM)
    
    Nm <- nReg - radPopSizes[radPopSize.index]
    NM <- nReg + radPopSizes[radPopSize.index]
    fv <- list('unif', xMin = Nm, xMax = NM)
    
    auxResults <- postNestimates(nMNO, nReg, fu, fv, flambda)
    results.Mean[radShare.index, radPopSize.index] <- auxResults[['postMean']] 
    results.Median[radShare.index, radPopSize.index] <- auxResults[['postMedian']]
    results.Mode[radShare.index, radPopSize.index] <- auxResults[['postMode']]
  }
}
rownames(results.Mean) <- round(2 * radShares, 2)
rownames(results.Median) <- round(2 * radShares, 2)
rownames(results.Mode) <- round(2 * radShares, 2)
colnames(results.Mean) <- 2 * radPopSizes
colnames(results.Median) <- 2 * radPopSizes
colnames(results.Mode) <- 2 * radPopSizes
relBias.Mean <- round((results.Mean - nReg) / nReg * 100, 1)
relBias.Median <- round((results.Median - nReg) / nReg * 100, 1)
relBias.Mode <- round((results.Mode - nReg) / nReg * 100, 1)
knitr::kable(relBias.Mean, 
             caption = 'Relative bias (%) for posterior mean estimates')
knitr::kable(relBias.Median, 
             caption = 'Relative bias (%) for posterior median estimates')
knitr::kable(relBias.Mode, 
             caption = 'Relative bias (%) for posterior mode estimates')
