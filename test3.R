library(pestim)
library(ggplot2)

test_serial<-function() {
nReg <- 97
nMNO <- 19
fu <- list('unif', xMin = 0, xMax = 0.50)
fv <- list('triang', xMin = 87, xMax = 107, xMode = 97)
alphaSeq <- c(1, 10, 100, 1000)
flambdaList <- list()
for (alpha in alphaSeq){
  flambdaList[[as.character(alpha)]] <- list('gamma', shape = 1 + alpha, scale = nReg / alpha)
}
nSim <- 2

results <- lapply(alphaSeq, function(alpha){

  flambda <- flambdaList[[as.character(alpha)]]
  output <- replicate(nSim, postNestimates(nMNO, nReg, fu, fv, flambda))
  output <- as.data.table(t(matrix(unlist(output), nrow = 3)))
  setnames(output, c('postMean', 'postMedian', 'postMode'))
  output[, sim := 1:nSim]
  output <- melt(output, id.vars = 'sim')
  output[, 'alpha' := alpha]
  return(output)
})

names(results) <- alphaSeq
results <- rbindlist(results)
ggplot(results, aes(x = variable, y = value)) +
  geom_boxplot() + facet_grid(. ~ alpha) +
  xlab('') + ylab('') +
  geom_hline(yintercept = nReg) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
}

system.time(test_serial())
