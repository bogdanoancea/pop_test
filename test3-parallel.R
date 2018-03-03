library(pestim)
library(parallel)
library(foreach)

test_parallel<-function() {
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

cl1 <- makeCluster(length(alphaSeq))
clusterExport(cl1, c("nSim", "nMNO", "nReg", "fu", "fv", "flambdaList"))
clusterEvalQ(cl1, library(foreach))
clusterEvalQ(cl1, library(doParallel))
clusterEvalQ(cl1, library(data.table))

results <- parLapply(cl1, alphaSeq, function(alpha){
  no_cores <- detectCores() - 1
  cl2<-makeCluster(no_cores)
  registerDoParallel(cl2)
  flambda <- flambdaList[[as.character(alpha)]]
  
  output = foreach(i=1:nSim, .packages = "pestim", .combine=c, .export=c("nMNO", "nReg", "fu", "fv", "flambda", "flambdaList")) %dopar%
    postN0(nMNO, nReg, fu, fv, flambda)
  stopCluster(cl2)
  
  output <- as.data.table(t(matrix(unlist(output), nrow = 3)))
  setnames(output, c('postMean', 'postMedian', 'postMode'))
  output[, sim := 1:nSim]
  output <- melt(output, id.vars = 'sim')
  output[, 'alpha' := alpha]
  return(output)
})
stopCluster(cl1)

names(results) <- alphaSeq
results <- rbindlist(results)
ggplot(results, aes(x = variable, y = value)) +
  geom_boxplot() + facet_grid(. ~ alpha) +
  xlab('') + ylab('') +
  geom_hline(yintercept = nReg) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
}

system.time(test_parallel())
