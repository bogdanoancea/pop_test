p = vector(length = 24)
for(i in 1:24) {
  p[i] = which(N_0$NIDT[i] == rownames(trans_N1N0))
}

# List of priors for v
v0 <- nReg /100
cv_v0 <- 0.10
fv <- lapply(v0, function(u){
  umin <- max(0, u - cv_v0 * u)
  umax <- u + cv_v0 * u
  output <- list('unif', xMin = umin, xMax = umax)
  return(output)
})

# List of priors for lambda
cv_lambda <- 0.6
alpha <- 1 / cv_lambda**2 - 1
flambda <- lapply(v0, function(v){list('gamma', shape = 1 + alpha, scale =  v / alpha)})

# Names and parameters of priors for the transition probabilities
distNames <- rep('unif', 24)
variation <- rep(list(list(cv = 0.20)), 24)

system.time({
  for(i in 1:40) {
    x<-transitions_all[[i]][p,p]
    x[1,]<-14
    nMNOmat = x
    nReg = N_0$pop_official
    u0 <- rowSums(nMNOmat) /nReg
    cv_u0 <- 0.15
    fu <- lapply(u0, function(u){
      umin <- max(0, u - cv_u0 * u)
      umax <- min(1, u + cv_u0 * u)
      output <- list('unif', xMin = umin, xMax = umax)
      return(output)
    })
    #print(i)
    postNt(nMNOmat/100,nReg/100, fu, fv, flambda, distNames, variation, scale = 100)
    #print(y)
  }
})
