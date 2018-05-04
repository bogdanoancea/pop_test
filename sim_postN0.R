#### Simulation for postN0

fu<-list()
fv<-list()
flambda<-list()


sc = 100
cv_lambda <- 0.3
alpha <- 1 / cv_lambda**2 - 1
for(i in 1:24) {
  fu[[i]] <- list('unif', xMin = 0.3, xMax = 0.5)
  fv[[i]]<- list('unif', xMin = 1, xMax = N_0$pop_official[i]/sc+10)
  flambda[[i]]<-list('gamma', shape = 1+alpha, scale = N_0$pop_official[i]/sc/alpha)

}

system.time({
  estim_t0 <- postN0(nMNO = N_0$pop_phone/sc, nReg = N_0$pop_official/sc,
                     fu = fu, fv = fv, flambda = flambda, scale=sc, nThreads = 8)
})

