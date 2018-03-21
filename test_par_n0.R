library(rbenchmark)
set.seed(123)

nMNO1 = 20
nReg1 = 115
fu1 = list('unif', xMin = 0.3, xMax = 0.5)
fv1 = list('unif', xMin = 100, xMax = 120)
flambda1 = list('gamma', shape = 11, scale = 12)

nMNO2 = 22
nReg2 = 119
fu2 = list('unif', xMin = 0.35, xMax = 0.55)
fv2 = list('unif', xMin = 85, xMax = 125)
flambda2 = list('gamma', shape = 10, scale = 10)

nMNO3 = 22
nReg3 = 119
fu3 = list('unif', xMin = 0.4, xMax = 0.6)
fv3 = list('unif', xMin = 75, xMax = 135)
flambda3 = list('gamma', shape = 9, scale = 12)

nMNO4 = 21
nReg4 = 110
fu4 = list('unif', xMin = 0.25, xMax = 0.45)
fv4 = list('unif', xMin = 100, xMax = 110)
flambda4 = list('gamma', shape = 12, scale = 15)



nMNO5 = 30
nReg5 = 180
fu5 = list('unif', xMin = 0.22, xMax = 0.35)
fv5 = list('unif', xMin = 90, xMax = 115)
flambda5 = list('gamma', shape = 11, scale = 10)

nMNO6 = 45
nReg6 = 165
fu6 = list('unif', xMin = 0.32, xMax = 0.52)
fv6 = list('unif', xMin = 75, xMax = 95)
flambda6 = list('gamma', shape = 9, scale = 8)

nMNO6 = 55
nReg6 = 265
fu6 = list('unif', xMin = 0.12, xMax = 0.72)
fv6 = list('unif', xMin = 35, xMax = 195)
flambda6 = list('gamma', shape = 10, scale = 18)

nMNO7 = 15
nReg7 = 65
fu7 = list('unif', xMin = 0.62, xMax = 0.72)
fv7 = list('unif', xMin = 50, xMax = 95)
flambda7 = list('gamma', shape = 9, scale = 11)


nMNO8 = 200
nReg8 = 500
fu8 = list('unif', xMin = 0.12, xMax = 0.72)
fv8 = list('unif', xMin = 5, xMax = 95)
flambda8 = list('gamma', shape = 19, scale = 8)





#postN0(nMNO = nMNO1 , nReg = nReg1, fu = fu1, fv = fv1, flambda = flambda1)


#postN0(nMNO = nMNO2 , nReg = nReg2, fu = fu2, fv = fv2, flambda = flambda2)

nm = c(nMNO1, nMNO2, nMNO3, nMNO4, nMNO5, nMNO6, nMNO7, nMNO8)
nr = c(nReg1, nReg2, nReg3, nReg4, nReg5, nReg6, nReg7, nReg8)

fuu = list(fu1, fu2, fu3, fu4, fu5, fu6, fu7, fu8)
fvv = list(fv1, fv2, fv3, fv4, fv5, fv6, fv7, fv8)
fl = list(flambda1, flambda2, flambda3, flambda4, flambda5, flambda6, flambda7, flambda8)

set.seed(1)
o1<-postN0(nm , nr, fuu, fvv, fl, nThreads = 1)
set.seed(1)
o2<-postN0(nm, nr, fuu, fvv, fl)

res <- benchmark(postN0(nm , nr, fuu, fvv, fl, nThreads = 1), postN0(nm , nr, fuu, fvv, fl, nThreads = 4), order="relative", replications = 4)
res[, 1:4]

Nvalues <- rN0(n=1e3, nm, nr, fuu, fvv, fl)[['N0']]

Nvalues2 <- rN0(n=1e3, c(nm[1],nm[2]), c(nr[1],nr[2]), list(fuu[[1]],fuu[[2]]), list(fvv[[1]],fvv[[2]]), list(fl[[1]], fl[[2]]))

o3<-postN0(c(nm[1],nm[2]), c(nr[1],nr[2]), list(fuu[[1]],fuu[[2]]), list(fvv[[1]],fvv[[2]]), list(fl[[1]], fl[[2]]))



nCells <-8
v<-c(1,2,3,4,5,6,7,8,9)
ptest<-function(v) {
  n<-length(v)
  if (n ==1)
  {
    output<-sqrt(v)
    return(output)
  }
  else {
    registerDoParallel(cores = nCells)
    output <- foreach(i=1:n, .combine = rbind, .options.multicore = list(preschedule = TRUE)) %dopar% {
      outlocal<-sqrt(v[i]^3)
    }
    return(output)
  }
}

stest<-function(v) {
  n<-length(v)
  if (n ==1)
  {
    output<-sqrt(v)
    return(output)
  }
  else {
    output <- lapply(seq(along = v), function(i){
      outlocal<-sqrt(v[i]^3)
    })
    output <- Reduce(rbind, lapply(output, rbind))
    return(output)
  }
}

res <- benchmark(stest(v), ptest(v), order="relative", replications = 4)
res[, 1:4]

a<-stest(v)
b<-ptest(v)




output <- lapply(seq(along = ret), function(i){

     ret[i]

   })
output <- Reduce(rbind, lapply(output, rbind))

nn<-c(1,2,3,3,2,2,1)
mode<-function(nn) {
  return (nn[which.max(names(table(nn)))])
}
mode(nn)


Mode = function(x){
  ta = table(x)
  tam = max(ta)
  if (all(ta == tam))
    mod = NA
  else
    if(is.numeric(x))
      mod = as.numeric(names(ta)[ta == tam])
  else
    mod = names(ta)[ta == tam]
  return(mod)
}
Mode(nn)
