# ==========================================
#     TEST EUROPEAN PACKAGE pestim
# ==========================================
library(dplyr)
library(ggplot2)
library(pestim)
#library(doSNOW)
library(foreach)
library(pestim)


rm(list = ls())

path <- "D:/W3CRK9/Mes Documents/2018-04-25 test package européen"

dfindiv <- readr::read_csv(paste0(path,"/output.csv"), col_names = c("CodeId","timestamp","NIDT"))

df_flux <- dfindiv %>% group_by(CodeId) %>%
  dplyr::mutate(n = row_number()) %>% filter(n()>1) %>%
  arrange(CodeId, n) %>% dplyr::mutate(origin = lag(NIDT,1), current = NIDT, destination = lead(NIDT,1))


# stat1 <- df_flux %>% ungroup() %>% arrange(CodeId,timestamp) %>%
#   group_by(CodeId) %>%
#   dplyr::mutate(diff = timestamp - lag(timestamp)) %>%
#   group_by(CodeId) %>% summarise(moyenne_diff = mean(diff, na.rm=T)/3600, mediane_diff = median(diff,na.rm = T)/3600)
#
# stat1 <- reshape2::melt(stat1, id.vars = 'CodeId') %>% tbl_df()
#
# ggplot(data = stat1) + geom_histogram(aes(x = value, fill = variable))
# ggplot(data = stat1) + stat_ecdf(aes(x = value, color = variable)) + scale_x_continuous(limits = c(0,80),
#                                                                                         breaks=seq(0,80,10))


# =====================
# CREATE FUNCTIONS
# =====================

#fluxdata = df_flux

#' Aggregate timestamps given into time intervals
#'
#' From a dataframe with timestamp measures, create time intervals of fixed length
#' @param fluxdata Dataframe or tibble
#' @param length.interval Length of interval in hour
#' @param time.var Timestamp variable name in \code{fluxdata}
#' @return \code{fluxdata} with two new columns

create_time <- function(fluxdata, length.interval = 6, time.var = "timestamp"){

  fluxdata %>% ungroup() %>%
    dplyr::mutate_(.dots = setNames(c(sprintf("lubridate::day(%s)",time.var),
                                      sprintf("1+floor(lubridate::hour(%s)/%d)",
                                              time.var,length.interval)),c("jour","tranche"))) %>%
    dplyr::mutate(tranche = 4*(jour - min(jour,na.rm=T)) + tranche)


}

#' Compute probability of being to a given antenna by time interval
#'
#' @param fluxdata Dataframe


proba_tranche <- function(fluxdata){

  fluxdata %>% group_by(tranche, CodeId, NIDT) %>%
    dplyr::summarise(N = n()) %>% ungroup() %>%
    group_by(tranche, CodeId) %>%
    dplyr::mutate(N = N/sum(N)) %>%
    rename(proba = N)

}

#' Compute population at antenna level for each
#'  time interval

population_tranche <- function(fluxdata){

  df_flux2 <- fluxdata %>% proba_tranche() %>%
    filter(!is.na(tranche)) %>%
    group_by(tranche, NIDT) %>% dplyr::summarise(pop = sum(proba))

}


#' Compute benchmark population (time t0 population)
#'
#' Compute voronoi where people live using max activity
#' during the night
#'
#' @param fluxdata Dataframe
#' @param hour.home Time for night to be considered
#' @return Dataframe where people are assigned a voronoi of leaving

population_t0 <- function(fluxdata, hour.home = c(21,8)){

  df_flux0 <- fluxdata %>%
    dplyr::filter(between(lubridate::hour(timestamp),max(hour.home),24)|between(lubridate::hour(timestamp),0,min(hour.home))) %>%
    group_by(CodeId,NIDT) %>% dplyr::mutate(N = n()) %>%
    group_by(CodeId) %>% filter(N == max(N, na.rm=T)) %>%
    sample_n(1) %>% select(CodeId,NIDT) %>%
    dplyr::mutate(proba = 1, tranche = 0)

}

#' Row: NIDT of origin, Column: destination NIDT

compute_transition <- function(dataflux,df_flux0, t = 1){

  if (!between(t,1,max(dataflux$tranche))) stop("Time considered should be between 1 and max t available in CDR")

  df_fluxt <- dataflux %>% group_by(tranche, CodeId, NIDT) %>%
    dplyr::summarise(N = n()) %>% ungroup() %>% filter(tranche == t) %>%
    group_by(CodeId) %>%
    dplyr::mutate(N = N/sum(N)) %>%
    rename(proba = N)


  df_flux_01 <- rbind(df_flux0,df_fluxt)

  df_flux_01_bis <- df_flux_01 %>% group_by(CodeId) %>% filter(n()>1) %>%
    group_by(CodeId,tranche) %>% mutate(listNIDT1 = paste(sort(unique(NIDT)),collapse=","), p = 1/length(listNIDT1))
  df_flux_01_bis1 = df_flux_01_bis %>% filter(tranche==0)
  df_flux_01_bis2 = df_flux_01_bis%>% filter(tranche==t) %>% filter(row_number()>1)
  df_flux_01_bis = rbind(df_flux_01_bis1,df_flux_01_bis2) %>% group_by(CodeId) %>% mutate(listNIDT_res = lag(NIDT)) %>%
    filter(tranche==1)

  mat <- matrix(0, nrow = length(unique(df_flux_01$NIDT)), ncol = length(unique(df_flux_01$NIDT)))
  rownames(mat) = unique(df_flux_01$NIDT)
  colnames(mat) = unique(df_flux_01$NIDT)

  # Matrice transition
  for (i in 1:nrow(df_flux_01_bis)){
    nam = as.character(stringr::str_split(df_flux_01_bis$listNIDT1[i],",",simplify = T))
    for (nidt in nam){
      mat[df_flux_01_bis$NIDT[i],nidt] <- mat[df_flux_01_bis$NIDT[i],nidt] + 1/df_flux_01_bis$p[i]
    }
  }

  return(mat)
}


compute_transition_all <- function(dataflux, df_flux0){

  tseq <- sort(unique(dataflux$tranche))

  cl <- parallel::makeCluster(parallel::detectCores()-1)
  doSNOW::registerDoSNOW(cl)

  pb <- txtProgressBar(max = length(tseq), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  # APPELS PARALLELISES DE LA FONCTION
  transitions_all <- foreach(t = tseq, .combine = "list",
                    .multicombine = TRUE,
                    .maxcombine = length(tseq),
                    .options.snow=opts,
                    .export = c('dataflux','compute_transition','df_flux0'),
                    .packages = c("dplyr","stringr")) %dopar% {
                                    return(
                                      compute_transition(dataflux, df_flux0, t=1)
                                    )
                                  }

  # STOP LES CLUSTERS
  parallel::stopCluster(cl)


  return(transitions_all)
}


compare_popN0 <- function(df_flux2){

  # Pop in filosofi
  df_pop <- readr::read_csv(paste0(path,"/popmarseille.csv")) %>% filter(NIDT %in% unique(df_flux2$NIDT)) %>%
    rename(pop_official = pop)

  df_flux2 <- df_flux2 %>% rename(pop_phone = pop)

  return(df_pop %>% left_join(.,df_flux2))

}



# CREATE TIME INTERVAL
df_flux <- df_flux %>% create_time()

# COMPUTE POPULATION AT ANTENNA LEVEL FOR EACH TIMEPERIOD
df_flux2 <- df_flux %>% population_tranche()

# COMPUTE WHERE PEOPLE LIVE
df_flux0 <- df_flux %>% population_t0()

# TRANSITION TO T=1
trans_N1N0 <- compute_transition(dataflux = df_flux, df_flux0, t=1)

# ALL TRANSITIONS FROM T=0 TO T=t
transitions_all <- df_flux %>% compute_transition_all(., df_flux0)

# COMPARE POPULATION AT N0 BETWEEN TWO SOURCES
N_0 <- compare_popN0(df_flux0 %>% group_by(NIDT) %>% summarise(pop = sum(proba)))



# =============================== PART II: TESTING PESTIM PACKAGE ====================================

# ======================
#   postNO
# ======================

system.time({
  estim_t0 <- lapply(1:nrow(N_0), function(i) ({
  try(postN0(nMNO = N_0$pop_phone[i], nReg = N_0$pop_official[i],
         fu = list('unif', xMin = 0.3, xMax = 0.5),
         fv = list('unif', xMin = 100, xMax = N_0$pop_official[1]+1000),
         flambda = list('gamma', shape = 11, scale = 12),
         verbose = F))
}))
})

modifyError <- function(estim){
    nam <- c("postMean", "StDev", "CV", "postMedian", "Median_CI_LB", "Median_CI_UB", "MedianQuantileCV",
             "postMode", "Mode_CI_LB", "Mode_CI_UB", "ModeQuantileCV")
    if (class(estim)=="try-error") setNames(rep(NA,11), nam) else estim
}

df.estim_t0 <- do.call(rbind, lapply(1:length(estim_t0), function(i) modifyError(estim_t0[[i]]))) %>%
  tbl_df(.)

df.estim_t0


save(N_0, trans_N1N0, transitions_all, df.estim_t0, file = "D:/W3CRK9/Mes Documents/2018-04-25 test package européen/objects.RData")


# ======================
#   postNt
# ======================
# t = 1

p = vector(length = 24)
for(i in 1:24) {
  p[i] = which(N_0$NIDT[i] == rownames(trans_N1N0))
}
x<-trans_N1N0[p,p]

x[1,]<-14

nMNOmat = x
#trans_N1N0
nReg = N_0$pop_official

u0 <- rowSums(nMNOmat) /nReg
cv_u0 <- 0.15
fu <- lapply(u0, function(u){
  umin <- max(0, u - cv_u0 * u)
  umax <- min(1, u + cv_u0 * u)
  output <- list('unif', xMin = umin, xMax = umax)
  return(output)
})

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
  post_N1 <- postNt(nMNOmat/100,nReg/100, fu, fv, flambda, distNames, variation, scale = 100)
})


# ======================
#   postNtcondN0
# ======================

nMNOmat = x
#trans_N1N0
N0 = N_0$pop_official


# Names and parameters of priors for the transition probabilities
distNames <- rep('unif', length(N0))
variation <- rep(list(list(cv = 0.20)), length(N0))


postNtcondN0(N0, nMNOmat, distNames, variation)

########################### my code ###################

fu<-list()
fv<-list()
flambda<-list()

sc = 100

for(i in 1:24) {
  fu[[i]] <- list('unif', xMin = 0.3, xMax = 0.5)
  fv[[i]]<- list('unif', xMin = 1, xMax = N_0$pop_official[i]/sc+10)
  flambda[[i]]<-list('gamma', shape = 11, scale = 12)

}

system.time({
  estim_t0 <- postN0(nMNO = N_0$pop_phone/sc, nReg = N_0$pop_official/sc,
               fu = fu, fv = fv, flambda = flambda, scale=sc)
})



p = vector(length = 24)
for(i in 1:24) {
  p[i] = which(N_0$NIDT[i] == rownames(trans_N1N0))
}
x<-trans_N1N0[p,p]

x[1,]<-14

nMNOmat = x
#trans_N1N0
nReg = N_0$pop_official

u0 <- rowSums(nMNOmat) /nReg
cv_u0 <- 0.15
fu <- lapply(u0, function(u){
  umin <- max(0, u - cv_u0 * u)
  umax <- min(1, u + cv_u0 * u)
  output <- list('unif', xMin = umin, xMax = umax)
  return(output)
})

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
  post_N1 <- postNt(nMNOmat/100,nReg/100, fu, fv, flambda, distNames, variation, scale = 100)
})
#######################################
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
  y<-postNt(nMNOmat/100,nReg/100, fu, fv, flambda, distNames, variation, scale = 100)
  print(y)
}
})

