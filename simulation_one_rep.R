# one replication simulation
# on USR manuscript simulation plan version
# included sd/se, CI, point estimate and inference
# and how to calculate C-inex (use Cindex function correctly)
# scenario: beta1=1, ei follows a standard normal distribution, n=200
# no G(x) transformation
# beta1 <- 1

rm(list = ls())
########### packages ##############
library(tidyverse)
library(devtools)
library(DEoptim)
library(rje)
library(MASS)
library(BB)
library(boot)
library(survival)
library(SurvMetrics)
library(parallel)
library(dfoptim)

######### parameters settings ############
# variables ei, n, seed need to be assigned values
# ei <- "norm" #, "lnorm", "exp"
# n <- 200 #, 500
# seed <- 1 #, 2, ..., 1000
set.seed(seed)

beta0 = 0
beta1 = 1
alpha1 = 10
alpha2 = -10

z1 <- rnorm(n, 0, 1)
z2 <- rbinom(n, 1, 0.5)

x <- runif(n, -50, 50)
t1 <- -x + alpha1 * z1
t2 <- beta0 + beta1 * x + alpha2 * z2

if (ei=="norm"){
  # epsilon follows the standard normal distribution
  epsilon = rnorm(n, 0, 1)
}else if(ei=="lnorm"){
  # epsilon follows the log standard normal distribution
  epsilon = rlnorm(n, 0, 1)
}else{
  # epsilon follows the exponential distribution
  epsilon = rexp(n, 1)
}

## control the censor rate within 10%-20%
censor = runif(n, 8.5, 9.4)

t = 10*plogis(pmax(t1,t2) + epsilon, scale=20)

y <- pmin(t, censor)
delta <- sapply(1:n, function(i) {
  return(as.numeric(t[i] < censor[i]))
})
censor_rate = sum(delta)/n

########## find optimal parameters  #############
#### get initial values ####
mymat = cbind(x, z1, z2, y, delta)
nr = nrow(mymat)

theta_true = c(beta0, log(beta1), alpha1, alpha2)

############ point estimate ###########
#### use DEoptim to find a reasonable initial value for optim
cindex_DE = function(theta){
  htemp <- apply(mymat[,1:4], 1, function(x){
    return(max(-x[1] + theta[3] * x[2], theta[1] + exp(theta[2]) * x[1] + theta[4] * x[3]))
  })
  d <- mymat[,5]
  y <- mymat[,4]
  return(-Cindex(Surv(y, d), htemp))
}

par_DE = c("cindex_DE", "mymat", "theta_true")
controlDE <- list(reltol=.000001, steptol=100, itermax = 1000, trace = 50,
                  parallelType = 1, parVar = par_DE)
fit_DE = DEoptim(fn = cindex_DE, lower = theta_true-5, upper = theta_true+5, control=controlDE)

est_par = fit_DE$optim$bestmem
est_cindex = -fit_DE$optim$bestval

## given z1=0.5,z2=1,estimate change/critical point, critical region
est_cpoint = (-est_par[1]+est_par[3]*0.5-est_par[4]*1)/(1+exp(est_par[2]))
est_cregin = c(est_par[3]*0.5-5, (5-est_par[1]-est_par[4]*1)/exp(est_par[2]))

cp_true = (10*0.5+10)/(1+1)
cr_true = c(5-5, (5+10)/1)
bias=est_par-theta_true
bias_cp = est_cpoint-cp_true
bias_cr = est_cregin-cr_true

############# inference ##############
nboot = 1000

mat_big= lapply(1:nboot, function(i){
  set.seed(i)
  mymat[sample(n,replace = TRUE),]
})

mylist=lapply(mat_big, function(mat){
  cindex_my = function(theta){
    htemp <- apply(mat[,1:4], 1, function(x){
      return(max(-x[1] + theta[3] * x[2], theta[1] + exp(theta[2]) * x[1] + theta[4] * x[3]))
    })
    d <- mat[,5]
    y <- mat[,4]
    return(-Cindex(Surv(y, d), htemp))
  }
  fit_nmk = nmk(par = est_par, fn = cindex_my)
  temp_cpoint = (-fit_nmk$par[1]+fit_nmk$par[3]*0.5-fit_nmk$par[4]*1)/(1+exp(fit_nmk$par[2]))
  temp_cregin = c(fit_nmk$par[3]*0.5-5, (5-fit_nmk$par[1]-fit_nmk$par[4]*1)/exp(fit_nmk$par[2]))
  return(c(fit_nmk$par, fit_nmk$value, fit_nmk$convergence, temp_cpoint, temp_cregin))
})

boot_result = do.call(rbind, mylist)
converg = mean(boot_result[,6])
sd_est = apply(boot_result, 2, sd)
ci_par = c(est_par - qnorm(0.975, mean=0, sd=1) * sd_est[1:4], 
           est_par + qnorm(0.975, mean=0, sd=1) * sd_est[1:4])
ci_cindex = c(est_cindex - qnorm(0.975, mean=0, sd=1) * sd_est[5], 
              est_cindex + qnorm(0.975, mean=0, sd=1) * sd_est[5])
ci_cpoint = c(est_cpoint - qnorm(0.975, mean=0, sd=1) * sd_est[7], 
              est_cpoint + qnorm(0.975, mean=0, sd=1) * sd_est[7])
ci_cr_up = c(est_cregin[2] - qnorm(0.975, mean=0, sd=1) * sd_est[9], 
             est_cregin[2] + qnorm(0.975, mean=0, sd=1) * sd_est[9])
ci_cr_lo = c(est_cregin[1] - qnorm(0.975, mean=0, sd=1) * sd_est[8], 
             est_cregin[1] + qnorm(0.975, mean=0, sd=1) * sd_est[8])

CI_index_par = theta_true>ci_par[1:4] & theta_true<ci_par[5:8]
CI_index_cp = cp_true>ci_cpoint[1] & cp_true<ci_cpoint[2]
CI_index_cr_up = cr_true[2]>ci_cr_up[1] & cr_true[2]<ci_cr_up[2]
CI_index_cr_lo = cr_true[1]>ci_cr_lo[1] & cr_true[1]<ci_cr_lo[2]

########### estimating St  ############
## St(5.1|0, c(0, 0)) ##calculated by R function?
t0 = 5.1
y_order = sort(unique(y[delta == 1]))
H_est = pmax(-x + est_par[3] * z1, 
             est_par[1] + exp(est_par[2]) * x + est_par[4] * z2)
H0 = pmax(0, est_par[1])

sort_x = sort(H_est-H0,index.return=TRUE)$x
sort_ind = sort(H_est-H0,index.return=TRUE)$ix
H0_ind = which(abs(sort_x)==min(abs(H_est-H0)))

if(H0_ind[1]<=5){
  neighbor = sort_ind[1:10]
}else if(tail(H0_ind,1)>=length(H_est)-5){
  neighbor = tail(sort_ind, 10)
}else{
  neighbor = sort_ind[(H0_ind[1]-5):(tail(H0_ind,1)+5)]
}
yj = y[neighbor]
dj = delta[neighbor]

relevant_y_order = y_order[y_order <= t0]

S = prod(1 - sapply(relevant_y_order, function(yi) sum(yj == yi & dj == 1) / sum(yj >= yi)))

main.dir <- "/users/zfu1/tree_project/final/"
pilot.dir <- "sim_fin"
if(file.exists(pilot.dir)) {
  setwd(file.path(main.dir, pilot.dir))
} else {
  dir.create(file.path(main.dir, pilot.dir))
  setwd(file.path(main.dir, pilot.dir))
}

sub.dir <- paste0("beta_", beta1*10, "_e_", name_ei, "_n_", n)
if(file.exists(sub.dir)) {
  setwd(file.path(main.dir, pilot.dir, sub.dir))
} else {
  dir.create(file.path(main.dir, pilot.dir, sub.dir))
  setwd(file.path(main.dir, pilot.dir, sub.dir))
}
save.image(file = paste0(sub.dir, "-k", seed, ".rda"))
setwd(file.path(main.dir))

