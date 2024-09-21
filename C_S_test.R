library(survival)
##### S and C on test dataset #####
mydir = ".../tree_project/sim_fin/"
setwd(mydir)
allfiles = list.files()
args = commandArgs(TRUE)
file_num = as.numeric(args[1])

senario  = allfiles[file_num]
subdir = paste0(mydir, "/", senario)
setwd(subdir)
print(subdir)
array_num <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
  
load(list.files()[array_num])
set.seed(seed)
print(seed)
n=10000
z1 <- rnorm(n, 0, 1)
z2 <- rbinom(n, 1, 0.5)
x <- runif(n, -50, 50)
t1 <- -x + alpha1 * z1
t2 <- beta0 + beta1 * x + alpha2 * z2
if (ei==1){
  epsilon = rnorm(n, 0, 1)
}else if(ei==2){
  epsilon = rlnorm(n, 0, 1)
}else{
  epsilon = rexp(n, 1)
}
    
if (beta1==1){
  censor = runif(n, 8.5, 9.4)
}else if(beta1==1.5){
  censor = runif(n, 9, 10)
}else{
  censor = runif(n, 8.5, 9.5)
}
    
t = 10*plogis(pmax(t1,t2) + epsilon, scale=20)
y <- pmin(t, censor)
delta <- sapply(1:n, function(i) {
  return(as.numeric(t[i] < censor[i]))
})
    
htemp <- pmax(-x + est_par[3] * z1, est_par[1] + exp(est_par[2]) * x + est_par[4] * z2)
c_test =  Cindex(Surv(y, delta), htemp)
    
## St(5|1, c(0, 0.5)) ##
t0 = 5
y_order = sort(unique(y[delta == 1]))
H_est = pmax(-x + est_par[3] * z1, 
              est_par[1] + exp(est_par[2]) * x + est_par[4] * z2)
H0 = pmax(-1, est_par[1] + exp(est_par[2])  + est_par[4]*0.5)
    
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
s_test = prod(1 - sapply(relevant_y_order, function(yi) sum(yj == yi & dj == 1) / sum(yj >= yi)))
    
main.dir <- ".../tree_project/"
pilot.dir <- "sim_fin_test_c_s"

if(file.exists(pilot.dir)) {
  setwd(file.path(main.dir, pilot.dir))
} else {
  dir.create(file.path(main.dir, pilot.dir))
  setwd(file.path(main.dir, pilot.dir))
}
    
sub.dir <- paste0("_e_", name_ei, "_n_", n)
    
if(file.exists(sub.dir)) {
  setwd(file.path(main.dir, pilot.dir, sub.dir))
} else {
  dir.create(file.path(main.dir, pilot.dir, sub.dir))
  setwd(file.path(main.dir, pilot.dir, sub.dir))
}
    
save(list = c("c_test", "s_test"), file = paste0(sub.dir, "-k", seed, "_", file_num, ".RData"))
    
  

