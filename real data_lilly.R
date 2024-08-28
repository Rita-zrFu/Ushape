## worked on Eli Lilly dataset
## one biomarker at a time, start from HDL first, then use BMI
# Y: only focus on the primary outcome
# Z: binary covariate
# Z1: Age>65 -> 1, Z2: smoke Yes->1

############## packages ################
rm(list = ls())
library(readxl)
library(tidyverse)
library(devtools)
library(survcomp)
library(DEoptim)
library(rje)
library(MASS)
library(BB)
library(boot)
library(survival)
library(survminer)## Nice Summary of a Survival Curve
library(SurvMetrics)
library(dfoptim)
library(pec)
library(ggplot2)
library(plot3D)


############## import data #################
setwd("D:/working/USR project/real data")
realdata = as.data.frame(read_excel("simulated_datav5.xlsx"))

mydat = realdata[which(!is.na(realdata$BMIBASE) & !is.na(realdata$SMOKE) & !is.na(realdata$Simulated_ID)),]
mydat=mydat%>%filter(BMIBASE<35)%>%filter(BMIBASE>18)%>%filter(Age>55)
Y = NULL; censor = NULL; delta = NULL; X = NULL; age = NULL; smoke = NULL; ID = NULL; gender=NULL

for (i in unique(mydat$Simulated_ID)) {# some IDs are NA like 997
  temp_dat = subset(mydat, Simulated_ID == as.numeric(i))
  X[i] = unique(temp_dat$BMIBASE)
  age[i] = unique(temp_dat$Age)
  #gender[i] = unique(temp_dat$Gender)
  smoke[i] = unique(temp_dat$SMOKE)
  ID[i] = i
  #if(sum(temp_dat$PRIOUTCM=="Y")>0){
  if("CV Death(other than Stroke/MI)"%in% temp_dat$OUTCOME|"Fatal MI" %in%temp_dat$OUTCOME
     |"Fatal Stroke"%in%temp_dat$OUTCOME|"Non-CV Death" %in% temp_dat$OUTCOME
     |"Non-Fatal Stroke"%in% temp_dat$OUTCOME|"Non-Fatal Silent MI"%in% temp_dat$OUTCOME 
     #|"Coronary Revascularization"%in% temp_dat$OUTCOME
     ){
    Y[i] = tail(temp_dat$time,1)
    # Y is the first time observed primary outcome
    censor[i] = subset(temp_dat, OUTCOME=="Last Observation for Patient")$time
    delta[i] = 1
  }else{
    Y[i] = censor[i] = subset(temp_dat, OUTCOME=="Last Observation for Patient")$time
    delta[i] = 0
  }
}
# do not identify censor like this, do not use censor below
id_dat = data.frame(ID=ID, smoke = smoke, age = age, Y = Y, censor = censor, 
                    delta = delta, X = X)
id_dat = id_dat%>% #z1 = ifelse(id_dat$age>65, 1, 0),
  mutate(
         z2 = ifelse(id_dat$smoke=="Y",1,0)  )%>% #  
  filter(!if_all(.fns = is.na))
#nrow(id_dat)
censor_rate = sum(id_dat$delta)/length(id_dat$X)

#new_z = cbind(id_dat$z1, id_dat$z2)
new_z=id_dat$z2
plot(id_dat$X, id_dat$Y, xlab = "BMI", ylab = "survival outcome", 
     main = "Scatter plot between survival outcome and biomarker")

############  Generate survival outcomes #######
set.seed(23)
x_pair = NULL
for (ix in 1:2){
  set.seed(ix)
  Zi <- unique(new_z)[ix]
  submat = as.data.frame(subset(id_dat, z2==Zi))#z1==Zi[1]&
  Xi <- submat$X
  ni <- nrow(submat)
  Yi <- submat$Y
  deltai<- submat$delta
  censori<- submat$censor
  ti <- tail(unique(sort(Yi[submat$delta==1])), n=14)[1:10]  ## ti is the sequence of observed failure time
  xl_t = NULL
  xr_t = NULL
  # i=0
  # num_event = NULL
  # num_censor = NULL
  # for (k in unique(sort(Yi[submat$delta==1]))) {
  #   y_bin = sapply(1:ni, function(i) ifelse(Yi[i]>=k,0,ifelse(deltai[i]==1,1,NA)))
  #   i=i+1
  #   num_event[i] = sum(y_bin==1, na.rm = TRUE)
  #   num_censor[i] = sum(y_bin==0, na.rm = TRUE)
  #   }
  
  for (j in 1:length(ti)) {
    #j=6
    t=ti[j]
    y_bin = sapply(1:ni, function(i) ifelse(Yi[i]>=t,0,ifelse(deltai[i]==1,1,NA)))
    
    ####### kernel smooth 
    set.seed(floor(t))
    x = Xi
    y = y_bin
    dat = data.frame(x, y)
    dat<-dat[complete.cases(dat),]
    plot(dat$x,dat$y)

    ### find xl xr pairs
    xl=NULL; xr=NULL
    mod = ksmooth(dat$x, dat$y, kernel = 'normal', bandwidth = 3)
    plot(mod,ylab="risk",xlab="BMI")
    
    x_min=min(mod$x); x_max=max(mod$x); y_min=min(mod$y); y_max=max(mod$y)
    if(ix==1){
      f1 = function(x0){
        ksmooth(dat$x, dat$y, kernel = 'normal', bandwidth = 5, x.points = x0)$y
      }
      f2 = function(x0){
        -ksmooth(dat$x, dat$y, kernel = 'normal', bandwidth = 5, x.points = x0)$y
      }
      x_mid_low = optimize(f1, c(25, 32),tol = 1e-04)$minimum
    }else{
      f1 = function(x0){
        ksmooth(dat$x, dat$y, kernel = 'normal', bandwidth = 3, x.points = x0)$y
      }
      f2 = function(x0){
        -ksmooth(dat$x, dat$y, kernel = 'normal', bandwidth = 3, x.points = x0)$y
      }
      x_mid_low = optimize(f1, c(25, 32),tol = 1e-04)$minimum
    }
    

    x_mid_high1 = optimize(f2, c(x_min,x_mid_low),tol = 1e-04)$minimum
    x_mid_high2 = optimize(f2, c(x_mid_low,x_max),tol = 1e-04)$minimum
    y_mid_low = f1(x_mid_low)
    y_mid_high1 = f1(x_mid_high1)
    y_mid_high2 = f1(x_mid_high2)

    y_right = f1(x_max)
    y_left = f1(x_min)
    a=seq(y_mid_low, min(y_mid_high1, y_mid_high2), 0.0005)
    jpeg(paste("D:/working/USR project/real data/u_smoke/",ix,"-",t,".jpeg",sep=''))
    plot(mod, xlab='x', ylab='y_hat', ylim=c(0,1))
    dev.off()
    
    for (i in 1:length(a)){
      f3 = function(x0){abs(f1(x0)-a[i])}
      if(ix==1){
        xl[i] = optimize(f3, c(x_mid_high1, x_mid_low), tol = 1e-06)$minimum
        xr[i] = optimize(f3, c(x_mid_low, x_mid_high2), tol = 1e-06)$minimum
      }else{
        xl[i] = optimize(f3, c(x_mid_high1, x_mid_low), tol = 1e-06)$minimum
        xr[i] = optimize(f3, c(x_mid_low, x_mid_high2), tol = 1e-06)$minimum
      }
    }
    x_pair_one = cbind(ix,j,t,xl,xr,rep(Zi,length(xl)))
    x_pair = rbind(x_pair, x_pair_one)
  }
}
est_par_list=NULL
set.seed(23)
for (j in 1:10) {
  group_t = x_pair[x_pair[,2]==j & (x_pair[,1]==1 | x_pair[,1]==2|x_pair[,1]==3),]
  beta0_est=NULL;beta1_est=NULL;alpha2_est=NULL;alpha1_est=NULL;
  npair=nrow(unique(group_t))
  R=1000
  for (i in 1:R) {
    new_pair = group_t[sample(npair,replace=TRUE),]
    fit_lm = lm(new_pair[,5]~new_pair[,4]+new_pair[,6])#+new_pair[,7]) # xr~xl+z
    beta0_est[i] = fit_lm$coefficients[1]/fit_lm$coefficients[2]
    beta1_est[i] = -1/fit_lm$coefficients[2]
    #alpha1_est[i] = -fit_lm$coefficients[3]/fit_lm$coefficients[2]
    alpha2_est[i] = fit_lm$coefficients[3]/fit_lm$coefficients[2]
  }
  est_par = c(mean(beta0_est), mean(beta1_est), #mean(alpha1_est,na.rm = TRUE), 
              mean(alpha2_est,na.rm = TRUE))
  est_par_list = rbind(est_par_list, est_par)
}


### get a pool of initial values
### In this simple linear regression model,use bootstrap to find the CI
saveRDS(est_par,"D:/working/USR project/real data/est_par.rds")
load('D:/working/USR project/real data/214.rdata')
############# find the optimized parameters  ############
est_par = c(mean(est_par_list[,1]),mean(est_par_list[,2],),
            mean(est_par_list[,3],na.rm=TRUE))
penalty <- function(theta) {
  if (theta[2]<=1) {
    return(1e6) # large penalty value if the constraint is violated
  }
  return(0) # no penalty if the constraint is satisfied
}
cindex_fun <- function(theta) {
  # htemp <- apply(id_dat[,7:9], 1,function(x){
  #   return(max(-x[1] + theta[3]*x[2], theta[1] + theta[2] * x[1] + theta[4]*x[3]))
  # })#x2 is z1, x3 is z2
  htemp <- apply(id_dat[,7:8], 1,function(x){
       return(max(-x[1], theta[1] + theta[2] * x[1] + theta[3]*x[2]))
     })#x2 is z2
  y = id_dat[,4]
  d = id_dat[,6]
  cindex_value = -Cindex(Surv(y, d), htemp)
  return(cindex_value+penalty(theta))
}
par_DE = c("cindex_fun", "id_dat", "est_par", "penalty")
controlDE <- list(reltol=.0001, steptol=100, itermax = 1000, trace = 50,
                  parallelType = 1, parVar = par_DE)
fit_DE = DEoptim(fn = cindex_fun, lower = est_par-10, upper = est_par+10, control=controlDE)
mypar = fit_DE$optim$bestmem

cindex = -fit_DE$optim$bestval
H_fun = pmax(-id_dat$X, mypar[1] + mypar[2] * id_dat$X + mypar[3]*id_dat$z2 )
plot(id_dat$X, H_fun, xlab = "BMI", ylab = "estimated H function",
     main = "scatter plot of estimated H function and BMI")

  #### calculate critical point ####
# given z1=0,z2=1,estimate change/critical point, critical region
est_cpoint = (-mypar[1]-mypar[3])/(1+mypar[2])#0 1
#est_cpoint_2 = (-mypar[1]+mypar[3]-mypar[4])/(1+mypar[2])#1 1
est_cpoint_2 = (-mypar[1])/(1+mypar[2])#0 0
#est_cpoint_4 = (-mypar[1]+mypar[3])/(1+mypar[2])#1 0

saveRDS(fit_DE, "D:/working/USR project/real data/fit_DE.rds")
saveRDS(mypar,"D:/working/USR project/real data/mypar.rds")
mypar = readRDS("D:/working/USR project/real data/mypar.rds")

est_cregin = c(25, (-25-mypar[1]-mypar[4]*1)/mypar[2])
#est_cregin_2 = c(mypar[3]+25, (-25-mypar[1]-mypar[4]*1)/mypar[2])
est_cregin_3 = c(25, (-25-mypar[1])/mypar[2])
#est_cregin_4 = c(mypar[3]+25, (-25-mypar[1])/mypar[2])


########## bootstrap to get variance  #########
set.seed(23)

ind_big= lapply(1:100, function(i){
  set.seed(i)
  sample(nrow(id_dat),replace = TRUE)
})#sample index
i=0
mylist=lapply(ind_big, function(ind){
  i=i+1
  print(i)
  penalty <- function(theta) {
    if (theta[1]>=0|theta[2]<=1|
        17*theta[1]+35*theta[3]>=0) {
      return(1e6) # large penalty value if the constraint is violated
    }
    return(0) # no penalty if the constraint is satisfied
  }
  
  cindex_my = function(theta){
    mat = id_dat[ind, ]
    htemp <- apply(mat[,7:8], 1, function(x){
      return(max(-x[1], theta[1] + theta[2] * x[1] + theta[3] * x[2]))
    })
    d <- mat[,6]
    y <- mat[,4]
    return(-Cindex(Surv(y, d), htemp)+penalty(theta))
  }
  fit_nmk = nmkb(par = mypar, fn = cindex_my,upper = c(0,4,20))
  temp_cpoint = c((-fit_nmk$par[1]-fit_nmk$par[3])/(1+fit_nmk$par[2]),#1
                  (-fit_nmk$par[1])/(1+fit_nmk$par[2])#0 
  )
  return(c(fit_nmk$par, fit_nmk$value, fit_nmk$convergence, temp_cpoint))
})

boot_result = do.call(rbind, mylist)

boot1 = boot_result[which(boot_result[,6]>18),]
boot2 = boot1[which(boot1[,7]<35),]
boot3 = boot2[which(boot2[,6]<35),]
boot4 = boot3[which(boot3[,7]>18),]
saveRDS(boot_result,"D:/working/USR project/real data/boot_result1.rds")
saveRDS(boot4,"D:/working/USR project/real data/boot_result1_flt.rds")

converg = mean(boot_result[,5])
sd_est = apply(boot4, 2, sd)
ci_par = c(mypar - qnorm(0.975, mean=0, sd=1) * sd_est[1:3], 
           mypar + qnorm(0.975, mean=0, sd=1) * sd_est[1:3])
ci_cindex = c(cindex - qnorm(0.975, mean=0, sd=1) * sd_est[4], 
              cindex + qnorm(0.975, mean=0, sd=1) * sd_est[4])
ci_cpoint = c(est_cpoint - qnorm(0.975, mean=0, sd=1) * sd_est[6], 
              est_cpoint + qnorm(0.975, mean=0, sd=1) * sd_est[6])
ci_cpoint2 = c(est_cpoint_2 - qnorm(0.975, mean=0, sd=1) * sd_est[7], 
               est_cpoint_2 + qnorm(0.975, mean=0, sd=1) * sd_est[7])
# ci_cpoint3 = c(est_cpoint_3 - qnorm(0.975, mean=0, sd=1) * sd_est[9], 
#                est_cpoint_3 + qnorm(0.975, mean=0, sd=1) * sd_est[9])
# ci_cpoint4 = c(est_cpoint_4 - qnorm(0.975, mean=0, sd=1) * sd_est[10], 
#                est_cpoint_4 + qnorm(0.975, mean=0, sd=1) * sd_est[10])
#CI_index_par = theta_true>ci_par[1:4] & theta_true<ci_par[5:8]
#CI_index_cp = cp_true>ci_cpoint[1] & cp_true<ci_cpoint[2]
#CI_index_cr_up = cr_true[2]>ci_cr_up[1] & cr_true[2]<ci_cr_up[2]
#CI_index_cr_lo = cr_true[1]>ci_cr_lo[1] & cr_true[1]<ci_cr_lo[2]

#### estimating St  ####
###### knn kernel ####
## use knn near H0
g1 = id_dat[id_dat$z2==0,]
g2 = id_dat[id_dat$z2==1,]

plot_t1 = sort(unique(g1$Y[g1$delta == 1]))
plot_t2 = sort(unique(g2$Y[g2$delta == 1]))

get_s_list = function(group_dat){
  iso_result=NULL
  y_order = sort(unique(group_dat$Y[group_dat$delta == 1]))
  H_est = pmax(-group_dat$X, 
               mypar[1] + mypar[2] * group_dat$X + mypar[3] * group_dat$z2)
  n = nrow(group_dat)
  num_times <- length(y_order)
  S_list_knn <- matrix(0, n, num_times)
  for(i in 1:n) { # for every patient
    if(i%%10==0) print(floor(i/10))
    for (j in 1:num_times) { # for every time point
      t0 = y_order[j]
      H0 = H_est[i]      
      H0_ind = which(sort(H_est-H0)==0)
      sort_ind = sort(H_est-H0,index.return=TRUE)$ix
      if(H0_ind[1]<=5){
        neighbor = sort_ind[1:10]
      }else if(tail(H0_ind,1)>=n-5){
        neighbor = tail(sort_ind, 10)
      }else{
        neighbor = sort_ind[(H0_ind[1]-5):(tail(H0_ind,1)+5)]
      }
      yj = group_dat$Y[neighbor]
      dj = group_dat$delta[neighbor]
      relevant_y_order = y_order[y_order <= t0]
      S = prod(1 - sapply(relevant_y_order, function(yi) sum(yj == yi & dj == 1) / sum(yj >= yi)))
      S_list_knn[i, j] <- S
    }
  }
    #na_row = unique(which(is.na(S_list_knn),arr.ind = TRUE)[,1])
    #S_list_knn=na.omit(S_list_knn)
  slist_rm_col = S_list_knn[, apply(S_list_knn, 2, function(y) all(!is.na(y)))]
    #H_new = H_est[-na_row]
    #x_new = group_dat$X[-na_row]

  # iso_result = apply(S_list_knn, 2, function(col){
  #   iso_fit = isoreg(H_new, col)#返回的yf是排序了的，返回的isofun$x[isofun$ord]是排序的H，不能再用原来的H了
  #   return(list(yf=iso_fit$yf,
  #               h_ord=iso_fit$x[iso_fit$ord],
  #               x_ord=x_new[iso_fit$ord]))})

  iso_result = apply(slist_rm_col, 2, function(col){
    iso_fit = isoreg(-H_est, col)#返回的yf是排序了的，
    #返回的isofun$x[isofun$ord]是排序的H，不能再用原来的H了
    return(list(yf=iso_fit$yf,
                #列都是yorder，只是行在排序
                h_ord=-iso_fit$x[iso_fit$ord],
                x_ord=group_dat$X[iso_fit$ord]))})
  return(iso_result)
}


s_list1 = get_s_list(g1)
yf1 = do.call(cbind,lapply(s_list1, function(list) return(list$yf)))
x_order1 = s_list1[[1]]$x_ord
#h_order1 = s_list1[[1]]$h_ord

s_list2 = get_s_list(g2)
yf2 = do.call(cbind,lapply(s_list2, function(list) return(list$yf)))
x_order2 = s_list2[[1]]$x_ord

km_fit = survfit(Surv(Y, delta)~1,data = g1)
avg_S_knn = colMeans(yf1)
#S_list_knn[is.na(S_list_knn)]=0
plot(km_fit, 
     main = "Kaplan-Meier Survival Curve", 
     xlab = "Time", 
     ylab = "Probability of Survival")
lines(km_fit$time[km_fit$n.event != 0][1:length(avg_S_knn)], 
      avg_S_knn, col = "blue", lwd = 2, lty = 1)


z_sm = apply(yf1,2,function(m){ksmooth(x_order1,m,bandwidth = 5)$y})
#yf1对x的图是0-1图，z_sm就是smooth过后的yf1
#ksmooth先把x增序排列了，然后smooth的y！！！我们的xorder1是降序排列的！！
#所以画出来的ushape是相反的
#所以我们在画图的时候要用ksmooth产生的x,或者用reverse的ksmooth产生的y
jpeg('D:/working/USR project/real data/g1fun_2.jpeg')
persp3D(x=ksmooth(x_order1,yf1[,1],bandwidth = 5)$x,
        y=plot_t1[1:length(avg_S_knn)],#变成[1:6]
        z=z_sm[,1:length(avg_S_knn)],
        main="estimated risk",
        xlab="biomarker",
        ylab="time",
        zlab="survival probability",
        clab = "survival probability",
        scale=10,
        theta =  0,
        phi=300,
        ticktype="detailed",
        d=0.8,
        ltheta =10,lphi=90,
        r=10)
dev.off()


##### plot figure 2#####
  get_s_list2 = function(group_dat){
    y_order = sort(unique(group_dat$Y[group_dat$delta == 1]))
    xplace = sort(group_dat$X,index.return=TRUE)$ix
    H_est = pmax(-group_dat$X, 
                 mypar[1] + mypar[2] * group_dat$X + mypar[3] * group_dat$z2)
    n = nrow(group_dat)
    num_times <- length(y_order)
    S_list_knn <- matrix(0, n, num_times)
    for(i in 1:n) { # for every patient
      if(i%%10==0) print(floor(i/10))
      for (j in 1:num_times) { # for every time point
        t0 = y_order[j]
        H0 = H_est[xplace[i]]
        H0_ind = which(sort(H_est-H0)==0)
        sort_ind = sort(H_est-H0,index.return=TRUE)$ix
        if(H0_ind[1]<=5){
          neighbor = sort_ind[1:10]
        }else if(tail(H0_ind,1)>=n-5){
          neighbor = tail(sort_ind, 10)
        }else{
          neighbor = sort_ind[(H0_ind[1]-5):(tail(H0_ind,1)+5)]
        }
        yj = group_dat$Y[neighbor]
        dj = group_dat$delta[neighbor]
        relevant_y_order = y_order[y_order <= t0]
        S = prod(1 - sapply(relevant_y_order, function(yi) sum(yj == yi & dj == 1) / sum(yj >= yi)))
        S_list_knn[i, j] <- S
      }
    }
    return(S_list_knn)
  }
s_order1 = get_s_list2(g1)
library(reshape2)
library(ggplot2)  


z_sm = apply(yf1,2,function(m){ksmooth(x_order1,m,bandwidth = 2)$y})
s_long1 = data.frame(t = rep(plot_t1,times=nrow(yf1)),
                     bmi = rep(ksmooth(x_order1,yf1[,10],bandwidth = 2)$x,
                               each=ncol(yf1)),
                     st = c(t(z_sm)))
# length(plot_t2)
# length(unique(plot_t2))
# length(ksmooth(x_order2,yf2[,10],bandwidth = 1.5)$x)
# length(unique(ksmooth(x_order2,yf2[,10],bandwidth = 1.5)$x))
slong2 = na.omit(s_long1)#%>%
  #filter(st<0.9&st>0.7)

#slong2=na.omit(slong2)
slong2$t=as.factor(slong2$t)
slong2$bmi=as.factor(slong2$bmi)
jpeg('D:/working/USR project/real data/fun3_new.jpeg',width = 800, height=1000,res=0.05)
ggplot(slong2, aes(x=t, y=bmi))+
  geom_tile(aes(fill=st))+
  theme(axis.text.x=element_text(angle=45,hjust=1, vjust=1))+
  theme_classic()+
  scale_fill_viridis_c()
dev.off()

##### pava ####
S_list_knn[is.na(S_list_knn)]=0
iso_fit = isoreg(id_dat$X, S_list_knn[,2])
iso_result = apply(S_list_knn,2,function(i){
  iso_fit = isoreg(id_dat$X, i)
  return(iso_fit$x[iso_fit$ord])})

plot(id_dat$X, S_list1[,2], pch = 4,xlab = "BMI in testing dataset", 
     ylab = "estimated survival probability",
     main = "scatter plot of estimated survival probability and BMI")
points(iso_fit$x[iso_fit$ord], iso_fit$yf, pch = 16, col = "blue")
legend("bottomright",                                   
       legend=c("before isotonation","after isotonation"),       
       col=c("black","blue"),               
       lty=1,lwd=2,pch = c(4,16), cex = 1,
       bty = "n",
       pt.cex = 1, pt.lwd = 2, seg.len = 1,
       y.intersp = 0.2,
       x.intersp = 0.5,
       text.width = 5) 


plot(km_fit, main = "Kaplan-Meier Survival Curve", xlab = "Time", ylab = "Probability of Survival")
lines(km_fit$time, avg_S_knn, col = "red", lwd = 2, lty = 1)
lines(km_fit$time[km_fit$n.event != 0], avg_S_knn, col = "blue", lwd = 2, lty = 1)
#lines(km_fit$time, avg_S_gauss2, col = "green", lwd = 2, lty = 1)

legend("bottomright", cex=1, bty="n",c("KNN kernel"), col=c("blue"), 
       lty=1, lwd=c(2), x.intersp = c(0.5), y.intersp = c(0.5),
       xjust = 0, yjust=0)

legend("bottomright", cex=1, bty="n",c("gaussian", "knn"), col=c("red","blue"), 
       lty=1:1, lwd=c(2,2), x.intersp = c(0.5,0.5), y.intersp = c(0.5,0.5),
       xjust = 0, yjust=0)


#### real C-index ####
# c-index calculated by other packages

myframe = data.frame(X=id_dat$X, Y=id_dat$Y, z2=id_dat$z2, delta = id_dat$delta)
mod1 = coxph(Surv(Y, delta)~X+z2, data=myframe,x=TRUE, ties = "efron")
c_true = summary(mod1)$concordance[1] #0.7441342
cindex(mod1, formula=Surv(Y,delta)~X+z2,data=myframe)
