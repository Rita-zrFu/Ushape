## worked on Eli Lilly dataset
## one biomarker: BMI
# Y: only focus on the primary outcome (or all cause mortality?)
# Z: binary covariate smoke Yes->1; age 65+ -> 1
#restrict age>55; BMI ~ (18,35)
#use gaussian kernel

############## packages ################
rm(list = ls())
library(readxl)
library(tidyverse)
library(devtools)
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
library(plotly)

############## import data #################
setwd(".../real data")
realdata = as.data.frame(read_excel("simulated_datav5.xlsx"))
mydat <- realdata %>%
  filter(!is.na(BMIBASE), !is.na(Age), !is.na(Simulated_ID), !is.na(SMOKE)) 

mydat=mydat%>%filter(BMIBASE<35, BMIBASE>18)%>%
  mutate(taroutcm =ifelse(OUTCOME %in% c("CV Death(other than Stroke/MI)", 
                                         "Fatal MI", "Fatal Stroke", "Non-CV Death",
                                         "Non-Fatal Stroke", "Non-Fatal Silent MI"),
                          "Y","N"))#%>%filter(Age>55)

Y = delta = X = age = smoke = ID = NULL

for (i in unique(mydat$Simulated_ID)) {# some IDs are NA like 997
  temp_dat = subset(mydat, Simulated_ID == as.numeric(i))
  X[i] = unique(temp_dat$BMIBASE)
  age[i] = unique(temp_dat$Age)
  smoke[i] = unique(temp_dat$SMOKE)
  ID[i] = i
  #if(sum(temp_dat$PRIOUTCM=="Y")>0){
  if(sum(temp_dat$taroutcm == "Y")>0){
    Y[i] = tail(subset(temp_dat, taroutcm=="Y")$time,1)
    delta[i] = 1
  }else{
    Y[i] = subset(temp_dat, OUTCOME=="Last Observation for Patient")$time
    delta[i] = 0
  }
}
id_dat = data.frame(ID=ID, age = age, smoke = smoke, Y = Y, delta = delta, X = X)
id_dat = id_dat%>% 
  mutate(z1 = ifelse(id_dat$smoke=="Y",1,0),
         z2 = ifelse(id_dat$age>65, 1, 0)
         )%>% 
  filter(!if_all(.fns = is.na))
#nrow(id_dat)
censor_rate = sum(id_dat$delta)/length(id_dat$X)

new_z = cbind(id_dat$z1, id_dat$z2)
plot(id_dat$X, id_dat$Y, xlab = "BMI", ylab = "survival outcome", 
     main = "Scatter plot between survival outcome and biomarker")

############  Generate survival outcomes #######
set.seed(23)
x_pair = NULL
for (ix in 1:4){
  set.seed(ix)
  Zi <- unique(new_z)[ix,]
  submat = as.data.frame(subset(id_dat, z1==Zi[1]&z2==Zi[2]))
  Xi <- submat$X
  ni <- nrow(submat)
  Yi <- submat$Y
  deltai<- submat$delta
  censori<- submat$censor
  ti <- tail(unique(sort(Yi[submat$delta==1])), n=14)[1:10]  ## ti is the sequence of observed failure time
  xl_t = NULL
  xr_t = NULL
  
  for (j in 1:length(ti)) {
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
    mod = ksmooth(dat$x, dat$y, kernel = 'normal', bandwidth = 4)
    plot(mod,ylab="risk",xlab="BMI")
    
    x_min=min(mod$x); x_max=max(mod$x); y_min=min(mod$y); y_max=max(mod$y)
 
      f1 = function(x0){
        ksmooth(dat$x, dat$y, kernel = 'normal', bandwidth = 4, x.points = x0)$y
      }
      f2 = function(x0){
        -ksmooth(dat$x, dat$y, kernel = 'normal', bandwidth = 4, x.points = x0)$y
      }
      x_mid_low = optimize(f1, c(22, 35),tol = 1e-04)$minimum
    
    x_mid_high1 = optimize(f2, c(x_min,x_mid_low),tol = 1e-04)$minimum
    x_mid_high2 = optimize(f2, c(x_mid_low,x_max),tol = 1e-04)$minimum
    y_mid_low = f1(x_mid_low)
    y_mid_high1 = f1(x_mid_high1)
    y_mid_high2 = f1(x_mid_high2)

    y_right = f1(x_max)
    y_left = f1(x_min)
    a=seq(y_mid_low, min(y_mid_high1, y_mid_high2), 0.0005)
    jpeg(paste("D:/working/USR project/real data/ushape_age_smoke/",ix,"-",t,".jpeg",sep=''))
    plot(mod, xlab='x', ylab='y_hat', ylim=c(0,1))
    dev.off()
    
    for (i in 1:length(a)){
      f3 = function(x0){abs(f1(x0)-a[i])}
        xl[i] = optimize(f3, c(x_mid_high1, x_mid_low), tol = 1e-06)$minimum
        xr[i] = optimize(f3, c(x_mid_low, x_mid_high2), tol = 1e-06)$minimum
    }
    x_pair_one = cbind(ix,j,t,xl,xr,matrix(rep(Zi,length(xl)), ncol=2,byrow = TRUE))
    x_pair = rbind(x_pair, x_pair_one)
  }
}
est_par_list=NULL
set.seed(23)
for (j in 1:10) {
  group_t = x_pair[x_pair[,2]==j,]
  beta0_est=NULL;beta1_est=NULL;alpha2_est=NULL;alpha1_est=NULL;
  npair=nrow(unique(group_t))
  R=1000
  for (i in 1:R) {
    fit_lm = lm(new_pair[,5]~new_pair[,4]+new_pair[,6]+new_pair[,7]) # xr~xl+z
    beta0_est[i] = fit_lm$coefficients[1]/fit_lm$coefficients[2]
    beta1_est[i] = -1/fit_lm$coefficients[2]
    alpha1_est[i] = -fit_lm$coefficients[3]/fit_lm$coefficients[2]
    alpha2_est[i] = fit_lm$coefficients[4]/fit_lm$coefficients[2]
  }
  est_par = c(mean(beta0_est), mean(beta1_est),mean(alpha1_est,na.rm = TRUE), 
              mean(alpha2_est,na.rm = TRUE))
  est_par_list = rbind(est_par_list, est_par)
}


### get a pool of initial values
### In this simple linear regression model,use bootstrap to find the CI
############# find the optimized parameters  ############
est_par = c(mean(est_par_list[,1]),mean(est_par_list[,2],),
            mean(est_par_list[,3],na.rm=TRUE), mean(est_par_list[,4]))

cindex_fun <- function(theta) {
  htemp <- apply(id_dat[,7:9], 1,function(x){
    return(max(-x[1] + theta[3]*x[2], theta[1] + theta[2] * x[1] + theta[4]*x[3]))
  })#x2 is z1, x3 is z2
  y = id_dat[,4]
  d = id_dat[,6]
  cindex_value = -Cindex(Surv(y, d), htemp)
  return(cindex_value)
}
controlDE <- list(reltol=.0001, steptol=100, itermax = 1000, trace = 50)
fit_DE = DEoptim(fn = cindex_fun, lower = est_par-10, upper = est_par+10, control=controlDE)
mypar = fit_DE$optim$bestmem

cindex = -fit_DE$optim$bestval#
H_fun = pmax(-id_dat$X+mypar[3]*id_dat$z1, mypar[1] + mypar[2] * id_dat$X + mypar[4]*id_dat$z2 )
plot(id_dat$X, H_fun, xlab = "BMI", ylab = "estimated H function",
     main = "scatter plot of estimated H function and BMI",
     xlim = c(15,40),ylim = c(-40,-10))

###### estimating critical point and critical region ####

cpoint1_fun = function(par){
  return((-par[1]-par[4])/(1+par[2]))#0 1 nonsmoking old
}
cpoint2_fun = function(par){
  return((-par[1])/(1+par[2]))#z1=0 z2=0 non-smoking young 
}
cpoint3_fun = function(par){
  return((-par[1]+par[3])/(1+par[2]))#1 0 smoking young
}
cpoint4_fun = function(par){
  return((-par[1]+par[3]-par[4])/(1+par[2]))#1 1 smoking old 
}
cpoint1_fun(mypar)
cpoint2_fun(mypar)
cpoint3_fun(mypar)
cpoint4_fun(mypar)

cregin1_fun = function(par){
  return(c(25, (-25-par[1]-par[3]*1)/par[2]))
}
cregin2_fun = function(par){
  return(c(25, (-25-par[1])/par[2]))
}
cregin1_fun(mypar)
cregin2_fun(mypar)

est_cregin = c(25, (-25-mypar[1]-mypar[4]*1)/mypar[2])
est_cregin_2 = c(mypar[3]+25, (-25-mypar[1]-mypar[4]*1)/mypar[2])
est_cregin_3 = c(25, (-25-mypar[1])/mypar[2])
est_cregin_4 = c(mypar[3]+25, (-25-mypar[1])/mypar[2])

####### bootstrap variance #######
temp = list.files(".../results z1 smk z2 age")
mypar_list = NULL
cindex_list = NULL
cp1=NULL
cp2=NULL
cp3=NULL
cp4=NULL

for (i in temp) {
  boot_one = readRDS(paste0(".../results smk ag/",i))
  mypar_list = rbind(mypar_list, boot_one[1:4])
  cindex_list = c(cindex_list, boot_one[5])
  cp1 = c(cp1, boot_one[6])
  cp2 = c(cp2, boot_one[7])
  cp3 = c(cp3, boot_one[8])
  cp4 = c(cp4, boot_one[9])
}
boot_par = apply(mypar_list, 2, mean)
apply(mypar_list, 2, sd)
mean(cindex_list)

cp_m = cbind(cp1,cp2,cp3,cp4)
sum(cp_m<35&cp_m>18)/length(cp_m)
sum(cp1>18&cp1<35&cp2<35&cp3>18&cp4>18)/length(temp)
mean(cp1)
mean(cp2)
mean(cp3)
mean(cp4)
sd(cp1)
sd(cp3)
sd(cp2)
sd(cp4)
cpoint1_fun = function(par){
  return((-par[1]-par[4])/(1+par[2]))#0 z2=1 old group
}
cpoint2_fun = function(par){
  return((-par[1])/(1+par[2]))#z1=0 z2=0
}
cpoint3_fun = function(par){
  return((par[3]-par[1])/(1+par[2]))#z1=1 z2=0
}
cpoint4_fun = function(par){
  return((par[3]-par[1]-par[4])/(1+par[2]))#z1=1 z2=1
}
ci_cpoint = c(cpoint1_fun(mypar) - qnorm(0.975, mean=0, sd=1) * sd(cp1), 
              cpoint1_fun(mypar) + qnorm(0.975, mean=0, sd=1) * sd(cp1))
ci_cpoint2 = c(cpoint2_fun(mypar) - qnorm(0.975, mean=0, sd=1) * sd(cp2), 
               cpoint2_fun(mypar) + qnorm(0.975, mean=0, sd=1) * sd(cp2))
ci_cpoint3 = c(cpoint3_fun(mypar) - qnorm(0.975, mean=0, sd=1) * sd(cp3), 
               cpoint3_fun(mypar) + qnorm(0.975, mean=0, sd=1) * sd(cp3))
ci_cpoint4 = c(cpoint4_fun(mypar) - qnorm(0.975, mean=0, sd=1) * sd(cp4), 
               cpoint3_fun(mypar) + qnorm(0.975, mean=0, sd=1) * sd(cp4))

setdiff(1:1000,as.numeric((do.call(c,strsplit(do.call(rbind,strsplit(list.files(),"_"))[,3],".rds")))))

#### estimating Survival probability  ####
###### gaussian kernel ####
## use H and Y
gaussian_kernel <- function(x, h=2){ # h is bandwidth
  (1/(sqrt(2*pi)*h)) * exp(-x^2/(2*h^2))
}
get_s_list = function(group_dat){
  iso_result=NULL
  y_order = sort(unique(group_dat$Y[group_dat$delta == 1]))
  H_est = pmax(-group_dat$X + mypar[3]*group_dat$z1, 
               mypar[1] + mypar[2] * group_dat$X + mypar[4] * group_dat$z2)
  n = nrow(group_dat)
  num_times <- length(y_order)
  gaussian_values <- gaussian_kernel(matrix(rep(H_est, each = n), nrow=n) - 
                                         matrix(rep(H_est, n), nrow=n))
  S_list_gauss <- matrix(0, n, num_times)
  for(i in 1:n)  { # for every patient
    for (j in 1:num_times) { # at every time point
      t0 = y_order[j]
      yi_values <- unique(group_dat$Y[group_dat$Y <= t0 & group_dat$delta == 1])
      S_numerators <- sapply(yi_values, function(yi) {
        sum(gaussian_values[i, ] * (group_dat$Y == yi) * (group_dat$delta == 1))
      })
      S_denominators <- sapply(yi_values, function(yi) {
        sum(gaussian_values[i, ] * (group_dat$Y >= yi))
      })
      S <- prod(1 - S_numerators / S_denominators)
      S_list_gauss[i, j] <- S # each row is the result of each patient
    }
  }
  iso_result = apply(S_list_gauss, 2, function(col){
    iso_fit = isoreg(-H_est, col)
    return(list(yf=iso_fit$yf,
                h_ord=-iso_fit$x[iso_fit$ord],
                x_ord=group_dat$X[iso_fit$ord]))})
  return(iso_result)
}

####### plot figure 1 ######
for (i in 1:4){
  group_dat = id_dat[id_dat$z1==nique(new_z)[i,1] & id_dat$z2==unique(new_z)[i,2],]
  plot_t = sort(unique(group_dat$Y[group_dat$delta == 1]))
  s_list = get_s_list(group_dat)
  yf = do.call(cbind,lapply(s_list, function(list) return(list$yf)))
  x_order = s_list[[1]]$x_ord
  h_order = s_list[[1]]$h_ord
  km_fit = survfit(Surv(Y, delta)~1,data = group_dat)
  avg_S_gauss = colMeans(yf) 
  jpeg(paste0('.../z1 smk z2 age pic/km', i, '.jpeg'))
plot(km_fit, 
     main = "Kaplan-Meier Survival Curve", 
     xlab = "Time", 
     ylab = "Probability of Survival")
lines(km_fit$time[km_fit$n.event != 0], 
      avg_S_gauss, col = "blue", lwd = 2, lty = 1)
legend("bottomleft",legend=c("Gaussian Kernel smoothed"), col=c("blue"),
       lty=1,lwd=2)
dev.off()

time_points = km_fit$time[km_fit$n.event != 0]
l2_distance <- sum((km_fit$surv[km_fit$n.event != 0] - avg_S_gauss)^2 * c(diff(time_points),0))
z_sm = apply(yf,2,function(m){1-ksmooth(x_order,m,bandwidth = 2)$y})
jpeg(paste0('.../z1 smk z2 age pic/fun1_g', i, '.jpeg'))
persp3D(x=ksmooth(x_order,yf[,1],bandwidth = 2)$x,
        y=plot_t,
        z=z_sm,
        clim = c(0,1), clab = "risk", main="estimated risk", 
        xlab="biomarker", ylab="time", zlab="risk",
        scale=10, theta =  10, phi=10,
        ticktype="detailed",
        d=0.8, ltheta =10,lphi=90, r=10)
dev.off()

##### plot figure 2 #####
#z_sm = apply(yf2,2,function(m){ksmooth(x_order2,m,bandwidth = 2)$y})
x_sm = ksmooth(x_order, yf[,5],bandwidth = 2)$x
s_long = data.frame(t = rep(plot_t, times=nrow(yf)),
                    bmi = rep(x_sm, each=ncol(yf)),
                    st = c(t(z_sm)))%>%
  mutate(t=as.factor(t))%>%
  mutate(bmi=as.factor(bmi))
jpeg('.../z1 smk z2 age pic/fun2_g', i, '.jpeg',width = 800, height=1000,res=0.05)
ggplot(s_long, aes(x=t, y=bmi))+
  geom_tile(aes(fill=st))+
  theme(axis.text.x=element_text(angle=45,hjust=1, vjust=1))+
  theme_classic()+
  scale_fill_viridis_c(limits=c(0,1))
dev.off()
}



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

myframe = data.frame(X=id_dat$X, Y=id_dat$Y, z1=id_dat$z1,z2=id_dat$z2, delta = id_dat$delta)
mod1 = coxph(Surv(Y, delta)~X+z1+z2, data=myframe,x=TRUE, ties = "efron")
c_true = summary(mod1)$concordance[1] #0.7441342
cindex(mod1, formula=Surv(Y,delta)~X+z1+z2,data=myframe)


