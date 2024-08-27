rm(list = ls())

# summarize rda files from simulations into data.frame

mydir = ".../simulation_fin"
setwd(mydir)
### parameter setting name
scenario_list = NULL
samplesize = NULL
true_para = NULL
sd_par_est = NULL
mean_sd_par = NULL
mean_point_est = NULL
mean_cindex_est = NULL
mean_cpoint_est = NULL
mean_cregin_est = NULL
covrg = NULL
sd_point_est = NULL
mean_S_est = NULL
mean_cind_est = NULL
all_conv = NULL

### point estimate and inference
allfiles = list.files()
for (i in 1:length(allfiles)) { ## for each scenario
  senario  = allfiles[i]
  subdir = paste0(mydir, "/", senario)
  setwd(subdir)
  print(subdir)
  point_est = NULL; 
  cindex_est = NULL; 
  cpoint_est = NULL; 
  cregin_est = NULL; 
  S_est = NULL; 
  S_true_est = NULL
  sd_par = NULL
  ci_par_mat = NULL; 
  ci_cpt_mat = NULL; 
  ci_crup_mat = NULL; 
  ci_crlo_mat = NULL
  conv = NULL
  CI.index = NULL
  CI.cp = NULL
  CI.crup = NULL
  CI.crlo = NULL
  bias_ls = NULL; bias_cp_ls = NULL; bias_cr_ls = NULL
  
  for (j in 1:length(list.files())) { ## for each replication
    load(list.files()[j])
    point_est = rbind(point_est, est_par)
    cindex_est  = append(cindex_est, est_cindex) 
    cpoint_est = append(cpoint_est, est_cpoint)
    cregin_est = rbind(cregin_est, est_cregin)
    sd_par = rbind(sd_par, sd_est)
    bias_ls = rbind(bias_ls,bias)
    bias_cp_ls = append(bias_cp_ls, bias_cp)
    bias_cr_ls = rbind(bias_cr_ls, bias_cr)
    ci_cpt_mat = rbind(ci_cpt_mat, ci_cpoint)
    ci_par_mat = rbind(ci_par_mat, ci_par)
    ci_crup_mat = rbind(ci_crup_mat, ci_cr_up)
    ci_crlo_mat = rbind(ci_crlo_mat, ci_cr_lo)
      
    conv = append(conv, converg)
    CI.index = rbind(CI.index, CI_index_par)
    CI.cp = append(CI.cp, CI_index_cp)
    CI.crup = append(CI.crup, CI_index_cr_up)
    CI.crlo = append(CI.crlo, CI_index_cr_lo)
    S_est = append(S_est, S)
  }
    samplesize = append(samplesize, n)
    true_para  = rbind(true_para, theta_true)
    mean_point_est = rbind(mean_point_est, apply(point_est, 2, mean))
    mean_cindex_est = rbind(mean_cindex_est, apply(cindex_est, 2, mean))
    mean_cpoint_est = rbind(mean_cpoint_est, apply(cpoint_est, 2, mean))
    mean_cregin_est = rbind(mean_cregin_est, apply(cregin_est, 2, mean))
    sd_point_est = rbind(sd_point_est, apply(point_est, 2, sd))
    mean_sd_par = rbind(mean_sd_par, apply(sd_par, 2, mean))
    covrg = rbind(covrg, apply(CI.index, 2, mean))
    all_conv = append(all_conv, mean(conv))
    mean_S_est = append(mean_S_est, mean(S_est,na.rm = TRUE))
    mean_cind_est = append(mean_cind_est, mean(cind_est))
    setwd(mydir)
}

bias = mean_point_est -true_para

mysheet = data.frame(scenario_list, true_para, mean_point_est, bias, sd_point_est, 
                     mean_sd_par[,1:4], mean_S_est, mean_cindex_est, all_conv, covrg, mean_cpoint_est)
colnames(mysheet) = c("scenario", 
                      "true value b0", "true value logb1", 
                      "true value a1", "true value a2", 
                      "mean est b0", "mean est logb1", "mean est a1", "mean est a2", 
                      "bias b0", "bias logb1", "bias a1", "bias a2",
                      "sd est b0", "sd est logb1", "sd est a1", "sd est a2", 
                      "mean sd b0", "mean sd b1", "mean sd a1", "mean sd a2",
                      "mean_est_S", "mean est C index", "convergence.index",
                      "95% coverage rate b0", "95% coverage rate logb1", 
                      "95% coverage rate a1", "95% coverage rate a2", "critical point")

# Rearrange the columns to be more intuitive
mysheet <- mysheet[ , c(1,2,6,10,14,18,3,7,11,15,19,4,8,12,16,20,5,9,13,17,21,22,
                        23,24,25,26,27,28,29)]
write.csv(mysheet, "sim_result.csv", row.names = FALSE)
